/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *
 * BP-v2RDM: a boundary-point semidefinite solver for variational 2-RDM
 *           computations.
 *
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 * 
 *@END LICENSE
 *
 * This code performs a semidefinite optimization of the electronic
 * energy according to the boundary-point algorithm described in
 * PRL 106, 083001 (2011).
 *
 */
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libqt/qt.h>

#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>

#include<libmints/wavefunction.h>
#include<libmints/mints.h>
#include<libmints/vector.h>
#include<libmints/matrix.h>
#include<../bin/fnocc/blas.h>
#include<time.h>

#include"v2rdm_solver.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace boost;
using namespace psi;
using namespace fnocc;

namespace psi{ namespace v2rdm_casscf{

void v2RDMSolver::TEI() {

    
    // one-electron integrals:  
    boost::shared_ptr<Matrix> K1 = GetOEI();

    // size of the tei buffer
    if ( is_df_ ) {

        // size of the 3-index integral buffer
        tei_full_dim_ = (long int) nQ_ * (long int) ( nmo_ - nfrzv_ ) * ( (long int) ( nmo_ - nfrzv_ ) + 1L ) / 2L ;

        // just point to 3-index integral buffer
        tei_full_sym_      = Qmo_;

    }else {

        // size of the 4-index integral buffer
        tei_full_dim_ = 0;
        for (int h = 0; h < nirrep_; h++) {
            tei_full_dim_ += (long int)gems_full[h] * ( (long int)gems_full[h] + 1L ) / 2L;
        }

        tei_full_sym_ = (double*)malloc(tei_full_dim_*sizeof(double));
        memset((void*)tei_full_sym_,'\0',tei_full_dim_*sizeof(double));

    }

    // size of d2, blocked by symmetry, including the core orbitals
    d2_plus_core_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim_ += (long int)gems_plus_core[h] * ( (long int)gems_plus_core[h] + 1L ) / 2L;
    }

    d2_plus_core_sym_  = (double*)malloc(d2_plus_core_dim_*sizeof(double));
    memset((void*)d2_plus_core_sym_,'\0',d2_plus_core_dim_*sizeof(double));


    // allocate memory for oei tensor, blocked by symmetry, excluding frozen virtuals
    oei_full_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        oei_full_dim_ += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // allocate memory for d1 tensor, blocked by symmetry, including the core orbitals
    d1_plus_core_dim_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_plus_core_dim_ += (rstcpi_[h] + frzcpi_[h] + amopi_[h]) * ( rstcpi_[h] + frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }

    oei_full_sym_ = (double*)malloc(oei_full_dim_*sizeof(double));
    memset((void*)oei_full_sym_,'\0',oei_full_dim_*sizeof(double));

    d1_plus_core_sym_ = (double*)malloc(d1_plus_core_dim_*sizeof(double));
    memset((void*)d1_plus_core_sym_,'\0',d1_plus_core_dim_*sizeof(double));

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                oei_full_sym_[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    if ( is_df_ ) {
        // build tei's from 3-index integrals 
        RepackIntegralsDF();
    }else {
        // read tei's from disk
        TEIFromDisk();
    }

}

void v2RDMSolver::TEIFromDisk() {

    double * temptei = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)temptei,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // read two-electron integrals from disk
    ReadIntegrals(temptei,nmo_);
  
    // load tei_full_sym_
    long int offset = 0;
    long int n2 = (long int)nmo_*(long int)nmo_;
    long int n3 = n2 * (long int)nmo_;
    for (int h = 0; h < nirrep_; h++) {
        for (long int ij = 0; ij < gems_full[h]; ij++) {
            long int i = bas_full_sym[h][ij][0];
            long int j = bas_full_sym[h][ij][1];
            for (long int kl = ij; kl < gems_full[h]; kl++) {
                long int k = bas_full_sym[h][kl][0];
                long int l = bas_full_sym[h][kl][1];
                tei_full_sym_[offset + INDEX(ij,kl)] = temptei[i*n3+j*n2+k*(long int)nmo_+l];
            }
        }
        offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
    }

    free(temptei);

    RepackIntegrals();
}

}}
