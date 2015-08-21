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


//build subset of tei's from 3-index integrals 
void v2RDMSolver::DF_TEI() {

    // one-electron integrals:  
    boost::shared_ptr<Matrix> K1 = GetOEI();

    // size of the 3-index integral buffer
    tei_full_dim_ = (long int) nQ_ * (long int) nmo_ * ( (long int) nmo_ + 1 ) /2 ;

    d2_plus_core_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim_ += gems_plus_core[h] * ( gems_plus_core[h] + 1 ) / 2;
    }

    // just point to 3-index integral buffer
    tei_full_sym_      = Qmo_;

    d2_plus_core_sym_  = (double*)malloc(d2_plus_core_dim_*sizeof(double));
    memset((void*)d2_plus_core_sym_,'\0',d2_plus_core_dim_*sizeof(double));

    // allocate memory for full OEI tensor blocked by symmetry for Greg
    oei_full_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        oei_full_dim_ += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
    d1_plus_core_dim_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_plus_core_dim_ += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }
    oei_full_sym_ = (double*)malloc(oei_full_dim_*sizeof(double));
    d1_plus_core_sym_  = (double*)malloc(d1_plus_core_dim_*sizeof(double));
    memset((void*)oei_full_sym_,'\0',oei_full_dim_*sizeof(double));
    memset((void*)d1_plus_core_sym_,'\0',d1_plus_core_dim_*sizeof(double));

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h]; i++) {
            for (long int j = i; j < nmopi_[h]; j++) {
                oei_full_sym_[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
    RepackIntegralsDF();

    enuc_ = Process::environment.molecule()->nuclear_repulsion_energy();

}

void v2RDMSolver::TEI() {

    double * temptei = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)temptei,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // read two-electron integrals from disk
    ReadIntegrals(temptei,nmo_);
   
    // read oei's from disk (this will simplify things when we add symmetry):
    boost::shared_ptr<Matrix> K1 = GetOEI();

    // allocate memory for full ERI tensor blocked by symmetry for Greg
    tei_full_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        tei_full_dim_ += (long int)gems_full[h] * ( (long int)gems_full[h] + 1L ) / 2L;
    }
    d2_plus_core_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim_ += (long int)gems_plus_core[h] * ( (long int)gems_plus_core[h] + 1L ) / 2L;
    }
   
    tei_full_sym_      = (double*)malloc(tei_full_dim_*sizeof(double));
    d2_plus_core_sym_  = (double*)malloc(d2_plus_core_dim_*sizeof(double));
    memset((void*)tei_full_sym_,'\0',tei_full_dim_*sizeof(double));
    memset((void*)d2_plus_core_sym_,'\0',d2_plus_core_dim_*sizeof(double));
    //for (int i = 0; i < d2_plus_core_dim_; i++) {
    //    d2_plus_core_sym_[i] = -9.0e99;
    //}
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
    // allocate memory for full OEI tensor blocked by symmetry for Greg
    oei_full_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        oei_full_dim_ += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
    d1_plus_core_dim_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_plus_core_dim_ += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }
    oei_full_sym_ = (double*)malloc(oei_full_dim_*sizeof(double));
    d1_plus_core_sym_  = (double*)malloc(d1_plus_core_dim_*sizeof(double));
    memset((void*)oei_full_sym_,'\0',oei_full_dim_*sizeof(double));
    memset((void*)d1_plus_core_sym_,'\0',d1_plus_core_dim_*sizeof(double));

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h]; i++) {
            for (long int j = i; j < nmopi_[h]; j++) {
                oei_full_sym_[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
   
    // if frozen core, adjust oei's and compute frozen core energy:
    efzc_ = 0.0;
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < frzcpi_[h]; i++) {
            efzc_ += 2.0 * K1->pointer(h)[i][i];

        long int offset2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
          for (long int j = 0; j < frzcpi_[h2]; j++) {
                    efzc_ += (2.0 * temptei[(i+offset)*n3+(i+offset) *n2+(j+offset2)*(long int)nmo_+(j+offset2)]
                                  - temptei[(i+offset)*n3+(j+offset2)*n2+(i+offset) *(long int)nmo_+(j+offset2)]);
                }
                offset2 += nmopi_[h2];
            }
        }
        offset += nmopi_[h];
    }

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = frzcpi_[h]; i < nmopi_[h] - frzvpi_[h]; i++) {
            for (long int j = frzcpi_[h]; j < nmopi_[h] - frzvpi_[h]; j++) {

                double dum = 0.0;

                int offset2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {

                    for (long int k = 0; k < frzcpi_[h2]; k++) {
                    dum += (2.0 * temptei[(i+offset)*n3+(j+offset) *n2+(k+offset2)*(long int)nmo_+(k+offset2)]
                                - temptei[(i+offset)*n3+(k+offset2)*n2+(k+offset2)*(long int)nmo_+(j+offset)]);
                    }

                    offset2 += nmopi_[h2];
                }
                K1->pointer(h)[i][j] += dum;
            }
        }
        offset += nmopi_[h];
    }

    double* c_p = c->pointer();
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < amopi_[h]; i++) {
            for (long int j = 0; j < amopi_[h]; j++) {
                //c_p[d1aoff[h] + INDEX(i-frzcpi_[h],j-frzcpi_[h])] = K1->pointer(h)[i][j];
                //c_p[d1boff[h] + INDEX(i-frzcpi_[h],j-frzcpi_[h])] = K1->pointer(h)[i][j];
                //c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = K1->pointer(h)[i][j];
                //c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = K1->pointer(h)[i][j];

                c_p[d1aoff[h] + i*amopi_[h] + j] = K1->pointer(h)[i+frzcpi_[h]][j+frzcpi_[h]];
                c_p[d1boff[h] + i*amopi_[h] + j] = K1->pointer(h)[i+frzcpi_[h]][j+frzcpi_[h]];
            }
        }
    }

    long int na = nalpha_ - nfrzc_;
    long int nb = nbeta_ - nfrzc_;
    for (int h = 0; h < nirrep_; h++) {
        for (long int ij = 0; ij < gems_ab[h]; ij++) {
            long int i = bas_ab_sym[h][ij][0];
            long int j = bas_ab_sym[h][ij][1];
            long int ii = full_basis[i];
            long int jj = full_basis[j];

            for (long int kl = 0; kl < gems_ab[h]; kl++) {
                long int k = bas_ab_sym[h][kl][0];
                long int l = bas_ab_sym[h][kl][1];

                long int kk = full_basis[k];
                long int ll = full_basis[l];

                c_p[d2aboff[h] + ij*gems_ab[h]+kl]    = temptei[ii*n3+kk*n2+jj*(long int)nmo_+ll];

                long int hi = symmetry[i];
                long int hj = symmetry[j];
            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (long int ij = 0; ij < gems_aa[h]; ij++) {
            long int i = bas_aa_sym[h][ij][0];
            long int j = bas_aa_sym[h][ij][1];

            long int ii = full_basis[i];
            long int jj = full_basis[j];

            for (long int kl = 0; kl < gems_aa[h]; kl++) {
                long int k = bas_aa_sym[h][kl][0];
                long int l = bas_aa_sym[h][kl][1];

                long int kk = full_basis[k];
                long int ll = full_basis[l];

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]    = 0.5 * temptei[ii*n3+kk*n2+jj*(long int)nmo_+ll];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * temptei[ii*n3+ll*n2+jj*(long int)nmo_+kk];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * temptei[jj*n3+kk*n2+ii*(long int)nmo_+ll];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   += 0.5 * temptei[jj*n3+ll*n2+ii*(long int)nmo_+kk];

                long int hi = symmetry[i];
                long int hj = symmetry[j];

                c_p[d2bboff[h] + ij*gems_aa[h]+kl]    = 0.5 * temptei[ii*n3+kk*n2+jj*(long int)nmo_+ll];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * temptei[ii*n3+ll*n2+jj*(long int)nmo_+kk];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * temptei[jj*n3+kk*n2+ii*(long int)nmo_+ll];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   += 0.5 * temptei[jj*n3+ll*n2+ii*(long int)nmo_+kk];

            }
        }
    }

    enuc_ = Process::environment.molecule()->nuclear_repulsion_energy();
    free(temptei);
}

}}
