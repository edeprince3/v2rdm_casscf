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

#include"v2rdm_solver.h"

#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libmints/sieve.h>
#include <psifiles.h>

#include <../bin/fnocc/blas.h>

#include <libtrans/integraltransform.h>

using namespace psi;
using namespace fnocc;

namespace psi { namespace v2rdm_casscf {

void v2RDMSolver::ThreeIndexIntegrals() {

    basisset_ = reference_wavefunction_->basisset();

    // get ntri from sieve
    boost::shared_ptr<ERISieve> sieve (new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
    const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
    long int ntri = function_pairs.size();

    // read integrals that were written to disk in the scf
    nQ_ = Process::environment.globals["NAUX (SCF)"];
    if ( options_.get_str("SCF_TYPE") == "DF" ) {
        boost::shared_ptr<BasisSet> primary = BasisSet::pyconstruct_orbital(molecule_,
            "BASIS", options_.get_str("BASIS"));

        boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(molecule_,
            "DF_BASIS_SCF", options_.get_str("DF_BASIS_SCF"), "JKFIT",
            options_.get_str("BASIS"), primary->has_puream());

        nQ_ = auxiliary->nbf();
        Process::environment.globals["NAUX (SCF)"] = nQ_;
    }

    // 100 mb extra to account for all mapping arrays already 
    // allocated. this should be WAY more than necessary.
    long int extra = 100 * 1024 * 1024;  
    long int ndoubles = (memory_-extra) / 8;

    // orbitals will end up in energy order.  
    // we will want them in pitzer.  for sorting: 
    int * reorder  = (int*)malloc(nso_*sizeof(int));
    int * sym      = (int*)malloc(nso_*sizeof(int));
    bool * skip    = (bool*)malloc(nso_*sizeof(bool));

    for (int i = 0; i < nso_; i++) {
        skip[i] = false;
    }
    for (int i = 0; i < nso_; i++) {
        double min   = 1.e99;
        int count    = 0;
        int minj     = -999;
        int minh     = -999;
        int mincount = -999;
        for (int h = 0; h < nirrep_; h++) {
            for (int j = 0; j < nsopi_[h]; j++) {
                if ( skip[count+j] ) continue;
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    min      = epsilon_a_->pointer(h)[j];
                    mincount = count;
                    minj     = j;
                    minh     = h;
                }
            }
            count += nsopi_[h];
        }
        skip[mincount + minj]     = true;
        reorder[i]                = minj;
        sym[i]                    = minh;
    }

    // how many rows of (Q|mn) can we read in at once?
    if ( ndoubles < nso_*nso_ ) {
        throw PsiException("holy moses, we can't fit nso^2 doubles in memory.  increase memory!",__FILE__,__LINE__);
    }

    long int nrows = 1;
    long int rowsize = nQ_;
    while ( rowsize*nso*nso*2 > ndoubles ) {
        nrows++;
        rowsize = nQ_ / nrows;
        if (nrows * rowsize < nQ_) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ_ - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    double * tmp1 = (double*)malloc(rowdims[0]*nso_*nso_*sizeof(double));
    double * tmp2 = (double*)malloc(rowdims[0]*nso_*nso_*sizeof(double));

    long int nn1 = nso_*(nso_+1)/2;

    boost::shared_ptr<PSIO> psio(new PSIO());

    psio->open(PSIF_DCC_QSO,PSIO_OPEN_NEW);
    psio->open(PSIF_DCC_QMO,PSIO_OPEN_NEW);
    psio_address addr  = PSIO_ZERO;
    psio_address addr2 = PSIO_ZERO;
    for (int row = 0; row < nrows; row++) {
        psio->write(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * rowdims[row] * nso_ * nso_,addr,&addr);
        psio->write(PSIF_DCC_QMO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * rowdims[row] * nn1,addr2,&addr2);
    }
    psio->close(PSIF_DCC_QSO,1);
    psio->close(PSIF_DCC_QMO,1);

    // read integrals from SCF and unpack them
    addr  = PSIO_ZERO;
    addr2 = PSIO_ZERO;
    psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    for (int row = 0; row < nrows; row++) {

        // read
        psio->read(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) tmp2, sizeof(double) * ntri * rowdims[row],addr,&addr);

        // unpack
        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < rowdims[row]; Q++) {
            for (long int mn = 0; mn < ntri; mn++) {

                long int m = function_pairs[mn].first;
                long int n = function_pairs[mn].second;

                tmp1[Q*nso_*nso_+m*nso_+n] = tmp2[Q*ntri+mn];
                tmp1[Q*nso_*nso_+n*nso_+m] = tmp2[Q*ntri+mn];
            }
        }

        // write
        psio->write(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nso_*nso_ * rowdims[row],addr2,&addr2);
    }
    psio->close(PSIF_DFSCF_BJ,1);

    // AO->MO transformation matrix:
    boost::shared_ptr<Matrix> myCa (new Matrix(reference_wavefunction_->Ca_subset("AO","ALL")));

    // transform first index:
    addr  = PSIO_ZERO;
    addr2 = PSIO_ZERO;
    for (int row = 0; row < nrows; row++) {
        // read
        psio->read(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nso_*nso_ * rowdims[row],addr,&addr);

        // transform first index:
        F_DGEMM('n','n',nso_,nso_*rowdims[row],nso_,1.0,&(myCa->pointer()[0][0]),nso_,tmp1,nso_,0.0,tmp2,nso_);

        // sort
        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < rowdims[row]; Q++) {
            for (long int i = 0; i < nso_; i++) {
                for (long int m = 0; m < nso_; m++) {
                    tmp1[Q*nso_*nso_+i*nso_+m] = tmp2[Q*nso_*nso_+m*nso_+i];
                }
            }
        }

        // write
        psio->write(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nso_*nso_ * rowdims[row],addr2,&addr2);
    }
    // transform second index:
    addr  = PSIO_ZERO;
    addr2 = PSIO_ZERO;
    psio->open(PSIF_DCC_QMO,PSIO_OPEN_OLD);
    for (int row = 0; row < nrows; row++) {
        // read
        psio->read(PSIF_DCC_QSO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nso_*nso_ * rowdims[row],addr,&addr);

        // transform second index:
        F_DGEMM('n','n',nso_,nso_*rowdims[row],nso_,1.0,&(myCa->pointer()[0][0]),nso_,tmp1,nso_,0.0,tmp2,nso_);

        // sort orbitals into pitzer order
        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < rowdims[row]; Q++) {
            for (int m = 0; m < nso_; m++) {
                int hm = sym[m];
                int offm = 0;
                for (int h = 0; h < hm; h++) {
                    offm += nsopi_[h];
                }
                int mm = reorder[m] + offm;
                for (int n = 0; n < nso_; n++) {
                    int hn = sym[n];
                    int offn = 0;
                    for (int h = 0; h < hn; h++) {
                        offn += nsopi_[h];
                    }
                    int nn = reorder[n] + offn;
                    tmp1[Q*nn1+INDEX(mm,nn)] = tmp2[Q*nso_*nso_+m*nso+n];
                }
            }
        }

        // write
        psio->write(PSIF_DCC_QMO, "(Q|mn) Integrals", (char*) tmp1, sizeof(double) * nn1 * rowdims[row],addr2,&addr2);
    }
    psio->close(PSIF_DCC_QMO,1);
    psio->close(PSIF_DCC_QSO,1);

    //F_DGEMM('t','t',nso_*nQ_,nso_,nso_,1.0,tmp1,nso_,&(myCa->pointer()[0][0]),nso_,0.0,tmp2,nso_*nQ_);
    //F_DGEMM('t','t',nso_*nQ_,nso_,nso_,1.0,tmp2,nso_,&(myCa->pointer()[0][0]),nso_,0.0,tmp1,nso_*nQ_);

    free(reorder);
    free(skip);
    free(sym);
    free(tmp2);
    free(tmp1);

    Qmo_ = (double*)malloc(nn1*nQ_*sizeof(double));
    memset((void*)Qmo_,'\0',nn1*nQ_*sizeof(double));

    // TODO: read larger chunks
    double * tmp3 = (double*)malloc(nn1*sizeof(double));
    memset((void*)tmp3,'\0',nn1*sizeof(double));
    addr = PSIO_ZERO;
    psio->open(PSIF_DCC_QMO,PSIO_OPEN_OLD);
    for (long int Q = 0; Q < nQ_; Q++) {
        psio->read(PSIF_DCC_QMO,"(Q|mn) Integrals",(char*)tmp3,sizeof(double)*nn1,addr,&addr);
        for (long int mn = 0; mn < nn1; mn++) {
            Qmo_[mn*nQ_+Q] = tmp3[mn];
        }
    }
    psio->close(PSIF_DCC_QMO,1);

    free(tmp3);
}


}}
