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
        boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(molecule_,
              "DF_BASIS_SCF", options_.get_str("DF_BASIS_SCF"), "JKFIT", options_.get_str("BASIS"));
        nQ_ = auxiliary->nbf();
        Process::environment.globals["NAUX (SCF)"] = nQ_;
    }

    double * tmp1 = (double*)malloc(nQ_*nso_*nso_*sizeof(double));
    double * tmp2 = (double*)malloc(nQ_*nso_*nso_*sizeof(double));

    boost::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) tmp2, sizeof(double) * ntri * nQ_);
    psio->close(PSIF_DFSCF_BJ,1);

    // unpack

    #pragma omp parallel for schedule (static)
    for (long int Q = 0; Q < nQ_; Q++) {
        for (long int mn = 0; mn < ntri; mn++) {

            long int m = function_pairs[mn].first;
            long int n = function_pairs[mn].second;

            tmp1[Q*nso_*nso_+m*nso_+n] = tmp2[Q*ntri+mn];
            tmp1[Q*nso_*nso_+n*nso_+m] = tmp2[Q*ntri+mn];
        }
    }

    boost::shared_ptr<Matrix> myCa_ (new Matrix(reference_wavefunction_->Cb_subset("AO","ALL")));

    F_DGEMM('t','t',nso_*nQ_,nso_,nso_,1.0,tmp1,nso_,&(myCa_->pointer()[0][0]),nso_,0.0,tmp2,nso_*nQ_);
    F_DGEMM('t','t',nso_*nQ_,nso_,nso_,1.0,tmp2,nso_,&(myCa_->pointer()[0][0]),nso_,0.0,tmp1,nso_*nQ_);

    free(tmp2);

    // orbitals are in energy order.  we want them in pitzer

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

    Qmo_ = (double*)malloc(nso_*nso_*nQ_*sizeof(double));
    memset((void*)Qmo_,'\0',nso_*nso_*nQ_*sizeof(double));

    // sort integrals: (Q|mn) -> (Q|m'n') mn are energy order, m'n' are pitzer order
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
            //C_DCOPY(nQ_,qmop[mm*nso_+nn], 1 ,tmp1 + m*nQ_*nso_+n*nQ_,1);
            C_DCOPY(nQ_,tmp1 + nQ_*(m*nso_+n), 1 ,Qmo_ + nQ_*INDEX(mm,nn),1);
        }
    }

    free(tmp1);
    free(reorder);
    free(skip);
    free(sym);
}


}}
