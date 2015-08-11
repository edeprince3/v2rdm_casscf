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

    // have three-index integrals in AO basis. now, transform to SO basis

    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral));
    boost::shared_ptr<Matrix> AO2USO_ (new Matrix(pet->aotoso()));

    long int off1 = 0;
    for (int h = 0; h < nirrep_; h++) {
        double ** cp = AO2USO_->pointer(h);

        F_DGEMM('n','n',nsopi_[h],nQ_*nso_,nso_,1.0,&cp[0][0],nsopi_[h],tmp1,nso_,0.0,tmp2+off1,nsopi_[h]);
        off1 += nQ_ * nsopi_[h] * nso_;

    }
    off1 = 0;
    long int off2 = 0;
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < nQ_; Q++) {
            for (long int i = 0; i < nsopi_[h]; i++) {
                long int ii = i + off2;
                for (long int m = 0; m < nso_; m++) {
                    tmp1[Q*nso_*nso_+ii*nso_+m] = tmp2[Q*nso_*nsopi_[h]+m*nsopi_[h] + i + off1];
                }
            }
        }
        off1 += nQ_ * nsopi_[h] * nso_;
        off2 += nsopi_[h];
    }

    off1 = 0;
    for (int h = 0; h < nirrep_; h++) {
        double ** cp = AO2USO_->pointer(h);

        F_DGEMM('n','n',nsopi_[h],nQ_*nso_,nso_,1.0,&cp[0][0],nsopi_[h],tmp1,nso_,0.0,tmp2+off1,nsopi_[h]);
        off1 += nQ_ * nsopi_[h] * nso_;

    }

    off1 = 0;
    off2 = 0;
    for (int h = 0; h < nirrep_; h++) {

        #pragma omp parallel for schedule (static)
        for (long int Q = 0; Q < nQ_; Q++) {
            for (long int m = 0; m < nso_; m++) {
                for (long int i = 0; i < nsopi_[h]; i++) {
                    long int ii = i + off2;
                    tmp1[Q*nso_*nso_+m*nso_+ii] = tmp2[Q*nso_*nsopi_[h]+m*nsopi_[h] + i + off1];
                }
            }
        }
        off1 += nQ_ * nsopi_[h] * nso_;
        off2 += nsopi_[h];
    }

    int off3 = 0;
    off1 = 0;
    for (int h = 0; h < nirrep_; h++) {
        double ** cp = Ca_->pointer(h);

        off2 = 0;
        for (int h2 = 0; h2 < nirrep_; h2++) {

            for (int Q = 0; Q < nQ_; Q++) {

                for (int m = 0; m < nsopi_[h2]; m++) {
                    int mm = m + off2;

                    double dum = 0.0;
                    for (int n = 0; n < nsopi_[h]; n++) {
                        int nn = n + off1;

                        tmp2[Q*nsopi_[h]*nsopi_[h2] + m*nsopi_[h] + n + off3] = tmp1[Q*nso_*nso_+mm*nso_+nn];

                    }
                }
            }
            off3 += nQ_*nsopi_[h]*nsopi_[h2];
            off2 += nsopi_[h2];
        }
        off1 += nsopi_[h];
    }

    off3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        double ** cp = Ca_->pointer(h);
        for (int h2 = 0; h2 < nirrep_; h2++) {

            F_DGEMM('n','n',nsopi_[h],nQ_*nsopi_[h2],nsopi_[h],1.0,&cp[0][0],nsopi_[h],tmp2+off3,nsopi_[h],0.0,tmp1+off3,nsopi_[h]);
            off3 += nQ_*nsopi_[h]*nsopi_[h2];

        }
    }

    off3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int h2 = 0; h2 < nirrep_; h2++) {

            for (int Q = 0; Q < nQ_; Q++) {
                for (int m = 0; m < nsopi_[h2]; m++) {
                    for (int i = 0; i < nsopi_[h]; i++) {
                        tmp2[Q*nsopi_[h]*nsopi_[h2] + i*nsopi_[h2] + m + off3] = tmp1[Q*nsopi_[h]*nsopi_[h2]+m*nsopi_[h] + i + off3];
                    }
                }
            }

            off3 += nQ_*nsopi_[h]*nsopi_[h2];
        }
    }

    // SO -> MO transformation:
    off3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int h2 = 0; h2 < nirrep_; h2++) {
            double ** cp = Ca_->pointer(h2);

            F_DGEMM('n','n',nsopi_[h2],nQ_*nsopi_[h],nsopi_[h2],1.0,&cp[0][0],nsopi_[h2],tmp2+off3,nsopi_[h2],0.0,tmp1+off3,nsopi_[h2]);

            off3 += nQ_*nsopi_[h]*nsopi_[h2];
        }
    }

    free(tmp2);

    Qmo_ = (boost::shared_ptr<Matrix>)(new Matrix(nso_*nso_,nQ_));
    double ** qmop = Qmo_->pointer();

    off3 = 0;
    off1 = 0;
    for (int h = 0; h < nirrep_; h++) {
        off2 = 0;
        for (int h2 = 0; h2 < nirrep_; h2++) {

            for (int Q = 0; Q < nQ_; Q++) {
                for (int i = 0; i < nmopi_[h]; i++) {
                    int ii = i + off1;

                    for (int j = 0; j < nsopi_[h2]; j++) {
                        int jj = j + off2;

                        qmop[ii*nso_+jj][Q] = tmp1[Q*nsopi_[h]*nsopi_[h2] + i * nsopi_[h2]+j + off3];
                    }
                }
            }
            off2 += nsopi_[h2];
            off3 += nQ_*nsopi_[h]*nsopi_[h2];
        }
        off1 += nsopi_[h];
    }

    free(tmp1);

    //boost::shared_ptr<MintsHelper> mints(new MintsHelper() );
    //boost::shared_ptr<Matrix> eri = mints->mo_eri(Ca_,Ca_);
    //for (int i = 0; i < nso_; i++) {
    //    for (int j = 0; j < nso_; j++) {
    //        for (int k = 0; k < nso_; k++) {
    //            for (int l = 0; l < nso_; l++) {
    //                double dum = 0.0;
    //                for (int Q = 0; Q < nQ_; Q++) {
    //                    dum += Qmo_->pointer()[i*nso_+j][Q] * Qmo_->pointer()[k*nso_+l][Q];
    //                }
    //                printf("%5i %5i %5i %5i %20.12lf %20.12lf\n",i,j,k,l,dum,eri->pointer()[i*nso_+j][k*nso_+l]);
    //            }
    //        }
    //    }
    //}
    //exit(0);

}


}}
