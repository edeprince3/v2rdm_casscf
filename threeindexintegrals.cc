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

namespace psi { namespace v2rdm_casscf {


void v2RDMSolver::ThreeIndexIntegrals() {

    basisset_ = reference_wavefunction_->basisset();

    // get ntri from sieve
    boost::shared_ptr<ERISieve> sieve (new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
    const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
    long int ntri = function_pairs.size();

    // read integrals that were written to disk in the scf

    boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(molecule_,
          "DF_BASIS_SCF", options_.get_str("DF_BASIS_SCF"), "JKFIT", options_.get_str("BASIS"));
    nQ_ = auxiliary->nbf();

    boost::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",nQ_,ntri));
    double** Qmnp = Qmn->pointer();
    boost::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * nQ_);
    psio->close(PSIF_DFSCF_BJ,1);

    // have three-index integrals in AO basis. now, transform to SO basis

    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral));
    boost::shared_ptr<Matrix> AO2USO_ (new Matrix(pet->aotoso()));

    boost::shared_ptr<Matrix> Qso (new Matrix(nQ_,nso_*nso_) );

    for (int Q = 0; Q < nQ_; Q++) {
        int offh = 0;
        for (int h = 0; h < nirrep_; h++) {
            double ** c1 = AO2USO_->pointer(h);
            for (int i = 0; i < nsopi_[h]; i++) {
                int ii = i + offh;
                int offh2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
                    double ** c2 = AO2USO_->pointer(h2);
                    for (int j = 0; j < nsopi_[h2]; j++) {
                        int jj = j + offh2;

                        double dum = 0.0;
                        for (int mn = 0; mn < ntri; mn++) {
                            long int m = function_pairs[mn].first;
                            long int n = function_pairs[mn].second;
                            dum += Qmnp[Q][mn] * c1[m][i] * c2[n][j];
                            if ( m != n ) {
                                dum += Qmnp[Q][mn] * c1[m][i] * c2[n][j];
                            }
                        }
                        Qso->pointer()[Q][ii*nso_+jj] = dum;

                    }
                    offh2 += nsopi_[h2];
                }
            }
            offh += nsopi_[h];
        }
    }
    Qso->print();

    // SO -> MO transformation:
    Qmo_ = (boost::shared_ptr<Matrix>)(new Matrix(nQ_,nso_*nso_));

    for (int Q = 0; Q < nQ_; Q++) {
        int offh = 0;
        for (int h = 0; h < nirrep_; h++) {
            double ** c1 = Ca_->pointer(h);
            for (int i = 0; i < nmopi_[h]; i++) {
                int ii = i + offh;
                int offh2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
                    double ** c2 = Ca_->pointer(h2);
                    for (int j = 0; j < nmopi_[h2]; j++) {
                        int jj = j + offh2;

                        double dum = 0.0;
                        for (int m = 0; m < nsopi_[h]; m++) {
                            int mm = m + offh;
                            for (int n = 0; n < nsopi_[h2]; n++) {
                                int nn = n + offh2;
                                dum += Qso->pointer()[Q][mm*nso_+nn] * c1[m][i] * c2[n][j];
                            }
                        }
                        Qmo_->pointer()[Q][ii*nso_+jj] = dum;
                    }
                    offh2 += nsopi_[h2];
                }
            }
            offh += nsopi_[h];
        }
    }
}


}}
