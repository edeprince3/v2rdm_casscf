/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
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
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 * 
 *@END LICENSE
 *
 */

#ifndef _python_api_h_
#define _python_api_h_

#include "v2rdm_solver.h"

//#include <pybind11/pybind11.h>
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>

#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/wavefunction.h"

using namespace psi;

//namespace py = pybind11;
//using namespace pybind11::literals;

namespace psi { namespace v2rdm_casscf {

/*
void export_v2RDMSolver(py::module& m) {
    py::class_<v2rdm_casscf::v2RDMSolver, std::shared_ptr<v2rdm_casscf::v2RDMSolver>, Wavefunction>(m, "v2RDMSolver")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("set_orbitals", &v2RDMSolver::set_orbitals)
        .def("get_orbitals", &v2RDMSolver::get_orbitals)
        .def("get_opdm", &v2RDMSolver::get_opdm)
        .def("get_tpdm", &v2RDMSolver::get_tpdm)
        .def("compute_energy", &v2RDMSolver::compute_energy);
}

PYBIND11_MODULE(v2rdm_casscf, m) {
    m.doc() = "Python API of v2rdm_casscf: a variational 2-RDM-driven CASSCF plugin to Psi4";
    export_v2RDMSolver(m);
}
*/

std::shared_ptr<Matrix> v2RDMSolver::get_opdm() {

    Dimension amopi = Dimension(nirrep_);
    for (int h = 0; h < nirrep_; h++) {
        amopi[h] = amopi_[h];
    }
    std::shared_ptr<Matrix> opdm(new Matrix("OPDM",amopi,amopi));

    double * x_p = x->pointer();
    for (int h = 0; h < nirrep_; h++) {
        double ** opdm_p = opdm->pointer(h);
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                opdm_p[i][j]  = x_p[d1aoff[h]+i*amopi_[h]+j];
                opdm_p[i][j] += x_p[d1aoff[h]+i*amopi_[h]+j];
            }
        }
    }

    return opdm;
}

std::shared_ptr<Matrix> v2RDMSolver::get_tpdm() {

    if ( nirrep_ > 1 ) {
        throw PsiException("not sure how CIWavefunction stores TPDM with symmetry",__FILE__,__LINE__);
    }
    //Dimension amopi = Dimension(nirrep_);
    //for (int h = 0; h < nirrep_; h++) {
    //    amopi[h] = amopi_[h];
    //}
    std::shared_ptr<Matrix> tpdm(new Matrix("TPDM",amo_*amo_,amo_*amo_));

    double ** tpdm_p = tpdm->pointer();
    double * x_p = x->pointer();
    for (int i = 0; i < amo_; i++) {
        for (int j = 0; j < amo_; j++) {
            int ij_ab = ibas_ab_sym[0][i][j];
            int ji_ab = ibas_ab_sym[0][j][i];
            int ij_aa = ibas_aa_sym[0][i][j];
            for (int k = 0; k < amo_; k++) {
                for (int l = 0; l < amo_; l++) {
                    int kl_ab = ibas_ab_sym[0][k][l];
                    int lk_ab = ibas_ab_sym[0][l][k];
                    int il_ab = ibas_ab_sym[0][i][l];
                    int kj_ab = ibas_ab_sym[0][k][j];
                    int li_ab = ibas_ab_sym[0][l][i];
                    int jk_ab = ibas_ab_sym[0][j][k];

                    int kl_aa = ibas_aa_sym[0][k][l];
                    int il_aa = ibas_aa_sym[0][i][l];
                    int kj_aa = ibas_aa_sym[0][k][j];

                    double dum_ikjl = x_p[d2aboff[0] + ij_ab * gems_ab[0] + kl_ab];
                    dum_ikjl       += x_p[d2aboff[0] + ji_ab * gems_ab[0] + lk_ab];

                    double dum_iklj = x_p[d2aboff[0] + il_ab * gems_ab[0] + kj_ab];
                    dum_iklj       += x_p[d2aboff[0] + li_ab * gems_ab[0] + jk_ab];

                    if ( i != j && k != l ) {
                        int sg = 1;
                        if ( i > j ) sg = -sg;
                        if ( k > l ) sg = -sg;
                        dum_ikjl += sg * x_p[d2aaoff[0] + ij_aa * gems_aa[0] + kl_aa];
                        dum_ikjl += sg * x_p[d2bboff[0] + ij_aa * gems_aa[0] + kl_aa];
                    }
                    if ( i != l && k != j ) {
                        int sg = 1;
                        if ( i > l ) sg = -sg;
                        if ( k > j ) sg = -sg;
                        dum_ikjl += sg * x_p[d2aaoff[0] + il_aa * gems_aa[0] + kj_aa];
                        dum_ikjl += sg * x_p[d2bboff[0] + il_aa * gems_aa[0] + kj_aa];
                    }

                    //tpdm_p[i*amo_+j][k*amo_+l] =  dum;
                    tpdm_p[i*amo_+k][j*amo_+l] =  0.5 * (dum_ikjl + dum_iklj);
                    //tpdm_p[i*amo_+j][k*amo_+l] =  0.5 * (dum_ijkl + dum_ijlk);
                }
            }
        }
    }
    return tpdm;
}

void v2RDMSolver::orbital_locations(const std::string& orbitals, int* start, int* end) {
    if (orbitals == "FROZEN_DOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = frzcpi_[h];
        }
    } else if (orbitals == "RESTRICTED_DOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h];
        }
    } else if (orbitals == "DOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = frzcpi_[h] + rstcpi_[h];
        }
    } else if (orbitals == "ACTIVE") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h] + rstcpi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h];
        }
    } else if (orbitals == "RESTRICTED_UOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h] + rstvpi_[h];
        }
    } else if (orbitals == "FROZEN_UOCC") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h] + rstvpi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h] + rstvpi_[h] + frzvpi_[h];
        }
    } else if (orbitals == "VIRTUAL") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h];
            end[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h] + rstvpi_[h] + frzvpi_[h];
        }

    } else if (orbitals == "ALL") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = 0;
            end[h] = nmopi_[h];
        }
    } else if (orbitals == "ROT") {
        for (int h = 0; h < nirrep_; h++) {
            start[h] = frzcpi_[h];
            end[h] = nmopi_[h] - frzvpi_[h];
        }
    } else {
        throw PSIEXCEPTION(
            "v2RDMSolver: Orbital subset is not defined, should be FROZEN_DOCC, "
            "RESTRICTED_DOCC, DOCC, ACTIVE, RESTRICTED_UOCC, FROZEN_UOCC, "
            "VIRTUAL, ROT, or ALL");
    }
}

std::shared_ptr<Matrix> v2RDMSolver::get_orbitals(const std::string& orbital_name) {
    /// Figure out orbital positions
    auto* start = new int[nirrep_];
    auto* end = new int[nirrep_];

    orbital_locations(orbital_name, start, end);

    auto* spread = new int[nirrep_];
    for (int h = 0; h < nirrep_; h++) {
        spread[h] = end[h] - start[h];
    }

    /// Fill desired orbitals
    auto retC = std::make_shared<Matrix>("C " + orbital_name, nirrep_, nsopi_, spread);
    for (int h = 0; h < nirrep_; h++) {
        for (int i = start[h], pos = 0; i < end[h]; i++, pos++) {
            C_DCOPY(nsopi_[h], &Ca_->pointer(h)[0][i], nmopi_[h], &retC->pointer(h)[0][pos], spread[h]);
        }
    }

    /// Cleanup
    delete[] start;
    delete[] end;
    delete[] spread;

    return retC;
}

void v2RDMSolver::set_orbitals(const std::string& orbital_name, SharedMatrix orbitals) {
    /// Figure out orbital positions
    auto* start = new int[nirrep_];
    auto* end = new int[nirrep_];

    orbital_locations(orbital_name, start, end);

    auto* spread = new int[nirrep_];
    for (int h = 0; h < nirrep_; h++) {
        spread[h] = end[h] - start[h];
    }

    /// Fill desired orbitals
    for (int h = 0; h < nirrep_; h++) {
        for (int i = start[h], pos = 0; i < end[h]; i++, pos++) {
            C_DCOPY(nsopi_[h], &orbitals->pointer(h)[0][pos], spread[h], &Ca_->pointer(h)[0][i], nmopi_[h]);
        }
    }

    /// Cleanup
    delete[] start;
    delete[] end;
    delete[] spread;
}


}} // End namespaces

#endif
