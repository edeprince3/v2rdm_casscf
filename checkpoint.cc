/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF by A. Eugene DePrince III, a plugin to:
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

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>

#include<psi4/libmints/wavefunction.h>
#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>
#include<time.h>

#include"v2rdm_solver.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace psi;

namespace psi{ namespace v2rdm_casscf{


void v2RDMSolver::InitializeCheckpointFile() {

    std::shared_ptr<PSIO> psio ( new PSIO() );
    psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_NEW);

    // scf energy
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"SCF ENERGY",(char*)(&escf_),sizeof(double));

    // number of irreps
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"NIRREP",(char*)(&nirrep_),sizeof(int));

    // is_df_
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"IS DF?",(char*)(&is_df_),sizeof(bool));

    // mu
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"MU",(char*)(&mu),sizeof(double));

    // x
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"PRIMAL",(char*)x->pointer(),dimx_*sizeof(double));

    // y
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"DUAL 1",(char*)y->pointer(),nconstraints_*sizeof(double));

    // z
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"DUAL 2",(char*)z->pointer(),dimx_*sizeof(double));

    // one-electron integrals
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"OEI",(char*)oei_full_sym_,oei_full_dim_*sizeof(double));

    // two-electron integrals (or DF/CD integrals)
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"TEI",(char*)tei_full_sym_,tei_full_dim_*sizeof(double));

    // orbital optimization transformation matrix
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"TRANSFORMATION MATRIX",
        (char*)orbopt_transformation_matrix_,(nmo_-nfrzc_-nfrzv_)*(nmo_-nfrzc_-nfrzv_)*sizeof(double));

    // energy order to pitzer order mapping array
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"ENERGY_TO_PITZER_ORDER",
        (char*)energy_to_pitzer_order,(nmo_-nfrzv_)*sizeof(int));

    // energy order to pitzer order mapping array
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"ENERGY_TO_PITZER_ORDER_REALLY_FULL",
        (char*)energy_to_pitzer_order_really_full,nmo_*sizeof(int));

    // orbital symmetries (energy order)
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"SYMMETRY_ENERGY_ORDER",
        (char*)symmetry_energy_order,(nmo_-nfrzv_)*sizeof(int));

    psio->close(PSIF_V2RDM_CHECKPOINT,1);
}
void v2RDMSolver::WriteCheckpointFile() {

    std::shared_ptr<PSIO> psio ( new PSIO() );
    psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_OLD);

    // mu
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"MU",(char*)(&mu),sizeof(double));

    // x
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"PRIMAL",(char*)x->pointer(),dimx_*sizeof(double));

    // y
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"DUAL 1",(char*)y->pointer(),nconstraints_*sizeof(double));

    // z
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"DUAL 2",(char*)z->pointer(),dimx_*sizeof(double));

    // one-electron integrals
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"OEI",(char*)oei_full_sym_,oei_full_dim_*sizeof(double));

    // two-electron integrals (or DF/CD integrals)
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"TEI",(char*)tei_full_sym_,tei_full_dim_*sizeof(double));

    // orbital optimization transformation matrix
    psio->write_entry(PSIF_V2RDM_CHECKPOINT,"TRANSFORMATION MATRIX",
        (char*)orbopt_transformation_matrix_,(nmo_-nfrzc_-nfrzv_)*(nmo_-nfrzc_-nfrzv_)*sizeof(double));

    psio->close(PSIF_V2RDM_CHECKPOINT,1);
}
void v2RDMSolver::ReadFromCheckpointFile() {

    std::shared_ptr<PSIO> psio ( new PSIO() );
    psio->open(PSIF_V2RDM_CHECKPOINT,PSIO_OPEN_OLD);

    double dume;

    // scf energy
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"SCF ENERGY",(char*)(&dume),sizeof(double));
    if ( fabs(dume - escf_) > 1e-8 ) {
        throw PsiException("CHECKPOINT and current SCF energies do not agree",__FILE__,__LINE__);
    }
    escf_ = dume;

    int dumn;

    // number of irreps
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"NIRREP",(char*)(&dumn),sizeof(int));
    if ( dumn != nirrep_) {
        throw PsiException("CHECKPOINT and current number of irreps do not agree",__FILE__,__LINE__);
    }

    // is_df_
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"IS DF?",(char*)(&is_df_),sizeof(bool));
    if ( dumn != nirrep_) {
        throw PsiException("CHECKPOINT and current integral type do not agree",__FILE__,__LINE__);
    }

    // mu
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"MU",(char*)(&mu),sizeof(double));

    // x
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"PRIMAL",(char*)x->pointer(),dimx_*sizeof(double));

    // y
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"DUAL 1",(char*)y->pointer(),nconstraints_*sizeof(double));

    // z
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"DUAL 2",(char*)z->pointer(),dimx_*sizeof(double));

    // one-electron integrals
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"OEI",(char*)oei_full_sym_,oei_full_dim_*sizeof(double));

    // two-electron integrals (or DF/CD integrals)
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"TEI",(char*)tei_full_sym_,tei_full_dim_*sizeof(double));

    // orbital optimization transformation matrix
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"TRANSFORMATION MATRIX",
        (char*)orbopt_transformation_matrix_,(nmo_-nfrzc_-nfrzv_)*(nmo_-nfrzc_-nfrzv_)*sizeof(double));

    // energy order to pitzer order mapping array
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"ENERGY_TO_PITZER_ORDER",
        (char*)energy_to_pitzer_order,(nmo_-nfrzv_)*sizeof(int));

    // energy order to pitzer order mapping array
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"ENERGY_TO_PITZER_ORDER_REALLY_FULL",
        (char*)energy_to_pitzer_order_really_full,nmo_*sizeof(int));

    // orbital symmetries (energy order)
    psio->read_entry(PSIF_V2RDM_CHECKPOINT,"SYMMETRY_ENERGY_ORDER",
        (char*)symmetry_energy_order,(nmo_-nfrzv_)*sizeof(int));

    psio->close(PSIF_V2RDM_CHECKPOINT,1);

    RepackIntegrals();
}

}}
