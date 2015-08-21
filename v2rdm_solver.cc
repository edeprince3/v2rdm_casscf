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

#define hey(n) printf("hey%i\n",n);fflush(stdout);

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <libmints/writer.h>
#include <libmints/writer_file_prefix.h>

#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libqt/qt.h>

#include<libpsio/psio.hpp>
#include<libmints/wavefunction.h>
#include<psifiles.h>
#include<libpsio/psio.hpp>
#include<libmints/mints.h>
#include<libmints/vector.h>
#include<libmints/matrix.h>
#include<../bin/fnocc/blas.h>
#include<time.h>
#include <../bin/fnocc/blas.h>

#include <libiwl/iwl.h>

#include"cg_solver.h"
#include"v2rdm_solver.h"
#include"oeprop.h"

// greg
#include "fortran.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace boost;
using namespace psi;
using namespace fnocc;

static void evaluate_Ap(long int n, SharedVector Ax, SharedVector x, void * data) {
  
    // reinterpret void * as an instance of v2RDMSolver
    v2rdm_casscf::v2RDMSolver* BPSDPcg = reinterpret_cast<v2rdm_casscf::v2RDMSolver*>(data);
    // call a function from class to evaluate Ax product:
    BPSDPcg->cg_Ax(n,Ax,x);

}
namespace psi{ namespace v2rdm_casscf{

v2RDMSolver::v2RDMSolver(boost::shared_ptr<Wavefunction> reference_wavefunction,Options & options):
    Wavefunction(options,_default_psio_lib_){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

v2RDMSolver::~v2RDMSolver()
{
}

void  v2RDMSolver::common_init(){

    is_df_ = false;
    if ( options_.get_str("SCF_TYPE") == "DF" || options_.get_str("SCF_TYPE") == "CD" ) {
        is_df_ = true;
    }

    // This function is a copy of common_init() from SDPSolver
    escf_      = reference_wavefunction_->reference_energy();
    nalpha_   = reference_wavefunction_->nalpha();
    nbeta_    = reference_wavefunction_->nbeta();
    nalphapi_ = reference_wavefunction_->nalphapi();
    nbetapi_  = reference_wavefunction_->nbetapi();
    doccpi_   = reference_wavefunction_->doccpi();
    soccpi_   = reference_wavefunction_->soccpi();
    frzcpi_   = reference_wavefunction_->frzcpi();
    frzvpi_   = reference_wavefunction_->frzvpi();
    nmopi_    = reference_wavefunction_->nmopi();
    nirrep_   = reference_wavefunction_->nirrep();
    nso_      = reference_wavefunction_->nso();
    nmo_      = reference_wavefunction_->nmo();
    nsopi_    = reference_wavefunction_->nsopi();

    // multiplicity:
    multiplicity_ = Process::environment.molecule()->multiplicity();

    if (options_["FROZEN_DOCC"].has_changed()) {
        if (options_["FROZEN_DOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzcpi_[h] = options_["FROZEN_DOCC"][h].to_double();
        }
    }
    if (options_["FROZEN_UOCC"].has_changed()) {
        if (options_["FROZEN_UOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzvpi_[h] = options_["FROZEN_UOCC"][h].to_double();
        }
    }


    // were there linear dependencies in the primary basis set?
    if ( nmo_ != nso_ ) {

        // which irreps lost orbitals?
        int * lost = (int*)malloc(nirrep_*sizeof(int));
        memset((void*)lost,'\0',nirrep_*sizeof(int));
        bool active_space_changed = false;
        for (int h = 0; h < factory_->nirrep(); h++){
            lost[h] = nsopi_[h] - nmopi_[h];
            if ( lost[h] > 0 ) {
                active_space_changed = true;
            }

            if ( frzvpi_[h] > 0 && lost[h] > 0 ) {
                frzvpi_[h] -= ( frzvpi_[h] < lost[h] ? frzvpi_[h] : lost[h] );
            }
        }
        if ( active_space_changed ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING!! >>>\n");
            outfile->Printf("\n");
            outfile->Printf("    Your basis set may have linear dependencies.\n");
            outfile->Printf("    The number of frozen virtual orbitals per irrep were adjusted accordingly:\n");
            outfile->Printf("\n");
            outfile->Printf("    No. orbitals removed per irrep: [");
            for (int h = 0; h < nirrep_; h++) 
                outfile->Printf("%4i",lost[h]);
            outfile->Printf(" ]\n");
            outfile->Printf("    No. frozen virtuals per irrep:  [");
            for (int h = 0; h < nirrep_; h++) 
                outfile->Printf("%4i",frzvpi_[h]);
            outfile->Printf(" ]\n");
            outfile->Printf("\n");
            outfile->Printf("    Check that your active space is still correct.\n");
            outfile->Printf("\n");
        }
       
        //if ( is_df_ ) { 
        //    throw PsiException("v2rdm doesn't work when scf_type = df and linear dependencies in basis",__FILE__,__LINE__);
        //}
    }



    // active orbitals per irrep:
    amopi_    = (int*)malloc(nirrep_*sizeof(int));
    
    S_  = SharedMatrix(reference_wavefunction_->S());
    Da_ = SharedMatrix(reference_wavefunction_->Da());
    Ca_ = SharedMatrix(reference_wavefunction_->Ca());
    Fa_ = SharedMatrix(reference_wavefunction_->Fa());
    Db_ = SharedMatrix(reference_wavefunction_->Db());
    Cb_ = SharedMatrix(reference_wavefunction_->Cb());
    Fb_ = SharedMatrix(reference_wavefunction_->Fb());
    
    //Ca_->print();

    epsilon_a_= boost::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= boost::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());
    
    amo_      = 0;
    nfrzc_    = 0;

    int ndocc = 0;
    int nvirt = 0;
    int nfrzv = 0;
    for (int h = 0; h < nirrep_; h++){
        nfrzc_   += frzcpi_[h];
        nfrzv    += frzvpi_[h];
        amo_   += nmopi_[h]-frzcpi_[h]-frzvpi_[h];
        ndocc    += doccpi_[h];
        amopi_[h] = nmopi_[h]-frzcpi_[h]-frzvpi_[h];
    }

    int ndoccact = ndocc - nfrzc_;
    nvirt    = amo_ - ndoccact;

    // sanity check for orbital occupancies:
    for (int h = 0; h < nirrep_; h++) {
        int tot = doccpi_[h] + soccpi_[h] + frzvpi_[h];
        if (doccpi_[h] + soccpi_[h] + frzvpi_[h] > nmopi_[h] ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING >>> irrep %5i has too many orbitals:\n",h);
            outfile->Printf("\n");
            outfile->Printf("                    docc = %5i\n",doccpi_[h]);
            outfile->Printf("                    socc = %5i\n",soccpi_[h]);
            outfile->Printf("                    frzv = %5i\n",frzvpi_[h]);
            outfile->Printf("                    tot  = %5i\n",doccpi_[h] + soccpi_[h] + frzvpi_[h]);
            outfile->Printf("\n");
            outfile->Printf("                    total no. orbitals should be %5i\n",nmopi_[h]);
            outfile->Printf("\n");
            throw PsiException("at least one irrep has too many orbitals",__FILE__,__LINE__);
        }
        if (frzcpi_[h] > doccpi_[h] ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING >>> irrep %5i has too many frozen core orbitals:\n",h);
            outfile->Printf("                    frzc = %5i\n",frzcpi_[h]);
            outfile->Printf("                    docc = %5i\n",doccpi_[h]);
            outfile->Printf("\n");
            throw PsiException("at least one irrep has too many frozen core orbitals",__FILE__,__LINE__);
        }
    }
    
    // memory is from process::environment        
    memory_ = Process::environment.get_memory();
    // set the wavefunction name
    name_ = "V2RDM CASSCF";

    // pick conditions.  default is dqg
    constrain_q2_ = true;
    constrain_g2_ = true;
    constrain_t1_ = false;
    constrain_t2_ = false;
    constrain_d3_ = false;
    if (options_.get_str("POSITIVITY")=="D") {
        constrain_q2_ = false;
        constrain_g2_ = false;
    }else if (options_.get_str("POSITIVITY")=="DQ") {
        constrain_q2_ = true;
        constrain_g2_ = false;
    }else if (options_.get_str("POSITIVITY")=="DG") {
        constrain_q2_ = false;
        constrain_g2_ = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT1") {
        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT2") {
        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t2_ = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT1T2") {
        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;
        constrain_t2_ = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT") {
        constrain_q2_ = true;
        constrain_g2_ = true;
        constrain_t1_ = true;
        constrain_t2_ = true;
    }

    if ( options_.get_bool("CONSTRAIN_D3") ) {
        constrain_d3_ = true;
    }

    spin_adapt_g2_  = options_.get_bool("SPIN_ADAPT_G2");
    spin_adapt_q2_  = options_.get_bool("SPIN_ADAPT_Q2");
    constrain_spin_ = options_.get_bool("CONSTRAIN_SPIN");

    if ( constrain_t1_ || constrain_t2_ ) {
        if (spin_adapt_g2_) {
            throw PsiException("If constraining T1/T2, G2 cannot currently be spin adapted.",__FILE__,__LINE__);
        }
        if (spin_adapt_q2_) {
            throw PsiException("If constraining T1/T2, Q2 cannot currently be spin adapted.",__FILE__,__LINE__);
        }
    }

    // build mapping arrays and determine the number of geminals per block
    BuildBasis();

    int ms = (multiplicity_ - 1)/2;
    if ( ms > 0 ) {
        if (spin_adapt_g2_) {
            throw PsiException("G2 not spin adapted for S = M != 0",__FILE__,__LINE__);
        }
        if (spin_adapt_q2_) {
            throw PsiException("Q2 not spin adapted for S = M != 0",__FILE__,__LINE__);
        }
    }

    // dimension of variable buffer (x) 
    dimx_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        dimx_ += gems_ab[h]*gems_ab[h]; // D2ab
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx_ += gems_aa[h]*gems_aa[h]; // D2aa
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx_ += gems_aa[h]*gems_aa[h]; // D2bb
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx_ += amopi_[h]*amopi_[h]; // D1a
        dimx_ += amopi_[h]*amopi_[h]; // D1b
        dimx_ += amopi_[h]*amopi_[h]; // Q1b
        dimx_ += amopi_[h]*amopi_[h]; // Q1a
    }
    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // Q2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2aa
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_00[h]*gems_00[h]; // Q2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_aa[h]*gems_aa[h]; // Q2t_m1
            }
        }
    }
    if ( constrain_g2_ ) {
        if ( !spin_adapt_g2_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2ba
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += 2*gems_ab[h]*2*gems_ab[h]; // G2aa/bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx_ += gems_ab[h]*gems_ab[h]; // G2t_m1
            }
        }
    }
    if ( constrain_t1_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aaa[h]*trip_aaa[h]; // T1aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aaa[h]*trip_aaa[h]; // T1bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // T1aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // T1bba
        }
    }
    if ( constrain_t2_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += (trip_aba[h]+trip_aab[h])*(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += (trip_aba[h]+trip_aab[h])*(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // T2aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // T2bba
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aba[h]*trip_aba[h]; // T2aba
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aba[h]*trip_aba[h]; // T2bab
        }
    }
    if ( constrain_d3_ ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aaa[h] * trip_aaa[h]; // D3aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aaa[h] * trip_aaa[h]; // D3bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // D3aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx_ += trip_aab[h]*trip_aab[h]; // D3bba
        }
    }

    // offsets in x
    offset = 0;

    d2aboff = (int*)malloc(nirrep_*sizeof(int));
    d2aaoff = (int*)malloc(nirrep_*sizeof(int));
    d2bboff = (int*)malloc(nirrep_*sizeof(int));
    for (int h = 0; h < nirrep_; h++) {
        d2aboff[h] = offset; offset += gems_ab[h]*gems_ab[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        d2aaoff[h] = offset; offset += gems_aa[h]*gems_aa[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        d2bboff[h] = offset; offset += gems_aa[h]*gems_aa[h];
    }

    d1aoff = (int*)malloc(nirrep_*sizeof(int));
    d1boff = (int*)malloc(nirrep_*sizeof(int));
    q1aoff = (int*)malloc(nirrep_*sizeof(int));
    q1boff = (int*)malloc(nirrep_*sizeof(int));
    for (int h = 0; h < nirrep_; h++) {
        d1aoff[h] = offset; offset += amopi_[h]*amopi_[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        d1boff[h] = offset; offset += amopi_[h]*amopi_[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        q1aoff[h] = offset; offset += amopi_[h]*amopi_[h];
    }
    for (int h = 0; h < nirrep_; h++) {
        q1boff[h] = offset; offset += amopi_[h]*amopi_[h];
    }

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            q2aboff = (int*)malloc(nirrep_*sizeof(int));
            q2aaoff = (int*)malloc(nirrep_*sizeof(int));
            q2bboff = (int*)malloc(nirrep_*sizeof(int));
            for (int h = 0; h < nirrep_; h++) {
                q2aboff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2aaoff[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2bboff[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
        }else {
            q2soff = (int*)malloc(nirrep_*sizeof(int));
            q2toff = (int*)malloc(nirrep_*sizeof(int));
            q2toff_p1 = (int*)malloc(nirrep_*sizeof(int));
            q2toff_m1 = (int*)malloc(nirrep_*sizeof(int));
            for (int h = 0; h < nirrep_; h++) {
                q2soff[h] = offset; offset += gems_00[h]*gems_00[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2toff[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2toff_p1[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                q2toff_m1[h] = offset; offset += gems_aa[h]*gems_aa[h];
            }
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            g2aboff = (int*)malloc(nirrep_*sizeof(int));
            g2baoff = (int*)malloc(nirrep_*sizeof(int));
            g2aaoff = (int*)malloc(nirrep_*sizeof(int));
            for (int h = 0; h < nirrep_; h++) {
                g2aboff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2baoff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2aaoff[h] = offset; offset += 2*gems_ab[h]*2*gems_ab[h];
            }
        }else {
            g2soff = (int*)malloc(nirrep_*sizeof(int));
            g2toff = (int*)malloc(nirrep_*sizeof(int));
            g2toff_p1 = (int*)malloc(nirrep_*sizeof(int));
            g2toff_m1 = (int*)malloc(nirrep_*sizeof(int));
            for (int h = 0; h < nirrep_; h++) {
                g2soff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2toff[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2toff_p1[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
            for (int h = 0; h < nirrep_; h++) {
                g2toff_m1[h] = offset; offset += gems_ab[h]*gems_ab[h];
            }
        }
    }

    if ( constrain_t1_ ) {
        t1aaboff = (int*)malloc(nirrep_*sizeof(int));
        t1bbaoff = (int*)malloc(nirrep_*sizeof(int));
        t1aaaoff = (int*)malloc(nirrep_*sizeof(int));
        t1bbboff = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            t1aaaoff[h] = offset; offset += trip_aaa[h]*trip_aaa[h]; // T1aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            t1bbboff[h] = offset; offset += trip_aaa[h]*trip_aaa[h]; // T1bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            t1aaboff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // T1aab
        }
        for (int h = 0; h < nirrep_; h++) {
            t1bbaoff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // T1bba
        }
    }

    if ( constrain_t2_ ) {
        t2aaboff = (int*)malloc(nirrep_*sizeof(int));
        t2bbaoff = (int*)malloc(nirrep_*sizeof(int));
        t2aaaoff = (int*)malloc(nirrep_*sizeof(int));
        t2bbboff = (int*)malloc(nirrep_*sizeof(int));
        t2abaoff = (int*)malloc(nirrep_*sizeof(int));
        t2baboff = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            t2aaaoff[h] = offset; offset += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            t2bbboff[h] = offset; offset += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            t2aaboff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // T2aab
        }
        for (int h = 0; h < nirrep_; h++) {
            t2bbaoff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // T2bba
        }
        for (int h = 0; h < nirrep_; h++) {
            t2abaoff[h] = offset; offset += trip_aba[h]*trip_aba[h]; // T2aba
        }
        for (int h = 0; h < nirrep_; h++) {
            t2baboff[h] = offset; offset += trip_aba[h]*trip_aba[h]; // T2bab
        }
    }
    if ( constrain_d3_ ) {
        d3aaaoff = (int*)malloc(nirrep_*sizeof(int));
        d3bbboff = (int*)malloc(nirrep_*sizeof(int));
        d3aaboff = (int*)malloc(nirrep_*sizeof(int));
        d3bbaoff = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            d3aaaoff[h] = offset; offset += trip_aaa[h]*trip_aaa[h]; // D3aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            d3bbboff[h] = offset; offset += trip_aaa[h]*trip_aaa[h]; // D3bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            d3aaboff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // D3aab
        }
        for (int h = 0; h < nirrep_; h++) {
            d3bbaoff[h] = offset; offset += trip_aab[h]*trip_aab[h]; // D3bba
        }
    }
    // constraints:
    nconstraints_ = 0;

    if ( constrain_spin_ ) {
        nconstraints_ += 1;               // spin
    }
    nconstraints_ += 1;                   // Tr(D2ab)
    nconstraints_ += 1;                   // Tr(D2aa)
    nconstraints_ += 1;                   // Tr(D2bb)
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints_ += amopi_[h]*amopi_[h]; // D1a <-> Q1a
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints_ += amopi_[h]*amopi_[h]; // D1b <-> Q1b
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints_ += amopi_[h]*amopi_[h]; // contract D2ab        -> D1 a
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints_ += amopi_[h]*amopi_[h]; // contract D2ab        -> D1 b
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints_ += amopi_[h]*amopi_[h]; // contract D2aa        -> D1 a
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints_ += amopi_[h]*amopi_[h]; // contract D2bb        -> D1 b
    }
    if ( constrain_q2_ ) {
        if ( ! spin_adapt_q2_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // Q2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2aa
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_00[h]*gems_00[h]; // Q2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_aa[h]*gems_aa[h]; // Q2t_m1
            }
        }
        
    }
    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2ba
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += 2*gems_ab[h]*2*gems_ab[h]; // G2aa
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints_ += gems_ab[h]*gems_ab[h]; // G2t_m1
            }
        }
    }
    if ( constrain_t1_ ) {
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aaa[h]*trip_aaa[h]; // T1aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aaa[h]*trip_aaa[h]; // T1bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aab[h]*trip_aab[h]; // T1aab
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aab[h]*trip_aab[h]; // T1bba
        }
    }
    if ( constrain_t2_ ) {
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aab[h]*trip_aab[h]; // T2aab
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aab[h]*trip_aab[h]; // T2bba
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aba[h]*trip_aba[h]; // T2aba
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += trip_aba[h]*trip_aba[h]; // T2bab
        }
    }
    if ( constrain_d3_ ) {
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += gems_aa[h]*gems_aa[h]; // D3aaa -> D2aa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += gems_aa[h]*gems_aa[h]; // D3bbb -> D2bb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += gems_aa[h]*gems_aa[h]; // D3aab -> D2aa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += gems_aa[h]*gems_aa[h]; // D3bba -> D2bb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += gems_ab[h]*gems_ab[h]; // D3aab -> D2ab
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints_ += gems_ab[h]*gems_ab[h]; // D3bba -> D2ab
        }
    }

    // list of dimensions_
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(gems_ab[h]); // D2ab
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(gems_aa[h]); // D2aa
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(gems_aa[h]); // D2bb
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(amopi_[h]); // D1a
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(amopi_[h]); // D1b
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(amopi_[h]); // Q1a
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions_.push_back(amopi_[h]); // Q1b
    }
    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // Q2ab
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2aa
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2bb
            }
        }else {
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_00[h]); // Q2s
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2t
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2t_p1
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_aa[h]); // Q2t_m1
            }
        }
    }
    if ( constrain_g2_ ) {
        if ( !spin_adapt_g2_ ) {
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2ab
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2ba
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(2*gems_ab[h]); // G2aa
            }
        }else {
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2s
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2t
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2t_p1
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions_.push_back(gems_ab[h]); // G2t_m1
            }
        }
    }
    if ( constrain_t1_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aaa[h]); // T1aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aaa[h]); // T1bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // T1aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // T1bba
        }
    }
    if ( constrain_t2_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // T2aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // T2bba
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aba[h]); // T2aba
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aba[h]); // T2bab
        }
    }
    if ( constrain_d3_ ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aaa[h]); // D3aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aaa[h]); // D3bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // D3aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions_.push_back(trip_aab[h]); // D3bba
        }
    }

    // v2rdm sdp convergence thresholds:
    r_convergence_  = options_.get_double("R_CONVERGENCE");
    e_convergence_  = options_.get_double("E_CONVERGENCE");
    maxiter_        = options_.get_int("MAXITER");

    // conjugate gradient solver thresholds:
    cg_convergence_ = options_.get_double("CG_CONVERGENCE");
    cg_maxiter_     = options_.get_double("CG_MAXITER");


    // memory check happens here

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    v-2RDM                                           *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Ground-state variational 2-RDM optimization      *\n");
    outfile->Printf( "        *    using a boundary-point semidefinite solver       *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Eugene DePrince                                  *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    outfile->Printf("\n");
    outfile->Printf("  ==> Input parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("        Freeze core orbitals?               %5s\n",nfrzc_ > 0 ? "yes" : "no");
    outfile->Printf("        Number of frozen core orbitals:     %5i\n",nfrzc_);
    outfile->Printf("        Number of active occupied orbitals: %5i\n",ndoccact);
    outfile->Printf("        Number of active virtual orbitals:  %5i\n",nvirt);
    outfile->Printf("        Number of frozen virtual orbitals:  %5i\n",nfrzv);
    outfile->Printf("        r_convergence:                  %5.3le\n",r_convergence_);
    outfile->Printf("        e_convergence:                  %5.3le\n",e_convergence_);
    outfile->Printf("        cg_convergence:                 %5.3le\n",cg_convergence_);
    outfile->Printf("        maxiter:                         %8i\n",maxiter_);
    outfile->Printf("        cg_maxiter:                      %8i\n",cg_maxiter_);
    outfile->Printf("\n");
    outfile->Printf("  ==> Memory requirements <==\n");
    outfile->Printf("\n");
    int nd2   = 0;
    int ng2    = 0;
    int nt1    = 0;
    int nt2    = 0;
    int maxgem = 0;
    for (int h = 0; h < nirrep_; h++) {
        nd2 +=     gems_ab[h]*gems_ab[h];
        nd2 += 2 * gems_aa[h]*gems_aa[h];

        ng2 +=     gems_ab[h] * gems_ab[h]; // G2ab
        ng2 +=     gems_ab[h] * gems_ab[h]; // G2ba
        ng2 += 4 * gems_ab[h] * gems_ab[h]; // G2aa

        if ( gems_ab[h] > maxgem ) {
            maxgem = gems_ab[h];
        }
        if ( constrain_g2_ ) {
            if ( 2*gems_ab[h] > maxgem ) {
                maxgem = 2*gems_ab[h];
            }
        }

        if ( constrain_t1_ ) {
            nt1 += trip_aaa[h] * trip_aaa[h]; // T1aaa
            nt1 += trip_aaa[h] * trip_aaa[h]; // T1bbb
            nt1 += trip_aab[h] * trip_aab[h]; // T1aab
            nt1 += trip_aab[h] * trip_aab[h]; // T1bba
            if ( trip_aab[h] > maxgem ) {
                maxgem = trip_aab[h];
            }
        }

        if ( constrain_t2_ ) {
            nt2 += (trip_aab[h]+trip_aba[h]) * (trip_aab[h]+trip_aba[h]); // T2aaa
            nt2 += (trip_aab[h]+trip_aba[h]) * (trip_aab[h]+trip_aba[h]); // T2bbb
            nt2 += trip_aab[h] * trip_aab[h]; // T2aab
            nt2 += trip_aab[h] * trip_aab[h]; // T2bba
            nt2 += trip_aba[h] * trip_aba[h]; // T2aba
            nt2 += trip_aba[h] * trip_aba[h]; // T2bab

            if ( trip_aab[h]+trip_aaa[h] > maxgem ) {
                maxgem = trip_aab[h]+trip_aaa[h];
            }
        }

    }

    outfile->Printf("        D2:                       %7.2lf mb\n",nd2 * 8.0 / 1024.0 / 1024.0);
    if ( constrain_q2_ ) {
        outfile->Printf("        Q2:                       %7.2lf mb\n",nd2 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_g2_ ) {
        outfile->Printf("        G2:                       %7.2lf mb\n",ng2 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_d3_ ) {
        outfile->Printf("        D3:                       %7.2lf mb\n",nt1 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_t1_ ) {
        outfile->Printf("        T1:                       %7.2lf mb\n",nt1 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_t2_ ) {
        outfile->Printf("        T2:                       %7.2lf mb\n",nt2 * 8.0 / 1024.0 / 1024.0);
    }
    outfile->Printf("\n");

    // we have 4 arrays the size of x and 4 the size of y
    // in addition, we need to store whatever the 2x largest block is
    // for the diagonalization step
    // integrals:
    //     K2a, K2b
    // casscf:
    //     4-index integrals (no permutational symmetry)
    //     3-index integrals 

    double tot = 4.0*dimx_ + 4.0*nconstraints_ + 2.0*maxgem*maxgem;
    tot += nd2; // for K2a, K2b

    // for casscf, need d2 and 3- or 4-index integrals

    // allocate memory for full ERI tensor blocked by symmetry for Greg
    tei_full_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        tei_full_dim_ += gems_full[h] * ( gems_full[h] + 1 ) / 2;
    }
    d2_plus_core_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim_ += gems_plus_core[h] * ( gems_plus_core[h] + 1 ) / 2;
    }
    tot += d2_plus_core_dim_;
    if ( is_df_ ) {
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
        tot += (long int)nQ_*(long int)nmo_*((long int)nmo_+1)/2;
    }else {
        tei_full_dim_ = 0;
        for (int h = 0; h < nirrep_; h++) {
            tei_full_dim_ += gems_full[h] * ( gems_full[h] + 1 ) / 2;
        }
        tot += tei_full_dim_;
        tot += nmo_*nmo_*nmo_*nmo_; // for four-index integrals stored stupidly 
    }
    
    outfile->Printf("        Total number of variables:     %10i\n",dimx_);
    outfile->Printf("        Total number of constraints:   %10i\n",nconstraints_);
    outfile->Printf("        Total memory requirements:     %7.2lf mb\n",tot * 8.0 / 1024.0 / 1024.0);
    outfile->Printf("\n");

    if ( tot * 8.0 > (double)memory_ ) {
        outfile->Printf("\n");
        outfile->Printf("        Not enough memory!\n");
        outfile->Printf("\n");
        if ( !is_df_ ) {
            outfile->Printf("        Either increase the available memory by %7.2lf mb\n",(8.0 * tot - memory_)/1024.0/1024.0);
            outfile->Printf("        or try scf_type = df or scf_type = cd\n");
        
        }else {
            outfile->Printf("        Increase the available memory by %7.2lf mb.\n",(8.0 * tot - memory_)/1024.0/1024.0);
        }
        outfile->Printf("\n");
        throw PsiException("Not enough memory",__FILE__,__LINE__);
    }

    // if using 3-index integrals, transform them before allocating any memory integrals, transform 
    if ( is_df_ ) {
        outfile->Printf("    ==> Transform three-electron integrals <==\n");
        outfile->Printf("\n");

        double start = omp_get_wtime();
        ThreeIndexIntegrals();
        double end = omp_get_wtime();

        outfile->Printf("\n");
        outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
        outfile->Printf("\n");
    } else {
        // transform integrals
        outfile->Printf("    ==> Transform two-electron integrals <==\n");
        outfile->Printf("\n");
        
        double start = omp_get_wtime();
        boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
        std::vector<shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::all);
        boost::shared_ptr<IntegralTransform> ints(new IntegralTransform(wfn, spaces, IntegralTransform::Restricted,
            				      IntegralTransform::IWLOnly, IntegralTransform::PitzerOrder, IntegralTransform::None, false));
        ints->set_dpd_id(0);
        ints->set_keep_iwl_so_ints(true);
        ints->set_keep_dpd_so_ints(true);
        ints->initialize();
        ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
        double end = omp_get_wtime();
        outfile->Printf("\n");
        outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
        outfile->Printf("\n");
    
    }

    // allocate vectors
    Ax     = SharedVector(new Vector("A . x",nconstraints_));
    ATy    = SharedVector(new Vector("A^T . y",dimx_));
    x      = SharedVector(new Vector("primal solution",dimx_));
    c      = SharedVector(new Vector("OEI and TEI",dimx_));
    y      = SharedVector(new Vector("dual solution",nconstraints_));
    z      = SharedVector(new Vector("dual solution 2",dimx_));
    b      = SharedVector(new Vector("constraints",nconstraints_));

    // input/output array for orbopt sweeps

    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    orbopt_data_    = (double*)malloc(13*sizeof(double));
    orbopt_data_[0] = (double)nthread;
    orbopt_data_[1] = (double)options_.get_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS");
    orbopt_data_[2] = (double)options_.get_int("ORBOPT_FROZEN_CORE");
    orbopt_data_[3] = (double)options_.get_double("ORBOPT_GRADIENT_TOLERANCE");
    orbopt_data_[4] = (double)options_.get_double("ORBOPT_ENERGY_CONVERGENCE");
    orbopt_data_[5] = (double)options_.get_bool("ORBOPT_WRITE");
    orbopt_data_[6] = (double)options_.get_int("ORBOPT_EXACT_DIAGONAL_HESSIAN");
    orbopt_data_[7] = (double)options_.get_int("ORBOPT_NUM_DIIS_VECTORS");
    orbopt_data_[8] = 0.0;
    if ( is_df_ ) {
      orbopt_data_[8] = 1.0;
    }
    orbopt_data_[9] = 0.0;  // number of iterations (output)
    orbopt_data_[10] = 0.0;  // gradient norm (output)
    orbopt_data_[11] = 0.0;  // change in energy (output)
    orbopt_data_[12] = 0.0;  // converged?
    orbopt_converged_ = false;

    orbopt_transformation_matrix_ = (double*)malloc(nmo_*nmo_*sizeof(double));
    memset((void*)orbopt_transformation_matrix_,'\0',nmo_*nmo_*sizeof(double));
    for (int i = 0; i < nmo_; i++) {
        orbopt_transformation_matrix_[i*nmo_+i] = 1.0;
    }

    // don't change the length of this filename
    orbopt_outfile_ = (char*)malloc(120*sizeof(char));
    std::string filename = get_writer_file_prefix() + ".orbopt";
    strcpy(orbopt_outfile_,filename.c_str());
    if ( options_.get_bool("ORBOPT_WRITE") ) { 
        FILE * fp = fopen(orbopt_outfile_,"w");
        fclose(fp);
    }

}

int v2RDMSolver::SymmetryPair(int i,int j) {
    return table[i*8+j];
}
int v2RDMSolver::TotalSym(int i,int j,int k, int l) {
    return SymmetryPair(SymmetryPair(symmetry[i],symmetry[j]),SymmetryPair(symmetry[k],symmetry[l]));
}

// compute the energy!
double v2RDMSolver::compute_energy() {

    // hartree-fock guess
    Guess();

    // get integrals
    if ( is_df_ ) {
        DF_TEI();
    }else {
        TEI();
    }

    // generate constraint vector
    BuildConstraints();

    // AATy = A(c-z)+tu(b-Ax) rearange w.r.t cg solver
    // Ax   = AATy and b=A(c-z)+tu(b-Ax)
    SharedVector B   = SharedVector(new Vector("compound B",nconstraints_));

    tau = 1.6;
    mu  = 1.0;

    // congugate gradient solver
    long int N = nconstraints_;
    shared_ptr<CGSolver> cg (new CGSolver(N));
    cg->set_max_iter(cg_maxiter_);
    cg->set_convergence(cg_convergence_);

    // evaluate guess energy (c.x):
    double energy_primal = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);

    outfile->Printf("\n");
    outfile->Printf("    reference energy:     %20.12lf\n",escf_);
    outfile->Printf("    frozen core energy:   %20.12lf\n",efzc_);
    outfile->Printf("    initial 2-RDM energy: %20.12lf\n",energy_primal + enuc_ + efzc_);
    outfile->Printf("\n");
    outfile->Printf("      oiter");
    outfile->Printf(" iiter");
    outfile->Printf("        E(p)");
    outfile->Printf("        E(d)");
    outfile->Printf("      E gap)");
    outfile->Printf("      mu");
    outfile->Printf("     eps(p)");
    outfile->Printf("     eps(d)\n");

    double energy_dual,egap;
    double denergy_primal = fabs(energy_primal);

    int oiter=0;
    do {

        double start = omp_get_wtime();

        // evaluate tau * mu * (b - Ax) for CG
        bpsdp_Au(Ax, x);
        Ax->subtract(b);
        Ax->scale(-tau*mu);
        
        // evaluate A(c-z) ( but don't overwrite c! )
        z->scale(-1.0);
        z->add(c);
        bpsdp_Au(B,z);
        
        // add tau*mu*(b-Ax) to A(c-z) and put result in B
        B->add(Ax);
        // solve CG problem (step 1 in table 1 of PRL 106 083001)
        if (oiter == 0) cg->set_convergence(0.01);
        else            cg->set_convergence( ( ep > ed ) ? 0.01 * ed : 0.01 * ep);
        cg->solve(N,Ax,y,B,evaluate_Ap,(void*)this);
        int iiter = cg->total_iterations();

        double end = omp_get_wtime();
        double cg_time = end - start;

        start = omp_get_wtime();
        // build and diagonalize U according to step 2 in PRL
        // and update z and z
        Update_xz();
        end = omp_get_wtime();
        double diag_time = end - start;

        // update mu (step 3)

        // evaluate || A^T y - c + z||
        bpsdp_ATu(ATy, y);
        ATy->add(z);
        ATy->subtract(c);
        ed = sqrt(ATy->norm());
        
        // evaluate || Ax - b ||
        bpsdp_Au(Ax, x);
        Ax->subtract(b);
        ep = sqrt(Ax->norm());

        //RotateOrbitals();
        // don't update mu every iteration
        if ( oiter % 200 == 0 && oiter > 0) {
            mu = mu*ep/ed;
        }
        if ( oiter % 200 == 0 && oiter > 0) {

            // gidofalvi debug
            //outfile->Printf("current nuclear repulsion energy is %20.13lf \n",enuc_);
            //outfile->Printf("             current core energy is %20.13lf \n",efzc_);
            //double current_energy_tmp = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);
            //outfile->Printf("           current active energy is %20.13lf \n",current_energy_tmp);
            //outfile->Printf("            current total energy is %20.13lf \n",current_energy_tmp+efzc_+enuc_);
            // end debug

            RotateOrbitals();
        }

        // compute current primal and dual energies

        //energy_primal = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);
        double current_energy = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);

        energy_dual   = C_DDOT(nconstraints_,b->pointer(),1,y->pointer(),1);

        outfile->Printf("      %5i %5i %11.6lf %11.6lf %11.6lf %7.3lf %10.5lf %10.5lf\n",
                    oiter,iiter,current_energy+enuc_+efzc_,energy_dual+efzc_+enuc_,fabs(current_energy-energy_dual),mu,ep,ed);
        oiter++;
    
        if (oiter == maxiter_) break;

        egap = fabs(current_energy-energy_dual);
        denergy_primal = fabs(energy_primal - current_energy);
        energy_primal = current_energy;

        if ( ep < r_convergence_ && ed < r_convergence_ && egap < e_convergence_ ) {
            RotateOrbitals();
            energy_primal = C_DDOT(dimx_,c->pointer(),1,x->pointer(),1);
        }

    }while( ep > r_convergence_ || ed > r_convergence_  || egap > e_convergence_ || !orbopt_converged_);

    if ( oiter == maxiter_ ) {
        throw PsiException("v2RDM did not converge.",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("      v2RDM iterations converged!\n");
    outfile->Printf("\n");

    // evaluate spin squared
    double s2 = 0.0;
    double * x_p = x->pointer();
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            int ij = ibas_ab_sym[h][i][j];
            int ji = ibas_ab_sym[h][j][i];
            s2 += x_p[d2aboff[h] + ij*gems_ab[h]+ji];
        }
    }
    int na = nalpha_ - nfrzc_;
    int nb = nbeta_ - nfrzc_;
    int ms = (multiplicity_ - 1)/2;
    outfile->Printf("      v2RDM total spin [S(S+1)]: %20.6lf\n", 0.5 * (na + nb) + ms*ms - s2);
    outfile->Printf("    * v2RDM total energy:        %20.12lf\n",energy_primal+enuc_+efzc_);
    outfile->Printf("\n");

    Process::environment.globals["CURRENT ENERGY"]     = energy_primal+enuc_+efzc_;
    Process::environment.globals["v2RDM TOTAL ENERGY"] = energy_primal+enuc_+efzc_;

    // compute and print natural orbital occupation numbers
    FinalTransformationMatrix();
    MullikenPopulations();
    NaturalOrbitals();

    return energy_primal + enuc_ + efzc_;
}

void v2RDMSolver::NaturalOrbitals() {
    boost::shared_ptr<Matrix> Da (new Matrix(nirrep_,nmopi_,nmopi_));
    boost::shared_ptr<Matrix> eigveca (new Matrix(nirrep_,nmopi_,nmopi_));
    boost::shared_ptr<Vector> eigvala (new Vector("Natural Orbital Occupation Numbers (alpha)",nirrep_,nmopi_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            Da->pointer(h)[i][i] = 1.0;
        }
        for (int i = frzcpi_[h]; i < nmopi_[h]-frzvpi_[h]; i++) {
            for (int j = frzcpi_[h]; j < nmopi_[h]-frzvpi_[h]; j++) {
                Da->pointer(h)[i][j] = x->pointer()[d1aoff[h]+(i-frzcpi_[h])*amopi_[h]+(j-frzcpi_[h])];
            }
        }
    }
    boost::shared_ptr<Matrix> saveda ( new Matrix(Da) );
    Da->diagonalize(eigveca,eigvala,descending);
    eigvala->print();
    //Ca_->print();
    // build AO/NO transformation matrix (Ca_)
    for (int h = 0; h < nirrep_; h++) {
        for (int mu = 0; mu < nsopi_[h]; mu++) {
            double *  temp = (double*)malloc(nmopi_[h]*sizeof(double));
            double ** cp   = Ca_->pointer(h);
            double ** ep   = eigveca->pointer(h);
            for (int i = 0; i < nmopi_[h]; i++) {
                double dum = 0.0;
                for (int j = 0; j < nmopi_[h]; j++) {
                    dum += cp[mu][j] * ep[j][i];
                }
                temp[i] = dum;
            }
            for (int i = 0; i < nmopi_[h]; i++) {
                cp[mu][i] = temp[i];
            }
            free(temp);
        }
    }
    //Ca_->print();

    boost::shared_ptr<Matrix> Db (new Matrix(nirrep_,nmopi_,nmopi_));
    boost::shared_ptr<Matrix> eigvecb (new Matrix(nirrep_,nmopi_,nmopi_));
    boost::shared_ptr<Vector> eigvalb (new Vector("Natural Orbital Occupation Numbers (beta)",nirrep_,nmopi_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            Db->pointer(h)[i][i] = 1.0;
        }
        for (int i = frzcpi_[h]; i < nmopi_[h]-frzvpi_[h]; i++) {
            for (int j = frzcpi_[h]; j < nmopi_[h]-frzvpi_[h]; j++) {
                Db->pointer(h)[i][j] = x->pointer()[d1boff[h]+(i-frzcpi_[h])*amopi_[h]+(j-frzcpi_[h])];
            }
        }
    }
    Db->diagonalize(eigvecb,eigvalb,descending);
    eigvalb->print();
    // build AO/NO transformation matrix (Cb_)
    for (int h = 0; h < nirrep_; h++) {
        for (int mu = 0; mu < nsopi_[h]; mu++) {
            double * temp  = (double*)malloc(nmopi_[h]*sizeof(double));
            double ** cp   = Cb_->pointer(h);
            double ** ep   = eigvecb->pointer(h);
            for (int i = 0; i < nmopi_[h]; i++) {
                double dum = 0.0;
                for (int j = 0; j < nmopi_[h]; j++) {
                    dum += cp[mu][j] * ep[j][i];
                }
                temp[i] = dum;
            }
            for (int i = 0; i < nmopi_[h]; i++) {
                cp[mu][i] = temp[i];
            }
            free(temp);
        }
    }

    // Print a molden file
    if ( options_.get_bool("MOLDEN_WRITE") ) {
       //boost::shared_ptr<MoldenWriter> molden(new MoldenWriter((boost::shared_ptr<Wavefunction>)this));
       boost::shared_ptr<MoldenWriter> molden(new MoldenWriter(reference_wavefunction_));
       boost::shared_ptr<Vector> zero (new Vector("",nirrep_,nmopi_));
       zero->zero();
       std::string filename = get_writer_file_prefix() + ".molden";
       molden->write(filename,Ca_,Cb_,zero, zero,eigvala,eigvalb);
    }
}

void v2RDMSolver::MullikenPopulations() {
    // nee

    std::stringstream ss;
    ss << "v-2RDM";
    std::stringstream ss_a;
    std::stringstream ss_b;
    ss_a << ss.str() << " alpha";
    ss_b << ss.str() << " beta";
    boost::shared_ptr<Matrix> opdm_a(new Matrix(ss_a.str(), Ca_->colspi(), Ca_->colspi()));
    boost::shared_ptr<Matrix> opdm_b(new Matrix(ss_b.str(), Ca_->colspi(), Ca_->colspi()));

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            opdm_a->pointer(h)[i][i] = 1.0;
        }
        for (int i = frzcpi_[h]; i < nmopi_[h]-frzvpi_[h]; i++) {
            for (int j = frzcpi_[h]; j < nmopi_[h]-frzvpi_[h]; j++) {
                opdm_a->pointer(h)[i][j] = x->pointer()[d1aoff[h]+(i-frzcpi_[h])*amopi_[h]+(j-frzcpi_[h])];
            }
        }
    }

    int symm = opdm_a->symmetry();

    double* temp = (double*)malloc(Ca_->max_ncol() * Ca_->max_nrow() * sizeof(double));

    Da_->zero();
    for (int h = 0; h < nirrep_; h++) {
        int nmol = Ca_->colspi()[h];
        int nmor = Ca_->colspi()[h^symm];
        int nsol = Ca_->rowspi()[h];
        int nsor = Ca_->rowspi()[h^symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Clp = Ca_->pointer(h);
        double** Crp = Ca_->pointer(h^symm);
        double** Dmop = opdm_a->pointer(h^symm);
        double** Dsop = Da_->pointer(h^symm);
        C_DGEMM('N','T',nmol,nsor,nmor,1.0,Dmop[0],nmor,Crp[0],nmor,0.0,temp,nsor);
        C_DGEMM('N','N',nsol,nsor,nmol,1.0,Clp[0],nmol,temp,nsor,0.0,Dsop[0],nsor);
    }

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            opdm_b->pointer(h)[i][i] = 1.0;
        }
        for (int i = frzcpi_[h]; i < nmopi_[h]-frzvpi_[h]; i++) {
            for (int j = frzcpi_[h]; j < nmopi_[h]-frzvpi_[h]; j++) {
                opdm_b->pointer(h)[i][j] = x->pointer()[d1boff[h]+(i-frzcpi_[h])*amopi_[h]+(j-frzcpi_[h])];
            }
        }
    }
    Db_->zero();
    // hmm... the IntegralTransform type is Restricted, so only use Ca_ here.
    for (int h = 0; h < nirrep_; h++) {
        int nmol = Ca_->colspi()[h];
        int nmor = Ca_->colspi()[h^symm];
        int nsol = Ca_->rowspi()[h];
        int nsor = Ca_->rowspi()[h^symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Clp = Ca_->pointer(h);
        double** Crp = Ca_->pointer(h^symm);
        double** Dmop = opdm_b->pointer(h^symm);
        double** Dsop = Db_->pointer(h^symm);
        C_DGEMM('N','T',nmol,nsor,nmor,1.0,Dmop[0],nmor,Crp[0],nmor,0.0,temp,nsor);
        C_DGEMM('N','N',nsol,nsor,nmol,1.0,Clp[0],nmol,temp,nsor,0.0,Dsop[0],nsor);
    }

    // compute mulliken charges:
    boost::shared_ptr<MyOEProp> oe(new MyOEProp());
    oe->compute_mulliken_charges_custom(Da_,Db_);
    
    free(temp);    
}

void v2RDMSolver::Guess(){

    double* x_p = x->pointer();
    double* z_p = z->pointer();

    memset((void*)x_p,'\0',dimx_*sizeof(double));
    memset((void*)z_p,'\0',dimx_*sizeof(double));

    // random guess
    srand(0);
    for (int i = 0; i < dimx_; i++) {
        x_p[i] = ( (double)rand()/RAND_MAX - 1.0 ) * 2.0;
    }
    return;

    // Hartree-Fock guess for D2, D1, Q1, Q2, and G2

    // D2ab
    int poff1 = 0;
    for (int h1 = 0; h1 < nirrep_; h1++) {
        for (int i = 0; i < soccpi_[h1] + doccpi_[h1] - frzcpi_[h1]; i++){
            int poff2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (int j = 0; j < doccpi_[h2] - frzcpi_[h2]; j++){
                    int ii = i + poff1;
                    int jj = j + poff2;
                    int h3 = SymmetryPair(symmetry[ii],symmetry[jj]);
                    int ij = ibas_ab_sym[h3][ii][jj];
                    x_p[d2aboff[h3] + ij*gems_ab[h3]+ij] = 1.0;
                }
                poff2   += nmopi_[h2] - frzcpi_[h2] - frzvpi_[h2];
            }
        }
        poff1   += nmopi_[h1] - frzcpi_[h1] - frzvpi_[h1];
    }

    // d2aa
    poff1 = 0;
    for (int h1 = 0; h1 < nirrep_; h1++) {
        for (int i = 0; i < soccpi_[h1] + doccpi_[h1] - frzcpi_[h1]; i++){
            int poff2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (int j = 0; j < soccpi_[h2] + doccpi_[h2] - frzcpi_[h2]; j++){
                    int ii = i + poff1;
                    int jj = j + poff2;
                    if ( jj >= ii ) continue;
                    int h3 = SymmetryPair(symmetry[ii],symmetry[jj]);
                    int ij = ibas_aa_sym[h3][ii][jj];
                    x_p[d2aaoff[h3] + ij*gems_aa[h3]+ij] = 1.0;
                }
                poff2   += nmopi_[h2] - frzcpi_[h2] - frzvpi_[h2];
            }
        }
        poff1   += nmopi_[h1] - frzcpi_[h1] - frzvpi_[h1];
    }

    // d2bb
    poff1 = 0;
    for (int h1 = 0; h1 < nirrep_; h1++) {
        for (int i = 0; i < doccpi_[h1] - frzcpi_[h1]; i++){
            int poff2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (int j = 0; j < doccpi_[h2] - frzcpi_[h2]; j++){
                    int ii = i + poff1;
                    int jj = j + poff2;
                    if ( jj >= ii ) continue;
                    int h3 = SymmetryPair(symmetry[ii],symmetry[jj]);
                    int ij = ibas_aa_sym[h3][ii][jj];
                    x_p[d2bboff[h3] + ij*gems_aa[h3]+ij] = 1.0;
                }
                poff2   += nmopi_[h2] - frzcpi_[h2] - frzvpi_[h2];
            }
        }
        poff1   += nmopi_[h1] - frzcpi_[h1] - frzvpi_[h1];
    }

    // D1
    for (int h = 0; h < nirrep_; h++) {
        for (int i = frzcpi_[h]; i < doccpi_[h]+soccpi_[h]; i++) {
            int ii = i - frzcpi_[h];
            x_p[d1aoff[h]+ii*amopi_[h]+ii] = 1.0;
        }
        for (int i = frzcpi_[h]; i < doccpi_[h]; i++) {
            int ii = i - frzcpi_[h];
            x_p[d1boff[h]+ii*amopi_[h]+ii] = 1.0;
        }
        // Q1
        for (int i = doccpi_[h]+soccpi_[h]; i < nmopi_[h]-frzvpi_[h]; i++) {
            int ii = i - frzcpi_[h];
            x_p[q1aoff[h]+ii*amopi_[h]+ii] = 1.0;
        }
        for (int i = doccpi_[h]; i < nmopi_[h]-frzvpi_[h]; i++) {
            int ii = i - frzcpi_[h];
            x_p[q1boff[h]+ii*amopi_[h]+ii] = 1.0;
        }
    }

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_) {
            Q2_constraints_guess(x);
        }else {
            Q2_constraints_guess_spin_adapted(x);
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            G2_constraints_guess(x);
        }else {
            G2_constraints_guess_spin_adapted(x);
        }
    }

    if ( constrain_t1_ ) {
        T1_constraints_guess(x);
    }
    if ( constrain_t2_ ) {
        T2_constraints_guess(x);
    }

}

void v2RDMSolver::BuildConstraints(){

    //constraint on the Trace of D2(s=0,ms=0)

    int na = nalpha_ - nfrzc_;
    int nb = nbeta_ - nfrzc_;
    double trdab = na * nb;

    //constraint on the Trace of D2(s=1,ms=0)
    double trdaa  = na*(na-1.0);
    double trdbb  = nb*(nb-1.0);

    b->zero();
    double* b_p = b->pointer();

    offset = 0;

    // funny ab trace with spin: N/2 + Ms^2 - S(S+1)
    if ( constrain_spin_ ) {
        double ms = (multiplicity_-1.0)/2.0;
        b_p[offset++] = (0.5 * (na + nb) + ms*ms - ms*(ms+1.0));
    }

    ///Trace of D2(s=0,ms=0) and D2(s=1,ms=0)
    b_p[offset++] = trdab;   
    b_p[offset++] = trdaa;
    b_p[offset++] = trdbb;


    // d1 / q1 a
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                b_p[offset + i*amopi_[h]+j] = (double)(i==j);
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    // d1 / q1 b
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                b_p[offset + i*amopi_[h]+j] = (double)(i==j);
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    //contract D2ab -> D1a
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                b_p[offset + i*amopi_[h]+j] = 0.0;
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    //contract D2ab -> D1b
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                b_p[offset + i*amopi_[h]+j] = 0.0;
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    //contract D2aa -> D1a
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                b_p[offset + i*amopi_[h]+j] = 0.0;
            }
        }
        offset += amopi_[h]*amopi_[h];
    }
    //contract D2bb -> D1b
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                b_p[offset + i*amopi_[h]+j] = 0.0;
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            // map d2ab to q2ab
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h]; 
            }

            // map d2aa to q2aa
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }

            // map d2bb to q2bb
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h]; 
            }
        }else {
            // map d2ab to q2s
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_00[h]; i++){
                    for(int j = 0; j < gems_00[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_00[h]*gems_00[h];
            }
            // map d2ab to q2t
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
            // map d2aa to q2t_p1
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
            // map d2bb to q2t_m1
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_aa[h]; i++){
                    for(int j = 0; j < gems_aa[h]; j++){
                        b_p[offset + INDEX(i,j)] = 0.0;
                    }
                }
                offset += gems_aa[h]*gems_aa[h];
            }
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            // map d2 and d1 to g2ab
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }

            // map d2 and d1 to g2ba
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }

            // map d2 and d1 to g2aa
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < 2*gems_ab[h]; i++){
                    for(int j = 0; j < 2*gems_ab[h]; j++){
                        b_p[offset + i*2*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += 2*gems_ab[h]*2*gems_ab[h]; 
            }
        }else {
            // map d2 and d1 to g2s
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h]; 
            }
            // map d2 and d1 to g2t
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h]; 
            }
            // map d2 and d1 to g2t_p1
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h]; 
            }
            // map d2 and d1 to g2t_m1
            for (int h = 0; h < nirrep_; h++) {
                for(int i = 0; i < gems_ab[h]; i++){
                    for(int j = 0; j < gems_ab[h]; j++){
                        b_p[offset + i*gems_ab[h]+j] = 0.0;
                    }
                }
                offset += gems_ab[h]*gems_ab[h]; 
            }
        }
    }

    if ( constrain_t1_ ) {
        // T1aaa
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aaa[h]; i++){
                for(int j = 0; j < trip_aaa[h]; j++){
                    b_p[offset + i*trip_aaa[h]+j] = 0.0;
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
        // T1bbb
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aaa[h]; i++){
                for(int j = 0; j < trip_aaa[h]; j++){
                    b_p[offset + i*trip_aaa[h]+j] = 0.0;
                }
            }
            offset += trip_aaa[h]*trip_aaa[h];
        }
        // T1aab
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aab[h]; i++){
                for(int j = 0; j < trip_aab[h]; j++){
                    b_p[offset + i*trip_aab[h]+j] = 0.0;
                }
            }
            offset += trip_aab[h]*trip_aab[h];
        }
        // T1bba
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aab[h]; i++){
                for(int j = 0; j < trip_aab[h]; j++){
                    b_p[offset + i*trip_aab[h]+j] = 0.0;
                }
            }
            offset += trip_aab[h]*trip_aab[h];
        }
    }

    if ( constrain_t2_ ) {
        // T2aaa
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aba[h]+trip_aab[h]; i++){
                for(int j = 0; j < trip_aba[h] + trip_aab[h]; j++){
                    b_p[offset + i*(trip_aab[h]+trip_aba[h])+j] = 0.0;
                }
            }
            offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
            //for(int i = 0; i < trip_aab[h]; i++){
            //    for(int j = 0; j < trip_aab[h]; j++){
            //        b_p[offset + i*trip_aab[h]+j] = 0.0;
            //    }
            //}
            //offset += trip_aab[h]*trip_aab[h];
        }
        // T2bbb
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aba[h]+trip_aab[h]; i++){
                for(int j = 0; j < trip_aba[h] + trip_aab[h]; j++){
                    b_p[offset + i*(trip_aab[h]+trip_aba[h])+j] = 0.0;
                }
            }
            offset += (trip_aba[h]+trip_aab[h])*(trip_aba[h]+trip_aab[h]);
            //for(int i = 0; i < trip_aaa[h]; i++){
            //    for(int j = 0; j < trip_aaa[h]; j++){
            //        b_p[offset + i*trip_aaa[h]+j] = 0.0;
            //    }
            //}
            //offset += trip_aaa[h]*trip_aaa[h];
        }
        // T2aab
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aab[h]; i++){
                for(int j = 0; j < trip_aab[h]; j++){
                    b_p[offset + i*trip_aab[h]+j] = 0.0;
                }
            }
            offset += trip_aab[h]*trip_aab[h];
        }
        // T2bba
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aab[h]; i++){
                for(int j = 0; j < trip_aab[h]; j++){
                    b_p[offset + i*trip_aab[h]+j] = 0.0;
                }
            }
            offset += trip_aab[h]*trip_aab[h];
        }
        // T2aba
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aba[h]; i++){
                for(int j = 0; j < trip_aba[h]; j++){
                    b_p[offset + i*trip_aba[h]+j] = 0.0;
                }
            }
            offset += trip_aba[h]*trip_aba[h];
        }
        // T2bab
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < trip_aba[h]; i++){
                for(int j = 0; j < trip_aba[h]; j++){
                    b_p[offset + i*trip_aba[h]+j] = 0.0;
                }
            }
            offset += trip_aba[h]*trip_aba[h];
        }
    }
    if ( constrain_d3_ ) {
        // D3aaa -> D2aa
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < gems_aa[h]; i++){
                for(int j = 0; j < gems_aa[h]; j++){
                    b_p[offset + i*gems_aa[h]+j] = 0.0;
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D3bbb -> D2bb
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < gems_aa[h]; i++){
                for(int j = 0; j < gems_aa[h]; j++){
                    b_p[offset + i*gems_aa[h]+j] = 0.0;
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D3aab -> D2aa
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < gems_aa[h]; i++){
                for(int j = 0; j < gems_aa[h]; j++){
                    b_p[offset + i*gems_aa[h]+j] = 0.0;
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D3bba -> D2bb
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < gems_aa[h]; i++){
                for(int j = 0; j < gems_aa[h]; j++){
                    b_p[offset + i*gems_aa[h]+j] = 0.0;
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D3aab -> D2ab
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < gems_ab[h]; i++){
                for(int j = 0; j < gems_ab[h]; j++){
                    b_p[offset + i*gems_ab[h]+j] = 0.0;
                }
            }
            offset += gems_ab[h]*gems_ab[h];
        }
        // D3bba -> D2ab
        for (int h = 0; h < nirrep_; h++) {
            for(int i = 0; i < gems_ab[h]; i++){
                for(int j = 0; j < gems_ab[h]; j++){
                    b_p[offset + i*gems_ab[h]+j] = 0.0;
                }
            }
            offset += gems_ab[h]*gems_ab[h];
        }
    }

}

///Build A dot u where u =[z,c]
void v2RDMSolver::bpsdp_Au(SharedVector A, SharedVector u){

    //A->zero();  
    memset((void*)A->pointer(),'\0',nconstraints_*sizeof(double));

    offset = 0;
    D2_constraints_Au(A,u);

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            Q2_constraints_Au(A,u);
        }else {
            Q2_constraints_Au_spin_adapted(A,u);
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            G2_constraints_Au(A,u);
        }else {
            G2_constraints_Au_spin_adapted(A,u);
        }
    }

    if ( constrain_t1_ ) {
        T1_constraints_Au(A,u);
    }

    if ( constrain_t2_ ) {
        //T2_constraints_Au(A,u);
        T2_constraints_Au_slow(A,u);
    }

    if ( constrain_d3_ ) {
        D3_constraints_Au(A,u);
    }

} // end Au

void v2RDMSolver::bpsdp_Au_slow(SharedVector A, SharedVector u){

    //A->zero();  
    memset((void*)A->pointer(),'\0',nconstraints_*sizeof(double));

    offset = 0;
    D2_constraints_Au(A,u);

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            Q2_constraints_Au(A,u);
        }else {
            Q2_constraints_Au_spin_adapted(A,u);
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            G2_constraints_Au(A,u);
        }else {
            G2_constraints_Au_spin_adapted(A,u);
        }
    }

    if ( constrain_t1_ ) {
        T1_constraints_Au(A,u);
    }

    if ( constrain_t2_ ) {
        //T2_constraints_Au(A,u);
        T2_constraints_Au_slow(A,u);
    }
    if ( constrain_d3_ ) {
        D3_constraints_Au(A,u);
    }

} // end Au

///Build AT dot u where u =[z,c]
void v2RDMSolver::bpsdp_ATu(SharedVector A, SharedVector u){

    //A->zero();
    memset((void*)A->pointer(),'\0',dimx_*sizeof(double));

    offset = 0;
    D2_constraints_ATu(A,u);

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            Q2_constraints_ATu(A,u);
        }else {
            Q2_constraints_ATu_spin_adapted(A,u);
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            G2_constraints_ATu(A,u);
        }else {
            G2_constraints_ATu_spin_adapted(A,u);
        }
    }

    if ( constrain_t1_ ) {
        T1_constraints_ATu(A,u);
    }

    if ( constrain_t2_ ) {
        //T2_constraints_ATu(A,u);
        T2_constraints_ATu_slow(A,u);
    }
    if ( constrain_d3_ ) {
        D3_constraints_ATu(A,u);
    }

}//end ATu

void v2RDMSolver::bpsdp_ATu_slow(SharedVector A, SharedVector u){

    //A->zero();
    memset((void*)A->pointer(),'\0',dimx_*sizeof(double));

    offset = 0;
    D2_constraints_ATu(A,u);

    if ( constrain_q2_ ) {
        if ( !spin_adapt_q2_ ) {
            Q2_constraints_ATu(A,u);
        }else {
            Q2_constraints_ATu_spin_adapted(A,u);
        }
    }

    if ( constrain_g2_ ) {
        if ( ! spin_adapt_g2_ ) {
            G2_constraints_ATu(A,u);
        }else {
            G2_constraints_ATu_spin_adapted(A,u);
        }
    }

    if ( constrain_t1_ ) {
        T1_constraints_ATu(A,u);
    }

    if ( constrain_t2_ ) {
        //T2_constraints_ATu(A,u);
        T2_constraints_ATu_slow(A,u);
    }

    if ( constrain_d3_ ) {
        D3_constraints_ATu(A,u);
    }

}//end ATu

void v2RDMSolver::cg_Ax(long int N,SharedVector A,SharedVector ux){

    A->zero();
    bpsdp_ATu(ATy,ux);
    bpsdp_Au(A,ATy);

}//end cg_Ax

void v2RDMSolver::Update_xz() {

    // evaluate M(mu*x + ATy - c)
    bpsdp_ATu(ATy,y);
    ATy->subtract(c);
    x->scale(mu);
    ATy->add(x);

    // loop over each block of x/z
    for (int i = 0; i < dimensions_.size(); i++) {
        if ( dimensions_[i] == 0 ) continue;
        int myoffset = 0;
        for (int j = 0; j < i; j++) {
            myoffset += dimensions_[j] * dimensions_[j];
        }

        boost::shared_ptr<Matrix> mat    (new Matrix(dimensions_[i],dimensions_[i]));
        boost::shared_ptr<Matrix> eigvec (new Matrix(dimensions_[i],dimensions_[i]));
        boost::shared_ptr<Vector> eigval (new Vector(dimensions_[i]));
        boost::shared_ptr<Vector> Up     (new Vector(dimensions_[i]));
        boost::shared_ptr<Vector> Um     (new Vector(dimensions_[i]));
        double ** mat_p = mat->pointer();
        double * A_p   = ATy->pointer();

        //C_DCOPY(dimensions_[i]*dimensions_[i],&A_p[myoffset],1,&mat_p[0][0],1);
        for (int p = 0; p < dimensions_[i]; p++) {
            for (int q = p; q < dimensions_[i]; q++) {
                double dum = 0.5 * ( A_p[myoffset + p * dimensions_[i] + q] +
                                     A_p[myoffset + q * dimensions_[i] + p] );
                mat_p[p][q] = mat_p[q][p] = dum;
                 
            }
        }
        mat->diagonalize(eigvec,eigval);

        // separate U+ and U-
        double * u_p    = Up->pointer();
        double * u_m    = Um->pointer();
        double * eval_p = eigval->pointer();
        for (int p = 0; p < dimensions_[i]; p++) {
            if ( eval_p[p] < 0.0 ) {
                u_m[p] = -eval_p[p];
                u_p[p] = 0.0;
            }else {
                u_m[p] = 0.0;
                u_p[p] = eval_p[p]/mu;
            }
        }

        // transform U+ and U- back to nondiagonal basis
        double ** evec_p = eigvec->pointer();
        double * x_p = x->pointer();
        double * z_p = z->pointer();
        #pragma omp parallel for schedule (dynamic)
        for (int pq = 0; pq < dimensions_[i] * dimensions_[i]; pq++) {

            int q = pq % dimensions_[i];
            int p = (pq-q) / dimensions_[i];
            if ( p > q ) continue;

            double sumx = 0.0;
            double sumz = 0.0;
            for (int j = 0; j < dimensions_[i]; j++) {
                sumx += u_p[j] * evec_p[p][j] * evec_p[q][j];
                sumz += u_m[j] * evec_p[p][j] * evec_p[q][j];
            }
            x_p[myoffset+p*dimensions_[i]+q] = sumx;
            z_p[myoffset+p*dimensions_[i]+q] = sumz;

            x_p[myoffset+q*dimensions_[i]+p] = sumx;
            z_p[myoffset+q*dimensions_[i]+p] = sumz;
        }
    }
}

void v2RDMSolver::UnpackDensityPlusCore() {

    memset((void*)d2_plus_core_sym_,'\0',d2_plus_core_dim_*sizeof(double));

    // D2 first
    double * x_p = x->pointer();
    // active active; active active
    int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i            = bas_ab_sym[h][ij][0];
            int j            = bas_ab_sym[h][ij][1];
            int hi           = symmetry[i];
            int hj           = symmetry[j];
            int ifull        = i - pitzer_offset[hi] + pitzer_offset_full[hi] + frzcpi_[hi];
            int jfull        = j - pitzer_offset[hj] + pitzer_offset_full[hj] + frzcpi_[hj];
            int ij_ab        = ibas_ab_sym[h][i][j];
            int ji_ab        = ibas_ab_sym[h][j][i];
            int ij_aa        = ibas_aa_sym[h][i][j];

            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k          = bas_ab_sym[h][kl][0];
                int l          = bas_ab_sym[h][kl][1];
                int hk         = symmetry[k];
                int hl         = symmetry[l];
                int kfull      = k - pitzer_offset[hk] + pitzer_offset_full[hk] + frzcpi_[hk];
                int lfull      = l - pitzer_offset[hl] + pitzer_offset_full[hl] + frzcpi_[hl];
                int kl_ab      = ibas_ab_sym[h][k][l];
                int lk_ab      = ibas_ab_sym[h][l][k];
                int kl_aa      = ibas_aa_sym[h][k][l];

                //if ( i > k ) continue;
                //if ( j > l ) continue;

                int hik = SymmetryPair(hi,hk);

                int ik_full      = ibas_full_sym[hik][ifull][kfull];
                int jl_full      = ibas_full_sym[hik][jfull][lfull];

                //if ( ik_plus_core > jl_plus_core ) continue;

                int hkj = SymmetryPair(hk,hj);

                int kj_ab = ibas_ab_sym[hkj][k][j];
                int il_ab = ibas_ab_sym[hkj][i][l];

                int jk_ab = ibas_ab_sym[hkj][j][k];
                int li_ab = ibas_ab_sym[hkj][l][i];

                int kj_aa = ibas_aa_sym[hkj][k][j];
                int il_aa = ibas_aa_sym[hkj][i][l];

                offset = 0;
                for (int myh = 0; myh < hik; myh++) {
                    offset += gems_plus_core[myh] * ( gems_plus_core[myh] + 1 ) / 2;
                }
                int id = offset + INDEX(ik_full,jl_full);

                double val = 0.0;

                val += 0.5 * x_p[d2aboff[h]   + ij_ab*gems_ab[h]  + kl_ab];
                val += 0.5 * x_p[d2aboff[hkj] + kj_ab*gems_ab[hkj]+ il_ab] * (1.0 - (double)(l==j));
                val += 0.5 * x_p[d2aboff[hkj] + il_ab*gems_ab[hkj]+ kj_ab] * (1.0 - (double)(i==k));
                val += 0.5 * x_p[d2aboff[h]   + kl_ab*gems_ab[h]  + ij_ab] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                val += 0.5 * x_p[d2aboff[h]   + ji_ab*gems_ab[h]  + lk_ab];
                val += 0.5 * x_p[d2aboff[hkj] + jk_ab*gems_ab[hkj]+ li_ab] * (1.0 - (double)(l==j));
                val += 0.5 * x_p[d2aboff[hkj] + li_ab*gems_ab[hkj]+ jk_ab] * (1.0 - (double)(i==k));
                val += 0.5 * x_p[d2aboff[h]   + lk_ab*gems_ab[h]  + ji_ab] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));

                // aa / bb
                if ( i != j && k != l ) {
                    int sij = ( i < j ? 1 : -1 );
                    int skl = ( k < l ? 1 : -1 );
                    val += 0.5 * sij * skl * x_p[d2aaoff[h]   + ij_aa*gems_aa[h]  + kl_aa];
                    val += 0.5 * sij * skl * x_p[d2aaoff[h]   + kl_aa*gems_aa[h]  + ij_aa] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                    val += 0.5 * sij * skl * x_p[d2bboff[h]   + ij_aa*gems_aa[h]  + kl_aa];
                    val += 0.5 * sij * skl * x_p[d2bboff[h]   + kl_aa*gems_aa[h]  + ij_aa] * (1.0 - (double)(l==j))*(1.0-(double)(i==k));
                }
                if ( k != j && i != l ) {
                    int skj = ( k < j ? 1 : -1 );
                    int sil = ( i < l ? 1 : -1 );
                    val += 0.5 * skj * sil * x_p[d2aaoff[hkj] + kj_aa*gems_aa[hkj]+ il_aa] * (1.0 - (double)(l==j));
                    val += 0.5 * skj * sil * x_p[d2aaoff[hkj] + il_aa*gems_aa[hkj]+ kj_aa] * (1.0 - (double)(i==k));
                    val += 0.5 * skj * sil * x_p[d2bboff[hkj] + kj_aa*gems_aa[hkj]+ il_aa] * (1.0 - (double)(l==j));
                    val += 0.5 * skj * sil * x_p[d2bboff[hkj] + il_aa*gems_aa[hkj]+ kj_aa] * (1.0 - (double)(i==k));
                }

                // scale the off-diagonal elements
                if ( ik_full != jl_full ) {
                    val *= 2.0;
                }
                d2_plus_core_sym_[id] = val;
            }
        }
    }

    // core core; core core
    double en = 0.0;
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int i = 0; i < frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];

            for (int hj = 0; hj < nirrep_; hj++) {
                for (int j = 0; j < frzcpi_[hj]; j++) {

                    int jfull      = j + pitzer_offset_full[hj];
                    int hij = SymmetryPair(hi,hj);

                    if ( ifull == jfull ) {

                        int iifull = ibas_full_sym[0][ifull][ifull];
                        int jjfull = ibas_full_sym[0][jfull][jfull];

                        //if ( fabs(d2_plus_core_sym_[INDEX(iifull,jjfull)]) < 1e-10 ) {
                        //    en += tei_full_sym_[INDEX(iifull,jjfull)];
                        //}

                        d2_plus_core_sym_[INDEX(iifull,jjfull)] =  1.0;

                    }else {

                        int iifull = ibas_full_sym[0][ifull][ifull];
                        int jjfull = ibas_full_sym[0][jfull][jfull];

                        //if ( fabs(d2_plus_core_sym_[INDEX(iifull,jjfull)]) < 1e-10 ) {
                        //    en += 4.0 * tei_full_sym_[INDEX(iifull,jjfull)];
                        //}

                        d2_plus_core_sym_[INDEX(iifull,jjfull)] =  4.0;

                        int offset2 = 0;
                        for (int myh = 0; myh < hij; myh++) {
                            offset2 += gems_full[myh] * ( gems_full[myh] + 1 ) / 2;
                        }
                        offset = 0;
                        for (int myh = 0; myh < hij; myh++) {
                            offset += gems_plus_core[myh] * ( gems_plus_core[myh] + 1 ) / 2;
                        }

                        int ijfull = ibas_full_sym[hij][ifull][jfull];
                        //if ( fabs(d2_plus_core_sym_[offset + INDEX(ijfull,ijfull)]) < 1e-10 ) {
                        //    en -= 2.0 * tei_full_sym_[offset2 + INDEX(ijfull,ijfull)];
                        //}

                        d2_plus_core_sym_[offset + INDEX(ijfull,ijfull)] = -2.0;
                    }
                }
            }
        }
    }

    // core active; core active
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int i = 0; i < frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];
            int iifull     = ibas_full_sym[0][ifull][ifull];

            // D2(il; ij) ab, ba, aa, bb
            for (int hj = 0; hj < nirrep_; hj++) {
                for (int j = 0; j < amopi_[hj]; j++) {

                    int jfull      = j + pitzer_offset_full[hj] + frzcpi_[hj];

                    for (int l = j; l < amopi_[hj]; l++) {

                        int lfull      = l + pitzer_offset_full[hj] + frzcpi_[hj];

                        int jlfull = ibas_full_sym[0][jfull][lfull];

                        int id = INDEX(iifull,jlfull);

                        double val = 0.0;

                        // aa and bb pieces
                        if ( j == l ) {
                            val += 1.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val += 1.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }else {
                            val += 2.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val += 2.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }

                        // ab and ba pieces
                        if ( j == l ) {
                            val += 1.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val += 1.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }else {
                            val += 2.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val += 2.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }

                        d2_plus_core_sym_[id] = val;

                        // also (il|ji) with a minus sign
                        int hil = SymmetryPair(symmetry_full[ifull],symmetry_full[lfull]);
                        int hji = SymmetryPair(symmetry_full[ifull],symmetry_full[jfull]);

                        int ilfull = ibas_full_sym[hil][ifull][lfull];
                        int jifull = ibas_full_sym[hji][jfull][ifull];

                        //if ( il > ji ) continue;

                        val = 0.0;

                        // aa and bb pieces
                        if ( j == l ) {
                            val -= 1.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val -= 1.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }else {
                            val -= 2.0 * x_p[d1aoff[hj]+j*amopi_[hj]+l];
                            val -= 2.0 * x_p[d1boff[hj]+j*amopi_[hj]+l];
                        }

                        offset = 0;
                        for (int myh = 0; myh < hil; myh++) {
                            offset += gems_plus_core[myh] * ( gems_plus_core[myh] + 1 ) / 2;
                        }

                        d2_plus_core_sym_[offset + INDEX(ilfull,jifull)] = val;

                    }
                }
            }
        }
    }

    // now D1
    // active; active
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {

            int iplus_core = i + frzcpi_[h];

            for (int j = i; j < amopi_[h]; j++) {

                int jplus_core = j + frzcpi_[h];

                int id = offset + INDEX(iplus_core,jplus_core);
                
                d1_plus_core_sym_[id]  = x_p[d1aoff[h] + i * amopi_[h] + j];
                d1_plus_core_sym_[id] += x_p[d1boff[h] + i * amopi_[h] + j];

                // scale off-diagonal elements
                if ( i != j ) {
                    d1_plus_core_sym_[id] *= 2.0;
                }

            }
        }
        offset += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }

    // core; core;
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            d1_plus_core_sym_[offset + INDEX(i,i)] = 2.0;
        }
        offset += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }
}

// repack rotated full-space integrals into active-space integrals
void v2RDMSolver::RepackIntegralsDF(){

    // if frozen core, adjust oei's and compute frozen core energy:
    efzc_ = 0.0;
    offset = 0;
    long int offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < frzcpi_[h]; i++) {
            long int ii = i + offset;
            efzc_ += 2.0 * oei_full_sym_[offset3 + INDEX(i,i)];

            long int offset2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (long int j = 0; j < frzcpi_[h2]; j++) {
                    long int jj = j + offset2;
                    double dum1 = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(ii,ii),1,Qmo_+nQ_*INDEX(jj,jj),1);
                    double dum2 = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(ii,jj),1,Qmo_+nQ_*INDEX(ii,jj),1);
                    efzc_ += 2.0 * dum1 - dum2;
                }
                offset2 += nmopi_[h2];
            }
        }
        offset += nmopi_[h];
        offset3 += nmopi_[h] * (nmopi_[h] + 1 ) / 2;
    }

    double* c_p = c->pointer();

    // adjust oeis
    offset = 0;
    offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = frzcpi_[h]; i < nmopi_[h] - frzvpi_[h]; i++) {
            long int ii = i + offset;
            for (long int j = frzcpi_[h]; j < nmopi_[h] - frzvpi_[h]; j++) {

                long int jj = j + offset;
                double dum1 = 0.0;
                double dum2 = 0.0;

                long int offset2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {

                    for (long int k = 0; k < frzcpi_[h2]; k++) {
                        long int kk = k + offset2;
                        dum1 += C_DDOT(nQ_,Qmo_ + nQ_*INDEX(ii,jj),1,Qmo_+nQ_*INDEX(kk,kk),1);
                        dum2 += C_DDOT(nQ_,Qmo_ + nQ_*INDEX(ii,kk),1,Qmo_+nQ_*INDEX(jj,kk),1);
                    }
                    offset2 += nmopi_[h2];
                }
                c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = oei_full_sym_[offset3+INDEX(i,j)];
                c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = oei_full_sym_[offset3+INDEX(i,j)];

                c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] += 2.0 * dum1 - dum2;
                c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] += 2.0 * dum1 - dum2;
            }
        }
        offset += nmopi_[h];
        offset3 += nmopi_[h] * (nmopi_[h]+1)/2;
    }

    // two-electron part
    long int na = nalpha_ - nfrzc_;
    long int nb = nbeta_ - nfrzc_;
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
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

                c_p[d2aboff[h] + ij*gems_ab[h]+kl] = C_DDOT(nQ_,Qmo_+nQ_*INDEX(ii,kk),1,Qmo_+nQ_*INDEX(jj,ll),1);
            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
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

                double dum1 = C_DDOT(nQ_,Qmo_+nQ_*INDEX(ii,kk),1,Qmo_+nQ_*INDEX(jj,ll),1);
                double dum2 = C_DDOT(nQ_,Qmo_+nQ_*INDEX(ii,ll),1,Qmo_+nQ_*INDEX(jj,kk),1);

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]    = dum1 - dum2;
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]    = dum1 - dum2;
            }
        }
    }

}

// repack rotated full-space integrals into active-space integrals
void v2RDMSolver::RepackIntegrals(){

    // if frozen core, adjust oei's and compute frozen core energy:
    efzc_ = 0.0;
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {

            int ifull = i + pitzer_offset_full[h];
            int ii    = ibas_full_sym[0][ifull][ifull];

            efzc_ += 2.0 * oei_full_sym_[offset + INDEX(i,i)]; 

            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (int j = 0; j < frzcpi_[h2]; j++) {

                    int jfull = j + pitzer_offset_full[h2];
                    int jj    = ibas_full_sym[0][jfull][jfull];

                    int hij = SymmetryPair(h,h2);

                    int ij    = ibas_full_sym[hij][ifull][jfull];

                    int myoff = 0;
                    for (int myh = 0; myh < hij; myh++) {
                        myoff += gems_full[myh] * ( gems_full[myh] + 1 ) / 2;
                    }
                    efzc_ += (2.0 * tei_full_sym_[INDEX(ii,jj)] - tei_full_sym_[myoff + INDEX(ij,ij)]);
                }
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }

    double* c_p = c->pointer();

    // adjust one-electron integrals for core repulsion contribution
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = frzcpi_[h]; i < nmopi_[h] - frzvpi_[h]; i++) {

            int ifull = i + pitzer_offset_full[h];

            for (int j = frzcpi_[h]; j < nmopi_[h] - frzvpi_[h]; j++) {

                int jfull = j + pitzer_offset_full[h];

                int ij = ibas_full_sym[0][ifull][jfull];

                double dum = 0.0;

                for (int h2 = 0; h2 < nirrep_; h2++) {
                    for (int k = 0; k < frzcpi_[h2]; k++) {

                        int kfull = k + pitzer_offset_full[h2];

                        int kk = ibas_full_sym[0][kfull][kfull];

                        int hik = SymmetryPair(h,h2);

                        int ik = ibas_full_sym[hik][ifull][kfull];
                        int jk = ibas_full_sym[hik][jfull][kfull];

                        int myoff = 0;
                        for (int myh = 0; myh < hik; myh++) {
                            myoff += gems_full[myh] * ( gems_full[myh] + 1 ) / 2;
                        }

                        dum += (2.0 * tei_full_sym_[INDEX(ij,kk)] - tei_full_sym_[myoff+INDEX(ik,jk)]);

                    }
                }
                c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = oei_full_sym_[offset+INDEX(i,j)];
                c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = oei_full_sym_[offset+INDEX(i,j)];
                c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] += dum;
                c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] += dum;
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }

    int na = nalpha_ - nfrzc_;
    int nb = nbeta_ - nfrzc_;
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            int hi = symmetry[i];
            int ifull = i - pitzer_offset[hi] + pitzer_offset_full[hi] + frzcpi_[hi];

            int hj = symmetry[j];
            int jfull = j - pitzer_offset[hj] + pitzer_offset_full[hj] + frzcpi_[hj];

            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];

                int hk = symmetry[k];
                int kfull = k - pitzer_offset[hk] + pitzer_offset_full[hk] + frzcpi_[hk];

                int hl = symmetry[l];
                int lfull = l - pitzer_offset[hl] + pitzer_offset_full[hl] + frzcpi_[hl];

                int hik = SymmetryPair(symmetry_full[ifull],symmetry_full[kfull]);
                int ik  = ibas_full_sym[hik][ifull][kfull];
                int jl  = ibas_full_sym[hik][jfull][lfull];

                int myoff = 0;
                for (int myh = 0; myh < hik; myh++) {
                    myoff += gems_full[myh] * ( gems_full[myh] + 1 ) / 2;
                }

                c_p[d2aboff[h] + ij*gems_ab[h]+kl]    = tei_full_sym_[myoff + INDEX(ik,jl)];

            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];

            int hi = symmetry[i];
            int ifull = i - pitzer_offset[hi] + pitzer_offset_full[hi] + frzcpi_[hi];

            int hj = symmetry[j];
            int jfull = j - pitzer_offset[hj] + pitzer_offset_full[hj] + frzcpi_[hj];

            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];

                int hk = symmetry[k];
                int kfull = k - pitzer_offset[hk] + pitzer_offset_full[hk] + frzcpi_[hk];

                int hl = symmetry[l];
                int lfull = l - pitzer_offset[hl] + pitzer_offset_full[hl] + frzcpi_[hl];

                int hik = SymmetryPair(symmetry_full[ifull],symmetry_full[kfull]);
                int ik  = ibas_full_sym[hik][ifull][kfull];
                int jl  = ibas_full_sym[hik][jfull][lfull];

                int myoff = 0;
                for (int myh = 0; myh < hik; myh++) {
                    myoff += gems_full[myh] * ( gems_full[myh] + 1 ) / 2;
                }

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]    = 0.5 * tei_full_sym_[myoff + INDEX(ik,jl)];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   += 0.5 * tei_full_sym_[myoff + INDEX(jl,ik)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]    = 0.5 * tei_full_sym_[myoff + INDEX(ik,jl)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   += 0.5 * tei_full_sym_[myoff + INDEX(jl,ik)];

                int hil = SymmetryPair(symmetry_full[ifull],symmetry_full[lfull]);
                int il  = ibas_full_sym[hil][ifull][lfull];
                int jk  = ibas_full_sym[hil][jfull][kfull];

                myoff = 0;
                for (int myh = 0; myh < hil; myh++) {
                    myoff += gems_full[myh] * ( gems_full[myh] + 1 ) / 2;
                }

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym_[myoff + INDEX(il,jk)];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym_[myoff + INDEX(jk,il)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym_[myoff + INDEX(il,jk)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym_[myoff + INDEX(jk,il)];

            }
        }
    }
}

void v2RDMSolver::FinalTransformationMatrix() {

    // test: transform oeis
    /*offset = 0;
    double diff = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ieo = 0; ieo < nmo_; ieo++) {
            int ifull = energy_to_pitzer_order[ieo];
            int hi    = symmetry_full[ifull];
            if ( h != hi ) continue;
            int i     = ifull - pitzer_offset_full[hi];
            for (int jeo = 0; jeo < nmo_; jeo++) {
                int jfull = energy_to_pitzer_order[jeo];
                int hj    = symmetry_full[jfull];
                if ( h != hj ) continue;
                int j     = jfull - pitzer_offset_full[hj];
                double dum = 0.0;
                for (int keo = 0; keo < nmo_; keo++) {
                    int kfull = energy_to_pitzer_order[keo];
                    int hk    = symmetry_full[kfull];
                    if ( h != hk ) continue;
                    int k     = kfull - pitzer_offset_full[hk];
                    for (int leo = 0; leo < nmo_; leo++) {
                        int lfull = energy_to_pitzer_order[leo];
                        int hl    = symmetry_full[lfull];
                        if ( h != hl ) continue;
                        int l     = lfull - pitzer_offset_full[hl];
                        //dum += saveOEI_->pointer(h)[k][l] * orbopt_transformation_matrix_[ieo*nmo_+keo]
                        //                                * orbopt_transformation_matrix_[jeo*nmo_+leo];
                        dum += saveOEI_->pointer(h)[k][l] * orbopt_transformation_matrix_[keo*nmo_+ieo]
                                                        * orbopt_transformation_matrix_[leo*nmo_+jeo];
                    }
                }
                diff += (dum - oei_full_sym_[offset+INDEX(i,j)]) * (dum - oei_full_sym_[offset+INDEX(i,j)]);
                printf("%5i %5i %20.12lf %20.12lf\n",i,j,dum,oei_full_sym_[offset+INDEX(i,j)]);
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
    printf("%20.12lf\n",diff);*/


    // update so/mo coefficient matrix (only need Ca_):
    for (int h = 0; h < nirrep_; h++) {
        double **ca_p = Ca_->pointer(h);
        double **cb_p = Cb_->pointer(h);
        for (int mu = 0; mu < nsopi_[h]; mu++) {
            double * temp = (double*)malloc(nmopi_[h] * sizeof(double));

            // new basis function i in energy order
            for (int ieo = 0; ieo < nmo_; ieo++) {
                int ifull = energy_to_pitzer_order[ieo];
                int hi    = symmetry_full[ifull];
                if ( h != hi ) continue;
                int i     = ifull - pitzer_offset_full[hi];

                double dum = 0.0;

                // old basis function j in energy order
                for (int jeo = 0; jeo < nmo_; jeo++) {
                    int jfull = energy_to_pitzer_order[jeo];
                    int hj    = symmetry_full[jfull];
                    if ( h != hj ) continue;
                    int j     = jfull - pitzer_offset_full[hj];

                    dum += ca_p[mu][j] * orbopt_transformation_matrix_[jeo*nmo_+ieo];
                }
                temp[i] = dum;
            }
            for (int i = 0; i < nmopi_[h]; i++) {
                ca_p[mu][i] = temp[i];
                cb_p[mu][i] = temp[i];
            }
            free(temp);
        }
    }
    //for (int i = 0; i < nmo_; i++) {
    //    for (int j = 0; j < nmo_; j++) {
    //        printf("%5i %5i %20.12lf\n",i,j,orbopt_transformation_matrix_[i*nmo_+j]);
    //    }
    //}
    //Ca_->print();
}

void v2RDMSolver::RotateOrbitals(){

/*   
    if ( is_df_ ) {
        throw PsiException("orbital optimization does not work with 3-index integrals yet",__FILE__,__LINE__);
    }
*/

    UnpackDensityPlusCore();

    outfile->Printf("\n");
    outfile->Printf("        ==> Orbital Optimization <==\n");
    outfile->Printf("\n");

    int nfrzv = nmo_-amo_-nfrzc_;
    OrbOpt(orbopt_transformation_matrix_,
          oei_full_sym_,oei_full_dim_,tei_full_sym_,tei_full_dim_,
          d1_plus_core_sym_,d1_plus_core_dim_,d2_plus_core_sym_,d2_plus_core_dim_,
          symmetry_energy_order,nfrzc_,amo_,nfrzv,nirrep_,
          orbopt_data_,orbopt_outfile_);

    outfile->Printf("            Orbital Optimization %s in %3i iterations \n",(int)orbopt_data_[12] ? "converged" : "did not converge",(int)orbopt_data_[9]);
    outfile->Printf("            Total energy change: %11.6le\n",orbopt_data_[11]);
    outfile->Printf("            Final gradient norm: %11.6le\n",orbopt_data_[10]);
    outfile->Printf("\n");

    if ( fabs(orbopt_data_[11]) < orbopt_data_[4] ) {
        orbopt_converged_ = true;
    }

    if ( is_df_ ) {
        RepackIntegralsDF();
    }else {
        RepackIntegrals();
    }
}

}} //end namespaces
