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
bool prints2 = false;
namespace psi{ namespace v2rdm_casscf{

v2RDMSolver::v2RDMSolver(boost::shared_ptr<Wavefunction> reference_wavefunction,Options & options):
    Wavefunction(options,_default_psio_lib_){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

v2RDMSolver::~v2RDMSolver()
{
}

int v2RDMSolver::SymmetryPair(int i,int j) {
    return table[i*8+j];
}
int v2RDMSolver::TotalSym(int i,int j,int k, int l) {
    return SymmetryPair(SymmetryPair(symmetry[i],symmetry[j]),SymmetryPair(symmetry[k],symmetry[l]));
}


void v2RDMSolver::BuildBasis() {

    // product table:
    table = (int*)malloc(64*sizeof(int));
    memset((void*)table,'\0',64*sizeof(int));
    table[0*8+1] = table[1*8+0] = 1;
    table[0*8+2] = table[2*8+0] = 2;
    table[0*8+3] = table[3*8+0] = 3;
    table[0*8+4] = table[4*8+0] = 4;
    table[0*8+5] = table[5*8+0] = 5;
    table[0*8+6] = table[6*8+0] = 6;
    table[0*8+7] = table[7*8+0] = 7;
    table[1*8+2] = table[2*8+1] = 3;
    table[1*8+3] = table[3*8+1] = 2;
    table[1*8+4] = table[4*8+1] = 5;
    table[1*8+5] = table[5*8+1] = 4;
    table[1*8+6] = table[6*8+1] = 7;
    table[1*8+7] = table[7*8+1] = 6;
    table[2*8+3] = table[3*8+2] = 1;
    table[2*8+4] = table[4*8+2] = 6;
    table[2*8+5] = table[5*8+2] = 7;
    table[2*8+6] = table[6*8+2] = 4;
    table[2*8+7] = table[7*8+2] = 5;
    table[3*8+4] = table[4*8+3] = 7;
    table[3*8+5] = table[5*8+3] = 6;
    table[3*8+6] = table[6*8+3] = 5;
    table[3*8+7] = table[7*8+3] = 4;
    table[4*8+5] = table[5*8+4] = 1;
    table[4*8+6] = table[6*8+4] = 2;
    table[4*8+7] = table[7*8+4] = 3;
    table[5*8+6] = table[6*8+5] = 3;
    table[5*8+7] = table[7*8+5] = 2;
    table[6*8+7] = table[7*8+6] = 1;

    // orbitals are in pitzer order:
    symmetry               = (int*)malloc((nmo+nfrzc+nfrzv)*sizeof(int));
    symmetry_full          = (int*)malloc((nmo+nfrzc+nfrzv)*sizeof(int));
    symmetry_plus_core     = (int*)malloc((nmo+nfrzc)*sizeof(int));
    symmetry_energy_order  = (myint*)malloc((nmo+nfrzc+nfrzv)*sizeof(myint));
    pitzer_to_energy_order = (myint*)malloc((nmo+nfrzc+nfrzv)*sizeof(myint));
    energy_to_pitzer_order = (myint*)malloc((nmo+nfrzc+nfrzv)*sizeof(myint));
    memset((void*)symmetry,'\0',(nmo+nfrzc+nfrzv)*sizeof(int));
    memset((void*)symmetry_full,'\0',(nmo+nfrzc+nfrzv)*sizeof(int));
    memset((void*)symmetry_plus_core,'\0',(nmo+nfrzc)*sizeof(int));
    memset((void*)symmetry_energy_order,'\0',(nmo+nfrzc+nfrzv)*sizeof(myint));
    memset((void*)pitzer_to_energy_order,'\0',(nmo+nfrzc+nfrzv)*sizeof(myint));
    memset((void*)energy_to_pitzer_order,'\0',(nmo+nfrzc+nfrzv)*sizeof(myint));
    full_basis = (int*)malloc((nmo+nfrzc+nfrzv)*sizeof(int));
    int count = 0;
    int count_full = 0;

    // symmetry of ACTIVE orbitals
    for (int h = 0; h < nirrep_; h++) {
        count_full += frzcpi_[h];
        for (int norb = frzcpi_[h]; norb < nmopi_[h] - frzvpi_[h]; norb++){
            full_basis[count] = count_full++;
            symmetry[count++] = h;
        }
        count_full += frzvpi_[h];
    }

    // symmetry of ALL orbitals
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int norb = 0; norb < nmopi_[h]; norb++){
            symmetry_full[count++] = h;
        }
    }
    // symmetry of ALL orbitals, except frozen virtual
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int norb = 0; norb < amopi_[h] + frzcpi_[h]; norb++){
            symmetry_plus_core[count++] = h;
        }
    }
    // symmetry of ALL orbitals in energy order
    double min = 1.0e99;
    int imin = -999;
    int isym = -999;
    int * skip = (int*)malloc((nmo+nfrzc+nfrzv)*sizeof(int));
    memset((void*)skip,'\0',(nmo+nfrzc+nfrzv)*sizeof(int));

    // warning to future eugene:  it is possible that
    // this ordering will differ from that printed at the
    // end of the SCF routine if orbitals are truly
    // degenerate.  past eugene hasn't convinved himself
    // of whether or not this is actually a problem.

    // TODO: the orbital ordering should be according to 
    // energy within each type of orbital

    // core
    for (int i = 0; i < nfrzc; i++){

        int me = 0;
        min = 1.0e99;
        for (int h = 0; h < nirrep_; h++) {
            for (int j = 0; j < frzcpi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += nmopi_[h] - frzcpi_[h];
        }
        skip[imin] = 1;
        symmetry_energy_order[i] = isym + 1;
        pitzer_to_energy_order[imin] = i;
        energy_to_pitzer_order[i] = imin;
    }
    // active
    for (int i = nfrzc; i < nmo + nfrzc; i++){

        int me = 0;
        min = 1.0e99;
        for (int h = 0; h < nirrep_; h++) {
            me += frzcpi_[h];
            for (int j = frzcpi_[h]; j < frzcpi_[h]+amopi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += frzvpi_[h];
        }
        skip[imin] = 1;
        symmetry_energy_order[i] = isym + 1;
        pitzer_to_energy_order[imin] = i;
        energy_to_pitzer_order[i] = imin;
    }
    // virtual
    for (int i = nmo + nfrzc; i < nfrzv + nmo + nfrzc; i++){

        int me = 0;
        min = 1.0e99;
        for (int h = 0; h < nirrep_; h++) {
            me += frzcpi_[h] + amopi_[h];
            for (int j = frzcpi_[h] + amopi_[h]; j < nmopi_[h]; j++){
                if ( epsilon_a_->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = epsilon_a_->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
        }
        skip[imin] = 1;
        symmetry_energy_order[i] = isym + 1;
        pitzer_to_energy_order[imin] = i;
        energy_to_pitzer_order[i] = imin;
    }

    pitzer_offset           = (int*)malloc(nirrep_*sizeof(int));
    pitzer_offset_full      = (int*)malloc(nirrep_*sizeof(int));
    pitzer_offset_plus_core = (int*)malloc(nirrep_*sizeof(int));
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        pitzer_offset[h] = count;
        count += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        pitzer_offset_full[h] = count;
        count += nmopi_[h];
    }
    count = 0;
    for (int h = 0; h < nirrep_; h++) {
        pitzer_offset_plus_core[h] = count;
        count += nmopi_[h] - frzvpi_[h];
    }

    // symmetry pairs:
    int ** sympairs = (int**)malloc(nirrep_*sizeof(int*));
    for (int h = 0; h < nirrep_; h++) {
        sympairs[h] = (int*)malloc(nmo*sizeof(int));
        memset((void*)sympairs[h],'\0',nmo*sizeof(int));
    }

    for (int h = 0; h < nirrep_; h++) {
        std::vector < std::pair<int,int> > mygems;
        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                int sym = SymmetryPair(symmetry[i],symmetry[j]);
                if (h==sym) {
                    mygems.push_back(std::make_pair(j,i));
                }

            }
        }
        gems.push_back(mygems);
    }
    for (int h = 0; h < nirrep_; h++) {
        std::vector < std::pair<int,int> > mygems;
        for (int i = 0; i < nmo + nfrzc + nfrzv; i++) {
            for (int j = 0; j <= i; j++) {
                int sym = SymmetryPair(symmetry_full[i],symmetry_full[j]);
                if (h==sym) {
                    mygems.push_back(std::make_pair(i,j));
                }

            }
        }
        gems_fullspace.push_back(mygems);
    }
    for (int h = 0; h < nirrep_; h++) {
        std::vector < std::pair<int,int> > mygems;
        for (int i = 0; i < nmo + nfrzc; i++) {
            for (int j = 0; j <= i; j++) {
                int sym = SymmetryPair(symmetry_plus_core[i],symmetry_plus_core[j]);
                if (h==sym) {
                    mygems.push_back(std::make_pair(i,j));
                }

            }
        }
        gems_plus_corespace.push_back(mygems);
    }

    bas_ab_sym         = (int***)malloc(nirrep_*sizeof(int**));
    bas_aa_sym         = (int***)malloc(nirrep_*sizeof(int**));
    bas_00_sym         = (int***)malloc(nirrep_*sizeof(int**));
    bas_full_sym       = (int***)malloc(nirrep_*sizeof(int**));

    ibas_ab_sym        = (int***)malloc(nirrep_*sizeof(int**));
    ibas_aa_sym        = (int***)malloc(nirrep_*sizeof(int**));
    ibas_00_sym        = (int***)malloc(nirrep_*sizeof(int**));
    ibas_full_sym      = (int***)malloc(nirrep_*sizeof(int**));

    gems_ab            = (int*)malloc(nirrep_*sizeof(int));
    gems_aa            = (int*)malloc(nirrep_*sizeof(int));
    gems_00            = (int*)malloc(nirrep_*sizeof(int));
    gems_full          = (int*)malloc(nirrep_*sizeof(int));
    gems_plus_core     = (int*)malloc(nirrep_*sizeof(int));

    for (int h = 0; h < nirrep_; h++) {

        ibas_ab_sym[h]        = (int**)malloc(nmo*sizeof(int*));
        ibas_aa_sym[h]        = (int**)malloc(nmo*sizeof(int*));
        ibas_00_sym[h]        = (int**)malloc(nmo*sizeof(int*));
        ibas_full_sym[h]      = (int**)malloc((nmo+nfrzc+nfrzv)*sizeof(int*));

        bas_ab_sym[h]         = (int**)malloc(nmo*nmo*sizeof(int*));
        bas_aa_sym[h]         = (int**)malloc(nmo*nmo*sizeof(int*));
        bas_00_sym[h]         = (int**)malloc(nmo*nmo*sizeof(int*));
        bas_full_sym[h]       = (int**)malloc((nmo+nfrzc+nfrzv)*(nmo+nfrzc+nfrzv)*sizeof(int*));

        // active space geminals
        for (int i = 0; i < nmo; i++) {
            ibas_ab_sym[h][i] = (int*)malloc(nmo*sizeof(int));
            ibas_aa_sym[h][i] = (int*)malloc(nmo*sizeof(int));
            ibas_00_sym[h][i] = (int*)malloc(nmo*sizeof(int));
            for (int j = 0; j < nmo; j++) {
                ibas_ab_sym[h][i][j] = -999;
                ibas_aa_sym[h][i][j] = -999;
                ibas_00_sym[h][i][j] = -999;
            }
        }
        for (int i = 0; i < nmo*nmo; i++) {
            bas_ab_sym[h][i] = (int*)malloc(2*sizeof(int));
            bas_aa_sym[h][i] = (int*)malloc(2*sizeof(int));
            bas_00_sym[h][i] = (int*)malloc(2*sizeof(int));
            for (int j = 0; j < 2; j++) {
                bas_ab_sym[h][i][j] = -999;
                bas_aa_sym[h][i][j] = -999;
                bas_00_sym[h][i][j] = -999;
            }
        }
        // full space geminals
        for (int i = 0; i < nmo+nfrzv+nfrzc; i++) {
            ibas_full_sym[h][i] = (int*)malloc((nmo+nfrzc+nfrzv)*sizeof(int));
            for (int j = 0; j < nfrzv+nfrzc+nmo; j++) {
                ibas_full_sym[h][i][j] = -999;
            }
        }
        for (int i = 0; i < (nfrzc+nfrzv+nmo)*(nfrzc+nfrzv+nmo); i++) {
            bas_full_sym[h][i] = (int*)malloc(2*sizeof(int));
            for (int j = 0; j < 2; j++) {
                bas_full_sym[h][i][j] = -999;
            }
        }

        // active space mappings:
        int count_ab = 0;
        int count_aa = 0;
        int count_00 = 0;
        for (int n = 0; n < gems[h].size(); n++) {
            int i = gems[h][n].first;
            int j = gems[h][n].second;

            ibas_ab_sym[h][i][j] = n;
            bas_ab_sym[h][n][0]  = i;
            bas_ab_sym[h][n][1]  = j;
            count_ab++;

            if ( i < j ) continue;

            ibas_00_sym[h][i][j] = count_00;
            ibas_00_sym[h][j][i] = count_00;
            bas_00_sym[h][count_00][0] = i;
            bas_00_sym[h][count_00][1] = j;
            count_00++;

            if ( i <= j ) continue;

            ibas_aa_sym[h][i][j] = count_aa;
            ibas_aa_sym[h][j][i] = count_aa;
            bas_aa_sym[h][count_aa][0] = i;
            bas_aa_sym[h][count_aa][1] = j;
            count_aa++;
        }
        gems_ab[h] = count_ab;
        gems_aa[h] = count_aa;
        gems_00[h] = count_00;

    }

    // new way:
    memset((void*)gems_full,'\0',nirrep_*sizeof(int));
    memset((void*)gems_plus_core,'\0',nirrep_*sizeof(int));

    for (int ieo = 0; ieo < nmo + nfrzc + nfrzv; ieo++) {
        int ifull = energy_to_pitzer_order[ieo];
        int hi    = symmetry_full[ifull];
        int i     = ifull - pitzer_offset_full[hi];
        for (int jeo = 0; jeo <= ieo; jeo++) {
            int jfull = energy_to_pitzer_order[jeo];
            int hj    = symmetry_full[jfull];
            int j     = jfull - pitzer_offset_full[hj];
           
            int hij = SymmetryPair(hi,hj); 
            ibas_full_sym[hij][ifull][jfull] = gems_full[hij];
            ibas_full_sym[hij][jfull][ifull] = gems_full[hij];
            bas_full_sym[hij][gems_full[hij]][0] = ifull;
            bas_full_sym[hij][gems_full[hij]][1] = jfull;
            gems_full[hij]++;
            if ( ieo < nmo + nfrzc && jeo < nmo + nfrzc ) {
                gems_plus_core[hij]++;
            }
        }
    }

    if ( constrain_t1 || constrain_t2 || constrain_d3 ) {
        // make all triplets
        for (int h = 0; h < nirrep_; h++) {
            std::vector < boost::tuple<int,int,int> > mytrip;
            for (int i = 0; i < nmo; i++) {
                for (int j = 0; j < nmo; j++) {
                    int s1 = SymmetryPair(symmetry[i],symmetry[j]);
                    for (int k = 0; k < nmo; k++) {
                        int s2 = SymmetryPair(s1,symmetry[k]);
                        if (h==s2) {
                            mytrip.push_back(boost::make_tuple(i,j,k));
                        }
                    }

                }
            }
            triplets.push_back(mytrip);
        }
        bas_aaa_sym  = (int***)malloc(nirrep_*sizeof(int**));
        bas_aab_sym  = (int***)malloc(nirrep_*sizeof(int**));
        bas_aba_sym  = (int***)malloc(nirrep_*sizeof(int**));
        ibas_aaa_sym = (int****)malloc(nirrep_*sizeof(int***));
        ibas_aab_sym = (int****)malloc(nirrep_*sizeof(int***));
        ibas_aba_sym = (int****)malloc(nirrep_*sizeof(int***));
        trip_aaa    = (int*)malloc(nirrep_*sizeof(int));
        trip_aab    = (int*)malloc(nirrep_*sizeof(int));
        trip_aba    = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            ibas_aaa_sym[h] = (int***)malloc(nmo*sizeof(int**));
            ibas_aab_sym[h] = (int***)malloc(nmo*sizeof(int**));
            ibas_aba_sym[h] = (int***)malloc(nmo*sizeof(int**));
            bas_aaa_sym[h]  = (int**)malloc(nmo*nmo*nmo*sizeof(int*));
            bas_aab_sym[h]  = (int**)malloc(nmo*nmo*nmo*sizeof(int*));
            bas_aba_sym[h]  = (int**)malloc(nmo*nmo*nmo*sizeof(int*));
            for (int i = 0; i < nmo; i++) {
                ibas_aaa_sym[h][i] = (int**)malloc(nmo*sizeof(int*));
                ibas_aab_sym[h][i] = (int**)malloc(nmo*sizeof(int*));
                ibas_aba_sym[h][i] = (int**)malloc(nmo*sizeof(int*));
                for (int j = 0; j < nmo; j++) {
                    ibas_aaa_sym[h][i][j] = (int*)malloc(nmo*sizeof(int));
                    ibas_aab_sym[h][i][j] = (int*)malloc(nmo*sizeof(int));
                    ibas_aba_sym[h][i][j] = (int*)malloc(nmo*sizeof(int));
                    for (int k = 0; k < nmo; k++) {
                        ibas_aaa_sym[h][i][j][k] = -999;
                        ibas_aab_sym[h][i][j][k] = -999;
                        ibas_aba_sym[h][i][j][k] = -999;
                    }
                }
            }
            for (int i = 0; i < nmo*nmo*nmo; i++) {
                bas_aaa_sym[h][i] = (int*)malloc(3*sizeof(int));
                bas_aab_sym[h][i] = (int*)malloc(3*sizeof(int));
                bas_aba_sym[h][i] = (int*)malloc(3*sizeof(int));
                for (int j = 0; j < 3; j++) {
                    bas_aaa_sym[h][i][j] = -999;
                    bas_aab_sym[h][i][j] = -999;
                    bas_aba_sym[h][i][j] = -999;
                }
            }

            // mappings:
            int count_aaa = 0;
            int count_aab = 0;
            int count_aba = 0;
            for (int n = 0; n < triplets[h].size(); n++) {
                int i = get<0>(triplets[h][n]);
                int j = get<1>(triplets[h][n]);
                int k = get<2>(triplets[h][n]);

                ibas_aba_sym[h][i][j][k] = count_aba;
                bas_aba_sym[h][count_aba][0]  = i;
                bas_aba_sym[h][count_aba][1]  = j;
                bas_aba_sym[h][count_aba][2]  = k;
                count_aba++;

                if ( i >= j ) continue;

                ibas_aab_sym[h][i][j][k] = count_aab;
                ibas_aab_sym[h][j][i][k] = count_aab;
                bas_aab_sym[h][count_aab][0]  = i;
                bas_aab_sym[h][count_aab][1]  = j;
                bas_aab_sym[h][count_aab][2]  = k;
                count_aab++;

                if ( j >= k ) continue;

                ibas_aaa_sym[h][i][j][k] = count_aaa;
                ibas_aaa_sym[h][i][k][j] = count_aaa;
                ibas_aaa_sym[h][j][i][k] = count_aaa;
                ibas_aaa_sym[h][j][k][i] = count_aaa;
                ibas_aaa_sym[h][k][i][j] = count_aaa;
                ibas_aaa_sym[h][k][j][i] = count_aaa;
                bas_aaa_sym[h][count_aaa][0]  = i;
                bas_aaa_sym[h][count_aaa][1]  = j;
                bas_aaa_sym[h][count_aaa][2]  = k;
                count_aaa++;

            }
            trip_aaa[h] = count_aaa;
            trip_aab[h] = count_aab;
            trip_aba[h] = count_aba;
        }
    }
}

// compute the energy!
double v2RDMSolver::compute_energy() {

    PrintHeader();
    Guess();

    if ( is_df_ ) {
        DFK2();
    }else {
        K2();
    }

    BuildConstraints();
    
    // AATy = A(c-z)+tu(b-Ax) rearange w.r.t cg solver
    // Ax   = AATy and b=A(c-z)+tu(b-Ax)
    SharedVector B   = SharedVector(new Vector("compound B",nconstraints));

    tau = 1.6;
    mu  = 1.0;

    // congugate gradient solver
    long int N = nconstraints;
    shared_ptr<CGSolver> cg (new CGSolver(N));
    cg->set_max_iter(cg_maxiter);
    cg->set_convergence(cg_conv);

    // evaluate guess energy (c.x):
    double energy_primal = C_DDOT(dimx,c->pointer(),1,x->pointer(),1);

    outfile->Printf("\n");
    outfile->Printf("    reference energy:     %20.12lf\n",escf);
    outfile->Printf("    frozen core energy:   %20.12lf\n",efrzc);
    outfile->Printf("    initial 2-RDM energy: %20.12lf\n",energy_primal + enuc + efrzc);
    outfile->Printf("\n");
    outfile->Printf("      oiter");
    outfile->Printf(" iiter");
    outfile->Printf("        E(p)");
    outfile->Printf("        E(d)");
    outfile->Printf("      E gap)");
    outfile->Printf("      mu");
    outfile->Printf("     eps(p)");
    outfile->Printf("     eps(d)\n");
    //outfile->Printf("    t(CG)");
    //outfile->Printf("  t(Diag)\n");

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
        //cg->preconditioned_solve(N,Ax,y,B,precon,evaluate_Ap,(void*)this);
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
            RotateOrbitals();
        }

        // compute current primal and dual energies

        //energy_primal = C_DDOT(dimx,c->pointer(),1,x->pointer(),1);
        double current_energy = C_DDOT(dimx,c->pointer(),1,x->pointer(),1);

        energy_dual   = C_DDOT(nconstraints,b->pointer(),1,y->pointer(),1);

        outfile->Printf("      %5i %5i %11.6lf %11.6lf %11.6lf %7.3lf %10.5lf %10.5lf\n",
                    oiter,iiter,current_energy+enuc+efrzc,energy_dual+efrzc+enuc,fabs(current_energy-energy_dual),mu,ep,ed);
        d2timeAu=q2timeAu=g2timeAu=d2timeATu=q2timeATu=g2timeATu=0.0;
        oiter++;
    
        if (oiter == maxiter) break;

        egap = fabs(current_energy-energy_dual);
        denergy_primal = fabs(energy_primal - current_energy);
        energy_primal = current_energy;

        if ( ep < r_conv && ed < r_conv && egap < e_conv ) {
            RotateOrbitals();
            energy_primal = C_DDOT(dimx,c->pointer(),1,x->pointer(),1);
        }

    }while( ep > r_conv || ed > r_conv  || egap > e_conv || !jacobi_converged_);

    if ( oiter == maxiter ) {
        throw PsiException("v2RDM did not converge.",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("      v2RDM iterations converged!\n");
    outfile->Printf("\n");


    // evaluate spin squared
    double s2 = 0.0;
    double * x_p = x->pointer();
    for (int i = 0; i < nmo; i++){
        for (int j = 0; j < nmo; j++){
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            int ij = ibas_ab_sym[h][i][j];
            int ji = ibas_ab_sym[h][j][i];
            s2 += x_p[d2aboff[h] + ij*gems_ab[h]+ji];
        }
    }
    int na = nalpha_ - nfrzc;
    int nb = nbeta_ - nfrzc;
    int ms = (multiplicity_ - 1)/2;
    outfile->Printf("      v2RDM total spin [S(S+1)]: %20.6lf\n", 0.5 * (na + nb) + ms*ms - s2);
    outfile->Printf("    * v2RDM total energy:        %20.12lf\n",energy_primal+enuc+efrzc);

    double en1 = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                int id = d1aoff[h] + i*amopi_[h] + j;
                en1 += x->pointer()[id] * c->pointer()[id];
                id = d1boff[h] + i*amopi_[h] + j;
                en1 += x->pointer()[id] * c->pointer()[id];
            }
        }
    }
    double enaa = 0.0;
    double enbb = 0.0;
    double enab = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int id = d2aboff[h] + ij*gems_ab[h] + kl;
                enab += x->pointer()[id] * c->pointer()[id];
            }
        }
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int id = d2aaoff[h] + ij*gems_aa[h] + kl;
                enaa += x->pointer()[id] * c->pointer()[id];
                id = d2bboff[h] + ij*gems_aa[h] + kl;
                enbb += x->pointer()[id] * c->pointer()[id];
            }
        }
    }
    //printf("correct answer: %20.12lf\n",energy_primal+efrzc);
    //printf("correct en1:    %20.12lf\n",en1);
    //printf("correct enaa:   %20.12lf\n",enaa);
    //printf("correct enbb:   %20.12lf\n",enbb);
    //printf("correct enab:   %20.12lf\n",enab);
    //printf("correct en2:    %20.12lf\n",enaa+enbb+enab);
    outfile->Printf("\n");

    //UnpackFullDensity();
    //UnpackDensityPlusCore();


    // maximal spin constraint:
    /*for (int h = 0; h < nirrep_; h++) {
        for (int klg = 0; klg < gems_ab[h]; klg++) {

            int k = bas_ab_sym[h][klg][0];
            int l = bas_ab_sym[h][klg][1];

            double maxs = 0.0;

            for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                int i = bas_ab_sym[h][ijg][0];
                int j = bas_ab_sym[h][ijg][1];
                if ( i != j ) continue;


                maxs += x_p[g2toff_p1[h] + ijg*gems_ab[h] + klg];

            }
            printf("%5i %5i %20.12lf \n",k,l,maxs);
        }
    }*/

    // compute and print natural orbital occupation numbers
    //Ca_->print();
    FinalTransformationMatrix();
    MullikenPopulations();
    NaturalOrbitals();

    return energy_primal + enuc + efrzc;
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


void  v2RDMSolver::common_init(){


    is_df_ = false;
    if ( options_.get_str("SCF_TYPE") == "DF" || options_.get_str("SCF_TYPE") == "CD" ) {
        is_df_ = true;
    }

    // This function is a copy of common_init() from SDPSolver
    escf      = reference_wavefunction_->reference_energy();
    nalpha_   = reference_wavefunction_->nalpha();
    nbeta_    = reference_wavefunction_->nbeta();
    nalphapi_ = reference_wavefunction_->nalphapi();
    nbetapi_  = reference_wavefunction_->nbetapi();
    doccpi_   = reference_wavefunction_->doccpi();
    soccpi_   = reference_wavefunction_->soccpi();
    frzcpi_   = reference_wavefunction_->frzcpi();
    nmopi_    = reference_wavefunction_->nmopi();
    nirrep_   = reference_wavefunction_->nirrep();
    nso_      = reference_wavefunction_->nso();
    nsopi_    = reference_wavefunction_->nsopi();

    // multiplicity:
    multiplicity_ = Process::environment.molecule()->multiplicity();
    if ( options_["MULTIPLICITY"].has_changed() ) {
        multiplicity_ = options_.get_int("MULTIPLICITY");
        int ms = (multiplicity_ - 1)/2;
        if ( nalpha_ != nbeta_ ) {
            throw PsiException("error: multiplicity keyword only works with na=nb",__FILE__,__LINE__);
        }
        nalpha_ += ms;
        nbeta_  -= ms;
    }

    if (options_["FROZEN_DOCC"].has_changed()) {
        if (options_["FROZEN_DOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_DOCC array has the wrong dimensions",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzcpi_[h] = options_["FROZEN_DOCC"][h].to_double();
        }
    }
    if (options_["FROZEN_UOCC"].has_changed()) {
        if (options_["FROZEN_UOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_UOCC array has the wrong dimensions",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzvpi_[h] = options_["FROZEN_UOCC"][h].to_double();
        }
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

    epsilon_a_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());
    
    nso = nmo = ndocc = nvirt = nfrzc = nfrzv = 0;
    for (int h = 0; h < nirrep_; h++){
        nfrzc    += frzcpi_[h];
        nfrzv    += frzvpi_[h];
        nso      += nsopi_[h];
        nmo      += nmopi_[h]-frzcpi_[h]-frzvpi_[h];
        ndocc    += doccpi_[h];
        amopi_[h] = nmopi_[h]-frzcpi_[h]-frzvpi_[h];
    }
    ndoccact = ndocc - nfrzc;
    nvirt    = nmo - ndoccact;
    //if (nfrzv > 0) {
    //    throw PsiException("bpsdp does not yet work with frozen virtuals",__FILE__,__LINE__);
    //}
    //if (nfrzc > 0) {
    //    throw PsiException("bpsdp does not yet work with frozen core",__FILE__,__LINE__);
    //}
    if (nmo + nfrzc + nfrzv != nso) {
        throw PsiException("bpsdp does not yet work when nmo!=nso",__FILE__,__LINE__);
    }
    // memory is from process::environment        
    memory_ = Process::environment.get_memory();
    // set the wavefunction name
    name_ = "BPSDP";

    // pick conditions.  default is dqg
    constrain_q2 = true;
    constrain_g2 = true;
    constrain_t1 = false;
    constrain_t2 = false;
    constrain_d3 = false;
    if (options_.get_str("POSITIVITY")=="D") {
        constrain_q2 = false;
        constrain_g2 = false;
    }else if (options_.get_str("POSITIVITY")=="DQ") {
        constrain_q2 = true;
        constrain_g2 = false;
    }else if (options_.get_str("POSITIVITY")=="DG") {
        constrain_q2 = false;
        constrain_g2 = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT1") {
        constrain_q2 = true;
        constrain_g2 = true;
        constrain_t1 = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT2") {
        constrain_q2 = true;
        constrain_g2 = true;
        constrain_t2 = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT1T2") {
        constrain_q2 = true;
        constrain_g2 = true;
        constrain_t1 = true;
        constrain_t2 = true;
    }else if (options_.get_str("POSITIVITY")=="DQGT") {
        constrain_q2 = true;
        constrain_g2 = true;
        constrain_t1 = true;
        constrain_t2 = true;
    }

    if ( options_.get_bool("CONSTRAIN_D3") ) {
        constrain_d3 = true;
    }

    spin_adapt_g2  = options_.get_bool("SPIN_ADAPT_G2");
    spin_adapt_q2  = options_.get_bool("SPIN_ADAPT_Q2");
    constrain_spin = options_.get_bool("CONSTRAIN_SPIN");

    if ( constrain_t1 || constrain_t2 ) {
        if (spin_adapt_g2) {
            throw PsiException("If constraining T1/T2, G2 cannot currently be spin adapted.",__FILE__,__LINE__);
        }
        if (spin_adapt_q2) {
            throw PsiException("If constraining T1/T2, Q2 cannot currently be spin adapted.",__FILE__,__LINE__);
        }
    }

    // build mapping arrays and determine the number of geminals per block
    BuildBasis();

    int ms = (multiplicity_ - 1)/2;
    if ( ms > 0 ) {
        if (spin_adapt_g2) {
            throw PsiException("G2 not spin adapted for S = M != 0",__FILE__,__LINE__);
        }
        if (spin_adapt_q2) {
            throw PsiException("Q2 not spin adapted for S = M != 0",__FILE__,__LINE__);
        }
    }

    // dimension of variable buffer (x) 
    dimx = 0;
    for ( int h = 0; h < nirrep_; h++) {
        dimx += gems_ab[h]*gems_ab[h]; // D2ab
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx += gems_aa[h]*gems_aa[h]; // D2aa
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx += gems_aa[h]*gems_aa[h]; // D2bb
    }
    for ( int h = 0; h < nirrep_; h++) {
        dimx += amopi_[h]*amopi_[h]; // D1a
        dimx += amopi_[h]*amopi_[h]; // D1b
        dimx += amopi_[h]*amopi_[h]; // Q1b
        dimx += amopi_[h]*amopi_[h]; // Q1a
    }
    if ( constrain_q2 ) {
        if ( !spin_adapt_q2 ) {
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_ab[h]*gems_ab[h]; // Q2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_aa[h]*gems_aa[h]; // Q2aa
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_aa[h]*gems_aa[h]; // Q2bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_00[h]*gems_00[h]; // Q2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_aa[h]*gems_aa[h]; // Q2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_aa[h]*gems_aa[h]; // Q2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_aa[h]*gems_aa[h]; // Q2t_m1
            }
        }
    }
    if ( constrain_g2 ) {
        if ( !spin_adapt_g2 ) {
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_ab[h]*gems_ab[h]; // G2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_ab[h]*gems_ab[h]; // G2ba
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += 2*gems_ab[h]*2*gems_ab[h]; // G2aa/bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_ab[h]*gems_ab[h]; // G2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_ab[h]*gems_ab[h]; // G2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_ab[h]*gems_ab[h]; // G2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                dimx += gems_ab[h]*gems_ab[h]; // G2t_m1
            }
        }
    }
    if ( constrain_t1 ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aaa[h]*trip_aaa[h]; // T1aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aaa[h]*trip_aaa[h]; // T1bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aab[h]*trip_aab[h]; // T1aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aab[h]*trip_aab[h]; // T1bba
        }
    }
    if ( constrain_t2 ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx += (trip_aba[h]+trip_aab[h])*(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += (trip_aba[h]+trip_aab[h])*(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aab[h]*trip_aab[h]; // T2aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aab[h]*trip_aab[h]; // T2bba
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aba[h]*trip_aba[h]; // T2aba
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aba[h]*trip_aba[h]; // T2bab
        }
    }
    if ( constrain_d3 ) {
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aaa[h] * trip_aaa[h]; // D3aaa
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aaa[h] * trip_aaa[h]; // D3bbb
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aab[h]*trip_aab[h]; // D3aab
        }
        for ( int h = 0; h < nirrep_; h++) {
            dimx += trip_aab[h]*trip_aab[h]; // D3bba
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

    if ( constrain_q2 ) {
        if ( !spin_adapt_q2 ) {
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

    if ( constrain_g2 ) {
        if ( ! spin_adapt_g2 ) {
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

    if ( constrain_t1 ) {
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

    if ( constrain_t2 ) {
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
    if ( constrain_d3 ) {
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
    nconstraints = 0;

    if ( constrain_spin ) {
        nconstraints += 1;               // spin
    }
    nconstraints += 1;                   // Tr(D2ab)
    nconstraints += 1;                   // Tr(D2aa)
    nconstraints += 1;                   // Tr(D2bb)
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints += amopi_[h]*amopi_[h]; // D1a <-> Q1a
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints += amopi_[h]*amopi_[h]; // D1b <-> Q1b
    }
    //for ( int h = 0; h < nirrep_; h++) {
    //    nconstraints += amopi_[h]*(amopi_[h]+1)/2; // contract D2ab + D2aa -> D1 a
    //}
    //for ( int h = 0; h < nirrep_; h++) {
    //    nconstraints += amopi_[h]*(amopi_[h]+1)/2; // contract D2ab + D2bb -> D1 b
    //}
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints += amopi_[h]*amopi_[h]; // contract D2ab        -> D1 a
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints += amopi_[h]*amopi_[h]; // contract D2ab        -> D1 b
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints += amopi_[h]*amopi_[h]; // contract D2aa        -> D1 a
    }
    for ( int h = 0; h < nirrep_; h++) {
        nconstraints += amopi_[h]*amopi_[h]; // contract D2bb        -> D1 b
    }
    if ( constrain_q2 ) {
        if ( ! spin_adapt_q2 ) {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_ab[h]*gems_ab[h]; // Q2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_aa[h]*gems_aa[h]; // Q2aa
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_aa[h]*gems_aa[h]; // Q2bb
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_00[h]*gems_00[h]; // Q2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_aa[h]*gems_aa[h]; // Q2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_aa[h]*gems_aa[h]; // Q2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_aa[h]*gems_aa[h]; // Q2t_m1
            }
        }
        
    }
    if ( constrain_g2 ) {
        if ( ! spin_adapt_g2 ) {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_ab[h]*gems_ab[h]; // G2ab
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_ab[h]*gems_ab[h]; // G2ba
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += 2*gems_ab[h]*2*gems_ab[h]; // G2aa
            }
        }else {
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_ab[h]*gems_ab[h]; // G2s
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_ab[h]*gems_ab[h]; // G2t
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_ab[h]*gems_ab[h]; // G2t_p1
            }
            for ( int h = 0; h < nirrep_; h++) {
                nconstraints += gems_ab[h]*gems_ab[h]; // G2t_m1
            }
            // maximal spin constraint
            //nconstraints += 2*nmo*nmo;
        }
    }
    if ( constrain_t1 ) {
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += trip_aaa[h]*trip_aaa[h]; // T1aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += trip_aaa[h]*trip_aaa[h]; // T1bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += trip_aab[h]*trip_aab[h]; // T1aab
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += trip_aab[h]*trip_aab[h]; // T1bba
        }
    }
    if ( constrain_t2 ) {
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += (trip_aab[h]+trip_aba[h])*(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += trip_aab[h]*trip_aab[h]; // T2aab
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += trip_aab[h]*trip_aab[h]; // T2bba
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += trip_aba[h]*trip_aba[h]; // T2aba
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += trip_aba[h]*trip_aba[h]; // T2bab
        }
    }
    if ( constrain_d3 ) {
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += gems_aa[h]*gems_aa[h]; // D3aaa -> D2aa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += gems_aa[h]*gems_aa[h]; // D3bbb -> D2bb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += gems_aa[h]*gems_aa[h]; // D3aab -> D2aa
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += gems_aa[h]*gems_aa[h]; // D3bba -> D2bb
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += gems_ab[h]*gems_ab[h]; // D3aab -> D2ab
        }
        for (int h = 0; h < nirrep_; h++) {
            nconstraints += gems_ab[h]*gems_ab[h]; // D3bba -> D2ab
        }
    }

    // list of dimensions
    for (int h = 0; h < nirrep_; h++) {
        dimensions.push_back(gems_ab[h]); // D2ab
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions.push_back(gems_aa[h]); // D2aa
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions.push_back(gems_aa[h]); // D2bb
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions.push_back(amopi_[h]); // D1a
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions.push_back(amopi_[h]); // D1b
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions.push_back(amopi_[h]); // Q1a
    }
    for (int h = 0; h < nirrep_; h++) {
        dimensions.push_back(amopi_[h]); // Q1b
    }
    if ( constrain_q2 ) {
        if ( !spin_adapt_q2 ) {
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_ab[h]); // Q2ab
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_aa[h]); // Q2aa
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_aa[h]); // Q2bb
            }
        }else {
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_00[h]); // Q2s
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_aa[h]); // Q2t
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_aa[h]); // Q2t_p1
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_aa[h]); // Q2t_m1
            }
        }
    }
    if ( constrain_g2 ) {
        if ( !spin_adapt_g2 ) {
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_ab[h]); // G2ab
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_ab[h]); // G2ba
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(2*gems_ab[h]); // G2aa
            }
        }else {
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_ab[h]); // G2s
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_ab[h]); // G2t
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_ab[h]); // G2t_p1
            }
            for (int h = 0; h < nirrep_; h++) {
                dimensions.push_back(gems_ab[h]); // G2t_m1
            }
        }
    }
    if ( constrain_t1 ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aaa[h]); // T1aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aaa[h]); // T1bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aab[h]); // T1aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aab[h]); // T1bba
        }
    }
    if ( constrain_t2 ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aab[h]+trip_aba[h]); // T2aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aab[h]+trip_aba[h]); // T2bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aab[h]); // T2aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aab[h]); // T2bba
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aba[h]); // T2aba
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aba[h]); // T2bab
        }
    }
    if ( constrain_d3 ) {
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aaa[h]); // D3aaa
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aaa[h]); // D3bbb
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aab[h]); // D3aab
        }
        for (int h = 0; h < nirrep_; h++) {
            dimensions.push_back(trip_aab[h]); // D3bba
        }
    }

    // allocate vectors
    Ax     = SharedVector(new Vector("A . x",nconstraints));
    ATy    = SharedVector(new Vector("A^T . y",dimx));
    x      = SharedVector(new Vector("primal solution",dimx));
    c      = SharedVector(new Vector("OEI and TEI",dimx));
    y      = SharedVector(new Vector("dual solution",nconstraints));
    z      = SharedVector(new Vector("dual solution 2",dimx));
    b      = SharedVector(new Vector("constraints",nconstraints));

    // bpsdp convergence thresholds:
    r_conv  = options_.get_double("R_CONVERGENCE");
    e_conv  = options_.get_double("E_CONVERGENCE");
    maxiter = options_.get_int("MAXITER");

    // conjugate gradient solver thresholds:
    cg_conv    = options_.get_double("CG_CONVERGENCE");
    cg_maxiter = options_.get_double("CG_MAXITER");

    // input/output array for jacobi sweeps
    jacobi_data_    = (double*)malloc(11*sizeof(double));
    jacobi_data_[0] = (double)options_.get_int("JACOBI_NTHREAD");
    jacobi_data_[1] = (double)options_.get_bool("JACOBI_ACTIVE_ACTIVE_ROTATIONS");
    jacobi_data_[2] = (double)options_.get_int("JACOBI_FROZEN_CORE");
    jacobi_data_[3] = (double)options_.get_double("JACOBI_ANGLE_TOLERANCE");
    jacobi_data_[4] = (double)options_.get_double("JACOBI_E_CONVERGENCE");
    jacobi_data_[5] = (double)options_.get_bool("JACOBI_WRITE");
    jacobi_data_[6] = 0.0;  // ntsweep: total sweeps (output)
    jacobi_data_[7] = 0.0;  // ntrot: total pairs rotated (output)
    jacobi_data_[8] = 0.0;  // delrot: change in energy
    jacobi_data_[9] = 0.0;  // converged?
    jacobi_data_[10] = 0.0;
    if ( is_df_ ) {
      jacobi_data_[10] = 1.0;
    }
    jacobi_converged_ = false;

    int full = nmo + nfrzc + nfrzv;
    jacobi_transformation_matrix_ = (double*)malloc(full*full*sizeof(double));
    memset((void*)jacobi_transformation_matrix_,'\0',full*full*sizeof(double));
    for (int i = 0; i < full; i++) {
        jacobi_transformation_matrix_[i*full+i] = 1.0;
    }

    // don't change the length of this filename
    jacobi_outfile_ = (char*)malloc(120*sizeof(char));
    std::string filename = get_writer_file_prefix() + ".jacobi";
    strcpy(jacobi_outfile_,filename.c_str());
    if ( options_.get_bool("JACOBI_WRITE") ) { 
        FILE * fp = fopen(jacobi_outfile_,"w");
        fclose(fp);
    }
}

boost::shared_ptr<Matrix> v2RDMSolver::GetOEI() {
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    boost::shared_ptr<Matrix> K1 (new Matrix(mints->so_potential()));
    K1->add(mints->so_kinetic());
    K1->transform(Ca_);
    return K1;
}

boost::shared_ptr<Matrix> v2RDMSolver::ReadOEI() {

    boost::shared_ptr<Matrix> K1 (new Matrix("one-electron integrals",nirrep_,nmopi_,nmopi_));

    long int full = nmo + nfrzc + nfrzv;
    double * tempoei = (double*)malloc(full*(full+1)/2*sizeof(double));
    memset((void*)tempoei,'\0',full*(full+1)/2*sizeof(double));
    boost::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_OEI,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_OEI,"MO-basis One-electron Ints",(char*)&tempoei[0],full*(full+1)/2*sizeof(double));
    psio->close(PSIF_OEI,1);
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h]; i++) {
            for (long int j = 0; j < nmopi_[h]; j++) {
                K1->pointer(h)[i][j] = tempoei[INDEX(i+offset,j+offset)];
            }
        }
        offset += nmopi_[h];
    }
    free(tempoei);

    return K1;
}

//build K2 using 3-index integrals 
void v2RDMSolver::DFK2() {

    // one-electron integrals:  
    boost::shared_ptr<Matrix> K1 = GetOEI();

    outfile->Printf("    ==> Transform three-electron integrals <==\n");
    outfile->Printf("\n");

    double start = omp_get_wtime();
    ThreeIndexIntegrals();

    double end = omp_get_wtime();

    outfile->Printf("\n");
    outfile->Printf("    Time for integral transformation:  %7.2lf s\n",end-start);
    outfile->Printf("\n");

    // size of the 3-index integral buffer
    tei_full_dim = nQ_*nso_*nso_;

    d2_plus_core_dim = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim += gems_plus_core[h] * ( gems_plus_core[h] + 1 ) / 2;
    }

    // just point to 3-index integral buffer
    tei_full_sym      = Qmo_;

    d2_plus_core_sym  = (double*)malloc(d2_plus_core_dim*sizeof(double));
    memset((void*)d2_plus_core_sym,'\0',d2_plus_core_dim*sizeof(double));

    // allocate memory for full OEI tensor blocked by symmetry for Greg
    oei_full_dim = 0;
    for (int h = 0; h < nirrep_; h++) {
        oei_full_dim += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
    d1_plus_core_dim = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_plus_core_dim += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }
    oei_full_sym = (double*)malloc(oei_full_dim*sizeof(double));
    d1_plus_core_sym  = (double*)malloc(d1_plus_core_dim*sizeof(double));
    memset((void*)oei_full_sym,'\0',oei_full_dim*sizeof(double));
    memset((void*)d1_plus_core_sym,'\0',d1_plus_core_dim*sizeof(double));

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h]; i++) {
            for (long int j = i; j < nmopi_[h]; j++) {
                oei_full_sym[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }

    RepackIntegralsDF();

    enuc = Process::environment.molecule()->nuclear_repulsion_energy();

}

//build K2 using 4-index integrals 
void v2RDMSolver::K2() {

    long int full = nmo + nfrzc + nfrzv;
    double * temptei = (double*)malloc(full*full*full*full*sizeof(double));
    memset((void*)temptei,'\0',full*full*full*full*sizeof(double));
      
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
    outfile->Printf("    Time for integral transformation:  %7.2lf s\n",end-start);
    outfile->Printf("\n");
    
    // read two-electron integrals from disk
    ReadIntegrals(temptei,full);
    
    // read oei's from disk (this will simplify things when we add symmetry):
    //boost::shared_ptr<Matrix> K1 = ReadOEI();
    boost::shared_ptr<Matrix> K1 = GetOEI();

    // allocate memory for full ERI tensor blocked by symmetry for Greg
    tei_full_dim = 0;
    for (int h = 0; h < nirrep_; h++) {
        tei_full_dim += gems_full[h] * ( gems_full[h] + 1 ) / 2;
    }
    d2_plus_core_dim = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim += gems_plus_core[h] * ( gems_plus_core[h] + 1 ) / 2;
    }
    
    tei_full_sym      = (double*)malloc(tei_full_dim*sizeof(double));
    d2_plus_core_sym  = (double*)malloc(d2_plus_core_dim*sizeof(double));
    memset((void*)tei_full_sym,'\0',tei_full_dim*sizeof(double));
    memset((void*)d2_plus_core_sym,'\0',d2_plus_core_dim*sizeof(double));
    //for (int i = 0; i < d2_plus_core_dim; i++) {
    //    d2_plus_core_sym[i] = -9.0e99;
    //}
    // load tei_full_sym
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int ij = 0; ij < gems_full[h]; ij++) {
            long int i = bas_full_sym[h][ij][0];
            long int j = bas_full_sym[h][ij][1];
            for (long int kl = ij; kl < gems_full[h]; kl++) {
                long int k = bas_full_sym[h][kl][0];
                long int l = bas_full_sym[h][kl][1];
                tei_full_sym[offset + INDEX(ij,kl)] = temptei[i*full*full*full+j*full*full+k*full+l];
            }
        }
        offset += gems_full[h] * ( gems_full[h] + 1 ) / 2;
    }
    // allocate memory for full OEI tensor blocked by symmetry for Greg
    oei_full_dim = 0;
    for (int h = 0; h < nirrep_; h++) {
        oei_full_dim += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
    d1_plus_core_dim = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_plus_core_dim += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }
    oei_full_sym = (double*)malloc(oei_full_dim*sizeof(double));
    d1_plus_core_sym  = (double*)malloc(d1_plus_core_dim*sizeof(double));
    memset((void*)oei_full_sym,'\0',oei_full_dim*sizeof(double));
    memset((void*)d1_plus_core_sym,'\0',d1_plus_core_dim*sizeof(double));
    
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h]; i++) {
            for (long int j = i; j < nmopi_[h]; j++) {
                oei_full_sym[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
    


    // if frozen core, adjust oei's and compute frozen core energy:
    efrzc1 = 0.0;
    efrzc2 = 0.0;
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < frzcpi_[h]; i++) {
            efrzc1 += 2.0 * K1->pointer(h)[i][i];

	    long int offset2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
	      for (long int j = 0; j < frzcpi_[h2]; j++) {
                    efrzc2 += (2.0 * temptei[(i+offset)*full*full*full+(i+offset)*full*full+(j+offset2)*full+(j+offset2)] 
                                  - temptei[(i+offset)*full*full*full+(j+offset2)*full*full+(i+offset)*full+(j+offset2)]);
                }
                offset2 += nmopi_[h2];
            }
        }
        offset += nmopi_[h];
    }
    efrzc = efrzc1 + efrzc2;
    printf("efrzc1       %20.12lf\n",efrzc1);
    printf("efrzc2       %20.12lf\n",efrzc2);
    
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = frzcpi_[h]; i < nmopi_[h] - frzvpi_[h]; i++) {
            for (long int j = frzcpi_[h]; j < nmopi_[h] - frzvpi_[h]; j++) {
	      
                double dum = 0.0;

		        int offset2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
		 
                    for (long int k = 0; k < frzcpi_[h2]; k++) {
		            dum += (2.0 * temptei[(i+offset)*full*full*full+(j+offset)*full*full+(k+offset2)*full+(k+offset2)] 
			            - temptei[(i+offset)*full*full*full+(k+offset2)*full*full+(k+offset2)*full+(j+offset)]);
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

    long int na = nalpha_ - nfrzc;
    long int nb = nbeta_ - nfrzc;
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

                c_p[d2aboff[h] + ij*gems_ab[h]+kl]    = temptei[ii*full*full*full+kk*full*full+jj*full+ll];

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

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]    = 0.5 * temptei[ii*full*full*full+kk*full*full+jj*full+ll];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * temptei[ii*full*full*full+ll*full*full+jj*full+kk];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * temptei[jj*full*full*full+kk*full*full+ii*full+ll];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   += 0.5 * temptei[jj*full*full*full+ll*full*full+ii*full+kk];

                long int hi = symmetry[i];
                long int hj = symmetry[j];

                c_p[d2bboff[h] + ij*gems_aa[h]+kl]    = 0.5 * temptei[ii*full*full*full+kk*full*full+jj*full+ll];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * temptei[ii*full*full*full+ll*full*full+jj*full+kk];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * temptei[jj*full*full*full+kk*full*full+ii*full+ll];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   += 0.5 * temptei[jj*full*full*full+ll*full*full+ii*full+kk];

            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (long int i = frzcpi_[h]; i < nmopi_[h] - frzvpi_[h]; i++) {
            for (long int j = frzcpi_[h]; j < nmopi_[h] - frzvpi_[h]; j++) {
                //c_p[d1aoff[h] + INDEX(i-frzcpi_[h],j-frzcpi_[h])] = K1->pointer(h)[i][j];
                //c_p[d1boff[h] + INDEX(i-frzcpi_[h],j-frzcpi_[h])] = K1->pointer(h)[i][j];
                //c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h]+(j-frzcpi_[h])] = 0.0;//K1->pointer(h)[i][j];
                //c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h]+(j-frzcpi_[h])] = 0.0;//K1->pointer(h)[i][j];
            }
        }
    }

    enuc = Process::environment.molecule()->nuclear_repulsion_energy();
    free(temptei);
}


void v2RDMSolver::Guess(){

    srand(0);//time(NULL));

    double* x_p = x->pointer();
    double* z_p = z->pointer();

    x->zero();
    z->zero();

    //for (int i = 0; i < dimx; i++) {
    //    x->pointer()[i] = 0.1;
    //    z->pointer()[i] = 0.1;
    //}
    //return;

    // Hartree-Fock guess for D2, D1, Q1, Q2, and G2

    // D2ab
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            if ( i < nalpha_ - nfrzc && j < nbeta_ - nfrzc ) {
                x_p[d2aboff[h] + ij*gems_ab[h]+ij] = 1.0;

            }
        }
    }
  
    // D2aa/D2bb
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            if ( i < nalpha_-nfrzc && j < nalpha_-nfrzc ) {
                x_p[d2aaoff[h] + ij*gems_aa[h]+ij] = 1.0;
            }
            if ( i < nbeta_-nfrzc && j < nbeta_-nfrzc ) {
                x_p[d2bboff[h] + ij*gems_aa[h]+ij] = 1.0;
            }
        }
    }

    // D1
    for (int h = 0; h < nirrep_; h++) {
        for (int i = frzcpi_[h]; i < doccpi_[h]+soccpi_[h]; i++) {
            int ii = i - frzcpi_[h];
            //x_p[d1aoff[h]+INDEX(ii,ii)] = 1.0;
            x_p[d1aoff[h]+ii*amopi_[h]+ii] = 1.0;
        }
        for (int i = frzcpi_[h]; i < doccpi_[h]; i++) {
            int ii = i - frzcpi_[h];
            //x_p[d1boff[h]+INDEX(ii,ii)] = 1.0;
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

    //z->copy(x.get());

    if ( constrain_q2 ) {
        if ( !spin_adapt_q2) {
            // map D2ab to Q2ab
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ij = 0; ij < gems_ab[h]; ij++) {
                    int i = bas_ab_sym[h][ij][0];
                    int j = bas_ab_sym[h][ij][1];
                    for (int kl = ij; kl < gems_ab[h]; kl++) {
                        int k = bas_ab_sym[h][kl][0];
                        int l = bas_ab_sym[h][kl][1];
                        double dum  =  x_p[d2aboff[h] + INDEX(kl,ij)];          // +D2(kl,ij)

                        if ( j==l ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        +=  x_p[q1aoff[h2] + INDEX(ii,kk)]; // +Q1(i,k) djl
                        }
                        if ( i==k ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            //dum        -=  x_p[d1boff[h2] + INDEX(ll,jj)]; // -D1(l,j) dik
                            dum        -=  x_p[d1boff[h2] + ll*amopi_[h2]+jj]; // -D1(l,j) dik
                        }

                        x_p[q2aboff[h] + INDEX(ij,kl)] = dum;              // -Q2(ij,kl)
                    }
                }
                offset += gems_ab[h]*(gems_ab[h]+1)/2;
            }
            // map D2aa to Q2aa
            for (int h = 0; h < nirrep_; h++) {
                for (int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for (int kl = ij; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum  =  x_p[d2aaoff[h] + INDEX(kl,ij)];    // +D2(kl,ij)

                        if ( j==l ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        +=  x_p[q1aoff[h2] + INDEX(ii,kk)];  // +Q1(i,k) djl
                        }
                        if ( j==k ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            //dum        +=  x_p[d1aoff[h2] + INDEX(ll,ii)];  // +D1(l,i) djk
                            dum        +=  x_p[d1aoff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                        }
                        if ( i==l ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        -=  x_p[q1aoff[h2] + INDEX(jj,kk)];  // -Q1(j,k) dil
                        }
                        if ( i==k ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        -=  x_p[d1aoff[h2] + INDEX(ll,jj)];  // -D1(l,j) dkl
                        }

                        x_p[q2aaoff[h] + INDEX(ij,kl)] = dum;    // -Q2(ij,kl)
                    }
                }
                offset += gems_aa[h]*(gems_aa[h]+1)/2;
            }
            // map D2bb to Q2bb
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for (int kl = ij; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum  =  x_p[d2bboff[h] + INDEX(kl,ij)];    // +D2(kl,ij)

                        if ( j==l ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        +=  x_p[q1boff[h2] + INDEX(ii,kk)];  // +Q1(i,k) djl
                        }
                        if ( j==k ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        +=  x_p[d1boff[h2] + INDEX(ll,ii)];  // +D1(l,i) djk
                        }
                        if ( i==l ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        -=  x_p[q1boff[h2] + INDEX(jj,kk)];  // -Q1(j,k) dil
                        }
                        if ( i==k ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        -=  x_p[d1boff[h2] + INDEX(ll,jj)];  // -D1(l,j) dkl
                        }

                        x_p[q2bboff[h] + INDEX(ij,kl)] = dum;    // -Q2(ij,kl)
                    }
                }
                offset += gems_aa[h]*(gems_aa[h]+1)/2;
            }
        }else {
            // map D2ab to Q2s
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ij = 0; ij < gems_00[h]; ij++) {
                    int i = bas_00_sym[h][ij][0];
                    int j = bas_00_sym[h][ij][1];
                    int ijd = ibas_ab_sym[h][i][j];
                    int jid = ibas_ab_sym[h][j][i];
                    for (int kl = ij; kl < gems_00[h]; kl++) {
                        int k = bas_00_sym[h][kl][0];
                        int l = bas_00_sym[h][kl][1];

                        double dum  = 0.0;

                        // not spin adapted
                        int kld = ibas_ab_sym[h][k][l];
                        int lkd = ibas_ab_sym[h][l][k];
                        dum        +=  0.5 * x_p[d2aboff[h] + INDEX(kld,ijd)];          // +D2(kl,ij)
                        dum        +=  0.5 * x_p[d2aboff[h] + INDEX(lkd,ijd)];          // +D2(lk,ij)
                        dum        +=  0.5 * x_p[d2aboff[h] + INDEX(kld,jid)];          // +D2(kl,ji)
                        dum        +=  0.5 * x_p[d2aboff[h] + INDEX(lkd,jid)];          // +D2(lk,ji)

                        // spin adapted
                        //dum        +=  x_p[d2soff[h] + INDEX(kl,ij)];          // +D2(kl,ij)

                        if ( j==l ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        +=  0.5 * x_p[q1aoff[h2] + INDEX(ii,kk)]; // +Q1(i,k) djl
                            dum        -=  0.5 * x_p[d1boff[h2] + INDEX(ii,kk)]; // -D1(i,k) djl
                        }
                        if ( i==k ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        +=  0.5 * x_p[q1aoff[h2] + INDEX(ll,jj)]; // +Q1(l,j) dik
                            dum        -=  0.5 * x_p[d1boff[h2] + INDEX(ll,jj)]; // -D1(l,j) dik
                        }
                        if ( j==k ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        +=  0.5 * x_p[q1aoff[h2] + INDEX(ll,ii)]; // +Q1(l,i) djk
                            dum        -=  0.5 * x_p[d1boff[h2] + INDEX(ll,ii)]; // -D1(l,i) djk
                        }
                        if ( i==l ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        +=  0.5 * x_p[q1aoff[h2] + INDEX(kk,jj)]; // +Q1(k,j) dil
                            dum        -=  0.5 * x_p[d1boff[h2] + INDEX(kk,jj)]; // -D1(k,j) dil
                        }

                        x_p[q2soff[h] + INDEX(ij,kl)] = dum;          // -Q2(ij,kl)

                    }
                }
            }
            // map D2ab to Q210
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ij = 0; ij < gems_aa[h]; ij++) {
                    int i   =  bas_aa_sym[h][ij][0];
                    int j   =  bas_aa_sym[h][ij][1];
                    int ijd = ibas_ab_sym[h][i][j];
                    int jid = ibas_ab_sym[h][j][i];
                    for (int kl = ij; kl < gems_aa[h]; kl++) {
                        int   k =  bas_aa_sym[h][kl][0];
                        int   l =  bas_aa_sym[h][kl][1];

                        double dum  = 0.0;

                        // not spin adapted
                        int kld = ibas_ab_sym[h][k][l];
                        int lkd = ibas_ab_sym[h][l][k];
                        dum        +=  0.5 * x_p[d2aboff[h] + INDEX(kld,ijd)];          // +D2(kl,ij)
                        dum        -=  0.5 * x_p[d2aboff[h] + INDEX(lkd,ijd)];          // -D2(lk,ij)
                        dum        -=  0.5 * x_p[d2aboff[h] + INDEX(kld,jid)];          // -D2(kl,ji)
                        dum        +=  0.5 * x_p[d2aboff[h] + INDEX(lkd,jid)];          // +D2(lk,ji)

                        // spin adapted
                        //dum        +=  x_p[d2toff[h] + INDEX(kl,ij)];          // +D2(kl,ij)

                        if ( j==l ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        +=  0.5 * x_p[q1aoff[h2] + INDEX(ii,kk)]; // +Q1(i,k) djl
                            dum        -=  0.5 * x_p[d1boff[h2] + INDEX(ii,kk)]; // -D1(i,k) djl
                        }
                        if ( i==k ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        +=  0.5 * x_p[q1aoff[h2] + INDEX(ll,jj)]; // +Q1(l,j) dik
                            dum        -=  0.5 * x_p[d1boff[h2] + INDEX(ll,jj)]; // -D1(l,j) dik
                        }
                        if ( j==k ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        -=  0.5 * x_p[q1aoff[h2] + INDEX(ll,ii)]; // -Q1(l,i) djk
                            dum        +=  0.5 * x_p[d1boff[h2] + INDEX(ll,ii)]; // +D1(l,i) djk
                        }
                        if ( i==l ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        -=  0.5 * x_p[q1aoff[h2] + INDEX(kk,jj)]; // -Q1(k,j) dil
                            dum        +=  0.5 * x_p[d1boff[h2] + INDEX(kk,jj)]; // +D1(k,j) dil
                        }

                        x_p[q2toff[h] + INDEX(ij,kl)] = dum;
                    }
                }
            }
            // map D2aa to Q211
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for (int kl = ij; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum  = 0.0;
                        //dum        +=  x_p[d2toff_p1[h] + INDEX(kl,ij)];    // +D2(kl,ij)
                        dum        +=  x_p[d2aaoff[h] + INDEX(kl,ij)];    // +D2(kl,ij)

                        if ( j==l ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        +=  x_p[q1aoff[h2] + INDEX(ii,kk)];  // +Q1(i,k) djl
                        }
                        if ( j==k ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        +=  x_p[d1aoff[h2] + INDEX(ll,ii)];  // +D1(l,i) djk
                        }
                        if ( i==l ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        -=  x_p[q1aoff[h2] + INDEX(jj,kk)];  // -Q1(j,k) dil
                        }
                        if ( i==k ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        -=  x_p[d1aoff[h2] + INDEX(ll,jj)];  // -D1(l,j) dkl
                        }

                        x_p[q2toff_p1[h] + INDEX(ij,kl)] = dum;
                    }
                }
            }
            // map D2bb to Q21-1
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ij = 0; ij < gems_aa[h]; ij++) {
                    int i = bas_aa_sym[h][ij][0];
                    int j = bas_aa_sym[h][ij][1];
                    for (int kl = ij; kl < gems_aa[h]; kl++) {
                        int k = bas_aa_sym[h][kl][0];
                        int l = bas_aa_sym[h][kl][1];
                        double dum  = 0.0;
                        //dum        +=  x_p[d2toff_m1[h] + INDEX(kl,ij)];    // +D2(kl,ij)
                        dum        +=  x_p[d2bboff[h] + INDEX(kl,ij)];    // +D2(kl,ij)

                        if ( j==l ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        +=  x_p[q1boff[h2] + INDEX(ii,kk)];  // +Q1(i,k) djl
                        }
                        if ( j==k ) {
                            int h2 = symmetry[i];
                            int ii = i - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        +=  x_p[d1boff[h2] + INDEX(ll,ii)];  // +D1(l,i) djk
                        }
                        if ( i==l ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int kk = k - pitzer_offset[h2];
                            dum        -=  x_p[q1boff[h2] + INDEX(jj,kk)];  // -Q1(j,k) dil
                        }
                        if ( i==k ) {
                            int h2 = symmetry[j];
                            int jj = j - pitzer_offset[h2];
                            int ll = l - pitzer_offset[h2];
                            dum        -=  x_p[d1boff[h2] + INDEX(ll,jj)];  // -D1(l,j) dkl
                        }

                        x_p[q2toff_m1[h] + INDEX(ij,kl)] = dum;
                    }
                }
            }
        }
    }

    if ( constrain_g2 ) {
        if ( ! spin_adapt_g2 ) {
            // G2ab constraints:
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = ijg; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];


                        double dum = 0.0;

                        if (j==l) {
                            int h3 = symmetry[i];
                            int ii = i - pitzer_offset[h3];
                            int kk = k - pitzer_offset[h3];
                            dum   +=  x_p[d1aoff[h3] + INDEX(ii,kk)];      //   D1(i,k) djl
                        }

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                        int ild = ibas_ab_sym[h2][i][l];
                        int kjd = ibas_ab_sym[h2][k][j];

                        dum       -=  x_p[d2aboff[h2] + INDEX(ild,kjd)];   // - D2ab(il,kj)

                        x_p[g2aboff[h] + INDEX(ijg,klg)] = dum;    // - G2ab(ij,kl)
                    }
                }
                offset += gems_ab[h]*(gems_ab[h]+1)/2;
            }
            // G2ba constraints:
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = ijg; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];

                        double dum = 0.0;

                        if (j==l) {
                            int h3 = symmetry[i];
                            int ii = i - pitzer_offset[h3];
                            int kk = k - pitzer_offset[h3];
                            dum       +=  x_p[d1boff[h3] + INDEX(ii,kk)];      //   D1(i,k) djl
                        }

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                        int lid = ibas_ab_sym[h2][l][i];
                        int jkd = ibas_ab_sym[h2][j][k];

                        dum       -=  x_p[d2aboff[h2] + INDEX(lid,jkd)];       //   D2ab(li,jk)

                        x_p[g2baoff[h] + INDEX(ijg,klg)] = dum;        // - G2ba(ij,kl)
                    }
                }
                offset += gems_ab[h]*(gems_ab[h]+1)/2;
            }
            // G2aaaa / G2aabb / G2bbaa / G2bbbb
            for (int h = 0; h < nirrep_; h++) {
                // G2aaaa
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = ijg; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                        double dum = 0.0;
                        if ( j == l ) {
                            int h3 = symmetry[i];
                            int ii = i - pitzer_offset[h3];
                            int kk = k - pitzer_offset[h3];
                            dum       +=  x_p[d1aoff[h3] + INDEX(ii,kk)];    //   D1(i,k) djl
                        }
                        
                        if ( i != l && k != j ) {

                            int sil = ( i < l ? 1 : -1 );
                            int skj = ( k < j ? 1 : -1 );

                            int ild = ibas_aa_sym[h2][i][l];
                            int kjd = ibas_aa_sym[h2][k][j];

                            dum       -=  x_p[d2aaoff[h2] + INDEX(ild,kjd)] * sil * skj; // -D2aa(il,kj)

                        }

                        x_p[g2aaoff[h] + INDEX(ijg,klg)] = dum;       // - G2aaaa(ij,kl)
                    }
                }
                // G2bbbb
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = ijg; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                        double dum = 0.0;
                        if ( j == l ) {
                            int h3 = symmetry[i];
                            int ii = i - pitzer_offset[h3];
                            int kk = k - pitzer_offset[h3];
                            dum       +=  x_p[d1boff[h3] + INDEX(ii,kk)];    //   D1(i,k) djl
                        }
                        
                        if ( i != l && k != j ) {

                            int sil = ( i < l ? 1 : -1 );
                            int skj = ( k < j ? 1 : -1 );

                            int ild = ibas_aa_sym[h2][i][l];
                            int kjd = ibas_aa_sym[h2][k][j];

                            dum       -=  x_p[d2bboff[h2] + INDEX(ild,kjd)] * sil * skj; // -D2bb(il,kj)

                        }

                        x_p[g2aaoff[h] + INDEX(gems_ab[h] + ijg,gems_ab[h] + klg)] = dum; // - G2bbbb(ij,kl)
                    }
                }
                // G2aabb
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = ijg; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                        double dum = 0.0;
                        
                        int ild = ibas_ab_sym[h2][i][l];
                        int jkd = ibas_ab_sym[h2][j][k];

                        dum       +=  x_p[d2aboff[h2] + INDEX(ild,jkd)]; // D2ab(il,jk)

                        x_p[g2aaoff[h] + INDEX(ijg,klg+gems_ab[h])] = dum;       // - G2aabb(ij,kl)
                    }
                }
                // G2bbaa
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = ijg; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);

                        double dum = 0.0;
                        
                        int lid = ibas_ab_sym[h2][l][i];
                        int kjd = ibas_ab_sym[h2][k][j];

                        dum       +=  x_p[d2aboff[h2] + INDEX(lid,kjd)]; // D2ab(li,kj)

                        x_p[g2aaoff[h] + INDEX(ijg+gems_ab[h],klg)] = dum;       // - G2bbaa(ij,kl)
                    }
                }
                offset += 2*gems_ab[h]*(2*gems_ab[h]+1)/2;
            }
        }else {
            // G200
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = 0; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];

                        double dum = 0.0;

                        if ( j == l ) {
                            int h3 = symmetry[i];
                            int ii = i - pitzer_offset[h3];
                            int kk = k - pitzer_offset[h3];
                            dum       +=  x_p[d1aoff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                            dum       +=  x_p[d1boff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                        }

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                        //int ils = ibas_00_sym[h2][i][l];
                        //int jks = ibas_00_sym[h2][j][k];

                        //dum       +=  x_p[d2soff[h2] + INDEX(ils,jks)] * 0.5; //   D2s(li,kj)

                        if ( i != l && k != j ) {

                            int sil = ( i < l ? 1 : -1 );
                            int skj = ( k < j ? 1 : -1 );


                            int ild = ibas_aa_sym[h2][i][l];
                            int kjd = ibas_aa_sym[h2][k][j];
                            dum       -=  x_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2aa(il,kj)
                            dum       -=  x_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2bb(il,kj)

                            //int ilt = ibas_aa_sym[h2][i][l];
                            //int kjt = ibas_aa_sym[h2][k][j];
                            //dum       -=  x_p[d2toff[h2]    + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D210(il,kj)
                            //dum       -=  x_p[d2toff_p1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D211(il,kj)
                            //dum       -=  x_p[d2toff_m1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D21-1(il,kj)

                        }

                        int ild = ibas_ab_sym[h2][i][l];
                        int jkd = ibas_ab_sym[h2][j][k];

                        dum       +=  x_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] * 0.5; // D2ab(il,jk)

                        int lid = ibas_ab_sym[h2][l][i];
                        int kjd = ibas_ab_sym[h2][k][j];

                        dum       +=  x_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] * 0.5; // D2ab(li,kj)

                        x_p[g2soff[h] + ijg*gems_ab[h]+klg] = dum;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
            // G210
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = 0; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];


                        double dum = 0.0;

                        if ( j == l ) {
                            int h3 = symmetry[i];
                            int ii = i - pitzer_offset[h3];
                            int kk = k - pitzer_offset[h3];
                            dum       +=  x_p[d1aoff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                            dum       +=  x_p[d1boff[h3] + ii*amopi_[h3]+kk] * 0.5; //   D1(i,k) djl
                        }

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                        //int ils = ibas_00_sym[h2][i][l];
                        //int jks = ibas_00_sym[h2][j][k];

                        //dum       -=  x_p[d2soff[h2] + INDEX(ils,jks)] * 0.5; //   D2s(li,kj)

                        if ( i != l && k != j ) {

                            int sil = ( i < l ? 1 : -1 );
                            int skj = ( k < j ? 1 : -1 );

                            int ild = ibas_aa_sym[h2][i][l];
                            int kjd = ibas_aa_sym[h2][k][j];
                            dum       -=  x_p[d2aaoff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2aa(il,kj)
                            dum       -=  x_p[d2bboff[h2] + ild*gems_aa[h2]+kjd] * sil * skj * 0.5; // -D2bb(il,kj)

                            //int ilt = ibas_aa_sym[h2][i][l];
                            //int kjt = ibas_aa_sym[h2][k][j];
                            //dum       +=  x_p[d2toff[h2]    + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D210(il,kj)
                            //dum       -=  x_p[d2toff_p1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D211(il,kj)
                            //dum       -=  x_p[d2toff_m1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D21-1(il,kj)

                        }

                        int ild = ibas_ab_sym[h2][i][l];
                        int jkd = ibas_ab_sym[h2][j][k];

                        dum       -=  x_p[d2aboff[h2] + ild*gems_ab[h2]+jkd] * 0.5; // D2ab(il,jk)

                        int lid = ibas_ab_sym[h2][l][i];
                        int kjd = ibas_ab_sym[h2][k][j];

                        dum       -=  x_p[d2aboff[h2] + lid*gems_ab[h2]+kjd] * 0.5; // D2ab(li,kj)

                        x_p[g2toff[h] + ijg*gems_ab[h]+klg] = dum;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
            // G211 constraints:
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = 0; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];


                        double dum = 0.0;

                        if ( j == l ) {
                            int h3 = symmetry[i];
                            int ii = i - pitzer_offset[h3];
                            int kk = k - pitzer_offset[h3];
                            dum   +=  x_p[d1aoff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                        }

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                        int ild = ibas_ab_sym[h2][l][i];
                        int kjd = ibas_ab_sym[h2][j][k];
                        dum       -=  x_p[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                        //int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                        //int ils = ibas_00_sym[h2][i][l];
                        //int kjs = ibas_00_sym[h2][k][j];
                        //dum       -=  x_p[d2soff[h2] + INDEX(ils,kjs)] * 0.5; //   D2s(li,kj)

                        //if ( i != l && k != j ) {

                        //    int sil = ( i < l ? 1 : -1 );
                        //    int skj = ( k < j ? 1 : -1 );

                        //    int ilt = ibas_aa_sym[h2][i][l];
                        //    int kjt = ibas_aa_sym[h2][k][j];

                        //    dum       -=  x_p[d2toff_p1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D211(il,kj)

                        //}

                        x_p[g2toff_p1[h] + ijg*gems_ab[h]+klg] = dum;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
            // G21-1 constraints:
            for (int h = 0; h < nirrep_; h++) {
                #pragma omp parallel for schedule (static)
                for (int ijg = 0; ijg < gems_ab[h]; ijg++) {

                    int i = bas_ab_sym[h][ijg][0];
                    int j = bas_ab_sym[h][ijg][1];

                    for (int klg = 0; klg < gems_ab[h]; klg++) {

                        int k = bas_ab_sym[h][klg][0];
                        int l = bas_ab_sym[h][klg][1];


                        double dum = 0.0;

                        if ( j == l ) {
                            int h3 = symmetry[i];
                            int ii = i - pitzer_offset[h3];
                            int kk = k - pitzer_offset[h3];
                            dum   +=  x_p[d1boff[h3] + ii*amopi_[h3]+kk];      //   D1(i,k) djl
                        }

                        int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                        int ild = ibas_ab_sym[h2][l][i];
                        int kjd = ibas_ab_sym[h2][j][k];
                        dum       -=  x_p[d2aboff[h2] + ild*gems_ab[h2]+kjd];   // - D2ab(il,kj)

                        //int h2 = SymmetryPair(symmetry[i],symmetry[l]);
                        //int ils = ibas_00_sym[h2][i][l];
                        //int kjs = ibas_00_sym[h2][k][j];
                        //dum       -=  x_p[d2soff[h2] + INDEX(ils,kjs)] * 0.5; //   D2s(li,kj)

                        //if ( i != l && k != j ) {

                        //    int sil = ( i < l ? 1 : -1 );
                        //    int skj = ( k < j ? 1 : -1 );

                        //    int ilt = ibas_aa_sym[h2][i][l];
                        //    int kjt = ibas_aa_sym[h2][k][j];

                        //    dum       -=  x_p[d2toff_m1[h2] + INDEX(ilt,kjt)] * sil * skj * 0.5; //   D211(il,kj)

                        //}

                        x_p[g2toff_m1[h] + ijg*gems_ab[h]+klg] = dum;
                    }
                }
                offset += gems_ab[h]*gems_ab[h];
            }
        }
    }

    if ( constrain_t1 ) {
        T1_constraints_guess(x);
    }
    if ( constrain_t2 ) {
        T2_constraints_guess(x);
    }

    // zero guess for z
    z->zero();
}

void v2RDMSolver::PrintHeader(){
       
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
    outfile->Printf("        Freeze core orbitals?               %5s\n",nfrzc > 0 ? "yes" : "no");
    outfile->Printf("        Number of frozen core orbitals:     %5i\n",nfrzc);
    outfile->Printf("        Number of active occupied orbitals: %5i\n",ndoccact);
    outfile->Printf("        Number of active virtual orbitals:  %5i\n",nvirt);
    outfile->Printf("        Number of frozen virtual orbitals:  %5i\n",nfrzv);
    outfile->Printf("        r_convergence:                  %5.3le\n",r_conv);
    outfile->Printf("        e_convergence:                  %5.3le\n",e_conv);
    outfile->Printf("        cg_convergence:                 %5.3le\n",cg_conv);
    outfile->Printf("        maxiter:                            %5i\n",maxiter);
    outfile->Printf("        cg_maxiter:                         %5i\n",cg_maxiter);
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
        if ( constrain_g2 ) {
            if ( 2*gems_ab[h] > maxgem ) {
                maxgem = 2*gems_ab[h];
            }
        }

        if ( constrain_t1 ) {
            nt1 += trip_aaa[h] * trip_aaa[h]; // T1aaa
            nt1 += trip_aaa[h] * trip_aaa[h]; // T1bbb
            nt1 += trip_aab[h] * trip_aab[h]; // T1aab
            nt1 += trip_aab[h] * trip_aab[h]; // T1bba
            if ( trip_aab[h] > maxgem ) {
                maxgem = trip_aab[h];
            }
        }

        if ( constrain_t2 ) {
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
    if ( constrain_q2 ) {
        outfile->Printf("        Q2:                       %7.2lf mb\n",nd2 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_g2 ) {
        outfile->Printf("        G2:                       %7.2lf mb\n",ng2 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_d3 ) {
        outfile->Printf("        D3:                       %7.2lf mb\n",nt1 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_t1 ) {
        outfile->Printf("        T1:                       %7.2lf mb\n",nt1 * 8.0 / 1024.0 / 1024.0);
    }
    if ( constrain_t2 ) {
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

    double tot = 4.0*dimx + 4.0*nconstraints + 2.0*maxgem*maxgem;
    tot += nd2; // for K2a, K2b

    // for casscf, need d2 and 3- or 4-index integrals

    // allocate memory for full ERI tensor blocked by symmetry for Greg
    tei_full_dim = 0;
    for (int h = 0; h < nirrep_; h++) {
        tei_full_dim += gems_full[h] * ( gems_full[h] + 1 ) / 2;
    }
    d2_plus_core_dim = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim += gems_plus_core[h] * ( gems_plus_core[h] + 1 ) / 2;
    }
    tot += d2_plus_core_dim;
    if ( is_df_ ) {
        nQ_ = Process::environment.globals["NAUX (SCF)"];
        // for 3-index integrals.  the factor of 2 comes from the AO/MO 
        // transformation.  TODO: transform these before allocating the
        // other memory so the factor of 2 doesn't affect the overall
        // memory footprint
        tot += nQ_*nso_*nso_*2;
    }else {
        tei_full_dim = 0;
        for (int h = 0; h < nirrep_; h++) {
            tei_full_dim += gems_full[h] * ( gems_full[h] + 1 ) / 2;
        }
        tot += tei_full_dim;
        tot += nso_*nso_*nso_*nso_; // for four-index integrals stored stupidly 
    }
    
    outfile->Printf("        Total number of variables:     %10i\n",dimx);
    outfile->Printf("        Total number of constraints:   %10i\n",nconstraints);
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

}

void v2RDMSolver::BuildConstraints(){

    //constraint on the Trace of D2(s=0,ms=0)

    int na = nalpha_ - nfrzc;
    int nb = nbeta_ - nfrzc;
    double trdab = na * nb;

    //constraint on the Trace of D2(s=1,ms=0)
    double trdaa  = na*(na-1.0);
    double trdbb  = nb*(nb-1.0);

    b->zero();
    double* b_p = b->pointer();

    offset = 0;

    // funny ab trace with spin: N/2 + Ms^2 - S(S+1)
    if ( constrain_spin ) {
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

    if ( constrain_q2 ) {
        if ( !spin_adapt_q2 ) {
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

    if ( constrain_g2 ) {
        if ( ! spin_adapt_g2 ) {
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
            // maximal spin constraint:
            //for(int i = 0; i < nmo; i++){
            //    for(int j = 0; j < nmo; j++){
            //        b_p[offset + i*nmo+j] = 0.0;
            //    }
            //}
            //offset += nmo*nmo;
            //for(int i = 0; i < nmo; i++){
            //    for(int j = 0; j < nmo; j++){
            //        b_p[offset + i*nmo+j] = 0.0;
            //    }
            //}
            //offset += nmo*nmo;
        }
    }

    if ( constrain_t1 ) {
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

    if ( constrain_t2 ) {
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
    if ( constrain_d3 ) {
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

// D2 portion of A.x (and D1/Q1)
void v2RDMSolver::D2_constraints_Au(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    if ( constrain_spin ) {
        // spin
        double s2 = 0.0;
        for (int i = 0; i < nmo; i++){
            for (int j = 0; j < nmo; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_ab[h] == 0 ) continue;
                int ij = ibas_ab_sym[h][i][j];
                int ji = ibas_ab_sym[h][j][i];
                s2 += u_p[d2aboff[h] + ij*gems_ab[h]+ji];
            }
        }
        A_p[offset] = s2;
        offset++;
    }
    if (prints2) {
        // spin
        double s2 = 0.0;
        for (int i = 0; i < nmo; i++){
            for (int j = 0; j < nmo; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_ab[h] == 0 ) continue;
                int ij = ibas_ab_sym[h][i][j];
                int ji = ibas_ab_sym[h][j][i];
                s2 += u_p[d2aboff[h] + ij*gems_ab[h]+ji];
            }
        }
        printf("%20.12lf\n",s2);
    }

    // Traces
    // Tr(D2ab)
    double sumab =0.0;
    for (int i = 0; i < nmo; i++){
        for (int j = 0; j < nmo; j++){
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            if ( gems_ab[h] == 0 ) continue;
            int ij = ibas_ab_sym[h][i][j];
            sumab += u_p[d2aboff[h] + ij*gems_ab[h]+ij];
        }
    }
    A_p[offset] = sumab;
    offset++;

    // Tr(D2aa)
    double sumaa =0.0;
    for (int i = 0; i < nmo; i++){
        for (int j = 0; j < nmo; j++){
            if ( i==j ) continue;
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            if ( gems_aa[h] == 0 ) continue;
            int ij = ibas_aa_sym[h][i][j];
            sumaa += u_p[d2aaoff[h] + ij*gems_aa[h]+ij];
        }

    }
    A_p[offset] = sumaa;
    offset++;

    // Tr(D2bb)
    double sumbb =0.0;
    for (int i = 0; i < nmo; i++){
        for (int j = 0; j < nmo; j++){
            if ( i==j ) continue;
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            if ( gems_aa[h] == 0 ) continue;
            int ij = ibas_aa_sym[h][i][j];
            sumbb += u_p[d2bboff[h] + ij*gems_aa[h]+ij];
        }

    }
    A_p[offset] = sumbb;
    offset++;

    // d1 / q1 a
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A_p[offset+i*amopi_[h]+j] = u_p[d1aoff[h]+j*amopi_[h]+i] + u_p[q1aoff[h]+i*amopi_[h]+j];
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    // d1 / q1 b
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A_p[offset+i*amopi_[h]+j] = u_p[d1boff[h]+j*amopi_[h]+i] + u_p[q1boff[h]+i*amopi_[h]+j];
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    int na = nalpha_ - nfrzc;
    int nb = nbeta_ - nfrzc;
    double n = na + nb;

    // contraction: D2aa + D2ab -> D1 a
    int poff = 0;

    // contraction: D2ab -> D1 a
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = nb * u_p[d1aoff[h] + i*amopi_[h]+j];
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < nmo; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][ii][k];
                    int jk = ibas_ab_sym[h2][jj][k];
                    sum -= u_p[d2aboff[h2] + ik*gems_ab[h2]+jk];
                }
                A_p[offset + i*amopi_[h]+j] = sum;
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }

    // contraction: D2ab -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = na * u_p[d1boff[h] + i*amopi_[h]+j];
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < nmo; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][k][ii];
                    int jk = ibas_ab_sym[h2][k][jj];
                    sum -= u_p[d2aboff[h2] + ik*gems_ab[h2]+jk];
                }
                A_p[offset + i*amopi_[h]+j] = sum;
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }

    //contract D2aa -> D1 a
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = (na - 1.0) * u_p[d1aoff[h] + i*amopi_[h]+j];
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < nmo; k++){
                    if( ii==k || jj==k ) continue;
                    int h2   = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik  = ibas_aa_sym[h2][ii][k];
                    int jk  = ibas_aa_sym[h2][jj][k];
                    int sik = ( ii < k ) ? 1 : -1;
                    int sjk = ( jj < k ) ? 1 : -1;
                    sum -= sik*sjk*u_p[d2aaoff[h2] + ik*gems_aa[h2]+jk];
                }
                A_p[offset+i*amopi_[h]+j] = sum;
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }

    //contract D2bb -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = (nb - 1.0) * u_p[d1boff[h] + i*amopi_[h]+j];
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < nmo; k++){
                    if( ii==k || jj==k ) continue;
                    int h2   = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik  = ibas_aa_sym[h2][ii][k];
                    int jk  = ibas_aa_sym[h2][jj][k];
                    int sik = ( ii < k ) ? 1 : -1;
                    int sjk = ( jj < k ) ? 1 : -1;
                    sum -= sik*sjk*u_p[d2bboff[h2] + ik*gems_aa[h2]+jk];
                }
                A_p[offset+i*amopi_[h]+j] = sum;
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }
}

///Build A dot u where u =[z,c]
void v2RDMSolver::bpsdp_Au(SharedVector A, SharedVector u){

    //A->zero();  
    memset((void*)A->pointer(),'\0',nconstraints*sizeof(double));

    offset = 0;
    double start = omp_get_wtime();
    D2_constraints_Au(A,u);
    double end = omp_get_wtime();
    d2timeAu += (end - start);

    if ( constrain_q2 ) {
        start = omp_get_wtime();
        if ( !spin_adapt_q2 ) {
            Q2_constraints_Au(A,u);
        }else {
            Q2_constraints_Au_spin_adapted(A,u);
        }
        end = omp_get_wtime();
        q2timeAu += (end - start);
    }

    if ( constrain_g2 ) {
        start = omp_get_wtime();
        if ( ! spin_adapt_g2 ) {
            G2_constraints_Au(A,u);
        }else {
            G2_constraints_Au_spin_adapted(A,u);
        }
        end = omp_get_wtime();
        g2timeAu += (end - start);
    }

    if ( constrain_t1 ) {
        T1_constraints_Au(A,u);
    }

    if ( constrain_t2 ) {
        double start = omp_get_wtime();
        //T2_constraints_Au(A,u);
        T2_constraints_Au_slow(A,u);
        double end = omp_get_wtime();
        t2timeAu += end - start;
    }

    if ( constrain_d3 ) {
        D3_constraints_Au(A,u);
    }

} // end Au
void v2RDMSolver::bpsdp_Au_slow(SharedVector A, SharedVector u){

    //A->zero();  
    memset((void*)A->pointer(),'\0',nconstraints*sizeof(double));

    offset = 0;
    double start = omp_get_wtime();
    D2_constraints_Au(A,u);
    double end = omp_get_wtime();
    d2timeAu += (end - start);

    if ( constrain_q2 ) {
        start = omp_get_wtime();
        if ( !spin_adapt_q2 ) {
            Q2_constraints_Au(A,u);
        }else {
            Q2_constraints_Au_spin_adapted(A,u);
        }
        end = omp_get_wtime();
        q2timeAu += (end - start);
    }

    if ( constrain_g2 ) {
        start = omp_get_wtime();
        if ( ! spin_adapt_g2 ) {
            G2_constraints_Au(A,u);
        }else {
            G2_constraints_Au_spin_adapted(A,u);
        }
        end = omp_get_wtime();
        g2timeAu += (end - start);
    }

    if ( constrain_t1 ) {
        T1_constraints_Au(A,u);
    }

    if ( constrain_t2 ) {
        double start = omp_get_wtime();
        //T2_constraints_Au(A,u);
        T2_constraints_Au_slow(A,u);
        double end = omp_get_wtime();
        t2timeAu += end - start;
    }
    if ( constrain_d3 ) {
        D3_constraints_Au(A,u);
    }

} // end Au

// D2 portion of A^T.y ( and D1 / Q1 )
void v2RDMSolver::D2_constraints_ATu(SharedVector A,SharedVector u){
    double* A_p = A->pointer();
    double* u_p = u->pointer();

    if ( constrain_spin ) {
        // spin
        for (int i = 0; i < nmo; i++){
            for (int j = 0; j < nmo; j++){
                int h = SymmetryPair(symmetry[i],symmetry[j]);
                if ( gems_ab[h] == 0 ) continue;
                int ij = ibas_ab_sym[h][i][j];
                int ji = ibas_ab_sym[h][j][i];
                A_p[d2aboff[h] + ij*gems_ab[h]+ji] += u_p[offset];
            }
        }
        offset++;
    }

    // Traces
    // Tr(D2ab)
    for (int i = 0; i < nmo; i++){
        for (int j = 0; j < nmo; j++){
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            int ij = ibas_ab_sym[h][i][j];
            if ( gems_ab[h] == 0 ) continue;
            A_p[d2aboff[h] + ij*gems_ab[h]+ij] += u_p[offset];
        }
    }
    offset++;

    // Tr(D2aa)
    for (int i = 0; i < nmo; i++){
        for (int j = 0; j < nmo; j++){
            if ( i==j ) continue;
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            if ( gems_aa[h] == 0 ) continue;
            int ij = ibas_aa_sym[h][i][j];
            A_p[d2aaoff[h]+ij*gems_aa[h]+ij] += u_p[offset];
        }
    }
    offset++;

    // Tr(D2bb)
    for (int i = 0; i < nmo; i++){
        for (int j = 0; j < nmo; j++){
            if ( i==j ) continue;
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            if ( gems_aa[h] == 0 ) continue;
            int ij = ibas_aa_sym[h][i][j];
            A_p[d2bboff[h]+ij*gems_aa[h]+ij] += u_p[offset];
        }
    }
    offset++;

    // d1 / q1 a
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                double dum = u_p[offset + i*amopi_[h]+j];
                A_p[d1aoff[h] + j*amopi_[h]+i] += dum;
                A_p[q1aoff[h] + i*amopi_[h]+j] += dum;
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    // d1 / q1 b
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                double dum = u_p[offset + i*amopi_[h]+j];
                A_p[d1boff[h] + j*amopi_[h]+i] += dum;
                A_p[q1boff[h] + i*amopi_[h]+j] += dum;
            }
        }
        offset += amopi_[h]*amopi_[h];
    }

    int na = nalpha_ - nfrzc;
    int nb = nbeta_ - nfrzc;

    int poff = 0;

    // contraction: D2ab -> D1 a
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                A_p[d1aoff[h] + i*amopi_[h]+j] += nb * u_p[offset + i*amopi_[h]+j];
                int ii = i + poff;
                int jj = j + poff;
                for (int k = 0; k < nmo; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][ii][k];
                    int jk = ibas_ab_sym[h2][jj][k];
                    A_p[d2aboff[h2] + ik*gems_ab[h2]+jk] -= u_p[offset + i*amopi_[h]+j];
                }
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }

    // contraction: D2ab -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A_p[d1boff[h] + i*amopi_[h]+j] += na * u_p[offset + i*amopi_[h]+j];
                int ii = i + poff;
                int jj = j + poff;
                for(int k = 0; k < nmo; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][k][ii];
                    int jk = ibas_ab_sym[h2][k][jj];
                    A_p[d2aboff[h2] + ik*gems_ab[h2]+jk] -= u_p[offset + i*amopi_[h]+j];
                }
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }

    //contract D2aa -> D1 a
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A_p[d1aoff[h] + i*amopi_[h]+j] += (na - 1.0) * u_p[offset + i*amopi_[h]+j];
                int ii = i + poff;
                int jj = j + poff;
                for(int k =0; k < nmo; k++){
                    if( ii==k || jj==k )continue;
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_aa_sym[h2][ii][k];
                    int jk = ibas_aa_sym[h2][jj][k];
                    int sik = ( ii < k ? 1 : -1);
                    int sjk = ( jj < k ? 1 : -1);
                    A_p[d2aaoff[h2] + ik*gems_aa[h2]+jk] -= sik*sjk*u_p[offset + i*amopi_[h]+j];
                }
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }

    //contract D2bb -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A_p[d1boff[h] + i*amopi_[h]+j] += (nb - 1.0) * u_p[offset + i*amopi_[h]+j];
                int ii = i + poff;
                int jj = j + poff;
                for(int k =0; k < nmo; k++){
                    if( ii==k || jj==k )continue;
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_aa_sym[h2][ii][k];
                    int jk = ibas_aa_sym[h2][jj][k];
                    int sik = ( ii < k ? 1 : -1);
                    int sjk = ( jj < k ? 1 : -1);
                    A_p[d2bboff[h2] + ik*gems_aa[h2]+jk] -= sik*sjk*u_p[offset + i*amopi_[h]+j];
                }
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
    }

}

///Build AT dot u where u =[z,c]
void v2RDMSolver::bpsdp_ATu(SharedVector A, SharedVector u){

    //A->zero();
    memset((void*)A->pointer(),'\0',dimx*sizeof(double));

    offset = 0;
    double start = omp_get_wtime();
    D2_constraints_ATu(A,u);
    double end = omp_get_wtime();
    d2timeATu += (end - start);

    if ( constrain_q2 ) {
        start = omp_get_wtime();
        if ( !spin_adapt_q2 ) {
            Q2_constraints_ATu(A,u);
        }else {
            Q2_constraints_ATu_spin_adapted(A,u);
        }
        end = omp_get_wtime();
        q2timeATu += (end - start);
    }

    if ( constrain_g2 ) {
        start = omp_get_wtime();
        if ( ! spin_adapt_g2 ) {
            G2_constraints_ATu(A,u);
        }else {
            G2_constraints_ATu_spin_adapted(A,u);
        }
        end = omp_get_wtime();
        g2timeATu += (end - start);
    }

    if ( constrain_t1 ) {
        T1_constraints_ATu(A,u);
    }

    if ( constrain_t2 ) {
        double start = omp_get_wtime();
        //T2_constraints_ATu(A,u);
        T2_constraints_ATu_slow(A,u);
        double end = omp_get_wtime();
        t2timeATu += end - start;
    }
    if ( constrain_d3 ) {
        D3_constraints_ATu(A,u);
    }

}//end ATu

void v2RDMSolver::bpsdp_ATu_slow(SharedVector A, SharedVector u){

    //A->zero();
    memset((void*)A->pointer(),'\0',dimx*sizeof(double));

    offset = 0;
    double start = omp_get_wtime();
    D2_constraints_ATu(A,u);
    double end = omp_get_wtime();
    d2timeATu += (end - start);

    if ( constrain_q2 ) {
        start = omp_get_wtime();
        if ( !spin_adapt_q2 ) {
            Q2_constraints_ATu(A,u);
        }else {
            Q2_constraints_ATu_spin_adapted(A,u);
        }
        end = omp_get_wtime();
        q2timeATu += (end - start);
    }

    if ( constrain_g2 ) {
        start = omp_get_wtime();
        if ( ! spin_adapt_g2 ) {
            G2_constraints_ATu(A,u);
        }else {
            G2_constraints_ATu_spin_adapted(A,u);
        }
        end = omp_get_wtime();
        g2timeATu += (end - start);
    }

    if ( constrain_t1 ) {
        T1_constraints_ATu(A,u);
    }

    if ( constrain_t2 ) {
        double start = omp_get_wtime();
        //T2_constraints_ATu(A,u);
        T2_constraints_ATu_slow(A,u);
        double end = omp_get_wtime();
        t2timeATu += end - start;
    }

    if ( constrain_d3 ) {
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
    for (int i = 0; i < dimensions.size(); i++) {
        if ( dimensions[i] == 0 ) continue;
        int myoffset = 0;
        for (int j = 0; j < i; j++) {
            myoffset += dimensions[j] * dimensions[j];
        }

        boost::shared_ptr<Matrix> mat    (new Matrix(dimensions[i],dimensions[i]));
        boost::shared_ptr<Matrix> eigvec (new Matrix(dimensions[i],dimensions[i]));
        boost::shared_ptr<Vector> eigval (new Vector(dimensions[i]));
        boost::shared_ptr<Vector> Up     (new Vector(dimensions[i]));
        boost::shared_ptr<Vector> Um     (new Vector(dimensions[i]));
        double ** mat_p = mat->pointer();
        double * A_p   = ATy->pointer();

        //C_DCOPY(dimensions[i]*dimensions[i],&A_p[myoffset],1,&mat_p[0][0],1);
        for (int p = 0; p < dimensions[i]; p++) {
            for (int q = p; q < dimensions[i]; q++) {
                double dum = 0.5 * ( A_p[myoffset + p * dimensions[i] + q] +
                                     A_p[myoffset + q * dimensions[i] + p] );
                mat_p[p][q] = mat_p[q][p] = dum;
                 
            }
        }
        mat->diagonalize(eigvec,eigval);

        // separate U+ and U-
        double * u_p    = Up->pointer();
        double * u_m    = Um->pointer();
        double * eval_p = eigval->pointer();
        for (int p = 0; p < dimensions[i]; p++) {
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
        for (int pq = 0; pq < dimensions[i] * dimensions[i]; pq++) {

            int q = pq % dimensions[i];
            int p = (pq-q) / dimensions[i];
            if ( p > q ) continue;

            double sumx = 0.0;
            double sumz = 0.0;
            for (int j = 0; j < dimensions[i]; j++) {
                sumx += u_p[j] * evec_p[p][j] * evec_p[q][j];
                sumz += u_m[j] * evec_p[p][j] * evec_p[q][j];
            }
            x_p[myoffset+p*dimensions[i]+q] = sumx;
            z_p[myoffset+p*dimensions[i]+q] = sumz;

            x_p[myoffset+q*dimensions[i]+p] = sumx;
            z_p[myoffset+q*dimensions[i]+p] = sumz;
        }
    }
}

void v2RDMSolver::UnpackDensityPlusCore() {

    memset((void*)d2_plus_core_sym,'\0',d2_plus_core_dim*sizeof(double));

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
                d2_plus_core_sym[id] = val;
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

                        //if ( fabs(d2_plus_core_sym[INDEX(iifull,jjfull)]) < 1e-10 ) {
                        //    en += tei_full_sym[INDEX(iifull,jjfull)];
                        //}

                        d2_plus_core_sym[INDEX(iifull,jjfull)] =  1.0;

                    }else {

                        int iifull = ibas_full_sym[0][ifull][ifull];
                        int jjfull = ibas_full_sym[0][jfull][jfull];

                        //if ( fabs(d2_plus_core_sym[INDEX(iifull,jjfull)]) < 1e-10 ) {
                        //    en += 4.0 * tei_full_sym[INDEX(iifull,jjfull)];
                        //}

                        d2_plus_core_sym[INDEX(iifull,jjfull)] =  4.0;

                        int offset2 = 0;
                        for (int myh = 0; myh < hij; myh++) {
                            offset2 += gems_full[myh] * ( gems_full[myh] + 1 ) / 2;
                        }
                        offset = 0;
                        for (int myh = 0; myh < hij; myh++) {
                            offset += gems_plus_core[myh] * ( gems_plus_core[myh] + 1 ) / 2;
                        }

                        int ijfull = ibas_full_sym[hij][ifull][jfull];
                        //if ( fabs(d2_plus_core_sym[offset + INDEX(ijfull,ijfull)]) < 1e-10 ) {
                        //    en -= 2.0 * tei_full_sym[offset2 + INDEX(ijfull,ijfull)];
                        //}

                        d2_plus_core_sym[offset + INDEX(ijfull,ijfull)] = -2.0;
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

                        d2_plus_core_sym[id] = val;

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

                        d2_plus_core_sym[offset + INDEX(ilfull,jifull)] = val;

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
                
                d1_plus_core_sym[id]  = x_p[d1aoff[h] + i * amopi_[h] + j];
                d1_plus_core_sym[id] += x_p[d1boff[h] + i * amopi_[h] + j];

                // scale off-diagonal elements
                if ( i != j ) {
                    d1_plus_core_sym[id] *= 2.0;
                }

            }
        }
        offset += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }

    // core; core;
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            d1_plus_core_sym[offset + INDEX(i,i)] = 2.0;
        }
        offset += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }

    // check energy:

    // two-electron part:
    //double en2 = 0.0;
    //int offset_full      = 0;
    //int offset_plus_core = 0;
    //for (int h = 0; h < nirrep_; h++) {
    //    en2 += C_DDOT(gems_plus_core[h]*(gems_plus_core[h]+1)/2,tei_full_sym+offset_full,1,d2_plus_core_sym+offset_plus_core,1);
    //    offset_full += gems_full[h] * ( gems_full[h] + 1 ) / 2;
    //    offset_plus_core += gems_plus_core[h] * ( gems_plus_core[h] + 1 ) / 2;
    //}

    // one-electron part:
    //double en1 = 0.0;
    //offset_full      = 0;
    //offset_plus_core = 0;
    //for (int h = 0; h < nirrep_; h++) {
    //    en1 += C_DDOT( ( frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2 ,oei_full_sym + offset_full,1,d1_plus_core_sym + offset_plus_core,1);
    //    offset_full      += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    //    offset_plus_core += (frzcpi_[h] + amopi_[h]) * ( frzcpi_[h] + amopi_[h] + 1 ) / 2;
    //}

    //printf("\n");
    //printf("efrzc1       %20.12lf\n",efrzc1);
    //printf("efrzc2       %20.12lf\n",efrzc2);
    //printf("\n");
    //printf("en1          %20.12lf\n",en1);
    //printf("en2          %20.12lf\n",en2);
    //printf("en1+en2      %20.12lf\n",en1+en2);
    //printf("en1+en2+enuc %20.12lf\n",en1+en2 + enuc);
}

// repack rotated full-space integrals into active-space integrals
void v2RDMSolver::RepackIntegralsDF(){

    long int full = nmo + nfrzc + nfrzv;

    // if frozen core, adjust oei's and compute frozen core energy:
    efrzc1 = 0.0;
    efrzc2 = 0.0;
    offset = 0;
    long int offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < frzcpi_[h]; i++) {
            long int ii = i + offset;
            efrzc1 += 2.0 * oei_full_sym[offset3 + INDEX(i,i)];

            long int offset2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (long int j = 0; j < frzcpi_[h2]; j++) {
                    long int jj = j + offset2;
                    double dum1 = C_DDOT(nQ_,Qmo_ + nQ_*(ii*nso_+ii),1,Qmo_+nQ_*(jj*nso_+jj),1);
                    double dum2 = C_DDOT(nQ_,Qmo_ + nQ_*(ii*nso_+jj),1,Qmo_+nQ_*(ii*nso_+jj),1);
                    efrzc2 += 2.0 * dum1 - dum2;
                }
                offset2 += nmopi_[h2];
            }
        }
        offset += nmopi_[h];
        offset3 += nmopi_[h] * (nmopi_[h] + 1 ) / 2;
    }
    efrzc = efrzc1 + efrzc2;
    //printf("%20.12lf\n",efrzc1);
    //printf("%20.12lf\n",efrzc2);
    //printf("%20.12lf\n",efrzc);
    //exit(0);

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
                        int kk = k + offset2;
                        dum1 += C_DDOT(nQ_,Qmo_ + nQ_*(ii*nso_+jj),1,Qmo_+nQ_*(kk*nso_+kk),1);
                        dum2 += C_DDOT(nQ_,Qmo_ + nQ_*(ii*nso_+kk),1,Qmo_+nQ_*(jj*nso_+kk),1);
                    }
                    offset2 += nmopi_[h2];
                }
                c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = oei_full_sym[offset3+INDEX(i,j)];
                c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = oei_full_sym[offset3+INDEX(i,j)];

                c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] += 2.0 * dum1 - dum2;
                c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] += 2.0 * dum1 - dum2;
            }
        }
        offset += nmopi_[h];
        offset3 += nmopi_[h] * (nmopi_[h]+1)/2;
    }

    // two-electron part
    long int na = nalpha_ - nfrzc;
    long int nb = nbeta_ - nfrzc;
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

                c_p[d2aboff[h] + ij*gems_ab[h]+kl] = C_DDOT(nQ_,Qmo_+nQ_*(ii*nso_+kk),1,Qmo_+nQ_*(jj*nso_+ll),1);
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

                double dum1 = C_DDOT(nQ_,Qmo_+nQ_*(ii*nso_+kk),1,Qmo_+nQ_*(jj*nso_+ll),1);
                double dum2 = C_DDOT(nQ_,Qmo_+nQ_*(ii*nso_+ll),1,Qmo_+nQ_*(jj*nso_+kk),1);

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]    = dum1 - dum2;
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]    = dum1 - dum2;
            }
        }
    }
}

// repack rotated full-space integrals into active-space integrals
void v2RDMSolver::RepackIntegrals(){

    long int full = nmo + nfrzc + nfrzv;

    // if frozen core, adjust oei's and compute frozen core energy:
    efrzc1 = 0.0;
    efrzc2 = 0.0;
    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {

            int ifull = i + pitzer_offset_full[h];
            int ii    = ibas_full_sym[0][ifull][ifull];

            efrzc1 += 2.0 * oei_full_sym[offset + INDEX(i,i)]; 

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
                    efrzc2 += (2.0 * tei_full_sym[INDEX(ii,jj)] - tei_full_sym[myoff + INDEX(ij,ij)]);
                }
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }
    efrzc = efrzc1 + efrzc2;

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

                        dum += (2.0 * tei_full_sym[INDEX(ij,kk)] - tei_full_sym[myoff+INDEX(ik,jk)]);

                    }
                }
                c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = oei_full_sym[offset+INDEX(i,j)];
                c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] = oei_full_sym[offset+INDEX(i,j)];
                c_p[d1aoff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] += dum;
                c_p[d1boff[h] + (i-frzcpi_[h])*amopi_[h] + (j-frzcpi_[h])] += dum;
            }
        }
        offset += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;
    }

    int na = nalpha_ - nfrzc;
    int nb = nbeta_ - nfrzc;
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

                c_p[d2aboff[h] + ij*gems_ab[h]+kl]    = tei_full_sym[myoff + INDEX(ik,jl)];

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

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]    = 0.5 * tei_full_sym[myoff + INDEX(ik,jl)];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   += 0.5 * tei_full_sym[myoff + INDEX(jl,ik)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]    = 0.5 * tei_full_sym[myoff + INDEX(ik,jl)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   += 0.5 * tei_full_sym[myoff + INDEX(jl,ik)];

                int hil = SymmetryPair(symmetry_full[ifull],symmetry_full[lfull]);
                int il  = ibas_full_sym[hil][ifull][lfull];
                int jk  = ibas_full_sym[hil][jfull][kfull];

                myoff = 0;
                for (int myh = 0; myh < hil; myh++) {
                    myoff += gems_full[myh] * ( gems_full[myh] + 1 ) / 2;
                }

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym[myoff + INDEX(il,jk)];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym[myoff + INDEX(jk,il)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym[myoff + INDEX(il,jk)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym[myoff + INDEX(jk,il)];

            }
        }
    }
}

void v2RDMSolver::FinalTransformationMatrix() {

    // test: transform oeis
    /*offset = 0;
    double diff = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ieo = 0; ieo < nmo+nfrzc+nfrzv; ieo++) {
            int ifull = energy_to_pitzer_order[ieo];
            int hi    = symmetry_full[ifull];
            if ( h != hi ) continue;
            int i     = ifull - pitzer_offset_full[hi];
            for (int jeo = 0; jeo < nmo+nfrzc+nfrzv; jeo++) {
                int jfull = energy_to_pitzer_order[jeo];
                int hj    = symmetry_full[jfull];
                if ( h != hj ) continue;
                int j     = jfull - pitzer_offset_full[hj];
                double dum = 0.0;
                for (int keo = 0; keo < nmo+nfrzc+nfrzv; keo++) {
                    int kfull = energy_to_pitzer_order[keo];
                    int hk    = symmetry_full[kfull];
                    if ( h != hk ) continue;
                    int k     = kfull - pitzer_offset_full[hk];
                    for (int leo = 0; leo < nmo+nfrzc+nfrzv; leo++) {
                        int lfull = energy_to_pitzer_order[leo];
                        int hl    = symmetry_full[lfull];
                        if ( h != hl ) continue;
                        int l     = lfull - pitzer_offset_full[hl];
                        //dum += saveK1->pointer(h)[k][l] * jacobi_transformation_matrix_[ieo*(nmo+nfrzc+nfrzv)+keo]
                        //                                * jacobi_transformation_matrix_[jeo*(nmo+nfrzc+nfrzv)+leo];
                        dum += saveK1->pointer(h)[k][l] * jacobi_transformation_matrix_[keo*(nmo+nfrzc+nfrzv)+ieo]
                                                        * jacobi_transformation_matrix_[leo*(nmo+nfrzc+nfrzv)+jeo];
                    }
                }
                diff += (dum - oei_full_sym[offset+INDEX(i,j)]) * (dum - oei_full_sym[offset+INDEX(i,j)]);
                printf("%5i %5i %20.12lf %20.12lf\n",i,j,dum,oei_full_sym[offset+INDEX(i,j)]);
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
            for (int ieo = 0; ieo < nmo+nfrzc+nfrzv; ieo++) {
                int ifull = energy_to_pitzer_order[ieo];
                int hi    = symmetry_full[ifull];
                if ( h != hi ) continue;
                int i     = ifull - pitzer_offset_full[hi];

                double dum = 0.0;

                // old basis function j in energy order
                for (int jeo = 0; jeo < nmo+nfrzc+nfrzv; jeo++) {
                    int jfull = energy_to_pitzer_order[jeo];
                    int hj    = symmetry_full[jfull];
                    if ( h != hj ) continue;
                    int j     = jfull - pitzer_offset_full[hj];

                    dum += ca_p[mu][j] * jacobi_transformation_matrix_[jeo*(nmo+nfrzc+nfrzv)+ieo];
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
    //for (int i = 0; i < nfrzc + nmo + nfrzv; i++) {
    //    for (int j = 0; j < nfrzc + nmo + nfrzv; j++) {
    //        printf("%5i %5i %20.12lf\n",i,j,jacobi_transformation_matrix_[i*(nfrzc+nmo+nfrzv)+j]);
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

    Jacobi(jacobi_transformation_matrix_,
          oei_full_sym,oei_full_dim,tei_full_sym,tei_full_dim,
          d1_plus_core_sym,d1_plus_core_dim,d2_plus_core_sym,d2_plus_core_dim,
          symmetry_energy_order,nfrzc,nmo,nfrzv,nirrep_,
          jacobi_data_,jacobi_outfile_);

    outfile->Printf("            Jacobi Optimization %s.\n",(int)jacobi_data_[9] ? "converged" : "did not converge");
    outfile->Printf("            Total energy change: %11.6le\n",jacobi_data_[8]);
    outfile->Printf("\n");

    if ( fabs(jacobi_data_[8]) < jacobi_data_[4] ) {
        jacobi_converged_ = true;
    }

    if ( is_df_ ) {
        RepackIntegralsDF();
    }else {
        RepackIntegrals();
    }
}

}} //end namespaces
