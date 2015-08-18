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
    // tei_full_dim = nQ_*nso_*(nso_+1)/2;

    // gidofalvi -- modified this so that number of nnz 2-e integrals is correct for large bases
    tei_full_dim = (long int) nQ_ * (long int) nso_ * ( (long int) nso_ + 1 ) /2 ;

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

void v2RDMSolver::TEI() {

    long int full = nmo + nfrzc + nfrzv;
    double * temptei = (double*)malloc(full*full*full*full*sizeof(double));
    memset((void*)temptei,'\0',full*full*full*full*sizeof(double));

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
    //printf("efrzc1       %20.12lf\n",efrzc1);
    //printf("efrzc2       %20.12lf\n",efrzc2);

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

}}
