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
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>
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

void v2RDMSolver::Q2_constraints_guess_spin_adapted(SharedVector u){

    double * u_p = u->pointer();

    // map D2ab to Q2s
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_00[h]; ij++) {
            int i = bas_00_sym[h][ij][0];
            int j = bas_00_sym[h][ij][1];
            int ijd = ibas_ab_sym[h][i][j];
            int jid = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_00[h]; kl++) {
                int k = bas_00_sym[h][kl][0];
                int l = bas_00_sym[h][kl][1];

                double dum  = 0.0;

                int kld = ibas_ab_sym[h][k][l];
                int lkd = ibas_ab_sym[h][l][k];
                dum        +=  0.5 * u_p[d2aboff[h] + kld*gems_ab[h]+ijd];          // +D2(kl,ij)
                dum        +=  0.5 * u_p[d2aboff[h] + lkd*gems_ab[h]+ijd];          // +D2(lk,ij)
                dum        +=  0.5 * u_p[d2aboff[h] + kld*gems_ab[h]+jid];          // +D2(kl,ji)
                dum        +=  0.5 * u_p[d2aboff[h] + lkd*gems_ab[h]+jid];          // +D2(lk,ji)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ii*amopi_[h2]+kk]; // +Q1(i,k) djl
                    dum        -=  0.5 * u_p[d1boff[h2] + ii*amopi_[h2]+kk]; // -D1(i,k) djl
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ll*amopi_[h2]+jj]; // +Q1(l,j) dik
                    dum        -=  0.5 * u_p[d1boff[h2] + ll*amopi_[h2]+jj]; // -D1(l,j) dik
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ll*amopi_[h2]+ii]; // +Q1(l,i) djk
                    dum        -=  0.5 * u_p[d1boff[h2] + ll*amopi_[h2]+ii]; // -D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + kk*amopi_[h2]+jj]; // +Q1(k,j) dil
                    dum        -=  0.5 * u_p[d1boff[h2] + kk*amopi_[h2]+jj]; // -D1(k,j) dil
                }
                u_p[q2soff[h] + ij*gems_00[h]+kl] = dum;
            }
        }
        offset += gems_00[h]*gems_00[h];
    }
    // map D2ab to Q210
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i   =  bas_aa_sym[h][ij][0];
            int j   =  bas_aa_sym[h][ij][1];
            int ijd = ibas_ab_sym[h][i][j];
            int jid = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int   k =  bas_aa_sym[h][kl][0];
                int   l =  bas_aa_sym[h][kl][1];

                double dum  = 0.0;

                // not spin adapted
                int kld = ibas_ab_sym[h][k][l];
                int lkd = ibas_ab_sym[h][l][k];
                dum        +=  0.5 * u_p[d2aboff[h] + kld*gems_ab[h]+ijd];          // +D2(kl,ij)
                dum        -=  0.5 * u_p[d2aboff[h] + lkd*gems_ab[h]+ijd];          // -D2(lk,ij)
                dum        -=  0.5 * u_p[d2aboff[h] + kld*gems_ab[h]+jid];          // -D2(kl,ji)
                dum        +=  0.5 * u_p[d2aboff[h] + lkd*gems_ab[h]+jid];          // +D2(lk,ji)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ii*amopi_[h2]+kk]; // +Q1(i,k) djl
                    dum        -=  0.5 * u_p[d1boff[h2] + ii*amopi_[h2]+kk]; // -D1(i,k) djl
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ll*amopi_[h2]+jj]; // +Q1(l,j) dik
                    dum        -=  0.5 * u_p[d1boff[h2] + ll*amopi_[h2]+jj]; // -D1(l,j) dik
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  0.5 * u_p[q1aoff[h2] + ll*amopi_[h2]+ii]; // -Q1(l,i) djk
                    dum        +=  0.5 * u_p[d1boff[h2] + ll*amopi_[h2]+ii]; // +D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  0.5 * u_p[q1aoff[h2] + kk*amopi_[h2]+jj]; // -Q1(k,j) dil
                    dum        +=  0.5 * u_p[d1boff[h2] + kk*amopi_[h2]+jj]; // +D1(k,j) dil
                }

                u_p[q2toff[h] + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
    // map D2aa to Q211
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];

                double dum  = 0.0;
                dum        +=  u_p[d2aaoff[h] + kl*gems_aa[h]+ij];    // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[q1aoff[h2] + ii*amopi_[h2]+kk];  // +Q1(i,k) djl
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  u_p[d1aoff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  u_p[q1aoff[h2] + jj*amopi_[h2]+kk];  // -Q1(j,k) dil
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1aoff[h2] + ll*amopi_[h2]+jj];  // -D1(l,j) dkl
                }

                u_p[q2toff_p1[h] + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
    // map D2bb to Q21-1
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];

                double dum  = 0.0;

                dum        +=  u_p[d2bboff[h] + kl*gems_aa[h]+ij];    // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[q1boff[h2] + ii*amopi_[h2]+kk];  // +Q1(i,k) djl
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  u_p[d1boff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  u_p[q1boff[h2] + jj*amopi_[h2]+kk];  // -Q1(j,k) dil
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1boff[h2] + ll*amopi_[h2]+jj];  // -D1(l,j) dkl
                }

                u_p[q2toff_m1[h] + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
}

void v2RDMSolver::Q2_constraints_Au_spin_adapted(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    // map D2ab to Q2s
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_00[h]; ij++) {
            int i = bas_00_sym[h][ij][0];
            int j = bas_00_sym[h][ij][1];
            int ijd = ibas_ab_sym[h][i][j];
            int jid = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_00[h]; kl++) {
                int k = bas_00_sym[h][kl][0];
                int l = bas_00_sym[h][kl][1];

                double dum  = -u_p[q2soff[h] + ij*gems_00[h]+kl];          // -Q2(ij,kl)

                // not spin adapted
                int kld = ibas_ab_sym[h][k][l];
                int lkd = ibas_ab_sym[h][l][k];
                dum        +=  0.5 * u_p[d2aboff[h] + kld*gems_ab[h]+ijd];          // +D2(kl,ij)
                dum        +=  0.5 * u_p[d2aboff[h] + lkd*gems_ab[h]+ijd];          // +D2(lk,ij)
                dum        +=  0.5 * u_p[d2aboff[h] + kld*gems_ab[h]+jid];          // +D2(kl,ji)
                dum        +=  0.5 * u_p[d2aboff[h] + lkd*gems_ab[h]+jid];          // +D2(lk,ji)

                // spin adapted
                //dum        +=  u_p[d2soff[h] + INDEX(kl,ij)];          // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ii*amopi_[h2]+kk]; // +Q1(i,k) djl
                    dum        -=  0.5 * u_p[d1boff[h2] + ii*amopi_[h2]+kk]; // -D1(i,k) djl
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ll*amopi_[h2]+jj]; // +Q1(l,j) dik
                    dum        -=  0.5 * u_p[d1boff[h2] + ll*amopi_[h2]+jj]; // -D1(l,j) dik
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ll*amopi_[h2]+ii]; // +Q1(l,i) djk
                    dum        -=  0.5 * u_p[d1boff[h2] + ll*amopi_[h2]+ii]; // -D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + kk*amopi_[h2]+jj]; // +Q1(k,j) dil
                    dum        -=  0.5 * u_p[d1boff[h2] + kk*amopi_[h2]+jj]; // -D1(k,j) dil
                }

                A_p[offset + ij*gems_00[h]+kl] = dum;
            }
        }
        offset += gems_00[h]*gems_00[h];
    }
    // map D2ab to Q210
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i   =  bas_aa_sym[h][ij][0];
            int j   =  bas_aa_sym[h][ij][1];
            int ijd = ibas_ab_sym[h][i][j];
            int jid = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int   k =  bas_aa_sym[h][kl][0];
                int   l =  bas_aa_sym[h][kl][1];

                double dum  = -u_p[q2toff[h] + ij*gems_aa[h]+kl];          // -Q2(ij,kl)

                // not spin adapted
                int kld = ibas_ab_sym[h][k][l];
                int lkd = ibas_ab_sym[h][l][k];
                dum        +=  0.5 * u_p[d2aboff[h] + kld*gems_ab[h]+ijd];          // +D2(kl,ij)
                dum        -=  0.5 * u_p[d2aboff[h] + lkd*gems_ab[h]+ijd];          // -D2(lk,ij)
                dum        -=  0.5 * u_p[d2aboff[h] + kld*gems_ab[h]+jid];          // -D2(kl,ji)
                dum        +=  0.5 * u_p[d2aboff[h] + lkd*gems_ab[h]+jid];          // +D2(lk,ji)

                // spin adapted
                //dum        +=  u_p[d2toff[h] + INDEX(kl,ij)];          // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ii*amopi_[h2]+kk]; // +Q1(i,k) djl
                    dum        -=  0.5 * u_p[d1boff[h2] + ii*amopi_[h2]+kk]; // -D1(i,k) djl
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  0.5 * u_p[q1aoff[h2] + ll*amopi_[h2]+jj]; // +Q1(l,j) dik
                    dum        -=  0.5 * u_p[d1boff[h2] + ll*amopi_[h2]+jj]; // -D1(l,j) dik
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  0.5 * u_p[q1aoff[h2] + ll*amopi_[h2]+ii]; // -Q1(l,i) djk
                    dum        +=  0.5 * u_p[d1boff[h2] + ll*amopi_[h2]+ii]; // +D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  0.5 * u_p[q1aoff[h2] + kk*amopi_[h2]+jj]; // -Q1(k,j) dil
                    dum        +=  0.5 * u_p[d1boff[h2] + kk*amopi_[h2]+jj]; // +D1(k,j) dil
                }

                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
    // map D2aa to Q211
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum  = -u_p[q2toff_p1[h] + ij*gems_aa[h]+kl];    // -Q2(ij,kl)
                dum        +=  u_p[d2aaoff[h] + kl*gems_aa[h]+ij];    // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[q1aoff[h2] + ii*amopi_[h2]+kk];  // +Q1(i,k) djl
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  u_p[d1aoff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  u_p[q1aoff[h2] + jj*amopi_[h2]+kk];  // -Q1(j,k) dil
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1aoff[h2] + ll*amopi_[h2]+jj];  // -D1(l,j) dkl
                }

                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
    // map D2bb to Q21-1
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum  = -u_p[q2toff_m1[h] + ij*gems_aa[h]+kl];    // -Q2(ij,kl)
                dum        +=  u_p[d2bboff[h] + kl*gems_aa[h]+ij];    // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[q1boff[h2] + ii*amopi_[h2]+kk];  // +Q1(i,k) djl
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  u_p[d1boff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  u_p[q1boff[h2] + jj*amopi_[h2]+kk];  // -Q1(j,k) dil
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1boff[h2] + ll*amopi_[h2]+jj];  // -D1(l,j) dkl
                }

                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
}

// Q2 portion of A^T.y (spin adapted)
void v2RDMSolver::Q2_constraints_ATu_spin_adapted(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    // map D2ab to Q2s
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_00[h]; ij++) {
            int i = bas_00_sym[h][ij][0];
            int j = bas_00_sym[h][ij][1];
            int ijd = ibas_ab_sym[h][i][j];
            int jid = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_00[h]; kl++) {
                int k = bas_00_sym[h][kl][0];
                int l = bas_00_sym[h][kl][1];

                double dum  = u_p[offset + ij*gems_00[h]+kl];

                A_p[q2soff[h] + ij*gems_00[h]+kl]    -= dum;          // -Q2(ij,kl)

                // not spin adapted
                int kld = ibas_ab_sym[h][k][l];
                int lkd = ibas_ab_sym[h][l][k];
                A_p[d2aboff[h] + kld*gems_ab[h]+ijd] += 0.5 * dum;          // +D2(kl,ij)
                A_p[d2aboff[h] + lkd*gems_ab[h]+ijd] += 0.5 * dum;          // +D2(lk,ij)
                A_p[d2aboff[h] + kld*gems_ab[h]+jid] += 0.5 * dum;          // +D2(kl,ji)
                A_p[d2aboff[h] + lkd*gems_ab[h]+jid] += 0.5 * dum;          // +D2(lk,ji)

                // spin adapted
                //A_p[d2soff[h] + INDEX(kl,ij)] += dum;          // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[q1aoff[h2] + ii*amopi_[h2]+kk] += 0.5 * dum; // +Q1(i,k) djl
                    A_p[d1boff[h2] + ii*amopi_[h2]+kk] -= 0.5 * dum; // -D1(i,k) djl
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[q1aoff[h2] + ll*amopi_[h2]+jj] += 0.5 * dum; // +Q1(l,j) dik
                    A_p[d1boff[h2] + ll*amopi_[h2]+jj] -= 0.5 * dum; // -D1(l,j) dik
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[q1aoff[h2] + ll*amopi_[h2]+ii] += 0.5 * dum; // +Q1(l,i) djk
                    A_p[d1boff[h2] + ll*amopi_[h2]+ii] -= 0.5 * dum; // -D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[q1aoff[h2] + kk*amopi_[h2]+jj] += 0.5 * dum; // +Q1(k,j) dil
                    A_p[d1boff[h2] + kk*amopi_[h2]+jj] -= 0.5 * dum; // -D1(k,j) dil
                }
            }
        }
        offset += gems_00[h]*gems_00[h];
    }
    // map D2ab to Q210
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i   =  bas_aa_sym[h][ij][0];
            int j   =  bas_aa_sym[h][ij][1];
            int ijd = ibas_ab_sym[h][i][j];
            int jid = ibas_ab_sym[h][j][i];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k   =  bas_aa_sym[h][kl][0];
                int l   =  bas_aa_sym[h][kl][1];

                double dum  = u_p[offset + ij*gems_aa[h]+kl];

                A_p[q2toff[h] + ij*gems_aa[h]+kl]    -= dum;          // -Q2(ij,kl)

                // not spin adapted
                int kld = ibas_ab_sym[h][k][l];
                int lkd = ibas_ab_sym[h][l][k];
                A_p[d2aboff[h] + kld*gems_ab[h]+ijd] += 0.5 * dum;          // +D2(kl,ij)
                A_p[d2aboff[h] + lkd*gems_ab[h]+ijd] -= 0.5 * dum;          // +D2(lk,ij)
                A_p[d2aboff[h] + kld*gems_ab[h]+jid] -= 0.5 * dum;          // +D2(kl,ji)
                A_p[d2aboff[h] + lkd*gems_ab[h]+jid] += 0.5 * dum;          // +D2(lk,ji)

                // spin adapted
                //A_p[d2toff[h] + INDEX(kl,ij)] += dum;          // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[q1aoff[h2] + ii*amopi_[h2]+kk] += 0.5 * dum; // +Q1(i,k) djl
                    A_p[d1boff[h2] + ii*amopi_[h2]+kk] -= 0.5 * dum; // -D1(i,k) djl
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[q1aoff[h2] + ll*amopi_[h2]+jj] += 0.5 * dum; // +Q1(l,j) dik
                    A_p[d1boff[h2] + ll*amopi_[h2]+jj] -= 0.5 * dum; // -D1(l,j) dik
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[q1aoff[h2] + ll*amopi_[h2]+ii] -= 0.5 * dum; // +Q1(l,i) djk
                    A_p[d1boff[h2] + ll*amopi_[h2]+ii] += 0.5 * dum; // -D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[q1aoff[h2] + kk*amopi_[h2]+jj] -= 0.5 * dum; // +Q1(k,j) dil
                    A_p[d1boff[h2] + kk*amopi_[h2]+jj] += 0.5 * dum; // -D1(k,j) dil
                }
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
    // map D2aa to Q211
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double val = u_p[offset + ij*gems_aa[h]+kl];
                A_p[q2toff_p1[h] + ij*gems_aa[h]+kl] -= val;
                A_p[d2aaoff[h] + kl*gems_aa[h]+ij]   += val;
                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[q1aoff[h2]  + ii*amopi_[h2]+kk]      += val;
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[d1aoff[h2]  + ll*amopi_[h2]+ii]      += val;
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[q1aoff[h2]  + jj*amopi_[h2]+kk]      -= val;
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[d1aoff[h2]  + ll*amopi_[h2]+jj]      -= val;
                }
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
    // map D2bb to Q21-1
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double val = u_p[offset + ij*gems_aa[h]+kl];
                A_p[q2toff_m1[h] + ij*gems_aa[h]+kl] -= val;
                //A_p[d2toff_m1[h] + INDEX(kl,ij)] += u_p[offset + INDEX(ij,kl)];
                A_p[d2bboff[h] + kl*gems_aa[h]+ij] += val;
                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[q1boff[h2]  + ii*amopi_[h2]+kk]      += val;
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[d1boff[h2]  + ll*amopi_[h2]+ii]      += val;
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[q1boff[h2]  + jj*amopi_[h2]+kk]      -= val;
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[d1boff[h2]  + ll*amopi_[h2]+jj]      -= val;
                }
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
}

// Q2 guess
void v2RDMSolver::Q2_constraints_guess(SharedVector u){

    double * u_p = u->pointer();

    // map D2ab to Q2ab

    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];

                double dum  = 0.0;

                dum        +=  u_p[d2aboff[h] + kl*gems_ab[h]+ij];          // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[q1aoff[h2] + ii*amopi_[h2]+kk]; // +Q1(i,k) djl
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1boff[h2] + ll*amopi_[h2]+jj]; // -D1(l,j) dik
                }
                u_p[q2aboff[h] + ij*gems_ab[h]+kl] = dum;
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }

    // map D2aa to Q2aa
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];

                double dum  = 0.0;

                dum        +=  u_p[d2aaoff[h] + kl*gems_aa[h]+ij];    // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[q1aoff[h2] + ii*amopi_[h2]+kk];  // +Q1(i,k) djl
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  u_p[d1aoff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  u_p[q1aoff[h2] + jj*amopi_[h2]+kk];  // -Q1(j,k) dil
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1aoff[h2] + ll*amopi_[h2]+jj];  // -D1(l,j) dkl
                }

                u_p[q2aaoff[h] + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }

    // map D2bb to Q2bb
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];

                double dum  = 0.0;

                dum        +=  u_p[d2bboff[h] + kl*gems_aa[h]+ij];    // +D2(kl,ij)

                if ( j==l ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[q1boff[h2] + ii*amopi_[h2]+kk];  // +Q1(i,k) djl
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  u_p[d1boff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                }
                if ( i==l ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  u_p[q1boff[h2] + jj*amopi_[h2]+kk];  // -Q1(j,k) dil
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1boff[h2] + ll*amopi_[h2]+jj];  // -D1(l,j) dkl
                }
                u_p[q2bboff[h] + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
}

// Q2 portion of A.x (with symmetry)
void v2RDMSolver::Q2_constraints_Au(SharedVector A,SharedVector u){
    double * A_p = A->pointer();
    double * u_p = u->pointer();

    long int blocksize_ab = 0;
    long int blocksize_aa = 0;
    for (int h = 0; h < nirrep_; h++) {
        blocksize_ab += gems_ab[h]*gems_ab[h];
        blocksize_aa += gems_aa[h]*gems_aa[h];
    }

    // map D2ab to Q2ab
    C_DCOPY(blocksize_ab,u_p + d2aboff[0],1,A_p + offset,1);      // + D2(kl,ij)
    C_DAXPY(blocksize_ab,-1.0,u_p + q2aboff[0],1,A_p + offset,1); // - Q2(kl,ij)
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            // +Q1(i,k) djl
            //int hi = symmetry[i];
            //int ii = i - pitzer_offset[hi];
            //for (int kk = 0; kk < amopi_[hi]; kk++) {
            //    int k  = kk + pitzer_offset[hi];
            //    int kj = ibas_ab_sym[h][k][j];
            //    A_p[offset + ij*gems_ab[h]+kj] += u_p[q1aoff[hi] + ii*amopi_[hi]+kk]; // +Q1(i,k) djl
            //}
            // -D1(k,i) djl
            int hi = symmetry[i];
            int ii = i - pitzer_offset[hi];
            for (int kk = 0; kk < amopi_[hi]; kk++) {
                int k  = kk + pitzer_offset[hi];
                int kj = ibas_ab_sym[h][k][j];
                A_p[offset + ij*gems_ab[h]+kj] -= u_p[d1aoff[hi] + kk*amopi_[hi]+ii]; // +Q1(k,i) djl
            }

            // -D1(l,j) dik
            int hj = symmetry[j];
            int jj = j - pitzer_offset[hj];
            for (int ll = 0; ll < amopi_[hj]; ll++) {
                int l  = ll + pitzer_offset[hj];
                int il = ibas_ab_sym[h][i][l];
                A_p[offset + ij*gems_ab[h]+il] -= u_p[d1boff[hj] + jj*amopi_[hj]+ll]; // -D1(l,j) dik
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }

    // map D2aa to Q2aa
    C_DCOPY(blocksize_aa,u_p + d2aaoff[0],1,A_p + offset,1);      // + D2(kl,ij)
    C_DAXPY(blocksize_aa,-1.0,u_p + q2aaoff[0],1,A_p + offset,1); // - Q2(kl,ij)
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum  = 0.0;
                if ( j==l ) {
                    //int h2 = symmetry[i];
                    //int ii = i - pitzer_offset[h2];
                    //int kk = k - pitzer_offset[h2];
                    //dum        +=  u_p[q1aoff[h2] + ii*amopi_[h2]+kk];  // +Q1(i,k) djl
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  u_p[d1aoff[h2] + kk*amopi_[h2]+ii];  // -D1(k,i) djl
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  u_p[d1aoff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                }
                if ( i==l ) {
                    //int h2 = symmetry[j];
                    //int jj = j - pitzer_offset[h2];
                    //int kk = k - pitzer_offset[h2];
                    //dum        -=  u_p[q1aoff[h2] + jj*amopi_[h2]+kk];  // -Q1(j,k) dil
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[d1aoff[h2] + kk*amopi_[h2]+jj];  // +D1(k,j) dil
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1aoff[h2] + ll*amopi_[h2]+jj];  // -D1(l,j) dkl
                }
                A_p[offset + ij*gems_aa[h]+kl] += dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }


    // map D2bb to Q2bb
    C_DCOPY(blocksize_aa,u_p + d2bboff[0],1,A_p + offset,1);      // + D2(kl,ij)
    C_DAXPY(blocksize_aa,-1.0,u_p + q2bboff[0],1,A_p + offset,1); // - Q2(kl,ij)
    for (int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum  = 0.0;
                if ( j==l ) {
                    //int h2 = symmetry[i];
                    //int ii = i - pitzer_offset[h2];
                    //int kk = k - pitzer_offset[h2];
                    //dum        +=  u_p[q1boff[h2] + ii*amopi_[h2]+kk];  // +Q1(i,k) djl
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        -=  u_p[d1boff[h2] + kk*amopi_[h2]+ii];  // -D1(k,i) djl
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        +=  u_p[d1boff[h2] + ll*amopi_[h2]+ii];  // +D1(l,i) djk
                }
                if ( i==l ) {
                    //int h2 = symmetry[j];
                    //int jj = j - pitzer_offset[h2];
                    //int kk = k - pitzer_offset[h2];
                    //dum        -=  u_p[q1boff[h2] + jj*amopi_[h2]+kk];  // -Q1(j,k) dil
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    dum        +=  u_p[d1boff[h2] + kk*amopi_[h2]+jj];  // +Q1(k,j) dil
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    dum        -=  u_p[d1boff[h2] + ll*amopi_[h2]+jj];  // -D1(l,j) dkl
                }
                A_p[offset + ij*gems_aa[h]+kl] += dum;
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }
}

// Q2 portion of A^T.y (with symmetry)
void v2RDMSolver::Q2_constraints_ATu(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    long int blocksize_ab = 0;
    long int blocksize_aa = 0;
    for (int h = 0; h < nirrep_; h++) {
        blocksize_ab += gems_ab[h]*gems_ab[h];
        blocksize_aa += gems_aa[h]*gems_aa[h];
    }

    // map D2ab to Q2ab
    C_DAXPY(blocksize_ab, 1.0,u_p + offset,1,A_p + d2aboff[0],1); // + D2(kl,ij)
    C_DAXPY(blocksize_ab,-1.0,u_p + offset,1,A_p + q2aboff[0],1); // - Q2(ij,kl)
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            // +Q1(i,k) djl
            //int hi = symmetry[i];
            //int ii = i - pitzer_offset[hi];
            //for (int kk = 0; kk < amopi_[hi]; kk++) {
            //    int k  = kk + pitzer_offset[hi];
            //    int kj = ibas_ab_sym[h][k][j];
            //    A_p[q1aoff[hi] + ii*amopi_[hi]+kk] += u_p[offset + ij*gems_ab[h]+kj]; // +Q1(i,k) djl
            //}
            // -D1(k,i) djl
            int hi = symmetry[i];
            int ii = i - pitzer_offset[hi];
            for (int kk = 0; kk < amopi_[hi]; kk++) {
                int k  = kk + pitzer_offset[hi];
                int kj = ibas_ab_sym[h][k][j];
                A_p[d1aoff[hi] + kk*amopi_[hi]+ii] -= u_p[offset + ij*gems_ab[h]+kj]; // -D1(k,i) djl
            }

            // -D1(l,j) dik
            int hj = symmetry[j];
            int jj = j - pitzer_offset[hj];
            for (int ll = 0; ll < amopi_[hj]; ll++) {
                int l  = ll + pitzer_offset[hj];
                int il = ibas_ab_sym[h][i][l];
                A_p[d1boff[hj] + jj*amopi_[hj]+ll] -= u_p[offset + ij*gems_ab[h]+il]; // -D1(l,j) dik
            }
        }
        offset += gems_ab[h]*gems_ab[h];
    }

    // map D2aa to Q2aa
    C_DAXPY(blocksize_aa, 1.0,u_p + offset,1,A_p + d2aaoff[0],1); // + D2(kl,ij)
    C_DAXPY(blocksize_aa,-1.0,u_p + offset,1,A_p + q2aaoff[0],1); // - Q2(ij,kl)
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double val = u_p[offset + ij*gems_aa[h]+kl];
                if ( j==l ) {
                    //int h2 = symmetry[i];
                    //int ii = i - pitzer_offset[h2];
                    //int kk = k - pitzer_offset[h2];
                    //A_p[q1aoff[h2]  + ii*amopi_[h2]+kk] += val;
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[d1aoff[h2]  + kk*amopi_[h2]+ii] -= val;
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[d1aoff[h2]  + ll*amopi_[h2]+ii] += val;
                }
                if ( i==l ) {
                    //int h2 = symmetry[j];
                    //int jj = j - pitzer_offset[h2];
                    //int kk = k - pitzer_offset[h2];
                    //A_p[q1aoff[h2]  + jj*amopi_[h2]+kk] -= val;
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[d1aoff[h2]  + kk*amopi_[h2]+jj] += val;
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[d1aoff[h2]  + ll*amopi_[h2]+jj] -= val;
                }
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }


    // map D2bb to Q2bb
    C_DAXPY(blocksize_aa, 1.0,u_p + offset,1,A_p + d2bboff[0],1); // + D2(kl,ij)
    C_DAXPY(blocksize_aa,-1.0,u_p + offset,1,A_p + q2bboff[0],1); // - Q2(ij,kl)
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double val = u_p[offset + ij*gems_aa[h]+kl];
                if ( j==l ) {
                    //int h2 = symmetry[i];
                    //int ii = i - pitzer_offset[h2];
                    //int kk = k - pitzer_offset[h2];
                    //A_p[q1boff[h2]  + ii*amopi_[h2]+kk] += val;
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[d1boff[h2]  + kk*amopi_[h2]+ii] -= val;
                }
                if ( j==k ) {
                    int h2 = symmetry[i];
                    int ii = i - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[d1boff[h2]  + ll*amopi_[h2]+ii] += val;
                }
                if ( i==l ) {
                    //int h2 = symmetry[j];
                    //int jj = j - pitzer_offset[h2];
                    //int kk = k - pitzer_offset[h2];
                    //A_p[q1boff[h2]  + jj*amopi_[h2]+kk] -= val;
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int kk = k - pitzer_offset[h2];
                    A_p[d1boff[h2]  + kk*amopi_[h2]+jj] += val;
                }
                if ( i==k ) {
                    int h2 = symmetry[j];
                    int jj = j - pitzer_offset[h2];
                    int ll = l - pitzer_offset[h2];
                    A_p[d1boff[h2]  + ll*amopi_[h2]+jj] -= val;
                }
            }
        }
        offset += gems_aa[h]*gems_aa[h];
    }

}

}} // end namespaces
