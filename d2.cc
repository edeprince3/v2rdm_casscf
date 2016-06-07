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


// D2 portion of A^T.y ( and D1 / Q1 )
void v2RDMSolver::D2_constraints_ATu(SharedVector A,SharedVector u){
    double* A_p = A->pointer();
    double* u_p = u->pointer();

    if ( constrain_spin_ ) {
        // spin
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
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
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            int ij = ibas_ab_sym[h][i][j];
            if ( gems_ab[h] == 0 ) continue;
            A_p[d2aboff[h] + ij*gems_ab[h]+ij] += u_p[offset];
        }
    }
    offset++;

    // Tr(D2aa)
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
            if ( i==j ) continue;
            int h = SymmetryPair(symmetry[i],symmetry[j]);
            if ( gems_aa[h] == 0 ) continue;
            int ij = ibas_aa_sym[h][i][j];
            A_p[d2aaoff[h]+ij*gems_aa[h]+ij] += u_p[offset];
        }
    }
    offset++;
    // Tr(D2bb)
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
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

    int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;

    int poff = 0;

    // contraction: D2ab -> D1 a
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                A_p[d1aoff[h] + i*amopi_[h]+j] += nb * u_p[offset + i*amopi_[h]+j];
                int ii = i + poff;
                int jj = j + poff;
                for (int k = 0; k < amo_; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][ii][k];
                    int jk = ibas_ab_sym[h2][jj][k];
                    A_p[d2aboff[h2] + ik*gems_ab[h2]+jk] -= u_p[offset + i*amopi_[h]+j];
                }
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    // contraction: D2ab -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A_p[d1boff[h] + i*amopi_[h]+j] += na * u_p[offset + i*amopi_[h]+j];
                int ii = i + poff;
                int jj = j + poff;
                for(int k = 0; k < amo_; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][k][ii];
                    int jk = ibas_ab_sym[h2][k][jj];
                    A_p[d2aboff[h2] + ik*gems_ab[h2]+jk] -= u_p[offset + i*amopi_[h]+j];
                }
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    //contract D2aa -> D1 a
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A_p[d1aoff[h] + i*amopi_[h]+j] += (na - 1.0) * u_p[offset + i*amopi_[h]+j];
                int ii = i + poff;
                int jj = j + poff;
                for(int k =0; k < amo_; k++){
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
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    //contract D2bb -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for(int i = 0; i < amopi_[h]; i++){
            for(int j = 0; j < amopi_[h]; j++){
                A_p[d1boff[h] + i*amopi_[h]+j] += (nb - 1.0) * u_p[offset + i*amopi_[h]+j];
                int ii = i + poff;
                int jj = j + poff;
                for(int k =0; k < amo_; k++){
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
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }


    // additional spin constraints for singlets:
    if ( constrain_spin_ && nalpha_ == nbeta_ ) {
        // D1a = D1b
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(amopi_[h]*amopi_[h], 1.0, u_p + offset, 1, A_p + d1aoff[h],1);
            C_DAXPY(amopi_[h]*amopi_[h],-1.0, u_p + offset, 1, A_p + d1boff[h],1);
            offset += amopi_[h]*amopi_[h];
        }
        // D2aa = D2bb
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_aa[h]*gems_aa[h], 1.0, u_p + offset, 1, A_p + d2aaoff[h],1);
            C_DAXPY(gems_aa[h]*gems_aa[h],-1.0, u_p + offset, 1, A_p + d2bboff[h],1);
            offset += gems_aa[h]*gems_aa[h];
        }
        // D2aa[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_aa[h]*gems_aa[h],1.0,u_p + offset,1,A_p + d2aaoff[h],1);
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0]; 
                int j = bas_aa_sym[h][ij][1];
                int ijb = ibas_ab_sym[h][i][j];
                int jib = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0]; 
                    int l = bas_aa_sym[h][kl][1];
                    int klb = ibas_ab_sym[h][k][l];
                    int lkb = ibas_ab_sym[h][l][k];
                    A_p[d2aboff[h] + ijb*gems_ab[h] + klb] -= 0.5 * u_p[offset + ij*gems_aa[h] + kl];
                    A_p[d2aboff[h] + jib*gems_ab[h] + klb] += 0.5 * u_p[offset + ij*gems_aa[h] + kl];
                    A_p[d2aboff[h] + ijb*gems_ab[h] + lkb] += 0.5 * u_p[offset + ij*gems_aa[h] + kl];
                    A_p[d2aboff[h] + jib*gems_ab[h] + lkb] -= 0.5 * u_p[offset + ij*gems_aa[h] + kl];
                }   
            }   
            offset += gems_aa[h]*gems_aa[h];
        }   
        // D2bb[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_aa[h]*gems_aa[h],1.0,u_p + offset,1,A_p + d2bboff[h],1);
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                int ijb = ibas_ab_sym[h][i][j];
                int jib = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    int klb = ibas_ab_sym[h][k][l];
                    int lkb = ibas_ab_sym[h][l][k];
                    A_p[d2aboff[h] + ijb*gems_ab[h] + klb] -= 0.5 * u_p[offset + ij*gems_aa[h] + kl];
                    A_p[d2aboff[h] + jib*gems_ab[h] + klb] += 0.5 * u_p[offset + ij*gems_aa[h] + kl];
                    A_p[d2aboff[h] + ijb*gems_ab[h] + lkb] += 0.5 * u_p[offset + ij*gems_aa[h] + kl];
                    A_p[d2aboff[h] + jib*gems_ab[h] + lkb] -= 0.5 * u_p[offset + ij*gems_aa[h] + kl];
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D200 = 1/(2 sqrt(1+dpq)sqrt(1+drs)) ( D2ab[pq][rs] + D2ab[pq][sr] + D2ab[qp][rs] + D2ab[qp][sr] )
        for ( int h = 0; h < nirrep_; h++) {
            C_DAXPY(gems_ab[h]*gems_ab[h],1.0,u_p + offset,1,A_p + d200off[h],1);
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A_p[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[offset + ij*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ji*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[offset + ij*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ij*gems_ab[h] + lk] -= 0.5 / ( dij * dkl ) * u_p[offset + ij*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ji*gems_ab[h] + lk] -= 0.5 / ( dij * dkl ) * u_p[offset + ij*gems_ab[h] + kl];
                }
            }
            offset += gems_ab[h]*gems_ab[h];
        }
    }else if ( constrain_spin_ ) { // nonsinglets ... big block

        for ( int h = 0; h < nirrep_; h++) {
            // D200
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A_p[d200off[h] + ij*2*gems_ab[h] + kl] += u_p[offset + ij*2*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[offset + ij*2*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ji*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[offset + ij*2*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ij*gems_ab[h] + lk] -= 0.5 / ( dij * dkl ) * u_p[offset + ij*2*gems_ab[h] + kl];
                    A_p[d2aboff[h] + ji*gems_ab[h] + lk] -= 0.5 / ( dij * dkl ) * u_p[offset + ij*2*gems_ab[h] + kl];
                }
            }
            // D201
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A_p[d200off[h] + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] += u_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 / dij * u_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[d2aboff[h] + ij*gems_ab[h] + lk] += 0.5 / dij * u_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[d2aboff[h] + ji*gems_ab[h] + kl] -= 0.5 / dij * u_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[d2aboff[h] + ji*gems_ab[h] + lk] += 0.5 / dij * u_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                }
            }
            // D210
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A_p[d200off[h] + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] += u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A_p[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 / dkl * u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A_p[d2aboff[h] + ij*gems_ab[h] + lk] -= 0.5 / dkl * u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A_p[d2aboff[h] + ji*gems_ab[h] + kl] += 0.5 / dkl * u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A_p[d2aboff[h] + ji*gems_ab[h] + lk] += 0.5 / dkl * u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                }
            }
            // D211
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A_p[d200off[h] + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] += u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[d2aboff[h] + ij*gems_ab[h] + kl] -= 0.5 * u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[d2aboff[h] + ji*gems_ab[h] + kl] += 0.5 * u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[d2aboff[h] + ij*gems_ab[h] + lk] += 0.5 * u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[d2aboff[h] + ji*gems_ab[h] + lk] -= 0.5 * u_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                }
            }
            offset += 4*gems_ab[h]*gems_ab[h];
        }
    }

}

// D2 portion of A.x (and D1/Q1)
void v2RDMSolver::D2_constraints_Au(SharedVector A,SharedVector u){

    double* A_p = A->pointer();
    double* u_p = u->pointer();

    if ( constrain_spin_ ) {
        // spin
        double s2 = 0.0;
        for (int i = 0; i < amo_; i++){
            for (int j = 0; j < amo_; j++){
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

    // Traces
    // Tr(D2ab)
    double sumab =0.0;
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
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
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
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
    for (int i = 0; i < amo_; i++){
        for (int j = 0; j < amo_; j++){
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

    int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;
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
                for(int k = 0; k < amo_; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][ii][k];
                    int jk = ibas_ab_sym[h2][jj][k];
                    sum -= u_p[d2aboff[h2] + ik*gems_ab[h2]+jk];
                }
                A_p[offset + i*amopi_[h]+j] = sum;
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    // contraction: D2ab -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = na * u_p[d1boff[h] + i*amopi_[h]+j];
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < amo_; k++){
                    int h2  = SymmetryPair(symmetry[ii],symmetry[k]);
                    int ik = ibas_ab_sym[h2][k][ii];
                    int jk = ibas_ab_sym[h2][k][jj];
                    sum -= u_p[d2aboff[h2] + ik*gems_ab[h2]+jk];
                }
                A_p[offset + i*amopi_[h]+j] = sum;
            }
        }
        offset += amopi_[h]*amopi_[h];
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    //contract D2aa -> D1 a
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = (na - 1.0) * u_p[d1aoff[h] + i*amopi_[h]+j];
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < amo_; k++){
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
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    //contract D2bb -> D1 b
    poff = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++){
            for (int j = 0; j < amopi_[h]; j++){
                double sum = (nb - 1.0) * u_p[d1boff[h] + i*amopi_[h]+j];
                int ii  = i + poff;
                int jj  = j + poff;
                for(int k = 0; k < amo_; k++){
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
        poff   += nmopi_[h] - rstcpi_[h] - frzcpi_[h] - rstvpi_[h] - frzvpi_[h];
    }

    // additional spin constraints for singlets:
    if ( constrain_spin_ && nalpha_ == nbeta_ ) {
        // D1a = D1b
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(amopi_[h]*amopi_[h],     u_p + d1aoff[h],1,A_p + offset,1);
            C_DAXPY(amopi_[h]*amopi_[h],-1.0,u_p + d1boff[h],1,A_p + offset,1);
            offset += amopi_[h]*amopi_[h]; 
        }
        // D2aa = D2bb
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(gems_aa[h]*gems_aa[h],     u_p + d2aaoff[h],1,A_p + offset,1);
            C_DAXPY(gems_aa[h]*gems_aa[h],-1.0,u_p + d2bboff[h],1,A_p + offset,1);
            offset += gems_aa[h]*gems_aa[h];
        }
        // D2aa[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(gems_aa[h]*gems_aa[h],u_p + d2aaoff[h],1,A_p + offset,1);
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                int ijb = ibas_ab_sym[h][i][j];
                int jib = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    int klb = ibas_ab_sym[h][k][l];
                    int lkb = ibas_ab_sym[h][l][k];
                    A_p[offset + ij*gems_aa[h] + kl] -= 0.5 * u_p[d2aboff[h] + ijb*gems_ab[h] + klb];
                    A_p[offset + ij*gems_aa[h] + kl] += 0.5 * u_p[d2aboff[h] + jib*gems_ab[h] + klb];
                    A_p[offset + ij*gems_aa[h] + kl] += 0.5 * u_p[d2aboff[h] + ijb*gems_ab[h] + lkb];
                    A_p[offset + ij*gems_aa[h] + kl] -= 0.5 * u_p[d2aboff[h] + jib*gems_ab[h] + lkb];
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D2bb[pq][rs] = 1/2(D2ab[pq][rs] - D2ab[pq][sr] - D2ab[qp][rs] + D2ab[qp][sr])
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(gems_aa[h]*gems_aa[h],u_p + d2bboff[h],1,A_p + offset,1);
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                int ijb = ibas_ab_sym[h][i][j];
                int jib = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    int klb = ibas_ab_sym[h][k][l];
                    int lkb = ibas_ab_sym[h][l][k];
                    A_p[offset + ij*gems_aa[h] + kl] -= 0.5 * u_p[d2aboff[h] + ijb*gems_ab[h] + klb];
                    A_p[offset + ij*gems_aa[h] + kl] += 0.5 * u_p[d2aboff[h] + jib*gems_ab[h] + klb];
                    A_p[offset + ij*gems_aa[h] + kl] += 0.5 * u_p[d2aboff[h] + ijb*gems_ab[h] + lkb];
                    A_p[offset + ij*gems_aa[h] + kl] -= 0.5 * u_p[d2aboff[h] + jib*gems_ab[h] + lkb];
                }
            }
            offset += gems_aa[h]*gems_aa[h];
        }
        // D200 = 1/(2 sqrt(1+dpq)sqrt(1+drs)) ( D2ab[pq][rs] + D2ab[pq][sr] + D2ab[qp][rs] + D2ab[qp][sr] )
        for ( int h = 0; h < nirrep_; h++) {
            C_DCOPY(gems_ab[h]*gems_ab[h],u_p + d200off[h],1,A_p + offset,1);
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A_p[offset + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                    A_p[offset + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[d2aboff[h] + ji*gems_ab[h] + kl];
                    A_p[offset + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[d2aboff[h] + ij*gems_ab[h] + lk];
                    A_p[offset + ij*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            offset += gems_ab[h]*gems_ab[h];
        }
    }else if ( constrain_spin_ ) { // nonsinglets ... big block

        for ( int h = 0; h < nirrep_; h++) {
            // D200
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A_p[offset + ij*2*gems_ab[h] + kl] += u_p[d200off[h] + ij*2*gems_ab[h] + kl];
                    A_p[offset + ij*2*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                    A_p[offset + ij*2*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[d2aboff[h] + ji*gems_ab[h] + kl];
                    A_p[offset + ij*2*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[d2aboff[h] + ij*gems_ab[h] + lk];
                    A_p[offset + ij*2*gems_ab[h] + kl] -= 0.5 / ( dij * dkl ) * u_p[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            // D201
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                double dij = ( i == j ) ? sqrt(2.0) : 1.0;
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] += u_p[d200off[h] + (ij)*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] -= 0.5 / dij * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                    A_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] += 0.5 / dij * u_p[d2aboff[h] + ij*gems_ab[h] + lk];
                    A_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] -= 0.5 / dij * u_p[d2aboff[h] + ji*gems_ab[h] + kl];
                    A_p[offset + (ij)*2*gems_ab[h] + (kl+gems_ab[h])] += 0.5 / dij * u_p[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            // D210
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    double dkl = ( k == l ) ? sqrt(2.0) : 1.0;
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] += u_p[d200off[h] + (ij+gems_ab[h])*2*gems_ab[h] + (kl)];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] -= 0.5 / dkl * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] -= 0.5 / dkl * u_p[d2aboff[h] + ij*gems_ab[h] + lk];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] += 0.5 / dkl * u_p[d2aboff[h] + ji*gems_ab[h] + kl];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl)] += 0.5 / dkl * u_p[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            // D211
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                int ji = ibas_ab_sym[h][j][i];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    int lk = ibas_ab_sym[h][l][k];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] += u_p[d200off[h] + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] -= 0.5 * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] += 0.5 * u_p[d2aboff[h] + ji*gems_ab[h] + kl];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] += 0.5 * u_p[d2aboff[h] + ij*gems_ab[h] + lk];
                    A_p[offset + (ij+gems_ab[h])*2*gems_ab[h] + (kl+gems_ab[h])] -= 0.5 * u_p[d2aboff[h] + ji*gems_ab[h] + lk];
                }
            }
            offset += 4*gems_ab[h]*gems_ab[h];
        }
    }

}

}}
