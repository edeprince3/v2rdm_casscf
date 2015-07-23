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

// D3 portion of A.u 
void v2RDMSolver::D3_constraints_Au(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    int na = nalpha_ - nfrzc;
    int nb = nbeta_ - nfrzc;

    // D3aaa -> D2aa
    for ( int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = (na - 2.0) * u_p[d2aaoff[h] + ij*gems_aa[h] + kl];
                for ( int p = 0; p < nmo; p++) {
                    if ( i == p || j == p ) continue;
                    if ( k == p || l == p ) continue;
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aaa_sym[h2][i][j][p];
                    int klp = ibas_aaa_sym[h2][k][l][p];
                    int s = 1;
                    if ( p < i ) s = -s;
                    if ( p < j ) s = -s;
                    if ( p < k ) s = -s;
                    if ( p < l ) s = -s;
                    dum -= s * u_p[d3aaaoff[h2] + ijp*trip_aaa[h2]+klp];
                }
                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3bbb -> D2bb
    for ( int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = (nb - 2.0) * u_p[d2bboff[h] + ij*gems_aa[h] + kl];
                for ( int p = 0; p < nmo; p++) {
                    if ( i == p || j == p ) continue;
                    if ( k == p || l == p ) continue;
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aaa_sym[h2][i][j][p];
                    int klp = ibas_aaa_sym[h2][k][l][p];
                    int s = 1;
                    if ( p < i ) s = -s;
                    if ( p < j ) s = -s;
                    if ( p < k ) s = -s;
                    if ( p < l ) s = -s;
                    dum -= s * u_p[d3bbboff[h2] + ijp*trip_aaa[h2]+klp];
                }
                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3aab -> D2aa
    for ( int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = nb * u_p[d2aaoff[h] + ij*gems_aa[h] + kl];
                for ( int p = 0; p < nmo; p++) {
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][j][p];
                    int klp = ibas_aab_sym[h2][k][l][p];
                    dum -= u_p[d3aaboff[h2] + ijp*trip_aab[h2]+klp];
                }
                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3bba -> D2bb
    for ( int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = na * u_p[d2bboff[h] + ij*gems_aa[h] + kl];
                for ( int p = 0; p < nmo; p++) {
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][j][p];
                    int klp = ibas_aab_sym[h2][k][l][p];
                    dum -= u_p[d3bbaoff[h2] + ijp*trip_aab[h2]+klp];
                }
                A_p[offset + ij*gems_aa[h]+kl] = dum;
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3aab -> D2ab
    for ( int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for ( int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            for ( int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = (na - 1.0) * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                for ( int p = 0; p < nmo; p++) {
                    if ( i == p) continue;
                    if ( k == p) continue;
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][p][j];
                    int klp = ibas_aab_sym[h2][k][p][l];
                    int s = 1;
                    if ( p < i ) s = -s;
                    if ( p < k ) s = -s;
                    dum -= s * u_p[d3aaboff[h2] + ijp*trip_aab[h2]+klp];
                }
                A_p[offset + ij*gems_ab[h]+kl] = dum;
            }
        }
        offset += gems_ab[h] * gems_ab[h];
    }
    // D3bba -> D2ab
    for ( int h = 0; h < nirrep_; h++) {
        #pragma omp parallel for schedule (static)
        for ( int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            for ( int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = (nb - 1.0) * u_p[d2aboff[h] + ij*gems_ab[h] + kl];
                for ( int p = 0; p < nmo; p++) {
                    if ( j == p) continue;
                    if ( l == p) continue;
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][j][p][i];
                    int klp = ibas_aab_sym[h2][l][p][k];
                    int s = 1;
                    if ( p < j ) s = -s;
                    if ( p < l ) s = -s;
                    dum -= s * u_p[d3bbaoff[h2] + ijp*trip_aab[h2]+klp];
                }
                A_p[offset + ij*gems_ab[h]+kl] = dum;
            }
        }
        offset += gems_ab[h] * gems_ab[h];
    }

}

// D3 portion of A^T.y 
void v2RDMSolver::D3_constraints_ATu(SharedVector A,SharedVector u){

    double * A_p = A->pointer();
    double * u_p = u->pointer();

    int na = nalpha_ - nfrzc;
    int nb = nbeta_ - nfrzc;

    // D3aaa -> D2aa
    for ( int h = 0; h < nirrep_; h++) {
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = u_p[offset + ij*gems_aa[h] + kl];
                A_p[d2aaoff[h] + ij*gems_aa[h] + kl] += (na - 2.0) * dum;
                for ( int p = 0; p < nmo; p++) {
                    if ( i == p || j == p ) continue;
                    if ( k == p || l == p ) continue;
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aaa_sym[h2][i][j][p];
                    int klp = ibas_aaa_sym[h2][k][l][p];
                    int s = 1;
                    if ( p < i ) s = -s;
                    if ( p < j ) s = -s;
                    if ( p < k ) s = -s;
                    if ( p < l ) s = -s;
                    A_p[d3aaaoff[h2] + ijp*trip_aaa[h2]+klp] -= s * dum;
                }
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3bbb -> D2bb
    for ( int h = 0; h < nirrep_; h++) {
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = u_p[offset + ij*gems_aa[h] + kl];
                A_p[d2bboff[h] + ij*gems_aa[h] + kl] += (nb - 2.0) * dum;
                for ( int p = 0; p < nmo; p++) {
                    if ( i == p || j == p ) continue;
                    if ( k == p || l == p ) continue;
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aaa_sym[h2][i][j][p];
                    int klp = ibas_aaa_sym[h2][k][l][p];
                    int s = 1;
                    if ( p < i ) s = -s;
                    if ( p < j ) s = -s;
                    if ( p < k ) s = -s;
                    if ( p < l ) s = -s;
                    A_p[d3bbboff[h2] + ijp*trip_aaa[h2]+klp] -= s * dum;
                }
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3aab -> D2aa
    for ( int h = 0; h < nirrep_; h++) {
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = u_p[offset + ij*gems_aa[h] + kl];
                A_p[d2aaoff[h] + ij*gems_aa[h] + kl] += nb * dum;
                for ( int p = 0; p < nmo; p++) {
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][j][p];
                    int klp = ibas_aab_sym[h2][k][l][p];
                    A_p[d3aaboff[h2] + ijp*trip_aab[h2]+klp] -= dum;
                }
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3bba -> D2bb
    for ( int h = 0; h < nirrep_; h++) {
        for ( int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for ( int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                double dum = u_p[offset + ij*gems_aa[h] + kl];
                A_p[d2bboff[h] + ij*gems_aa[h] + kl] += na * dum;
                for ( int p = 0; p < nmo; p++) {
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][j][p];
                    int klp = ibas_aab_sym[h2][k][l][p];
                    A_p[d3bbaoff[h2] + ijp*trip_aab[h2]+klp] -= dum;
                }
            }
        }
        offset += gems_aa[h] * gems_aa[h];
    }
    // D3aab -> D2ab
    for ( int h = 0; h < nirrep_; h++) {
        for ( int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            for ( int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = u_p[offset + ij*gems_ab[h] + kl];
                A_p[d2aboff[h] + ij*gems_ab[h] + kl] += (na - 1.0) * dum;
                for ( int p = 0; p < nmo; p++) {
                    if ( i == p) continue;
                    if ( k == p) continue;
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][i][p][j];
                    int klp = ibas_aab_sym[h2][k][p][l];
                    int s = 1;
                    if ( p < i ) s = -s;
                    if ( p < k ) s = -s;
                    A_p[d3aaboff[h2] + ijp*trip_aab[h2]+klp] -= s * dum;
                }
            }
        }
        offset += gems_ab[h] * gems_ab[h];
    }
    // D3bba -> D2ab
    for ( int h = 0; h < nirrep_; h++) {
        for ( int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            for ( int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = u_p[offset + ij*gems_ab[h] + kl];
                A_p[d2aboff[h] + ij*gems_ab[h] + kl] += (nb - 1.0) * dum;
                for ( int p = 0; p < nmo; p++) {
                    if ( j == p) continue;
                    if ( l == p) continue;
                    int h2 = SymmetryPair(h,symmetry[p]);
                    int ijp = ibas_aab_sym[h2][j][p][i];
                    int klp = ibas_aab_sym[h2][l][p][k];
                    int s = 1;
                    if ( p < j ) s = -s;
                    if ( p < l ) s = -s;
                    A_p[d3bbaoff[h2] + ijp*trip_aab[h2]+klp] -= s * dum;
                }
            }
        }
        offset += gems_ab[h] * gems_ab[h];
    }
}

}} // end namespaces
