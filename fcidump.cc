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

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/mintshelper.h>

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>

#include<psi4/libmints/wavefunction.h>
//#include<psi4/libmints/mints.h>
#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>
//#include<../bin/fnocc/blas.h>
#include<time.h>

#include"v2rdm_solver.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace psi;
//using namespace fnocc;

namespace psi{ namespace v2rdm_casscf{


// 
// dump 1-/2-electron integrals and 1-/2-RDM to disk
// 
// notes:
//
// - all quantities should be in the NO basis
// - all 1-/2-electron integrals are required
// - only active 1-/2-RDM elements are required
// - as a first pass, this will only work with DF integrals
//
void v2RDMSolver::FCIDUMP() {

    if ( nfrzv_ > 0 ) {
        throw PsiException("FCIDUMP does not work with frozen virtual orbitals",__FILE__,__LINE__);
    }

    FILE * int_fp = fopen("int.dump","wb");
    FILE * rdm_fp = fopen("rdm.dump","wb");
    //FILE * int_fp = fopen("int.txt","w");
    //FILE * rdm_fp = fopen("rdm.txt","w");

    int zero = 0;

    // two-electron integrals
    for (int p = 0; p < nmo_; p++) {
        for (int q = p; q < nmo_; q++) {
            long int pq = INDEX(p,q);
            for (int r = 0; r < nmo_; r++) {
                for (int s = r; s < nmo_; s++) {
                    long int rs = INDEX(r,s);
                    if ( pq > rs ) continue;
                    double dum = TEI(p,q,r,s,0);
                    //fprintf(int_fp,"%20.12lf %5i %5i %5i %5i\n",dum,p+1,q+1,r+1,s+1);
                    int pp = p+1;
                    int qq = q+1;
                    int rr = r+1;
                    int ss = s+1;
                    fwrite (&dum , sizeof(double), 1, int_fp);
                    fwrite (&pp , sizeof(int), 1, int_fp);
                    fwrite (&qq , sizeof(int), 1, int_fp);
                    fwrite (&rr , sizeof(int), 1, int_fp);
                    fwrite (&ss , sizeof(int), 1, int_fp);
                }
            }
        }
    }
    // one-electron integrals

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    std::shared_ptr<Matrix> T (new Matrix(mints->so_kinetic()));
    std::shared_ptr<Matrix> V (new Matrix(mints->so_potential()));

    T->transform(Ca_);
    V->transform(Ca_);

    for (int h = 0; h < nirrep_; h++) {
        double ** Tp = T_->pointer(h);
        double ** Vp = V_->pointer(h);
        for (int p = 0; p < nmopi_[h]; p++) {
            for (int q = p; q < nmopi_[h]; q++) {
                double dum = Tp[p][q] + Vp[p][q];
                //fprintf(int_fp,"%20.12lf %5i %5i %5i %5i\n",dum,p + pitzer_offset[h]+1,q + pitzer_offset[h]+1,0,0);
                int pp = p+1;
                int qq = q+1;
                fwrite (&dum , sizeof(double), 1, int_fp);
                fwrite (&pp , sizeof(int), 1, int_fp);
                fwrite (&qq , sizeof(int), 1, int_fp);
                fwrite (&zero , sizeof(int), 1, int_fp);
                fwrite (&zero , sizeof(int), 1, int_fp);
            }
        }
    }

    // two-electron RDM
    double * x_p = x->pointer();
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                double dum = 2.0 * x_p[d2aboff[h] + ij * gems_ab[h] + kl];
                if ( i != j && k != l ) {

                    int ija = ibas_aa_sym[h][i][j];
                    int kla = ibas_aa_sym[h][k][l];
                    int sg = 1;
                    if ( i > j ) sg = -sg;
                    if ( k > l ) sg = -sg;

                    dum += sg * x_p[d2aaoff[h] + ija * gems_aa[h] + kla];
                    dum += sg * x_p[d2bboff[h] + ija * gems_aa[h] + kla];

                }
                //fprintf(rdm_fp,"%20.12lf %5i %5i %5i %5i\n",dum,i+1,j+1,k+1,l+1);
                int ii = i+1;
                int jj = j+1;
                int kk = k+1;
                int ll = l+1;
                fwrite (&dum , sizeof(double), 1, rdm_fp);
                fwrite (&ii , sizeof(int), 1, rdm_fp);
                fwrite (&jj , sizeof(int), 1, rdm_fp);
                fwrite (&kk , sizeof(int), 1, rdm_fp);
                fwrite (&ll , sizeof(int), 1, rdm_fp);
            }
        }
    }

    // one-electron integrals
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = i; j < amopi_[h]; j++) {
                double dum = x_p[d1aoff[h] + i * amopi_[h] + j] + x_p[d1boff[h] + i * amopi_[h] + j];
                //fprintf(rdm_fp,"%20.12lf %5i %5i %5i %5i\n",dum,i + pitzer_offset[h]+1,j + pitzer_offset[h]+1,0,0);
                int ii = i+1;
                int jj = j+1;
                fwrite (&dum , sizeof(double), 1, rdm_fp);
                fwrite (&ii , sizeof(int), 1, rdm_fp);
                fwrite (&jj , sizeof(int), 1, rdm_fp);
                fwrite (&zero , sizeof(int), 1, rdm_fp);
                fwrite (&zero , sizeof(int), 1, rdm_fp);
            }
        }
    }

    fclose(int_fp);
    fclose(rdm_fp);

}


}}
