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
    FILE * int_fp_txt = fopen("int.txt","w");
    FILE * rdm_fp_txt = fopen("rdm.txt","w");

    int zero = 0;

    // orbitals are ordered by irrep and by space within each irrep. 
    // need a map to order by space and by irrep within each space
    // (plus 1)
    int * map = (int*)malloc(nmo_*sizeof(int));
    int count = 0;
    // core
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
            //printf("core    %5i -> %5i\n",i + pitzer_offset_full[h] + frzcpi_[h],count);
            map[i + pitzer_offset_full[h]] = count + 1;
            count++;
        }
    }
    // active
    for (int h = 0; h < nirrep_; h++) {
        for (int t = 0; t < amopi_[h]; t++) {
            //printf("active  %5i -> %5i\n",t + pitzer_offset_full[h] + frzcpi_[h] + rstcpi_[h],count);
            map[t + pitzer_offset_full[h] + frzcpi_[h] + rstcpi_[h]] = count + 1;
            count++;
        }
    }
    // virtual
    for (int h = 0; h < nirrep_; h++) {
        for (int a = 0; a < nmopi_[h] - ( frzcpi_[h] + rstcpi_[h] + amopi_[h] ); a++) {
            //printf("virtual %5i -> %5i\n",a + pitzer_offset_full[h] + frzcpi_[h] + rstcpi_[h] + amopi_[h],count);
            map[a + pitzer_offset_full[h] + frzcpi_[h] + rstcpi_[h] + amopi_[h]] = count + 1;
            count++;
        }
    }
    //for (int p = 0; p < nmo_; p++) {
    //    printf("%5i %5i %5i\n",symmetry_full[p],p,map[p]);
    //}

    // two-electron integrals
    for (int p = 0; p < nmo_; p++) {
        int hp = symmetry_really_full[p];
        for (int q = p; q < nmo_; q++) {
            int hq = symmetry_really_full[q];
            int hpq = hp ^ hq;
            long int pq = INDEX(p,q);
            for (int r = 0; r < nmo_; r++) {
                int hr = symmetry_really_full[r];
                for (int s = r; s < nmo_; s++) {
                    int hs = symmetry_really_full[s];
                    int hrs = hr ^ hs;
                    if ( hpq != hrs ) continue;
                    long int rs = INDEX(r,s);
                    if ( pq > rs ) continue;
                    double dum = TEI(p,q,r,s,0);
                    //if ( fabs(dum) < 1e-12 ) continue;
                    fprintf(int_fp_txt,"%20.12lf %5i %5i %5i %5i\n",dum,map[p],map[q],map[r],map[s]);
                    int pp = map[p];
                    int qq = map[q];
                    int rr = map[r];
                    int ss = map[s];
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
            for (int q = 0; q < nmopi_[h]; q++) {
                double dum = Tp[p][q] + Vp[p][q];
                //if ( fabs(dum) < 1e-12 ) continue;
                int pp = map[p + pitzer_offset_full[h]];
                int qq = map[q + pitzer_offset_full[h]];
                fprintf(int_fp_txt,"%20.12lf %5i %5i %5i %5i\n",dum,pp,qq,0,0);
                fwrite (&dum , sizeof(double), 1, int_fp);
                fwrite (&pp , sizeof(int), 1, int_fp);
                fwrite (&qq , sizeof(int), 1, int_fp);
                fwrite (&zero , sizeof(int), 1, int_fp);
                fwrite (&zero , sizeof(int), 1, int_fp);
            }
        }
    }

    // two-electron rdm
    double * x_p = x->pointer();
    double e2 = 0.0;
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

                e2 += 0.5 * dum * TEI(full_basis[i],full_basis[k],full_basis[j],full_basis[l],0);

                //if ( fabs(dum) < 1e-12 ) continue;
                int ii = map[full_basis[i]];
                int jj = map[full_basis[j]];
                int kk = map[full_basis[k]];
                int ll = map[full_basis[l]];
                fprintf(rdm_fp_txt,"%20.12lf %5i %5i %5i %5i\n",dum,ii,jj,kk,ll);
                fwrite (&dum , sizeof(double), 1, rdm_fp);
                fwrite (&ii , sizeof(int), 1, rdm_fp);
                fwrite (&jj , sizeof(int), 1, rdm_fp);
                fwrite (&kk , sizeof(int), 1, rdm_fp);
                fwrite (&ll , sizeof(int), 1, rdm_fp);
            }
        }
    }
    //printf("%20.12lf\n",e2);

    // one-electron rdm
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            int ii = i - pitzer_offset[h];
            int jj = j - pitzer_offset[h];
            double dum = x_p[d1aoff[h] + ii * amopi_[h] + jj] + x_p[d1boff[h] + ii * amopi_[h] + jj];

            ii = map[full_basis[i]];
            jj = map[full_basis[j]];
            fprintf(rdm_fp_txt,"%20.12lf %5i %5i %5i %5i\n",dum,ii,jj,0,0);
            fwrite (&dum , sizeof(double), 1, rdm_fp);
            fwrite (&ii , sizeof(int), 1, rdm_fp);
            fwrite (&jj , sizeof(int), 1, rdm_fp);
            fwrite (&zero , sizeof(int), 1, rdm_fp);
            fwrite (&zero , sizeof(int), 1, rdm_fp);

        }
    }

    fclose(int_fp);
    fclose(rdm_fp);
    fclose(int_fp_txt);
    fclose(rdm_fp_txt);

}


}}
