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

void v2RDMSolver::TEI() {

    
    // one-electron integrals:  
    boost::shared_ptr<Matrix> K1 = GetOEI();

    // size of the tei buffer
    if ( is_df_ ) {

        // size of the 3-index integral buffer
        tei_full_dim_ = (long int) nQ_ * (long int) ( nmo_ - nfrzv_ ) * ( (long int) ( nmo_ - nfrzv_ ) + 1L ) / 2L ;

        // just point to 3-index integral buffer
        tei_full_sym_      = Qmo_;

    }else {

        // size of the 4-index integral buffer
        tei_full_dim_ = 0;
        for (int h = 0; h < nirrep_; h++) {
            tei_full_dim_ += (long int)gems_full[h] * ( (long int)gems_full[h] + 1L ) / 2L;
        }

        tei_full_sym_ = (double*)malloc(tei_full_dim_*sizeof(double));
        memset((void*)tei_full_sym_,'\0',tei_full_dim_*sizeof(double));

    }

    // size of d2, blocked by symmetry, including the core orbitals
    d2_plus_core_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        d2_plus_core_dim_ += (long int)gems_plus_core[h] * ( (long int)gems_plus_core[h] + 1L ) / 2L;
    }

    d2_plus_core_sym_  = (double*)malloc(d2_plus_core_dim_*sizeof(double));
    memset((void*)d2_plus_core_sym_,'\0',d2_plus_core_dim_*sizeof(double));


    // allocate memory for oei tensor, blocked by symmetry, excluding frozen virtuals
    oei_full_dim_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        oei_full_dim_ += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // allocate memory for d1 tensor, blocked by symmetry, including the core orbitals
    d1_plus_core_dim_ = 0;
    for ( int h = 0; h < nirrep_; h++) {
        d1_plus_core_dim_ += (rstcpi_[h] + frzcpi_[h] + amopi_[h]) * ( rstcpi_[h] + frzcpi_[h] + amopi_[h] + 1 ) / 2;
    }

    oei_full_sym_ = (double*)malloc(oei_full_dim_*sizeof(double));
    memset((void*)oei_full_sym_,'\0',oei_full_dim_*sizeof(double));

    d1_plus_core_sym_ = (double*)malloc(d1_plus_core_dim_*sizeof(double));
    memset((void*)d1_plus_core_sym_,'\0',d1_plus_core_dim_*sizeof(double));

    offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                oei_full_sym_[offset + INDEX(i,j)] = K1->pointer(h)[i][j];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    if ( is_df_ ) {
        // build tei's from 3-index integrals 
        RepackIntegralsDF();
    }else {
        // read tei's from disk
        TEIFromDisk();
    }

}

void v2RDMSolver::TEIFromDisk() {

    double * temptei = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)temptei,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // read two-electron integrals from disk
    ReadIntegrals(temptei,nmo_);
  
    // load tei_full_sym_
    long int offset = 0;
    long int n2 = (long int)nmo_*(long int)nmo_;
    long int n3 = n2 * (long int)nmo_;
    for (int h = 0; h < nirrep_; h++) {
        for (long int ij = 0; ij < gems_full[h]; ij++) {
            long int i = bas_full_sym[h][ij][0];
            long int j = bas_full_sym[h][ij][1];
            for (long int kl = ij; kl < gems_full[h]; kl++) {
                long int k = bas_full_sym[h][kl][0];
                long int l = bas_full_sym[h][kl][1];
                tei_full_sym_[offset + INDEX(ij,kl)] = temptei[i*n3+j*n2+k*(long int)nmo_+l];
            }
        }
        offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
    }

    free(temptei);

    RepackIntegrals();
}

// repack rotated full-space integrals into active-space integrals
void v2RDMSolver::RepackIntegralsDF(){

    FrozenCoreEnergy();

    double * c_p = c->pointer();

    // two-electron part
    long int na = nalpha_ - nrstc_ - nfrzc_;
    long int nb = nbeta_ - nrstc_ - nfrzc_;
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

    FrozenCoreEnergy();

    double* c_p = c->pointer();

    int na = nalpha_ - nrstc_ - nfrzc_;
    int nb = nbeta_ - nrstc_ - nfrzc_;
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            int hi = symmetry[i];
            int ifull = i - pitzer_offset[hi] + pitzer_offset_full[hi] + rstcpi_[hi] + frzcpi_[hi];

            int hj = symmetry[j];
            int jfull = j - pitzer_offset[hj] + pitzer_offset_full[hj] + rstcpi_[hj] + frzcpi_[hj];

            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];

                int hk = symmetry[k];
                int kfull = k - pitzer_offset[hk] + pitzer_offset_full[hk] + rstcpi_[hk] + frzcpi_[hk];

                int hl = symmetry[l];
                int lfull = l - pitzer_offset[hl] + pitzer_offset_full[hl] + rstcpi_[hl] + frzcpi_[hl];

                int hik = SymmetryPair(symmetry_full[ifull],symmetry_full[kfull]);
                int ik  = ibas_full_sym[hik][ifull][kfull];
                int jl  = ibas_full_sym[hik][jfull][lfull];

                int myoff = 0;
                for (int myh = 0; myh < hik; myh++) {
                    myoff += (long int)gems_full[myh] * ( (long int)gems_full[myh] + 1L ) / 2L;
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
            int ifull = i - pitzer_offset[hi] + pitzer_offset_full[hi] + rstcpi_[hi] + frzcpi_[hi];

            int hj = symmetry[j];
            int jfull = j - pitzer_offset[hj] + pitzer_offset_full[hj] + rstcpi_[hj] + frzcpi_[hj];

            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];

                int hk = symmetry[k];
                int kfull = k - pitzer_offset[hk] + pitzer_offset_full[hk] + rstcpi_[hk] + frzcpi_[hk];

                int hl = symmetry[l];
                int lfull = l - pitzer_offset[hl] + pitzer_offset_full[hl] + rstcpi_[hl] + frzcpi_[hl];

                int hik = SymmetryPair(symmetry_full[ifull],symmetry_full[kfull]);
                int ik  = ibas_full_sym[hik][ifull][kfull];
                int jl  = ibas_full_sym[hik][jfull][lfull];

                int myoff = 0;
                for (int myh = 0; myh < hik; myh++) {
                    myoff += (long int)gems_full[myh] * ( (long int)gems_full[myh] + 1L ) / 2L;
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
                    myoff += (long int)gems_full[myh] * ( (long int)gems_full[myh] + 1L ) / 2L;
                }

                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym_[myoff + INDEX(il,jk)];
                c_p[d2aaoff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym_[myoff + INDEX(jk,il)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym_[myoff + INDEX(il,jk)];
                c_p[d2bboff[h] + ij*gems_aa[h]+kl]   -= 0.5 * tei_full_sym_[myoff + INDEX(jk,il)];

            }
        }
    }
}

void v2RDMSolver::FrozenCoreEnergy() {


    // if frozen core, adjust oei's and compute frozen core energy:
    efzc_ = 0.0;
    offset = 0;
    long int offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {

            int ifull = i + offset;
            int ii    = ibas_full_sym[0][ifull][ifull];

            efzc_ += 2.0 * oei_full_sym_[offset3 + INDEX(i,i)];

            long int offset2 = 0;
            for (int h2 = 0; h2 < nirrep_; h2++) {
                for (int j = 0; j < rstcpi_[h2] + frzcpi_[h2]; j++) {

                    int jfull = j + offset2;
                    int jj    = ibas_full_sym[0][jfull][jfull];
                    int hij = SymmetryPair(h,h2);

                    int ij    = ibas_full_sym[hij][ifull][jfull];

                    if ( ! is_df_ ) {
                        int myoff = 0;
                        for (int myh = 0; myh < hij; myh++) {
                            myoff += (long int)gems_full[myh] * ( (long int)gems_full[myh] + 1L ) / 2L;
                        }
                        efzc_ += (2.0 * tei_full_sym_[INDEX(ii,jj)] - tei_full_sym_[myoff + INDEX(ij,ij)]);
                    }else {
                        double dum1 = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(ifull,ifull),1,Qmo_+nQ_*INDEX(jfull,jfull),1);
                        double dum2 = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(ifull,jfull),1,Qmo_+nQ_*INDEX(ifull,jfull),1);
                        efzc_ += 2.0 * dum1 - dum2;
                    }

                }
                offset2 += nmopi_[h2] - frzvpi_[h2];
            }
        }
        offset += nmopi_[h] - frzvpi_[h];
        offset3 += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    double * c_p = c->pointer();

    // adjust one-electron integrals for core repulsion contribution
    offset = 0;
    offset3 = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {

            int ifull = i + offset;

            for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; j++) {

                int jfull = j + offset;

                int ij = ibas_full_sym[0][ifull][jfull];

                double dum = 0.0;

                long int offset2 = 0;
                for (int h2 = 0; h2 < nirrep_; h2++) {
                    for (int k = 0; k < rstcpi_[h2] + frzcpi_[h2]; k++) {

                        int kfull = k + offset2;

                        int kk = ibas_full_sym[0][kfull][kfull];

                        int hik = SymmetryPair(h,h2);

                        int ik = ibas_full_sym[hik][ifull][kfull];
                        int jk = ibas_full_sym[hik][jfull][kfull];



                        if ( ! is_df_ ) {
                            int myoff = 0;
                            for (int myh = 0; myh < hik; myh++) {
                                myoff += (long int)gems_full[myh] * ( (long int)gems_full[myh] + 1L ) / 2L;
                            }
                            dum += (2.0 * tei_full_sym_[INDEX(ij,kk)] - tei_full_sym_[myoff+INDEX(ik,jk)]);
                        }else {
                            double dum1 = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(ifull,jfull),1,Qmo_+nQ_*INDEX(kfull,kfull),1);
                            double dum2 = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(ifull,kfull),1,Qmo_+nQ_*INDEX(jfull,kfull),1);
                            dum += 2.0 * dum1 - dum2;
                        }


                    }
                    offset2 += nmopi_[h2] - frzvpi_[h2];
                }
                c_p[d1aoff[h] + (i-rstcpi_[h]-frzcpi_[h])*amopi_[h] + (j-rstcpi_[h]-frzcpi_[h])] = oei_full_sym_[offset3+INDEX(i,j)];
                c_p[d1boff[h] + (i-rstcpi_[h]-frzcpi_[h])*amopi_[h] + (j-rstcpi_[h]-frzcpi_[h])] = oei_full_sym_[offset3+INDEX(i,j)];
                c_p[d1aoff[h] + (i-rstcpi_[h]-frzcpi_[h])*amopi_[h] + (j-rstcpi_[h]-frzcpi_[h])] += dum;
                c_p[d1boff[h] + (i-rstcpi_[h]-frzcpi_[h])*amopi_[h] + (j-rstcpi_[h]-frzcpi_[h])] += dum;
            }
        }
        offset += nmopi_[h] - frzvpi_[h];
        offset3 += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

}



}}
