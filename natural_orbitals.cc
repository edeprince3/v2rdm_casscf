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

#include "v2rdm_solver.h"


using namespace psi;
//using namespace fnocc;

namespace psi{ namespace v2rdm_casscf {

void v2RDMSolver::ComputeNaturalOrbitals() {

    if ( nalpha_ != nbeta_ || !constrain_spin_ ) {
        throw PsiException("natural orbital transformation only implemented for closed-shell systems with spin symmetry enforced",__FILE__,__LINE__);
    }

    SharedMatrix D (new Matrix(nirrep_,amopi_,amopi_));
    SharedMatrix eigvec (new Matrix(nirrep_,amopi_,amopi_));
    SharedVector eigval (new Vector("Natural Orbital Occupation Numbers",nirrep_,amopi_));

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                D->pointer(h)[i][j]  = x->pointer()[d1aoff[h]+i*amopi_[h]+j];
                D->pointer(h)[i][j] += x->pointer()[d1boff[h]+i*amopi_[h]+j];
            }
        }
    }
    D->diagonalize(eigvec,eigval,descending);

    // build AO/NO transformation matrix (both Ca_ and Cb_)
    for (int h = 0; h < nirrep_; h++) {

        for (int mu = 0; mu < nsopi_[h]; mu++) {

            double *  temp = (double*)malloc(nmopi_[h]*sizeof(double));
            double ** ep   = eigvec->pointer(h);

            // Ca_

            double ** cp = Ca_->pointer(h);
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                double dum = 0.0;
                for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                    dum += cp[mu][j] * ep[j-rstcpi_[h]-frzcpi_[h]][i-rstcpi_[h]-frzcpi_[h]];
                }
                temp[i] = dum;
            }
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                cp[mu][i] = temp[i];
            }


            // Cb_
            cp = Cb_->pointer(h);
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                double dum = 0.0;
                for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                    dum += cp[mu][j] * ep[j-rstcpi_[h]-frzcpi_[h]][i-rstcpi_[h]-frzcpi_[h]];
                }
                temp[i] = dum;
            }
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                cp[mu][i] = temp[i];
            }

            free(temp);
        }
    }

    // transform the alpha 1-RDM to natural orbital basis
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                D->pointer(h)[i][j]  = x->pointer()[d1aoff[h]+i*amopi_[h]+j];
            }
        }
    }
    D->transform(eigvec);
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                x->pointer()[d1aoff[h]+i*amopi_[h]+j] = D->pointer(h)[i][j];
            }
        }
    }

    // transform the beta 1-RDM to natural orbital basis
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                D->pointer(h)[i][j]  = x->pointer()[d1boff[h]+i*amopi_[h]+j];
            }
        }
    }
    D->transform(eigvec);
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                x->pointer()[d1boff[h]+i*amopi_[h]+j] = D->pointer(h)[i][j];
            }
        }
    }


    // transform the ab block of the 2-RDM to natural orbital basis
    std::shared_ptr<Vector> tempx (new Vector(dimx_));

    // D2ab(ij,kl): transform index 4
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double dum = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int kp = ibas_ab_sym[h][k][p];

                                dum += x->pointer()[d2aboff[h]+ij*gems_ab[h]+kp] * ep[pp][ll];

                            }

                            int kl = ibas_ab_sym[h][k][l];

                            tempx->pointer()[d2aboff[h]+ij*gems_ab[h]+kl] = dum;

                        }
                    }
                }
            }
        }
    }

    // D2ab(ij,lk): transform index 3
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double dum = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int pk = ibas_ab_sym[h][p][k];

                                dum += tempx->pointer()[d2aboff[h]+ij*gems_ab[h]+pk] * ep[pp][ll];

                            }

                            int lk = ibas_ab_sym[h][l][k];

                            x->pointer()[d2aboff[h]+ij*gems_ab[h]+lk] = dum;

                        }
                    }
                }
            }
        }
    }

    // D2ab(kl,ij): transform index 2
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double dum = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int kp = ibas_ab_sym[h][k][p];

                                dum += x->pointer()[d2aboff[h]+kp*gems_ab[h]+ij] * ep[pp][ll];

                            }

                            int kl = ibas_ab_sym[h][k][l];

                            tempx->pointer()[d2aboff[h]+kl*gems_ab[h]+ij] = dum;

                        }
                    }
                }
            }
        }
    }

    // D2ab(lk,ij): transform index 1
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double dum = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int pk = ibas_ab_sym[h][p][k];

                                dum += tempx->pointer()[d2aboff[h]+pk*gems_ab[h]+ij] * ep[pp][ll];

                            }

                            int lk = ibas_ab_sym[h][l][k];

                            x->pointer()[d2aboff[h]+lk*gems_ab[h]+ij] = dum;

                        }
                    }
                }
            }
        }
    }

    // now transform the aa and bb blocks of the 2-RDM to natural orbital basis

    // D2aa(ij,kl): transform index 4
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];

                            if ( k <= l ) continue; // TODO: check basis.cc to be sure we shouldn't be skipping k >= l
    
                            double duma = 0.0;
                            double dumb = 0.0;

                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];
                                if ( p == k ) continue;

                                int kp = ibas_aa_sym[h][k][p];

                                int sg = 1;
                                if ( p > k ) sg = -1;

                                duma += sg * x->pointer()[d2aaoff[h]+ij*gems_aa[h]+kp] * ep[pp][ll];
                                dumb += sg * x->pointer()[d2bboff[h]+ij*gems_aa[h]+kp] * ep[pp][ll];

                            }

                            int kl = ibas_aa_sym[h][k][l];

                            tempx->pointer()[d2aaoff[h]+ij*gems_aa[h]+kl] = duma;
                            tempx->pointer()[d2bboff[h]+ij*gems_aa[h]+kl] = dumb;

                        }
                    }
                }
            }
        }
    }

    // D2aa(ij,lk): transform index 3
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];

                            if ( l <= k ) continue; 
    
                            double duma = 0.0;
                            double dumb = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];
                                if ( p == k ) continue;

                                int pk = ibas_aa_sym[h][p][k];

                                int sg = 1;
                                if ( k > p ) sg = -1;

                                duma += sg * tempx->pointer()[d2aaoff[h]+ij*gems_aa[h]+pk] * ep[pp][ll];
                                dumb += sg * tempx->pointer()[d2bboff[h]+ij*gems_aa[h]+pk] * ep[pp][ll];

                            }

                            int lk = ibas_aa_sym[h][l][k];

                            x->pointer()[d2aaoff[h]+ij*gems_aa[h]+lk] = duma;
                            x->pointer()[d2bboff[h]+ij*gems_aa[h]+lk] = dumb;

                        }
                    }
                }
            }
        }
    }

    // D2aa(kl,ij): transform index 2
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];

                            if ( k <= l ) continue; // TODO: check basis.cc to be sure we shouldn't be skipping k >= l
    
                            double duma = 0.0;
                            double dumb = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                if ( p == k ) continue;

                                int kp = ibas_aa_sym[h][k][p];

                                int sg = 1;
                                if ( p > k ) sg = -1;

                                duma += sg * x->pointer()[d2aaoff[h]+kp*gems_aa[h]+ij] * ep[pp][ll];
                                dumb += sg * x->pointer()[d2bboff[h]+kp*gems_aa[h]+ij] * ep[pp][ll];

                            }

                            int kl = ibas_aa_sym[h][k][l];

                            tempx->pointer()[d2aaoff[h]+kl*gems_aa[h]+ij] = duma;
                            tempx->pointer()[d2bboff[h]+kl*gems_aa[h]+ij] = dumb;

                        }
                    }
                }
            }
        }
    }

    // D2aa(lk,ij): transform index 1
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];

                            if ( l <= k ) continue; 
    
                            double duma = 0.0;
                            double dumb = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                if ( p == k ) continue;

                                int pk = ibas_aa_sym[h][p][k];

                                int sg = 1;
                                if ( k > p ) sg = -1;

                                duma += sg * tempx->pointer()[d2aaoff[h]+pk*gems_aa[h]+ij] * ep[pp][ll];
                                dumb += sg * tempx->pointer()[d2bboff[h]+pk*gems_aa[h]+ij] * ep[pp][ll];

                            }

                            int lk = ibas_aa_sym[h][l][k];

                            x->pointer()[d2aaoff[h]+lk*gems_aa[h]+ij] = duma;
                            x->pointer()[d2bboff[h]+lk*gems_aa[h]+ij] = dumb;

                        }
                    }
                }
            }
        }
    }

}

void v2RDMSolver::PrintNaturalOrbitalOccupations() {

    SharedMatrix Da (new Matrix(nirrep_,nmopi_,nmopi_));
    SharedMatrix eigveca (new Matrix(nirrep_,nmopi_,nmopi_));
    SharedVector eigvala (new Vector("Natural Orbital Occupation Numbers (alpha)",nirrep_,nmopi_));

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
            Da->pointer(h)[i][i] = 1.0;
        }
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {
            for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                Da->pointer(h)[i][j] = x->pointer()[d1aoff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
            }
        }
    }
    SharedMatrix saveda ( new Matrix(Da) );
    Da->diagonalize(eigveca,eigvala,descending);
    eigvala->print();

    SharedMatrix Db (new Matrix(nirrep_,nmopi_,nmopi_));
    SharedMatrix eigvecb (new Matrix(nirrep_,nmopi_,nmopi_));
    SharedVector eigvalb (new Vector("Natural Orbital Occupation Numbers (beta)",nirrep_,nmopi_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {
            Db->pointer(h)[i][i] = 1.0;
        }
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
            for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                Db->pointer(h)[i][j] = x->pointer()[d1boff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
            }
        }
    }
    Db->diagonalize(eigvecb,eigvalb,descending);
    eigvalb->print();

}


}} // End namespaces

