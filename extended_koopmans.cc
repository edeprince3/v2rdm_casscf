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
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/mintshelper.h>

using namespace psi;

namespace psi{ namespace v2rdm_casscf {

extern "C" {
void F77NAME(dgeev)(char &jobvl,char &jobvr,long int &n,double *a,long int &lda,
           double *wr,double *wi, double *vl,long int &ldvl,double *vr,
           long int &ldvr,double * work,long int &lwork, long int &info);
};
inline void DGEEV(char &jobvl,char &jobvr,long int &n,double*a,long int &lda,
           double *wr,double *wi, double *vl,long int &ldvl,double *vr,
           long int &ldvr,double * work,long int &lwork,long int &info)
{
   F77NAME(dgeev)(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info);
};

void v2RDMSolver::ExtendedKoopmans() {

    Dimension noccpi(nirrep_,"Number of occupied / partially occupied orbitals per irrep");
    for (int h = 0; h < nirrep_; h++) {
        noccpi[h] = frzcpi_[h] + rstcpi_[h] + amopi_[h];
    }

    // OPDM
    std::shared_ptr<Matrix> Da (new Matrix(nirrep_,noccpi,noccpi));
    std::shared_ptr<Matrix> Db (new Matrix(nirrep_,noccpi,noccpi));

    double * x_p = x->pointer();

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
            Da->pointer(h)[i][i] = 1.0;
        }
        int off = rstcpi_[h] + frzcpi_[h];
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                Da->pointer(h)[i + off][j + off] = x_p[d1aoff[h] + i * amopi_[h] + j];
                Db->pointer(h)[i + off][j + off] = x_p[d1boff[h] + i * amopi_[h] + j];
            }
        }
    }
    std::shared_ptr<Matrix> Da_save (new Matrix(Da));


    // one-electron integrals

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    std::shared_ptr<Matrix> hcore (new Matrix(mints->so_kinetic()));
    std::shared_ptr<Matrix> potential (new Matrix(mints->so_potential()));
    hcore->add(potential);
    hcore->transform(Ca_);

    // V = < p*[H,q] >
    std::shared_ptr<Matrix> Va (new Matrix(nirrep_,noccpi,noccpi));

    // core-core Vij = < i*[H,j] >
    int ioff = 0;
    int oei_off = 0;
    for (int h = 0; h < nirrep_; h++) {

        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {

            int ifull = i + ioff;

            for (int j = 0; j < rstcpi_[h] + frzcpi_[h]; j++) {

                int jfull = j + ioff;

                double dum = 0.0;

                long int koff = 0;
                for (int hk = 0; hk < nirrep_; hk++) {

                    // sum over core
                    for (int k = 0; k < rstcpi_[hk] + frzcpi_[hk]; k++) {

                        int kfull = k + koff;

                        int hik = SymmetryPair(h,hk);

                        double dum1 = TEI(jfull,ifull,kfull,kfull,0);
                        double dum2 = TEI(jfull,kfull,kfull,ifull,hik);
                        dum += 2.0 * dum1 - dum2;

                    }
                    koff += nmopi_[hk] - frzvpi_[hk];
                }

                long int toff = 0;
                for (int ht = 0; ht < nirrep_; ht++) {
                    // sum over active
                    for (int t = 0; t < amopi_[ht]; t++) {

                        int tfull = t + rstcpi_[ht] + frzcpi_[ht] + toff;

                        int hit = h ^ ht;

                        double dum1 = TEI(jfull,ifull,tfull,tfull,0);
                        double dum2 = TEI(jfull,tfull,tfull,ifull,hit);
                        dum += x_p[d1aoff[ht] + t * amopi_[ht] + t] * (dum1 - dum2);
                        dum += x_p[d1boff[ht] + t * amopi_[ht] + t] * dum1;

                    }
                    toff += nmopi_[ht] - frzvpi_[ht];
                }

                //Va->pointer(h)[i][j] = -oei_full_sym_[oei_off+INDEX(j,i)] - dum;
                Va->pointer(h)[i][j] = -hcore->pointer(h)[j][i] - dum;
            }
        }
        ioff += nmopi_[h] - frzvpi_[h];
        oei_off += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // core-active Vit = < i*[H,t] >
    ioff = 0;
    oei_off = 0;
    for (int h = 0; h < nirrep_; h++) {

        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {

            int ifull = i + ioff;

            for (int t = 0; t < amopi_[h]; t++) {

                int tfull = t + rstcpi_[h] + frzcpi_[h] + ioff;

                double dum = 0.0;

                long int koff = 0;
                for (int hk = 0; hk < nirrep_; hk++) {

                    // sum over core
                    for (int k = 0; k < rstcpi_[hk] + frzcpi_[hk]; k++) {

                        int kfull = k + koff;

                        int hki = hk ^ h;

                        double dum1 = TEI(tfull,ifull,kfull,kfull,0);
                        double dum2 = TEI(tfull,kfull,kfull,ifull,hki);
                        dum += 2.0 * dum1 - dum2;

                    }
                    koff += nmopi_[hk] - frzvpi_[hk];
                }

                long int uoff = 0;
                for (int hu = 0; hu < nirrep_; hu++) {

                    // sum over active
                    for (int u = 0; u < amopi_[hu]; u++) {

                        int ufull = u + rstcpi_[hu] + frzcpi_[hu] + uoff;

                        int hui = hu ^ h;

                        double dum1 = TEI(tfull,ifull,ufull,ufull,0);
                        double dum2 = TEI(tfull,ufull,ufull,ifull,hui);
                        dum += x_p[d1aoff[hu] + u * amopi_[hu] + u] * (dum1 - dum2);
                        dum += x_p[d1boff[hu] + u * amopi_[hu] + u] * dum1;

                    }
                    uoff += nmopi_[hu] - frzvpi_[hu];
                }

                //Va->pointer(h)[i][t + rstcpi_[h] + frzcpi_[h]] = -oei_full_sym_[oei_off+INDEX(t + rstcpi_[h] + frzcpi_[h],i)] - dum;
                Va->pointer(h)[i][t + rstcpi_[h] + frzcpi_[h]] = -hcore->pointer(h)[t + rstcpi_[h] + frzcpi_[h]][i] - dum;
            }
        }
        ioff += nmopi_[h] - frzvpi_[h];
        oei_off += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // active-core Vti = < t*[H,i] >
    int toff = 0;
    oei_off = 0;
    for (int h = 0; h < nirrep_; h++) {

        for (int t = 0; t < amopi_[h]; t++) {

            int tfull = t + rstcpi_[h] + frzcpi_[h] + toff;

            for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {
            
                int ifull = i + toff;

                double dum = 0.0;

                // one-electron part (sum over active) 
                // - 1D(tu) h(iu)
                for (int u = 0; u < amopi_[h]; u++) {
                    //dum += x_p[d1aoff[h] + t * amopi_[h] + u]  * oei_full_sym_[oei_off+INDEX(i,u + rstcpi_[h] + frzcpi_[h])];
                    dum += x_p[d1aoff[h] + t * amopi_[h] + u]  * hcore->pointer(h)[i][u + rstcpi_[h] + frzcpi_[h]];
                }

                // sum over core + active
                // - 1D(tu) [ 2(iu|jj) - (ij|ju) ]
                for (int u = 0; u < amopi_[h]; u++) {

                    int ufull = u + rstcpi_[h] + frzcpi_[h] + toff;

                    int joff = 0;

                    for (int hj = 0; hj < nirrep_; hj++) {

                        for (int j = 0; j < rstcpi_[hj] + frzcpi_[hj]; j++) {

                            int jfull = j + joff;

                            int hju = hj ^ h;

                            double dum1 = TEI(ifull,ufull,jfull,jfull,0);
                            double dum2 = TEI(ifull,jfull,jfull,ufull,hju);
                            dum += x_p[d1aoff[h] + t * amopi_[h] + u] * ( 2.0 * dum1 - dum2 );

                        }

                        joff += nmopi_[hj] - frzvpi_[hj];
                    }
                }

                int uoff = 0;
                for (int hu = 0; hu < nirrep_; hu++) {

                    for (int u = 0; u < amopi_[hu]; u++) {
                        // sum over active
                        // D2(tu;vw) (iv|uw)

                        int ufull = u + rstcpi_[hu] + frzcpi_[hu] + uoff;

                        int htu = h ^ hu;
                        int hvw = htu;

                        int tt = t + pitzer_offset[h];
                        int uu = u + pitzer_offset[hu];

                        int tu_ab = ibas_ab_sym[htu][tt][uu];

                        for (int vw_ab = 0; vw_ab < gems_ab[hvw]; vw_ab++) {
                            int vv = bas_ab_sym[hvw][vw_ab][0];
                            int ww = bas_ab_sym[hvw][vw_ab][1];

                            int hv = symmetry[vv];
                            int hw = symmetry[ww];

                            int w = ww - pitzer_offset[hw];
                            int v = vv - pitzer_offset[hv];

                            int woff = 0;
                            for (int myh = 0; myh < hw; myh++) {
                                woff += nmopi_[myh] - frzvpi_[myh];
                            }
                            int wfull = w + rstcpi_[hw] + frzcpi_[hw] + woff;

                            int voff = 0;
                            for (int myh = 0; myh < hv; myh++) {
                                voff += nmopi_[myh] - frzvpi_[myh];
                            }
                            int vfull = v + rstcpi_[hv] + frzcpi_[hv] + voff;

                            int hiv = h ^ hv;
                            double eri = TEI(ifull,vfull,ufull,wfull,hiv);

                            dum += eri * x_p[d2aboff[htu] + tu_ab * gems_ab[htu] + vw_ab];
                            
                            if ( tt != uu && vv != ww ) {
                                int tu_aa = ibas_aa_sym[htu][tt][uu];
                                int vw_aa = ibas_aa_sym[hvw][vv][ww];
                                int sg = 1;
                                if ( tt > uu ) sg = -sg;
                                if ( ww > vv ) sg = -sg;
                                dum += eri * sg * x_p[d2aaoff[htu] + tu_aa * gems_aa[htu] + vw_aa];
                            }

                        }
                    }
                    uoff += nmopi_[hu] - frzvpi_[hu];
                }

                Va->pointer(h)[t + rstcpi_[h] + frzcpi_[h]][i] = -dum;
            }
        }
        toff += nmopi_[h] - frzvpi_[h];
        oei_off += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // active-active Vtu = < t*[H,u] >
    toff = 0;
    oei_off = 0;
    for (int h = 0; h < nirrep_; h++) {

        for (int t = 0; t < amopi_[h]; t++) {

            int tfull = t + rstcpi_[h] + frzcpi_[h] + toff;

            for (int u = 0; u < amopi_[h]; u++) {
            
                int ufull = u + rstcpi_[h] + frzcpi_[h] + toff;

                double dum = 0.0;

                // one-electron part (sum over active) 
                // - 1D(tv) h(uv)
                for (int v = 0; v < amopi_[h]; v++) {
                    //dum += x_p[d1aoff[h] + t * amopi_[h] + v]  * oei_full_sym_[oei_off+INDEX(u + rstcpi_[h] + frzcpi_[h],v + rstcpi_[h] + frzcpi_[h])];
                    dum += x_p[d1aoff[h] + t * amopi_[h] + v]  * hcore->pointer(h)[u + rstcpi_[h] + frzcpi_[h]][v + rstcpi_[h] + frzcpi_[h]];
                }

                // sum over core + active
                // - 1D(tv) [ 2(uv|ii) - (ui|iv) ]
                for (int v = 0; v < amopi_[h]; v++) {
                    int vfull = v + rstcpi_[h] + frzcpi_[h] + toff;

                    int ioff = 0;
                    for (int hi = 0; hi < nirrep_; hi++) {

                        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {

                            int ifull = i + ioff;

                            int hui = h ^ hi;

                            double dum1 = TEI(ufull,vfull,ifull,ifull,0);
                            double dum2 = TEI(ufull,ifull,ifull,vfull,hui);
                            dum += x_p[d1aoff[h] + t * amopi_[h] + v] * ( 2.0 * dum1 - dum2 );

                        }

                        ioff += nmopi_[hi] - frzvpi_[hi];
                    }
                }

                int voff = 0;
                for (int hv = 0; hv < nirrep_; hv++) {

                    for (int v = 0; v < amopi_[hv]; v++) {
                        // sum over active
                        // D2(tv;wx) (uw|vx)

                        int htv = h ^ hv;
                        int hwx = htv;

                        int vfull = v + rstcpi_[hv] + frzcpi_[hv] + voff;

                        int tt = t + pitzer_offset[h];
                        int vv = v + pitzer_offset[hv];

                        int tv_ab = ibas_ab_sym[htv][tt][vv];

                        for (int wx_ab = 0; wx_ab < gems_ab[hwx]; wx_ab++) {
                            int ww = bas_ab_sym[hwx][wx_ab][0];
                            int xx = bas_ab_sym[hwx][wx_ab][1];

                            int hw = symmetry[ww];
                            int hx = symmetry[xx];

                            int w = ww - pitzer_offset[hw];
                            int x = xx - pitzer_offset[hx];

                            int woff = 0;
                            for (int myh = 0; myh < hw; myh++) {
                                woff += nmopi_[myh] - frzvpi_[myh];
                            }
                            int wfull = w + rstcpi_[hw] + frzcpi_[hw] + woff;

                            int xoff = 0;
                            for (int myh = 0; myh < hx; myh++) {
                                xoff += nmopi_[myh] - frzvpi_[myh];
                            }
                            int xfull = x + rstcpi_[hx] + frzcpi_[hx] + xoff;

                            int huw = h ^ hw;
                            double eri = TEI(ufull,wfull,vfull,xfull,huw);

                            dum += eri * x_p[d2aboff[htv] + tv_ab * gems_ab[htv] + wx_ab];
                            
                            if ( tt != vv && ww != xx ) {
                                int wx_aa = ibas_aa_sym[hwx][ww][xx];
                                int tv_aa = ibas_aa_sym[htv][tt][vv];
                                int sg = 1;
                                if ( tt > vv ) sg = -sg;
                                if ( ww > xx ) sg = -sg;
                                dum += eri * sg * x_p[d2aaoff[htv] + tv_aa * gems_aa[htv] + wx_aa];
                            }

                        }
                    }
                    voff += nmopi_[hv] - frzvpi_[hv];
                }

                Va->pointer(h)[t + rstcpi_[h] + frzcpi_[h]][u + rstcpi_[h] + frzcpi_[h]] = -dum;
            }
        }
        toff += nmopi_[h] - frzvpi_[h];
        oei_off += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // print
    //Va->print();

    outfile->Printf("\n");
    outfile->Printf("    ==> Extended Koopman's Theorem <==\n");
    outfile->Printf("\n");

    // now ... diagonalize each block using nonsymmetric eigensolver 
    for (int h = 0; h < nirrep_; h++) {

        long int N = noccpi[h];
        if ( N == 0 ) continue;

        // symmetrize
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                double dum = 0.5 * ( Va->pointer(h)[i][j] + Va->pointer(h)[j][i] );
                Va->pointer(h)[i][j] = dum;
                Va->pointer(h)[j][i] = dum;
            }
        }
        // 1 / D * V
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Va->pointer(h)[i][j] /= Da->pointer(h)[i][i];
            }
        }

        outfile->Printf("    ");
        outfile->Printf("    Irrep");
        outfile->Printf("    State");
        outfile->Printf("            IP");
        outfile->Printf("     Orb. Occ.");
        outfile->Printf("\n");
        //printf("    Irrep: %2i\n",h);


/*
        char JOBVR = 'V';
        char JOBVL = 'V';

        long int LDA  = N;
        long int LDB  = N;
        long int LDVR = N;
        long int LDVL = N;
        long int LWORK = 8*N;

        double*WORK=(double*)malloc(LWORK*sizeof(double));
        double*ALPHAR=(double*)malloc(N*sizeof(double));
        double*ALPHAI=(double*)malloc(N*sizeof(double));
        double*BETA=(double*)malloc(N*sizeof(double));
        double*VR=(double*)malloc(N*N*sizeof(double));
        double*VL=(double*)malloc(N*N*sizeof(double));

        long int INFO = 0;
        INFO = C_DGGEV(JOBVL,JOBVR,N,Va->pointer(h)[0],LDA,Da->pointer(h)[0],LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK);

        // check orthonormality
        //for (int i = 0; i < N; i++) {
        //    for (int j = 0; j < N; j++) {
        //        double dum1 = 0.0;
        //        double dum2 = 0.0;
        //        for (int k = 0; k < N; k++) {
        //            dum1 += VL[i*N+k] * VR[k*N+j];
        //            dum2 += VL[k*N+i] * VR[j*N+k];
        //        }
        //        printf("ortho %5i %5i %20.12lf %20.12lf\n",i,j,dum1,dum2);
        //    }
        //}

        for (int i = 0; i < N; i++) {
            double max = 0.0;
            int jmax = -999;
            for (int j = 0; j < N; j++) {
printf("%5i %5i %20.12lf %20.12lf\n",i,j,VL[i*N+j],VR[i*N+j]);
                if ( fabs(VL[i*N+j]) > max ) {
                    max = fabs(VL[i*N+j]);
                    jmax = j;
                }
            }
            double jocc = Da_save->pointer(h)[jmax][jmax];
            outfile->Printf("    ");
            outfile->Printf("    %5i",h);
            outfile->Printf("    %5i",i);
            outfile->Printf("%14.6lf",ALPHAR[i]/BETA[i]);
            outfile->Printf("%14.6lf",jocc);
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
  
        free(VR);
        free(VL);
        free(WORK);
        free(ALPHAR);
        free(ALPHAI);
        free(BETA); 
*/
        long int info;
        char jobl = 'V';
        char jobr = 'V';

        long int lwork = 4*N;
        double * work  = (double*)malloc(lwork*sizeof(double));
        double * wi = (double*)malloc(N*sizeof(double));
        double * eigval = (double*)malloc(N*sizeof(double));
        double * vl = (double*)malloc(N*N*sizeof(double));
        double * vr = (double*)malloc(N*N*sizeof(double));

        DGEEV(jobl,jobr,N,Va->pointer(h)[0],N,eigval,wi,vl,N,vr,N,work,lwork,info);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("coeffs %5i %5i %20.12lf %20.12lf %20.12lf\n",i,j,vl[i*N+j],vr[i*N+j],Da_save->pointer(h)[j][j]);
            }
        }
        // check orthonormality
        //for (int i = 0; i < N; i++) {
        //    for (int j = 0; j < N; j++) {
        //        double dum1 = 0.0;
        //        double dum2 = 0.0;
        //        for (int k = 0; k < N; k++) {
        //            //dum1 += vl[i*N+k] * vr[j*N+k];
        //            //dum2 += vl[k*N+i] * vr[k*N+j];
        //            dum1 += vl[i*N+k] * vl[j*N+k];
        //            dum2 += vr[i*N+k] * vr[j*N+k];
        //        }
        //        printf("ortho %5i %5i %20.12lf %20.12lf\n",i,j,dum1,dum2);
        //    }
        //}

        for (int i = 0; i < N; i++) {
            double max = 0.0;
            int jmax = -999;
            for (int j = 0; j < N; j++) {
                if ( fabs(vr[j*N+i]) > max ) {
                    max = fabs(vr[j*N+i]);
                    jmax = j;
                }
            }
            double jocc = Da_save->pointer(h)[jmax][jmax];
            outfile->Printf("    ");
            outfile->Printf("    %5i",h);
            outfile->Printf("    %5i",i);
            outfile->Printf("%14.6lf",eigval[i]);
            outfile->Printf("%14.6lf",jocc);
            outfile->Printf("\n");
        }
        outfile->Printf("\n");

        free(wi);
        free(work);
        free(vr);
        free(vl);
        free(eigval);
    }

}

}} // End namespaces

