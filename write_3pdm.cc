/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
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
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 * 
 *@END LICENSE
 *
 */

#include "psi4-dec.h"
#include <psifiles.h>
#include <libiwl/iwl.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>

#include "v2rdm_solver.h"

using namespace psi;

namespace psi{namespace v2rdm_casscf{

struct dm3 {
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double val;
};


void v2RDMSolver::WriteActive3PDM(){

    double * x_p = x->pointer();

    boost::shared_ptr<PSIO> psio (new PSIO());

    psio->open(PSIF_V2RDM_D3AAA,PSIO_OPEN_NEW);
    psio->open(PSIF_V2RDM_D3AAB,PSIO_OPEN_NEW);
    psio->open(PSIF_V2RDM_D3BBA,PSIO_OPEN_NEW);
    psio->open(PSIF_V2RDM_D3BBB,PSIO_OPEN_NEW);

    psio_address addr_aaa = PSIO_ZERO;
    psio_address addr_aab = PSIO_ZERO;
    psio_address addr_bba = PSIO_ZERO;
    psio_address addr_bbb = PSIO_ZERO;

    long int countaaa = 0;
    long int countaab = 0;
    long int countbba = 0;
    long int countbbb = 0;

    // active-active part

    // aab, bba
    for (int h = 0; h < nirrep_; h++) {
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i     = bas_aba_sym[h][ijk][0];
            int j     = bas_aba_sym[h][ijk][1];
            int k     = bas_aba_sym[h][ijk][2];

            if ( i == j ) continue;

            int ifull = full_basis[i];
            int jfull = full_basis[j];
            int kfull = full_basis[k];

            int ijk_aab = ibas_aab_sym[h][i][j][k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l     = bas_aba_sym[h][lmn][0];
                int m     = bas_aba_sym[h][lmn][1];
                int n     = bas_aba_sym[h][lmn][2];

                if ( l == m ) continue;

                int lfull = full_basis[l];
                int mfull = full_basis[m];
                int nfull = full_basis[n];

                int lmn_aab = ibas_aab_sym[h][l][m][n];

                double valaab = x_p[d3aaboff[h] + ijk_aab*trip_aab[h] + lmn_aab];
                double valbba = x_p[d3bbaoff[h] + ijk_aab*trip_aab[h] + lmn_aab];

                dm3 d3;
                
                int sijk = 1;
                if ( i > j ) sijk = -sijk;
                int slmn = 1;
                if ( l > m ) slmn = -slmn;

                // ijk; lmn
                d3.i   = ifull;
                d3.j   = jfull;
                d3.k   = kfull;
                d3.l   = lfull;
                d3.m   = mfull;
                d3.n   = nfull;
                d3.val = sijk*slmn*valaab;
                psio->write(PSIF_V2RDM_D3AAB,"D3aab",(char*)&d3,sizeof(dm3),addr_aab,&addr_aab);
                countaab++;

                d3.val = sijk*slmn*valbba;
                psio->write(PSIF_V2RDM_D3BBA,"D3bba",(char*)&d3,sizeof(dm3),addr_bba,&addr_bba);
                countbba++;

            }
        }
    }
    // aaa, bbb
    for (int h = 0; h < nirrep_; h++) {
        for (int ijk = 0; ijk < trip_aba[h]; ijk++) {

            int i     = bas_aba_sym[h][ijk][0];
            int j     = bas_aba_sym[h][ijk][1];
            int k     = bas_aba_sym[h][ijk][2];

            if ( i == j || i == k || j == k ) continue;

            int ifull = full_basis[i];
            int jfull = full_basis[j];
            int kfull = full_basis[k];

            int ijk_aaa = ibas_aaa_sym[h][i][j][k];

            for (int lmn = 0; lmn < trip_aba[h]; lmn++) {

                int l     = bas_aba_sym[h][lmn][0];
                int m     = bas_aba_sym[h][lmn][1];
                int n     = bas_aba_sym[h][lmn][2];

                if ( l == m || l == n || m == n ) continue;

                int lfull = full_basis[l];
                int mfull = full_basis[m];
                int nfull = full_basis[n];

                int lmn_aaa = ibas_aaa_sym[h][l][m][n];

                double valaaa = x_p[d3aaaoff[h] + ijk_aaa*trip_aaa[h] + lmn_aaa];
                double valbbb = x_p[d3bbboff[h] + ijk_aaa*trip_aaa[h] + lmn_aaa];

                dm3 d3;
                
                int sijk = 1;
                if ( i > j ) sijk = -sijk;
                if ( i > k ) sijk = -sijk;
                if ( j > k ) sijk = -sijk;
                int slmn = 1;
                if ( l > m ) slmn = -slmn;
                if ( l > n ) slmn = -slmn;
                if ( m > n ) slmn = -slmn;

                // ijk; lmn
                d3.i   = ifull;
                d3.j   = jfull;
                d3.k   = kfull;
                d3.l   = lfull;
                d3.m   = mfull;
                d3.n   = nfull;
                d3.val = sijk*slmn*valaaa;
                psio->write(PSIF_V2RDM_D3AAA,"D3aaa",(char*)&d3,sizeof(dm3),addr_aaa,&addr_aaa);
                countaaa++;

                d3.val = sijk*slmn*valbbb;
                psio->write(PSIF_V2RDM_D3BBB,"D3bbb",(char*)&d3,sizeof(dm3),addr_bbb,&addr_bbb);
                countbbb++;

            }
        }
    }

    // write the number of entries in each file
    psio->write_entry(PSIF_V2RDM_D3AAA,"length",(char*)&countaaa,sizeof(long int));
    psio->write_entry(PSIF_V2RDM_D3AAB,"length",(char*)&countaab,sizeof(long int));
    psio->write_entry(PSIF_V2RDM_D3BBA,"length",(char*)&countbba,sizeof(long int));
    psio->write_entry(PSIF_V2RDM_D3BBB,"length",(char*)&countbbb,sizeof(long int));

    // close files
    psio->close(PSIF_V2RDM_D3AAA,1);
    psio->close(PSIF_V2RDM_D3AAB,1);
    psio->close(PSIF_V2RDM_D3BBA,1);
    psio->close(PSIF_V2RDM_D3BBB,1);

}

void v2RDMSolver::Read3PDM(){

    boost::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D3AAA) ) return;
    if ( !psio->exists(PSIF_V2RDM_D3AAB) ) return;
    if ( !psio->exists(PSIF_V2RDM_D3AAB) ) return;
    if ( !psio->exists(PSIF_V2RDM_D3BBB) ) return;

    //Ca_->print();

    double * D3aaa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D3aab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D3bba = (double*)malloc(nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D3bbb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)D3aaa,'\0',nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D3aab,'\0',nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D3bba,'\0',nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D3bbb,'\0',nmo_*nmo_*nmo_*nmo_*nmo_*nmo_*sizeof(double));

    psio_address addr_aaa = PSIO_ZERO;
    psio_address addr_bbb = PSIO_ZERO;
    psio_address addr_aab = PSIO_ZERO;
    psio_address addr_bba = PSIO_ZERO;

    // aab
    psio->open(PSIF_V2RDM_D3AAB,PSIO_OPEN_OLD);

    long int naab;
    psio->read_entry(PSIF_V2RDM_D3AAB,"length",(char*)&naab,sizeof(long int));

    for (int count = 0; count < naab; count++) {
        dm3 d3;
        psio->read(PSIF_V2RDM_D3AAB,"D3aab",(char*)&d3,sizeof(dm3),addr_aab,&addr_aab);
        int i = d3.i;
        int j = d3.j;
        int k = d3.k;
        int l = d3.l;
        int m = d3.m;
        int n = d3.n;
        long int id = i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n;
        D3aab[id] = d3.val;
    }
    psio->close(PSIF_V2RDM_D3AAB,1);

    // bba
    psio->open(PSIF_V2RDM_D3BBA,PSIO_OPEN_OLD);

    long int nbba;
    psio->read_entry(PSIF_V2RDM_D3BBA,"length",(char*)&nbba,sizeof(long int));

    for (int count = 0; count < nbba; count++) {
        dm3 d3;
        psio->read(PSIF_V2RDM_D3BBA,"D3bba",(char*)&d3,sizeof(dm3),addr_bba,&addr_bba);
        int i = d3.i;
        int j = d3.j;
        int k = d3.k;
        int l = d3.l;
        int m = d3.m;
        int n = d3.n;
        long int id = i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n;
        D3bba[id] = d3.val;
    }
    psio->close(PSIF_V2RDM_D3BBA,1);

    // aaa
    psio->open(PSIF_V2RDM_D3AAA,PSIO_OPEN_OLD);

    long int naaa;
    psio->read_entry(PSIF_V2RDM_D3AAA,"length",(char*)&naaa,sizeof(long int));

    for (int count = 0; count < naaa; count++) {
        dm3 d3;
        psio->read(PSIF_V2RDM_D3AAA,"D3aaa",(char*)&d3,sizeof(dm3),addr_aaa,&addr_aaa);
        int i = d3.i;
        int j = d3.j;
        int k = d3.k;
        int l = d3.l;
        int m = d3.m;
        int n = d3.n;
        long int id = i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n;
        D3aaa[id] = d3.val;
    }
    psio->close(PSIF_V2RDM_D3AAA,1);

    // bbb
    psio->open(PSIF_V2RDM_D3BBB,PSIO_OPEN_OLD);

    long int nbbb;
    psio->read_entry(PSIF_V2RDM_D3BBB,"length",(char*)&nbbb,sizeof(long int));

    for (int count = 0; count < nbbb; count++) {
        dm3 d3;
        psio->read(PSIF_V2RDM_D3BBB,"D3bbb",(char*)&d3,sizeof(dm3),addr_bbb,&addr_bbb);
        int i = d3.i;
        int j = d3.j;
        int k = d3.k;
        int l = d3.l;
        int m = d3.m;
        int n = d3.n;
        long int id = i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n;
        D3bbb[id] = d3.val;
    }
    psio->close(PSIF_V2RDM_D3BBB,1);

    // check traces:
    /*double traaa = 0.0;
    double trbbb = 0.0;
    double traab = 0.0;
    double trbba = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                traaa += D3aaa[i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n];
                traab += D3aab[i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n];
                trbba += D3bba[i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n];
                trbbb += D3bbb[i*nmo_*nmo_*nmo_*nmo_*nmo_ + j*nmo_*nmo_*nmo_*nmo_ + k*nmo_*nmo_*nmo_ + l*nmo_*nmo_ + m*nmo_ + n];
            }
        }
    }
    printf("  tr(d3aaa) = %20.12lf\n",traaa);
    printf("  tr(d3aab) = %20.12lf\n",traab);
    printf("  tr(d3bba) = %20.12lf\n",trbba);
    printf("  tr(d3bbb) = %20.12lf\n",trbbb);*/

    // check energy from D2 = L(D3):
    /*double en2 = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {

                    double eri = C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,k),1,Qmo_+nQ_*INDEX(j,l),1);
                    
                    en2 +=       eri * D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                    en2 += 0.5 * eri * D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                    en2 += 0.5 * eri * D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];

                }
            }
        }
    }

    boost::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    boost::shared_ptr<Matrix> K1 (new Matrix(mints->so_potential()));
    K1->add(mints->so_kinetic());
    K1->transform(Ca_);

    double en1 = 0.0;

    long int offset = 0;
    long int offset2 = 0;
    for (int h = 0; h < nirrep_; h++) {

        for (int i = 0; i < nmopi_[h]; i++) {

            int ifull = i + offset;

            for (int j = 0; j < nmopi_[h]; j++) {

                int jfull = j + offset;


                en1 += oei_full_sym_[offset2 + INDEX(i,j)] * Da[ifull*nmo_+jfull];
                en1 += oei_full_sym_[offset2 + INDEX(i,j)] * Db[ifull*nmo_+jfull];
                //en1 += K1->pointer(h)[i][j] * Da[ifull*nmo_+jfull];
                //en1 += K1->pointer(h)[i][j] * Db[ifull*nmo_+jfull];
            }
        }

        offset  += nmopi_[h] - frzvpi_[h];
        offset2 += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;

    }

    printf("%20.12lf %20.12lf %20.12lf %20.12lf\n",en1,en2,enuc_,en1+en2+enuc_);

    free(D2aa);
    free(D2bb);
    free(D2ab);
    free(Da);
    free(Db);*/

    free(D3aaa);
    free(D3aab);
    free(D3bba);
    free(D3bbb);
}


}} //end namespaces


