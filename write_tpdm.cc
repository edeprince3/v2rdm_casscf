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

struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double val;
};


void v2RDMSolver::WriteTPDM(){

    double * x_p = x->pointer();

    boost::shared_ptr<PSIO> psio (new PSIO());

    psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_NEW);
    psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_NEW);
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_NEW);

    psio_address addr_aa = PSIO_ZERO;
    psio_address addr_bb = PSIO_ZERO;
    psio_address addr_ab = PSIO_ZERO;
    // active-active part
    for (int h = 0; h < nirrep_; h++) {

        for (int ij = 0; ij < gems_ab[h]; ij++) {

            int i     = bas_ab_sym[h][ij][0];
            int j     = bas_ab_sym[h][ij][1];
            int ifull = full_basis[i];
            int jfull = full_basis[j];

            for (int kl = 0; kl < gems_ab[h]; kl++) {

                int k     = bas_ab_sym[h][kl][0];
                int l     = bas_ab_sym[h][kl][1];
                int kfull = full_basis[k];
                int lfull = full_basis[l];

                double valab = x_p[d2aboff[h] + ij*gems_ab[h] + kl];

                tpdm d2;
                d2.i   = ifull;
                d2.j   = jfull;
                d2.k   = kfull;
                d2.l   = lfull;
                d2.val = valab;
                psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);

                if ( i != j && k != l ) {

                    int ija = ibas_aa_sym[h][i][j];
                    int kla = ibas_aa_sym[h][k][l];

                    int sij = i < j ? 1 : -1;
                    int skl = k < l ? 1 : -1;

                    double valaa = sij * skl * x_p[d2aaoff[h] + ija*gems_aa[h] + kla];
                    double valbb = sij * skl * x_p[d2bboff[h] + ija*gems_aa[h] + kla];

                    d2.val = valaa;
                    psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);

                    d2.val = valbb;
                    psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);

                }
            }
        }
    }

    // core-core
    for (int hi = 0; hi < nirrep_; hi++) {

        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];

            for (int hj = 0; hj < nirrep_; hj++) {

                for (int j = 0; j < rstcpi_[hj] + frzcpi_[hj]; j++) {

                    int jfull      = j + pitzer_offset_full[hj];

                    tpdm d2;
                    d2.i   = ifull;
                    d2.j   = jfull;
                    d2.k   = ifull;
                    d2.l   = jfull;
                    d2.val = 1.0;
                    psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);

                    if ( i != j ) {

                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);

                        // ij;ji

                        d2.val = -1.0;
                        d2.k   = jfull;
                        d2.l   = ifull;

                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);

                    }
                }
            }
        }
    }

    // core active; core active
    for (int hi = 0; hi < nirrep_; hi++) {

        for (int i = 0; i < rstcpi_[hi] + frzcpi_[hi]; i++) {

            int ifull      = i + pitzer_offset_full[hi];

            // D2(ij; il) ab, ba, aa, bb
            for (int hj = 0; hj < nirrep_; hj++) {

                for (int j = 0; j < amopi_[hj]; j++) {

                    int jfull      = full_basis[j+pitzer_offset[hj]];

                    for (int l = 0; l < amopi_[hj]; l++) {

                        int lfull      = full_basis[l+pitzer_offset[hj]];

                        // aa and bb pieces
                        double valaa = x_p[d1aoff[hj]+j*amopi_[hj]+l];
                        double valbb = x_p[d1boff[hj]+j*amopi_[hj]+l];

                        tpdm d2;

                        // ij;il
                        d2.i   = ifull;
                        d2.j   = jfull;
                        d2.k   = ifull;
                        d2.l   = lfull;

                        d2.val = valaa;
                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);

                        d2.val = valbb;
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);

                        // ij;li
                        d2.k   = lfull;
                        d2.l   = ifull;

                        d2.val = -valaa;
                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);

                        d2.val = -valbb;
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);

                        // ji;li
                        d2.i   = jfull;
                        d2.j   = ifull;

                        d2.val = valaa;
                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);

                        d2.val = valbb;
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);

                        // ji;il
                        d2.k   = ifull;
                        d2.l   = lfull;

                        d2.val = -valaa;
                        psio->write(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);

                        d2.val = -valbb;
                        psio->write(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);


                        // ab (ij;il) and ba (ji;li) pieces
                        double valab = x_p[d1boff[hj]+j*amopi_[hj]+l];
                        double valba = x_p[d1aoff[hj]+j*amopi_[hj]+l];

                        // ij;il
                        d2.i   = ifull;
                        d2.j   = jfull;

                        d2.val = valab;
                        psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);

                        // ji;li
                        d2.i   = jfull;
                        d2.j   = ifull;
                        d2.k   = lfull;
                        d2.l   = ifull;

                        d2.val = valba;
                        psio->write(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);

                    }
                }
            }
        }
    }

    // close files
    psio->close(PSIF_V2RDM_D2AA,1);
    psio->close(PSIF_V2RDM_D2BB,1);
    psio->close(PSIF_V2RDM_D2AB,1);
}

/*void v2RDMSolver::ReadTPDM(){

    int full = nmo + nfrzc;

    double * D2aa = (double*)malloc(full*full*full*full*sizeof(double));
    double * D2bb = (double*)malloc(full*full*full*full*sizeof(double));
    double * D2ab = (double*)malloc(full*full*full*full*sizeof(double));

    memset((void*)D2aa,'\0',full*full*full*full*sizeof(double));
    memset((void*)D2bb,'\0',full*full*full*full*sizeof(double));
    memset((void*)D2ab,'\0',full*full*full*full*sizeof(double));

    boost::shared_ptr<PSIO> psio (new PSIO());

    psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);
    psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    psio_address addr_aa = PSIO_ZERO;
    psio_address addr_bb = PSIO_ZERO;
    psio_address addr_ab = PSIO_ZERO;

    psio_address end_aa = psio_get_entry_end(PSIF_V2RDM_D2AA,"D2aa");
    psio_address end_bb = psio_get_entry_end(PSIF_V2RDM_D2BB,"D2bb");
    psio_address end_ab = psio_get_entry_end(PSIF_V2RDM_D2AB,"D2ab");

    do { 
        int i,j,k,l;
        double val;
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&i,sizeof(int),addr_ab,&addr_ab);
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&j,sizeof(int),addr_ab,&addr_ab);
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&k,sizeof(int),addr_ab,&addr_ab);
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&l,sizeof(int),addr_ab,&addr_ab);
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&val,sizeof(double),addr_ab,&addr_ab);
        long int id = i*full*full*full+j*full*full+k*full+l;
        D2ab[id] = val;
    }while(addr_ab != end_ab);
    do { 
        int i,j,k,l;
        double val;
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&i,sizeof(int),addr_aa,&addr_aa);
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&j,sizeof(int),addr_aa,&addr_aa);
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&k,sizeof(int),addr_aa,&addr_aa);
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&l,sizeof(int),addr_aa,&addr_aa);
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&val,sizeof(double),addr_aa,&addr_aa);
        long int id = i*full*full*full+j*full*full+k*full+l;
        D2aa[id] = val;
    }while(addr_aa != end_aa);
    do { 
        int i,j,k,l;
        double val;
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&i,sizeof(int),addr_bb,&addr_bb);
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&j,sizeof(int),addr_bb,&addr_bb);
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&k,sizeof(int),addr_bb,&addr_bb);
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&l,sizeof(int),addr_bb,&addr_bb);
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&val,sizeof(double),addr_bb,&addr_bb);
        long int id = i*full*full*full+j*full*full+k*full+l;
        D2bb[id] = val;
    }while(addr_bb != end_bb);

    // close files
    psio->close(PSIF_V2RDM_D2AA,1);
    psio->close(PSIF_V2RDM_D2BB,1);
    psio->close(PSIF_V2RDM_D2AB,1);
}*/

}} //end namespaces


