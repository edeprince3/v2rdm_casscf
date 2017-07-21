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

#include"v2rdm_solver.h"
#include<psi4/libqt/qt.h>
#include<psi4/psifiles.h>
#include "blas.h"

using namespace psi;
using namespace fnocc;

/*================================================================

   diis functions

================================================================*/

namespace psi{ namespace v2rdm_casscf{

// TODO: replace junk1/junk2
void v2RDMSolver::DIIS(double*c,long int nvec,int replace_diis_iter){
    long int nvar      = nvec+1;
    long int * ipiv    = (long int*)malloc(nvar*sizeof(long int));
    double * temp      = (double*)malloc(sizeof(double)*maxdiis_*maxdiis_);
    double * A         = (double*)malloc(sizeof(double)*nvar*nvar);
    double * B         = (double*)malloc(sizeof(double)*nvar);
    memset((void*)A,'\0',nvar*nvar*sizeof(double));
    memset((void*)B,'\0',nvar*sizeof(double));
    B[nvec] = -1.;

    char*evector=(char*)malloc(1000*sizeof(char));

    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_EVEC,PSIO_OPEN_OLD);

    // add row to matrix, don't build the whole thing.
    psio->read_entry(PSIF_DCC_EVEC,"error matrix",(char*)&temp[0],maxdiis_*maxdiis_*sizeof(double));
    for (long int i = 0; i < nvec; i++){
        for (long int j = 0; j < nvec; j++){
            A[i*nvar+j] = temp[i*maxdiis_+j];
        }
    }

    if (nvec <= 3) {
        for (long int i = 0; i < nvec; i++) {
            sprintf(evector,"evector%li",i+1);
            psio->read_entry(PSIF_DCC_EVEC,evector,(char*)&junk1[0],dimdiis_*sizeof(double));
            for (long int j = i; j < nvec; j++){
                sprintf(evector,"evector%li",j+1);
                psio->read_entry(PSIF_DCC_EVEC,evector,(char*)&junk2[0],dimdiis_*sizeof(double));
                double sum  = C_DDOT(dimdiis_,junk1,1,junk2,1);
                A[i*nvar+j] = sum;
                A[j*nvar+i] = sum;
            }
        }
    }else {
        long int i;
        if (nvec<=maxdiis_ && diis_oiter_<=maxdiis_){
            i = nvec - 1;
        }
        else{
            i = replace_diis_iter - 1;
        }
        sprintf(evector,"evector%li",i+1);
        psio->read_entry(PSIF_DCC_EVEC,evector,(char*)&junk1[0],dimdiis_*sizeof(double));
        for (long int j = 0; j < nvec; j++){
            sprintf(evector,"evector%li",j+1);
            psio->read_entry(PSIF_DCC_EVEC,evector,(char*)&junk2[0],dimdiis_*sizeof(double));
            double sum  = C_DDOT(dimdiis_,junk1,1,junk2,1);
            A[i*nvar+j] = sum;
            A[j*nvar+i] = sum;
        }
    }

    long int j = nvec;
    for (long int i = 0; i < nvar; i++){
        A[j*nvar+i] = -1.0;
        A[i*nvar+j] = -1.0;
    }
    A[nvar*nvar-1] = 0.;

    // save matrix for next iteration
    for (long int i = 0; i < nvec; i++){
        for (long int j = 0; j < nvec; j++){
            temp[i*maxdiis_+j] = A[i*nvar+j];
        }
    }
    psio->write_entry(PSIF_DCC_EVEC,"error matrix",(char*)&temp[0],maxdiis_*maxdiis_*sizeof(double));
    free(temp);
    psio->close(PSIF_DCC_EVEC,1);
    free(evector);

    long int nrhs,lda,ldb,info;
    nrhs = 1;
    lda = ldb = nvar;
    info = 0;
    DGESV(nvar,nrhs,A,lda,ipiv,B,ldb, info);
    C_DCOPY(nvec,B,1,c,1);

    free(A);
    free(B);
    free(ipiv);
    psio.reset();
}

void v2RDMSolver::DIIS_WriteOldVector(long int iter,int diis_iter,int replace_diis_iter){

    char*oldvector=(char*)malloc(1000*sizeof(char));

    if (diis_iter<=maxdiis_ && iter<=maxdiis_){
       sprintf(oldvector,"oldvector%i",diis_iter);
    }
    else{
       sprintf(oldvector,"oldvector%i",replace_diis_iter);
    }

    std::shared_ptr<PSIO> psio(new PSIO());
    if (diis_iter==0) {
       psio->open(PSIF_DCC_OVEC,PSIO_OPEN_NEW);
    }else {
       psio->open(PSIF_DCC_OVEC,PSIO_OPEN_OLD);
    }

    psio_address addr;
    addr = PSIO_ZERO;

    psio->write(PSIF_DCC_OVEC,oldvector,(char*)&rx->pointer()[0],dimdiis_*sizeof(double),addr,&addr);
    psio->write(PSIF_DCC_OVEC,oldvector,(char*)&rz->pointer()[0],dimdiis_*sizeof(double),addr,&addr);
    psio->close(PSIF_DCC_OVEC,1);
    psio.reset();

    free(oldvector);
}
void v2RDMSolver::DIIS_WriteErrorVector(int diis_iter,int replace_diis_iter,int iter){

    char*evector   = (char*)malloc(1000*sizeof(char));
    if (diis_iter<=maxdiis_ && iter<=maxdiis_){
       sprintf(evector,"evector%i",diis_iter);
    }
    else{
       sprintf(evector,"evector%i",replace_diis_iter);
    }

    std::shared_ptr<PSIO> psio(new PSIO());
    if (diis_iter==0) {
       psio->open(PSIF_DCC_EVEC,PSIO_OPEN_NEW);
       double * temp = (double*)malloc(maxdiis_*maxdiis_*sizeof(double));
       memset((void*)temp,'\0',maxdiis_*maxdiis_*sizeof(double));
       psio->write_entry(PSIF_DCC_EVEC,"error matrix",(char*)&temp[0],maxdiis_*maxdiis_*sizeof(double));
       free(temp);
    }
    else {
       psio->open(PSIF_DCC_EVEC,PSIO_OPEN_OLD);
    }

    psio_address addr;
    addr = PSIO_ZERO;

    psio->write(PSIF_DCC_EVEC,evector,(char*)&rx_error->pointer()[0],dimdiis_*sizeof(double),addr,&addr);
    psio->write(PSIF_DCC_EVEC,evector,(char*)&rz_error->pointer()[0],dimdiis_*sizeof(double),addr,&addr);

    psio->close(PSIF_DCC_EVEC,1);
    psio.reset();

    free(evector);
}
void v2RDMSolver::DIIS_Extrapolate(int diis_iter,int&replace_diis_iter){

    char*oldvector;
    oldvector=(char*)malloc(1000*sizeof(char));

    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_OVEC,PSIO_OPEN_OLD);

    psio_address addr;

    memset((void*)rx->pointer(),'\0',dimdiis_*sizeof(double));
    memset((void*)rz->pointer(),'\0',dimdiis_*sizeof(double));

    long int max = diis_iter;
    if (max > maxdiis_) max = maxdiis_;

    double min = 1.e9;
    for (long int j=1; j<=max; j++){
        addr = PSIO_ZERO;
        sprintf(oldvector,"oldvector%li",j);

        psio->read(PSIF_DCC_OVEC,oldvector,(char*)&junk1[0],dimdiis_*sizeof(double),addr,&addr);
        C_DAXPY(dimdiis_,diisvec_[j-1],junk1,1,rx->pointer(),1);

        psio->read(PSIF_DCC_OVEC,oldvector,(char*)&junk1[0],dimdiis_*sizeof(double),addr,&addr);
        C_DAXPY(dimdiis_,diisvec_[j-1],junk1,1,rz->pointer(),1);

        //if ( fabs( diisvec[j-1] ) < min ) {
        //    min = fabs( diisvec[j-1] );
        //    replace_diis_iter = j;
        //}
    }
    psio->close(PSIF_DCC_OVEC,1);
    free(oldvector);

    // now, build x = rx^2, z = rz^2
    double * x_p  = x->pointer();
    double * z_p  = z->pointer();
    double * rx_p = rx->pointer();
    double * rz_p = rz->pointer();

    // loop over each block of x/z
    for (int i = 0; i < dimensions_.size(); i++) {
        if ( dimensions_[i] == 0 ) continue;
        int myoffset = 0;
        for (int j = 0; j < i; j++) {
            myoffset += dimensions_[j] * dimensions_[j];
        }
        F_DGEMM('n','n',dimensions_[i],dimensions_[i],dimensions_[i],1.0,rx_p+myoffset,dimensions_[i],rx_p+myoffset,dimensions_[i],0.0,x_p+myoffset,dimensions_[i]);
        F_DGEMM('n','n',dimensions_[i],dimensions_[i],dimensions_[i],1.0,rz_p+myoffset,dimensions_[i],rz_p+myoffset,dimensions_[i],0.0,z_p+myoffset,dimensions_[i]);
    }

    psio.reset();
}

}}
