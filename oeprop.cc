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

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <utility>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libciomr/libciomr.h>
#include <libmints/vector3.h>
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <physconst.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "oeprop.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi {
MyOEProp::MyOEProp() : OEProp()
{
    common_init();
}
MyOEProp::~MyOEProp()
{
}

SharedMatrix MyOEProp::Da_ao_custom(boost::shared_ptr<Matrix>Da_so)
{
    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    SharedMatrix D = SharedMatrix(new Matrix("Da (AO basis)", basisset_->nbf(), basisset_->nbf()));
    int symm = Da_so->symmetry();
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h^symm];
        if (!nsol || !nsor) continue;
        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h^symm);
        double** DSOp = Da_so->pointer(h^symm);
        double** DAOp = D->pointer();
        C_DGEMM('N','T',nsol,nao,nsor,1.0,DSOp[0],nsor,Urp[0],nsor,0.0,temp,nao);
        C_DGEMM('N','N',nao,nao,nsol,1.0,Ulp[0],nsol,temp,nao,1.0,DAOp[0],nao);
    }
    delete[] temp;
    return D;
}

void MyOEProp::compute_mulliken_charges_custom(boost::shared_ptr<Matrix>Da_so,boost::shared_ptr<Matrix>Db_so)
{

    outfile->Printf( "\n");
    outfile->Printf( "  Mulliken Charges: (a.u.)\n");
    outfile->Printf( "\n");

    boost::shared_ptr<Molecule> mol = basisset_->molecule();

    double* PSa = new double[basisset_->nbf()];
    double suma = 0.0;

    double* PSb = new double[basisset_->nbf()];
    double sumb = 0.0;

    double* Qa = new double[mol->natom()];
    double* Qb = new double[mol->natom()];

    ::memset(Qa, '\0', mol->natom()*sizeof(double));
    ::memset(Qb, '\0', mol->natom()*sizeof(double));

//    Compute the overlap matrix

    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);


    SharedMatrix Da = Da_ao_custom(Da_so);
    SharedMatrix Db = Da_ao_custom(Db_so);

//    Form the idempotent D*S matrix (well, its not idempotent any more)

    SharedMatrix PSam(new Matrix("PSa",basisset_->nbf(),basisset_->nbf()));
    PSam->gemm(false,false,1.0,Da,S,0.0);
    SharedMatrix PSbm(new Matrix("PSb",basisset_->nbf(),basisset_->nbf()));
    PSbm->gemm(false,false,1.0,Db,S,0.0);

//     Accumulate registers

    for (int mu = 0; mu < basisset_->nbf(); mu++) {
        PSa[mu] = PSam->get(0,mu,mu);
        PSb[mu] = PSbm->get(0,mu,mu);

        int shell = basisset_->function_to_shell(mu);
        int A = basisset_->shell_to_center(shell);

        Qa[A] += PSa[mu];
        Qb[A] += PSb[mu];

        suma += PSa[mu];
        sumb += PSb[mu];
    }

//    Print out the Mulliken populations and charges

    outfile->Printf( "   Center  Symbol    Alpha    Beta     Spin     Total\n");
    double nuc = 0.0;
    for (int A = 0; A < mol->natom(); A++) {
        double Qs = Qa[A] - Qb[A];
        double Qt = mol->Z(A) - (Qa[A] + Qb[A]);
        outfile->Printf("   %5d    %2s    %8.5f %8.5f %8.5f %8.5f\n", A+1,mol->label(A).c_str(), \
            Qa[A], Qb[A], Qs, Qt);
        nuc += (double) mol->Z(A);
   }

    outfile->Printf( "\n   Total alpha = %8.5f, Total beta = %8.5f, Total charge = %8.5f\n", \
        suma, sumb, nuc - suma - sumb);

//    Free memory
    delete[] PSa;
    delete[] PSb;
    delete[] Qa;
    delete[] Qb;
}

}
