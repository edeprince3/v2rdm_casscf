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

#include <psi4/psi4-dec.h>
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>

#include "v2rdm_solver.h"

using namespace psi;



namespace psi{namespace v2rdm_casscf{

void ReadAllIntegrals(iwlbuf *Buf,double*tei,long int nmo) {

  unsigned long int lastbuf;
  Label *lblptr;
  Value *valptr;

  unsigned long int idx, p, q, r, s, pq, rs, pqrs;

  lblptr = Buf->labels;
  valptr = Buf->values;
  lastbuf = Buf->lastbuf;

  outfile->Printf("\n");
  outfile->Printf("        Read integrals......");
  /**
    * first buffer (read in when Buf was initialized)
    */
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (unsigned long int) lblptr[idx++];
      q = (unsigned long int) lblptr[idx++];
      r = (unsigned long int) lblptr[idx++];
      s = (unsigned long int) lblptr[idx++];

      pq   = INDEX(p,q);
      rs   = INDEX(r,s);
      pqrs = INDEX(pq,rs);

      double val = (double)valptr[Buf->idx];

      tei[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] = val;
      tei[p*nmo*nmo*nmo+q*nmo*nmo+s*nmo+r] = val;
      tei[q*nmo*nmo*nmo+p*nmo*nmo+r*nmo+s] = val;
      tei[q*nmo*nmo*nmo+p*nmo*nmo+s*nmo+r] = val;

      tei[r*nmo*nmo*nmo+s*nmo*nmo+p*nmo+q] = val;
      tei[s*nmo*nmo*nmo+r*nmo*nmo+p*nmo+q] = val;
      tei[r*nmo*nmo*nmo+s*nmo*nmo+q*nmo+p] = val;
      tei[s*nmo*nmo*nmo+r*nmo*nmo+q*nmo+p] = val;
  }
  /**
    * now do the same for the rest of the buffers
    */
  while(!lastbuf){
      iwl_buf_fetch(Buf);
      lastbuf = Buf->lastbuf;
      for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {

          p = (unsigned long int) lblptr[idx++];
          q = (unsigned long int) lblptr[idx++];
          r = (unsigned long int) lblptr[idx++];
          s = (unsigned long int) lblptr[idx++];

          pq   = INDEX(p,q);
          rs   = INDEX(r,s);
          pqrs = INDEX(pq,rs);

          double val = (double)valptr[Buf->idx];

          tei[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] = val;
          tei[p*nmo*nmo*nmo+q*nmo*nmo+s*nmo+r] = val;
          tei[q*nmo*nmo*nmo+p*nmo*nmo+r*nmo+s] = val;
          tei[q*nmo*nmo*nmo+p*nmo*nmo+s*nmo+r] = val;

          tei[r*nmo*nmo*nmo+s*nmo*nmo+p*nmo+q] = val;
          tei[s*nmo*nmo*nmo+r*nmo*nmo+p*nmo+q] = val;
          tei[r*nmo*nmo*nmo+s*nmo*nmo+q*nmo+p] = val;
          tei[s*nmo*nmo*nmo+r*nmo*nmo+q*nmo+p] = val;

      }

  }
  outfile->Printf("done.\n\n");
}


void v2RDMSolver::ReadIntegrals(double * tei,long int nmo){
  struct iwlbuf Buf;
  iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);
  ReadAllIntegrals(&Buf,tei,nmo);
  iwl_buf_close(&Buf,1);
}


}}
