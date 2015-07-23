/*
 *@BEGIN LICENSE
 *
 * v2rdm_casscf by 
 *
 *   Jacob Fosso-Tande
 *   Greg Gidofalvi
 *   Eugene DePrince
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
 *@END LICENSE
 */

/*                                                                      
 *@BEGIN LICENSE

  Copyright (c) 2014, The Florida State University. All rights reserved.

 *@END LICENSE
 */


#include "v2rdm_solver.h"

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include<libciomr/libciomr.h>

#include<../bin/fnocc/frozen_natural_orbitals.h>

INIT_PLUGIN

using namespace boost;
using namespace psi;
using namespace fnocc;

namespace psi{ namespace v2rdm_casscf {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "V2RDM_CASSCF"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- The multiplicity.  This option will override the multiplicity 
            from the molecule group -*/
        options.add_int("MULTIPLICITY", 1);
        /*- The type of 2-positivity computation -*/
        options.add_str("POSITIVITY", "DQG", "DQG D DQ DG DQGT1 DQGT2 DQGT1T2");
        /*- Do constrain D3 to D2 mapping? -*/
        options.add_bool("CONSTRAIN_D3",false);
        /*- Do spin adapt G2 condition? -*/
        options.add_bool("SPIN_ADAPT_G2", false);
        /*- Do spin adapt Q2 condition? -*/
        options.add_bool("SPIN_ADAPT_Q2", false);
        /*- Do spin adapt D2? -*/
        //options.add_bool("SPIN_ADAPT_D2", false);
        /*- Do constrain spin squared? -*/
        options.add_bool("CONSTRAIN_SPIN", true);
        /*- convergence in the primal/dual energy gap -*/
        options.add_double("E_CONVERGENCE", 1e-4);
        /*- convergence in the primal error -*/
        options.add_double("R_CONVERGENCE", 1e-3);
        /*- convergence for conjugate gradient solver -*/
        options.add_double("CG_CONVERGENCE", 1e-5);
        /*- maximum number of bpsdp outer iterations -*/
        options.add_int("MAXITER", 10000);
        /*- maximum number of conjugate gradient iterations -*/
        options.add_int("CG_MAXITER", 10000);

        /*- SUBSECTION JACOBI -*/

        /*- number of threads to use for jacobi rotations -*/
        options.add_int("JACOBI_NTHREAD",1);
        /*- number of truly frozen orbitals (not optimized by jacobi) -*/
        options.add_int("JACOBI_FROZEN_CORE",0);
        /*- do rotate active/active orbital pairs? -*/
        options.add_bool("JACOBI_ACTIVE_ACTIVE_ROTATIONS",false);
        /*- tolerance for energy change for a given pair of orbitals -*/
        options.add_double("JACOBI_ANGLE_TOLERANCE",1.0e-6);
        /*- convergence in energy for rotations -*/
        options.add_double("JACOBI_E_CONVERGENCE",1.0e-8);
        /*- Do write a MOLDEN output file?  If so, the filename will end in
        .molden, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("MOLDEN_WRITE", false);
        /*- Do write a JACOBI output file?  If so, the filename will end in
        .molden, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("JACOBI_WRITE", false);
        /*- Base filename for text files written by PSI, such as the
        MOLDEN output file, the Hessian file, the internal coordinate file,
        etc. Use the add_str_i function to make this string case sensitive. -*/
        options.add_str_i("WRITER_FILE_LABEL", "v2rdm_casscf");
    }

    return true;
}

extern "C" 
PsiReturnType v2rdm_casscf(Options& options)
{

    tstart();

    boost::shared_ptr<Wavefunction> wfn;

    if ( options.get_str("SCF_TYPE") == "DF" || options.get_str("SCF_TYPE") == "CD") { 

        // get three-index integrals in usable form
        boost::shared_ptr<DFFrozenNO> fno(new DFFrozenNO(Process::environment.wavefunction(),options));
        fno->ThreeIndexIntegrals();
        wfn = (boost::shared_ptr<Wavefunction>)fno;

    }else {
        wfn = Process::environment.wavefunction();
    }

    boost::shared_ptr<v2RDMSolver > v2rdm (new v2RDMSolver(wfn,options));
    double energy = v2rdm->compute_energy();
    Process::environment.globals["CURRENT ENERGY"] = energy;

    tstop();

    return Success;
}

}} // End namespaces

