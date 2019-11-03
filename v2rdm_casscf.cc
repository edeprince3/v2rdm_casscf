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

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include<psi4/libciomr/libciomr.h>
#include <psi4/libpsi4util/process.h>

#include"backtransform_tpdm.h"

INIT_PLUGIN

using namespace psi;
//using namespace fnocc;

namespace psi{ namespace v2rdm_casscf {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "V2RDM_CASSCF"|| options.read_globals()) {
        /*- fractional charge -*/
        options.add_double("FRACTIONAL_CHARGE", 0.0);
        /*- do extended koopmans theorem computation? -*/
        options.add_bool("EXTENDED_KOOPMANS",false);
        /*- !expert -*/
        options.add_bool("DOCI",false);
        /*- !expert -*/
        options.add_double("DOCI_ALPHA", 0.0);
        /*- Do v2RDM-CASSCF gradient? !expert -*/
        /* Do write fcidump files? -*/
        options.add_bool("FCIDUMP", false);
        /* Do v2RDM-CASSCF gradient? !expert -*/
        options.add_str("DERTYPE", "NONE", "NONE FIRST");
        /*- Do localize orbitals? -*/
        options.add_bool("LOCALIZE_ORBITALS",false);
        /*- Do optimize orbitals? -*/
        options.add_bool("OPTIMIZE_ORBITALS",true);
        /*- Rotate guess orbitals -*/
        options.add("MCSCF_ROTATE", new ArrayType());
		/*- Do semicanonicalize orbitals? -*/
        options.add_bool("SEMICANONICALIZE_ORBITALS",false);
        /*- Type of guess -*/
        options.add_str("TPDM_GUESS","RANDOM", "RANDOM HF");
        /*- Do compute natural orbitals and transform 1- and 2-RDM to the natural orbital basis? 
        The OPDM and Ca/Cb matrices pushed onto the wavefunction will correspond to the natural orbital basis -*/
        options.add_bool("NAT_ORBS",false);
        /*- Do write the 1-RDM to disk? All nonzero elements of the 1-RDM will be written.  -*/
        options.add_bool("OPDM_WRITE_FULL",false);
        /*- Do write the 2-RDM to disk? All nonzero elements of the 2-RDM will be written.  -*/
        options.add_bool("TPDM_WRITE_FULL",false);
        /*- Do write the 2-RDM to disk? Only the nonzero elements of the active 2-RDM will be written. -*/
        options.add_bool("TPDM_WRITE",false);
        /*- Do write the 3-RDM to disk? -*/
        options.add_bool("3PDM_WRITE",false);
        /*- Do save progress in a checkpoint file? -*/
        options.add_bool("WRITE_CHECKPOINT_FILE",false);
        /*- Frequency of checkpoint file generation.  The checkpoint file is 
        updated every CHECKPOINT_FREQUENCY iterations.  The default frequency
        will be ORBOPT_FREQUENCY. -*/
        options.add_int("CHECKPOINT_FREQUENCY",500);
        /*- File containing previous primal/dual solutions and integrals. -*/
        options.add_str("RESTART_FROM_CHECKPOINT_FILE","");
        /*- A parameter introduced by Mazziotti [PRL 106, 083001 (2011)] to "increase the
        sensitivity of y on the deviation of x from primal feasibility."  Should 
        lie on the interval [1.0, 1.6]. -*/
        options.add_double("TAU_PARAMETER",1.0);
        /*- Frequency with which the pentalty-parameter, mu, is updated. mu is
        updated every MU_UPDATE_FREQUENCY iterations.   -*/
        options.add_int("MU_UPDATE_FREQUENCY",500);
        /*- The type of 2-positivity computation -*/
        options.add_str("POSITIVITY", "DQG", "DQG D DQ DG DQGT1 DQGT2 DQGT1T2");
        /*- Do constrain D4 to D3 mapping? -*/
        options.add_bool("CONSTRAIN_D4",false);
        /*- Do constrain D3 to D2 mapping? -*/
        options.add_bool("CONSTRAIN_D3",false);
        /*- Do spin adapt G2 condition? -*/
        options.add_bool("SPIN_ADAPT_G2", false);
        /*- Do spin adapt Q2 condition? -*/
        options.add_bool("SPIN_ADAPT_Q2", false);
        /*- Do constrain spin squared? -*/
        options.add_bool("CONSTRAIN_SPIN", true);
        /*- Do constrain sz? -*/
        options.add_bool("CONSTRAIN_SZ", true);
        /*- convergence in the primal/dual energy gap -*/
        options.add_double("E_CONVERGENCE", 1e-4);
        /*- convergence in the primal error -*/
        options.add_double("R_CONVERGENCE", 1e-4);
        /*- convergence for conjugate gradient solver. currently not used. -*/
        options.add_double("CG_CONVERGENCE", 1e-9);
        /*- maximum number of bpsdp outer iterations -*/
        options.add_int("MAXITER", 10000);
        /*- maximum number of conjugate gradient iterations -*/
        options.add_int("CG_MAXITER", 10000);
        /*- maximum number of diis vectors -*/
        options.add_int("DIIS_MAX_VECS", 8);
        /*- Frequency of DIIS extrapolation steps -*/
        options.add_int("DIIS_UPDATE_FREQUENCY",50);

        /*- Auxiliary basis set for SCF density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a JKFIT basis. -*/
        options.add_str("DF_BASIS_SCF", "");
        /*- What algorithm to use for the SCF computation. See Table :ref:`SCF
        Convergence & Algorithm <table:conv_scf>` for default algorithm for
        different calculation types. -*/
        options.add_str("SCF_TYPE", "DF", "DF CD PK OUT_OF_CORE DIRECT");
        /*- Tolerance for Cholesky decomposition of the ERI tensor -*/
        options.add_double("CHOLESKY_TOLERANCE",1e-4);

        /*- SUBSECTION ORBITAL OPTIMIZATION -*/

        /*- algorithm for orbital optimization -*/
        options.add_str("ORBOPT_ALGORITHM","QUASI_NEWTON", "QUASI_NEWTON CONJUGATE_GRADIENT NEWTON_RAPHSON");
        /*- flag to optimize orbitals using a one-step type approach -*/
        options.add_bool("ORBOPT_ONE_STEP",true);
        /*- do rotate active/active orbital pairs? -*/
        options.add_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS",false);
        /*- convergence in gradient norm -*/
        options.add_double("ORBOPT_GRADIENT_CONVERGENCE",1.0e-6);
        /*- convergence in energy for rotations -*/
        options.add_double("ORBOPT_ENERGY_CONVERGENCE",1.0e-8);
        /*- flag for using exact expresions for diagonal Hessian element -*/
        options.add_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN",false);
        /*- number of DIIS vectors to keep in orbital optimization -*/ 
        options.add_int("ORBOPT_NUM_DIIS_VECTORS",0);
        /*- frequency of orbital optimization.  optimization occurs every 
        orbopt_frequency iterations -*/
        options.add_int("ORBOPT_FREQUENCY",500);
        /*- maximum number of iterations for orbital optimization -*/
        options.add_int("ORBOPT_MAXITER",20);
        /*- Do write a MOLDEN output file?  If so, the filename will end in
        .molden, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("MOLDEN_WRITE", false);
        /*- Do write a MOLDEN file for guess orbitals?  If so, the filename will
        end in .guess.molden, and the prefix is determined by 
        |globals__writer_file_label| (if set), or else by the name of the output
        file plus the name of the current molecule. -*/
        options.add_bool("GUESS_ORBITALS_WRITE", false);
        /*- Do write a ORBOPT output file?  If so, the filename will end in
        .molden, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("ORBOPT_WRITE", false);
        /*- Base filename for text files written by PSI, such as the
        MOLDEN output file, the Hessian file, the internal coordinate file,
        etc. Use the add_str_i function to make this string case sensitive. -*/
        options.add_str_i("WRITER_FILE_LABEL", "v2rdm_casscf");
    }

    return true;
}

extern "C" PSI_API
SharedWavefunction v2rdm_casscf(SharedWavefunction ref_wfn, Options& options)
{
    tstart();

    std::shared_ptr<v2RDMSolver > v2rdm (new v2RDMSolver(ref_wfn,options));

    double energy = v2rdm->compute_energy();

    Process::environment.globals["CURRENT ENERGY"] = energy;

    if ( options.get_str("DERTYPE") == "FIRST" ) {
        // backtransform the tpdm
        std::vector<std::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::all);
        std::shared_ptr<TPDMBackTransform> transform = std::shared_ptr<TPDMBackTransform>(
        new TPDMBackTransform(ref_wfn,
                        spaces,
                        IntegralTransform::TransformationType::Unrestricted, // Transformation type
                        IntegralTransform::OutputType::DPDOnly,              // Output buffer
                        IntegralTransform::MOOrdering::QTOrder,              // MO ordering
                        IntegralTransform::FrozenOrbitals::None));           // Frozen orbitals?
        transform->backtransform_density();
        transform.reset();
    }

    tstop();

    return (SharedWavefunction)v2rdm;
}

}} // End namespaces

