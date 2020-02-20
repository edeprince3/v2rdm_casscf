#
#@BEGIN LICENSE
#
# v2rdm_casscf by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

import re
import os
import math
import warnings
import numpy

import psi4
from psi4.driver.procrouting import proc_util
import psi4.driver.p4util as p4util

from psi4.driver.p4util import solvers

import v2rdm_casscf

def run_v2rdm_casscf(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    v2rdm_casscf can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('v2rdm_casscf')

    """

    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # Your plugin's psi4 run sequence goes here
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # if restarting from a checkpoint file, this file
    # needs to be in scratch with the correct name
    filename = psi4.core.get_option("V2RDM_CASSCF","RESTART_FROM_CHECKPOINT_FILE")

    # Ensure IWL files have been written when not using DF/CD
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # reorder wavefuntions based on user input
    # apply a list of 2x2 rotation matrices to the orbitals in the form of [irrep, orbital1, orbital2, theta]
    # where an angle of 0 would do nothing and an angle of 90 would switch the two orbitals.
    # the indices of irreps and orbitals start from 0
    reorder_orbitals = psi4.core.get_option("V2RDM_CASSCF","MCSCF_ROTATE")
    for orbord in reorder_orbitals:
        if type(orbord) != list :
            raise psi4.p4util.PsiException("Each element of the orbtial rotate list requires 4 arguements (irrep, orb1, orb2, theta).")
        if len(orbord) != 4:
            raise psi4.p4util.PsiException("Each element of the orbtial rotate list requires 4 arguements (irrep, orb1, orb2, theta).")

        irrep, orb1, orb2, theta = orbord

        if irrep > ref_wfn.Ca().nirrep():
            raise psi4.p4util.PsiException("REORDER_ORBITALS: Expression %s irrep number is larger than the number of irreps" %
                                    (str(orbord)))

        if max(orb1, orb2) > ref_wfn.Ca().coldim()[irrep]:
            raise psi4.p4util.PsiException("REORDER_ORBITALS: Expression %s orbital number exceeds number of orbitals in irrep" %
                                    (str(orbord)))

        theta = numpy.deg2rad(theta)

        x_a = ref_wfn.Ca().nph[irrep][:, orb1].copy()
        y_a = ref_wfn.Ca().nph[irrep][:, orb2].copy()

        xp_a = numpy.cos(theta) * x_a - numpy.sin(theta) * y_a
        yp_a = numpy.sin(theta) * x_a + numpy.cos(theta) * y_a

        ref_wfn.Ca().nph[irrep][:, orb1] = xp_a
        ref_wfn.Ca().nph[irrep][:, orb2] = yp_a

        x_b = ref_wfn.Ca().nph[irrep][:, orb1].copy()
        y_b = ref_wfn.Ca().nph[irrep][:, orb2].copy()

        xp_b = numpy.cos(theta) * x_b - numpy.sin(theta) * y_b
        yp_b = numpy.sin(theta) * x_b + numpy.cos(theta) * y_b

        ref_wfn.Ca().nph[irrep][:, orb1] = xp_b
        ref_wfn.Ca().nph[irrep][:, orb2] = yp_b


    returnvalue = psi4.core.plugin('v2rdm_casscf.so', ref_wfn)

    optstash.restore()

    return returnvalue

def run_v2rdm_casscf_gradient(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    v2rdm_casscf can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> gradient('v2rdm_casscf')

    """

    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'],
        ['V2RDM_CASSCF', 'OPTIMIZE_ORBITALS'],
        ['V2RDM_CASSCF', 'SEMICANONICALIZE_ORBITALS'],
        ['V2RDM_CASSCF', 'ORBOPT_ACTIVE_ACTIVE_ROTATIONS'],
        ['V2RDM_CASSCF', 'RESTART_FROM_CHECKPOINT_FILE'],
        ['V2RDM_CASSCF', 'WRITE_CHECKPOINT_FILE'])

    psi4.core.set_global_option('DERTYPE', 'FIRST')
    psi4.core.set_local_option("V2RDM_CASSCF","OPTIMIZE_ORBITALS",True)
    psi4.core.set_local_option("V2RDM_CASSCF","ORBOPT_ACTIVE_ACTIVE_ROTATIONS",True)
    psi4.core.set_local_option("V2RDM_CASSCF","SEMICANONICALIZE_ORBITALS",False)
    psi4.core.set_local_option("V2RDM_CASSCF","RESTART_FROM_CHECKPOINT_FILE","DUMMY")
    psi4.core.set_local_option("V2RDM_CASSCF","WRITE_CHECKPOINT_FILE",True)

    # analytic derivatives do not work with scf_type df/cd
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'CD' or scf_type == 'DF' ):
        raise ValidationError("""Error: analytic v2RDM-CASSCF gradients not implemented for scf_type %s.""" % scf_type)


    v2rdm_wfn = run_v2rdm_casscf(name,**kwargs)

    derivobj = psi4.core.Deriv(v2rdm_wfn)
    derivobj.set_deriv_density_backtransformed(True)
    derivobj.set_ignore_reference(True)
    grad = derivobj.compute()

    v2rdm_wfn.set_gradient(grad)

    optstash.restore()

    return v2rdm_wfn

# Integration with driver routines
psi4.driver.procedures['energy']['v2rdm-casscf'] = run_v2rdm_casscf
psi4.driver.procedures['gradient']['v2rdm-casscf'] = run_v2rdm_casscf_gradient

def exampleFN():
    # Your Python code goes here
    pass

def test_sparse(ref_wfn):

    options = psi4.core.get_options()
    options.set_current_module('V2RDM_CASSCF')

    v2rdm = v2rdm_casscf.v2RDMHelper(ref_wfn,options)
    current_energy = v2rdm.compute_energy()

    d1 = v2rdm.get_opdm_sparse("SUM")
    d2 = v2rdm.get_tpdm_sparse("SUM")

    for i in range(0,len(d1)):
        print("D1(%i,%i) = %f" % (d1[i].i, d1[i].j, d1[i].value))
    for i in range(0,len(d2)):
        print("D2(%i, %i; %i ,%i) = %f" % (d2[i].i, d2[i].j, d2[i].k, d2[i].l, d2[i].value))



def print_iteration(mtype, niter, energy, de, orb_rms, ci_rms, nci, norb, stype):
    psi4.core.print_out("%s %2d:  % 18.12f   % 1.4e  %1.2e  %1.2e  %3d  %3d  %s\n" %
                    (mtype, niter, energy, de, orb_rms, ci_rms, nci, norb, stype))

def v2rdm_scf_solver(ref_wfn):

    # AED
    psi4.core.set_local_option('DETCI', 'WFN', 'CASSCF')

    # Build CIWavefunction
    psi4.core.prepare_options_for_module("DETCI")
    ciwfn = psi4.core.CIWavefunction(ref_wfn)

    # Hush a lot of CI output
    ciwfn.set_print(0)

    # Begin with a normal two-step
    step_type = 'Initial CI'
    total_step = psi4.core.Matrix("Total step", ciwfn.get_dimension('OA'), ciwfn.get_dimension('AV'))
    start_orbs = ciwfn.get_orbitals("ROT").clone()
    ciwfn.set_orbitals("ROT", start_orbs)

    # Grab da options
    mcscf_orb_grad_conv = psi4.core.get_option("DETCI", "MCSCF_R_CONVERGENCE")
    mcscf_e_conv = psi4.core.get_option("DETCI", "MCSCF_E_CONVERGENCE")
    mcscf_max_macroiteration = psi4.core.get_option("DETCI", "MCSCF_MAXITER")
    mcscf_type = psi4.core.get_option("DETCI", "MCSCF_TYPE")
    mcscf_d_file = psi4.core.get_option("DETCI", "CI_FILE_START") + 3
    mcscf_nroots = psi4.core.get_option("DETCI", "NUM_ROOTS")
    mcscf_wavefunction_type = psi4.core.get_option("DETCI", "WFN")
    mcscf_ndet = ciwfn.ndet()
    mcscf_nuclear_energy = ciwfn.molecule().nuclear_repulsion_energy()
    mcscf_steplimit = psi4.core.get_option("DETCI", "MCSCF_MAX_ROT")
    mcscf_rotate = psi4.core.get_option("DETCI", "MCSCF_ROTATE")

    # DIIS info
    mcscf_diis_start = psi4.core.get_option("DETCI", "MCSCF_DIIS_START")
    mcscf_diis_freq = psi4.core.get_option("DETCI", "MCSCF_DIIS_FREQ")
    mcscf_diis_error_type = psi4.core.get_option("DETCI", "MCSCF_DIIS_ERROR_TYPE")
    mcscf_diis_max_vecs = psi4.core.get_option("DETCI", "MCSCF_DIIS_MAX_VECS")

    # One-step info
    mcscf_target_conv_type = psi4.core.get_option("DETCI", "MCSCF_ALGORITHM")
    mcscf_so_start_grad = psi4.core.get_option("DETCI", "MCSCF_SO_START_GRAD")
    mcscf_so_start_e = psi4.core.get_option("DETCI", "MCSCF_SO_START_E")
    mcscf_current_step_type = 'Initial CI'

    # Start with SCF energy and other params
    scf_energy = ciwfn.variable("HF TOTAL ENERGY")
    eold = scf_energy
    norb_iter = 1
    converged = False
    ah_step = False
    qc_step = False
    approx_integrals_only = True

    # Fake info to start with the initial diagonalization
    ediff = 1.e-4
    orb_grad_rms = 1.e-3

    # Grab needed objects
    diis_obj = solvers.DIIS(mcscf_diis_max_vecs)
    mcscf_obj = ciwfn.mcscf_object()

    # Execute the rotate command
    for rot in mcscf_rotate:
        if len(rot) != 4:
            raise p4util.PsiException("Each element of the MCSCF rotate command requires 4 arguements (irrep, orb1, orb2, theta).")

        irrep, orb1, orb2, theta = rot
        if irrep > ciwfn.Ca().nirrep():
            raise p4util.PsiException("MCSCF_ROTATE: Expression %s irrep number is larger than the number of irreps" %
                                    (str(rot)))

        if max(orb1, orb2) > ciwfn.Ca().coldim()[irrep]:
            raise p4util.PsiException("MCSCF_ROTATE: Expression %s orbital number exceeds number of orbitals in irrep" %
                                    (str(rot)))

        theta = np.deg2rad(theta)

        x = ciwfn.Ca().nph[irrep][:, orb1].copy()
        y = ciwfn.Ca().nph[irrep][:, orb2].copy()

        xp = np.cos(theta) * x - np.sin(theta) * y
        yp = np.sin(theta) * x + np.cos(theta) * y

        ciwfn.Ca().nph[irrep][:, orb1] = xp
        ciwfn.Ca().nph[irrep][:, orb2] = yp

    # Limited RAS functionality
    if psi4.core.get_local_option("DETCI", "WFN") == "RASSCF" and mcscf_target_conv_type != "TS":
        psi4.core.print_out("\n  Warning! Only the TS algorithm for RASSCF wavefunction is currently supported.\n")
        psi4.core.print_out("             Switching to the TS algorithm.\n\n")
        mcscf_target_conv_type = "TS"

    # Print out headers
    if mcscf_type == "CONV":
        mtype = "   @MCSCF"
        psi4.core.print_out("\n   ==> Starting MCSCF iterations <==\n\n")
        psi4.core.print_out("        Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")
    elif mcscf_type == "DF":
        mtype = "   @DF-MCSCF"
        psi4.core.print_out("\n   ==> Starting DF-MCSCF iterations <==\n\n")
        psi4.core.print_out("           Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")
    else:
        mtype = "   @AO-MCSCF"
        psi4.core.print_out("\n   ==> Starting AO-MCSCF iterations <==\n\n")
        psi4.core.print_out("           Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")

    # Iterate !
    for mcscf_iter in range(1, mcscf_max_macroiteration + 1):

        ## Transform integrals, diagonalize H
        ciwfn.transform_mcscf_integrals(approx_integrals_only)
        #nci_iter = ciwfn.diag_h(abs(ediff) * 1.e-2, orb_grad_rms * 1.e-3)
        nci_iter = 0 #ciwfn.diag_h(abs(ediff) * 1.e-2, orb_grad_rms * 1.e-3)

        ## After the first diag we need to switch to READ
        #ciwfn.set_ci_guess("DFILE")

        #ciwfn.form_opdm()
        #ciwfn.form_tpdm()
        #ci_grad_rms = ciwfn.variable("DETCI AVG DVEC NORM")
        ci_grad_rms = 0.0

        # set options for v2RDM module (TODO: verify this is working correctly)
        psi4.core.set_local_option('V2RDM_CASSCF', 'OPTIMIZE_ORBITALS','FALSE')
        options = psi4.core.get_options()
        options.set_current_module('V2RDM_CASSCF')

        v2rdm = v2rdm_casscf.v2RDMHelper(ref_wfn,options)
        current_energy = v2rdm.compute_energy()
        opdm = v2rdm.get_opdm()
        tpdm = v2rdm.get_tpdm()

        Cocc = v2rdm.get_orbitals("DOCC")
        Cact = v2rdm.get_orbitals("ACTIVE")
        Cvir = v2rdm.get_orbitals("VIRTUAL")

        # END AED
        
        # Update MCSCF object
        #Cocc = ciwfn.get_orbitals("DOCC")
        #Cact = ciwfn.get_orbitals("ACT")
        #Cvir = ciwfn.get_orbitals("VIR")

        #opdm = ciwfn.get_opdm(-1, -1, "SUM", False)
        #tpdm = ciwfn.get_tpdm("SUM", True)

        Cact.print_out()
        mcscf_obj.update(Cocc, Cact, Cvir, opdm, tpdm)
        Cact.print_out()

        #current_energy = ciwfn.variable("v2RDM TOTAL ENERGY") #v2rdm.variable("v2RDM TOTAL ENERGY")

        orb_grad_rms = mcscf_obj.gradient_rms()
        ediff = current_energy - eold

        # Print iterations
        print_iteration(mtype, mcscf_iter, current_energy, ediff, orb_grad_rms, ci_grad_rms,
                        nci_iter, norb_iter, mcscf_current_step_type)
        eold = current_energy

        if mcscf_current_step_type == 'Initial CI':
            mcscf_current_step_type = 'TS'

        # Check convergence
        if (orb_grad_rms < mcscf_orb_grad_conv) and (abs(ediff) < abs(mcscf_e_conv)) and\
            (mcscf_iter > 3) and not qc_step:

            psi4.core.print_out("\n       %s has converged!\n\n" % mtype);
            converged = True
            break

        # Which orbital convergence are we doing?
        if ah_step:
            converged, norb_iter, step = ah_iteration(mcscf_obj, print_micro=False)
            norb_iter += 1

            if converged:
                mcscf_current_step_type = 'AH'
            else:
                psi4.core.print_out("      !Warning. Augmented Hessian did not converge. Taking an approx step.\n")
                step = mcscf_obj.approx_solve()
                mcscf_current_step_type = 'TS, AH failure'

        else:
            step = mcscf_obj.approx_solve()
            step_type = 'TS'

        maxstep = step.absmax()
        if maxstep > mcscf_steplimit:
            psi4.core.print_out('      Warning! Maxstep = %4.2f, scaling to %4.2f\n' % (maxstep, mcscf_steplimit))
            step.scale(mcscf_steplimit / maxstep)

        xstep = total_step.clone()
        total_step.add(step)

        # Do or add DIIS
        if (mcscf_iter >= mcscf_diis_start) and ("TS" in mcscf_current_step_type):

            # Figure out DIIS error vector
            if mcscf_diis_error_type == "GRAD":
                error = psi4.core.triplet(ciwfn.get_orbitals("OA"),
                                            mcscf_obj.gradient(),
                                            ciwfn.get_orbitals("AV"),
                                            False, False, True)
            else:
                error = step

            diis_obj.add(total_step, error)

            if not (mcscf_iter % mcscf_diis_freq):
                total_step = diis_obj.extrapolate()
                mcscf_current_step_type = 'TS, DIIS'

        # Build the rotation by continuous updates
        if mcscf_iter == 1:
            totalU = mcscf_obj.form_rotation_matrix(total_step)
        else:
            xstep.axpy(-1.0, total_step)
            xstep.scale(-1.0)
            Ustep = mcscf_obj.form_rotation_matrix(xstep)
            totalU = psi4.core.doublet(totalU, Ustep, False, False)

        # Build the rotation directly (not recommended)
        # orbs_mat = mcscf_obj.Ck(start_orbs, total_step)

        # Finally rotate and set orbitals in both ciwfn and v2rdm
        orbs_mat = psi4.core.doublet(start_orbs, totalU, False, False)
        ciwfn.set_orbitals("ROT", orbs_mat)
        v2rdm.set_orbitals("ROT", orbs_mat)

        # Figure out what the next step should be
        if (orb_grad_rms < mcscf_so_start_grad) and (abs(ediff) < abs(mcscf_so_start_e)) and\
                (mcscf_iter >= 2):

            if mcscf_target_conv_type == 'AH':
                approx_integrals_only = False
                ah_step = True
            elif mcscf_target_conv_type == 'OS':
                approx_integrals_only = False
                mcscf_current_step_type = 'OS, Prep'
                break
            else:
                continue
        #raise p4util.PsiException("")

    return current_energy
