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

import psi4
import re
import os
#import inputparser
import math
import warnings
#import driver
#from molutil import *
#from psi4.driver.procedures import proc_util
from psi4.driver.procrouting import proc_util

import psi4.driver.p4util as p4util
#from p4util.exceptions import *
#from procedures import *

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
