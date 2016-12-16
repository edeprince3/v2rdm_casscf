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
from psi4 import core
import psi4.driver.p4util as p4util
from psi4.driver.procedures import proc_util


def run_v2rdm_casscf(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    v2rdm_casscf can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('v2rdm_casscf')

    """

    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # Your plugin's psi4 run sequence goes here
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # if restarting from a checkpoint file, this file
    # needs to be in scratch with the correct name
    filename = core.get_option("V2RDM_CASSCF","RESTART_FROM_CHECKPOINT_FILE")

    # todo PSIF_V2RDM_CHECKPOINT should be definied in psifiles.h
    if ( filename != "" ):
        molname = ref_wfn.molecule().name()
        p4util.copy_file_to_scratch(filename,'psi',molname,269,False)

    # Ensure IWL files have been written when not using DF/CD
    scf_type = core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    returnvalue = core.plugin('v2rdm_casscf.so', ref_wfn)

    #core.set_variable('CURRENT ENERGY', returnvalue)

    #return core.get_variable('CURRENT ENERGY')
    return returnvalue


# Integration with driver routines
psi4.driver.procedures['energy']['v2rdm-casscf'] = run_v2rdm_casscf

def exampleFN():
    # Your Python code goes here
    pass
