import os
import sys

import pytest


base_dir = os.path.dirname(os.path.abspath(__file__))
plugin_sys_path = os.path.dirname(os.path.dirname(base_dir))


def exe_py(workspace, tdir, py):
    script = base_dir + '/' + tdir + '/' + py + '.py'
    workspace.run('python ' + script)


def exe_psi(workspace, tdir, p4):
    import psi4
    script = base_dir + '/' + tdir + '/' + p4
    workspace.run('PYTHONPATH={}:$PYTHONPATH '.format(plugin_sys_path) + psi4.executable + ' ' + script)


def test_v2rdm1(workspace):
    exe_psi(workspace, 'v2rdm1', 'input.dat')


@pytest.mark.quick
def test_v2rdm2(workspace):
    exe_psi(workspace, 'v2rdm2', 'input.dat')


def test_v2rdm3(workspace):
    exe_psi(workspace, 'v2rdm3', 'input.dat')


def test_v2rdm4(workspace):
    exe_psi(workspace, 'v2rdm4', 'input.dat')


@pytest.mark.long
def test_v2rdm5(workspace):
    exe_psi(workspace, 'v2rdm5', 'input.dat')


def test_v2rdm6(workspace):
    exe_psi(workspace, 'v2rdm6', 'input.dat')


def test_v2rdm7(workspace):
    exe_psi(workspace, 'v2rdm7', 'input.dat')


def test_v2rdm8(workspace):
    exe_psi(workspace, 'v2rdm8', 'input.dat')
