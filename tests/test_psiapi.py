import pytest


@pytest.mark.quick
def test_v2rdm1():
    """v2rdm_casscf/tests/v2rdm1"""
    #! cc-pvdz N2 (6,6) active space Test DQG

    print('        N2 / cc-pVDZ / DQG(6,6), scf_type = CD / 1e-12, rNN = 0.5 A')

    import psi4

    n2 = psi4.geometry("""
    0 1
    n
    n 1 r
    """)

    interloper = psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 90.0
    """)

    psi4.set_options({
      'basis': 'cc-pvdz',
      'scf_type': 'cd',
      'cholesky_tolerance': 1e-12,
      'd_convergence': 1e-10,
      'maxiter': 500,
      'restricted_docc': [ 2, 0, 0, 0, 0, 2, 0, 0 ],
      'active': [ 1, 0, 1, 1, 0, 1, 1, 1 ],
    })
    psi4.set_module_options('v2rdm_casscf', {
      'positivity': 'dqg',
      'r_convergence': 1e-5,
      'e_convergence': 1e-6,
      'maxiter': 20000,
      #'orbopt_frequency': 1000,
      #'mu_update_frequency': 1000,
    })

    psi4.activate(n2)

    n2.r     = 0.5 * 0.52917721067 / 0.52917720859
    refscf   = -103.04337420425350
    refv2rdm = -103.086205379481

    psi4.energy('v2rdm-casscf', molecule=n2)

    assert psi4.compare_values(refscf, psi4.variable("SCF TOTAL ENERGY"), 8, "SCF total energy")
    assert psi4.compare_values(refv2rdm, psi4.variable("CURRENT ENERGY"), 5, "v2RDM-CASSCF total energy")


def test_v2rdm2():
    #! cc-pvdz N2 (6,6) active space Test DQG

    print('        N2 / cc-pVDZ / DQG(6,6), scf_type = DF, rNN = 1.1 A')

    import psi4

    n2 = psi4.geometry("""
    0 1
    n
    n 1 r
    """)

    psi4.set_options({
      'basis':           'cc-pvdz',
      'scf_type':        'df',
      'd_convergence':   1e-10,
      'maxiter':         500,
      'restricted_docc': [ 2, 0, 0, 0, 0, 2, 0, 0 ],
      'active':          [ 1, 0, 1, 1, 0, 1, 1, 1 ],
    })
    psi4.set_module_options('v2rdm_casscf', {
      'positivity':      'dqg',
      'r_convergence':   1e-5,
      'e_convergence':   1e-6,
      'maxiter':         20000,
    })

    psi4.activate(n2)

    n2.r     = 1.1
    refscf   = -108.95348837831371
    refv2rdm = -109.094404909477

    psi4.energy('v2rdm-casscf')

    assert psi4.compare_values(refscf, psi4.variable("SCF TOTAL ENERGY"), 8, "SCF total energy")
    assert psi4.compare_values(refv2rdm, psi4.variable("CURRENT ENERGY"), 5, "v2RDM-CASSCF total energy")


@pytest.mark.quick
def test_v2rdm3():
    #! H3 / cc-pvdz / D+D3 vs full CI, scf_type = PK

    print('        H3 / cc-pvdz / D+D3 vs full CI, scf_type = PK')

    import psi4

    h3 = psi4.geometry("""
    0 2
    H
    H 1 1.0
    H 1 2.0 2 90.0
    """)

    psi4.set_options({
      'basis': 'sto-3g',
      'scf_type': 'pk',
      'reference': 'rohf',
      'guess': 'sad',
      'd_convergence': 1e-10,
      'maxiter': 500,
    })

    psi4.set_module_options('v2rdm_casscf', {
      'positivity': 'd',
      'constrain_d3': True,
      'r_convergence': 1e-5,
      'e_convergence': 1e-6,
      'maxiter': 20000,
      'optimize_orbitals': False,
      'semicanonicalize_orbitals': False,
    })

    psi4.activate(h3)
    v2rdm = psi4.energy('v2rdm-casscf')
    fci   = psi4.energy('fci')

    assert psi4.compare_values(v2rdm, fci, 5, "v2RDM vs full CI")


@pytest.mark.long
def test_v2rdm4():

    print('        1,4-phenylenedinitrene/(10,10)/cc-pVDZ')

    import psi4

    singlet = psi4.geometry("""
    0 1
    C      0.00000000    -1.43825841     0.00000000
    C      1.26321637    -0.67149862     0.00000000
    C      1.26321637     0.67149862     0.00000000
    C      0.00000000     1.43825841     0.00000000
    C     -1.26321637     0.67149862     0.00000000
    C     -1.26321637    -0.67149862     0.00000000
    H      2.18678709    -1.24095300     0.00000000
    H      2.18678709     1.24095300     0.00000000
    H     -2.18678709     1.24095300     0.00000000
    H     -2.18678709    -1.24095300     0.00000000
    N      0.00000000     2.71003942     0.00000000
    N      0.00000000    -2.71003942     0.00000000
    """)
    triplet = psi4.geometry("""
    0 3
    C      0.00000000    -1.43825841     0.00000000
    C      1.26321637    -0.67149862     0.00000000
    C      1.26321637     0.67149862     0.00000000
    C      0.00000000     1.43825841     0.00000000
    C     -1.26321637     0.67149862     0.00000000
    C     -1.26321637    -0.67149862     0.00000000
    H      2.18678709    -1.24095300     0.00000000
    H      2.18678709     1.24095300     0.00000000
    H     -2.18678709     1.24095300     0.00000000
    H     -2.18678709    -1.24095300     0.00000000
    N      0.00000000     2.71003942     0.00000000
    N      0.00000000    -2.71003942     0.00000000
    """)

    psi4.set_options({
      'basis': 'cc-pvdz',
      'scf_type': 'df',
      'd_convergence': 1e-8,
      'maxiter': 500,
      'reference': 'rohf',
      'restricted_docc': [   8,   3,   0,   0,   0,   0,   7,   4 ],
      'active':          [   0,   1,   1,   3,   1,   3,   0,   1 ],
    })
    psi4.set_module_options('v2rdm_casscf', {
      'positivity': 'dqg',
      'constrain_spin': True,

      'r_convergence': 1e-4,
      'e_convergence': 1e-4,
      'cg_convergence': 1e-8,

      'mu_update_frequency': 250,
      'orbopt_frequency': 500,

      'maxiter': 100000,
    })

    # singlet
    # from JCTC, with r_convergence = e_convergence = 1e-5
    #    * v2RDM total energy:           -338.502362198038

    psi4.activate(singlet)
    ref_singlet = -338.502342367694
    e_singlet = psi4.energy('v2rdm-casscf')

    assert psi4.compare_values(e_singlet, ref_singlet, 4, "singlet")

    # triplet
    # from JCTC, with r_convergence = e_convergence = 1e-5
    #    * v2RDM total energy:           -338.495146818081

    psi4.activate(triplet)
    ref_triplet = -338.495170913420
    e_triplet = psi4.energy('v2rdm-casscf')

    assert psi4.compare_values(e_triplet, ref_triplet, 4, "triplet")


def test_v2rdm5():
    #! cc-pvdz N2 (6,6) active space Test DQG

    print('        N2 / cc-pVDZ / DQG+T2(6,6), scf_type = PK, rNN = 1.1 A')

    import psi4

    n2 = psi4.geometry("""
    0 1
    n
    n 1 r
    """)

    psi4.set_options({
      'basis': 'cc-pvdz',
      'scf_type': 'pk',
      'd_convergence': 1e-10,
      'maxiter': 500,
      'restricted_docc': [ 2, 0, 0, 0, 0, 2, 0, 0 ],
      'active':          [ 1, 0, 1, 1, 0, 1, 1, 1 ],
    })
    psi4.set_module_options('v2rdm_casscf', {
      'positivity': 'dqgt2',
      'r_convergence':  1e-4,
      'e_convergence':  5e-4,
      'maxiter': 20000,
    })

    psi4.activate(n2)

    n2.r     = 1.1
    refscf   = -108.95379624015767
    refv2rdm = -109.091487394061

    psi4.energy('v2rdm-casscf')

    assert psi4.compare_values(refscf, psi4.variable("SCF TOTAL ENERGY"), 8, "SCF total energy")
    assert psi4.compare_values(refv2rdm, psi4.variable("CURRENT ENERGY"), 4, "v2RDM-CASSCF total energy")


def test_v2rdm6():
    #! cc-pvdz N2 (6,6) active space Test DQG

    print('        N2 / cc-pVDZ / DQG(6,6), geometry optimization')

    import psi4

    n2 = psi4.geometry("""
    0 1
    n
    n 1 1.1
    """)

    psi4.set_options({
      'basis': 'cc-pvdz',
      'scf_type': 'pk',
      'd_convergence': 1e-10,
      'maxiter': 500,
      'restricted_docc': [ 2, 0, 0, 0, 0, 2, 0, 0 ],
      'active':          [ 1, 0, 1, 1, 0, 1, 1, 1 ],
    })
    psi4.set_module_options('v2rdm_casscf', {
      'positivity': 'dqg',
      #'r_convergence': 1e-7,
      'r_convergence': 1e-6,
      'e_convergence': 1e-5,
      'orbopt_gradient_convergence': 1e-8,
      'maxiter': 20000,
    })

    psi4.activate(n2)

    psi4.optimize('v2rdm-casscf')

    refnuc   =   23.1968666562054260
    refscf   = -108.95016246035139
    refv2rdm = -109.095505119442

    assert psi4.compare_values(refnuc, n2.nuclear_repulsion_energy(),  4, "Nuclear repulsion energy")
    assert psi4.compare_values(refscf, psi4.variable("SCF TOTAL ENERGY"), 5, "SCF total energy")
    assert psi4.compare_values(refv2rdm, psi4.variable("CURRENT ENERGY"), 4, "v2RDM-CASSCF total energy")


def test_v2rdm7():
    # STO-3g benzene (6,6) guess orbital rotation test DQG

    print('        benzene (6,6), scf_type = PK')

    import psi4

    # orbital rotation needed before MCSCF calculation
    benzene_c1 = psi4.geometry("""
    0 1
    symmetry c1
     C                  0.00000000    1.38980400    0.00000000
     C                  1.20360500    0.69490200    0.00000000
     C                  1.20360500   -0.69490200    0.00000000
     C                  0.00000000   -1.38980400    0.00000000
     C                 -1.20360500   -0.69490200    0.00000000
     C                 -1.20360500    0.69490200    0.00000000
     H                  0.00000000    2.47523400    0.00000000
     H                  2.14361600    1.23761700    0.00000000
     H                  2.14361600   -1.23761700    0.00000000
     H                  0.00000000   -2.47523400    0.00000000
     H                 -2.14361600   -1.23761700    0.00000000
     H                 -2.14361600    1.23761700    0.00000000
    """)

    psi4.set_options({
      'basis': 'sto-3g',
      'scf_type': 'pk',
      'd_convergence': 1e-10,
      'maxiter': 500,
      'restricted_docc': [ 18 ],
      'active': [ 6 ],
    })

    psi4.energy ('hf')

    psi4.set_module_options('v2rdm_casscf', {
    # Switch the 17th (index 16) and the 19th (index 18) orbitals of the 1st irrep (index 0)
    # If more than one set of orbitals need to be rotated, use the following syntex
    # mcscf_rotate [[irrep_1, orb1_1, orb2_1, theta_1], [irrep_2, orb1_2, orb2_2, theta2],...]
    # Setting theta to 90 would switch the orbitals, setting it to 0 does nothing.
      'mcscf_rotate': [[ 0, 16, 18, 90 ]],
      'positivity': 'dq',
      'r_convergence': 1e-5,
      'e_convergence': 1e-6,
      'maxiter': 20000,
      'guess_orbitals_write': False,
      'molden_write': False,
    })

    psi4.activate(benzene_c1)

    E_c1 = psi4.energy('v2rdm-casscf')

    psi4.core.clean()

    # reference calculation, no need to rotate orbitals for MCSCF
    benzene_d2h = psi4.geometry("""
    0 1
    symmetry d2h
     C                  0.00000000    1.38980400    0.00000000
     C                  1.20360500    0.69490200    0.00000000
     C                  1.20360500   -0.69490200    0.00000000
     C                  0.00000000   -1.38980400    0.00000000
     C                 -1.20360500   -0.69490200    0.00000000
     C                 -1.20360500    0.69490200    0.00000000
     H                  0.00000000    2.47523400    0.00000000
     H                  2.14361600    1.23761700    0.00000000
     H                  2.14361600   -1.23761700    0.00000000
     H                  0.00000000   -2.47523400    0.00000000
     H                 -2.14361600   -1.23761700    0.00000000
     H                 -2.14361600    1.23761700    0.00000000
    """)

    psi4.set_options({
      'basis': 'sto-3g',
      'scf_type': 'pk',
      'd_convergence': 1e-10,
      'maxiter': 500,
      'restricted_docc': [ 6, 3, 0, 0, 0, 0, 5, 4 ],
      'active':          [ 0, 0, 1, 2, 1, 2, 0, 0 ],
    })

    psi4.energy('hf')

    psi4.set_module_options('v2rdm_casscf', {
    # Note that when mcscf_rotate is set for the previous molecule, the calculations afterwards
    #   also use this input, unless you overwrite it. A safe way to unset it is to overwrite it
    #   with mcscf_rotate [[0, 0, 0, 0]], which should work as long as your molecule has at
    #   least 1 orbital in the 1st irrep. This problem can also be simply avoided by calculating
    #   the molecules that do not require orbital rotations first.
      'mcscf_rotate': [[ 0, 0, 0, 0 ]],
      'positivity': 'dq',
      'r_convergence': 1e-5,
      'e_convergence': 1e-6,
      'maxiter': 20000,
      'guess_orbitals_write': False,
      'molden_write': False,
    })

    psi4.activate(benzene_d2h)

    E_d2h = psi4.energy('v2rdm-casscf')

    assert psi4.compare_values(E_c1, E_d2h, 4, "v2RDM-CASSCF total energy")


def test_v2rdm8():
    # H2 / cc-pVDZ / D(2,2), scf_type = DF, rHH = 1.0 A

    print('        H2 / cc-pVDZ / D(2,2), scf_type = DF, rHH = 1.0 A')

    import v2rdm_casscf
    import psi4

    h2 = psi4.geometry("""
    0 1
    h
    h 1 r
    symmetry c1
    """)

    psi4.set_options({
      "basis": "cc-pvdz",
      "mcscf_type": "df",
      "scf_type": "df",
      "d_convergence": 1e-10,
      "maxiter": 500,
      "restricted_docc": [  0 ],
      "restricted_uocc": [  8 ],
      "active":          [  2 ],
    })
    psi4.set_module_options('v2rdm_casscf', {
      "positivity": "d",
      "r_convergence": 1e-6,
      "e_convergence": 1e-8,
      "maxiter": 20000,
    })

    psi4.activate(h2)

    h2.r = 1.0

    psi4.set_options({"df_ints_io": "save"})
    en,wfn = psi4.energy('scf',return_wfn = True)

    en1 = v2rdm_casscf.v2rdm_scf_solver(wfn)
    en2 = psi4.energy('casscf')

    assert psi4.compare_values(en1, en2, 6, "v2RDM-CASSCF total energy")
