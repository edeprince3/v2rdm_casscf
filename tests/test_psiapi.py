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

    n2.r     = 0.5
    refscf   = -103.04337420425350
    refv2rdm = -103.086205379481

    psi4.energy('v2rdm-casscf', molecule=n2)

    assert psi4.compare_values(refscf, psi4.get_variable("SCF TOTAL ENERGY"), 8, "SCF total energy")
    assert psi4.compare_values(refv2rdm, psi4.get_variable("CURRENT ENERGY"), 5, "v2RDM-CASSCF total energy")


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

    assert psi4.compare_values(refscf, psi4.get_variable("SCF TOTAL ENERGY"), 8, "SCF total energy")
    assert psi4.compare_values(refv2rdm, psi4.get_variable("CURRENT ENERGY"), 5, "v2RDM-CASSCF total energy")


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

    assert psi4.compare_values(refscf, psi4.get_variable("SCF TOTAL ENERGY"), 8, "SCF total energy")
    assert psi4.compare_values(refv2rdm, psi4.get_variable("CURRENT ENERGY"), 4, "v2RDM-CASSCF total energy")


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
    assert psi4.compare_values(refscf, psi4.get_variable("SCF TOTAL ENERGY"), 5, "SCF total energy")
    assert psi4.compare_values(refv2rdm, psi4.get_variable("CURRENT ENERGY"), 4, "v2RDM-CASSCF total energy")
