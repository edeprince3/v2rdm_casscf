# v2rdm_casscf
A variational 2-RDM-driven CASSCF plugin to Psi4

##OVERVIEW

This plugin to Psi4[1] performs variational two-electron reduced-density-matrix (2-RDM)-driven complete active space self consistent field (CASSCF) computations.  In principle, because  variational 2-RDM (v-2RDM) methods scale only polynomially with system size, v-2RDM-driven CASSCF computations can be performed using active spaces that are larger than can be used within conventional configuration-interaction-driven CASSCF methods.  For more information regarding the performance of the method, see Refs. 2-3.

##INSTALLATION

To run the Psi4 plugin v2rdm_casscf:

* Download Psi4 (1.1a2.dev200 or later) from github.com: https://github.com/psi4/psi4, and follow the installation instructions given here: http://psicode.org/psi4manual/master/build_planning.html . Make sure to keep the name of the plugin directory v2rdm_casscf.

*  Configure with CMake to generate a Makefile. Run `psi4 --plugin-compile` to get a CMake command. Modify it as needed with `-D` for compiler, libraries, and options.

* Note that, if you configured Psi4 with a fortran compiler, you shouldn't have to specify these things here. If the configure shows no errors, compile the plugin:

  > make

* If the plugin compiles without any errors, you can run a few tests:

  > cd tests

  > make

* The test directories (tests/v2rdm1, etc.) contain input files that can help you get started using v2rdm-casscf.

##INPUT OPTIONS

###N-representability conditions

* **POSITIVITY** (string):

    The positivity conditions enforced in the computation.  Allowed values
    include DQG, DQ, DG, D, DQGT1, DQGT2, and DQGT1T2.  The default value
    is DQG.

* **CONSTRAIN_D3** (bool):

    Enforce the additional condition that D3 be possitive and correctly
    contract to D2?  Default false.

* **CONSTRAIN_SPIN** (bool):

    Do constrain the expectation value of spin squared? Default true.

###Convergence

* **E_CONVERGENCE** (double):

    The convergence in the primal/dual energy gap.  Default 1e-4.

* **R_CONVERGENCE** (double):

    The convergence in the primal and dual errors. Default 1e-4.

* **MAXITER** (int):

    The maximum number of outer iterations.  Default 10000.

###Active space specification

* **FROZEN_DOCC** (array):

    An array containing the number of frozen doubly-occupied orbitals per
    irrep.  These orbitals are not included in the active space, nor are
    they optimized during the v2RDM-CASSCF procedure.  This option trumps
    Psi4's **NUM_FROZEN_DOCC** and **FREEZE_CORE** options.

* **FROZEN_UOCC** (array):

    An array containing the number of frozen unoccupied orbitals per
    irrep.  These orbitals are not included in the active space, nor are
    they optimized during the v2RDM-CASSCF procedure.  This option trumps
    Psi4's **NUM_FROZEN_UOCC** option.

* **RESTRICTED_DOCC** (array):

    An array giving the number of restricted doubly-occupied orbitals per 
    irrep.  These orbitals are not included in the active space, but they are
    optimized by the v2RDM-CASSCF procedure.

* **RESTRICTED_UOCC** (array):

    An array giving the number of restricted unoccupied orbitals per 
    irrep.  These orbitals are not included in the active space, but they are
    optimized by the v2RDM-CASSCF procedure.

* **ACTIVE** (array):

    An array giving the number of active orbitals (occupied plus
    unoccupied) per irrep.  This option provides a more intuitive way of
    specifying the active space than the
    **FROZEN_DOCC**/**RESTRICTED_DOCC**/**RESTRICTED_UOCC**/**FROZEN_UOCC**
    keywords.  The simplest specification of the active space would
    involve this keyword and at the **RESTRICTED_DOCC** keyword.  This
    option trumps the **RESTRICTED_UOCC** option, which will be determined
    from the **ACTIVE**, **RESTRICTED_DOCC**, **FROZEN_DOCC**, and
    **FROZEN_UOCC** arrays.

###Restarting jobs

* **WRITE_CHECKPOINT_FILE** (bool):

    Do save progress in a checkpoint file?  Default false.

* **CHECKPOINT_FREQUENCY** (bool):

    Frequency of checkpoint file generation.  The checkpoint file is 
    updated every **CHECKPOINT_FREQUENCY** iterations.  The default frequency
    will be **ORBOPT_FREQUENCY**.

* **RESTART_FROM_CHECKPOINT_FILE** (string):

    File containing previous primal/dual solutions and integrals.

###Integrals and SCF type

* **DF_BASIS_SCF** (string):

    Auxiliary basis set for SCF density fitting computations.  Defaults
    to a JKFIT basis.

* **SCF_TYPE** (string):

    What algorithm to use for the initial SCF computation.  Default DF.

* **CHOLESKY_TOLERANCE** (double):

    Tolerance for Cholesky decomposition of the ERI tensor.  Default 1e-4.

###Orbital optimization

* **ORBOPT_ONE_STEP** (int):

    Flag to optimize orbitals using a quasi one-step type approach. Default 1.

* **ORBOPT_GRADIENT_CONVERGENCE** (double):

    Convergence in the orbital gradient norm.  Default 1e-4.

* **ORBOPT_ENERGY_CONVERGENCE** (double):

    Convergence in the energy for orbital rotations. Default 1e-8.

* **ORBOPT_EXACT_DIAGONAL_HESSIAN** (int):

    Flag for using exact expresions for diagonal Hessian elements.  Default 0.

* **ORBOPT_FREQUENCY** (int):

    Frequency of orbital optimization.  Optimization occurs every 
    **ORBOPT_FREQUENCY** iterations.  Default 200.

* **ORBOPT_ACTIVE_ACTIVE_ROTATIONS** (bool):

    Do rotate active/active orbital pairs? Default false.

###Additional files

* **MOLDEN_WRITE** (bool):

    Do write a MOLDEN output file containing the natural orbitals?  If
    yes, the filename will end in .molden, and the prefix is determined by
    **WRITER_FILE_LABEL** (if set), or else by the name of the output file
    plus the name of the current molecule.  Default false.

* **WRITER_FILE_LABEL** (string):

    Base filename for text files written by PSI, such as the MOLDEN output
    file, the Hessian file, the internal coordinate file, etc. Use the
    add_str_i function to make this string case sensitive.

* **ORBOPT_WRITE** (bool):

    Do write a ORBOPT output file?  If so, the filename will end in .orbopt,
    and the prefix is determined by **WRITER_FILE_LABEL** (if set), or else by
    the name of the output file plus the name of the current molecule.


##KNOWN ISSUES

* For large jobs, when running with multiple threads, sometimes a thread will hang and the job will stall.
* For large jobs, add "ulimit -s unlimited" to .bashrc to avoid segfault when calling the fortran orbital optimization routines.

##REFERENCES

[1] J. M. Turney, A. C. Simmonett, R. M. Parrish, E. G. Hohenstein, F. A. Evangelista, J. T. Fermann, B. J.  Mintz, L. A. Burns, J. J. Wilke, M. L. Abrams, N. J. Russ, M. L. Leininger, C. L. Janssen, E. T. Seidl, W. D. Allen, H. F. Schaefer, R. A. King, E. F. Valeev, C. D. Sherrill, and T. D. Crawford, *WIREs: Comp. Molec. Sci.* **2**, 556 (2012). "Psi4: an open-source ab initio electronic structure program"

[2] J. Fosso-Tande, D. R. Nascimento, and A. E. DePrince III, *Mol. Phys.* **114**, 423-430 (2015). "Accuracy of two-particle N-representability conditions for describing different spin states and the singlet-triplet gap in the linear acene series." http://dx.doi.org/10.1080/00268976.2015.1078008

[3] J. Fosso-Tande, T.-S. Nguyen, G. Gidofalvi, and A. E. DePrince III, *J. Chem. Theory Comput.*, accepted (2016).  "Large-scale v2RDM-driven CASSCF methods."  http://dx.doi.org/10.1021/acs.jctc.6b00190
