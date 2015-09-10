# v2rdm_casscf
A variational 2-RDM-driven CASSCF plugin to Psi4

OVERVIEW
---

This plugin to Psi4[1] performs variational two-electron reduced-density-matrix (2-RDM)-driven complete active space self consistent field (CASSCF) computations.  In principle, because  variational 2-RDM (v-2RDM) methods scale only polynomially with system size, v-2RDM-driven CASSCF computations can be performed using active spaces that are larger than can be used within conventional configuration-interaction-driven CASSCF methods.  For more information regarding the performance of the method, see Refs. 2-3.

INSTALLATION
---

To run the psi4 plugin v2rdm_casscf:

* Download psi4public from github.com: https://github.com/psi4/psi4public, and follow the installation instructions given here: http://www.psicode.org/.

*  Modify the setup script by specifying your fortran compiler, options, and libraries

  > cd v2rdm_casscf
  
  > vi setup

*  Modify F90, F90_FLAGS, and F90_LIB, as needed.

*  Compile the plugin:

  > ./setup

INPUT OPTIONS
---
* **POSITIVITY** (string):

    The positivity conditions enforced in the computation.  Allowed values include DQG, DQ, DG, D, DQGT1, DQGT2, and DQGT1T2.  The default value is DQG.

* **E_CONVERGENCE** (double):

    The convergence in the primal/dual energy gap.  Default 1e-4.

* **R_CONVERGENCE** (double):

    The convergence in the primal and dual errors. Default 1e-3.

* **MAXITER** (int):

    The maximum number of outer iterations.  Default 10000.

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

* **MOLDEN_WRITE** (bool):

    Do write a MOLDEN output file containing the natural orbitals?  If yes, the filename will end in .molden, and the prefix is determined by **WRITER_FILE_LABEL** (if set), or else by the name of the output file plus the name of the current molecule.  Default false.

* **WRITER_FILE_LABEL** (string):

    Base filename for text files written by PSI, such as the MOLDEN output file, the Hessian file, the internal coordinate file, etc. Use the add_str_i function to make this string case sensitive.

KNOWN ISSUES
---

* For large jobs, when running with multiple threads, sometimes a thread will hang and the job will stall.
* For large jobs, add "ulimit -s unlimited" to .bashrc to avoid segfault when calling the fortran orbital optimization routines.

REFERENCES
---

[1] J. M. Turney, A. C. Simmonett, R. M. Parrish, E. G. Hohenstein, F. A. Evangelista, J. T. Fermann, B. J.  Mintz, L. A. Burns, J. J. Wilke, M. L. Abrams, N. J. Russ, M. L. Leininger, C. L. Janssen, E. T. Seidl, W. D. Allen, H. F. Schaefer, R. A. King, E. F. Valeev, C. D. Sherrill, and T. D. Crawford, *WIREs: Comp. Molec. Sci.* **2**, 556 (2012). "Psi4: an open-source ab initio electronic structure program"

[2] J. Fosso-Tande, D. R. Nascimento, and A. E. DePrince III, *Mol. Phys.* accepted (2015). "Accuracy of two-particle N-representability conditions for describing different spin states and the singlet-triplet gap in the linear acene series"

[3] J. Fosso-Tande, G. Gidofalvi, and A. E. DePrince III, in preparation.
