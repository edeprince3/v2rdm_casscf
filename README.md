# v2rdm-casscf
A variational 2-RDM-driven CASSCF plugin to Psi4

OVERVIEW
---

This plugin to Psi4[1] performs variational two-electron reduced-density-matrix (2-RDM)-driven complete active space self consistent field (CASSCF) computations.  In principle, because  variational 2-RDM (v-2RDM) methods scales only polynomially with system size, v-2RDM-driven CASSCF computations can be performed using active spaces that are larger than can be used within conventional configuration-interaction-driven CASSCF methods.  For more information regarding the performance of the method, see Refs. 2-3.

INSTALLATION
---

To run the psi4 plugin v2rdm_casscf:

* Download and install psi4public from github.com:
https://github.com/psi4/psi4public.  You can obtain the source using git:

    > git clone git@github.com:psi4/psi4public.git

    Install psi4 as described on http://www.psicode.org/.

*  Make a fresh makefile for your plugin

  > psi4 --new-plugin-makefile
  
  > vi Makefile

* search for BINOBJ; delete everything BELOW that line
* add what is below (with the appropriate fortran compiler and flags) to the Makefile:

```
F90       = gfortran-mp-4.8
F90SRC    = $(notdir $(wildcard *.F90))
F90BINOBJ = $(F90SRC:%.F90=%.o)
F90FLAGS  = -O2 -fPIC
LDFLAGS  += -L/opt/local/lib/gcc48/ -lgfortran

fortran:
    $(F90) $(F90FLAGS) jacobi_data.F90 -c
    $(F90) $(F90FLAGS) jacobi_maxind_mod.F90 jacobi_data.o -c
    $(F90) $(F90FLAGS) jacobi_mod.F90 jacobi_data.o jacobi_maxind_mod.o -c
    $(F90) $(F90FLAGS) jacobi_interface.F90 jacobi_mod.o jacobi_data.o jacobi_maxind_mod.o -c
    rm *.mod

%.o: %.F90
    $(F90) $(F90FLAGS) -c $<

%.o: %.cc
    $(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<

$(PSITARGET): $(BINOBJ) $(F90BINOBJ)
    $(CXX) $(LDFLAGS) -o $@ $^ $(CXXDEFS) $(PSIPLUGIN)

clean:
    rm -f $(F90BINOBJ) $(BINOBJ) $(PSITARGET) *.d *.pyc *.test output.dat psi.timer.dat
```

* Compile the plugin

> make

INPUT OPTIONS
---
* **POSITIVITY** (string):

    The positivity conditions enforced in the computation.  Allowed values include DQG, DQ, DG, D, DQGT1, DQGT2, and DQGT1T2.  The default value is DQG.

* **CONSTRAIN_D3** (bool):

    Do constrain the 3-RDM to be positive?  Default no.
