# v2rdm-casscf test suite

SHELL := /bin/bash

# add new tests here
all-tests := $(addsuffix .test, v2rdm1 v2rdm2 v2rdm3)

quick-tests := $(addsuffix .test, v2rdm1)

.PHONY : test all %.test

test: $(all-tests)

quick: $(quick-tests)

%.test : 
	@echo ""
	@echo "    $(basename $@):"
	@echo ""
	@cd $(basename $@); psi4 -d
	@echo ""

