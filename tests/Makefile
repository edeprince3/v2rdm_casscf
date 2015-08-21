# v2rdm-casscf test suite
default: test

test:
	@echo "\033[0;36m"
	@echo ""
	@echo "    $(tput setaf 6)Running tests..."
	@echo ""
	@echo "\033[0;0m"
	@echo "    v2rdm1: N2 / cc-pVDZ / DQG(6,6), scf_type = CD / 1e-12, rNN = 0.5 A"
	@echo ""
	@cd v2rdm1; psi4 -d
	@echo ""
	@echo "    v2rdm2: N2 / cc-pVDZ / DQG(6,6), scf_type = PK, rNN = 1.1 A"
	@echo ""
	@cd v2rdm2; psi4 -d
	@echo ""
	@echo "    v2rdm3: N2 / cc-pVDZ / DQG(6,6), scf_type = DF, rNN = 3.0 A"
	@echo ""
	@cd v2rdm3; psi4 -d
	@echo ""
