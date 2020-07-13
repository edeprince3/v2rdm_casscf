/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 *
 *@END LICENSE
 *
 */

#ifndef V2RDM_SOLVER_H
#define V2RDM_SOLVER_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <psi4/libiwl/iwl.h>
#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libqt/qt.h>

#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>

// greg
#include"fortran.h"

// TODO: move to psifiles.h
#define PSIF_DCC_QMO          268
#define PSIF_V2RDM_CHECKPOINT 269
#define PSIF_V2RDM_D2AA       270
#define PSIF_V2RDM_D2AB       271
#define PSIF_V2RDM_D2BB       272
#define PSIF_V2RDM_D3AAA      273
#define PSIF_V2RDM_D3AAB      274
#define PSIF_V2RDM_D3BBA      275
#define PSIF_V2RDM_D3BBB      276
#define PSIF_V2RDM_D1A        277
#define PSIF_V2RDM_D1B        278
#define PSIF_V2RDM_D2_SPIN_FREE 279

namespace psi{ namespace v2rdm_casscf{

struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double value;
};

struct opdm {
    int i;
    int j;
    double value;
};

class v2RDMSolver: public Wavefunction{
  public:
    v2RDMSolver(SharedWavefunction reference_wavefunction,Options & options);
    ~v2RDMSolver();
    void common_init();

    double compute_energy();

    /// return spin-free one-particle density matrix. full space. sparse
    std::vector<opdm> get_opdm_sparse(std::string type);

    /// return spin-free two-particle density matrix. full space. sparse.
    std::vector<tpdm> get_tpdm_sparse(std::string type);

    /// return spin-free one-particle density matrix as shared matrix for python API
    std::shared_ptr<Matrix> get_opdm();

    /// return spin-free two-particle density matrix as shared matrix for python API
    std::shared_ptr<Matrix> get_tpdm();

    /// return subset of orbitals for python API
    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the definitions used by DETCI for now.
     * @param  orbital_name fzc, drc, docc, act, ras1, ras2, ras3, ras4, pop, vir, fzv, drv, or all
     * @return C            Returns the appropriate orbitals in the SO basis.
     */
    std::shared_ptr<Matrix> get_orbitals(const std::string &orbital_name);

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name FZC, DRC, DOCC, ACT, RAS1, RAS2, RAS3, RAS4, POP, VIR, FZV, DRV, or ALL
     * @param  orbitals     SharedMatrix to set
     */
    void set_orbitals(const std::string &orbital_name, SharedMatrix orbitals);

    /// Find out which orbitals belong where
    void orbital_locations(const std::string &orbital_name, int *start, int *end);


    // public methods
    void cg_Ax(long int n,SharedVector A, SharedVector u);

  protected:

    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

    double nalpha_;
    double nbeta_;

    /// constrain Q2 to be positive semidefinite?
    bool constrain_q2_;

    /// constrain G2 to be positive semidefinite?
    bool constrain_g2_;

    /// spin adapt g2 constraint?
    bool spin_adapt_g2_;

    /// spin adapt q2 constraint?
    bool spin_adapt_q2_;

    /// constraint T1 = D3 + Q3 to be positive semidefinite?
    bool constrain_t1_;

    /// constraint T2 = E3 + F3 to be positive semidefinite?
    bool constrain_t2_;

    /// keep d3 positive semidefinite and constrain D3->D2 mapping?
    bool constrain_d3_;

    /// keep d4 positive semidefinite and constrain D4->D3 mapping?
    bool constrain_d4_;

    /// constrain spin?
    bool constrain_spin_;

    /// constrain sz?
    bool constrain_sz_;

    /// symmetry product table:
    int * table;

    /// returns symmetry product for two orbitals
    int SymmetryPair(int i, int j);

    /// returns symmetry product for four orbitals
    int TotalSym(int i, int j,int k, int l);

    int * symmetry;
    int * symmetry_full;
    int * symmetry_really_full;
    int * energy_to_pitzer_order;
    int * energy_to_pitzer_order_really_full;
    int * symmetry_energy_order;
    int * pitzer_offset;
    int * pitzer_offset_full;

    /// active-space geminals for each symmetry:
    std::vector < std::vector < std::pair<int,int> > > gems;

    /// total number of active molecular orbitals
    int amo_;

    /// total number of frozen core orbitals
    int nfrzc_;

    /// total number of frozen virtual orbitals
    int nfrzv_;

    /// total number of restricted doubly occupied orbitals
    int nrstc_;

    /// total number of restricted unoccupied orbitals
    int nrstv_;

    /// active molecular orbitals per irrep
    int * amopi_;

    /// restricted core orbitals per irrep.  these will be optimized optimized
    int * rstcpi_;

    /// restricted virtual orbitals per irrep.  these will be optimized optimized
    int * rstvpi_;

    /// total number of constraints (dimension of dual solution vector)
    long int nconstraints_;

    /// total number of variables (dimension of primal solution vector)
    long int dimx_;

    /// number of auxilliary basis functions
    int nQ_;

    /// read three-index integrals and transform them to MO basis
    void ThreeIndexIntegrals();

    /// three-index integral buffer
    double * Qmo_;

    /// grab one-electron integrals (T+V) in MO basis
    SharedMatrix GetOEI();

    /// DIIS stuff
    void DIIS(double*c,long int nvec,int replace_diis_iter);
    void DIIS_WriteOldVector(long int iter,int diis_iter,int replace_diis_iter);
    void DIIS_WriteErrorVector(int diis_iter,int replace_diis_iter,int iter);
    void DIIS_Extrapolate(int diis_iter,int&replace_diis_iter);
    long int maxdiis_;
    double * diisvec_;
    double * junk1;
    double * junk2;
    long int diis_oiter_;
    long int dimdiis_;

    /// offsets
    int * d1aoff;
    int * d1boff;
    int * q1aoff;
    int * q1boff;
    int * d2aboff;
    int * d2aaoff;
    int * d2bboff;
    int * d200off;
    int * q2aboff;
    int * q2aaoff;
    int * q2bboff;
    int * g2aboff;
    int * g2baoff;
    int * g2aaoff;
    int * g2soff;
    int * g2toff;
    int * g2toff_p1;
    int * g2toff_m1;
    int * q2soff;
    int * q2toff;
    int * q2toff_p1;
    int * q2toff_m1;
    int * t1aaaoff;
    int * t1bbboff;
    int * t1aaboff;
    int * t1bbaoff;
    int * t2aaaoff;
    int * t2bbboff;
    int * t2aaboff;
    int * t2bbaoff;
    int * d3aaaoff;
    int * d3bbboff;
    int * d3aaboff;
    int * d3bbaoff;
    int * d4aaaaoff;
    int * d4aaaboff;
    int * d4aabboff;
    int * d4bbbaoff;
    int * d4bbbboff;

    /// convergence in primal energy
    double e_convergence_;

    /// convergence in primal and dual error
    double r_convergence_;

    /// convergence in conjugate gradient solver
    double cg_convergence_;

    /// maximum number of boundary-point SDP (outer) iterations
    int maxiter_;

    /// maximum number of conjugate gradient (inner) iterations
    int cg_maxiter_;

    /// standard vector of dimensions of each block of primal solution vector
    std::vector<int> dimensions_;

    int offset;

    // mapping arrays with abelian symmetry
    void BuildBasis();
    int * full_basis;

    /// mapping arrays with symmetry:
    int * gems_ab;
    int * gems_aa;
    int * gems_00;
    int * gems_full;
    int * gems_plus_core;
    int *** bas_ab_sym;
    int *** bas_aa_sym;
    int *** bas_00_sym;
    int *** bas_full_sym;
    int *** bas_really_full_sym;
    int *** ibas_ab_sym;
    int *** ibas_aa_sym;
    int *** ibas_00_sym;
    int *** ibas_full_sym;
    int *** ibas_really_full_sym;

    /// triplets for each irrep:
    std::vector < std::vector < std::tuple<int,int,int> > > triplets;
    int * trip_aaa;
    int * trip_aab;
    int * trip_aba;
    int  ***  bas_aaa_sym;
    int  ***  bas_aab_sym;
    int  ***  bas_aba_sym;
    int **** ibas_aaa_sym;
    int **** ibas_aab_sym;
    int **** ibas_aba_sym;

    /// quartets for each irrep:
    std::vector < std::vector < std::tuple<int,int,int,int> > > quartets;
    int * quartet_aaaa;
    int * quartet_aaab;
    int * quartet_aabb;
    int  ***  bas_aaaa_sym;
    int  ***  bas_aaab_sym;
    int  ***  bas_aabb_sym;
    int ***** ibas_aaaa_sym;
    int ***** ibas_aaab_sym;
    int ***** ibas_aabb_sym;

    void PrintHeader();

    /// grab one- and two-electron integrals
    void GetIntegrals();

    /// read two-electron integrals from disk
    void GetTEIFromDisk();

    // read teis from disk:
    void ReadAllIntegrals(iwlbuf *Buf);

    /// grab a specific two-electron integral
    double TEI(int i, int j, int k, int l, int h);

    void BuildConstraints();

    void Guess();
    void T1_constraints_guess(SharedVector u);
    void T2_constraints_guess(SharedVector u);
    void Q2_constraints_guess(SharedVector u);
    void Q2_constraints_guess_spin_adapted(SharedVector u);
    void G2_constraints_guess(SharedVector u);
    void G2_constraints_guess_spin_adapted(SharedVector u);

    void bpsdp_Au(SharedVector A, SharedVector u);
    void bpsdp_Au_slow(SharedVector A, SharedVector u);
    void Spin_constraints_Au(SharedVector A,SharedVector u);
    void D2_constraints_Au(SharedVector A,SharedVector u);
    void Q2_constraints_Au(SharedVector A,SharedVector u);
    void Q2_constraints_Au_spin_adapted(SharedVector A,SharedVector u);
    void G2_constraints_Au(SharedVector A,SharedVector u);
    void G2_constraints_Au_spin_adapted(SharedVector A,SharedVector u);
    void T1_constraints_Au(SharedVector A,SharedVector u);
    void T2_constraints_Au(SharedVector A,SharedVector u);
    void T2_constraints_Au_slow(SharedVector A,SharedVector u);
    void T2_tilde_constraints_Au(SharedVector A,SharedVector u);
    void D3_constraints_Au(SharedVector A,SharedVector u);
    void D4_constraints_Au(SharedVector A,SharedVector u);

    void bpsdp_ATu(SharedVector A, SharedVector u);
    void bpsdp_ATu_slow(SharedVector A, SharedVector u);
    void Spin_constraints_ATu(SharedVector A,SharedVector u);
    void D2_constraints_ATu(SharedVector A,SharedVector u);
    void Q2_constraints_ATu(SharedVector A,SharedVector u);
    void Q2_constraints_ATu_spin_adapted(SharedVector A,SharedVector u);
    void G2_constraints_ATu(SharedVector A,SharedVector u);
    void G2_constraints_ATu_spin_adapted(SharedVector A,SharedVector u);
    void T1_constraints_ATu(SharedVector A,SharedVector u);
    void T2_constraints_ATu(SharedVector A,SharedVector u);
    void T2_constraints_ATu_slow(SharedVector A,SharedVector u);
    void T2_tilde_constraints_ATu(SharedVector A,SharedVector u);
    void D3_constraints_ATu(SharedVector A,SharedVector u);
    void D4_constraints_ATu(SharedVector A,SharedVector u);

    /// SCF energy
    double escf_;

    /// nuclear repulsion energy
    double enuc_;

    double tau, mu, ed, ep;

    //vectors
    SharedVector Ax;     // vector to hold A . x
    SharedVector ATy;    // vector to hold A^T . y
    SharedVector c;      // 1ei and 2ei of bpsdp
    SharedVector y;      // dual solution
    SharedVector b;      // constraint vector
    SharedVector x;      // primal solution
    SharedVector z;      // second dual solution
    SharedVector rx;       // square root of x (for diis)
    SharedVector rz;       // square root of z (for diis)
    SharedVector rx_error; // error vector for x (for diis)
    SharedVector rz_error; // error vector for z (for diis)

    void Update_xz();
    void Update_xz_nonsymmetric();

    /// extended koopman's theorem
    void ExtendedKoopmans();
    void EKTEigensolver(std::shared_ptr<Matrix> V, std::shared_ptr<Matrix> D, std::shared_ptr<Vector> epsilon, bool use_dggev,std::string spin);


    /// compute natural orbitals and transform OPDM and TPDM to natural orbital basis
    void ComputeNaturalOrbitals();

    /// compute and print natural orbital occupation numbers
    void PrintNaturalOrbitalOccupations();

    /// push OPDM onto the wavefunction
    void FinalizeOPDM();

    // multiplicity
    int multiplicity_;

    /// full space of integrals for MO gradient / Hessian, blocked by symmetry
    double * tei_full_sym_;
    double * oei_full_sym_;
    // gidofalvi -- modified the type of tei_full_dim_ so that it is correct for large bases
    long int tei_full_dim_;
    int oei_full_dim_;
    std::shared_ptr<Matrix> T_;
    std::shared_ptr<Matrix> V_;

    /// full space D2, blocked by symmetry
    double * d2_plus_core_sym_;
    int d2_plus_core_dim_;

    /// active space spatial D2, blocked by symmetry
    double * d2_act_spatial_sym_;
    int d2_act_spatial_dim_;

    /// active space spatial D1, blocked by symmetry
    double * d1_act_spatial_sym_;
    int d1_act_spatial_dim_;

    /// unpack active-space density into full-space density
    void UnpackDensityPlusCore();

    /// pack active-space spin-blocked density into spatial density
    void PackSpatialDensity();

    /// repack rotated full-space integrals into active-space integrals
    void RepackIntegrals();
    void RepackIntegralsDF();

    /// compute frozen core energy and adjust oeis
    void FrozenCoreEnergy();

    /// function to rotate orbitals
    void RotateOrbitals();

    /// function to exponentiate step vector
    void exponentiate_step(double * X);

    double * orbopt_transformation_matrix_;
    double * orbopt_data_;
    char * orbopt_outfile_;
    bool orbopt_converged_;

    /// is this a v2RDM-DOCI computation?
    bool is_doci_;

    double doci_alpha_;
    double doci_ref_;

    /// are we using 3-index integrals?
    bool is_df_;

    /// write primal, dual, and orbitals to a checkpoint file
    void WriteCheckpointFile();

    /// read primal, dual, and orbitals from a checkpoint file
    void ReadFromCheckpointFile();

    /// read orbitals from a checkpoint file
    void ReadOrbitalsFromCheckpointFile();

    /// wall time for microiterations
    double iiter_time_;

    /// wall time for macroiterations
    double oiter_time_;

    /// wall time for orbital optimization
    double orbopt_time_;

    /// total number of microiterations
    long int iiter_total_;

    /// total number of macroiterations
    long int oiter_total_;

    /// total number of orbital optimization
    long int orbopt_iter_total_;

    /// write full 2RDM to disk
    void WriteTPDM();

    /// write full spin-free 2RDM to disk
    void WriteTPDMSpinFree();

    /// write full 1RDM to disk
    void WriteOPDM();

    /// write full 2RDM to disk in IWL format
    void WriteTPDM_IWL();

    /// write active-active-active-active 2RDM to disk
    void WriteActiveTPDM();

    /// read 2RDM from disk
    void ReadTPDM();

    /// write active 3RDM to disk
    void WriteActive3PDM();

    /// write molden file
    void WriteMoldenFile();

    /// read 3RDM from disk
    void Read3PDM();

    /// check spin structure of 1- and 2-RDM
    void CheckSpinStructure();

    /// orbital lagrangian
    double * X_;

    void OrbitalLagrangian();
    void DualD1Q1();

    /// memory available beyond what is allocated for v2RDM-CASSCF
    long int available_memory_;

    /// update primal solution after semicanonicalization
    void UpdatePrimal();

    /// transform a four-index quantity from one basis to another
    void TransformFourIndex(double * inout, double * tmp, SharedMatrix trans);

    /// update ao/mo transformation matrix after orbital optimization
    void UpdateTransformationMatrix();

    /// mo-mo transformation matrix
    SharedMatrix newMO_;

    /// FCIDUMP: dump integrals and RDMs to disk
    void FCIDUMP();

    /// break down energy into components
    void EnergyByComponent(double doci_ref, double doci_alpha, double &kinetic, double &potential, double &two_electron_energy);

};

}}
#endif

