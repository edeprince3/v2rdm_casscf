/*                                                                      
 *@BEGIN LICENSE

  Copyright (c) 2014, The Florida State University. All rights reserved.

 *@END LICENSE
 */

#ifndef BPSDP_H
#define BPSDP_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>

#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <../bin/fnocc/blas.h>
#include <libqt/qt.h>

#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#include<libmints/vector.h>

//analytical gradient
#include "GRAD_tensors.h"
#include "GRAD_arrays.h"

// greg
#include"fortran.h"

// TODO: move to psifiles.h
#define PSIF_DCC_QMO          268
#define PSIF_V2RDM_CHECKPOINT 269

namespace boost {
  template<class T> class shared_ptr;
}
namespace psi{ namespace v2rdm_grad{

class v2RDMSolver: public Wavefunction{
  public: 
    v2RDMSolver(boost::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~v2RDMSolver();
    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; } 

    // public methods
    void cg_Ax(long int n,SharedVector A, SharedVector u);

  protected:

    boost::shared_ptr<Matrix> saveOEI_;

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

    /// constrain spin?
    bool constrain_spin_;

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
    
    /// number of doubly occupied MO
    int ndocc;

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
    boost::shared_ptr<Matrix> GetOEI();

    /// offsets
    int * d1aoff;  
    int * d1boff;  
    int * q1aoff;  
    int * q1boff;  
    int * d2aboff; 
    int * d2aaoff; 
    int * d2bboff; 
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
    int * t2abaoff;
    int * t2baboff;
    int * d3aaaoff;
    int * d3bbboff;
    int * d3aaboff;
    int * d3bbaoff;

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

    /// triplets for each irrep:
    std::vector < std::vector < boost::tuple<int,int,int> > > triplets;
    int * trip_aaa;
    int * trip_aab;
    int * trip_aba;
    int  ***  bas_aaa_sym;
    int  ***  bas_aab_sym;
    int  ***  bas_aba_sym;
    int **** ibas_aaa_sym;
    int **** ibas_aab_sym;
    int **** ibas_aba_sym;

    void PrintHeader();

    /// grab one- and two-electron integrals
    void GetIntegrals();

    /// read two-electron integrals from disk
    void GetTEIFromDisk();

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

    void bpsdp_ATu(SharedVector A, SharedVector u);
    void bpsdp_ATu_slow(SharedVector A, SharedVector u);
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
    
    void Update_xz();
    void Update_xz_nonsymmetric();

    void NaturalOrbitals();
    void MullikenPopulations();
    void FinalTransformationMatrix();

    // read teis from disk:
    void ReadIntegrals(double * tei,long int nmo);

    // multiplicity
    int multiplicity_;

    /// full space of integrals for MO gradient / Hessian, blocked by symmetry
    double * tei_full_sym_;
    double * oei_full_sym_;
    // gidofalvi -- modified the type of tei_full_dim_ so that it is correct for large bases 
    long int tei_full_dim_;
    int oei_full_dim_;

    /// full space D2, blocked by symmetry
    double * d2_plus_core_sym_;
    int d2_plus_core_dim_;

    /// full space D1, blocked by symmetry
    double * d1_plus_core_sym_;
    int d1_plus_core_dim_;

    /// unpack active-space density into full-space density
    void UnpackDensityPlusCore();

    /// repack rotated full-space integrals into active-space integrals 
    void RepackIntegrals();
    void RepackIntegralsDF();

    /// compute frozen core energy and adjust oeis
    void FrozenCoreEnergy();

    /// function to rotate orbitals
    void RotateOrbitals();

    double * orbopt_transformation_matrix_;
    double * orbopt_data_;
    char * orbopt_outfile_;
    bool orbopt_converged_;

    /// are we using 3-index integrals?
    bool is_df_;

    /// initialize a checkpoint file
    void InitializeCheckpointFile();

    /// write current solution and integrals to a checkpoint file
    void WriteCheckpointFile();

    /// read solution and integrals from a checkpoint file
    void ReadFromCheckpointFile();

    //EVERYTHING ANALYTICAL GRADIENTS STARTS HERE
    void grad_common_init();
    void compute_grad();
    SharedMatrix AO2SO_;
    void jk_common_init();
    SharedMatrix jk_compute_gradient();
    // Gradients
    std::map<std::string, SharedMatrix> gradients;
    std::map<std::string, SharedMatrix> ref_gradients;
    std::vector<std::string> gradient_terms;
    int natom;
    int nmo;           // Number of MOs
    int nao;           // Number of AOs
    int nao_nz;        // Number of non-zero AOs
    int ndf_nz;        // Number of non-zero DF ints in AO-basis
    int nso;           // Number of SOs
    int noccA;         // Number of alpha occupied orbitals
    int noccB;         // Number of beta occupied orbitals
    int nvirA;         // Number of alpha virtual orbitals
    int nvirB;         // Number of beta virtual orbitals
    int naoccA;        // Number of active alpha occupied orbitals                               
    int naoccB;        // Number of active beta occupied orbitals
    int namo;          // Number of active SOs
    int navirA;        // Number of active alpha virtual orbitals
    int navirB;        // Number of active beta virtual orbitals
    int nirreps;       // Number of irreducible representations
    int nshell;        // Number of shells
    int nfrzc;         // Number of frozen cores
    int nfrzv;         // Number of frozen virtuals
    int npop;          // Number of populated orbitals: npop=nmo-nfrzv
    int dimtei;        // dimension of tei in pitzer order for all integrals
    int ntri;          // square matrix dimension (nmo) -> pitzer order 
    int ntri_so;       // square matrix dimension (nso) -> pitzer order
    int ntri_ijAA;
    int ntri_ijBB;
    int ntri_abAA;
    int ntri_abBB;
    
    int nQ;          // numer of aux-basis
    int nQ_ref;      // numer of aux-basis for DF_BASIS_SCF
    int nso2_;       // nso * nso
    int naocc2AA;     // # of active OO pairs
    int naocc2AB;     // # of active OO pairs
    int naocc2BB;     // # of active OO pairs
    int navir2AA;     // # of active VV pairs
    int navir2AB;     // # of active VV pairs
    int navir2BB;     // # of active VV pairs 
    int navoAA;       // # of active VO pairs
    int navoAB;       // # of active VO pairs
    int navoBA;       // # of active VO pairs
    int navoBB;       // # of active VO pairs
    int nvoAA;        // # of all VO pairs
    int nvoAB;        // # of all VO pairs
    int nvoBA;        // # of all VO pairs 
    int nvoBB;        // # of all VO pairs
    int nocc2AA;      // # of all OO pairs
    int nocc2AB;     // # of all OO pairs
    int nocc2BB;     // # of all OO pairs
    int nvir2AA;     // # of all VV pairs
    int nvir2AB;     // # of all VV pairs
    int nvir2BB;     // # of all VV pairs  
    int nidpA;       // # of alpha independent pairs
    int nidpB;       // # of beta independent pairs
    int nidp;
    int nidp_tot;
    int idp_returnA;
    int idp_returnB;
    
    
    // Common
    double Enuc;
    double Eelec;
    double Escf;
    double Eref;
    double Etotal;
    double msd_oo_scale;
    int grad_oo_pair_idxAA(int i, int j);
    double **J_mhalf;
    
    int nvar;
    int pcg_conver;
    void grad_get_moinfo();
    void grad_corr_dens();
    void grad_c_vv_ref();
    void grad_c_ov_ref();
    void grad_c_oo_ref();
    void grad_b_vv_ref();
    void grad_b_ov_ref();
    void grad_fock_so();
    void grad_df_ref();
    void grad_trans_ref();
    void ref_int_grad();
    void tei_grad_ref();
    void oei_grad();
    void title_grad();
    void grad_G_ref();
    void grad_F_ref();
    void grad_Gaux_ref();
    void grad_so_Deriv();
    SharedMatrix symmetrize_gradient(SharedMatrix grad);

    void mo_coeff_blocks();
    void grad_formJ_ref(boost::shared_ptr<BasisSet> auxiliary_, boost::shared_ptr<BasisSet> zero);
    void grad_b_so_ref(boost::shared_ptr<BasisSet> primary_, boost::shared_ptr<BasisSet> auxiliary_, boost::shared_ptr<BasisSet> zero);
    void grad_b_oo_ref();
    void grad_remove_binary_file(int fileno);
      
    string reference;
    string reference_;
    string jobtype;
    string orb_opt_;
    string qchf_;
    string dertype;
    string opt_method;
    string hess_type;
    
    // Common
    SharedMatrix Ca_ref;              // ref_mo_vector
    SharedTensor2d CmoA;              // C(mu, p)
    SharedTensor2d CmoB;              // C(mu, p)
    SharedTensor2d Cmo_refA;          // Reference (initial) MOs
    SharedTensor2d Cmo_refB;          // Reference (initial) MOs
    SharedTensor2d CaoccA;            // C(mu, i) active
    SharedTensor2d CaoccB;            // C(mu, i) active
    SharedTensor2d CavirA;            // C(mu, a) active
    SharedTensor2d CavirB;            // C(mu, a) active
    SharedTensor2d CoccA;             // C(mu, i) all
    SharedTensor2d CoccB;             // C(mu, i) all
    SharedTensor2d CvirA;             // C(mu, a) all
    SharedTensor2d CvirB;             // C(mu, a) all
    SharedTensor2d HmoA;
    SharedTensor2d HmoB;
    SharedTensor2d FockA;
    SharedTensor2d FockB;
    
    SharedTensor2d Hso;
    SharedTensor2d Sso;
    SharedTensor2d Dso;
    SharedTensor2d DsoA;
    SharedTensor2d FsoA;
    SharedTensor2d FsoB;
    SharedTensor2d Wso;
    SharedTensor2d DQmatA;
    SharedTensor2d HooA;
    SharedTensor2d HooB;
    SharedTensor2d HovA;
    SharedTensor2d HovB;
    SharedTensor2d HvoA;
    SharedTensor2d HvoB;
    SharedTensor2d HvvA;
    SharedTensor2d HvvB;
    
    // DF Integrals                                  
    SharedTensor2d Jmhalf;             // J Metric DF_BASIS_CC (RI)
    SharedTensor2d bQso;               // b(Q|mu nu) from DF_BASIS_CC (RI)
    SharedTensor2d bQnoA;              // b(Q|mu i)
    SharedTensor2d cQso;               // c(Q|mu nu) from DF_BASIS_CC (RI) 
    SharedTensor2d cQnoA;              // c(Q|mu i)
    SharedTensor2d F_ref;              // Fock ref
    
    //Definitions not needed now but may need in future
    SharedTensor2d bQnoB;              // b(Q|mu i)        
    SharedTensor2d bQnvA;              // b(Q|mu a)
    SharedTensor2d bQnvB;              // b(Q|mu a) 
    SharedTensor2d bQooA;              // b(Q|i j) : all
    SharedTensor2d bQooB;              // b(Q|i j) : all 
    SharedTensor2d bQovA;              // b(Q|i a) : all 
    SharedTensor2d bQovB;              // b(Q|i a) : all
    SharedTensor2d bQvvA;              // b(Q|a b) : all
    SharedTensor2d bQvvB;              // b(Q|a b) : all
    SharedTensor2d bQijA;              // b(Q|i j) : active
    SharedTensor2d bQijB;              // b(Q|i j) : active
    SharedTensor2d bQiaA;              // b(Q|i a) : active
    SharedTensor2d bQiaB;              // b(Q|i a) : active
    SharedTensor2d bQabA;              // b(Q|a b) : active 
    SharedTensor2d bQabB;              // b(Q|a b) : active
    
      //still not needed for now but could be in the future
    SharedTensor2d cQnoB;              // c(Q|mu i)    
    SharedTensor2d cQnvA;              // c(Q|mu a)
    SharedTensor2d cQnvB;              // c(Q|mu a)
    SharedTensor2d cQooA;              // c(Q|i j) : all
    SharedTensor2d cQooB;              // c(Q|i j) : all  
    SharedTensor2d cQovA;              // c(Q|i a) : all
    SharedTensor2d cQovB;              // c(Q|i a) : all 
    SharedTensor2d cQvvA;              // c(Q|a b) : all
    SharedTensor2d cQvvB;              // c(Q|a b) : all
    SharedTensor2d cQijA;              // c(Q|i j) : active
    SharedTensor2d cQijB;              // c(Q|i j) : active
    SharedTensor2d cQiaA;              // c(Q|i a) : active
    SharedTensor2d cQiaB;              // c(Q|i a) : active 
    SharedTensor2d cQabA;              // c(Q|a b) : active
    SharedTensor2d cQabB;              // c(Q|a b) : active   
    //Not needed for now
    SharedTensor2d FooA;          // Fock OO block
    SharedTensor2d FooB;          // Fock oo block
    SharedTensor2d FovA;          // Fock OV block
    SharedTensor2d FovB;          // Fock ov block
    SharedTensor2d FvoA;          // Fock VO block
    SharedTensor2d FvoB;          // Fock vo block
    SharedTensor2d FvvA;          // Fock VV block
    SharedTensor2d FvvB;          // Fock vv block
    
    SharedMatrix Tso_;
    SharedMatrix Vso_;
    SharedMatrix Hso_;
    SharedMatrix Sso_;
    
    // SO basis                                   
    SharedTensor2d gQso;              // Gamma(Q|mu nu): 3-index TPDM
    SharedTensor2d gQso_ref;          // Gamma(Q|mu nu): 3-index TPDM
    SharedTensor2d gQoo;              // Gamma(Q|i i): 3-index TPDM
    SharedTensor2d gQoo_ref;          // Gamma(Q|i i): 3-index TPDM
    SharedTensor2d gQon_ref;          // Gamma(Q|i nu): 3-index TPDM
    SharedTensor2d Gaux;              // Gamma(P,Q): 2-index TPDM
    SharedTensor2d Gaux_ref;          // Gamma(P,Q): 2-index TPDM 
    SharedTensor2d dQso;              // D(Q|mu nu): 3-index OPDM for REF WFN 
    
    SharedTensor1d eps_orbA;         //alpha orbital energy
    SharedTensor1d eps_orbB;
    
    // 1D-Tensors         
    SharedTensor1d DQvecA;
    
    
    
};

}}
#endif	
      
