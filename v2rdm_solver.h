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


// greg
#include"fortran.h"

// TODO: move to psifiles.h
#define PSIF_DCC_QMO 268

namespace boost {
  template<class T> class shared_ptr;
}
namespace psi{ namespace v2rdm_casscf{

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

    boost::shared_ptr<Matrix> saveK1;

    /// S(S+1)
    double spin_squared;

    /// constrain Q2 to be positive semidefinite?
    bool constrain_q2;

    /// constrain G2 to be positive semidefinite?
    bool constrain_g2;

    /// spin adapt constraints?
    bool spin_adapt_g2;
    bool spin_adapt_q2;

    /// constraint T1 = D3 + Q3 to be positive semidefinite?
    bool constrain_t1;

    /// constraint T2 = E3 + F3 to be positive semidefinite?
    bool constrain_t2;

    /// keep d3 positive semidefinite and constrain D3->D2 mapping?
    bool constrain_d3;

    /// constrain spin?
    bool constrain_spin;

    /// symmetry product table:
    int * table;

    /// returns symmetry product for two orbitals
    int SymmetryPair(int i, int j);
    int TotalSym(int i, int j,int k, int l);
    int * symmetry;
    int * symmetry_full;
    int * symmetry_plus_core;
    int * pitzer_to_energy_order;
    int * energy_to_pitzer_order;
    int * symmetry_energy_order;
    int * pitzer_offset;           // for orbital indexing within an irrep
    int * pitzer_offset_full;      // for orbital indexing within an irrep
    int * pitzer_offset_plus_core; // for orbital indexing within an irrep

    // geminals for each symmetry:
    std::vector < std::vector < std::pair<int,int> > > gems;
    std::vector < std::vector < std::pair<int,int> > > gems_fullspace;
    std::vector < std::vector < std::pair<int,int> > > gems_plus_corespace;

    /// total number of active molecular orbitals
    int amo_;
    /// active molecular orbitals per irrep
    int * amopi_;

    int ndocc,nso,nvirt,ndoccact,nfrzc,nfrzv;

    /// total number of constraints (dimension of dual solution vector)
    long int nconstraints;

    /// total number of variables (dimension of primal solution vector)
    long int dimx;  

    /// number of auxilliary basis functions
    int nQ_;
    void ThreeIndexIntegrals(); 
    double * Qmo_;

    boost::shared_ptr<Matrix> ReadOEI();
    boost::shared_ptr<Matrix> GetOEI();

    // offsets
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

    /// maximum number of outer bpsdp iterations
    int maxiter;

    /// maximum number of outer conjugate gradient iterations
    int cg_maxiter;

    // standard vector of dimensions of each block of x
    std::vector<int> dimensions;

    int offset;

    // mapping arrays with abelian symmetry
    void BuildBasis();
    int * full_basis;

    // mapping arrays with symmetry:
    int * gems_ab;
    int * gems_aa;
    int * gems_00;
    int * gems_full;
    int * gems_plus_core;
    int *** bas_ab_sym;
    int *** bas_aa_sym;
    int *** bas_00_sym;
    int *** bas_full_sym;
    int *** ibas_ab_sym;
    int *** ibas_aa_sym;
    int *** ibas_00_sym;
    int *** ibas_full_sym;

    // triplets for each symmetry:
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
    void TEI();
    void DF_TEI();
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


    double g2timeAu,q2timeAu,d2timeAu;
    double g2timeATu,q2timeATu,d2timeATu;
    double t2timeAu,t1timeAu;
    double t2timeATu,t1timeATu;

    double escf, enuc, efrzc, efrzc1, efrzc2, tau, mu, ed, ep;

    //vectors
    SharedVector Ax;     // vector to hold A . x
    SharedVector ATy;    // vector to hold A^T . y
    SharedVector c;      // 1ei and 2ei of bpsdp
    SharedVector y;      // dual solution
    SharedVector b;      // constraint vector
    SharedVector x;      // primal solution
    SharedVector z;      // second dual solution
    
    void Update_xz();
    void NaturalOrbitals();
    void MullikenPopulations();
    void FinalTransformationMatrix();

    // read teis from disk:
    void ReadIntegrals(double * tei,int nmo);

    // multiplicity
    int multiplicity_;

    /// full space of integrals for MO gradient / Hessian, blocked by symmetry
    double * tei_full_sym;
    double * oei_full_sym;
    // gidofalvi -- modified the type of tei_full_dim so that it is correct for large bases 
    long int tei_full_dim;
    int oei_full_dim;

    /// full space D2, blocked by symmetry
    double * d2_plus_core_sym;
    int d2_plus_core_dim;

    /// full space D1, blocked by symmetry
    double * d1_plus_core_sym;
    int d1_plus_core_dim;

    /// unpack active-space density into full-space density
    void UnpackFullDensity();
    void UnpackDensityPlusCore();

    /// repack rotated full-space integrals into active-space integrals 
    void RepackIntegrals();
    void RepackIntegralsDF();

    /// function to rotate orbitals
    void RotateOrbitals();

    double * jacobi_transformation_matrix_;
    double * jacobi_data_;
    char * jacobi_outfile_;
    bool jacobi_converged_;

    /// are we using 3-index integrals?
    bool is_df_;
  
};

}}
#endif	
      
