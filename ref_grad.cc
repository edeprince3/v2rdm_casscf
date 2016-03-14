/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/** Standard library includes */
#include <fstream>
#include <psifiles.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <libmints/sieve.h>

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>

#include <lib3index/3index.h>
#include <libciomr/libciomr.h>

#include<libpsio/psio.hpp>
#include<libmints/wavefunction.h>

#include "GRAD_tensors.h"
#include "GRAD_arrays.h"
#include "v2rdm_solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace std;

namespace psi{
  namespace v2rdm_casscf{


    void v2RDMSolver::title_grad()
    {
      outfile->Printf("\n");
      outfile->Printf(" ****************************************************************************** \n");
      outfile->Printf(" ****************************************************************************** \n");
      outfile->Printf(" ****************************************************************************** \n");
      outfile->Printf(" ********                                                              ******** \n");
      outfile->Printf(" ********                      DF-v2RDM-GRAD                           ******** \n");
      outfile->Printf(" ********           A General Analytic Gradients Code                  ******** \n");
      outfile->Printf(" ********              for Density-Fitted Methods                      ******** \n");
      outfile->Printf(" ********                  by Fosso Tande                              ******** \n") ;
      outfile->Printf(" ********                                                              ******** \n");
      outfile->Printf(" ****************************************************************************** \n");
      outfile->Printf(" ****************************************************************************** \n");
      outfile->Printf(" ****************************************************************************** \n");
      outfile->Printf("\n");

    }//       

    
    
    void v2RDMSolver::compute_dfgradient()
    {
      title_grad();
      df_ref();
      oei_grad_ref();
      form_ref_Udens();
      form_ref_dens();
      metric_derivative();
      df3index_derivative();
      reduce_gradient();
      
    }//end compute energy
    
    void v2RDMSolver::df_ref()
    {
      //Get maximum number of threads for DF

      
      // get ntri from sieve
      boost::shared_ptr<ERISieve> sieve (new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
      const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
      long int ntri_cd = function_pairs.size(); //In lowetr triangular form
      
      // read integrals from disk if they were generated in the SCF
      if ( options_.get_str("SCF_TYPE") == "DF") {
	outfile->Printf("\tReading DF integrals from disk ...\n");
	
	boost::shared_ptr<BasisSet> primary = BasisSet::pyconstruct_orbital(molecule(), "BASIS", options_.get_str("BASIS"));
	boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(molecule(),"DF_BASIS_SCF", 
					       options_.get_str("DF_BASIS_SCF"), "JKFIT", options_.get_str("BASIS"), primary->has_puream());
	boost::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
	naux = auxiliary->nbf();
		
	// ntri comes from sieve above
	boost::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",naux,ntri_cd));
	double** Qmnp = Qmn->pointer();
	//PSIF_DFSCF_BJ = B Matrix containing 3-index tensor in AOs with J^-1/2 for use with DF-SCF
	psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
	psio_->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri_cd * naux);
	psio_->close(PSIF_DFSCF_BJ,1);
	
	cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", naux, nso_*nso_));
	bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", naux, nso_, nso_));
	
	for (long int mn = 0; mn < ntri_cd; mn++) {
	  long int m = function_pairs[mn].first;
	  long int n = function_pairs[mn].second;
	  for (long int P = 0; P < naux; P++) {
 	    bQso->set(P, (m*nso_) + n, Qmnp[P][mn]);
	    bQso->set(P, (n*nso_) + m, Qmnp[P][mn]);
	  }
	}
	//bQso->write(psio_, PSIF_DFOCC_INTS, true, true);
	
	// Form J^-1/2
	timer_on("Form J");
	formJ(auxiliary, zero);
	timer_off("Form J");
	cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0);
      }// end if ( options_.get_str("SCF_TYPE") == "DF" ) 
      /*
      // read integrals from disk if they were generated in the SCF
      else if ( options_.get_str("SCF_TYPE") == "CD") {
	//Get number of threads for DF
	outfile->Printf("\tReading Cholesky vectors from disk ...\n");
	naux = Process::environment.globals["NAUX (SCF)"];
	boost::shared_ptr<BasisSet> primary = BasisSet::pyconstruct_orbital(molecule(), "BASIS", options_.get_str("BASIS"));
	boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(molecule(),"DF_BASIS_SCF", 
					        options_.get_str("DF_BASIS_SCF"), "JKFIT", options_.get_str("BASIS"), primary->has_puream());
	boost::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
	outfile->Printf("\tCholesky decomposition threshold: %8.2le\n", options_.get_double("CHOLESKY_TOLERANCE"));
	
	// ntri comes from sieve above
	boost::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",naux,ntri_cd));
	double** Qmnp = Qmn->pointer();
	psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
	psio_->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri_cd * naux);
	psio_->close(PSIF_DFSCF_BJ,1);
	
	cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", naux, nso_*nso_));
	bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", naux, nso_, nso_));
	
	for (long int mn = 0; mn < ntri_cd; mn++) {
	  long int m = function_pairs[mn].first;
	  long int n = function_pairs[mn].second;
	  for (long int P = 0; P < naux; P++) {
	    bQso->set(P, (m*nso_) + n, Qmnp[P][mn]);
	    bQso->set(P, (n*nso_) + m, Qmnp[P][mn]);
	  }
	}
	outfile->Printf("\t DONE READING Cholesky decomposition vector\n");
	// Form J^-1/2
	timer_on("Form J");
	formJ(auxiliary, zero);
	outfile->Printf("\t DONE computing Jmhalf\n");
	timer_off("Form J");
	cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0); //Jacob added
	outfile->Printf("\t DONE computing CQoo\n");
      }// end else if ( options_.get_str("SCF_TYPE") == "CD" ) */
      
      else {
	//Get number of threads for DF
	df_ints_num_threads_ = options_.get_int("DF_INTS_NUM_THREADS");
	// Read in the basis set informations
	boost::shared_ptr<BasisSet> auxiliary_ = BasisSet::pyconstruct_auxiliary(reference_wavefunction_->molecule(),
						 "DF_BASIS_SCF", options_.get_str("DF_BASIS_SCF"), "JKFIT", options_.get_str("BASIS"));
	boost::shared_ptr<BasisSet> primary_ = BasisSet::pyconstruct_orbital(reference_wavefunction_->molecule(),
                                               "BASIS", options_.get_str("BASIS"));
	boost::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());

	// Read number of auxilary basis
	naux = auxiliary_->nbf();
	
	// Form J^-1/2
	cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", naux, nso_*nso_));
	bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", naux, nso_,nso_));
	
	timer_on("Form J");
	formJ(auxiliary_, zero);
	timer_off("Form J");
	
	// Form B(Q,mu nu)
	timer_on("Form B(Q,munu)");
	formB(primary_, auxiliary_, zero);
	timer_off("Form B(Q,munu)");
	cQso->gemm(true, false, Jmhalf, bQso, 1.0, 0.0); 
	outfile->Printf("\t rms CQmn %12.5lf \n",cQso->rms()); 
      }
    }

    //One electron integral derivative begins
    void v2RDMSolver::oei_grad_ref()
    {      
      gradient_terms.push_back("Nuclear");
      gradient_terms.push_back("Kinetic");
      gradient_terms.push_back("Potential");
      gradient_terms.push_back("Overlap");
      gradient_terms.push_back("Metric");
      gradient_terms.push_back("3index");
      
      outfile->Printf("ref_grad: OPDM %12.5lf\n",Dt->rms());
      // => Sizings <= //
      int natom = molecule_->natom();
      int nalpha = nalpha_;
      int nbeta = nbeta_;
      
      // => Nuclear Gradient <= //
      gradients["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1().clone());
      gradients["Nuclear"]->set_name("Nuclear Gradient");
      gradients["Nuclear"]->print_atom_vector(); 
      // => Kinetic Gradient <= //
      timer_on("Grad: T");
      {
	double** Dp = Dt->pointer();
	
	gradients["Kinetic"] = SharedMatrix(gradients["Nuclear"]->clone());
	gradients["Kinetic"]->set_name("Kinetic Gradient");
	gradients["Kinetic"]->zero();
	double** Tp = gradients["Kinetic"]->pointer();
	
	// Kinetic derivatives
	boost::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));
	const double* buffer = Tint->buffer();   
	
	for (int P = 0; P < basisset_->nshell(); P++) {
	  for (int Q = 0; Q <= P; Q++) {
	    
	    Tint->compute_shell_deriv1(P,Q);
            
	    int nP = basisset_->shell(P).nfunction();
	    int oP = basisset_->shell(P).function_index();
	    int aP = basisset_->shell(P).ncenter();
	    
	    int nQ = basisset_->shell(Q).nfunction();
	    int oQ = basisset_->shell(Q).function_index();
	    int aQ = basisset_->shell(Q).ncenter();
	      
	    int offset = nP * nQ;
	    const double* ref = buffer;
	    double perm = (P == Q ? 1.0 : 2.0);
            
	    // Px 
	    for (int p = 0; p < nP; p++) {
	      for (int q = 0; q < nQ; q++) {
		Tp[aP][0] += perm * Dp[p + oP][q + oQ] * (*ref++);
	      }
	    }
              
	    // Py 
	    for (int p = 0; p < nP; p++) {
	      for (int q = 0; q < nQ; q++) {
		Tp[aP][1] += perm * Dp[p + oP][q + oQ] * (*ref++);
	      }
	    }
            
	    // Pz 
	    for (int p = 0; p < nP; p++) {
	      for (int q = 0; q < nQ; q++) {
		Tp[aP][2] += perm * Dp[p + oP][q + oQ] * (*ref++);
	      }
	    }
	    
	    // Qx 
	    for (int p = 0; p < nP; p++) {
	      for (int q = 0; q < nQ; q++) {
		  Tp[aQ][0] += perm * Dp[p + oP][q + oQ] * (*ref++);
	      }
	    }
            
	    // Qy 
	    for (int p = 0; p < nP; p++) {
	      for (int q = 0; q < nQ; q++) {
		Tp[aQ][1] += perm * Dp[p + oP][q + oQ] * (*ref++);
	      }
	    }
            
	    // Qz 
	    for (int p = 0; p < nP; p++) {
	      for (int q = 0; q < nQ; q++) {
		Tp[aQ][2] += perm * Dp[p + oP][q + oQ] * (*ref++);
	      }
	    }
	  }
	} 
      }
      timer_off("Grad: T");
      gradients["Kinetic"]->print_atom_vector(); 
      // => Potential Gradient <= //
      timer_on("Grad: V");
      {
	double** Dp = Dt->pointer();
	
	gradients["Potential"] = SharedMatrix(gradients["Nuclear"]->clone());
	gradients["Potential"]->set_name("Potential Gradient");
	gradients["Potential"]->zero();
	
	// Thread count
	int threads = 1;
#ifdef _OPENMP
	threads = omp_get_max_threads();
#endif
	
	// Potential derivatives
	std::vector<boost::shared_ptr<OneBodyAOInt> > Vint;
	std::vector<SharedMatrix> Vtemps;
	  for (int t = 0; t < threads; t++) { 
	    Vint.push_back(boost::shared_ptr<OneBodyAOInt>(integral_->ao_potential(1)));
	    Vtemps.push_back(SharedMatrix(gradients["Potential"]->clone()));
	  }
	  
	  // Lower Triangle
	  std::vector<std::pair<int,int> > PQ_pairs;
	  for (int P = 0; P < basisset_->nshell(); P++) {
	    for (int Q = 0; Q <= P; Q++) {
	      PQ_pairs.push_back(std::pair<int,int>(P,Q));
	    }
	  }
	  
#pragma omp parallel for schedule(dynamic) num_threads(threads)
	  
	  for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
	    
	    int P = PQ_pairs[PQ].first;
	    int Q = PQ_pairs[PQ].second;
	    
	    int thread = 0;
#ifdef _OPENMP
	    thread = omp_get_thread_num();
#endif
	    
	    Vint[thread]->compute_shell_deriv1(P,Q);
	    const double* buffer = Vint[thread]->buffer();
	    
	    int nP = basisset_->shell(P).nfunction();
	    int oP = basisset_->shell(P).function_index();
	    int aP = basisset_->shell(P).ncenter();
	    
	    int nQ = basisset_->shell(Q).nfunction();
	    int oQ = basisset_->shell(Q).function_index();
	    int aQ = basisset_->shell(Q).ncenter();
	    
	    double perm = (P == Q ? 1.0 : 2.0);
	    
	    double** Vp = Vtemps[thread]->pointer();
	    
	    for (int A = 0; A < natom; A++) {
	      
	      const double* ref0 = &buffer[3 * A * nP * nQ + 0 * nP * nQ];
	      const double* ref1 = &buffer[3 * A * nP * nQ + 1 * nP * nQ];
	      const double* ref2 = &buffer[3 * A * nP * nQ + 2 * nP * nQ];
	      for (int p = 0; p < nP; p++) {
		for (int q = 0; q < nQ; q++) {
		  double Vval = perm * Dp[p + oP][q + oQ];
		  Vp[A][0] += Vval * (*ref0++);
		  Vp[A][1] += Vval * (*ref1++);
		  Vp[A][2] += Vval * (*ref2++);
		}
	      }
            }
	  } 
	  
	  for (int t = 0; t < threads; t++) { 
            gradients["Potential"]->add(Vtemps[t]);
	  }
	}
	timer_off("Grad: V");
	gradients["Potential"]->print_atom_vector(); 
	// => Overlap Gradient <= //
	timer_on("Grad: S");
	{
	  SharedMatrix W(Dt->clone());
	  W->zero();
	  SharedMatrix Wso  = reference_wavefunction_->Lagrangian();
	  MintsHelper helper(options_, 0);
	  SharedMatrix sotoao = helper.petite_list()->sotoao();
	  W->remove_symmetry(Wso,sotoao);
	  W->scale(2.0);
	  W->set_name("W");
	  double** Wp = W->pointer();
	  	  
	    
	  gradients["Overlap"] = SharedMatrix(gradients["Nuclear"]->clone());
	  gradients["Overlap"]->set_name("Overlap Gradient");
	  gradients["Overlap"]->zero();
	  double** Sp = gradients["Overlap"]->pointer();
	  
	  // Overlap derivatives
	  boost::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));
	  const double* buffer = Sint->buffer();   
	  
	  for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {
	      
	      Sint->compute_shell_deriv1(P,Q);
              
	      int nP = basisset_->shell(P).nfunction();
	      int oP = basisset_->shell(P).function_index();
	      int aP = basisset_->shell(P).ncenter();
	      
	      int nQ = basisset_->shell(Q).nfunction();
	      int oQ = basisset_->shell(Q).function_index();
	      int aQ = basisset_->shell(Q).ncenter();
	      
	      int offset = nP * nQ;
	      const double* ref = buffer;
	      double perm = (P == Q ? 1.0 : 2.0);
              
	      // Px 
	      for (int p = 0; p < nP; p++) {
		for (int q = 0; q < nQ; q++) {
		  Sp[aP][0] -= perm * Wp[p + oP][q + oQ] * (*ref++);
		}
	      }
              
	      // Py 
	      for (int p = 0; p < nP; p++) {
		for (int q = 0; q < nQ; q++) {
		  Sp[aP][1] -= perm * Wp[p + oP][q + oQ] * (*ref++);
		}
	      }
              
	      // Pz 
	      for (int p = 0; p < nP; p++) {
		for (int q = 0; q < nQ; q++) {
		  Sp[aP][2] -= perm * Wp[p + oP][q + oQ] * (*ref++);
		}
	      }
	      
	      // Qx 
	      for (int p = 0; p < nP; p++) {
		for (int q = 0; q < nQ; q++) {
		  Sp[aQ][0] -= perm * Wp[p + oP][q + oQ] * (*ref++);
		}
	      }
              
	      // Qy 
	      for (int p = 0; p < nP; p++) {
		for (int q = 0; q < nQ; q++) {
		  Sp[aQ][1] -= perm * Wp[p + oP][q + oQ] * (*ref++);
		}
	      }
              
	      // Qz 
	      for (int p = 0; p < nP; p++) {
		for (int q = 0; q < nQ; q++) {
		  Sp[aQ][2] -= perm * Wp[p + oP][q + oQ] * (*ref++);
		}
	      }
            }
	  } 
	}
	timer_off("Grad: S");
	gradients["Overlap"]->print_atom_vector();
    }



    void v2RDMSolver::reduce_gradient()
    {
      // => Total Gradient <= //
      SharedMatrix total = SharedMatrix(gradients["Nuclear"]->clone());
      total->zero();
      
      for (int i = 0; i < gradient_terms.size(); i++) {
	if (gradients.count(gradient_terms[i])) {
	  total->add(gradients[gradient_terms[i]]); 
	}
      }
      
      gradients["Total"] = total; 
      gradients["Total"]->set_name("Total Gradient");
      
      // => Final Printing <= //
      if (print_ > 1) {
	for (int i = 0; i < gradient_terms.size(); i++) {
	  if (gradients.count(gradient_terms[i])) {
	    gradients[gradient_terms[i]]->print_atom_vector(); 
	  }
	}
      } else {
	gradients["Total"]->print_atom_vector();
      }
    }
    
    //Reference TEI derivative
    void v2RDMSolver::metric_derivative()
      {  
	gradients["Metric"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1().clone());
	gradients["Metric"]->set_name("Metric Gradient");
	gradients["Metric"]->zero();
	//Metric derivative J^x
	// => Sizing <= //
	boost::shared_ptr<BasisSet> auxiliary_ = BasisSet::pyconstruct_auxiliary(reference_wavefunction_->molecule(),
						 "DF_BASIS_SCF", options_.get_str("DF_BASIS_SCF"), "JKFIT", options_.get_str("BASIS"));
	boost::shared_ptr<BasisSet> primary_ = BasisSet::pyconstruct_orbital(reference_wavefunction_->molecule(),
	                                       "BASIS", options_.get_str("BASIS"));
	boost::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
	
	int natom = primary_->molecule()->natom();
	
	double** Vp= GPQ->to_block_matrix();
	
	// => Integrals <= //
	
	boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),auxiliary_,
									 BasisSet::zero_ao_basis_set()));
	std::vector<boost::shared_ptr<TwoBodyAOInt> > Jint;

#ifdef _OPENMP
      df_ints_num_threads_ = omp_get_max_threads();
#endif

	for (int t = 0; t < df_ints_num_threads_; t++) 
	  Jint.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
	
	
	// => Temporary Gradients <= //
	
	std::vector<SharedMatrix> Ktemps;
	for (int t = 0; t < df_ints_num_threads_; t++) 
	  Ktemps.push_back(SharedMatrix(new Matrix("Ktemp", natom, 3)));
	
	std::vector<std::pair<int,int> > PQ_pairs;
	for (int P = 0; P < auxiliary_->nshell(); P++) {
	  for (int Q = 0; Q <= P; Q++) {
	    PQ_pairs.push_back(std::pair<int,int>(P,Q));
	  }
	}
	
	int nthread_df = df_ints_num_threads_;
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
	
	for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
	  
	  int P = PQ_pairs[PQ].first;
	  int Q = PQ_pairs[PQ].second;
	  
	  int thread = 0;
#ifdef _OPENMP
	  thread = omp_get_thread_num();
#endif
	  
	  Jint[thread]->compute_shell_deriv1(P,0,Q,0);
	  const double* buffer = Jint[thread]->buffer();
	  
	  int nP = auxiliary_->shell(P).nfunction();
	  int cP = auxiliary_->shell(P).ncartesian();
	  int aP = auxiliary_->shell(P).ncenter();
	  int oP = auxiliary_->shell(P).function_index();
	  
	  int nQ = auxiliary_->shell(Q).nfunction();
	  int cQ = auxiliary_->shell(Q).ncartesian();
	  int aQ = auxiliary_->shell(Q).ncenter();
	  int oQ = auxiliary_->shell(Q).function_index();
	  
	  int ncart = cP * cQ;
	  const double *Px = buffer + 0*ncart;
	  const double *Py = buffer + 1*ncart;
	  const double *Pz = buffer + 2*ncart;
	  const double *Qx = buffer + 3*ncart;
	  const double *Qy = buffer + 4*ncart;
	  const double *Qz = buffer + 5*ncart;
	  
	  double perm = (P == Q ? 1.0 : 2.0);
	  
	  double** grad_Kp;
	  
	  grad_Kp = Ktemps[thread]->pointer();
	  
	  for(int p = 0; p < nP; p++){
	    for(int q = 0; q < nQ; q++) {
	      double Vval = 0.5 * perm * Vp[p + oP][q + oQ];  //0.5 from breaking down the DF derivative
	      grad_Kp[aP][0] -= Vval * (*Px);
	      grad_Kp[aP][1] -= Vval * (*Py);
	      grad_Kp[aP][2] -= Vval * (*Pz);
	      grad_Kp[aQ][0] -= Vval * (*Qx);
	      grad_Kp[aQ][1] -= Vval * (*Qy);
	      grad_Kp[aQ][2] -= Vval * (*Qz);
	      
	      Px++;
	      Py++;
	      Pz++;
	      Qx++;
	      Qy++;
	      Qz++;
	    }
	  }
	}
	
	// => Temporary Gradient Reduction <= //
	for (int t = 0; t < df_ints_num_threads_; t++) {
	  gradients["Metric"]->add(Ktemps[t]);
	}
	gradients["Metric"]->print_atom_vector(); 
      }
	
	/*===============================================================
	 * 3-index derivative
	 ===============================================================*/

    void v2RDMSolver::df3index_derivative()
    {      
      //Metric derivative C(Q|pq)^x
      gradients["3index"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1().clone());
      gradients["3index"]->set_name("3index Gradient");
      gradients["3index"]->zero();
      // => Sizing <= //
      boost::shared_ptr<BasisSet> auxiliary_ = BasisSet::pyconstruct_auxiliary(reference_wavefunction_->molecule(),
					      "DF_BASIS_SCF", options_.get_str("DF_BASIS_SCF"), "JKFIT", options_.get_str("BASIS"));
      boost::shared_ptr<BasisSet> primary_ = BasisSet::pyconstruct_orbital(reference_wavefunction_->molecule(),
									   "BASIS", options_.get_str("BASIS"));
      //	boost::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
      
      int natom = Process::environment.molecule()->natom();
      double** GQpq_p = GQso->to_block_matrix();

      
      boost::shared_ptr<ERISieve> sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, 0.0));
      const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
      int npairs = shell_pairs.size();
      
      // => Memory Constraints <= //
      
      int max_rows;
      int maxP = auxiliary_->max_function_per_shell();
      ULI row_cost = 0L;
	row_cost = nso_ * (ULI) nso_;
	
	ULI rows = memory_ / row_cost;
	rows = (rows > naux ? naux : rows);
	rows = (rows < maxP ? maxP : rows);
	max_rows = (int) rows;
	
	// => Block Sizing <= //
	
	std::vector<int> Pstarts;
	int counter = 0;
	Pstarts.push_back(0);
	for (int P = 0; P < auxiliary_->nshell(); P++) {
	  int nP = auxiliary_->shell(P).nfunction();
	  if (counter + nP > max_rows) {
	    counter = 0;
	    Pstarts.push_back(P);
	  }
	  counter += nP;
	}
	Pstarts.push_back(auxiliary_->nshell());

#ifdef _OPENMP
      df_ints_num_threads_ = omp_get_max_threads();
#endif
	
	// => Integrals <= //
	
	boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
	std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
	for (int t = 0; t < df_ints_num_threads_; t++) 
	  eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
	
	
	// => Temporary Gradients <= //
	
	std::vector<SharedMatrix> Ktemps;
	/// std::vector<SharedMatrix> wKtemps;
	for (int t = 0; t < df_ints_num_threads_; t++) 
	  Ktemps.push_back(SharedMatrix(new Matrix("Ktemp", natom, 3)));
	
	// => Master Loop <= //
	for (int block = 0; block < Pstarts.size() - 1; block++) {
	  
	  // > Sizing < //
	  
	  int Pstart = Pstarts[block];
	  int Pstop  = Pstarts[block+1];
	  int NP = Pstop - Pstart;
	  
	  int pstart = auxiliary_->shell(Pstart).function_index();
	  int pstop  = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop ).function_index());
	  int np = pstop - pstart;
	  
	  // > Integrals < //
	  int nthread_df = df_ints_num_threads_;
	  
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
	  
	  for (long int PMN = 0L; PMN < NP * npairs; PMN++) {
	    
	    int thread = 0;
#ifdef _OPENMP
	    thread = omp_get_thread_num();
#endif
	    
	    int P =  PMN / npairs + Pstart;
	    int MN = PMN % npairs;
	    int M = shell_pairs[MN].first;
	    int N = shell_pairs[MN].second;
	    
	    eri[thread]->compute_shell_deriv1(P,0,M,N);
	    
	    const double* buffer = eri[thread]->buffer();
	    
	    int nP = auxiliary_->shell(P).nfunction();
	    int cP = auxiliary_->shell(P).ncartesian();
	    int aP = auxiliary_->shell(P).ncenter();
	    int oP = auxiliary_->shell(P).function_index() - pstart;
	    
	    int nM = primary_->shell(M).nfunction();
	    int cM = primary_->shell(M).ncartesian();
	    int aM = primary_->shell(M).ncenter();
	    int oM = primary_->shell(M).function_index();
	    
	    int nN = primary_->shell(N).nfunction();
	    int cN = primary_->shell(N).ncartesian();
	    int aN = primary_->shell(N).ncenter();
	    int oN = primary_->shell(N).function_index();
	    
	    int ncart = cP * cM * cN;
	    const double *Px = buffer + 0*ncart;
	    const double *Py = buffer + 1*ncart;
	    const double *Pz = buffer + 2*ncart;
	    const double *Mx = buffer + 3*ncart;
	    const double *My = buffer + 4*ncart;
	    const double *Mz = buffer + 5*ncart;
	    const double *Nx = buffer + 6*ncart;
	    const double *Ny = buffer + 7*ncart;
	    const double *Nz = buffer + 8*ncart;
	    
	    double perm = (M == N ? 1.0 : 2.0);
	    
	    double** grad_Kp;
	    
	    grad_Kp = Ktemps[thread]->pointer();
	    
	    for (int p = 0; p < nP; p++) {
	      for (int m = 0; m < nM; m++) {
		for (int n = 0; n < nN; n++) {
		  double Jval = perm * GQpq_p[p + oP][(m + oM) * nso_ + (n + oN)];
		  grad_Kp[aP][0] += Jval * (*Px);
		  grad_Kp[aP][1] += Jval * (*Py);
		  grad_Kp[aP][2] += Jval * (*Pz);
		  grad_Kp[aM][0] += Jval * (*Mx);
		  grad_Kp[aM][1] += Jval * (*My);
		  grad_Kp[aM][2] += Jval * (*Mz);
		  grad_Kp[aN][0] += Jval * (*Nx);
		  grad_Kp[aN][1] += Jval * (*Ny);
		  grad_Kp[aN][2] += Jval * (*Nz);
		  
		  Px++;
		  Py++;
		  Pz++;
		  Mx++;
		  My++;
		  Mz++;
		  Nx++;
		  Ny++;
		  Nz++;
		}
	      }
	    }
	  }
	}
	
	// => Temporary Gradient Reduction <= //
	
	for (int t = 0; t < df_ints_num_threads_; t++) 
	  gradients["3index"]->add(Ktemps[t]);
	
	gradients["3index"]->print_atom_vector();
	
      }// end 
      
      /*=======================================================
       *          form J(P,Q)^-1/2
       =======================================================*/          
      void v2RDMSolver::formJ(boost::shared_ptr<BasisSet> auxiliary_, boost::shared_ptr<BasisSet> zero)
      {
	//TODO for now Jmhalf is in AO basis convert to MO?
	int nthreads = 1;
#ifdef _OPENMP
	nthreads = omp_get_max_threads();
#endif
	double **J =  block_matrix(naux, naux);
	double **J_mhalfp = block_matrix(naux, naux);
	// => Integrals <= //
	boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_,zero,auxiliary_,zero));
	std::vector<boost::shared_ptr<TwoBodyAOInt> > Jint;
	std::vector<const double*> buffer;
	for (int t = 0; t < nthreads; t++) {
	  Jint.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
	  buffer.push_back(Jint[t]->buffer());
	}
	
	std::vector<std::pair<int,int> > PQ_pairs;
	for (int P = 0; P < auxiliary_->nshell(); P++) {
	  for (int Q = 0; Q <= P; Q++) {
	    PQ_pairs.push_back(std::pair<int,int>(P,Q));
	  }
	}
	
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
	
	for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
	  
	  int P = PQ_pairs[PQ].first;
	  int Q = PQ_pairs[PQ].second;
	  
	  int thread = 0;
#ifdef _OPENMP
	  thread = omp_get_thread_num();
#endif
	  
	  Jint[thread]->compute_shell(P,0,Q,0);
	  //const double* buffer = Jint[thread]->buffer();
	  
	  int nP = auxiliary_->shell(P).nfunction();
	  int oP = auxiliary_->shell(P).function_index();
	  
	  int nQ = auxiliary_->shell(Q).nfunction();
	  int oQ = auxiliary_->shell(Q).function_index();
	  
	  int index = 0;
	  for (int p = 0; p < nP; p++) {
	    for (int q = 0; q < nQ; q++, ++index) {
	      J[p + oP][q + oQ] = buffer[thread][index];
	    }
	  }
	}
	
	
	// First, diagonalize J
	// the C_DSYEV call replaces the original matrix J with its eigenvectors
	int lwork = naux * 3;
	SharedVector eigval= SharedVector(new Vector(naux));
	double* eigvalp = eigval->pointer();
	SharedVector work = SharedVector(new Vector(lwork));
	double* workp= work->pointer();
	int status = C_DSYEV('v', 'u', naux, J[0], naux, eigvalp, workp, lwork);
	if(status){
	  throw PsiException("Diagonalization of J failed", __FILE__, __LINE__);
	}
	
	// Now J contains the eigenvectors of the original J
	//SharedTensor2d JPQ_copy = SharedTensor2d(new Tensor2d("Jmetric copy <P|Q>", naux, naux));
	double **J_copy = block_matrix(naux, naux);//JPQ_copy->to_block_matrix();
	C_DCOPY(naux*naux, J[0], 1, J_copy[0], 1);
	
	// Now form J^{-1/2} = U(T)*j^{-1/2}*U,
	// where j^{-1/2} is the diagonal matrix of the inverse square roots
	// of the eigenvalues, and U is the matrix of eigenvectors of J
	for(int i=0; i<naux; ++i){
	  eigvalp[i] = (eigvalp[i] < 1.0E-10) ? 0.0 : 1.0 / sqrt(eigvalp[i]);
	  // scale one set of eigenvectors by the diagonal elements j^{-1/2}
	  C_DSCAL(naux, eigvalp[i], J[i], 1);
	}
	
	// J_mhalf = J_copy(T) * J
	C_DGEMM('t','n', naux, naux, naux, 1.0, J_copy[0], naux, J[0], naux, 0.0, J_mhalfp[0], naux);
	free_block(J);
	free_block(J_copy);
	// write J
	Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", naux, naux));
	Jmhalf->set(J_mhalfp);
	outfile->Printf("Jmhalf norm %12.5lf ",Jmhalf->rms());
	//Jmhalf->write(psio_, PSIF_DFOCC_INTS);
	//Jmhalf.reset();
	
      } // end formJ
      
      /*=======================================================
	form B^Q_pq
	=======================================================*/ 
      void v2RDMSolver::formB(boost::shared_ptr<BasisSet> primary_ 
			      ,boost::shared_ptr<BasisSet> auxiliary_
			      , boost::shared_ptr<BasisSet> zero)
      { 
	// => Sizing <= //
	int na = Ca_->colspi()[0];
	int nb = Cb_->colspi()[0];
	double** Ap = block_matrix(naux, nso_*nso_); 
	double** Bp = block_matrix(naux, nso_*nso_);
	double** Cap =Ca_->pointer();
	double** Cbp =Cb_->pointer();
	double ** BJmhalf = block_matrix(naux, nso_*nso_); 
	SharedTensor2d Amn = SharedTensor2d(new Tensor2d("B(Q|mn) ",naux, nso_*nso_));
	bQso = SharedTensor2d(new Tensor2d("C(Q|mn) ",naux, nso_*nso_)); //so
	//CQpq = SharedTensor2d(new Tensor2d("C(Q|pq) ",naux, nso_*nso_)); //mo
	
      int nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif
      
      boost::shared_ptr<ERISieve> sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, 0.0));
      const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
      int npairs = shell_pairs.size();
      
      // => Memory Constraints <= //
      int max_rows;
      max_rows = auxiliary_->nshell();
      // => Block Sizing <= //
      std::vector<int> Pstarts;
      int counter = 0;
      Pstarts.push_back(0);
      for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P).nfunction();
        if (counter + nP > max_rows) {
	  counter = 0;
	  Pstarts.push_back(P);
        }
        counter += nP;
      }
      Pstarts.push_back(auxiliary_->nshell());
      
      // => Integrals <= //
      boost::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary_, zero, primary_, primary_));
      std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
      std::vector<const double*> buffer;
      for (int t = 0; t < nthreads; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory2->eri()));
        buffer.push_back(eri[t]->buffer());
      }
      
      outfile->Printf("\t begin Master loop  \n"); 
      // => Master Loop <= //
      
      for (int block = 0; block < Pstarts.size() - 1; block++) {
	
        // > Sizing < //
	
        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;
	
        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop ).function_index());
        int np = pstop - pstart;
	
        // > Integrals < //
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
	
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {
	  
	  int thread = 0;
#ifdef _OPENMP
	  thread = omp_get_thread_num();
#endif
	  
	  int P =  PMN / npairs + Pstart;
	  int MN = PMN % npairs;
	  int M = shell_pairs[MN].first;
	  int N = shell_pairs[MN].second;
	  
	  eri[thread]->compute_shell(P,0,M,N);
	  
	  int nP = auxiliary_->shell(P).nfunction();
	  int oP = auxiliary_->shell(P).function_index();
	  
	  int nM = primary_->shell(M).nfunction();
	  int oM = primary_->shell(M).function_index();
	  
	  int nN = primary_->shell(N).nfunction();
	  int oN = primary_->shell(N).function_index();
	  
	  int index = 0;
	  for (int p = 0; p < nP; p++) {
	    for (int m = 0; m < nM; m++) {
	      for (int n = 0; n < nN; n++, index++) {
		Bp[p + oP][(m + oM) * nso_ + (n + oN)] = buffer[thread][index]; 
		Bp[p + oP][(n + oN) * nso_ + (m + oM)] = buffer[thread][index]; 
	      }
	    }
	  }
	}
      }
      double** J_mhalf= Jmhalf->to_block_matrix(); 
      //Amn->set(Bp);
      C_DGEMM('N','N', naux, nso_*nso_, naux, 1.0, J_mhalf[0], naux, Bp[0], nso_*nso_, 0.0, Ap[0], nso_*nso_);
      bQso->set(Ap);
      // C_DGEMM('N','N', naux, nso*nso, naux, 1.0, J_mhalf[0], naux, Ap[0], nso*nso, 0.0, BJmhalf[0], nso*nso);
      //CQmn->set(BJmhalf);
      free_block(Bp);
      free_block(J_mhalf);
      free_block(BJmhalf);
      free_block(Ap);
      
      
      }
      
      /*=======================================================
       *  form 3-Index  and 2-index refenrence densities
       * Ugur's paper JCP 141, 124108 (2014); doi:10.1063/1.4896235
       =======================================================*/
      void v2RDMSolver::form_ref_Udens(){
	// => Sizing <= //
	int na = Process::environment.wavefunction()->nalpha();
	int nb = Process::environment.wavefunction()->nbeta();
	SharedMatrix Ca = Ca_subset("AO");
	//	SharedMatrix Cb = Process::environment.wavefunction()->Cb();
	
	
	//transform CQmn from AO to MO
	SharedTensor2d cQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mO)", naux, nso_ * na)); 
	SharedTensor2d cQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|OO)", naux, na * na));
	SharedTensor2d CmoA = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso_, na)); //need to be removed
        
	// RHF
	// Build CmoA 
	for (int mu = 0; mu < nso_; mu++) {
	  for (int i = 0; i < na; i++) {
	    CmoA->set(mu, i, Ca->get(mu, i));
	  }
	}
	cQnoA->contract(false, false, naux * nso_, na, nso_, cQso, CmoA, 1.0, 0.0);
	cQooA->contract233(true, false, na, na, CmoA, cQnoA, 1.0, 0.0);
	
	//Constructing Urgur's reference density. How the hell did he get this
	SharedTensor2d GQpq_ref = SharedTensor2d(new Tensor2d("DF_BASIS_SCF 3-Index TPDM <O|O>", naux, na* na));
	GQpq_ref->copy(cQooA);
	GQpq_ref->scale(-2.0);
	for (int Q = 0; Q < naux; Q++) {
	  for(int i = 0; i < na; ++i) {
	    double summ = 0.0;
	    int ii = i*na +i;
	    for(int j = 0; j < na; ++j) {
	      int jj = j*na+j;
	      summ +=  cQooA->get(Q,jj);
	    }
	    GQpq_ref->add(Q,ii,4.0*summ);
	  }
	}
	
	GQso = SharedTensor2d(new Tensor2d(" G (Q|mn)", naux, nso_*nso_));
	GPQ = SharedTensor2d(new Tensor2d("2-index TPDM G(P|Q)", naux, naux));
	
	SharedTensor2d gQon_ref = SharedTensor2d(new Tensor2d("DF_BASIS_SCF G_imu^Q", naux, nso_ * na));
	gQon_ref->contract(false, true, naux * na, nso_, na, GQpq_ref, CmoA, 1.0, 0.0);
	GQso->contract233(false, false, nso_, nso_, CmoA, gQon_ref, 1.0, 0.0);
	
	GPQ->gemm(false,true, cQooA, GQpq_ref, 1.0, 0.0);// MO basis 
	//GPQ->gemm(false,true, CQmn, GQmn, 1.0, 0.0);// MO basis 
	outfile->Printf("rms GQpq %12.5lf\n",GQpq_ref->rms());
	outfile->Printf("rms GQmn %12.5lf\n",GQso->rms());
	
	outfile->Printf("rms GPQ %12.5lf\n",GPQ->rms());
      }
      /*=======================================================
       *          form G^Q_pq 3-Index density
       =======================================================*/
      void v2RDMSolver::form_ref_dens(){
	// => Sizing <= //
	int na = Process::environment.wavefunction()->nalpha();
	int nb = Process::environment.wavefunction()->nbeta();
	
	
	GQso = SharedTensor2d(new Tensor2d(" G (Q|mn)", naux, nso_*nso_));
	GPQ = SharedTensor2d(new Tensor2d("2-index TPDM G(P|Q)", naux, naux));
	SharedTensor2d G_ref = SharedTensor2d(new Tensor2d(" G temp", nso_*nso_, nso_*nso_));
	SharedMatrix Da =Da_subset("AO");
	SharedMatrix Db =Db_subset("AO");
	double** Dap = Da->pointer();
	double** Dbp = Db->pointer();
	
	double** G_refp = block_matrix(nso_*nso_, nso_*nso_);
	
       	for(int p=0;p<nso_;p++){      
	  for(int q=0;q<nso_;q++){
	    for(int r=0;r<nso_;r++){      
	      for(int s=0;s<nso_;s++){      
		double dab = Dap[p][q]*Dbp[r][s];
		double daa = Dap[p][q]*Dap[r][s] - Dap[p][s]*Dap[r][q];
		double dbb = Dbp[p][q]*Dbp[r][s] - Dbp[p][s]*Dbp[r][q];
		G_ref->set(p*nso_+q,r*nso_+s, daa + dbb+ 2.0*dab) ;
	      }
	    }
	  }
	}
	
	//3-index density
	GQso->gemm(false, false,cQso, G_ref,1.0,0.0);
	
	//2-index density
	GPQ->gemm(false,true, cQso, GQso, 1.0, 0.0);
	
      }
          
  } 
} // End Namespaces




