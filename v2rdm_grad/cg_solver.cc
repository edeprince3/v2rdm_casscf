/*                                                                      
 *@BEGIN LICENSE

  Copyright (c) 2014, The Florida State University. All rights reserved.

 *@END LICENSE
 */


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
#include "cg_solver.h"

using namespace boost;

namespace psi{ 

CGSolver::CGSolver(long int n) {
    n_              = n;
    iter_           = 0;
    cg_max_iter_    = 10000;
    cg_convergence_ = 1e-6;
    p = boost::shared_ptr<Vector>(new Vector(n));
    r = boost::shared_ptr<Vector>(new Vector(n));
    //z = boost::shared_ptr<Vector>(new Vector(n));
}
CGSolver::~CGSolver(){
}
void CGSolver::set_max_iter(int iter) {
    cg_max_iter_ = iter;
}
void CGSolver::set_convergence(double conv) {
    cg_convergence_ = conv;
}

void CGSolver::preconditioned_solve(long int n,
                    boost::shared_ptr<Vector> Ap, 
                    boost::shared_ptr<Vector>  x, 
                    boost::shared_ptr<Vector>  b, 
                    boost::shared_ptr<Vector>  precon, 
                    CallbackType function, void * data) {

    if ( n != n_ ) {
        throw PsiException("Warning: dimension does not match dimension from initialization",__FILE__,__LINE__);
    }

    double * p_p = p->pointer();
    double * r_p = r->pointer();
    double * z_p = z->pointer();

    double alpha = 0.0;
    double beta  = 0.0;

    // call some function to evaluate A.x.  Result in Ap
    function(n,Ap,x,data);

    double * b_p      = b->pointer();
    double * x_p      = x->pointer();
    double * Ap_p     = Ap->pointer();
    double * precon_p = precon->pointer();

    for (int i = 0; i < n; i++) {
        r_p[i] = b_p[i] - Ap_p[i];
        z_p[i] = precon_p[i] * r_p[i];
    }
    C_DCOPY(n,z_p,1,p_p,1);

    iter_ = 0;
    do {

        // call some function to evaluate A.p.  Result in Ap
        function(n,Ap,p,data);

        double rz  = C_DDOT(n_,r_p,1,z_p,1);
        double pap = C_DDOT(n_,p_p,1,Ap_p,1);
        double alpha = rz / pap;
        C_DAXPY(n_,alpha,p_p,1,x_p,1);
        C_DAXPY(n_,-alpha,Ap_p,1,r_p,1);

        // if r is sufficiently small, then exit loop
        double rrnew = C_DDOT(n_,r_p,1,r_p,1);
        double nrm = sqrt(rrnew);
        if ( nrm < cg_convergence_ ) break;

        for (int i = 0; i < n; i++) {
            z_p[i] = precon_p[i] * r_p[i];
        }
        double rznew  = C_DDOT(n_,r_p,1,z_p,1);
        double beta = rznew/rz;

        C_DSCAL(n_,beta,p_p,1);
        C_DAXPY(n_,1.0,z_p,1,p_p,1);

        iter_++;

    }while(iter_ < cg_max_iter_ );
}

void CGSolver::solve(long int n,
                    boost::shared_ptr<Vector> Ap, 
                    boost::shared_ptr<Vector>  x, 
                    boost::shared_ptr<Vector>  b, 
                    CallbackType function, void * data) {

    if ( n != n_ ) {
        throw PsiException("Warning: dimension does not match dimension from initialization",__FILE__,__LINE__);
    }

    double * p_p = p->pointer();
    double * r_p = r->pointer();

    double alpha = 0.0;
    double beta  = 0.0;

    // call some function to evaluate A.x.  Result in Ap
    function(n,Ap,x,data);

    double * b_p  = b->pointer();
    double * x_p  = x->pointer();
    double * Ap_p = Ap->pointer();

    for (int i = 0; i < n; i++) {
        r_p[i] = b_p[i] - Ap_p[i];
    }
    C_DCOPY(n,r_p,1,p_p,1);

    iter_ = 0;
    do {

        // call some function to evaluate A.p.  Result in Ap
        function(n,Ap,p,data);

        double rr  = C_DDOT(n_,r_p,1,r_p,1);
        double pap = C_DDOT(n_,p_p,1,Ap_p,1);
        double alpha = rr / pap;
        C_DAXPY(n_,alpha,p_p,1,x_p,1);
        C_DAXPY(n_,-alpha,Ap_p,1,r_p,1);

        // if r is sufficiently small, then exit loop
        double rrnew = C_DDOT(n_,r_p,1,r_p,1);
        double nrm = sqrt(rrnew);
        double beta = rrnew/rr;
        if ( nrm < cg_convergence_ ) break;

        C_DSCAL(n_,beta,p_p,1);
        C_DAXPY(n_,1.0,r_p,1,p_p,1);

        iter_++;

    }while(iter_ < cg_max_iter_ );
}

int CGSolver::total_iterations() {
    return iter_;
}


}// end of namespace
