/*                                                                      
 *@BEGIN LICENSE

  Copyright (c) 2014, The Florida State University. All rights reserved.

 *@END LICENSE
 */


#ifndef CG_SOLVER_H
#define CG_SOLVER_H

#include<libmints/vector.h>

using namespace boost;

namespace psi{ 

typedef void (*CallbackType)(long int,SharedVector,SharedVector,void *);  

class CGSolver {
public:

    CGSolver(long int n);
    ~CGSolver();
    void preconditioned_solve(long int n,
               boost::shared_ptr<Vector> Ap,
               boost::shared_ptr<Vector>  x,
               boost::shared_ptr<Vector>  b,
               boost::shared_ptr<Vector>  precon,
               CallbackType function, void * data);
    void solve(long int n,
               boost::shared_ptr<Vector> Ap,
               boost::shared_ptr<Vector>  x,
               boost::shared_ptr<Vector>  b,
               CallbackType function, void * data);

    int total_iterations();
    void set_max_iter(int iter);
    void set_convergence(double conv);

private:

    int    n_;
    int    iter_;
    int    cg_max_iter_;
    double cg_convergence_;
    boost::shared_ptr<Vector> p;
    boost::shared_ptr<Vector> r;
    boost::shared_ptr<Vector> z;

};

} // end of namespace

#endif
