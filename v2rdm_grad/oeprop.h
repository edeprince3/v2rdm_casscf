#ifndef _libcim_oeprop_h
#define _libcim_oeprop_h

#include <set>
#include <vector>
#include <map>
#include <string>

#include <libmints/typedefs.h>
#include <libmints/vector3.h>

namespace boost {
template<class T> class shared_ptr;
// Forward declarations for boost.python used in the extract_subsets
}

namespace psi {

class Matrix;
class Wavefunction;
class IntegralFactory;
class MatrixFactory;
class BasisSet;

/**
* Modified OEProp object, computes arbitrary expectation values (scalars)
* analyses (typically vectors).  For CIM, I just want to be able to access
* Mulliken charges.
**/
class MyOEProp : public OEProp {
protected:
    SharedMatrix Da_ao_custom(boost::shared_ptr<Matrix>Da_so);
public:
    /// Constructor, uses globals and Process::environment::reference wavefunction
    MyOEProp();
    /// Destructor
    virtual ~MyOEProp();
    /// Compute Mulliken Charges using my Da and Db
    void compute_mulliken_charges_custom(boost::shared_ptr<Matrix>Da_so,boost::shared_ptr<Matrix>Db_so);
};

}
#endif
