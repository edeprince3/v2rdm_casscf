/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
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
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 * 
 *@END LICENSE
 *
 */

#ifndef _v2rdmsolver_oeprop_h
#define _v2rdmsolver_oeprop_h

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
