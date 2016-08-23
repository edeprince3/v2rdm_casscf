#
# @BEGIN LICENSE
#
# . by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(v2rdm_casscf CXX Fortran)

find_package(Psi4Plugins 1.0 REQUIRED)

add_psi4_plugin(v2rdm_casscf
    basis.cc
    blas.cc
    cg_solver.cc
    checkpoint.cc
    d2.cc
    d3.cc
    d4.cc
    diis.cc
    g2.cc
    oei.cc
    q2.cc
    sortintegrals.cc
    t1.cc
    t2.cc
    tei.cc
    threeindexintegrals.cc
    update_primal.cc
    update_transformation_matrix.cc
    v2rdm_casscf.cc
    v2rdm_solver.cc
    write_3pdm.cc
    write_tpdm.cc
    focas_data.F90             focas_hessian.F90
    focas_diis.F90             focas_interface.F90
    focas_driver.F90           focas_redundant.F90
    focas_energy.F90           focas_semicanonical.F90
    focas_exponential.F90      focas_transform_driver.F90
    focas_genfock.F90          focas_transform_oeints.F90
    focas_gradient.F90         focas_transform_teints.F90
)

