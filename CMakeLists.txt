#
# @BEGIN LICENSE
#
# v2rdm_casscf by A. Eugene DePrince III, a plugin to:
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

cmake_policy(SET CMP0048 NEW)  # project_VERSION* variables populated from project(... VERSION x.x.x) string
project(v2rdm_casscf
        VERSION 0.1
        LANGUAGES CXX Fortran)
set(v2rdm_casscf_AUTHORS      "A. Eugene DePrince III")
set(v2rdm_casscf_DESCRIPTION  "Variational 2-RDM-driven CASSCF plugin to Psi4")
set(v2rdm_casscf_URL          "github.com/edeprince3/v2rdm_casscf")
set(v2rdm_casscf_LICENSE      "GPL 3.0")

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

find_package(psi4 1.0 REQUIRED)

# <<<  Options  >>>

include(psi4OptionsTools)

option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)
option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON
                  "-xHost" "-march=native")

#   install alongside psi4 Py module by default, but overrideable
get_filename_component(_tmp ${psi4_LIBRARY} DIRECTORY)
get_filename_component(psi4_PYMOD_INSTALL_FULL_LIBDIR ${_tmp} DIRECTORY)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${psi4_PYMOD_INSTALL_FULL_LIBDIR} CACHE PATH "Install path" FORCE)
endif()
message(STATUS "Plugin install: ${CMAKE_INSTALL_PREFIX}")

include(custom_cxxstandard)
include(autocmake_static_library)

# <<<  Build  >>>

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -heap-arrays")
endif()

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

target_link_libraries(v2rdm_casscf PRIVATE ${LIBC_INTERJECT})

# <<<  Install  >>>

install(TARGETS v2rdm_casscf
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME})

install(FILES __init__.py pymodule.py LICENSE README.md
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME})

install(DIRECTORY tests/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/tests
        FILES_MATCHING PATTERN "v2rdm*/input.dat"
        PATTERN "old-tests" EXCLUDE)
