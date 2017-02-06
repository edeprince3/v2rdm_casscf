cmake_policy(SET CMP0048 NEW)  # project_VERSION* variables populated from project(... VERSION x.x.x) string
project(v2rdm_casscf
        VERSION 0.1
        LANGUAGES CXX Fortran)
set(v2rdm_casscf_AUTHORS      "A. Eugene DePrince III")
set(v2rdm_casscf_DESCRIPTION  "Variational 2-RDM-driven CASSCF plugin to Psi4")
set(v2rdm_casscf_URL          "github.com/edeprince3/v2rdm_casscf")
set(v2rdm_casscf_LICENSE      "GPL 2+")

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(psi4 1.0 REQUIRED)

# <<<  Options  >>>

include(psi4OptionsTools)

option_with_print(ENABLE_OPENMP "Enable threadsafe linking to OpenMP parallelized programs." ON)
option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)
if (APPLE AND (CMAKE_CXX_COMPILER_ID MATCHES GNU))
    option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" OFF)
else ()
    option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON
                      "-xHost" "-march=native")
endif()

#   install alongside psi4 Py module by default, but overrideable
get_filename_component(_tmp ${psi4_LIBRARY} DIRECTORY)
get_filename_component(psi4_PYMOD_INSTALL_FULL_LIBDIR ${_tmp} DIRECTORY)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${psi4_PYMOD_INSTALL_FULL_LIBDIR} CACHE PATH "Install path" FORCE)
endif()
message(STATUS "Plugin install: ${CMAKE_INSTALL_PREFIX}")

include(custom_cxxstandard)
include(autocmake_omp)
include(custom_static_library)

set(PN ${PROJECT_NAME})

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
    focas_data.F90
    focas_diis.F90
    focas_driver.F90
    focas_energy.F90
    focas_exponential.F90
    focas_genfock.F90
    focas_gradient.F90
    focas_hessian.F90
    focas_interface.F90
    focas_redundant.F90
    focas_semicanonical.F90
    focas_transform_driver.F90
    focas_transform_oeints.F90
    focas_transform_teints.F90
)

target_link_libraries(v2rdm_casscf PRIVATE ${LIBC_INTERJECT})

# <<<  Install  >>>

install(TARGETS v2rdm_casscf
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${PN})

install(FILES __init__.py pymodule.py LICENSE README.md
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${PN})

install(DIRECTORY tests/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${PN}/tests
        FILES_MATCHING PATTERN "v2rdm*/input.dat"
        PATTERN "old-tests" EXCLUDE)
