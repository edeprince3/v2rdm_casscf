#ifndef FORTRAN_H
#define FORTRAN_H

/**
 * wrappers to Greg's Fortran functions
 */

#ifndef FC_SYMBOL
#define FC_SYMBOL 2
#endif

#if   FC_SYMBOL==1
#define F77NAME(x) x
#elif FC_SYMBOL==2
#define F77NAME(x) x##_
#endif

// for greg
extern "C" {
    void F77NAME(focas_interface)(double*jacobi_transformation_matrix,
                                   double*oei_full_sym_,
                                   int &oei_full_sym__dim,
                                   double*tei_full_sym_,
                                   long int &tei_full_sym__dim,
                                   double*d1_full_sym,
                                   int &d1_full_sym_dim,
                                   double*d2_full_sym,
                                   int &d2_full_sym_dim,
                                   int *symmetry_energy_order,
                                   int &nfrzc,
                                   int &amo_,
                                   int &nfrzv,
                                   int &nirrep_,
                                   double*jacobi_data,
                                   char*jacobi_file);
};
inline void Jacobi(double*jacobi_transformation_matrix,
                   double*oei_full_sym_,
                   int &oei_full_sym__dim,
                   double*tei_full_sym_,
                   long int &tei_full_sym__dim,
                   double*d1_full_sym,
                   int &d1_full_sym_dim,
                   double*d2_full_sym,
                   int &d2_full_sym_dim,
                   int *symmetry_energy_order,
                   int &nfrzc,
                   int &amo_,
                   int &nfrzv,
                   int &nirrep_,
                   double*jacobi_data,
                   char*jacobi_file){
    F77NAME(focas_interface)(jacobi_transformation_matrix,oei_full_sym_,oei_full_sym__dim,tei_full_sym_,tei_full_sym__dim,d1_full_sym,d1_full_sym_dim,d2_full_sym,d2_full_sym_dim,symmetry_energy_order,nfrzc,amo_,nfrzv,nirrep_,jacobi_data,jacobi_file);
};

#endif
