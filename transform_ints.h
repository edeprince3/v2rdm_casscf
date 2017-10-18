#ifndef transform_ints_h
#define transform_ints_h

#include<stdlib.h>
#include<stdio.h>


void transform_ints_driver(double * int1, double *int2, double *U, int *nmopi, int *frzvpi, int nirrep, int nQ);

void transform_3index_tei(double *int2, double *U, int *nmopi, int* U_offset,
                          int* nmo_offset, int nirrep, int nQ, int max_nmopi,
                          int nmo_tot, int ngem_tot_lt, int max_num_threads);

void transform_oei(double *int1, double *U, int *nmopi, int * frzvpi, int* U_offset,
                   int* nmo_offset, int nirrep, int max_nmopi,
                   int nmo_tot, int max_num_threads);

#endif
