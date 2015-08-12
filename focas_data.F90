module focas_data
  implicit none

  ! *** parameters

  integer, parameter :: wp = selected_real_kind(10)                ! general working precision kind value
  integer, parameter :: ip = selected_int_kind(8)                  ! 64-bit integers (integral addressing)
  integer, parameter :: fid = 12345                                ! file identifier for output file
  integer, parameter :: max_nirrep_=8                              ! maximum number of irreps
  integer, parameter :: group_mult_tab_(max_nirrep_,max_nirrep_) & ! irrep multiplication table
      & = reshape( (/    &  
      & 1,2,3,4,5,6,7,8, &
      & 2,1,4,3,6,5,8,7, &
      & 3,4,1,2,7,8,5,6, &
      & 4,3,2,1,8,7,6,5, &
      & 5,6,7,8,1,2,3,4, &
      & 6,5,8,7,2,1,4,3, &
      & 7,8,5,6,3,4,1,2, &
      & 8,7,6,5,4,3,2,1  /), (/8,8/) )

  ! *** allocatable real arrays

  real(wp), allocatable :: fock_i_(:)                              ! inactive Fock matrix (symmetric, LT nmo*(nmo+1)/2 storage)
  real(wp), allocatable :: fock_a_(:)                              ! active Fock matrix (symmetric, LT nmo*(nmo+1)/2 storage)
  real(wp), allocatable :: q_(:,:)                                 ! auxiliary matrix (asymmetric, nact*nmo storage) 
  real(wp), allocatable :: z_(:,:)                                 ! auxiliary matrix that contains cotractions of the fock_i with den1 ( nact*nmo storage)
  real(wp), allocatable :: fock_gen_(:,:)                          ! generalized Fock matrix (nmo*nmo storage)
  real(wp), allocatable :: orbital_gradient_(:)                    ! orbital gradient
  real(wp), allocatable :: diagonal_orbital_hessian_(:)            ! diagonal elements of the orbital hessian
  real(wp), allocatable :: kappa_(:)                               ! orbital rotation parameters (lt elements of skew-symmetric matrix, npair_ storage)

  ! *** symmetry data for integrals and densities

  type sym_info
    integer, allocatable :: ngempi(:)                              ! number of geminals per irrep
    integer, allocatable :: nnzpi(:)                               ! number of nnz matrix elements
    integer, allocatable :: offset(: )                             ! offset for first matrix element in this irrep
    integer, allocatable :: gemind(:,:)                            ! symmetry-reduced index of a geminal 
  end type sym_info

  ! symmetry data for transformation

  type matrix_block
    real(wp), allocatable :: val(:,:)                       
  end type matrix_block

  type vector_block
    real(wp), allocatable :: val(:)
  end type vector_block

  type trans_info
    integer, allocatable :: npairpi(:)                             ! number of orbital rotation pairs per irrep
    integer, allocatable :: nmopi(:)                               ! number of orbitals per irrep
    integer, allocatable :: offset(:)                              ! first index of orbital (with a given symmetry) in the irrep_to_class_map array
    integer, allocatable :: irrep_to_class_map(:)                  ! mapping array to map symmetry-reduced index to class-index
    integer, allocatable :: class_to_irrep_map(:)                  ! mapping array to map class-index to symmetry-reduced index
    type(matrix_block), allocatable :: u_irrep_block(:)            ! transformation matrix for a symmetry block
  end type trans_info

  type rot_info
    integer :: act_doc_type                                        ! index for active-doubly occupied orbital pair
    integer :: ext_doc_type                                        ! index for external-doubly occupied orbital pair
    integer :: act_act_type                                        ! index for active-active orbital pair
    integer :: ext_act_type                                        ! index for external-active orbital pair 
    integer :: n_tot                                               ! total number of rotation pairs
    integer :: n_ad                                                ! total number of active-doubly occupied rotation pairs
    integer :: n_aa                                                ! total number of active-active pairs
    integer :: n_ed                                                ! total number of external-doubly occupied rotation pairs
    integer :: n_ea                                                ! total number of external-active rotation pairs
    integer, allocatable :: pair_offset(:,:)
  end type rot_info

  type df_info 
    integer :: nQ                                                  !  number of auxiliary function for density-fitted integrals
    integer :: use_df_teints                                       ! flag to use density-fitted 2-e integrals
    integer, allocatable :: class_to_df_map(:)                     ! mapping array to map orbital indeces from class order to df order
  end type df_info

  ! *** allocatable derived types

  type(sym_info)   :: dens_                                        ! density symmetry data
  type(sym_info)   :: ints_                                        ! integral symmetry data
  type(trans_info) :: trans_

  ! indexing derived types
  
  type(rot_info)   :: rot_pair_                                    ! info for rotation pair indexing
  type(df_info)    :: df_vars_

  ! *** indexing arrays

  integer :: ndocpi_(max_nirrep_)                                  ! number of doubly occupied orbitals per irrep
  integer :: nactpi_(max_nirrep_)                                  ! number of active orbitals per irrep
  integer :: nextpi_(max_nirrep_)                                  ! number of external orbitals per irrep

  ! *** allocatable orbital index arrays

  integer, allocatable :: first_index_(:,:)                        ! index of first orbital in this class and irrep (nirrep,3)
  integer, allocatable :: last_index_(:,:)                         ! index of last orbital in this class and irrep (nirrep,3)
  integer, allocatable :: orb_sym_scr_(:)                          ! scratch array to store the symmetries of orbitals

  ! *** integers

  integer :: nirrep_                                               ! number of irreps in point group
  integer :: ndoc_tot_                                             ! total number of doubly occupied orbitals
  integer :: nact_tot_                                             ! total number of active orbitals
  integer :: next_tot_                                             ! total number of external orbitals
  integer :: nmo_tot_                                              ! total number of orbitals
  integer :: include_aa_rot_                                       ! 1/0 = include/do not include rotations between active-active orbtials
  integer :: nthread_use_                                          ! number of threads to use in parallel parts of the code (this is the actuaal number of threads used)
  integer :: nthread_want_                                         ! number of threads to use in parallel parts of the code ( as specified by user )

  ! *** doubles
  real(wp) :: e1_c_                                                ! core contribution to 1-e energy
  real(wp) :: e2_cc_                                               ! core contribution to 2-e energy all indeces in g(pq|rs) in \D
  real(wp) :: e1_a_                                                ! active contribution to 1-e energy 
  real(wp) :: e2_aa_                                               ! active-active contribution to 2-e energy all indeces in g(pq|rs) in \A
  real(wp) :: e2_ca_                                               ! core-active contribution to 2-e energy only 2 indeces in g(pq|rs) in \A
  real(wp) :: e_nuc_rep_                                           ! nuclear repulsion energy
  real(wp) :: e_frozen_core_                                       ! frozen core energy
  real(wp) :: e_total_                                             ! total energy
  real(wp) :: e1_total_                                            ! total 1-e energy
  real(wp) :: e2_total_                                            ! total 2-e energy
  real(wp) :: e_active_                                            ! active space energy
  real(wp) :: grad_norm_                                           ! norm of the gradient ddot(g,g)
 
  contains

    pure function pq_index(i,j)
! this function computes the two-electron index index (lower triangular reference)
! index = ii*(ii-1)/2+jj where ii=max(i,j) and jj=min(i,j)
! the ishft(k,-1) divides the value of the integer k by 2 and seems to be somewhat
! faster than the regular human-readable expression
      implicit none
      integer, intent(in) ::i,j
      integer(ip) :: pq_index
      if (i.ge.j) then
        pq_index=ishft(i*(i-1),-1)+j
        return
      else
        pq_index=ishft(j*(j-1),-1)+i
        return
      end if
    end function pq_index

    pure function df_pq_index(i,j)
! this function computes the two-electron index index (lower triangular reference)
! index = ii*(ii+1)/2+jj where ii=max(i,j) and jj=min(i,j)
! the ishft(k,-1) divides the value of the integer k by 2 and seems to be somewhat
! faster than the regular human-readable expression
      implicit none
      integer, intent(in) ::i,j
      integer(ip) :: df_pq_index
      if (i.ge.j) then
        df_pq_index= i + j * nmo_tot_
        return
      else
        df_pq_index= j + i * nmo_tot_
        return
      end if
    end function df_pq_index

    function timer()
      real(wp) :: timer,omp_get_wtime
#ifdef OMP
      timer = omp_get_wtime()
#else
      call cpu_time(timer)
#endif
    end function timer

    integer function get_nthread()
      ! simple function to set the number of threads to use in the parallel sections of the code
      integer :: max_thread,omp_get_max_threads
      
#ifdef OMP
      max_thread = omp_get_max_threads()
#else
      max_thread = 1
#endif
      get_nthread = min( max_thread , nthread_want_ )
   
      return 

    end function get_nthread

end module focas_data
