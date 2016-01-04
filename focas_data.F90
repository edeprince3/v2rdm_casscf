module focas_data
  implicit none

  ! *** parameters

  integer, parameter :: wp = selected_real_kind(10)                ! general working precision kind value
  integer, parameter :: ip = selected_int_kind(16)                  ! 64-bit integers (integral addressing)
  integer, parameter :: fid_ = 12345                               ! file identifier for output file
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

  real(wp), allocatable :: q_(:,:)                                 ! auxiliary matrix (asymmetric, nact*nmo storage) 
  real(wp), allocatable :: z_(:,:)                                 ! auxiliary matrix that contains cotractions of the fock_i with den1 ( nact*nmo storage)
!  real(wp), allocatable :: fock_gen_(:,:)                          ! generalized Fock matrix (nmo*nmo storage)
  real(wp), allocatable :: orbital_gradient_(:)                    ! orbital gradient
  real(wp), allocatable :: diagonal_orbital_hessian_(:)            ! diagonal elements of the orbital hessian
  real(wp), allocatable :: kappa_(:)                               ! orbital rotation parameters (lt elements of skew-symmetric matrix, npair_ storage)

  ! *** symmetry data for integrals and densities

  type sym_info
    integer, allocatable     :: ngempi(:)                          ! number of geminals per irrep
    integer, allocatable     :: nnzpi(:)                           ! number of nnz matrix elements
    integer, allocatable     :: offset(:)                          ! offset for first matrix element in this irrep
    integer, allocatable     :: gemind(:,:)                        ! symmetry-reduced index of a geminal 
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

  type diis_info
    integer :: do_diis                                             ! flag for performing DIIS updates (set internally based on max_num_diis)
    integer :: error                                               ! return code from dgesv()
    integer :: update                                              ! internal flag for performing update
    integer :: max_num_diis                                        ! maximum number of diis vectors stored
    integer :: current_index                                       ! current diis index
    real(wp), allocatable :: B(:,:)                                ! DIIS B matrix (max_num_diis+1,max_num_diis+1)
    real(wp), allocatable :: c(:)                                  ! coefficent vector for DIIS interpolation
    integer, allocatable  :: ip(:)                                 ! temporary matrix used during solution of A * x = c
    real(wp), allocatable :: dP(:,:)                               ! dP vectors (npair,max_num_diis)
    real(wp), allocatable :: P(:,:)
  end type diis_info

  type fock_info
    type(matrix_block), allocatable :: occ(:)
    type(vector_block), allocatable :: ext(:)
  end type fock_info

  type qint_info
    type(matrix_block), allocatable :: tuQ(:) 
  end type qint_info

  ! *** allocatable derived types

  type(sym_info)   :: dens_                                        ! density symmetry data
  type(sym_info)   :: ints_                                        ! integral symmetry data
  type(trans_info) :: trans_
  type(diis_info)  :: diis_
  type(fock_info)  :: fock_i_ 
  type(fock_info)  :: fock_a_
  type(qint_info)  :: qint_
 
  ! indexing derived types
  
  type(rot_info)   :: rot_pair_                                    ! info for rotation pair indexing
  type(df_info)    :: df_vars_

  ! *** allocatable orbital index arrays

  integer, allocatable :: ndocpi_(:)                               ! number of doubly occupied orbitals per irrep
  integer, allocatable :: nactpi_(:)                               ! number of active orbitals per irrep
  integer, allocatable :: nextpi_(:)                               ! number of external orbitals per irrep  
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
  integer :: log_print_                                            ! 1/0 = flag for printing iteration/info for orbtial optimization
  integer :: num_negative_diagonal_hessian_                        ! number of negative diagonal Hessian matrix elements
  integer :: use_exact_hessian_diagonal_                           ! flag to use exact expressions for the diagonal elements of the Hessian
  integer :: num_diis_vectors_
 
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
  real(wp) :: min_diag_hessian_                                    ! smallest diagonal Hessian element
 
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

    pure function df_aa_index(g,a,a_sym)
! function to return the column index of df(:,ga) where 
! both g & a are an active orbitals (LT storage)
      integer, intent(in)  :: g,a,a_sym
      integer  :: a_i,g_i,df_aa_index

      ! adjust for the number of doubly-ococcupied orbital in this irrep
      a_i = trans_%class_to_irrep_map(a)-ndocpi_(a_sym)

      ! orbital index within irrep
      g_i = trans_%class_to_irrep_map(g)-ndocpi_(a_sym)

      if (a_i.ge.g_i) then
        df_aa_index=ishft(a_i*(a_i-1),-1)+g_i
        return
      else
        df_aa_index=ishft(g_i*(g_i-1),-1)+a_i
        return
      end if

      return
    end function df_aa_index

    pure function df_ga_index(g,a,a_sym)
! function to return the column index of df(:,ga) where 
! g is a general index and a is an active index
! assumes that for each general index g, all the a indeces are stored in contiguous order
      integer, intent(in)  :: g,a,a_sym
      integer  :: a_i,g_i,df_ga_index

      ! adjust for the number of doubly-ococcupied orbital in this irrep
      a_i = trans_%class_to_irrep_map(a)-ndocpi_(a_sym)

      ! orbital index within irrep
      g_i = trans_%class_to_irrep_map(g)


      df_ga_index = ( g_i - 1 ) * nactpi_(a_sym) + a_i

      return
    end function df_ga_index

    pure function df_gd_index(g,d,d_sym)
! function to return the column index of df(:,ga) where 
! g is a general index and d is an doubly-occupied index
! assumes that for each general index g, all the d indeces are stored in contiguous order
      integer, intent(in)  :: g,d,d_sym
      integer  :: d_i,g_i,df_gd_index

      ! adjust for the number of doubly-ococcupied orbital in this irrep
      d_i = trans_%class_to_irrep_map(d)

      ! orbital index within irrep
      g_i = trans_%class_to_irrep_map(g)


      df_gd_index = ( g_i - 1 ) * ndocpi_(d_sym) + d_i

      return
    end function df_gd_index

    pure function df_pq_index(i,j)
! this function computes the two-electron index index (lower triangular reference)
! index = ii*(ii+1)/2+jj where ii=max(i,j) and jj=min(i,j)
! the ishft(k,-1) divides the value of the integer k by 2 and seems to be somewhat
! faster than the regular human-readable expression
      implicit none
      integer, intent(in) ::i,j
      integer(ip) :: i_long,j_long
      integer(ip) :: df_pq_index
      i_long = int(i,kind=ip)
      j_long = int(j,kind=ip)
      if (i_long.ge.j_long) then
        df_pq_index = i_long * ( i_long + 1 ) / 2 + j_long
        df_pq_index = df_pq_index * int(df_vars_%nQ,kind=ip)
        return
      else
        df_pq_index= j_long * ( j_long + 1 ) / 2 + i_long
        df_pq_index = df_pq_index * int(df_vars_%nQ,kind=ip)
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

end module focas_data
