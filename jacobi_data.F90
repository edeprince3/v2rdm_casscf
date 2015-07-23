module jacobi_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! intfull -- 1e and 2e integrals in chemist order sorted according to irrep
!         -- the array intull%redind maps the index for orbital pair (ie,je) {energy order}
!            to the symmetry-reduced index (is,js) {symmetry order}
! denfull -- 1-e and 2e spatial scaled density matrices; off diagonals are scaled by 2 such that E=h.d1+g.d2
! symmetries    -- symmetry labels (Cotton ordering) [ntot]
! rotinds -- list of rotation orbital pairs [2*npair]
! group_mult_tab -- group multiplication table [8,8]
! symdim -- number of geminals per irrep [nirrep]
! symnnz -- number of nonzer geminal overlaps per irrep [nirrep]

! ncore -- total number of core orbitals
! nact --  total number of active orbitals
! nvirt -- total number of virtual orbitals
! npair -- total number of orbital pairs
! ntsweep -- total number of sweeps performed
! ntrot -- total number of orbital pairs updated
! converged -- integer flag for convergence criteria met
! delrot -- total energy change due to orbital rotations
! angtol -- tolerance for rotating angles 
! orbital  1 2 3 4 5
! symmetry a b a a b
! redind array = (
! 1_a
! 1_b 2_a
! 3_a 2_b 4_a    
! 5_a 3_b 6_a 7_a
! 4_b 8_a 5_b 6_b 9_a)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  implicit none

  integer, parameter :: wp = selected_real_kind(10)  ! general working precision kind value
  integer, parameter :: fid=123456

  real(wp) :: mcetol
  real(wp) :: angtol
  real(wp) :: dele
  real(wp) :: delrot
  integer :: nproc
  integer :: ntrot
  integer :: ntsweep
  integer :: converged
  integer :: nirrep
  integer :: nfrozen
  integer :: ncore
  integer :: nact
  integer :: nvirt
  integer :: ntot
  integer :: npair
  integer :: maxind_nnz
  integer :: aarot
  logical :: print_

!---------------
! integral lists
!---------------

  type int2_sub
    real(wp), allocatable :: val(:)      ! 2-E integrals for a given irrep                      (nnz)
    integer :: nnz                   ! number of nonnzero LT 2-E integrals for this irrep
  end type int2_sub
  type int_info
    real(wp), allocatable:: e1(:)       ! 1-E integrals                                        (nnz1)
    type(int2_sub), allocatable :: e2(:) ! 2-E integrals for each irrep                         (nirrep)
  end type int_info

  type(int_info) :: intfull,denfull

  integer, allocatable :: redind(:,:),symdim_i(:),symnnz_i(:),symdim_d(:),symnnz_d(:)
  integer, allocatable :: symmetries(:)
  integer, allocatable :: rotinds(:)
  real(wp), allocatable :: mo_coeff(:,:)


!---------------------------
! group multiplication table
!---------------------------
  integer, parameter :: group_mult_tab(8,8) = reshape( (/ &  ! irrep multiplication table
      & 1,2,3,4,5,6,7,8, &
      & 2,1,4,3,6,5,8,7, &
      & 3,4,1,2,7,8,5,6, &
      & 4,3,2,1,8,7,6,5, &
      & 5,6,7,8,1,2,3,4, &
      & 6,5,8,7,2,1,4,3, &
      & 7,8,5,6,3,4,1,2, &
      & 8,7,6,5,4,3,2,1  /), (/8,8/) )

end module jacobi_data
