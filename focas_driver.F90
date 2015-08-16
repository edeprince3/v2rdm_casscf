module focas_driver
  use focas_data
  use focas_energy
  use focas_gradient
  use focas_hessian
  use focas_transform_driver
  use focas_exponential

  implicit none

  contains

  subroutine focas_optimize(int1,nnz_int1,int2,nnz_int2,den1,nnz_den1,den2,nnz_den2,ndocpi,nactpi,nextpi,nirrep, &
                           & gnorm_tol,dele_tol,gnorm_last,dele_last,converged,df_ints_in,nthread,fname,print_fl)
 
    ! integer input
    integer, intent(in)     :: nthread         ! number of threads to use
    integer, intent(inout)  :: converged       ! flag for convergence
    integer, intent(in)     :: df_ints_in      ! flag for using density-fitted 2-e integrals (0--> 4-index integrals, 1--> 3-index DF integrals)
    integer, intent(in)     :: nirrep          ! number of irreps in point group
    integer, intent(in)     :: nnz_int1        ! total number of nonzero 1-e integrals
    integer(ip), intent(in)     :: nnz_int2        ! total number of nonzero 2-e integrals
    integer, intent(in)     :: nnz_den1        ! total number of nonzero 1-e density elements
    integer, intent(in)     :: nnz_den2        ! total number of nonzero 2-e density elements
    integer, intent(in)     :: ndocpi(nirrep)  ! number of doubly occupied orbitals per irrep (includes frozen doubly occupied orbitals)
    integer, intent(in)     :: nactpi(nirrep)  ! number of active orbitals per irrep
    integer, intent(in)     :: nextpi(nirrep)  ! number of virtual orbitals per irrep (excluding forzen virtual orbitals) 
    integer, intent(in)     :: print_fl        ! flag for printing output information
    ! real input
    real(wp), intent(inout) :: dele_last       ! final energy change
    real(wp), intent(inout) :: gnorm_last      ! final gradient norm
    real(wp), intent(in)    :: dele_tol        ! final energy change threshold
    real(wp), intent(in)    :: gnorm_tol       ! final gradient norm threshold
    real(wp), intent(in)    :: int1(nnz_int1)  ! nonzero 1-e integral matrix elements
    real(wp), intent(in)    :: int2(nnz_int2)  ! nonzero 2-e integral matrix elements 
    real(wp), intent(in)    :: den1(nnz_den1)  ! nonzero 1-e density matrix elements
    real(wp), intent(in)    :: den2(nnz_den2)  ! nonzero 2-e density matrix elements
    ! character
    character(*), intent(in) :: fname          ! name of file to print output

    ! timing variables
    real(wp) :: t0,t1,t_ene,t_gh,t_exp,t_trans

    ! iteration variables
    real(wp) :: current_energy,last_energy,delta_energy,gradient_norm_tolerance,delta_energy_tolerance
    real(wp) :: step_size,initial_energy
    integer  :: iter,max_iter,error
  
    ! other variables
    logical :: fexist

    ! orbitals should be sorted accoring to 
    ! 1) class (doubly occupied, active, virtual)
    ! 2) within each class, orbitals should be sorted accoring to irrep
    ! 3) within each class and irrep, orbitals are sorted according to energy

    log_print_ = print_fl
    if ( log_print_ == 1 ) then
      inquire(file=fname,exist=fexist)
      if (fexist) then
        open(fid_,file=fname,status='old',position='append')
      else
        open(fid_,file=fname,status='new')
      endif
    endif


    ! set df integral flag
    df_vars_%use_df_teints  = df_ints_in

    ! calculate the total number of orbitals in space
    ndoc_tot_ = sum(ndocpi)
    nact_tot_ = sum(nactpi)
    next_tot_ = sum(nextpi)
    nmo_tot_  = ndoc_tot_+nact_tot_+next_tot_

    include_aa_rot_ = 0

    ! figure out maxmimum number of threads to use
    nthread_want_ = nthread
    nthread_use_  =get_nthread()

    ! allocate indexing arrays
    call allocate_indexing_arrays(nirrep)

    ! determine integral/density addressing arrays
    call setup_indexing_arrays(ndocpi,nactpi,nextpi)

    ! allocate transformation matrices
    call allocate_transformation_matrices() 
   
    ! determine valid orbital rotation pairs (orbital_gradient allocated upon return)
    call setup_rotation_indeces()

    ! allocate temporary Fock matrices
    call allocate_temporary_fock_matrices()

    ! allocate remaining arrays/matrices
    call allocate_initial()

    ! determine indexing arrays (needed for sorts in the integral transformation step)
    call determine_transformation_maps()

    ! set up df mapping arrays if density-fitted 2-e integrals are used
    error = 0
    if ( df_vars_%use_df_teints == 1 ) error = df_map_setup(nnz_int2)
    if ( error /= 0 ) stop

    ! print information
    if ( log_print_ == 1 ) call print_info()    

    ! **********************************************************************************
    ! *** at this point, everything is allocated and we are ready to do the optimization
    ! **********************************************************************************


    last_energy             = 0.0_wp
    gradient_norm_tolerance = gnorm_tol
    delta_energy_tolerance  = dele_tol
    max_iter                = 20
    iter                    = 0
    kappa_                  = 0.0_wp
    converged               = 0

    do 
      
      ! transform integrals
      t0 = timer()

      if ( iter /= 0 ) call transform_driver(int1,int2)      

      t1 = timer()
      t_trans = t1 - t0

      ! calculate the current energy
      t0=timer()

      call compute_energy(int1,int2,den1,den2)

      current_energy = e_total_

      if ( iter == 0 ) initial_energy = current_energy

      delta_energy = last_energy - current_energy

      t1=timer()
      t_ene = t1 - t0

      t0 = timer()
      ! construct gradient (temporary Fock matrices allocated/deallocated within routine)
      call orbital_gradient(int1,int2,den1,den2)

      ! precondition gradient with the diagonal Hessian
      call diagonal_inverse_hessian_preconditioner(orbital_gradient_,fock_i_,fock_a_,q_,z_,den1)

      ! update kappa_
      step_size=1.0_wp
      kappa_ = - step_size * orbital_gradient_

      t1 = timer()
      t_gh = t1 - t0

      ! from kappa, construct unitary transformation matrix ( U = exp(kappa) )
      t0 = timer()

      call compute_exponential(kappa_)

      t1 = timer()
      t_exp = t1 - t0

      if ( log_print_ == 1 ) then  
        write(fid_,'(a,1x,i5,1x,1(a,1x,f15.8,3x),2(a,1x,es11.3,2x),4(a,1x,f6.2,1x))')'iter:',iter,&
                   & 'E(k):',current_energy,'E(k)-E(k-1)',-delta_energy,'||g||',grad_norm_,       &
                   & 't(g+h):',t_gh,'t(E):',t_ene,'t(e^U):',t_exp,'t(tran):',t_trans
      endif

      last_energy = current_energy

      iter = iter + 1

      if ( iter == max_iter ) exit

      if ( ( delta_energy > delta_energy_tolerance ) .or. (grad_norm_ > gradient_norm_tolerance) ) cycle

      converged = 1

      exit

    end do

    if ( log_print_ == 1 ) then
      if ( converged == 1 ) then
        write(fid_,'(a)')'gradient descent converged'
      else
        write(fid_,'(a)')'gradient descent did not converge'
      endif
    endif

    dele_last      = last_energy - initial_energy
    gnorm_last     = grad_norm_

    ! debug
    ! call compute_energy(int1,int2,den1,den2)

    ! deallocate indexing arrays
    call deallocate_indexing_arrays()

    ! final deallocation
    call deallocate_final()

    if ( log_print_ == 1 ) close(fid_)

    return
  end subroutine focas_optimize

  subroutine deallocate_final()
    implicit none
    call deallocate_temporary_fock_matrices()
    if (allocated(orbital_gradient_))        deallocate(orbital_gradient_)
    if (allocated(kappa_))                   deallocate(kappa_)
    if (allocated(rot_pair_%pair_offset))    deallocate(rot_pair_%pair_offset)
    if (allocated(df_vars_%class_to_df_map)) deallocate(df_vars_%class_to_df_map)
    call deallocate_transformation_matrices()
    call deallocate_hessian_data()
    return
  end subroutine deallocate_final

  subroutine allocate_initial()
    implicit none

    ! allocate orbital gradient/kappa
    allocate(orbital_gradient_(rot_pair_%n_tot))
    allocate(kappa_(rot_pair_%n_tot))
    orbital_gradient_ = 0.0_wp
    kappa_            = 0.0_wp
 
    ! hessian data
    call allocate_hessian_data()

    return
  end subroutine allocate_initial

  integer function df_map_setup(nnz_int2)
    implicit none
    integer(ip), intent(in) :: nnz_int2
    integer(ip) :: num
    integer :: npair,i_sym,i_class,ic,idf

    df_map_setup = 1

    ! check to make sure that input integral array is of reasonable size
    npair        = nmo_tot_* ( nmo_tot_ + 1 ) / 2

    num = nnz_int2/int(npair,kind=ip)
    num = num * int(npair,kind=ip)

    if ( num /= nnz_int2 ) then
      if ( log_print_ == 1) write(fid_,'(a)')'mod(nnz_int2,nmo_tot_*(nmo_tot_+1)/2) /= 0'
      return 
    endif
 
    ! determine the number of auxiliary function 
    df_vars_%nQ  = nnz_int2/npair

    ! allocate and determine mapping array from class order to df order
    allocate(df_vars_%class_to_df_map(nmo_tot_))

    idf = 0

    do i_sym = 1 , nirrep_

      do i_class = 1 , 3

        do ic = first_index_(i_sym,i_class) , last_index_(i_sym,i_class)

          df_vars_%class_to_df_map(ic) = idf
          idf = idf + 1          

        end do

      end do

    end do

    if ( minval(df_vars_%class_to_df_map) < 0 ) then

      if ( log_print_ == 1 ) write(fid_,'(a)')'error ... min(class_to_df_map(:)) < 0 )'
      return

    end if

    df_map_setup = 0

    return
  end function df_map_setup

  subroutine setup_rotation_indeces()
    implicit none
    ! subroutine to determine nonredundant orbital pairs
    integer :: i_sym,i,j,i_class,j_class,j_class_start,j_start,npair,npair_type(4),pair_ind

    ! initialize orbital pair type indeices
    rot_pair_%act_doc_type = 1
    rot_pair_%ext_doc_type = 2
    rot_pair_%act_act_type = 3
    rot_pair_%ext_act_type = 4

    ! initialize number of rotation paors per irrep
    trans_%npairpi = 0 

    ! initialize counters for each pair type
    npair_type     = 0

    ! initialize counter for total pairs
    npair          = 0

    ! allocate rotation index offset matrix
    allocate(rot_pair_%pair_offset(nirrep_,4))

    ! loop over symmetries for i

    do i_sym = 1 , nirrep_

      ! loop over orbital classes for i

      do i_class = 1 , 3

        ! determine first class for j

        j_class_start = i_class + 1
        if ( ( include_aa_rot_ == 1 ) .and. ( i_class == 2 ) ) j_class_start = i_class

        ! loop over classes for j

        do j_class = j_class_start , 3

          ! figure out the type of pair counter

          if ( i_class == 1 ) then

            pair_ind = rot_pair_%act_doc_type
  
            if ( j_class == 3 ) pair_ind = rot_pair_%ext_doc_type

          else

            pair_ind = rot_pair_%act_act_type

            if ( j_class == 3 ) pair_ind = rot_pair_%ext_act_type

          endif

          ! save offset for this rotation type

          rot_pair_%pair_offset(i_sym,pair_ind) = npair

          ! loop over i indeces

          do i = first_index_(i_sym,i_class) , last_index_(i_sym,i_class)

            ! determine first index for j
 
            j_start = first_index_(i_sym,j_class)
            if ( i_class == j_class ) j_start = i + 1

            ! loop over j indeces

            do j = j_start , last_index_(i_sym,j_class)            
  
              ! update rotation pair count

              npair                 = npair + 1 

              ! update pair count for this symmetry

              trans_%npairpi(i_sym) = trans_%npairpi(i_sym) + 1

              ! update pair counter for this pair

              npair_type(pair_ind)  = npair_type(pair_ind) + 1

            end do ! end j loop

          end do ! end i loop
        
        end do ! end j_class loop 

      end do ! end i_class loop 

    end do ! end i_sym loop

    ! save number of rotation pairs
    rot_pair_%n_tot   = sum(trans_%npairpi)

    ! save type of pair counters
    rot_pair_%n_ad    = npair_type(rot_pair_%act_doc_type)
    rot_pair_%n_ed    = npair_type(rot_pair_%ext_doc_type)
    rot_pair_%n_aa    = npair_type(rot_pair_%act_act_type)
    rot_pair_%n_ea    = npair_type(rot_pair_%ext_act_type)

    return
  end subroutine setup_rotation_indeces

  subroutine allocate_indexing_arrays(nirrep)
    implicit none
    integer, intent(in) :: nirrep
    nirrep_ = nirrep
    call allocate_indexing_array_help(dens_)
    call allocate_indexing_array_help(ints_)
    allocate(first_index_(nirrep_,3))
    allocate(last_index_(nirrep_,3))
    allocate(orb_sym_scr_(nmo_tot_))
    first_index_ = 0
    last_index_  = 0
    orb_sym_scr_ = 0
    return 
  end subroutine allocate_indexing_arrays

  subroutine allocate_indexing_array_help(styp)
    implicit none
    type(sym_info) :: styp
    allocate(styp%ngempi(nirrep_),styp%nnzpi(nirrep_),styp%offset(nirrep_),styp%gemind(nmo_tot_,nmo_tot_))
    styp%ngempi = 0
    styp%nnzpi  = 0
    styp%offset = 0
    styp%gemind = 0
    return
  end subroutine allocate_indexing_array_help

  subroutine deallocate_indexing_arrays()
    implicit none
    if (allocated(first_index_)) deallocate(first_index_)
    if (allocated(last_index_))  deallocate(last_index_)
    if (allocated(orb_sym_scr_)) deallocate(orb_sym_scr_)
    call deallocate_indexing_arrays_help(dens_)
    call deallocate_indexing_arrays_help(ints_)
    return
  end subroutine deallocate_indexing_arrays

  subroutine deallocate_indexing_arrays_help(styp)
    implicit none
    type(sym_info) :: styp
    if ( allocated(styp%ngempi) ) deallocate(styp%ngempi)
    if ( allocated(styp%nnzpi) ) deallocate(styp%nnzpi)
    if ( allocated(styp%offset) ) deallocate(styp%offset)
    if ( allocated(styp%gemind) ) deallocate(styp%gemind)
    return
  end subroutine deallocate_indexing_arrays_help

  subroutine setup_indexing_arrays(ndocpi,nactpi,nextpi)
    implicit none
    integer, intent(in) :: ndocpi(nirrep_),nactpi(nirrep_),nextpi(nirrep_)
    integer :: irrep,oclass,i,j,i_sym,j_sym,ij_sym
    ! ** Figure out index of the first orbital in each class and each irrep
    ! doubly occupied orbitals
    first_index_(1,1)=1
    do irrep=2,nirrep_
      first_index_(irrep,1) = first_index_(irrep-1,1) + ndocpi(irrep-1)
    end do
    ! active orbitals
    first_index_(1,2) = first_index_(nirrep_,1) + ndocpi(nirrep_)
    do irrep=2,nirrep_
      first_index_(irrep,2) = first_index_(irrep-1,2) + nactpi(irrep-1)
    end do
    ! external orbitals
    first_index_(1,3) = first_index_(nirrep_,2) + nactpi(nirrep_)
    do irrep=2,nirrep_
      first_index_(irrep,3) = first_index_(irrep-1,3) + nextpi(irrep-1)
    end do
    ! ** Figure out indexof the last orbital in each class and each irrep
    ! doubly occupied orbitals
    do irrep=1,nirrep_
      last_index_(irrep,1) = first_index_(irrep,1) + ndocpi(irrep)-1
    end do
    ! active orbitals
    do irrep=1,nirrep_
      last_index_(irrep,2) = first_index_(irrep,2) + nactpi(irrep)-1
    end do
    ! active orbitals
    do irrep=1,nirrep_
      last_index_(irrep,3) = first_index_(irrep,3) + nextpi(irrep)-1
    end do
    ! ** save orbital symmetries so that we can set up geminal indices **
    do oclass=1,3
      do irrep=1,nirrep_
        do i=first_index_(irrep,oclass),last_index_(irrep,oclass)
          orb_sym_scr_(i) = irrep
        end do
      end do
    end do
    ! figure out reduced geminal indeces for integral addressing
    do i=1,nmo_tot_
      i_sym = orb_sym_scr_(i)
      do j=1,i
        j_sym = orb_sym_scr_(j)
        ij_sym = group_mult_tab_(i_sym,j_sym)
        ints_%ngempi(ij_sym) = ints_%ngempi(ij_sym) + 1
        ints_%gemind(i,j) = ints_%ngempi(ij_sym)
        ints_%gemind(j,i) = ints_%gemind(i,j)
      end do
    end do
    do irrep=1,nirrep_
      ints_%nnzpi(irrep) = ints_%ngempi(irrep) * ( ints_%ngempi(irrep) + 1 ) /2
    end do
    do irrep=2,nirrep_
      ints_%offset(irrep) = ints_%offset(irrep-1) + ints_%nnzpi(irrep-1)
    end do
    ! figure out reduced geminal indeces for integral addressing
    do i=ndoc_tot_+1,ndoc_tot_+nact_tot_
      i_sym = orb_sym_scr_(i)
      do j=ndoc_tot_+1,i
        j_sym = orb_sym_scr_(j)
        ij_sym = group_mult_tab_(i_sym,j_sym)
        dens_%ngempi(ij_sym) = dens_%ngempi(ij_sym) + 1
        dens_%gemind(i,j) = dens_%ngempi(ij_sym)
        dens_%gemind(j,i) = dens_%gemind(i,j)
      end do
    end do
    do irrep=1,nirrep_
      dens_%nnzpi(irrep) = dens_%ngempi(irrep) * ( dens_%ngempi(irrep) + 1 ) /2
    end do
    do irrep=2,nirrep_
      dens_%offset(irrep) = dens_%offset(irrep-1) + dens_%nnzpi(irrep-1)
    end do
    deallocate(orb_sym_scr_)
    return
  end subroutine setup_indexing_arrays

  subroutine print_info()
    integer :: irrep,i
    ! print the information gathered so far
    write(fid_,'(a5,2x,3(3(a3,1x),5x))')'irrep','d_f','d_l','n_d','a_f','a_l','n_a','e_f','e_l','n_e'
    do irrep=1,nirrep_
      write(fid_,'(i5,2x,3(3(i3,1x),5x))')irrep,(first_index_(irrep,i),last_index_(irrep,i),&
               & last_index_(irrep,i)-first_index_(irrep,i)+1,i=1,3)
    end do
    write(fid_,'(a)')'density information'
    write(fid_,'(a,8(i9,1x))')'ngempi(:)=',dens_%ngempi
    write(fid_,'(a,8(i9,1x))')' nnzpi(:)=',dens_%nnzpi
    write(fid_,'(a,8(i9,1x))')'offset(:)=',dens_%offset
    write(fid_,'(a)')'integral information'
    write(fid_,'(a,8(i9,1x))')'ngempi(:)=',ints_%ngempi
    write(fid_,'(a,8(i9,1x))')' nnzpi(:)=',ints_%nnzpi
    write(fid_,'(a,8(i9,1x))')'offset(:)=',ints_%offset
    write(fid_,'(a)')'rotation pair information:'
    write(fid_,'(4(a,1x,i6,2x),a,1x,i9)')'act-doc pairs:',rot_pair_%n_ad,'ext-doc pairs:',rot_pair_%n_ed,&
                             &'act-act pairs:',rot_pair_%n_aa,'ext-act pairs:',rot_pair_%n_ea,&
                             &'total orbital pairs:',rot_pair_%n_tot
    write(fid_,*)
!
    return
  end subroutine print_info

end module focas_driver
