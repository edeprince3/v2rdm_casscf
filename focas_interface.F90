subroutine focas_interface(mo_coeff_out,integrals_1,nnz_i1,integrals_2,nnz_i2,density_1,nnz_d1,&
           &density_2,nnz_d2,syms,ncore_in,nact_in,nvirt_in,nirrep_in,jacobi_data_io, &
           &jacobi_log_file)

  use focas_driver, only : focas_optimize

  implicit none
  integer, parameter :: fid=99
  integer, parameter :: wp = selected_real_kind(10)
  integer, parameter :: mult_tab(8,8) = reshape( (/ &  ! irrep multiplication table
      & 1,2,3,4,5,6,7,8, &
      & 2,1,4,3,6,5,8,7, &
      & 3,4,1,2,7,8,5,6, &
      & 4,3,2,1,8,7,6,5, &
      & 5,6,7,8,1,2,3,4, &
      & 6,5,8,7,2,1,4,3, &
      & 7,8,5,6,3,4,1,2, &
      & 8,7,6,5,4,3,2,1  /), (/8,8/) )

  integer :: nirrep_in,ncore_in,nact_in,nvirt_in,aarot_in,nnz_d1,nnz_d2,nnz_i1,nnz_i2
  integer :: nproc,aarot,nfrozen,angtol,detol,mcetol,print
  real(wp) :: integrals_1(nnz_i1),integrals_2(nnz_i2),density_1(nnz_d1),density_2(nnz_d2)
  real(wp) :: mo_coeff_out(ncore_in+nact_in+nvirt_in,ncore_in+nact_in+nvirt_in)
  real(wp) :: jacobi_data_io(11) 
  character(120) :: jacobi_log_file
  integer  :: syms(ncore_in+nact_in+nvirt_in)
  logical :: fexist

  real(wp) :: delrot

  integer :: nactpi(nirrep_in),ndocpi(nirrep_in),nextpi(nirrep_in)
  integer :: ndoc,nact,next,nmo,nirrep,converged
  integer :: nnz_int1,nnz_int2,nnz_den1,nnz_den2,df_ints
  integer :: gemind_int(ncore_in+nact_in+nvirt_in,ncore_in+nact_in+nvirt_in)
  integer :: gemind_den_new(ncore_in+nact_in+nvirt_in,ncore_in+nact_in+nvirt_in)
  integer :: gemind_int_new(ncore_in+nact_in+nvirt_in,ncore_in+nact_in+nvirt_in)
  integer :: energy_to_class_map(ncore_in+nact_in+nvirt_in)
  integer :: class_to_energy_map(ncore_in+nact_in+nvirt_in)
  integer :: irrep_to_class_map(ncore_in+nact_in+nvirt_in)
  integer :: class_to_irrep_map(ncore_in+nact_in+nvirt_in)
  integer :: first_index(nirrep_in,3)
  integer :: last_index(nirrep_in,3)
  integer :: nnz_den_psi4(nirrep_in)
  integer :: nnz_den_new(nirrep_in)
  integer :: nnz_int(nirrep_in)
  integer :: offset_den_psi4(nirrep_in)
  integer :: offset_den_new(nirrep_in)
  integer :: offset_int(nirrep_in)
  integer :: offset_irrep_int1(nirrep_in)
  integer :: offset_irrep_den1(nirrep_in)
  real(wp) :: dele_tol,gnorm_tol,gnorm
 

!  jacobi_data_io:
!  1) nproc
!  2) aarot
!  3) nfrozen
!  4) angtol
!  5) mcetol
!  6) jacobi_print
!  7) ntsweep
!  8) ntrot
!  9) delrot
! 10) converged (1=yes/0=no)
! 11) df integral flag 
  ! copy some variables
 
  ndoc=ncore_in 
  nact=nact_in
  next=nvirt_in
  nmo=ncore_in+nact_in+nvirt_in
  nirrep=nirrep_in
  ! aarot=0/1 exclue/include active-active rotations
  nproc=int(jacobi_data_io(1))
  ! aarot=0/1 exclue/include active-active rotations
  aarot=int(jacobi_data_io(2))
  ! nfrozen is the number of core orbitals not optimized
  nfrozen=int(jacobi_data_io(3))
  ! set convergence variables
  gnorm_tol=jacobi_data_io(4)
  ! energy tolerance
  dele_tol=jacobi_data_io(5)
  ! set density-fitted integral flag
  df_ints = int(jacobi_data_io(11))
  ! print flag for jacobi routine
  print=int(jacobi_data_io(6))
  if (print == 1) then
    inquire(file=jacobi_log_file,exist=fexist)
    if (fexist) then
      open(fid,file=jacobi_log_file,status='old',position='append')
    else
      open(fid,file=jacobi_log_file,status='new')
    endif
    ! print some information
    write(fid,'(3(a,1x,i4,4x))')'ncore=',ndoc,'nact=',nact,'nvirt=',next
    write(fid,'(a,1x,i3)')'number of core orbitals not optimized=',nfrozen
    write(fid,'(a,1x,i3)')'include active-active orbital rotations=',aarot
    write(fid,'(a,1x,es12.3)')'energy convergence tolerance=',mcetol
    write(fid,'(a,1x,es12.3)')'threshold for including orbital pair=',angtol  
    !call flush(fid)
  endif
  if (size(syms,dim=1).ne.nmo) then
    if (print==1) then
      write(fid,'(a)')'size mismatch between symmetry label array and nmo'
      !call flush(fid)
    endif
    stop
  endif
  call setup_symmetry_arrays(syms)
  call initial_sort()
  nnz_int1 = gemind_int_new(last_index(nirrep,3),last_index(nirrep,3))
  nnz_den1 = gemind_den_new(last_index(nirrep,2),last_index(nirrep,2))
  if (df_ints == 0 ) then
    nnz_int2 = sum(nnz_int)  
  else
    nnz_int2 = size(integrals_2,dim=1)
  endif
  nnz_den2 = sum(nnz_den_new)
  call focas_optimize(integrals_1,nnz_int1,integrals_2,nnz_int2,                    &
                    & density_1(1:nnz_den1),nnz_den1,density_2(1:nnz_den2),nnz_den2,&
                    & ndocpi,nactpi,nextpi,nirrep,gnorm_tol,dele_tol,gnorm,delrot,converged,df_ints)
  call final_sort()
  jacobi_data_io(7)=real(1.0_wp,kind=wp)
  jacobi_data_io(8)=real(1.0_wp,kind=wp)
  jacobi_data_io(9)=delrot
  jacobi_data_io(10)=converged
  ! copy local copy into the array to be returned
!  mo_coeff_out=mo_coeff
!  mo_coeff_out=matmul(mo_coeff,mo_coeff_out)
  close(fid)
  contains

    subroutine final_sort()
      implicit none

      integer :: p_sym,q_sym,r_sym,s_sym,pq_sym
      integer :: pq_c,rs_c,pq_i,rs_i,pqrs_c,pqrs_i,pp_c,qq_c
      integer :: p_class,q_class,r_class,s_class,q_max,s_max
      integer :: p_c,q_c,r_c,s_c
      integer :: p_i,q_i,r_i,s_i
      integer :: pq_off,max_dim
      real(wp), allocatable :: block(:)

      max_dim = max(size(integrals_1,dim=1),maxval(nnz_int))
      allocate(block(max_dim))

      ! copy 1-e integrals

      block=huge(1.0_wp)
      do p_sym = 1 , nirrep
        pq_off = offset_irrep_int1(p_sym)
        do p_class = 1 , 3
          do q_class = 1 , p_class
            do p_c = first_index(p_sym,p_class) , last_index(p_sym,p_class)
              q_max=last_index(p_sym,q_class)
              p_i = class_to_irrep_map(p_c)
              if ( p_class == q_class ) q_max=p_c
              do q_c = first_index(p_sym,q_class),q_max
                q_i = class_to_irrep_map(q_c)
                pq_c = gemind_int_new(p_c,q_c)
                pq_i = pq_ind(p_i,q_i)+pq_off
                block(pq_i)=integrals_1(pq_c)
              end do
            end do
          end do
        end do
      end do
      integrals_1 = block(1:size(integrals_1,dim=1))

      ! copy 2-e integrals

      if ( df_ints == 0 ) then
        do pq_sym = 1 , nirrep
          block=huge(1.0_wp)
          pq_off=offset_int(pq_sym)
          do p_class = 1 , 3
            do q_class = 1 , p_class
              do p_sym = 1 , nirrep
                q_sym = mult_tab(pq_sym,p_sym)
                if ( ( q_sym > p_sym ) .and. ( p_class == q_class ) ) cycle
                do r_class = 1 , 3
                  do s_class = 1 , r_class
                    do r_sym = 1 , nirrep
                      s_sym = mult_tab(pq_sym,r_sym)
                      if ( ( s_sym > r_sym ) .and. ( r_class==q_class ) ) cycle
                      do p_c = first_index(p_sym,p_class),last_index(p_sym,p_class)
                        p_i = class_to_energy_map(p_c)
                        q_max = last_index(q_sym,q_class)
                        if ( ( p_class == q_class ) .and. ( p_sym == q_sym ) ) q_max = p_c
                        do q_c = first_index(q_sym,q_class),q_max
                          q_i = class_to_energy_map(q_c)
                          pq_c = gemind_int_new(p_c,q_c)
                          pq_i = gemind_int(p_i,q_i)
                          do r_c=first_index(r_sym,r_class),last_index(r_sym,r_class)
                            r_i = class_to_energy_map(r_c)
                            s_max = last_index(s_sym,s_class)
                            if ( ( s_sym == r_sym ) .and. ( s_class == r_class ) ) s_max = r_c
                            do s_c = first_index(s_sym,s_class),s_max
                              s_i = class_to_energy_map(s_c)
                              rs_c = gemind_int_new(r_c,s_c)
                              rs_i = gemind_int(r_i,s_i)
                              pqrs_c = pq_ind(pq_c,rs_c)
                              pqrs_i = pq_ind(pq_i,rs_i)
                              block(pqrs_i) = integrals_2(pqrs_c+pq_off)
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
          integrals_2(pq_off+1:pq_off+nnz_int(pq_sym))=block(1:nnz_int(pq_sym))
        end do
      endif

      deallocate(block)

      return

    end subroutine final_sort

    subroutine initial_sort()

      implicit none

      integer :: p_sym,q_sym,r_sym,s_sym,pq_sym
      integer :: pq_c,rs_c,pq_i,rs_i,pqrs_c,pqrs_i,pp_c,qq_c
      integer :: p_class,q_class,r_class,s_class,q_max,s_max
      integer :: p_c,q_c,r_c,s_c
      integer :: p_i,q_i,r_i,s_i
      integer :: pq_off,max_dim
      real(wp), allocatable :: block(:)
      real(wp) :: tr_d2

      max_dim = max(size(integrals_1,dim=1),maxval(nnz_int))
      allocate(block(max_dim))

      ! copy 1-e integrals

      block=huge(1.0_wp)
      do p_sym = 1 , nirrep
        pq_off = offset_irrep_int1(p_sym)
        do p_class = 1 , 3
          do q_class = 1 , p_class
            do p_c = first_index(p_sym,p_class) , last_index(p_sym,p_class)
              q_max=last_index(p_sym,q_class)
              p_i = class_to_irrep_map(p_c)
              if ( p_class == q_class ) q_max=p_c
              do q_c = first_index(p_sym,q_class),q_max
                q_i = class_to_irrep_map(q_c)
                pq_c = gemind_int_new(p_c,q_c)
                pq_i = pq_ind(p_i,q_i)+pq_off
                block(pq_c)=integrals_1(pq_i)
              end do
            end do
          end do 
        end do
      end do
      integrals_1 = block(1:size(integrals_1,dim=1))

      ! copy 1-e density
      tr_d2=0.0_wp
      block=huge(1.0_wp)
      do p_sym = 1 , nirrep
        pq_off = offset_irrep_den1(p_sym)
        do p_c = first_index(p_sym,2) , last_index(p_sym,2)
          p_i = class_to_irrep_map(p_c)
          do q_c = first_index(p_sym,2),p_c-1
            q_i = class_to_irrep_map(q_c)
            pq_c = gemind_den_new(p_c,q_c)
            pq_i = pq_ind(p_i,q_i)+pq_off
            block(pq_c)=0.5_wp*density_1(pq_i)
          end do
          pq_c = gemind_den_new(p_c,p_c)
          pq_i = pq_ind(p_i,p_i)+pq_off
          block(pq_c)=density_1(pq_i)
          tr_d2=tr_d2+block(pq_c)
        end do
      end do
      p_sym = gemind_den_new(last_index(nirrep,2),last_index(nirrep,2))
      density_1(1:p_sym) = block(1:p_sym)

      ! copy 2-e integrals
      if ( df_ints == 0 ) then
        do pq_sym = 1 , nirrep
          block=huge(1.0_wp)
          pq_off=offset_int(pq_sym)
          do p_class = 1 , 3
            do q_class = 1 , p_class
              do p_sym = 1 , nirrep
                q_sym = mult_tab(pq_sym,p_sym)
                if ( ( q_sym > p_sym ) .and. ( p_class == q_class ) ) cycle
                do r_class = 1 , 3
                  do s_class = 1 , r_class
                    do r_sym = 1 , nirrep
                      s_sym = mult_tab(pq_sym,r_sym)
                      if ( ( s_sym > r_sym ) .and. ( r_class==q_class ) ) cycle
                      do p_c = first_index(p_sym,p_class),last_index(p_sym,p_class)
                        p_i = class_to_energy_map(p_c)
                        q_max = last_index(q_sym,q_class)
                        if ( ( p_class == q_class ) .and. ( p_sym == q_sym ) ) q_max = p_c
                        do q_c = first_index(q_sym,q_class),q_max
                          q_i = class_to_energy_map(q_c)
                          pq_c = gemind_int_new(p_c,q_c)
                          pq_i = gemind_int(p_i,q_i)
                          do r_c=first_index(r_sym,r_class),last_index(r_sym,r_class)
                            r_i = class_to_energy_map(r_c)
                            s_max = last_index(s_sym,s_class)
                            if ( ( s_sym == r_sym ) .and. ( s_class == r_class ) ) s_max = r_c                        
                            do s_c = first_index(s_sym,s_class),s_max
                              s_i = class_to_energy_map(s_c)
                              rs_c = gemind_int_new(r_c,s_c)
                              rs_i = gemind_int(r_i,s_i)
                              pqrs_c = pq_ind(pq_c,rs_c)
                              pqrs_i = pq_ind(pq_i,rs_i)
                            block(pqrs_c) = integrals_2(pqrs_i+pq_off)
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
          integrals_2(pq_off+1:pq_off+nnz_int(pq_sym))=block(1:nnz_int(pq_sym))
        end do
      end if

      ! copy/scale 2-e active density
      do pq_sym = 1 , nirrep
        block=huge(1.0_wp)
        pq_off=offset_den_psi4(pq_sym)
        do p_sym = 1 , nirrep
          q_sym = mult_tab(pq_sym,p_sym)
          if ( q_sym > p_sym ) cycle
          do r_sym = 1 , nirrep
            s_sym = mult_tab(pq_sym,r_sym)
            if ( s_sym > r_sym ) cycle
            do p_c = first_index(p_sym,2),last_index(p_sym,2)
              p_i = class_to_energy_map(p_c)
              q_max = last_index(q_sym,2)
              if ( p_sym == q_sym ) q_max = p_c
              do q_c = first_index(q_sym,2),q_max
                q_i = class_to_energy_map(q_c)
                pq_c = gemind_den_new(p_c,q_c)
                pq_i = gemind_int(p_i,q_i)
                do r_c=first_index(r_sym,2),last_index(r_sym,2)
                  r_i = class_to_energy_map(r_c)
                  s_max = last_index(s_sym,2)
                  if ( s_sym == r_sym ) s_max = r_c
                  do s_c = first_index(s_sym,2),s_max
                    s_i = class_to_energy_map(s_c)
                    rs_c = gemind_den_new(r_c,s_c)
                    rs_i = gemind_int(r_i,s_i)
                    pqrs_c = pq_ind(pq_c,rs_c)
                    pqrs_i = pq_ind(pq_i,rs_i)
                    block(pqrs_c) = den_fac(p_c,q_c,r_c,s_c) * density_2(pqrs_i+pq_off)
                  end do
                end do
              end do
            end do
          end do
        end do
        pq_off = offset_den_new(pq_sym)
        density_2(pq_off+1:pq_off+nnz_den_new(pq_sym)) = block(1:nnz_den_new(pq_sym))
      end do

      density_2(1:sum(nnz_den_new)) = 2.0_wp * density_2(1:sum(nnz_den_new))

      ! check trace
      tr_d2=0.0_wp
      do p_sym = 1, nirrep
        do q_sym = 1,nirrep
          do p_c=first_index(p_sym,2),last_index(p_sym,2)
            do q_c=first_index(q_sym,2),last_index(q_sym,2)
              pp_c=gemind_den_new(p_c,p_c)
              qq_c=gemind_den_new(q_c,q_c)
              tr_d2=tr_d2+density_2(pq_ind(pp_c,qq_c))
            end do
          end do
        end do
      end do

      deallocate(block)

      return

    end subroutine initial_sort

    subroutine setup_symmetry_arrays(syms)
      implicit none
      integer :: syms(:)
      integer :: p,q,p_class,q_class,p_sym,pq_sym
      integer :: p_c,q_c,p_i,q_i
      integer :: dims(nirrep),sym_class(nmo)

      ! determine the number of mos per irrep for each class

      nactpi = 0
      ndocpi = 0
      nextpi = 0

      do p = 1 , ndoc
        ndocpi(syms(p)) = ndocpi(syms(p)) + 1
      end do
      do p = ndoc + 1 , ndoc + nact
        nactpi(syms(p)) = nactpi(syms(p)) + 1
      end do
      do p = ndoc + nact + 1 , nmo
        nextpi(syms(p)) = nextpi(syms(p)) + 1
      end do
      ! addressing for the density/integrals is the same for psi4
      ! the offset/nnz arrays for new order are also the same

      dims=0
      do p = 1 , nmo
        p_sym = syms(p)
        do q = 1 , p
          pq_sym = mult_tab(syms(q),p_sym)
          dims(pq_sym)=dims(pq_sym)+1
          gemind_int(p,q)=dims(pq_sym)
          gemind_int(q,p)=dims(pq_sym)
        end do
      end do

      nnz_int = 0
      do pq_sym = 1 , nirrep
        nnz_int(pq_sym) = dims(pq_sym)*(dims(pq_sym)+1)/2
      end do
      offset_int = 0
      do pq_sym = 2 , nirrep
        offset_int(pq_sym) = nnz_int(pq_sym-1) + offset_int(pq_sym-1)
      end do

      ! determine dimensions for density arrays

      dims=0
      do p = 1 , ndoc+nact
        p_sym = syms(p)
        do q = 1 , p
          pq_sym = mult_tab(syms(q),p_sym)
          dims(pq_sym)=dims(pq_sym)+1
        end do
      end do
      nnz_den_psi4 = 0
      do pq_sym = 1 , nirrep
        nnz_den_psi4(pq_sym) = dims(pq_sym)*(dims(pq_sym)+1)/2
      end do
      offset_den_psi4 = 0
      do pq_sym = 2 , nirrep
        offset_den_psi4(pq_sym) = nnz_den_psi4(pq_sym-1) + offset_den_psi4(pq_sym-1)
      end do

      ! determine irrep->class & class->irrep map arrays

      p_class=0
      do p_sym = 1 , nirrep
        do p = 1 , ndoc
          if ( syms(p) /= p_sym ) cycle
          p_class = p_class + 1
          energy_to_class_map(p) = p_class
          class_to_energy_map(p_class) = p
        end do    
      end do
      do p_sym = 1 , nirrep
        do p = ndoc + 1 , ndoc + nact
          if ( syms(p) /= p_sym ) cycle
          p_class = p_class + 1
          energy_to_class_map(p) = p_class
          class_to_energy_map(p_class) = p
        end do
      end do
      do p_sym = 1 , nirrep
        do p = ndoc + nact + 1 , nmo
          if ( syms(p) /= p_sym ) cycle
          p_class = p_class + 1
          energy_to_class_map(p) = p_class
          class_to_energy_map(p_class) = p
        end do
      end do

      ! set up orbital symmetry arrays in new order

      do p = 1 , nmo
        p_c = energy_to_class_map(p)
        sym_class(p_c) = syms(p)
      end do

      ! set up geminal addressing arrays for the integrals

      dims=0
      do p = 1 , nmo
        p_sym = sym_class(p)
        do q = 1 , p
          pq_sym = mult_tab(sym_class(q),p_sym)
          dims(pq_sym)=dims(pq_sym)+1
          gemind_int_new(p,q)=dims(pq_sym)
          gemind_int_new(q,p)=dims(pq_sym)
        end do
      end do

      ! set up geminal addressing arrays for the densities

      dims=0
      do p = ndoc + 1 , ndoc + nact
        p_sym = sym_class(p)
        do q = ndoc + 1 , p
          pq_sym = mult_tab(sym_class(q),p_sym)
          dims(pq_sym)=dims(pq_sym)+1
          gemind_den_new(p,q)=dims(pq_sym)
          gemind_den_new(q,p)=dims(pq_sym)
        end do
      end do

      nnz_den_new = 0
      do pq_sym = 1 , nirrep
        nnz_den_new(pq_sym) = dims(pq_sym)*(dims(pq_sym)+1)/2
      end do
      offset_den_new = 0
      do pq_sym = 2 , nirrep
        offset_den_new(pq_sym) = nnz_den_new(pq_sym-1) + offset_den_new(pq_sym-1)
      end do

      first_index = 0
      last_index = 0

      first_index(1,1)=1
      last_index(1,1)=first_index(1,1) + ndocpi(1) - 1
      do p_sym = 2,nirrep
        first_index(p_sym,1) = first_index(p_sym-1,1) + ndocpi(p_sym-1)
        last_index(p_sym,1) = first_index(p_sym,1) + ndocpi(p_sym) - 1
      end do

      first_index(1,2)=first_index(nirrep,1)+ndocpi(nirrep)
      last_index(1,2)=first_index(1,2) + nactpi(1) - 1
      do p_sym = 2,nirrep
        first_index(p_sym,2) = first_index(p_sym-1,2) + nactpi(p_sym-1)
        last_index(p_sym,2) = first_index(p_sym,2) + nactpi(p_sym) - 1
      end do

      first_index(1,3)=first_index(nirrep,2)+nactpi(nirrep)
      last_index(1,3)=first_index(1,3) + nextpi(1) - 1
      do p_sym = 2,nirrep
        first_index(p_sym,3) = first_index(p_sym-1,3) + nextpi(p_sym-1)
        last_index(p_sym,3) = first_index(p_sym,3) + nextpi(p_sym) - 1
      end do

      dims=0
      do p_sym = 1 , nirrep
        dims(p_sym) = 0
        do p_class = 1 , 3
          do p = first_index(p_sym,p_class),last_index(p_sym,p_class)
            dims(p_sym) = dims(p_sym) + 1
            class_to_irrep_map(p) = dims(p_sym)
          end do
        end do
      end do
      offset_irrep_int1 = 0
      do p_sym = 2 , nirrep
        offset_irrep_int1(p_sym) = offset_irrep_int1(p_sym-1) + dims(p_sym-1)*(dims(p_sym-1)+1)/2
      end do

      dims=0
      do p_sym = 1 , nirrep
        do p_class = 1 , 2
          do p = first_index(p_sym,p_class),last_index(p_sym,p_class)
            dims(p_sym) = dims(p_sym) + 1
          end do
        end do
      end do
      offset_irrep_den1 = 0
      do p_sym = 2 , nirrep
        offset_irrep_den1(p_sym) = offset_irrep_den1(p_sym-1) + dims(p_sym-1)*(dims(p_sym-1)+1)/2
      end do
      return
    end subroutine setup_symmetry_arrays

    pure function pq_ind(i,j)
! this function computes the two-electron index index (lower triangular reference)
! index = ii*(ii-1)/2+jj where ii=max(i,j) and jj=min(i,j)
! the ishft(k,-1) divides the value of the integer k by 2 and seems to be somewhat
! faster than the regular human-readable expression
      implicit none
      integer, intent(in) ::i,j
      integer :: pq_ind
      if (i.ge.j) then
        pq_ind=ishft(i*(i-1),-1)+j
        return
      else
        pq_ind=ishft(j*(j-1),-1)+i
        return
      end if
    end function pq_ind

    function den_fac(i,j,k,l)
      integer, intent(in) :: i,j,k,l
      real(wp) :: den_fac
      integer :: fac
      fac=1
      if ( i /= j ) fac = 2
      if ( k /= l ) fac = 2 * fac
      if (pq_ind(i,j) /= pq_ind(k,l)) fac = 2 * fac
      den_fac = 1.0_wp / real(fac,kind=wp)
    end function den_fac

end subroutine focas_interface
