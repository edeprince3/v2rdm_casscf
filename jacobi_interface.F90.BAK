subroutine jacobi_interface(mo_coeff_out,integrals_1,nnz_i1,integrals_2,nnz_i2,density_1,nnz_d1,&
           &density_2,nnz_d2,syms,ncore_in,nact_in,nvirt_in,nirrep_in,jacobi_data_io, &
           &jacobi_log_file)

  use jacobi_data, only : intfull,denfull,ncore,nact,nvirt,nfrozen,nirrep,aarot,ntot,symdim_d,&
                         & symnnz_d,symdim_i,symnnz_i,symmetries,redind,wp,group_mult_tab, &
                         & angtol,fid,delrot,ntrot,ntsweep,nproc,converged,mo_coeff,print_
  use jacobi_mod

  implicit none
  integer :: nirrep_in,ncore_in,nact_in,nvirt_in,aarot_in,nden,nnz_d1,nnz_d2,nnz_i1,nnz_i2
  real(wp) :: integrals_1(nnz_i1),integrals_2(nnz_i2),density_1(nnz_d1),density_2(nnz_d2)
  real(wp) :: mo_coeff_out(ncore_in+nact_in+nvirt_in,ncore_in+nact_in+nvirt_in)
  real(wp) :: jacobi_data_io(10) 
  character(120) :: jacobi_log_file
  integer  :: syms(ncore_in+nact_in+nvirt_in),i
  logical :: fexist

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
 
  ! copy some variables
  ncore=ncore_in
  nact=nact_in
  nvirt=nvirt_in
  ntot=ncore+nact+nvirt
  nirrep=nirrep_in
  ! maximum number of orbitals in the density arrays
  nden=ncore+nact
  ! aarot=0/1 exclue/include active-active rotations
  nproc=int(jacobi_data_io(1))
  ! aarot=0/1 exclue/include active-active rotations
  aarot=int(jacobi_data_io(2))
  ! nfrozen is the number of core orbitals not optimized
  nfrozen=int(jacobi_data_io(3))
  ! set convergence variables
  angtol=jacobi_data_io(4)
  ! energy tolerance
  mcetol=jacobi_data_io(5)
  ! print flag for jacobi routine
  print_=.false.
  if (jacobi_data_io(6).eq.1) print_=.true.
  if (print_) then
    inquire(file=jacobi_log_file,exist=fexist)
    if (fexist) then
      open(fid,file=jacobi_log_file,status='old',position='append')
    else
      open(fid,file=jacobi_log_file,status='new')
    endif
    ! print some information
    write(fid,'(3(a,1x,i4,4x))')'ncore=',ncore,'nact=',nact,'nvirt=',nvirt
    write(fid,'(a,1x,i3)')'number of core orbitals not optimized=',nfrozen
    write(fid,'(a,1x,i3)')'include active-active orbital rotations=',aarot
    write(fid,'(a,1x,es12.3)')'energy convergence tolerance=',mcetol
    write(fid,'(a,1x,es12.3)')'threshold for including orbital pair=',angtol  
    !call flush(fid)
  endif
  if (size(syms,dim=1).ne.ntot) then
    if (print_) then
      write(fid,'(a)')'size mismatch between symmetry label array and ntot'
      !call flush(fid)
    endif
    stop
  endif
  call setup_symmetry_arrays(syms)
  call allocate_arrays()
  ! save local copy of current MO coefficient array
  mo_coeff=mo_coeff_out
  call copy_arrays(integrals_1,integrals_2,density_1,density_2)
  call jacobi_optimize()
  call update_arrays(integrals_1,integrals_2,density_1,density_2)
  jacobi_data_io(7)=real(ntsweep,kind=wp)
  jacobi_data_io(8)=real(ntrot,kind=wp)
  jacobi_data_io(9)=delrot
  jacobi_data_io(10)=converged
  ! copy local copy into the array to be returned
  mo_coeff_out=mo_coeff
!  mo_coeff_out=matmul(mo_coeff,mo_coeff_out)
!  ! debug print, get rid of this if needed
  call deallocate_arrays()
  close(fid)
  contains

    subroutine update_arrays(integrals_1,integrals_2,density_1,density_2)
      implicit none
      real(wp) :: integrals_1(:),integrals_2(:),density_1(:),density_2(:)
      integer :: nnz1,start,irrep,p,q,p_red,q_red,pq_full,pq_red,p_sym,off(nirrep),pred(ntot),pdim(nirrep)
      ! 2-e integrals
      start=0
      do irrep=1,nirrep
        if (.not.allocated(intfull%e2(irrep)%val)) cycle
        integrals_2(start+1:start+symnnz_i(irrep))=intfull%e2(irrep)%val
        start=start+symnnz_i(irrep)
      end do
      ! figure out reduced index for each orbtal
      pred=0
      pdim=0
      ! number of orbitals in each irrep and for each orbital p save the reduced symmetry index
      do p=1,ntot
        p_sym=symmetries(p)
        pdim(p_sym)=pdim(p_sym)+1
        pred(p)=pdim(p_sym)
      end do
      ! total number of nonzeros
      nnz1=0
      do irrep=1,nirrep
        nnz1=nnz1+pdim(irrep)*(pdim(irrep)+1)/2
      end do
      ! offsets
      off(1)=0
      do irrep=2,nirrep
        off(irrep)=off(irrep-1)+pdim(irrep-1)*(pdim(irrep-1)+1)/2
      end do
      ! some checks
      if (nnz1.ne.size(intfull%e1,dim=1)) then
        if (print_) then
          write(fid,'(a)')'nnz1 based on number of orbitals in each irrep does not agree with size of 1-e integral array'
          !call flush(fid)
        endif
        stop
      end if
      if (nnz1.ne.symdim_i(1)) then
        if (print_) then
          write(fid,'(a)')'nnz1 based on number of orbitals in each irrep does not agree with nnz stored symdim_i(1)'
          call flush(fid)
        endif
        stop
      endif
      ! now copy the integrals
      pq_full=0
      do p=1,ntot
        p_sym=symmetries(p)
        p_red=pred(p)
        do q=1,p
          if (symmetries(q).ne.p_sym) cycle
          pq_full=pq_full+1
          q_red=pred(q)
          pq_red=off(p_sym)+p_red*(p_red-1)/2+q_red
          integrals_1(pq_red)=intfull%e1(pq_full)
        end do
      end do
      start=0
      do irrep=1,nirrep
        if (.not.allocated(denfull%e2(irrep)%val)) cycle
        density_2(start+1:start+symnnz_d(irrep))=denfull%e2(irrep)%val
        start=start+symnnz_d(irrep)
      end do
      ! figure out reduced index for each orbtal
      pred=0
      pdim=0
      ! number of orbitals in each irrep and for each orbital p save the reduced symmetry index
      do p=1,nden
        p_sym=symmetries(p)
        pdim(p_sym)=pdim(p_sym)+1
        pred(p)=pdim(p_sym)
      end do
      ! total number of nonzeros
      nnz1=0
      do irrep=1,nirrep
        nnz1=nnz1+pdim(irrep)*(pdim(irrep)+1)/2
      end do
      ! offsets
      off(1)=0
      do irrep=2,nirrep
        off(irrep)=off(irrep-1)+pdim(irrep-1)*(pdim(irrep-1)+1)/2
      end do
      ! some checks
      if (nnz1.ne.size(denfull%e1,dim=1)) then
        if (print_) then
          write(fid,'(a,1x,i6)')'nnz1 based on number of orbitals in each irrep does not agree with size of 1-e density array',&
            & nnz1
          !call flush(fid)
        endif
        stop
      end if
      if (nnz1.ne.symdim_d(1)) then
        if (print_) then
          write(fid,'(a,1x,i6)')'nnz1 based on number of orbitals in each irrep does not agree with nnz stored symdim_d(1)',&
           & symdim_d(1)
          !call flush(fid)
        endif
        stop
      endif
      ! now copy the integrals
      pq_full=0
      do p=1,nden
        p_sym=symmetries(p)
        p_red=pred(p)
        do q=1,p
          if (symmetries(q).ne.p_sym) cycle
          pq_full=pq_full+1
          q_red=pred(q)
          pq_red=off(p_sym)+p_red*(p_red-1)/2+q_red
          density_1(pq_red)=denfull%e1(pq_full)
        end do
      end do
      return
    end subroutine update_arrays

    subroutine copy_arrays(integrals_1,integrals_2,density_1,density_2)
      implicit none
      real(wp) :: integrals_1(:),integrals_2(:),density_1(:),density_2(:)
      integer :: nnz1,start,irrep,p,q,p_red,q_red,pq_full,pq_red,p_sym,off(nirrep),pred(ntot),pdim(nirrep)
      if (sum(symnnz_i).ne.size(integrals_2,dim=1)) then 
        if (print_) then
          write(fid,'(a)')'total number of nnz in integrals array does not match total nnz elements allocated in intfull%e2'
          !call flush(fid)
        endif
        stop
      else
        start=0
        do irrep=1,nirrep
          if (.not.allocated(intfull%e2(irrep)%val)) cycle
          intfull%e2(irrep)%val=integrals_2(start+1:start+symnnz_i(irrep))
          start=start+symnnz_i(irrep)
        end do
        ! figure out reduced index for each orbtal
        pred=0
        pdim=0
        ! number of orbitals in each irrep and for each orbital p save the reduced symmetry index
        do p=1,ntot
          p_sym=symmetries(p)
          pdim(p_sym)=pdim(p_sym)+1
          pred(p)=pdim(p_sym)
        end do
        nnz1=0
        do irrep=1,nirrep
          nnz1=nnz1+pdim(irrep)*(pdim(irrep)+1)/2
        end do
        off(1)=0
        do irrep=2,nirrep
          off(irrep)=off(irrep-1)+pdim(irrep-1)*(pdim(irrep-1)+1)/2
        end do
        if (nnz1.ne.size(intfull%e1,dim=1)) then
          if (print_) then
            write(fid,'(a)')'nnz1 based on number of orbitals in each irrep does not agree with size of 1-e integral array'
            !call flush(fid)
          endif
          stop
        end if
        if (nnz1.ne.symdim_i(1)) then
          if (print_) then
            write(fid,'(a)')'nnz1 based on number of orbitals in each irrep does not agree with nnz stored symdim_i(1)'
            !call flush(fid)
          endif
          stop
        endif
        ! now copy the integrals
        pq_full=0
        do p=1,ntot
          p_sym=symmetries(p)
          p_red=pred(p)
          do q=1,p
            if (symmetries(q).ne.p_sym) cycle
            pq_full=pq_full+1
            q_red=pred(q)
            pq_red=off(p_sym)+p_red*(p_red-1)/2+q_red
            intfull%e1(pq_full)=integrals_1(pq_red)           
          end do
        end do
      endif
      if (sum(symnnz_d).ne.size(density_2,dim=1)) then
        if (print_) then
          write(fid,'(a)')'total number of nnz in density array does not match total nnz elements allocated in denfull%e2'
          !call flush(fid)
        endif
        stop
      else
        start=0
        do irrep=1,nirrep
          if (.not.allocated(denfull%e2(irrep)%val)) cycle
          denfull%e2(irrep)%val=density_2(start+1:start+symnnz_d(irrep))
          start=start+symnnz_d(irrep)
        end do
        ! figure out reduced index for each orbtal
        pred=0
        pdim=0
        ! number of orbitals in each irrep and for each orbital p save the reduced symmetry index
        do p=1,nden
          p_sym=symmetries(p)
          pdim(p_sym)=pdim(p_sym)+1
          pred(p)=pdim(p_sym)
        end do
        nnz1=0
        do irrep=1,nirrep
          nnz1=nnz1+pdim(irrep)*(pdim(irrep)+1)/2
        end do
        off(1)=0
        do irrep=2,nirrep
          off(irrep)=off(irrep-1)+pdim(irrep-1)*(pdim(irrep-1)+1)/2
        end do
        if (nnz1.ne.size(denfull%e1,dim=1)) then
          if (print_) then
            write(fid,'(a,1x,i6)')'nnz1 based on number of orbitals in each irrep does not agree with size of 1-e density array', &
              & nnz1
            !call flush(fid)
          endif
          stop
        end if
        if (nnz1.ne.symdim_d(1)) then
          if (print_) then
            write(fid,'(a,1x,i6)')'nnz1 based on number of orbitals in each irrep does not agree with nnz stored symdim_d(1)', &
              & symdim_d(1)
            !call flush(fid)
          endif
          stop
        endif
        ! now copy the integrals
        pq_full=0
        do p=1,nden
          p_sym=symmetries(p)
          p_red=pred(p)
          do q=1,p
            if (symmetries(q).ne.p_sym) cycle
            pq_full=pq_full+1
            q_red=pred(q)
            pq_red=off(p_sym)+p_red*(p_red-1)/2+q_red
            denfull%e1(pq_full)=density_1(pq_red)
          end do
        end do
      endif
      return
    end subroutine copy_arrays

    subroutine deallocate_arrays()
      implicit none
      integer :: irrep
      if (allocated(mo_coeff)) deallocate(mo_coeff)
      if (allocated(intfull%e1)) deallocate(intfull%e1)
      if (allocated(intfull%e2)) then
        do irrep=1,nirrep
          if (allocated(intfull%e2(irrep)%val)) deallocate(intfull%e2(irrep)%val)
        end do
        deallocate(intfull%e2)
      endif
      if (allocated(denfull%e1)) deallocate(denfull%e1)
      if (allocated(denfull%e2)) then
        do irrep=1,nirrep
          if (allocated(denfull%e2(irrep)%val)) deallocate(denfull%e2(irrep)%val)
        end do
        deallocate(denfull%e2)
      endif
      deallocate(redind,symdim_d,symnnz_d,symdim_i,symnnz_i,symmetries)
      return
    end subroutine deallocate_arrays  

    subroutine allocate_arrays()
      implicit none
      integer :: irrep,i
      allocate(mo_coeff(ntot,ntot))
      do irrep=1,nirrep
        symnnz_d(irrep)=symdim_d(irrep)*(symdim_d(irrep)+1)/2
        symnnz_i(irrep)=symdim_i(irrep)*(symdim_i(irrep)+1)/2
      end do
      allocate(intfull%e1(symdim_i(1)))
      intfull%e1=0.0_wp
      allocate(intfull%e2(nirrep))
      do irrep=1,nirrep
        allocate(intfull%e2(irrep)%val(symnnz_i(irrep)))
        intfull%e2(irrep)%val=0.0_wp 
      end do
      allocate(denfull%e1(symdim_d(1)))
      denfull%e1=0.0_wp
      allocate(denfull%e2(nirrep))
      do irrep=1,nirrep
        allocate(denfull%e2(irrep)%val(symnnz_d(irrep)))
        denfull%e2(irrep)%val=0.0_wp
      end do
      return
    end subroutine allocate_arrays

    subroutine setup_symmetry_arrays(syms)
      implicit none
      integer :: syms(:)
      integer :: p,q,p_sym,pq_sym
      allocate(redind(ntot,ntot),symdim_i(nirrep),symnnz_i(nirrep),symnnz_d(nirrep),symdim_d(nirrep),symmetries(ntot))
      symmetries=syms
      symdim_i=0
      symnnz_i=0
      symdim_d=0
      symnnz_d=0
      redind=0
      do p=1,ntot
        p_sym=symmetries(p)
        do q=1,p
          pq_sym=group_mult_tab(p_sym,symmetries(q))
          symdim_i(pq_sym)=symdim_i(pq_sym)+1
          redind(p,q)=symdim_i(pq_sym)
          redind(q,p)=redind(p,q)
        end do
      end do
      do p=1,nden
        p_sym=symmetries(p)
        do q=1,p
          pq_sym=group_mult_tab(p_sym,symmetries(q))
          symdim_d(pq_sym)=symdim_d(pq_sym)+1
        end do
      end do
      return
    end subroutine setup_symmetry_arrays

end subroutine jacobi_interface
