module jacobi_mod
  use jacobi_data, only : intfull,denfull,redind,rotinds,ncore,nact,nvirt,ntot,npair,maxind_nnz, &
                        & symmetries,wp,nirrep,nproc,group_mult_tab,aarot,nfrozen,mcetol,angtol,delrot,&
                        & converged,ntrot,ntsweep,fid,mo_coeff,print_
  use jacobi_maxind_mod, only : jacobi_maxind
  implicit none

  contains

    subroutine jacobi_optimize()
      implicit none
      integer :: isweep,maxsweep,nrot,i,j
      integer, allocatable :: maxind(:,:)
      real(wp) :: xpi(65),ypi(65)
      real(wp) :: et_last,et_first,dele
      real(wp) :: e1_first,e2_first,e1_last,e2_last
! debug variables
      call jacobi_rotinds()
      if (print_) then
        write(fid,'(a)')'orbital symmetries'
        write(fid,'(30(i3))')(i,i=1,ntot)
        write(fid,'(30(i3))')symmetries
        write(fid,'(a,3x,i4)')'rotation pairs',npair
        write(fid,'(10(2(i3,1x),6x))')(rotinds(i+1:i+2),i=0,2*npair-2,2)
      endif
      !call flush(fid)
      call jacobi_maxind()
      allocate(maxind(maxind_nnz,nirrep))
      call jacobi_getpi(xpi,ypi)
      e1_first=0.0_wp
      e2_first=0.0_wp
      do i=1,nirrep
        j=0
        if (allocated(denfull%e2(i)%val)) j=size(denfull%e2(i)%val,dim=1)
        if (j.eq.0) cycle
        
        e2_first=e2_first+dot_product(denfull%e2(i)%val(1:j),intfull%e2(i)%val(1:j))
      end do
      j=size(denfull%e1)
      e1_first=dot_product(denfull%e1(1:j),intfull%e1(1:j))
      et_first=e1_first+e2_first
      if (print_) write(fid,'(3(a,1x,f12.6,1x))')'E(first)',et_first,'E1(first)',e1_first,'E2(first)',e2_first ! debug
      !call flush(fid)
!      call initialize_to_identity_matrix(mo_coeff)
      ntsweep=0
      ntrot=0
      converged=0
      delrot=0.0_wp
      dele=0.0_wp
      nrot=0
      do isweep=1,20
        call jacobi_sweep(maxind,xpi,ypi,nrot,angtol,dele)
        if (print_) write(fid,'(a,1x,i3,5x,a,1x,i4,5x,a,1x,es12.3)')'sweep=',isweep,'nrot=',nrot,'dE=',dele ! debug
        !call flush(fid)
        if (abs(dele).le.mcetol) converged=1
        delrot=delrot+dele
        ntsweep=ntsweep+1
        ntrot=ntrot+nrot
        if (converged.eq.1) exit
      end do
      e1_last=0.0_wp
      e2_last=0.0_wp
      do i=1,nirrep
        j=0
        if (allocated(denfull%e2(i)%val)) j=size(denfull%e2(i)%val,dim=1)
        if (j.eq.0) cycle
        e2_last=e2_last+dot_product(denfull%e2(i)%val(1:j),intfull%e2(i)%val(1:j))
      end do
      j=size(denfull%e1)
      e1_last=dot_product(denfull%e1(1:j),intfull%e1(1:j))
      et_last=e1_last+e2_last
      if (print_) write(fid,'(3(a,1x,f12.6,1x))')'E(last)',et_last,'E1(last)',e1_last,'E2(last)',e2_last ! debug 
      !call flush(fid)
      deallocate(maxind,rotinds)
      return
    end subroutine jacobi_optimize

    subroutine print_square_matrix(matrix)
      implicit none
      real(wp), intent(in) :: matrix(:,:)
      integer, parameter :: maxprint=8
      integer :: ntot_l,i,j,jprint,js,k,kmax
      ntot_l=size(matrix,dim=1)
      jprint=ntot_l/maxprint
      js=0
      if (.not.print_) return
      write(fid,'(a)')'printing current mo transformation matrix'
      do j=1,jprint
        write(fid,'(a,1x,8(2x,i3,2x,1x))')'MO ',(k,k=js*maxprint+1,min((js+1)*maxprint,ntot_l))
        kmax=min((js+1)*maxprint,ntot_l)
        do i=1,ntot_l
          write(fid,'(i3,1x)',advance='no')i
          do k=js*maxprint+1,kmax
            write(fid,'(8(f7.3,1x))',advance='no')matrix(k,i)
          end do
          write(fid,*)
        end do
        js=js+1
      end do
      if (mod(ntot_l,maxprint).ne.0) then
        write(fid,'(a,1x,8(2x,i3,2x,1x))')'MO ',(k,k=js*maxprint+1,ntot_l)
        do i=1,ntot_l
          write(fid,'(i3,1x)',advance='no')i
          do k=js*maxprint+1,ntot_l
            write(fid,'(8(f7.3,1x))',advance='no')matrix(k,i)
          end do
          write(fid,*)
        end do
      end if
      return
    end subroutine print_square_matrix

    subroutine initialize_to_identity_matrix(matrix)
      real(wp) :: matrix(:,:)
      integer :: i
      matrix=0.0_wp
      do i=1,size(matrix,dim=1)
        matrix(i,i)=1.0_wp
      end do
      return
    end subroutine initialize_to_identity_matrix

    subroutine jacobi_rotinds()
      implicit none
      integer :: p,q,psym,qsym,nind
      npair=0
      do p=nfrozen+1,ncore
        psym=symmetries(p)
        do q=ntot,ncore+1,-1 ! core-active and core-virtual indeces
          if (symmetries(q).ne.psym) cycle
          npair=npair+1
        end do
      end do
      do p=ncore+1,nact+ncore
        psym=symmetries(p)
        do q=ntot,p+1,-1     ! active-active and active-virtual indeces
          if ((aarot.eq.0).and.(q.le.nact+ncore)) cycle ! exclude active-active orbital rotations
          if (symmetries(q).ne.psym) cycle
          npair=npair+1
        end do
      end do
      allocate(rotinds(2*npair))
      nind=0
      do p=nfrozen+1,ncore
        psym=symmetries(p)
        do q=ntot,ncore+1,-1 ! core-active and core-virtual indeces
          if (symmetries(q).ne.psym) cycle
          rotinds(nind+1)=p
          rotinds(nind+2)=q
          nind=nind+2
        end do
      end do
      do p=ncore+1,nact+ncore
        psym=symmetries(p)
        do q=ntot,p+1,-1     ! active-active and active-virtual indeces
          if ((aarot.eq.0).and.(q.le.nact+ncore)) cycle ! exclude active-active orbital rotations
          if (symmetries(q).ne.psym) cycle
          rotinds(nind+1)=p
          rotinds(nind+2)=q
          nind=nind+2
        end do
      end do
      return
    end subroutine jacobi_rotinds

    subroutine jacobi_sweep(indx,xpi,ypi,nrot,angtol,dele)
      implicit none
      real(wp) :: xpi(:),ypi(:)
      real(wp), intent(out) :: dele
      real(wp), intent(in) :: angtol
      integer :: indz(nirrep)
      integer :: i,j,nrot,iang,minpos
      real(wp) :: pairde,xmin,ymin
      logical :: found
      integer, dimension(:) :: indx(maxind_nnz,nirrep),inz(14,nirrep),ibv(4) ! integer sratch used by routines
      real(wp), dimension(:) :: a(14),b(9),c(9)             ! real scratch used by routines
      real(wp), dimension(:) :: atemp(14,nirrep)
      integer :: irrep,procwant,procuse
      nrot=0
      dele=0.0_wp
      do iang=0,2*npair-2,2 ! loop over rotation pairs
        i=rotinds(iang+1) ! orbital i
        j=rotinds(iang+2) ! orbital j (j>i)
        inz=0
        if (j.le.nact+ncore) then
          call jacobi_makeinds_aa(indx,inz,i,j,indz)
          procwant=count(indz.ne.0)
          procuse=min(nproc,procwant)
# ifdef OMP
!$omp parallel num_threads(procuse)
!$omp do private(irrep)
# endif
          do irrep=1,nirrep
            call jacobi_a_aa(intfull%e1(:),intfull%e2(irrep)%val(:),       &
               & denfull%e1(:),denfull%e2(irrep)%val(:),atemp(:,irrep),indx(:,irrep),inz(:,irrep))
          end do
# ifdef OMP
!$omp end do
!$omp end parallel
# endif
          a=0.0_wp
          do irrep=1,nirrep
            a=a+atemp(:,irrep)
          end do
          call jacobi_get_bc_aa(a,b,c)
        else
          call jacobi_makeinds_av(indx,inz,i,j,indz)
          procwant=count(indz.ne.0)
          procuse=min(nproc,procwant)
# ifdef OMP
!$omp parallel num_threads(procuse)
!$omp do private(irrep)
# endif
          do irrep=1,nirrep
            call jacobi_a_av(intfull%e1(:),intfull%e2(irrep)%val(:),       &
               & denfull%e1(:),denfull%e2(irrep)%val(:),atemp(:,irrep),indx(:,irrep),inz(:,irrep))
          end do
# ifdef OMP
!$omp end do
!$omp end parallel
# endif
          a=0.0_wp
          do irrep=1,nirrep
            a=a+atemp(:,irrep)
          end do
          call jacobi_get_bc_av(a,b,c) 
        end if
        call jacobi_bracketmin(ncore,i,a,b,c,xpi,ypi,ibv,minpos,pairde,found) ! bracket the angle that minimizes the energy
        if (.not.found) cycle
        call jacobi_refinemin(b,c,xpi(minpos),xpi(minpos+1),xmin,ymin,pairde) ! refine value of angle that minimizes the energy
        if (pairde.gt.angtol) cycle
        dele=dele+pairde
# ifdef OMP
!$omp parallel num_threads(procuse) shared(xmin,ymin)
!$omp do private(irrep)
# endif
        do irrep=1,nirrep
          call jacobi_rotateints(intfull%e1,intfull%e2(irrep)%val(:),indx(:,irrep),inz(:,irrep),xmin,ymin)
        end do
# ifdef OMP
!$omp end do
!$omp end parallel
# endif
        call update_mo_coeff(i,j,xmin,ymin)
        nrot=nrot+1
      end do
      return
    end subroutine jacobi_sweep

    subroutine update_mo_coeff(i,j,x,y)
      implicit none
      integer, intent(in) :: i,j
      real(wp), intent(in) :: x,y
      real(wp) :: phi_i(ntot),phi_j(ntot)
      ! save original vectors
      phi_i=mo_coeff(:,i)
      phi_j=mo_coeff(:,j)
      ! update vectors y = cos(theta) && x = sin(theta)
      mo_coeff(:,i)=  y*phi_i + x*phi_j
      mo_coeff(:,j)= -x*phi_i + y*phi_j
      return
    end subroutine update_mo_coeff

    subroutine update_mo_coeff_full(i,j,x,y,matrix)
      implicit none
      integer :: i,j
      real(wp) :: x,y
      real(wp) :: matrix(ntot,ntot),cmat(ntot,ntot)
      call initialize_to_identity_matrix(cmat)
      cmat(i,i)=y
      cmat(j,j)=y
      cmat(i,j)=x
      cmat(j,i)=-x
      matrix=matmul(cmat,matrix)
      return
    end subroutine update_mo_coeff_full

    subroutine jacobi_makeinds_av(indx,inz,i,j,indz)
      implicit none
      integer, intent(in) :: i,j
      integer, intent(inout) :: indx(:,:),inz(:,:)
      integer :: indz(nirrep),isave
      integer :: p,q,r,iis,ijs,ips,r1,is1,is2,gli,glj,gri,grj,ioff,na,nt
      iis=symmetries(i)
      ijs=symmetries(j)
      inz=0
      na=ncore+nact
      nt=ntot
      r1=group_mult_tab(iis,ijs)
      indz=0
      indx(1,1)=redind(i,i)
      indx(2,1)=redind(i,j)
      indx(3,1)=redind(j,j)
      indz(1) = 3
      inz(1,:)=indz
      do p=1,i-1
         if (iis.ne.symmetries(p)) cycle
         isave=indz(1)
         indx(isave+1,1)=redind(i,p) ! (i|p)
         indx(isave+2,1)=redind(j,p) ! (j|p)
         indz(1) = indz(1) + 2
      end do
      do p=i+1,na
        if (iis.ne.symmetries(p)) cycle
         isave=indz(1)
         indx(isave+1,1)=redind(i,p) ! (i|p)
         indx(isave+2,1)=redind(j,p) ! (j|p)
         indz(1) = indz(1) + 2
      end do
      inz(2,:) = indz
      do r=1,i-1
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=1,i-1
          ips=symmetries(p)
          do q=1,p
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,na
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=1,r-1
          if (p.eq.i) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is2.ne.is1) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
          indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=i+1,na
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          gri=redind(p,q)
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            gli=redind(i,r)
            glj=redind(j,r)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,na
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do q=i+1,r
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|qr)
          indx(isave+2,is1)=index(glj,gri) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,na
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is2.ne.is1) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri)  ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri)  ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      inz(3,:) = indz
      do p=1,i-1
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,ips)
        do q=1,p-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri)  ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri)  ! (jp|iq)
          gri=redind(j,q)        
          indx(isave+3,is1)=index(gli,gri)  ! (ip|jq)
          indx(isave+4,is1)=index(glj,gri)  ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=i+1,na
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,ips)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri)  ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri)  ! (jp|iq)
          gri=redind(j,q)
          indx(isave+3,is1)=index(gli,gri)  ! (ip|jq)
          indx(isave+4,is1)=index(glj,gri)  ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=i+1,na
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,ips)
        do q=i+1,p-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri)  ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri)  ! (jp|iq)
          gri=redind(j,q)
          indx(isave+3,is1)=index(gli,gri)  ! (ip|jq)
          indx(isave+4,is1)=index(glj,gri)  ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      inz(4,:) = indz
      do p=1,i-1
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(symmetries(p),iis)
        isave=indz(is1)
        indx(isave+1,is1)=index(gli,gli) ! (ip|ip)
        indx(isave+2,is1)=index(gli,glj) ! (ip|jp)
        indx(isave+3,is1)=index(glj,glj) ! (jp|jp)
        indz(is1) = indz(is1) + 3
      end do
      do p=i+1,na
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(symmetries(p),iis)
        isave=indz(is1)
        indx(isave+1,is1)=index(gli,gli) ! (ip|ip)
        indx(isave+2,is1)=index(gli,glj) ! (ip|jp)
        indx(isave+3,is1)=index(glj,glj) ! (jp|jp)
        indz(is1) = indz(is1) + 3
      end do
      inz(5,:) = indz
      gli=redind(i,i)
      glj=redind(j,j)
      grj=redind(i,j)
      do p=1,i-1
        ips=symmetries(p)
        do q=1,p
          if (ips.ne.symmetries(q)) cycle
          gri=redind(p,q)
          isave=indz(1)
          indx(isave+1,1)=index(gli,gri) ! (ii|pq)
          indx(isave+2,1)=index(grj,gri) ! (ij|pq) 
          indx(isave+3,1)=index(glj,gri) ! (jj|pq)
          indz(1) = indz(1) + 3
        end do
      end do
      do p=i+1,na
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          if (ips.ne.symmetries(q)) cycle
          gri=redind(p,q)
          isave=indz(1)
          indx(isave+1,1)=index(gli,gri) ! (ii|pq)
          indx(isave+2,1)=index(grj,gri) ! (ij|pq)
          indx(isave+3,1)=index(glj,gri) ! (jj|pq)
          indz(1) = indz(1) + 3
        end do
      end do
      inz(6,:) = indz
      do p=1,i-1
        if (iis.ne.symmetries(p)) cycle
        gri=redind(i,p)
        isave=indz(1)
        indx(isave+1,1)=index(gli,gri) ! (ii|ip)
        indx(isave+2,1)=index(grj,gri) ! (ij|ip)
        indx(isave+3,1)=index(glj,gri) ! (jj|ip)
        gri=redind(j,p)
        indx(isave+4,1)=index(gli,gri) ! (ii|ip)
        indx(isave+5,1)=index(grj,gri) ! (ij|ip)
        indx(isave+6,1)=index(glj,gri) ! (jj|ip)
        indz(1) = indz(1) + 6
      end do
      do p=i+1,na
        if (iis.ne.symmetries(p)) cycle
        gri=redind(i,p)
        isave=indz(1)
        indx(isave+1,1)=index(gli,gri) ! (ii|ip)
        indx(isave+2,1)=index(grj,gri) ! (ij|ip)
        indx(isave+3,1)=index(glj,gri) ! (jj|ip)
        gri=redind(j,p)
        indx(isave+4,1)=index(gli,gri) ! (ii|ip)
        indx(isave+5,1)=index(grj,gri) ! (ij|ip)
        indx(isave+6,1)=index(glj,gri) ! (jj|ip)
        indz(1) = indz(1) + 6
      end do
      inz(7,:) = indz
      isave=indz(1)
      indx(isave+1,1)=index(gli,gli) ! (ii|ii)
      indx(isave+2,1)=index(grj,gli) ! (ij|ii)
      indx(isave+3,1)=index(gli,glj) ! (ii|jj)
      indx(isave+4,1)=index(grj,grj) ! (ij|ij)
      indx(isave+5,1)=index(grj,glj) ! (ij|jj)
      indx(isave+6,1)=index(glj,glj) ! (jj|jj)
      indz(1) = indz(1) + 6
      inz(8,:) = indz
      do p=na+1,j-1
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        indx(isave+1,1)=redind(i,p) ! (i|p)
        indx(isave+2,1)=redind(j,p) ! (j,p)
        indz(1) = indz(1) + 2
      end do
      do p=j+1,nt
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        indx(isave+1,1)=redind(i,p) ! (i|p)
        indx(isave+2,1)=redind(j,p) ! (j,p)
        indz(1) = indz(1) + 2
      end do
      inz(9,:) = indz
      do r=na+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=1,r-1
          if (p.eq.i) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(r,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|rq)
          indx(isave+2,is1)=index(glj,gri) ! (jr|rq)
          indz(is1) = indz(is1) + 2
        end do
      end do
      do r=j+1,nt
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=1,r-1
          if (p.eq.i.or.p.eq.j) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(r,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|rq)
          indx(isave+2,is1)=index(glj,gri) ! (jr|rq)
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=na+1,j-1
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          gri=redind(p,q)
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            gli=redind(i,r)
            glj=redind(j,r)
            isave=indz(is1)
            indx(isave+1,is1) =index(gli,gri) ! (ir,pq)
            indx(isave+2,is1) =index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          gri=redind(p,q)
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            gli=redind(i,r)
            glj=redind(j,r)
            isave=indz(is1)
            indx(isave+1,is1) =index(gli,gri) ! (ir,pq)
            indx(isave+2,is1) =index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=j+1,nt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,nt
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do q=j+1,r
          if (iis.ne.symmetries(q)) cycle
          gri=redind(r,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|rq)
          indx(isave+2,is1)=index(glj,gri) ! (jr|rq)
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,nt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,nt
        gli=redind(i,r)
        glj=redind(j,r)
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,j-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(r,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|rq)
          indx(isave+2,is1)=index(glj,gri) ! (jr|rq)
          indz(is1) = indz(is1) + 2
        end do
      end do
      do r=i+1,na
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=na+1,j-1
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=na+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do q=i+1,r
          if (iis.ne.symmetries(q)) cycle
          gri=redind(r,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|rq)
          indx(isave+2,is1)=index(glj,gri) ! (jr|rq)
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,j-1
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      inz(10,:) = indz
      do p=na+1,j-1
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(ips,iis)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          gri=redind(j,q)
          indx(isave+3,is1)=index(gli,gri) ! (ip,jq)
          indx(isave+4,is1)=index(glj,gri) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,ips)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          isave=indz(is1)
          gri=redind(i,q)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          gri=redind(j,q)
          indx(isave+3,is1)=index(gli,gri) ! (ip,jq)
          indx(isave+4,is1)=index(glj,gri) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=na+1,j-1
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,ips)
        do q=i+1,p-1
          if (ips.ne.symmetries(q)) cycle
          isave=indz(is1)
          gri=redind(i,q)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          gri=redind(j,q)
          indx(isave+3,is1)=index(gli,gri) ! (ip,jq)
          indx(isave+4,is1)=index(glj,gri) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,ips)
        do q=i+1,j-1
          if (ips.ne.symmetries(q)) cycle
          isave=indz(is1)
          gri=redind(i,q)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          gri=redind(j,q)
          indx(isave+3,is1)=index(gli,gri) ! (ip,jq)
          indx(isave+4,is1)=index(glj,gri) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,ips)
        do q=j+1,p-1
          if (ips.ne.symmetries(q)) cycle
          isave=indz(is1)
          gri=redind(i,q)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          gri=redind(j,q)
          indx(isave+3,is1)=index(gli,gri) ! (ip,jq)
          indx(isave+4,is1)=index(glj,gri) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      inz(11,:) = indz
      do p=na+1,j-1
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,symmetries(p))
        isave=indz(is1)
        indx(isave+1,is1)=index(gli,gli) ! (ip|ip)
        indx(isave+2,is1)=index(gli,glj) ! (ip|jp)
        indx(isave+3,is1)=index(glj,glj) ! (jp|jp)
        indz(is1) = indz(is1) + 3
      end do
      do p=j+1,nt
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,symmetries(p))
        isave=indz(is1)
        indx(isave+1,is1)=index(gli,gli) ! (ip|ip)
        indx(isave+2,is1)=index(gli,glj) ! (ip|jp)
        indx(isave+3,is1)=index(glj,glj) ! (jp|jp)
        indz(is1) = indz(is1) + 3
      end do
      inz(12,:) = indz
      gli=redind(i,i)
      glj=redind(j,j)
      grj=redind(i,j)
      do p=na+1,j-1
        ips=symmetries(p)
        do q=1,p
          iF (q.eq.i) cycle
          if (ips.ne.symmetries(q)) cycle
          gri=redind(p,q)
          isave=indz(1)
          indx(isave+1,1)=index(gli,gri) ! (ii|pq)
          indx(isave+2,1)=index(grj,gri) ! (ij|pq)
          indx(isave+3,1)=index(glj,gri) ! (jj|pq)
          indz(1) = indz(1) + 3
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        do q=1,p
          iF (q.eq.i.or.q.eq.j) cycle
          if (ips.ne.symmetries(q)) cycle
          gri=redind(p,q)
          isave=indz(1)
          indx(isave+1,1)=index(gli,gri) ! (ii|pq)
          indx(isave+2,1)=index(grj,gri) ! (ij|pq)
          indx(isave+3,1)=index(glj,gri) ! (jj|pq)
          indz(1) = indz(1) + 3
        end do
      end do
      inz(13,:) = indz
      do p=na+1,j-1
        if (iis.ne.symmetries(p)) cycle
        gri=redind(i,p)
        isave=indz(1)
        indx(isave+1,1)=index(gli,gri) ! (ii|ip)
        indx(isave+2,1)=index(grj,gri) ! (ij|ip)
        indx(isave+3,1)=index(glj,gri) ! (jj|ip)
        gri=redind(j,p)
        indx(isave+4,1)=index(gli,gri) ! (ii|jp)
        indx(isave+5,1)=index(grj,gri) ! (ij|jp)
        indx(isave+6,1)=index(glj,gri) ! (jj|jp)
        indz(1) = indz(1) + 6
      end do
      do p=j+1,nt
        if (iis.ne.symmetries(p)) cycle
        gri=redind(i,p)
        isave=indz(1)
        indx(isave+1,1)=index(gli,gri) ! (ii|ip)
        indx(isave+2,1)=index(grj,gri) ! (ij|ip)
        indx(isave+3,1)=index(glj,gri) ! (jj|ip)
        gri=redind(j,p)
        indx(isave+4,1)=index(gli,gri) ! (ii|jp)
        indx(isave+5,1)=index(grj,gri) ! (ij|jp)
        indx(isave+6,1)=index(glj,gri) ! (jj|jp)
        indz(1) = indz(1) + 6
      end do
      inz(14,:) = indz
      return
    end subroutine jacobi_makeinds_av

    subroutine jacobi_makeinds_aa(indx,inz,i,j,indz)
      implicit none
      integer, intent(inout) :: indx(:,:),inz(:,:)
      integer, intent(in) :: i,j
      integer :: indz(nirrep)
      integer :: p,q,r,iis,ijs,ips,r1,is1,is2,gli,glj,gri,grj,ioff,na,nt,isave
      iis=symmetries(i)
      ijs=symmetries(j)
      inz=0
      indz=0
      na=ncore+nact
      nt=ntot
      r1=group_mult_tab(iis,ijs)
! this routine is the jacobi source from GaMESS (mcjac.src) with minor modifications
! should only be called if both i and j belong to the active space
! 2) replaced GOtO statements with cycle's
! the routine computes two types of 1-E indeces (i|p) where p<=na (type 1 ... inz(1)) or p>na (type 2 ... inz(8))
!                   and two types of 2-E indeces
! (ij|pq) (ip|qr) (ir|jr) where two cases are possible
! case 1 ... all non-active indeces are within the active space inz(2:7)
! case 2 ... at least one of the non-active indeces is in the virtual space inz(9:13)
      indx(1,1)=redind(i,i) ! (i|i)
      indx(2,1)=redind(i,j) ! (j,j)
      indx(3,1)=redind(j,j) ! (i,j)
      indz(1)=3
      inz(1,:)=indz
      do p=1,i-1
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        indx(isave+1,1)=redind(i,p) ! (i|p)
        indx(isave+2,1)=redind(j,p) ! (j|p)
        indz(1) = indz(1) + 2
      end do
      do p=i+1,j-1
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        indx(isave+1,1)=redind(i,p) ! (i|p)
        indx(isave+2,1)=redind(j,p) ! (j|p)
        indz(1) = indz(1) + 2
      end do
      do p=j+1,na
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        indx(isave+1,1)=redind(i,p) ! (i|p)
        indx(isave+2,1)=redind(j,p) ! (j|p)
        indz(1) = indz(1) + 2
      end do
      inz(2,:) = indz
      do r=1,i-1
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(r,i)
        glj=redind(r,j)
        do p=1,i-1
          ips=symmetries(p)
          do q=1,p
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=1,r-1
          if (p.eq.i) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gri,gli) ! (ir|pq)
            indx(isave+2,is1)=index(gri,glj) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ir|qr)
          indx(isave+2,is1)=index(gri,glj) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
      end do
      do r=j+1,na
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(r,i)
        glj=redind(r,j)
        do p=1,r-1
          if (p.eq.i.or.p.eq.j) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1) = index(gli,gri) ! (ir,qr)
          indx(isave+2,is1) = index(glj,gri) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=i+1,j-1
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          gri=redind(p,q)
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            gli=redind(i,r)
            glj=redind(j,r)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        do q=1,p
          gri=redind(p,q)
          if (q.eq.i.or.q.eq.j) cycle
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            gli=redind(r,i)
            glj=redind(r,j)
            isave=indz(is1)
            indx(isave+1,is1) = index(gri,gli) ! (ir|pq)
            indx(isave+2,is1) = index(gri,glj) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(r,i)
        glj=redind(r,j)
        do p=j+1,na
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gri,gli) ! (ir|pq)
            indx(isave+2,is1)=index(gri,glj) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,na
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do q=j+1,r
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ir|qr)
          indx(isave+2,is1)=index(gri,glj) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,na
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (ir|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,na
        gli=redind(i,r)
        glj=redind(j,r)
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,j-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ir|qr)
          indx(isave+2,is1)=index(gri,glj) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
      end do
      do r=i+1,j-1
        gli=redind(i,r)
        glj=redind(j,r)
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,r
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ir|qr)
          indx(isave+2,is1)=index(gri,glj) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,j-1
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is2.ne.is1) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gri,gli) ! (ir|pq)
            indx(isave+2,is1)=index(gri,glj) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      inz(3,:) = indz
      do p=1,i-1
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(ips,iis)
        do q=1,p-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          indx(isave+3,is1)=index(gli,grj) ! (ip|jq)
          indx(isave+4,is1)=index(glj,grj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=i+1,j-1
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        gli=redind(i,p)
        glj=redind(j,p)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ip|iq)
          indx(isave+2,is1)=index(gri,glj) ! (ip|jq)
          indx(isave+3,is1)=index(grj,gli) ! (jp|iq)
          indx(isave+4,is1)=index(grj,glj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis) 
        gli=redind(i,p)
        glj=redind(j,p)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ip|iq)
          indx(isave+2,is1)=index(gri,glj) ! (ip|jq)
          indx(isave+3,is1)=index(grj,gli) ! (jp|iq)
          indx(isave+4,is1)=index(grj,glj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=i+1,j-1
        gli=redind(i,p)
        glj=redind(j,p)
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        do q=i+1,p-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ip|iq)
          indx(isave+2,is1)=index(gri,glj) ! (ip|jq)
          indx(isave+3,is1)=index(grj,gli) ! (jp|iq)
          indx(isave+4,is1)=index(grj,glj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        gli=redind(i,p)
        glj=redind(j,p)
        do q=i+1,j-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ip|iq)
          indx(isave+2,is1)=index(gri,glj) ! (ip|jq)
          indx(isave+3,is1)=index(grj,gli) ! (jp|iq)
          indx(isave+4,is1)=index(grj,glj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,na
        gli=redind(i,p)
        glj=redind(j,p)
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=j+1,p-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gri,gli) ! (ip|iq)
          indx(isave+2,is1)=index(gri,glj) ! (ip|jq)
          indx(isave+3,is1)=index(grj,gli) ! (jp|iq)
          indx(isave+4,is1)=index(grj,glj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      inz(4,:) = indz
      do p=1,i-1
        is1=group_mult_tab(symmetries(p),iis)
        gli=redind(i,p)
        glj=redind(j,p)
        isave=indz(is1)
        indx(isave+1,is1)=index(gli,gli) ! (ip|ip)
        indx(isave+2,is1)=index(gli,glj) ! (ip|jp)
        indx(isave+3,is1)=index(glj,glj) ! (jp|jp)
        indz(is1) = indz(is1) + 3
      end do
      do p=i+1,j-1
        is1=group_mult_tab(symmetries(p),iis)
        gli=redind(i,p)
        glj=redind(j,p)
        isave=indz(is1)
        indx(isave+1,is1)=index(gli,gli) ! (ip|ip)
        indx(isave+2,is1)=index(gli,glj) ! (ip|jp)
        indx(isave+3,is1)=index(glj,glj) ! (jp|jp)
        indz(is1) = indz(is1) + 3
      end do
      do p=j+1,na
        is1=group_mult_tab(symmetries(p),iis)
        gli=redind(i,p)
        glj=redind(j,p)
        isave=indz(is1)
        indx(isave+1,is1)=index(gli,gli) ! (ip|ip)
        indx(isave+2,is1)=index(gli,glj) ! (ip|jp)
        indx(isave+3,is1)=index(glj,glj) ! (jp|jp)
        indz(is1) = indz(is1) + 3
      end do
      inz(5,:) = indz
      gli=redind(i,i)
      glj=redind(j,j)
      grj=redind(i,j)
      do p=1,i-1
        ips=symmetries(p)
        do q=1,p
          if (ips.ne.symmetries(q)) cycle
          gri=redind(p,q)
          isave=indz(1)
          indx(isave+1,1)=index(gri,gli) ! (ii|pq)
          indx(isave+2,1)=index(gri,grj) ! (ij|pq)
          indx(isave+3,1)=index(gri,glj) ! (jj|pq)
          indz(1) = indz(1) + 3
        end do
      end do
      do p=i+1,j-1
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          if (ips.ne.symmetries(q)) cycle
          gri=redind(p,q)
          isave=indz(1)
          indx(isave+1,1)=index(gri,gli) ! (ii|pq)
          indx(isave+2,1)=index(gri,grj) ! (ij|pq)
          indx(isave+3,1)=index(gri,glj) ! (jj|pq)
          indz(1) = indz(1) + 3
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          if (ips.ne.symmetries(q)) cycle
          gri=redind(p,q)
          isave=indz(1)
          indx(isave+1,1)=index(gri,gli) ! (ii|pq)
          indx(isave+2,1)=index(gri,grj) ! (ij|pq)
          indx(isave+3,1)=index(gri,glj) ! (jj|pq)
          indz(1) = indz(1) + 3
        end do
      end do
      inz(6,:) = indz
      do p=1,i-1
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        gri=redind(i,p)
        indx(isave+1,1)=index(gri,gli) ! (ip|ii)
        indx(isave+2,1)=index(gri,grj) ! (ip|ij)
        indx(isave+3,1)=index(gri,glj) ! (ip|jj)
        gri=redind(j,p)
        indx(isave+4,1)=index(gli,gri) ! (jp|ii)
        indx(isave+5,1)=index(grj,gri) ! (jp|ij)
        indx(isave+6,1)=index(glj,gri) ! (jp|jj)
        indz(1) = indz(1) + 6
      end do
      do p=i+1,j-1
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        gri=redind(i,p)
        indx(isave+1,1)=index(gri,gli) ! (ip|ii)
        indx(isave+2,1)=index(gri,grj) ! (ip|ij)
        indx(isave+3,1)=index(gri,glj) ! (ip|jj)
        gri=redind(j,p)
        indx(isave+4,1)=index(gli,gri) ! (jp|ii)
        indx(isave+5,1)=index(grj,gri) ! (jp|ij)
        indx(isave+6,1)=index(glj,gri) ! (jp|jj)
        indz(1) = indz(1) + 6
      end do
      do p=j+1,na
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        gri=redind(i,p)
        indx(isave+1,1)=index(gri,gli) ! (ip|ii)
        indx(isave+2,1)=index(gri,grj) ! (ip|ij)
        indx(isave+3,1)=index(gri,glj) ! (ip|jj)
        gri=redind(j,p)
        indx(isave+4,1)=index(gli,gri) ! (jp|ii)
        indx(isave+5,1)=index(grj,gri) ! (jp|ij)
        indx(isave+6,1)=index(glj,gri) ! (jp|jj)
        indz(1) = indz(1) + 6
      end do
      inz(7,:) = indz
      gli=redind(i,i)
      glj=redind(j,j)
      gri=redind(i,j)
      isave=indz(1)
      indx(isave+1,1)=index(gli,gli) ! (ii|ii)
      indx(isave+2,1)=index(gli,gri) ! (ii|ij)
      indx(isave+3,1)=index(gli,glj) ! (ii|jj)
      indx(isave+4,1)=index(gri,gri) ! (ij|ij)
      indx(isave+5,1)=index(gri,glj) ! (ij,jj)
      indx(isave+6,1)=index(glj,glj) ! (jj,jj)
      indz(1) = indz(1) + 6
      inz(8,:) = indz
      do p=na+1,nt
        if (iis.ne.symmetries(p)) cycle
        isave=indz(1)
        indx(isave+1,1)=redind(i,p) ! (i|p)
        indx(isave+2,1)=redind(j,p) ! (j,p)
        indz(1) = indz(1) + 2
      end do
      inz(9,:) = indz
      do r=na+1,nt
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=1,r-1
          if (p.eq.i.or.p.eq.j) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(r,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|qr)
          indx(isave+2,is1)=index(glj,gri) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=na+1,nt
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          is2=group_mult_tab(ips,symmetries(q))
          gri=redind(p,q)
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            gli=redind(i,r)
            glj=redind(j,r)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=na+1,nt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,na
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do p=na+1,nt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
           end do
        end do
      end do
      do r=na+1,nt
        is1=group_mult_tab(iis,symmetries(r))
        gli=redind(i,r)
        glj=redind(j,r)
        do q=j+1,r
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|qr)
          indx(isave+2,is1)=index(glj,gri) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,nt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            gri=redind(p,q)
            isave=indz(is1)
            indx(isave+1,is1)=index(gli,gri) ! (ir|pq)
            indx(isave+2,is1)=index(glj,gri) ! (jr|pq)
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=na+1,nt
        gli=redind(i,r)
        glj=redind(j,r)
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,j-1
          if (iis.ne.symmetries(q)) cycle
          gri=redind(q,r)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ir|qr)
          indx(isave+2,is1)=index(glj,gri) ! (jr|qr)
          indz(is1) = indz(is1) + 2
        end do
      end do
      inz(10,:) = indz
      do p=na+1,nt
        gli=redind(i,p)
        glj=redind(j,p)
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          indx(isave+3,is1)=index(gli,grj) ! (ip|jq)
          indx(isave+4,is1)=index(glj,grj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=na+1,nt
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(ips,iis)
        do q=i+1,j-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          indx(isave+3,is1)=index(gli,grj) ! (ip|jq)
          indx(isave+4,is1)=index(glj,grj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=na+1,nt
        ips=symmetries(p)
        gli=redind(i,p)
        glj=redind(j,p)
        is1=group_mult_tab(iis,ips)
        do q=j+1,p-1
          if (ips.ne.symmetries(q)) cycle
          gri=redind(i,q)
          grj=redind(j,q)
          isave=indz(is1)
          indx(isave+1,is1)=index(gli,gri) ! (ip|iq)
          indx(isave+2,is1)=index(glj,gri) ! (jp|iq)
          indx(isave+3,is1)=index(gli,grj) ! (ip|jq)
          indx(isave+4,is1)=index(glj,grj) ! (jp|jq)
          indz(is1) = indz(is1) + 4
        end do
      end do
      inz(11,:) = indz
      do p=na+1,nt
        is1=group_mult_tab(iis,symmetries(p))
        gri=redind(i,p)
        grj=redind(j,p)
        isave=indz(is1)
        indx(isave+1,is1)=index(gri,gri) ! (ip|ip)
        indx(isave+2,is1)=index(gri,grj) ! (ip|jp)
        indx(isave+3,is1)=index(grj,grj) ! (jp|jp)
        indz(is1) = indz(is1) + 3
      end do
      inz(12,:) = indz
      gli=redind(i,i)
      glj=redind(j,j)
      grj=redind(i,j)
      do p=na+1,nt
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          if (ips.ne.symmetries(q)) cycle
          gri=redind(p,q)
          isave=indz(1)
          indx(isave+1,1)=index(gli,gri) ! (ii|pq)
          indx(isave+2,1)=index(grj,gri) ! (ij|pq)
          indx(isave+3,1)=index(glj,gri) ! (jj|pq)
          indz(1) = indz(1) + 3
        end do
      end do
      inz(13,:) = indz
      do p=na+1,nt
        if (iis.ne.symmetries(p)) cycle
        gri=redind(i,p)
        isave=indz(1)
        indx(isave+1,1)=index(gli,gri) ! (ii|ip)
        indx(isave+2,1)=index(grj,gri) ! (ij|ip)
        indx(isave+3,1)=index(glj,gri) ! (jj,ip)
        gri=redind(j,p)
        indx(isave+4,1)=index(gli,gri) ! (jj|jp)
        indx(isave+5,1)=index(grj,gri) ! (ij|jp)
        indx(isave+6,1)=index(glj,gri) ! (jj|jp)
        indz(1) = indz(1) + 6
      end do
      inz(14,:) = indz
      return
    end subroutine jacobi_makeinds_aa

    subroutine jacobi_rotateints(h1,g,indx,inz,x,y)
      implicit none
      real(wp), intent(in) :: x,y
      integer, intent(in)  :: indx(:),inz(:)
      real(wp), intent(inout) :: h1(:),g(:)
      integer :: i,i1,i2,i3,i4,i5,i6
      real(wp) :: x2,xy,t1,t2,y3,x3,xy2,t3,x2y,t4,t5,t6,x4,y4,x2y2
      real(wp) :: xy3,x3y,u1,u2,u3,u4,u5,u6,u7,u8,c1,c2,c3,c4,c5,c6,y2
      x2 = x*x
      y2 = y*y
      xy = x*y
      t1 = xy*2.0_wp
      t2 = (y2-x2)
      y3 = y2*y
      x3 = x2*x
      xy2 = xy*y
      t3 = 2.0_wp*xy2
      x2y = xy*x
      t4 = 2.0_wp*x2y
      t5 = (y3 - x2y)
      t6 = (x3 - xy2)
      x4 = x3*x
      y4 = y3*y
      x2y2 = x2*y2
      xy3 = xy2*y
      x3y = x2y*x
      u1 = 4.0_wp*xy3
      u2 = 2.0_wp*x2y2
      u3 = 2.0_wp*u2
      u4 = 4.0_wp*x3y
      u5 = y4 - 3.0_wp*x2y2
      u6 = x3y - xy3
      u7 = 2.0_wp*u6
      u8 = 3.0_wp*x2y2 - x4
      do i=1,inz(1),3
        i1 = indx(1)
        i2 = indx(2)
        i3 = indx(3)
        c1 = h1(i1)
        c2 = h1(i2)
        c3 = h1(i3)
        h1(i1) = y2*c1 - t1*c2 + x2*c3
        h1(i2) = xy*(c1-c3) + t2*c2
        h1(i3) = x2*c1 + t1*c2 + y2*c3
      end do
      do i=inz(1)+1,inz(2),2
         i1 = indx(i)
         i2 = indx(i+1)
         c1 = h1(i1)
         c2 = h1(i2)
         h1(i1) = y*c1 - x*c2
         h1(i2) = x*c1 + y*c2
      end do
      do i=inz(8)+1,inz(9),2
         i1 = indx(i)
         i2 = indx(i+1)
         c1 = h1(i1)
         c2 = h1(i2)
         h1(i1) = y*c1 - x*c2
         h1(i2) = x*c1 + y*c2
      end do
      do i=inz(2)+1,inz(3),2
         i1 = indx(i)
         i2 = indx(i+1)
         c1 = g(i1)
         c2 = g(i2)
         g(i1) = y*c1 - x*c2
         g(i2) = x*c1 + y*c2
      end do
      do i=inz(9)+1,inz(10),2
         i1 = indx(i)
         i2 = indx(i+1)
         c1 = g(i1)
         c2 = g(i2)
         g(i1) = y*c1 - x*c2
         g(i2) = x*c1 + y*c2
      end do
      do i=inz(3)+1,inz(4),4
         i1 = indx(i)
         i2 = indx(i+1)
         i3 = indx(i+2)
         i4 = indx(i+3)
         c1 = g(i1)
         c2 = g(i2)
         c3 = g(i3)
         c4 = g(i4)
         g(i1) = y2*c1 - xy*(c2+c3) + x2*c4
         g(i2) = xy*(c1-c4) + y2*c2 - x2*c3
         g(i3) = xy*(c1-c4) - x2*c2 + y2*c3
         g(i4) = x2*c1 + xy*(c2+c3) + y2*c4
      end do
      do i=inz(10)+1,inz(11),4
         i1 = indx(i)
         i2 = indx(i+1)
         i3 = indx(i+2)
         i4 = indx(i+3)
         c1 = g(i1)
         c2 = g(i2)
         c3 = g(i3)
         c4 = g(i4)
         g(i1) = y2*c1 - xy*(c2+c3) + x2*c4
         g(i2) = xy*(c1-c4) + y2*c2 - x2*c3
         g(i3) = xy*(c1-c4) - x2*c2 + y2*c3
         g(i4) = x2*c1 + xy*(c2+c3) + y2*c4
      end do
      do i=inz(4)+1,inz(5),3
         i1 = indx(i)
         i2 = indx(i+1)
         i3 = indx(i+2)
         c1 = g(i1)
         c2 = g(i2)
         c3 = g(i3)
         g(i1) = y2*c1 - t1*c2 + x2*c3
         g(i2) = xy*(c1-c3) + t2*c2
         g(i3) = x2*c1 + t1*c2 + y2*c3
      end do
      do i=inz(11)+1,inz(12),3
         i1 = indx(i)
         i2 = indx(i+1)
         i3 = indx(i+2)
         c1 = g(i1)
         c2 = g(i2)
         c3 = g(i3)
         g(i1) = y2*c1 - t1*c2 + x2*c3
         g(i2) = xy*(c1-c3) + t2*c2
         g(i3) = x2*c1 + t1*c2 + y2*c3
      end do
      do i=inz(5)+1,inz(6),3
         i1 = indx(i)
         i2 = indx(i+1)
         i3 = indx(i+2)
         c1 = g(i1)
         c2 = g(i2)
         c3 = g(i3)
         g(i1) = y2*c1 - t1*c2 + x2*c3
         g(i2) = xy*(c1-c3) + t2*c2
         g(i3) = x2*c1 + t1*c2 + y2*c3
      end do
      do i=inz(12)+1,inz(13),3
         i1 = indx(i)
         i2 = indx(i+1)
         i3 = indx(i+2)
         c1 = g(i1)
         c2 = g(i2)
         c3 = g(i3)
         g(i1) = y2*c1 - t1*c2 + x2*c3
         g(i2) = xy*(c1-c3) + t2*c2
         g(i3) = x2*c1 + t1*c2 + y2*c3
      end do
      do i=inz(6)+1,inz(7),6
         i1 = indx(i)
         i2 = indx(i+1)
         i3 = indx(i+2)
         i4 = indx(i+3)
         i5 = indx(i+4)
         i6 = indx(i+5)
         c1 = g(i1)
         c2 = g(i2)
         c3 = g(i3)
         c4 = g(i4)
         c5 = g(i5)
         c6 = g(i6)
         g(i1) = y3*c1-xy2*c4-t3*c2+x2y*c3+t4*c5-x3*c6
         g(i2) = xy2*c1-x2y*c4+t5*c2-xy2*c3+t6*c5+x2y*c6
         g(i3) = x2y*c1-x3*c4+t3*c2+y3*c3-t4*c5-xy2*c6
         g(i4) = xy2*c1+y3*c4-t4*c2+x3*c3-t3*c5+x2y*c6
         g(i5) = x2y*c1+xy2*c4-t6*c2-x2y*c3+t5*c5-xy2*c6
         g(i6) = x3*c1+x2y*c4+t4*c2+xy2*c3+t3*c5+y3*c6
      end do
      do i=inz(13)+1,inz(14),6
         i1 = indx(i)
         i2 = indx(i+1)
         i3 = indx(i+2)
         i4 = indx(i+3)
         i5 = indx(i+4)
         i6 = indx(i+5)
         c1 = g(i1)
         c2 = g(i2)
         c3 = g(i3)
         c4 = g(i4)
         c5 = g(i5)
         c6 = g(i6)
         g(i1) = y3*c1-xy2*c4-t3*c2+x2y*c3+t4*c5-x3*c6
         g(i2) = xy2*c1-x2y*c4+t5*c2-xy2*c3+t6*c5+x2y*c6
         g(i3) = x2y*c1-x3*c4+t3*c2+y3*c3-t4*c5-xy2*c6
         g(i4) = xy2*c1+y3*c4-t4*c2+x3*c3-t3*c5+x2y*c6
         g(i5) = x2y*c1+xy2*c4-t6*c2-x2y*c3+t5*c5-xy2*c6
         g(i6) = x3*c1+x2y*c4+t4*c2+xy2*c3+t3*c5+y3*c6
      end do
      do i=inz(7)+1,inz(8),6
        i1 = indx(i)
        i2 = indx(i+1)
        i3 = indx(i+2)
        i4 = indx(i+3)
        i5 = indx(i+4)
        i6 = indx(i+5)
        c1 = g(i1)
        c2 = g(i2)
        c3 = g(i3)
        c4 = g(i4)
        c5 = g(i5)
        c6 = g(i6)
        g(i1) = y4*c1-u1*c2+u2*c3+u3*c4-u4*c5+x4*c6
        g(i2) = xy3*c1+u5*c2+u6*c3+u7*c4+u8*c5-x3y*c6
        g(i3) = x2y2*c1-u7*c2+(x4+y4)*c3-u3*c4+u7*c5+x2y2*c6
        g(i4) = x2y2*c1-u7*c2-u2*c3+(x4+y4-u2)*c4+u7*c5+x2y2*c6
        g(i5) = x3y*c1+u8*c2-u6*c3-u7*c4+u5*c5-xy3*c6
        g(i6) = x4*c1+u4*c2+u2*c3+u3*c4+u1*c5+y4*c6
      end do
      return
    end subroutine jacobi_rotateints

    subroutine jacobi_a_aa(h1,g,dm1,dm2,a,indx,inz)
      implicit none
      real(wp),intent(in)     :: h1(:),g(:),dm1(:),dm2(:)
      real(wp), intent(inout) :: a(:)
      integer,intent(in)      :: indx(:),inz(:)
      integer :: i,i1,i2,i3,i4,i5,i6
      a=0.0_wp
      do i=1,inz(1),3
        i1 = indx(1)
        i2 = indx(2)
        i3 = indx(3)
        a(3)=a(3)+dm1(i1)*h1(i3)-dm1(i2)*h1(i2)+dm1(i3)*h1(i1)
        a(4)=a(4)+2.0_wp*h1(i2)*(dm1(i3)-dm1(i1))+dm1(i2)*(h1(i1)-h1(i3))
        a(5)=a(5)+dm1(i1)*h1(i1)+dm1(i2)*h1(i2)+dm1(i3)*h1(i3)
      end do
      do i=inz(1)+1,inz(2),2
         i1=indx(i)
         i2=indx(i+1)
         a(1)=a(1)+ dm1(i2)*h1(i1)-dm1(i1)*h1(i2)
         a(2)=a(2)+ dm1(i1)*h1(i1)+dm1(i2)*h1(i2)
      end do
      do i=inz(2)+1,inz(3),2
         i1=indx(i)
         i2=indx(i+1)
         a(1)=a(1)+ dm2(i2)*g(i1)-dm2(i1)*g(i2)
         a(2)=a(2)+ dm2(i1)*g(i1)+dm2(i2)*g(i2)
      end do
      do i=inz(3)+1,inz(4),4
         i1=indx(i)
         i2=indx(i+1)
         i3=indx(i+2)
         i4=indx(i+3)
         a(3)=a(3)+ dm2(i1)*g(i4)-dm2(i3)*g(i2)-dm2(i2)*g(i3)+&
                    dm2(i4)*g(i1)
         a(4)=a(4)+ (dm2(i4)-dm2(i1))*(g(i3)+g(i2))+(dm2(i3)+&
                    dm2(i2))*(g(i1)-g(i4))
         a(5)=a(5)+ dm2(i1)*g(i1)+dm2(i2)*g(i2)+dm2(i3)*g(i3)+&
                    dm2(i4)*g(i4)
      end do
      do i=inz(4)+1,inz(5),3
         i1=indx(i)
         i2=indx(i+1)
         i3=indx(i+2)
         a(3)=a(3)+ dm2(i1)*g(i3)-dm2(i2)*g(i2)+dm2(i3)*g(i1)
         a(4)=a(4)+ 2.0_wp*g(i2)*(dm2(i3)-dm2(i1))+&
                    dm2(i2)*(g(i1)-g(i3))
         a(5)=a(5)+ dm2(i1)*g(i1)+dm2(i2)*g(i2)+&
                    dm2(i3)*g(i3)
      end do
      do i=inz(5)+1,inz(6),3
         i1=indx(i)
         i2=indx(i+1)
         i3=indx(i+2)
         a(3)=a(3)+ dm2(i1)*g(i3)-dm2(i2)*g(i2)+&
                    dm2(i3)*g(i1)
         a(4)=a(4)+ 2.0_wp*g(i2)*(dm2(i3)-dm2(i1))+&
                    dm2(i2)*(g(i1)-g(i3))
         a(5)=a(5) + dm2(i1)*g(i1)+dm2(i2)*g(i2)+&
                    dm2(i3)*g(i3)
      end do
      do i=inz(6)+1,inz(7),6
         i1=indx(i)
         i2=indx(i+1)
         i3=indx(i+2)
         i4=indx(i+3)
         i5=indx(i+4)
         i6=indx(i+5)
         a(6)=a(6)- dm2(i1)*g(i6)+dm2(i2)*g(i5)-dm2(i3)*g(i4)+&
                    dm2(i4)*g(i3)-dm2(i5)*g(i2)+dm2(i6)*g(i1)
         a(7)=a(7)+ dm2(i1)*(g(i3)+2.0_wp*g(i5))+&
                    dm2(i2)*(g(i6)-g(i2)-g(i4))+&
                    dm2(i3)*(g(i1)-2.0_wp*g(i5))+&
                    dm2(i4)*(g(i6)-2.0_wp*g(i2))+&
                    dm2(i5)*(g(i1)-g(i3)-g(i5))+&
                    dm2(i6)*(g(i4)+2.0_wp*g(i2))
         a(8)=a(8)- dm2(i1)*(g(i4)+2.0_wp*g(i2))+&
                    dm2(i2)*(g(i1)-g(i3)-g(i5))-&
                    dm2(i3)*(g(i6)-2.0_wp*g(i2))+&
                    dm2(i4)*(g(i1)-2.0_wp*g(i5))-&
                    dm2(i5)*(g(i6)-g(i2)-g(i4))+&
                    dm2(i6)*(g(i3)+2.0_wp*g(i5))
         a(9)=a(9)+ dm2(i1)*g(i1)+dm2(i2)*g(i2)+dm2(i3)*g(i3)+&
                    dm2(i4)*g(i4)+dm2(i5)*g(i5)+dm2(i6)*g(i6)
      end do
      do i=inz(7)+1,inz(8),6
        i1=indx(i)
        i2=indx(i+1)
        i3=indx(i+2)
        i4=indx(i+3)
        i5=indx(i+4)
        i6=indx(i+5)
        a(10)= dm2(i1)*g(i6)-dm2(i2)*g(i5)+dm2(i3)*g(i3)+&
               dm2(i4)*g(i4)-dm2(i5)*g(i2)+dm2(i6)*g(i1)
        a(11)= dm2(i1)*(-4.0_wp*g(i5))+&
               dm2(i2)*(g(i3)+2.0_wp*g(i4)-g(i6))+&
               dm2(i3)*(-2.0_wp*g(i2)+2.0_wp*g(i5))+&
               dm2(i4)*(-2.0_wp*g(i2)+2.0_wp*g(i5))+&
               dm2(i5)*(g(i1)-2.0_wp*g(i4)-g(i3))+&
               dm2(i6)*(4.0_wp*g(i2))
        a(12)= dm2(i1)*(2.0_wp*g(i3)+4.0_wp*g(i4))+&
               dm2(i2)*(-3.0_wp*g(i2)+3.0_wp*g(i5))+&
               dm2(i3)*(g(i1)-4.0_wp*g(i4)+g(i6))+&
               dm2(i4)*(g(i1)+g(i6)-2.0_wp*g(i3)-&
               2.0_wp*g(i4))+&
               dm2(i5)*(3.0_wp*g(i2)-3.0_wp*g(i5))+&
               dm2(i6)*(2.0_wp*g(i3)+4.0_wp*g(i4))
        a(13)= dm2(i1)*(-4.0D+00*g(i2))+&
               dm2(i2)*(g(i1)-g(i3)-2.0D+00*g(i4))+&
               dm2(i3)*(2.0D+00*g(i2)-2.0D+00*g(i5))+&
               dm2(i4)*(2.0D+00*g(i2)-2.0D+00*g(i5))+&
               dm2(i5)*(g(i3)+2.0D+00*g(i4)-g(i6))+&
               dm2(i6)*(4.0D+00*g(i5))
        a(14)= dm2(i1)*g(i1)+dm2(i2)*g(i2)+dm2(i3)*g(i3)+&
               dm2(i4)*g(i4)+dm2(i5)*g(i5)+dm2(i6)*g(i6)
      end do
      return
    end subroutine jacobi_a_aa

    subroutine jacobi_get_bc_aa(a,b,c)
      implicit none
      real(wp), intent(in) :: a(:)
      real(wp), intent(out) :: b(:),c(:)      
      b(1) = a(5) + a(14)
      b(2) = a(1) + a(8)
      b(3) = a(3) - a(5) + a(12) - 2.0_wp*a(14)
      b(4) = a(6) - a(8)
      b(5) = a(10) - a(12) + a(14)
      b(6) = a(2) + a(9)
      b(7) = a(4) + a(13)
      b(8) = a(7) - a(9)
      b(9) = a(11) - a(13)
      c(1) = b(7)
      c(2) = 2.0_wp*b(8) - b(6)
      c(3) = 3.0_wp*b(9) - 2.0_wp*b(7)
      c(4) = -3.0_wp*b(8)
      c(5) = -4.0_wp*b(9)
      c(6) = b(2)
      c(7) = 2.0_wp*b(3)
      c(8) = 3.0_wp*b(4)
      c(9) = 4.0_wp*b(5)
      return
    end subroutine jacobi_get_bc_aa

    subroutine jacobi_a_av(h1,g,dm1,dm2,a,indx,inz)
      implicit none
      real(wp), intent(in)  :: h1(:),g(:),dm1(:),dm2(:)
      integer,  intent(in)  :: indx(:),inz(:)
      real(wp), intent(out) :: a(:)
      integer :: i,ii,in,i1,i2,i3,i4,i5,i6
      a=0.0_wp
      do i=1,inz(1),3
        i1 = indx(1)
        i2 = indx(2)
        i3 = indx(3)
        a(3) = a(3)+dm1(i1)*h1(i3)
        a(4) = a(4)-2.0_wp*dm1(i1)*h1(i2)
        a(5) = a(5)+dm1(i1)*h1(i1)
      end do
      do ii=inz(1)+1,inz(2),2
         i1 = indx(ii)
         i2 = indx(ii+1)
         a(1) = a(1) - dm1(i1)*h1(i2)
         a(2) = a(2) + dm1(i1)*h1(i1)
      end do
      do ii=inz(2)+1,inz(3),2
         i1 = indx(ii)
         i2 = indx(ii+1)
         a(1) = a(1) - dm2(i1)*g(i2)
         a(2) = a(2) + dm2(i1)*g(i1)
      end do
      do ii=inz(3)+1,inz(4),4
         i1 = indx(ii)
         i2 = indx(ii+1)
         i3 = indx(ii+2)
         i4 = indx(ii+3)
         a(3) = a(3) + dm2(i1)*g(i4)
         a(4) = a(4) + dm2(i1)*(-g(i2)-g(i3))
         a(5) = a(5) + dm2(i1)*g(i1)
      end do
      do ii=inz(4)+1,inz(5),3
         i1 = indx(ii)
         i2 = indx(ii+1)
         i3 = indx(ii+2)
         a(3) = a(3) + dm2(i1)*g(i3)
         a(4) = a(4) - 2.0_wp*dm2(i1)*g(i2)
         a(5) = a(5) + dm2(i1)*g(i1)
      end do
      do ii=inz(5)+1,inz(6),3
         i1 = indx(ii)
         i2 = indx(ii+1)
         i3 = indx(ii+2)
         a(3) = a(3) + dm2(i1)*g(i3)
         a(4) = a(4) -2.0_wp*dm2(i1)*g(i2)
         a(5) = a(5) + dm2(i1)*g(i1)
      end do
      do ii=inz(6)+1,inz(7),6
         i1 = indx(ii)
         i2 = indx(ii+1)
         i3 = indx(ii+2)
         i4 = indx(ii+3)
         i5 = indx(ii+4)
         i6 = indx(ii+5)
         a(6) = a(6) - dm2(i1)*g(i6)
         a(7) = a(7) + dm2(i1)*(g(i3)+2.0_wp*g(i5))
         a(8) = a(8) + dm2(i1)*(-2.0_wp*g(i2)-g(i4))
         a(9) = a(9) + dm2(i1)*g(i1)
      end do
      do i=inz(7)+1,inz(8),6
        i1 = indx(i)
        i2 = indx(i+1)
        i3 = indx(i+2)
        i4 = indx(i+3)
        i5 = indx(i+4)
        i6 = indx(i+5)
        a(10) = dm2(i1)*g(i6)
        a(11) = -4.0_wp*dm2(i1)*g(i5)
        a(12) = dm2(i1)*(2.0_wp*g(i3)+4.0_wp*g(i4))
        a(13) = -4.0_wp*dm2(i1)*g(i2)
        a(14) = dm2(i1)*g(i1)
      end do
      return
    end subroutine jacobi_a_av

    subroutine jacobi_get_bc_av(a,b,c)
      implicit none
      real(wp), intent(in) :: a(:)
      real(wp), intent(out) :: b(:),c(:)
      b=0.0_wp
      c=0.0_wp
      b(1) = a(5) + a(14)
      b(2) = a(1) + a(8)
      b(3) = a(3) - a(5) + a(12) - 2.0_wp*a(14)
      b(4) = a(6) - a(8)
      b(5) = a(10) - a(12) + a(14)
      b(6) = a(2) + a(9)
      b(7) = a(4) + a(13)
      b(8) = a(7) - a(9)
      b(9) = a(11) - a(13)
      c(1) = b(7)
      c(2) = 2.0_wp*b(8) - b(6)
      c(3) = 3.0_wp*b(9) - 2.0_wp*b(7)
      c(4) = -3.0_wp*b(8)
      c(5) = -4.0_wp*b(9)
      c(6) = b(2)
      c(7) = 2.0_wp*b(3)
      c(8) = 3.0_wp*b(4)
      c(9) = 4.0_wp*b(5)
      return
    end subroutine jacobi_get_bc_av

    subroutine jacobi_getpi(xpi,ypi)
      implicit none
      integer :: i
      real(wp) :: pi,gamma
      real(wp) :: xpi(:),ypi(:)
! this subroutine returns the values of sin(t) and cos(t) over the interval -pi <= t <=pi
! 64 increments
      pi=3.141592653589793238462643383279502884197169399D+00
      xpi(33) = 0.0_wp
      ypi(33) = 1.0_wp
      xpi(1) = 0.0_wp
      ypi(1) = -1.0_wp
      xpi(65) = 0.0_wp
      ypi(65) = -1.0_wp
      do i=1,31
        gamma = (dble(i)/32.0_wp)*pi
        xpi(i+33) = sin(gamma)
        ypi(i+33) = cos(gamma)
        xpi(33-i) = -xpi(i+33)
        ypi(33-i) = ypi(i+33)
      end do
      return
    end subroutine jacobi_getpi

    subroutine jacobi_bracketmin(ncor,i,a,b,c,xpi,ypi,ibv,xmin,pairde,found)
      implicit none
      logical :: found
      real(wp), intent(in)    :: a(:),b(:),c(:),xpi(:),ypi(:)
      real(wp), intent(inout) :: pairde
      integer, intent(in)     :: i,ncor
      integer, intent(inout)  :: xmin,ibv(4)

      integer :: ico,ii,ind,ismall,ping
      real(wp) :: dene,dene1,dene2,dene3

      found=.false.
      dene=jacobi_dedt(0.0_wp,1.0_wp,b)
      ico=0
      if (i.gt.ncor) then
        dene2=jacobi_dedt(xpi(1),ypi(1),c)
        do ii=5,65,4
          dene1=dene2
          dene2=jacobi_dedt(xpi(ii),ypi(ii),c)
          if (dene1.ge.0.0_wp) cycle
          if (dene2.ge.0.0_wp) then
            ico=ico+1
            ibv(ico)=ii-4
          end if
        end do
        if (ico.eq.0) return
      else
        dene2=jacobi_dedt(xpi(17),ypi(17),c)
        do ii=25,45,4
          dene1=dene2
          dene2=jacobi_dedt(xpi(ii),ypi(ii),c)
          if (dene1.ge.0.0_wp) cycle
          if (dene2.ge.0.0_wp) then
            ico=ico+1
            ibv(ico)=ii-4
          end if
        end do
        if (ico.eq.0) return
      end if
      do ii=1,ico
        ind=ibv(ii)+2
        dene3=jacobi_dedt(xpi(ind),ypi(ind),c)
        if (dene3.lt.0.0_wp) ibv(ii)=ind
      end do
      do ii=1,ico
        ind=ibv(ii)+1
        dene3=jacobi_dedt(xpi(ind),ypi(ind),c)
        if (dene3.lt.0.0_wp) ibv(ii)=ind
      end do
      ping=1
      ismall=iabs(ibv(1)-33)
      do ii=2,ico
        if (iabs(ibv(ii)-33).lt.ismall) ping=ii
      end do
      xmin=ibv(ping)
      pairde=a(ping)-dene
      found=.true.
      return
    end subroutine jacobi_bracketmin

    subroutine jacobi_refinemin(b,c,xs0,xs1,x2,y2,pairde)
      implicit none
      real(wp), intent(in) :: b(:),c(:),xs0,xs1
      real(wp), intent(inout) :: x2,y2,pairde
      integer, parameter :: mxit=50
      real(wp) :: e0,x0,x1,x02,x12,x22,y0,y1,fd0,fd1,fd2,f2
      integer :: isign,ii
      e0=jacobi_dedt(0.0_wp,1.0_wp,b)
      x0=xs0
      x1=xs0
      isign = +1
      if (xs0.lt.xs1) x1 = xs1
      if (xs0.Gt.xs1) then
         x0 = xs1
         isign = -1
      end if
      x02 = x0*x0 ! sin(t0)
      x12 = x1*x1 ! sin(t1)
      y0 = sqrt(abs(1.0_wp-x02)) ! cos(t0)
      y1 = sqrt(abs(1.0_wp-x12)) ! cos(t1)
      if (y0.lt.1.0E-16_wp.or.y1.lt.1.0E-16_wp) then
          x2 = xs0
          y2 = sqrt(abs(1.0_wp-xs0))
          pairde = jacobi_dedt(x2,y2,b) - e0
          return
      end if
      y0 = isign*y0
      y1 = isign*y1
      fd0 = jacobi_dedx(x0,y0,c) ! derivative of energy at t0
      fd1 = jacobi_dedx(x1,y1,c) ! derivative of energy at t1
      do ii=1,mxit
        x2 = x0 - ((x1-x0)*fd0)/(fd1 - fd0)
        x22 = x2*x2
        y2 = isign*sqrt(abs(1.0_wp-x22))
        fd2 = jacobi_dedx(x2,y2,c)
        if (abs(fd2).le.1.0E-16_wp)  then ! used to be 10E-16
          f2 = jacobi_dedt(x2,y2,b)
          pairde = f2 - e0 ! energy change with respect to starting angles
          if (pairde.ge.0.0_wp.or.abs(x2).lt.1.0E-16_wp) then
             x2 = 0.0_wp
             y2 = 1.0_wp
             pairde = 0.0_wp
          end if
          return
        end if
        if (fd2.lt.0.0_wp)  then
          fd0 = fd2
          x0 = x2
          y0 = y2
        else
          fd1 = fd2
          x1 = x2
          y1 = y2
        end if
      end do
      x2 = 0.0_wp
      y2 = 1.0_wp
      pairde = 0.0_wp
      return
    end subroutine jacobi_refinemin

    function jacobi_dedt(x,y,c)
! function to evaluate the derivative of the energy with respect to the rotation angle
! equation 3.32 in J. comp. chem. 24, 1250 (2003)
      implicit none
      real(wp), intent(in) :: x,y
      real(wp), intent(in) :: c(:)
      real(wp) :: jacobi_dedt
      real(wp) :: x2,x3,x4
      x2 = x*x
      x3 = x2*x
      x4 = x3*x
      jacobi_dedt = c(1)+c(2)*x+c(3)*x2+c(4)*x3+c(5)*x4+&
            y*(c(6)+c(7)*x+c(8)*x2+c(9)*x3)
      return
    end function jacobi_dedt

    function jacobi_dedx(x,y,c)
! function to evaluate the derivative of the energy with respect to the x=sin(t)
! equation 3.33 in J. comp. chem. 24, 1250 (2003)
      implicit none
      real(wp), intent(in) :: x,y
      real(wp), intent(in) :: c(:)
      real(wp) :: jacobi_dedx
      real(wp) :: x2,x3,x4
      x2 = x*x
      x3 = x2*x
      x4 = x3*x
      jacobi_dedx = (c(1)+c(2)*x+c(3)*x2+c(4)*x3+c(5)*x4)/y +&
            c(6)+c(7)*x+c(8)*x2+c(9)*x3
      return
    end function jacobi_dedx

    pure function index(i,j)
! this function computes the two-electron index index (lower triangular reference)
! index = ii*(ii-1)/2+jj where ii=max(i,j) and jj=min(i,j)
! the ishft(k,-1) divides the value of the integer k by 2 and seems to be somewhat
! faster than the regular human-readable expression
      implicit none
      integer, intent(in) ::i,j
      integer :: index
      if (i.ge.j) then
        index=ishft(i*(i-1),-1)+j
        return
      else
        index=ishft(j*(j-1),-1)+i
        return
      end if
    end function index

end module jacobi_mod
