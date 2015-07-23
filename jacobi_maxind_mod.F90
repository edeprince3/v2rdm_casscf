module jacobi_maxind_mod
  use jacobi_data, only : symmetries,rotinds,group_mult_tab,ncore,nact,nvirt,ntot,npair,maxind_nnz

  implicit none

  contains

    subroutine jacobi_maxind()
      implicit none
      integer :: indz(8)
      integer :: i,j,iang
      maxind_nnz=0
      do iang=0,2*npair-2,2  ! loop over rotation pairs
        i=rotinds(iang+1)    ! orbital i
        j=rotinds(iang+2)    ! orbital j (j>i)
        if (j.le.nact+ncore) then
          call jacobi_maxinds_aa(i,j,ncore,nact,ntot,indz)
          maxind_nnz=max(maxind_nnz,maxval(indz))
        else
          call jacobi_maxinds_av(i,j,ncore,nact,ntot,indz)
          maxind_nnz=max(maxind_nnz,maxval(indz))
        end if
      end do
    end subroutine jacobi_maxind

    subroutine jacobi_maxinds_av(i,j,rc,ra,rt,indz)
      implicit none
      integer, intent(in) :: rc,ra,rt,i,j
      integer, intent(out):: indz(:)
      integer :: p,q,r,iis,ijs,ips,r1,is1,is2,na,nt
      iis=symmetries(i)
      ijs=symmetries(j)
      indz=0
      na=rc+ra
      nt=rt
      r1=group_mult_tab(iis,ijs)
      indz=0
      indz(1) = 3
      do p=1,i-1
         if (iis.ne.symmetries(p)) cycle
         indz(1) = indz(1) + 2
      end do
      do p=i+1,na
        if (iis.ne.symmetries(p)) cycle
         indz(1) = indz(1) + 2
      end do
      do r=1,i-1
        is1=group_mult_tab(iis,symmetries(r))
        do p=1,i-1
          ips=symmetries(p)
          do q=1,p
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,na
        is1=group_mult_tab(iis,symmetries(r))
        do p=1,r-1
          if (p.eq.i) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is2.ne.is1) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=i+1,na
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,na
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,r
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,na
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is2.ne.is1) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do p=1,i-1
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=1,p-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=i+1,na
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=i+1,na
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=i+1,p-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=1,i-1
        is1=group_mult_tab(symmetries(p),iis)
        indz(is1) = indz(is1) + 3
      end do
      do p=i+1,na
        is1=group_mult_tab(symmetries(p),iis)
        indz(is1) = indz(is1) + 3
      end do
      do p=1,i-1
        ips=symmetries(p)
        do q=1,p
          if (ips.ne.symmetries(q)) cycle
          indz(1) = indz(1) + 3
        end do
      end do
      do p=i+1,na
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          if (ips.ne.symmetries(q)) cycle
          indz(1) = indz(1) + 3
        end do
      end do
      do p=1,i-1
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 6
      end do
      do p=i+1,na
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 6
      end do
      indz(1) = indz(1) + 6
      do p=na+1,j-1
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 2
      end do
      do p=j+1,nt
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 2
      end do
      do r=na+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        do p=1,r-1
          if (p.eq.i) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do r=j+1,nt
        is1=group_mult_tab(iis,symmetries(r))
        do p=1,r-1
          if (p.eq.i.or.p.eq.j) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=na+1,j-1
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        do p=j+1,nt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,nt
        is1=group_mult_tab(iis,symmetries(r))
        do q=j+1,r
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,nt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,nt
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,j-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do r=i+1,na
        is1=group_mult_tab(iis,symmetries(r))
        do p=na+1,j-1
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=na+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,r
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,j-1
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do p=na+1,j-1
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=na+1,j-1
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=i+1,p-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=i+1,j-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=j+1,p-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=na+1,j-1
        is1=group_mult_tab(iis,symmetries(p))
        indz(is1) = indz(is1) + 3
      end do
      do p=j+1,nt
        is1=group_mult_tab(iis,symmetries(p))
        indz(is1) = indz(is1) + 3
      end do
      do p=na+1,j-1
        ips=symmetries(p)
        do q=1,p
          IF (q.eq.i) cycle
          if (ips.ne.symmetries(q)) cycle
          indz(1) = indz(1) + 3
        end do
      end do
      do p=j+1,nt
        ips=symmetries(p)
        do q=1,p
          IF (q.eq.i.or.q.eq.j) cycle
          if (ips.ne.symmetries(q)) cycle
          indz(1) = indz(1) + 3
        end do
      end do
      do p=na+1,j-1
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 6
      end do
      do p=j+1,nt
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 6
      end do
      return
    end subroutine jacobi_maxinds_av

    subroutine jacobi_maxinds_aa(i,j,rc,ra,rt,indz)
      implicit none
      integer, intent(in) :: rc,ra,rt
      integer, intent(out) :: indz(:)
      integer :: i,j
      integer :: p,q,r,iis,ijs,ips,r1,is1,is2,na
      iis=symmetries(i)
      ijs=symmetries(j)
      indz=0
      na=rc+ra
      r1=group_mult_tab(iis,ijs)
      indz(1)=3
      do p=1,i-1
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 2
      end do
      do p=i+1,j-1
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 2
      end do
      do p=j+1,na
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 2
      end do
      do r=1,i-1
        is1=group_mult_tab(iis,symmetries(r))
        do p=1,i-1
          ips=symmetries(p)
          do q=1,p
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        do p=1,r-1
          if (p.eq.i) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do r=j+1,na
        is1=group_mult_tab(iis,symmetries(r))
        do p=1,r-1
          if (p.eq.i.or.p.eq.j) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=i+1,j-1
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        do p=j+1,na
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,na
        is1=group_mult_tab(iis,symmetries(r))
        do q=j+1,r
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,na
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,na
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,j-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,r
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,j-1
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is2.ne.is1) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do p=1,i-1
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        do q=1,p-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=i+1,j-1
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        do q=1,I-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis) 
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=i+1,j-1
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        do q=i+1,p-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=i+1,j-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=j+1,p-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=1,i-1
        is1=group_mult_tab(symmetries(p),iis)
        indz(is1) = indz(is1) + 3
      end do
      do p=i+1,j-1
        is1=group_mult_tab(symmetries(p),iis)
        indz(is1) = indz(is1) + 3
      end do
      do p=j+1,na
        is1=group_mult_tab(symmetries(p),iis)
        indz(is1) = indz(is1) + 3
      end do
      do p=1,i-1
        ips=symmetries(p)
        do q=1,p
          if (ips.ne.symmetries(q)) cycle
          indz(1) = indz(1) + 3
        end do
      end do
      do p=i+1,j-1
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i) cycle
          if (ips.ne.symmetries(q)) cycle
          indz(1) = indz(1) + 3
        end do
      end do
      do p=j+1,na
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          if (ips.ne.symmetries(q)) cycle
          indz(1) = indz(1) + 3
        end do
      end do
      do p=1,i-1
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 6
      end do
      do p=i+1,j-1
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 6
      end do
      do p=j+1,na
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 6
      end do
      indz(1) = indz(1) + 6
      do p=na+1,rt
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 2
      end do
      do r=na+1,rt
        is1=group_mult_tab(iis,symmetries(r))
        do p=1,r-1
          if (p.eq.i.or.p.eq.j) cycle
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
        do q=1,i-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=na+1,rt
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          is2=group_mult_tab(ips,symmetries(q))
          do r=1,i-1
            is1=group_mult_tab(iis,symmetries(r))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=i+1,j-1
        is1=group_mult_tab(iis,symmetries(r))
        do p=na+1,rt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=j+1,na
        is1=group_mult_tab(iis,symmetries(r))
        do p=na+1,rt
          ips=symmetries(p)
          do q=1,p
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
           end do
        end do
      end do
      do r=na+1,rt
        is1=group_mult_tab(iis,symmetries(r))
        do q=j+1,r
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
        do p=r+1,rt
          ips=symmetries(p)
          do q=1,P
            if (q.eq.i.or.q.eq.j) cycle
            is2=group_mult_tab(ips,symmetries(q))
            if (is1.ne.is2) cycle
            indz(is1) = indz(is1) + 2
          end do
        end do
      end do
      do r=na+1,rt
        is1=group_mult_tab(iis,symmetries(r))
        do q=i+1,j-1
          if (iis.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 2
        end do
      end do
      do p=na+1,rt
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        do q=1,i-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=na+1,rt
        ips=symmetries(p)
        is1=group_mult_tab(ips,iis)
        do q=i+1,j-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=na+1,rt
        ips=symmetries(p)
        is1=group_mult_tab(iis,ips)
        do q=j+1,p-1
          if (ips.ne.symmetries(q)) cycle
          indz(is1) = indz(is1) + 4
        end do
      end do
      do p=na+1,rt
        is1=group_mult_tab(iis,symmetries(p))
        indz(is1) = indz(is1) + 3
      end do
      do p=na+1,rt
        ips=symmetries(p)
        do q=1,p
          if (q.eq.i.or.q.eq.j) cycle
          if (ips.ne.symmetries(q)) cycle
          indz(1) = indz(1) + 3
        end do
      end do
      do p=na+1,rt
        if (iis.ne.symmetries(p)) cycle
        indz(1) = indz(1) + 6
      end do
      return
    end subroutine jacobi_maxinds_aa

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

end module jacobi_maxind_mod
