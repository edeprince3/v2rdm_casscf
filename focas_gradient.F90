module focas_gradient
  use focas_data
 
  implicit none

  contains

  subroutine orbital_gradient(int1,int2,den1,den2)
    implicit none
    real(wp), intent(in) :: int1(:),int2(:),den1(:),den2(:)
    real(wp) :: t0,t1

    ! calculate inactive Fock matrix
    t0=timer()
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_f_i(int1,int2)
    else
      call compute_f_i_df(int1,int2)
    endif
    t1=timer()
!    write(*,*)'f_i',t1-t0

    ! calculate active Fock matrix
    t0=timer()
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_f_a(den1,int2)
    else
      call compute_f_a_df(den1,int2)
    endif
    t1=timer()
!    write(*,*)'f_a',t1-t0


    ! calculate auxiliary q matrix
    t0=timer()
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_q(den2,int2)
    else
      call compute_q_df(den2,int2)
    endif
    t1=timer()
!    write(*,*)'q',t1-t0

 
    ! calculate auxiliary z matrix
    t0=timer()
    call compute_z(den1)
    t1=timer()
!    write(*,*)'z',t1-t0 

  !  stop

    ! compute gradient
    call compute_orbital_gradient()

!    ! print the gradient 
!    if ( log_print_ == 1 ) call print_orbital_gradient()

    return
  end subroutine orbital_gradient

  subroutine print_orbital_gradient()
    implicit none
    integer :: ij_pair,i,j,newline,i_class,j_class,j_class_start,j_start,i_sym
    character(1) :: ityp,jtyp
    ! loop over rotation pairs

    write(fid_,'(a)')'orbital gradient:'

    newline=0

    ! the loop structure below cycles through the possible rotation pairs
    ! rotation pairs are sorted according to symmetry and for each irrep,
    ! the rotation pairs are sorted according to orbital classes: ad,ed,aa,ea
    ! for each pair j>i; for more details, see subroutine setup_rotation_indeces in focas_main.F90

    ij_pair = 0

    do i_sym = 1 , nirrep_

      do i_class = 1 , 3

        j_class_start = i_class + 1

        if ( ( include_aa_rot_ == 1 ) .and. ( i_class == 2 ) ) j_class_start = i_class

        do j_class = j_class_start , 3

          do i = first_index_(i_sym,i_class) , last_index_(i_sym,i_class)

            j_start = first_index_(i_sym,j_class)

            if ( i_class == j_class ) j_start = i + 1

            do j = j_start , last_index_(i_sym,j_class)

              if ( i_class == 1 ) then

                ityp = 'd'

                jtyp = 'a'

                if ( j_class == 3 ) jtyp = 'e'

              else

                ityp = 'a'

                jtyp = 'a'

                if ( j_class == 3 ) jtyp = 'e'

              end if

              ij_pair = ij_pair + 1

              write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')jtyp,'-',ityp,' (',j,',',i,')',orbital_gradient_(ij_pair)

              newline = newline + 1

              if ( mod(newline,4) == 0 ) write(fid_,*)

            end do

          end do
 
        end do
 
      end do 

    end do
    if ( mod(newline,4) /= 0 ) write(fid_,*)

    write(fid_,'(a,1x,es10.3)')'gradient norm:',grad_norm_

    return

  end subroutine print_orbital_gradient

!  subroutine compute_orbital_gradient_old()
!    implicit none
!    ! subroutine to compute the gradient vector according to Eq. 12.5.5 in Helgaker on page 622
!    ! orbital_gradient(i|j) = 2 ( f_gen(i,j) - f_gen(j,i) )
!    integer :: i,j,ij_pair,j_start,i_sym,i_class,j_class,j_class_start
!
!    ! the loop structure below cycles through the possible rotation pairs
!    ! rotation pairs are sorted according to symmetry and for each irrep,
!    ! the rotation pairs are sorted according to orbital classes: ad,ed,aa,ea
!    ! for each pair j>i; for more details, see subroutine setup_rotation_indeces in focas_main.F90
!
!    ij_pair = 0
!
!    do i_sym = 1 , nirrep_
!
!      do i_class = 1 , 3
!
!        j_class_start = i_class + 1
! 
!        if ( ( include_aa_rot_ == 1 ) .and. ( i_class == 2 ) ) j_class_start = i_class
!
!        do j_class = j_class_start , 3
!
!          do i = first_index_(i_sym,i_class) , last_index_(i_sym,i_class)
!
!            j_start = first_index_(i_sym,j_class)
!
!            if ( i_class == j_class ) j_start = i + 1 
!
!            do j = j_start , last_index_(i_sym,j_class)
!
!              ij_pair = ij_pair + 1
!
!              orbital_gradient_(ij_pair) = 2.0_wp * ( fock_gen_(i,j) - fock_gen_(j,i) )
!          
!            end do ! end j loop
!          
!          end do ! end i loop
! 
!        end do ! end j_class loop
!
!      end do ! end i_class loop
!
!    end do ! end i_sym loop 
!
!    ! calculate orbital gradient norm
!    grad_norm_ = dot_product(orbital_gradient_,orbital_gradient_) 
!
!    return
!  end subroutine compute_orbital_gradient_old

  subroutine compute_orbital_gradient()

    ! subroutine to compute the orbital gradient without explicit storage
    ! of the generalized Fock matrix

    implicit none

    integer :: grad_ind
    integer :: a_sym,t_sym
    integer :: i,t,u,a
    integer :: ia,it
    real(wp) :: ddot

    orbital_gradient_ = 0.0_wp

    ! ********************************
    ! external - doubly-occupied pairs
    ! ********************************

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)
  
          grad_ind = grad_ind + 1

          ia = ints_%gemind(i,a)

          orbital_gradient_(grad_ind) = 4.0_wp*(fock_i_(ia)+fock_a_(ia))

        end do

      end do
 
    end do  

    ! ******************************
    ! acitve - doubly-occupied pairs
    ! ******************************

    do t_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

      do i = first_index_(t_sym,1) , last_index_(t_sym,1)

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

          grad_ind = grad_ind + 1

          it = ints_%gemind(i,t)

          orbital_gradient_(grad_ind) = 4.0_wp*(fock_i_(it)+fock_a_(it)) - 2.0_wp * ( q_(t - ndoc_tot_,i) + z_(t - ndoc_tot_,i))

        end do

      end do

    end do

    ! ***********************
    ! external - active pairs
    ! ***********************

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_act_type)

      do t = first_index_(a_sym,2) , last_index_(a_sym,2)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          orbital_gradient_(grad_ind) = 2.0_wp * ( q_(t - ndoc_tot_,a) + z_(t - ndoc_tot_,a))

        end do
    
      end do

    end do

    ! *********************
    ! active - active pairs
    ! *********************

    if ( include_aa_rot_ == 1 ) then

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1

            orbital_gradient_(grad_ind) = 2.0_wp * ( q_(u - ndoc_tot_,t) + &
                 & z_(u - ndoc_tot_,t) - q_(t - ndoc_tot_,u) - z_(t - ndoc_tot_,u))

          end do
 
        end do

      end do

    end if

    grad_norm_ = sqrt(ddot(rot_pair_%n_tot,orbital_gradient_,1,orbital_gradient_,1))

    return 

  end subroutine compute_orbital_gradient

!  subroutine compute_f_gen()
!    implicit none
!    ! subroutine to compute the generalized fock matrix accoring to Eqs. 12.5.-12.9.10 in Helgaker on page 622
!    ! f_gen(i|n) = 2 [ f_i(i|n) + f_a(i|n) ]               --> i \in D && n \in D,A,E
!    ! f_gen(v|n) = Q(v|n) + z(v|n)                         --> v \in A && w \in A && n \in D,A,E
!    ! f_gen(e,n) = 0                                       --> e \in E && n \in D,A,E
!    ! here :: z(v|n) = SUM(u) { f_i(v,u) * d1(v,u) }       --> v,u \in A && n \in D,A,E
!    !      :: Q(v|n) = SUM(xyz) { d(vxyz) * g(nxyz) }      --> v,x,y,z \in A && n \in D,A,E
!    integer :: i,v,n,ni,n_class,n_sym
!
!    ! initialize
!    fock_gen_ = 0.0_wp
!
!    ! loop over irreps for n
!
!    do n_sym = 1 , nirrep_
!
!      ! loop over class for n
!
!      do n_class = 1 , 3
!
!        ! loop over n indeces
!
!        do n = first_index_(n_sym,n_class) , last_index_(n_sym,n_class)
!
!          ! loop over i \in D
!
!          do i = first_index_(n_sym,1) , last_index_(n_sym,1)
!
!             ! ni geminal index
!             ni             = ints_%gemind(i,n)
!             ! update F(i|n) = 2 * [ f_a(i|n) + f_i(i|n) ]
!             fock_gen_(i,n) = 2.0_wp * ( fock_a_(ni) + fock_i_(ni) )
!
!          end do ! end i loop
!
!          ! loop over v in \A
!
!          do v = first_index_(n_sym,2) , last_index_(n_sym,2)
!
!            ! initialize F(v|n) = q_(v|n)
!            fock_gen_(v,n) = q_( v - ndoc_tot_ , n ) + z_( v - ndoc_tot_ , n )
!
!          end do ! end v loop
!
!        end do ! end n loop
!
!      end do ! end n_class loop
!
!    end do ! end n_sym loop
!
!    return 
!  end subroutine compute_f_gen

  subroutine compute_z(den1)
    implicit none
    ! function to compute contraction of the density with the inactive Fock matrix according to
    ! z(m,t) = sum_u { den1(tu) * f_i(mu) } t,u \in A m \in D,A,E
    real(wp), intent(in) :: den1(:)
    integer :: t,u,tu_den,mu_int,t_sym,m,m_class
    real(wp) :: val

    ! initialize output

    z_ = 0.0_wp

    ! loop over symmetries for t

    do t_sym = 1 , nirrep_

      ! loop over u indeces

      do t = first_index_(t_sym,2) , last_index_(t_sym,2)

        ! loop over classes for m

        do m_class = 1 , 3

          ! loop over m indeces
 
          do m = first_index_(t_sym,m_class) , last_index_(t_sym,m_class)

            ! initialize contraction value

            val = 0.0_wp

            do u = first_index_(t_sym,2) , last_index_(t_sym,2)

              ! integral/density addressing

              tu_den = dens_%gemind(t,u)
              mu_int = ints_%gemind(m,u)

              val = val + den1(tu_den) * fock_i_(mu_int)

            end do ! end u loop

            z_( t - ndoc_tot_ , m ) = val 

          end do ! end m loop

        end do ! end m_class loop

      end do ! end t loop

    end do ! end t_sym loop

    return

  end subroutine compute_z

  subroutine compute_q_df(den2,int2)
    implicit none
    ! this subroutine computes the "auxiliary Q" matrix according to Eq. 12.5.14 in Helgaker on page 622
    ! Q(v,m) = \SUM[w,x,y \in A] { d2(vw|xy) * g(mw|xy) }
    ! Because only those D2 elements with all four indeces \in A are nonzero, we have v \in A
    ! Also, m is a general orbital index && sym_m == sym_v
    real(wp), intent(in) :: den2(:),int2(:)
    integer :: mw_sym,x_sym,y_sym,w_sym,m_sym,m_class,m,v,w,x,y,mdf,wdf,xdf,ydf
    integer :: xy_den,vw_den
    integer(ip) :: mw,xy,nQ
    integer :: den_ind,den_sym_offset,int_col,int_row,y_last
    integer :: imw,ixy
    integer :: max_nmopi,max_ngempi,error
    real(wp) :: val,int_val,ddot
    real(wp) :: v_mw(df_vars_%nQ)

    type sym_block
      real(wp), allocatable :: mw_df(:,:)
      real(wp), allocatable :: xy_df(:,:)
      real(wp), allocatable :: ints(:,:)
    end type sym_block

    type(sym_block) :: gaaa(nirrep_)

    nQ = int(df_vars_%nQ,kind=ip)

    max_nmopi  = maxval(nactpi_)
    max_ngempi = max_nmopi * max_nmopi 

    error = allocate_gaaa()

    ! initialize
    q_ = 0.0_wp

    ! loop over symmetries for m

    do m_sym = 1 , nirrep_

      ! loop over classes for m

      do m_class = 1 , 3

        ! loop over orbital index m

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! m index in df order

          mdf = df_vars_%class_to_df_map(m)

          ! *************************************
          ! gather df integrals (:|mw) and (:|xy)
          ! *************************************

          do w_sym = 1 , nirrep_

            if ( nactpi_(w_sym) == 0 ) cycle
    
            imw = 0    

            mw_sym = group_mult_tab_(m_sym,w_sym)

            do w = first_index_(w_sym,2) , last_index_(w_sym,2)
 
              wdf = df_vars_%class_to_df_map(w)
 
              mw  = df_pq_index(wdf,mdf)

              imw = imw + 1

              call dcopy(nQ,int2(mw+1:mw+nQ),1,gaaa(w_sym)%mw_df(:,imw),1)

            end do ! end w loop

            ixy = 0

            do x_sym = 1 , nirrep_

              y_sym = group_mult_tab_(x_sym,mw_sym)
 
              if ( y_sym > x_sym ) cycle

              do x = first_index_(x_sym,2) , last_index_(x_sym,2) 

                xdf = df_vars_%class_to_df_map(x)

                y_last = last_index_(y_sym,2)
                if ( y_sym == x_sym ) y_last = x

                do y = first_index_(y_sym,2) , y_last

                  ydf     = df_vars_%class_to_df_map(y)

                  xy      = df_pq_index(xdf,ydf)

                  ixy     = ixy + 1

                  call dcopy(df_vars_%nQ,int2(xy+1:xy+nQ),1,gaaa(w_sym)%xy_df(:,ixy),1)

                end do ! end y loop  

              end do ! end x loop

            end do ! end x_sym loop

            if ( ixy == 0 ) cycle

            ! *******************************************************************************************
            ! multiply 3-index integral matrices to get the 4-index integrals (xy|mw) for all w,x,y \in A
            ! *******************************************************************************************

            call dgemm('t','n',ixy,imw,df_vars_%nQ,         &
                     & 1.0_wp,gaaa(w_sym)%xy_df(:,1:ixy),df_vars_%nQ,&
                     &        gaaa(w_sym)%mw_df(:,1:imw),df_vars_%nQ,&
                     & 0.0_wp,gaaa(w_sym)%ints,nirrep_*max_ngempi)

          end do ! end w_sym loop

          ! loop over orbital index v

#ifdef OMP
!$omp parallel shared(m_class,nQ,m_sym,m,mdf,first_index_,last_index_,df_vars_,den2,int2,dens_,q_) num_threads(nthread_use_)
!$omp do private(ixy,imw,val,v,w_sym,mw_sym,den_sym_offset,w,vw_den,x_sym,y_sym,x,y,xy_den,den_ind,int_val)
#endif

          do v = first_index_(m_sym,2) , last_index_(m_sym,2)

            ! initialize matrix element
  
            val = 0.0_wp

            ! loop over symmetries for w

            do w_sym = 1 , nirrep_

              imw = 0

              ! symmetry of mw geminal

              mw_sym = group_mult_tab_(m_sym,w_sym)

              ! vw density offset 
             
              den_sym_offset = dens_%offset(mw_sym)

              ! loop over w orbital index
 
              do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                ixy = 0   

                imw = imw + 1

                ! vw geminal index for density

                vw_den = dens_%gemind(v,w)

                ! loop over x_sym

                do x_sym = 1 , nirrep_

                  ! figure out symmetry of y

                  y_sym = group_mult_tab_(mw_sym,x_sym)

                  if ( y_sym == x_sym ) then
                 
                    ! loop over x orbital index
 
                    do x = first_index_(x_sym,2) , last_index_(x_sym,2)

                      ! ***********************************************
                      ! loop over y orbital index (x>y) --> factor of 2
                      ! ***********************************************

                      do y = first_index_(y_sym,2) , x - 1

                        ixy = ixy + 1

                        ! xy geminal index and density index

                        xy_den  = dens_%gemind(x,y)
                        den_ind = pq_index(vw_den,xy_den) + den_sym_offset
                    
                        ! xy geminal index in df order
                        xy  = df_pq_index(xdf,ydf)   

                        ! update matrix element

                        val = val + 2.0_wp * gaaa(w_sym)%ints(ixy,imw) * den2(den_ind)

                      end do ! end y loop

                      ! **********************
                      ! x == y --> factor of 1
                      ! **********************

                      xy_den  = dens_%gemind(x,x)
                      den_ind = pq_index(vw_den,xy_den) + den_sym_offset

                      ixy = ixy + 1

                      ! xy geminal index in df order
                      xy  = df_pq_index(xdf,xdf) 

                      ! update matrix element

                      val = val + gaaa(w_sym)%ints(ixy,imw) * den2(den_ind)

                    end do ! end x loop

                  elseif ( y_sym < x_sym) then

                    ! loop over x orbital index

                    do x = first_index_(x_sym,2) , last_index_(x_sym,2)


                      ! ***********************************************
                      ! loop over y orbital index (x>y) --> factor of 2
                      ! ***********************************************

                      do y = first_index_(y_sym,2) , last_index_(y_sym,2)

                        ! xy geminal index and density index

                        xy_den  = dens_%gemind(x,y)
                        den_ind = pq_index(vw_den,xy_den) + den_sym_offset

                        ! xy geminal index in df order
                        xy  = df_pq_index(xdf,ydf) 

                        ixy = ixy + 1

                        ! update matrix element

                        val = val + 2.0_wp * gaaa(w_sym)%ints(ixy,imw) * den2(den_ind)

                      end do ! end y loop

                    end do ! end x loop

                  endif 

                end do ! end x_sym loop

              end do ! end w loop

            end do ! end w_sym loop

            ! save matrix element

            q_( v - ndoc_tot_ , m ) = val
           
          end do ! end v loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

        end do ! end m loop

      end do ! end m_class loop 

    end do ! end m_sym loop
 
    error = deallocate_gaaa()

    return

    contains

      integer function allocate_gaaa()

        implicit none

        integer :: w_sym

        do w_sym = 1 , nirrep_

          if ( nactpi_(w_sym) == 0 ) cycle
 
          allocate(gaaa(w_sym)%mw_df(df_vars_%nQ,nactpi_(w_sym)))
          allocate(gaaa(w_sym)%xy_df(df_vars_%nQ,nirrep_*max_ngempi))
          allocate(gaaa(w_sym)%ints(nirrep_*max_ngempi,max_nmopi)) 

        end do
 
        allocate_gaaa = 0

      end function allocate_gaaa

      integer function deallocate_gaaa

        integer :: w_sym

        do w_sym = 1 , nirrep_       

          if (allocated(gaaa(w_sym)%mw_df)) deallocate(gaaa(w_sym)%mw_df)
          if (allocated(gaaa(w_sym)%xy_df)) deallocate(gaaa(w_sym)%xy_df)
          if (allocated(gaaa(w_sym)%ints))  deallocate(gaaa(w_sym)%ints)

        end do

        deallocate_gaaa = 0

      end function deallocate_gaaa 

  end subroutine compute_q_df

  subroutine compute_q(den2,int2)
    implicit none
    ! this subroutine computes the "auxiliary Q" matrix according to Eq. 12.5.14 in Helgaker on page 622
    ! Q(v,m) = \SUM[w,x,y \in A] { d2(vw|xy) * g(mw|xy) }
    ! Because only those D2 elements with all four indeces \in A are nonzero, we have v \in A
    ! Also, m is a general orbital index && sym_m == sym_v
    real(wp), intent(in) :: den2(:),int2(:)   
    integer :: mw_sym,x_sym,y_sym,w_sym,m_sym,m_class,m,v,w,x,y
    integer :: xy_den,xy_int,mw,vw
    integer :: den_ind,den_sym_offset,int_ind,int_sym_offset
    real(wp) :: val

    ! initialize
    q_ = 0.0_wp

    ! loop irreps for m and v

    do m_sym = 1 , nirrep_
 
      ! loop over irreps for w

      do w_sym = 1, nirrep_

        ! mw//vw symmetry
        mw_sym = group_mult_tab_(m_sym,w_sym)
        ! offsets for integral/density addressing
        den_sym_offset = dens_%offset(mw_sym)
        int_sym_offset = ints_%offset(mw_sym)

        ! loop over irreps for x

        do x_sym = 1 , nirrep_

          ! correspoing irrep for y
          y_sym = group_mult_tab_(x_sym,mw_sym)

          ! at this point, we have mw_sym == xy_sym && m_sym == v_sym
          ! loop over m_class

          do m_class = 1 , 3

            ! loop over v indeces
   
            do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

              ! loop over v \in A

              do v = first_index_(m_sym,2) , last_index_(m_sym,2)

                ! initialize q matrix element
                val     = 0.0_wp

                ! loop over w indeces

                do w = first_index_(w_sym,2) , last_index_(w_sym,2)
                  
                  ! save geminal indeces for integral/density addressing
                  mw = ints_%gemind(m,w)
                  vw = dens_%gemind(v,w)

                  ! loop over x indeces
 
                  do x = first_index_(x_sym,2) , last_index_(x_sym,2)

                    ! loop over y indeces

                    do y = first_index_(y_sym,2) , last_index_(y_sym,2) 

                      ! save geminal indeces for integral/density addressing
                      xy_int  = ints_%gemind(x,y)
                      xy_den  = dens_%gemind(x,y)
                      ! integral/density addresses
                      int_ind = pq_index(xy_int,mw) + int_sym_offset                      
                      den_ind = pq_index(xy_den,vw) + den_sym_offset
                      ! update temporary value
                      val     = val + den2(den_ind) * int2(int_ind)
 
                    end do ! end y loop

                  end do ! end x loop
                   
                end do ! end w loop

                ! update q matrix element
                q_( v - ndoc_tot_ , m ) = q_( v - ndoc_tot_ , m ) + val

              end do ! end v loop

            end do ! end m loop 

          end do ! end m_class loop

        end do ! end x_sym loop

      end do ! end w_sym loop

    end do ! end m_sym loop
 

    return
  end subroutine compute_q

  subroutine compute_f_a_df(den1,int2)
    implicit none
    ! this subroutine computes the "active fock matrix" according to Eq. 12.5.13 in Helgaker on page 622
    ! F_a(m|n) = h(m|n) + SUM[v,w \in A] { g(mn|vw) - 0.5 g(mw|vn) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes.
    ! v and w belong to the active orbital class
    real(wp), intent(in) :: den1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,v,w,w_sym,mdf,ndf,wdf,vdf
    integer :: den_ind,mn_fock
    integer(ip) :: mn,vw,mw,vn,wn,ww,nQ 
    real(wp) :: val,ival,dval,ddot
    real(wp) :: v_mn(df_vars_%nQ),v_mw(df_vars_%nQ)

    nQ = int(df_vars_%nQ,kind=ip)

    ! initialize
    fock_a_ = 0.0_wp

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! orbital index in df order
          mdf = df_vars_%class_to_df_map(m)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

#ifdef OMP
!$omp parallel shared(ints_,nQ,dens_,m_class,m_sym,m,mdf,first_index_,last_index_,df_vars_,den1,int2) num_threads(nthread_use_)
!$omp do private(n,ndf,mn,val,v_mn,v_mw,w_sym,w,wdf,v,den_ind,dval,vdf,vw,ival,ww,mw,wn,vn,mn_fock)
#endif

          do n = first_index_(m_sym,m_class) , m

            ! orbital index in df order
            ndf = df_vars_%class_to_df_map(n)

            ! mn-geminal index in df order
            mn  = df_pq_index(mdf,ndf) 

            call dcopy(df_vars_%nQ,int2(mn+1:mn+nQ),1,v_mn,1)

            ! initialize Fock matrix element
            val = 0.0_wp

            ! loop over irreps for v

            do w_sym = 1 , nirrep_

              ! loop over w indeces

              do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                ! orbital index in df order
                wdf = df_vars_%class_to_df_map(w)

                mw      = df_pq_index(mdf,wdf)

                call dcopy(df_vars_%nQ,int2(mw+1:mw+nQ),1,v_mw,1)

                ! *** w > v --> factor of 2

                ! loop over v indeces

                do v = first_index_(w_sym,2) , w - 1 

                  ! 1-e density element
                  den_ind = dens_%gemind(v,w)
                  dval    = den1(den_ind)

                  ! orbital index in df order
                  vdf     = df_vars_%class_to_df_map(v)

                  ! 2-e coulomb contribution 2 g(mn|vw)
                  vw      = df_pq_index(wdf,vdf) 

                  ival    = 2.0_wp * ddot(df_vars_%nQ,v_mn,1,int2(vw+1:vw+nQ),1)

                  ! 2-e exchange contribution - g(mw|vn)
                  vn      = df_pq_index(vdf,ndf) 

                  ival    = ival - ddot(df_vars_%nQ,v_mw,1,int2(vn+1:vn+nQ),1)

                  ! contract with integral/density matrix elements
                  val     = val + dval * ival

                end do ! end v loop

                ! **** w == v --> factor of 1

                ! 1-e density element
                den_ind = dens_%gemind(w,w)
                dval    = den1(den_ind)

                ! 2-e coulomb contribution 2 g(mn|ww)
                ww      = df_pq_index(wdf,wdf)

                ival    = ddot(df_vars_%nQ,v_mn,1,int2(ww+1:ww+nQ),1)

                ! 2-e exchange contribution - g(mw|wn)
                wn      = df_pq_index(wdf,ndf) 

                ival    = ival - 0.5_wp * ddot(df_vars_%nQ,v_mw,1,int2(wn+1:wn+nQ),1)

                ! contract with integral/density matrix elements
                val     = val + dval * ival

              end do ! end w loop

            end do ! end w_sym loop

            ! save Fock matrix element
            mn_fock          = ints_%gemind(m,n)
            fock_a_(mn_fock) = val

          end do ! end n loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

#ifdef OMP
!$omp parallel shared(n_class,nQ,ints_,dens_,m_class,m_sym,m,mdf,first_index_,last_index_,df_vars_,den1,int2) num_threads(nthread_use_)
!$omp do private(n,ndf,mn,val,v_mn,v_mw,w_sym,w,wdf,v,den_ind,dval,vdf,vw,ival,mw,wn,mn_fock)
#endif

            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              ! orbital index in df order
              ndf = df_vars_%class_to_df_map(n)

              ! mn-geminal index in df order
              mn  = df_pq_index(mdf,ndf) 

              call dcopy(df_vars_%nQ,int2(mn+1:mn+nQ),1,v_mn,1)

              ! initialize Fock matrix element
              val = 0.0_wp

              ! loop over irreps for v

              do w_sym = 1 , nirrep_

                ! loop over w indeces

                do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                  ! orbital index in df order
                  wdf = df_vars_%class_to_df_map(w)

                  mw      = df_pq_index(mdf,wdf)                  

                  call dcopy(df_vars_%nQ,int2(mw+1:mw+nQ),1,v_mw,1)

                  ! *** w > v --> factor of 2 

                  ! loop over v indeces

                  do v = first_index_(w_sym,2) , w -1 

                    ! 1-e density element
                    den_ind = dens_%gemind(v,w)
                    dval    = den1(den_ind)

                    ! orbital index in df order
                    vdf     = df_vars_%class_to_df_map(v)

                    ! 2-e coulomb contribution 2 g(mn|vw)
                    vw      = df_pq_index(wdf,vdf) 

                    ival    = 2.0_wp * ddot(df_vars_%nQ,v_mn,1,int2(vw+1:vw+nQ),1)

                    ! 2-e exchange contribution - g(mw|vn)
                    vn      = df_pq_index(vdf,ndf) 

                    ival    = ival - ddot(df_vars_%nQ,v_mw,1,int2(vn+1:vn+nQ),1)

                    ! contract with integral/density matrix elements
                    val     = val + dval * ival

                  end do ! end v loop

                  ! *** w == v --> factor of 1

                  ! 1-e density element
                  den_ind = dens_%gemind(w,w)
                  dval    = den1(den_ind)

                  ! 2-e coulomb contribution 2 g(mn|ww)
                  ww      = df_pq_index(wdf,wdf) 

                  ival    = ddot(df_vars_%nQ,v_mn,1,int2(ww+1:ww+nQ),1)

                  ! 2-e exchange contribution - g(mw|wn)
                  wn      = df_pq_index(wdf,ndf) 

                  ival    = ival - 0.5_wp * ddot(df_vars_%nQ,v_mw,1,int2(wn+1:wn+nQ),1)

                  ! contract with integral/density matrix elements
                  val     = val + dval * ival

                end do ! end w loop

              end do ! end w_sym loop

              ! save Fock matrix element
              mn_fock          = ints_%gemind(m,n)
              fock_a_(mn_fock) = val

            end do ! end n loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

          end do ! and n_class loop

        end do ! end m loop

      end do ! end m_sym loop

    end do ! end m_class loop

    return
  end subroutine compute_f_a_df

  subroutine compute_f_a(den1,int2)
    implicit none
    ! this subroutine computes the "active fock matrix" according to Eq. 12.5.13 in Helgaker on page 622
    ! F_a(m|n) = h(m|n) + SUM[v,w \in A] { g(mn|vw) - 0.5 g(mw|vn) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes.
    ! v and w belong to the active orbital class
    real(wp), intent(in) :: den1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,v,w,w_sym,mw_sym
    integer :: mn,vw,mw,vn,den_ind
    integer :: int_ind,sym_offset
    real(wp) :: val,ival,dval

    ! initialize
    fock_a_ = 0.0_wp

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          do n = first_index_(m_sym,m_class) , m

            ! mn-geminal index
            mn  = ints_%gemind(m,n)
 
            ! initialize Fock matrix element
            val = 0.0_wp
 
            ! loop over irreps for v

            do w_sym = 1 , nirrep_
  
              ! loop over w indeces

              do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                ! loop over v indeces

                do v = first_index_(w_sym,2) , last_index_(w_sym,2)

                  ! 1-e density element
                  den_ind    = dens_%gemind(v,w)
                  dval       = den1(den_ind)

                  ! 2-e coulomb contribution 2 g(mn|vw)
                  vw         = ints_%gemind(v,w)
                  int_ind    = pq_index(mn,vw)
                  ival       = int2(int_ind)

                  ! 2-e exchange contribution - g(mw|vn)
                  mw         = ints_%gemind(m,w)
                  vn         = ints_%gemind(v,n)
                  mw_sym     = group_mult_tab_(m_sym,w_sym)
                  sym_offset = ints_%offset(mw_sym)
                  int_ind    = sym_offset + pq_index(mw,vn)
                  ival       = ival - 0.5_wp * int2(int_ind)

                  ! contract with integral/density matrix elements
                  val        = val + dval * ival

                end do ! end v loop

              end do ! end w loop

            end do ! end w_sym loop

            ! save Fock matrix element
            fock_a_(mn) = val

          end do ! end n loop

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              ! mn-geminal index
              mn  = ints_%gemind(m,n)

              ! initialize Fock matrix element
              val = 0.0_wp

              ! loop over irreps for v

              do w_sym = 1 , nirrep_
  
                ! loop over w indeces

                do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                  ! loop over v indeces

                  do v = first_index_(w_sym,2) , last_index_(w_sym,2)

                    ! 1-e density element
                    den_ind    = dens_%gemind(v,w)
                    dval       = den1(den_ind)

                    ! 2-e coulomb contribution 2 g(mn|vw)
                    vw         = ints_%gemind(v,w)
                    int_ind    = pq_index(mn,vw)
                    ival       = int2(int_ind)

                    ! 2-e exchange contribution - g(mw|vn)
                    mw         = ints_%gemind(m,w)
                    vn         = ints_%gemind(v,n)
                    mw_sym     = group_mult_tab_(m_sym,w_sym)
                    sym_offset = ints_%offset(mw_sym)
                    int_ind    = sym_offset + pq_index(mw,vn)
                    ival       = ival - 0.5_wp * int2(int_ind)

                    ! contract with integral/density matrix elements
                    val        = val + dval * ival

                  end do ! end v loop

                end do ! end w loop
  
              end do ! end w_sym loop

              ! save Fock matrix element
              fock_a_(mn) = val

            end do ! end n loop

          end do ! and n_class loop

        end do ! end m loop

      end do ! end m_sym loop

    end do ! end m_class loop
   
    return
  end subroutine compute_f_a

  subroutine compute_f_i_df(int1,int2)
    ! this subroutine computes the "inactive fock matrix" according to Eq. 12.5.12 in Helgaker on page 622
    ! F_i(m|n) = h(m|n) + SUM[i \in D] { 2 g(mn|ii) - g(mi|in) }
    ! using 3-index 2-e integrals
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes. 
    implicit none
    real(wp), intent(in) :: int1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,i,i_sym,mdf,ndf,idf
    integer :: mn_fock,mn_int1
    integer(ip) :: mn,ii,mi,in,nQ
    real(wp) :: val,int_val,ddot
    real(wp) :: v_mn(df_vars_%nQ)

    nQ = int(df_vars_%nQ,kind=ip) 

    ! initialize
    fock_i_ = 0.0_wp

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

!#ifdef OMP
!!$omp parallel shared(nQ,first_index_,last_index_,fock_i_,int2,int1,ints_,df_vars_,m_class,m_sym) num_threads(nthread_use_)
!!$omp do private(m,mdf,n,mn_int1,val,ndf,mn,idf,ii,int_val,mi,in,mn_fock,n_class)
!#endif
        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! orbital index in df order
          mdf = df_vars_%class_to_df_map(m)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

#ifdef OMP
!$omp parallel shared(m,mdf,first_index_,last_index_,fock_i_,int2,int1,ints_,df_vars_,m_class,m_sym) num_threads(nthread_use_)
!$omp do private(n,v_mn,mn_int1,val,ndf,mn,idf,ii,int_val,mi,in,mn_fock,n_class)
#endif

          do n = first_index_(m_sym,m_class) , m

            ! mn-geminal index
            mn_int1 = ints_%gemind(m,n)

            ! 1-e contribution h(m|n)
            val     = int1(mn_int1)

            ! orbital index in df order
            ndf = df_vars_%class_to_df_map(n)

            ! mn-geminal index in df order
            mn      = df_pq_index(mdf,ndf) 

            call dcopy(df_vars_%nQ,int2(mn+1:mn+nQ),1,v_mn,1)

            ! loop over irreps for i

            do i_sym = 1 , nirrep_

              ! loop over i indeces

              do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                ! orbital index in df order
                idf        = df_vars_%class_to_df_map(i) 

                ! ii-geminal index in df order
                ii         = df_pq_index(idf,idf)

                ! 2-e coulomb contribution 2 g(mn|ii)

                int_val    = ddot(df_vars_%nQ,v_mn,1,int2(ii+1:ii+nQ),1)     

                val        = val + 2.0_wp * int_val

                ! geminal indeces in df order
                mi         = df_pq_index(mdf,idf) 
                in         = df_pq_index(idf,ndf)

                ! 2-e exchange contribution - g(mi|in)

                int_val    = ddot(df_vars_%nQ,int2(mi+1:mi+nQ),1,int2(in+1:in+nQ),1)

                val        = val - int_val

              end do ! end i loop

            end do ! end i_sym loop

            ! save Fock matrix element

            mn_fock          = ints_%gemind(m,n)
            fock_i_(mn_fock) = val

          end do ! end n loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

#ifdef OMP
!$omp parallel shared(n_class,m,mdf,first_index_,last_index_,fock_i_,int2,int1,ints_,df_vars_,m_class,m_sym) num_threads(nthread_use_)
!$omp do private(n,v_mn,mn_int1,val,ndf,mn,idf,ii,int_val,mi,in,mn_fock)
#endif
            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              mn_int1 = ints_%gemind(m,n)
              ! 1-e contribution h(m|n)
              val     = int1(mn_int1)

              ! orbital index in df order
              ndf     = df_vars_%class_to_df_map(n)

              ! mn geminal index in df order
              mn      = df_pq_index(mdf,ndf)

              call dcopy(df_vars_%nQ,int2(mn+1:mn+nQ),1,v_mn,1)

              ! loop over irreps for i

              do i_sym = 1 , nirrep_

                ! loop over i indeces

                do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                  ! orbital index in df order
                  idf        = df_vars_%class_to_df_map(i)

                  ! ii-geminal index in df order
                  ii         = df_pq_index(idf,idf) 

                  ! 2-e coulomb contribution 2 g(mn|ii)

                  int_val    = ddot(df_vars_%nQ,v_mn,1,int2(ii+1:ii+nQ),1)

                  val        = val + 2.0_wp * int_val

                  ! geminal indeces in df order
                  mi         = df_pq_index(mdf,idf) 
                  in         = df_pq_index(idf,ndf) 
 
                  ! 2-e exchange contribution - g(mi|in)
                  int_val    = ddot(df_vars_%nQ,int2(mi+1:mi+nQ),1,int2(in+1:in+nQ),1)

                  val        = val - int_val

                end do ! end i loop

              end do ! end i_sym loop

              ! save Fock matrix element
              mn_fock          = ints_%gemind(m,n)
              fock_i_(mn_fock) = val

            end do ! end n_loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

          end do ! end n_class loop

        end do ! end m_loop
!#ifdef OMP
!!$omp end do nowait
!!$omp end parallel
!#endif

      end do ! end m_sym loop

    end do ! end m_class loop

    return

  end subroutine compute_f_i_df

  subroutine compute_f_i(int1,int2)
    ! this subroutine computes the "inactive fock matrix" according to Eq. 12.5.12 in Helgaker on page 622
    ! F_i(m|n) = h(m|n) + SUM[i \in D] { 2 g(mn|ii) - g(mi|in) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes. 
    implicit none
    real(wp), intent(in) :: int1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,i,i_sym,mi_sym
    integer :: mn,ii,mi,in
    integer :: int_ind,sym_offset
    real(wp) :: val

    ! initialize
    fock_i_ = 0.0_wp

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          do n = first_index_(m_sym,m_class) , m

            ! mn-geminal index
            mn = ints_%gemind(m,n)
            ! 1-e contribution h(m|n)
            val = int1(mn)

            ! loop over irreps for i

            do i_sym = 1 , nirrep_

              ! loop over i indeces

              do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                ! 2-e coulomb contribution 2 g(mn|ii)
                ii         = ints_%gemind(i,i)
                int_ind    = pq_index(ii,mn)
                val        = val + 2.0_wp * int2(int_ind)

                ! 2-e exchange contribution - g(mi|in)
                mi         = ints_%gemind(m,i)
                in         = ints_%gemind(i,n)
                mi_sym     = group_mult_tab_(m_sym,i_sym)
                sym_offset = ints_%offset(mi_sym)
                int_ind    = sym_offset + pq_index(mi,in)
                val        = val - int2(int_ind) 

              end do ! end i loop

            end do ! end i_sym loop

            ! save Fock matrix element
            fock_i_(mn) = val

          end do ! end n loop

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              mn = ints_%gemind(m,n)
              ! 1-e contribution h(m|n)
              val = int1(mn) 

              ! loop over irreps for i

              do i_sym = 1 , nirrep_

                ! loop over i indeces

                do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                  ! 2-e coulomb contribution 2 g(mn|ii)
                  ii         = ints_%gemind(i,i)
                  int_ind    = pq_index(ii,mn)
                  val        = val + 2.0_wp * int2(int_ind)

                  ! 2-e exchange contribution - g(mi|in)
                  mi         = ints_%gemind(m,i)
                  in         = ints_%gemind(i,n)
                  mi_sym     = group_mult_tab_(m_sym,i_sym)
                  sym_offset = ints_%offset(mi_sym)
                  int_ind    = sym_offset + pq_index(mi,in)
                  val        = val - int2(int_ind)

                end do ! end i loop
 
              end do ! end i_sym loop

              ! save Fock matrix element
              fock_i_(mn) = val

            end do ! end n_loop

          end do ! end n_class loop

        end do ! end m_loop

      end do ! end m_sym loop
  
    end do ! end m_class loop

    return
  end subroutine compute_f_i

  integer function precondition_gradient(precond,gradient)
    implicit none
    ! simple subroutine to precondition the gradient according to
    ! g(i) = g(i)/precond(i) if (precond(i) >  0
    !      = g(i)            if (precond(i) <= 0
    real(wp) :: gradient(:)
    real(wp) :: precond(:)

    integer :: n_ij,ij_pair,nij_negative

    precondition_gradient = 1
  
    n_ij=size(gradient)
     
    if (n_ij /= size(precond) ) return

    nij_negative = 0

    do ij_pair = 1 , n_ij

      if ( precond(ij_pair) > 0.0_wp ) then
        gradient(ij_pair) = gradient(ij_pair) / precond(ij_pair)
      else
        gradient(ij_pair) = 0.0_wp
        nij_negative = nij_negative + 1
      endif    

    end do

    precondition_gradient = 0

    if ( ( nij_negative /= 0 ) .and. ( log_print_ == 1 ) ) then
      write(fid_,'(a,1x,i9,1x,a)')'!warning!',nij_negative,'diagonal Hessian elements were negative'
    end if

    return

  end function precondition_gradient
 
  subroutine allocate_temporary_fock_matrices()
    implicit none
    ! the total number of nonzero LT elements in the active/inactive Fock matrix 
    ! is equal to the number of geminals in the totally symmetric irrep
    allocate(fock_i_(ints_%ngempi(1)))
    allocate(fock_a_(ints_%ngempi(1)))
    ! the row index here will be between 1-nact, while the column index is arbitrary
    allocate(q_(nact_tot_,nmo_tot_))
    allocate(z_(nact_tot_,nmo_tot_))
    ! generalized Fock matrix
    !
    !        | F_dd F_da F_dv |   | F_dd F_da F_dv |
    !F_gen = | F_ad F_aa F_av | = | F_ad F_aa F_av |
    !        | F_vd F_va F_vv |   |  0    0    0   |
    !
    ! thus here we allocate more than needed, but that is not too much memory
    !allocate(fock_gen_(nmo_tot_,nmo_tot_))
    return
  end subroutine allocate_temporary_fock_matrices

  subroutine deallocate_temporary_fock_matrices()
    implicit none
    if (allocated(fock_i_))   deallocate(fock_i_)
    if (allocated(fock_a_))   deallocate(fock_a_)
    if (allocated(q_))        deallocate(q_)
    if (allocated(z_))        deallocate(z_)
    !if (allocated(fock_gen_)) deallocate(fock_gen_)
    return
  end subroutine deallocate_temporary_fock_matrices

end module focas_gradient
