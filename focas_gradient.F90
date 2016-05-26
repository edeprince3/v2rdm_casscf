!!
 !@BEGIN LICENSE
 !
 ! v2RDM-CASSCF, a plugin to:
 !
 ! PSI4: an ab initio quantum chemistry software package
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !@END LICENSE
 !
 !!

module focas_gradient
  use focas_data
 
  implicit none

  contains

  subroutine orbital_gradient(int1,int2,den1,den2)
    implicit none
    real(wp), intent(in) :: int1(:),int2(:),den1(:),den2(:)
    real(wp) :: t0,t1
    type(fock_info) :: fock
    integer :: i
    real(wp), allocatable :: tq(:,:)
   

    ! calculate inactive Fock matrix
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_f_i(int1,int2)
    else
      call compute_f_i_df_coulomb(int1,int2)
      call compute_f_i_df_exchange(int1,int2)
    endif
    call transpose_matrix(fock_i_)


    ! calculate active Fock matrix
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_f_a(den1,int2)
    else
      call compute_f_a_df_coulomb(den1,int2)
      call compute_f_a_df_exchange(den1,int2)
    endif
    call transpose_matrix(fock_a_)

    ! calculate auxiliary q matrix
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_q(den2,int2)
    else
      call compute_q_df(den2,int2)
    endif

!    if ( log_print_ == 1) then
!
!      write(fid_,*)'*** INACTIVE ***'
!      call print_f_matrix(fock_i_)
!
!      write(fid_,*)'*** ACTIVE ***'
!      call print_f_matrix(fock_a_)
!
!      write(fid_,*)'*** AUXILIARY Q ***'
!      call print_q_matrix(q_)
!
!    end if

    ! calculate auxiliary z matrix
    call compute_z(den1)

    ! compute gradient
    call compute_orbital_gradient()

    ! determine value,type, and orbital indices for largest gradient element
    call check_max_gradient()

    ! print the gradient 
    !if ( log_print_ == 1 ) call print_vector(orbital_gradient_)

    return
  end subroutine orbital_gradient

  subroutine check_q(q1,q2)
   real(wp), intent(in) :: q1(:,:),q2(:,:)
   integer :: p_sym,p,u_class,u

   do p_sym =1,nirrep_
     do p=1,nactpi_(p_sym)
       do u_class=1,3
         do u=first_index_(p_Sym,u_class),last_index_(p_Sym,u_class)
           write(*,'(i1,5x,2(i3,1x),10x,2(f10.6,1x),5x,es12.4)')u_class,p,u,q1(p,u),q2(p,u),q1(p,u)-q2(p,u)
         end do
       end do
     end do
   end do

  end subroutine 

  subroutine print_vector(vector)

    implicit none

    real(wp), intent(in) :: vector(:)

    integer :: a_sym,t_sym,i,a,t,u,grad_ind,newline
    character(1) :: a_typ,i_typ,t_typ,u_typ

    newline = 0

    ! ******************************
    ! doubly occupied - active pairs
    ! ******************************

    i_typ='d'
    t_typ='a'

    do t_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

      do i = first_index_(t_sym,1) , last_index_(t_sym,1)

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

          grad_ind = grad_ind + 1

          write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')t_typ,'-',i_typ,&
               & ' (',t,',',i,')',vector(grad_ind)

          newline = newline + 1

          if ( mod(newline,4) == 0 ) write(fid_,*)

        end do

      end do

    end do
 
    ! ********************************
    ! doubly occupied - external pairs
    ! ********************************

    i_typ='d'
    a_typ='e'
    
    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')a_typ,'-',i_typ,&
               & ' (',a,',',i,')',vector(grad_ind)

          newline = newline + 1

          if ( mod(newline,4) == 0 ) write(fid_,*)

        end do

      end do

    end do

    ! *********************
    ! active - active pairs
    ! *********************

    t_typ='a'
    u_typ='a'

    if ( include_aa_rot_ == 1 ) then

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1

            write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')t_typ,'-',u_typ,&
                 & ' (',t,',',u,')',vector(grad_ind)

            newline = newline + 1

            if ( mod(newline,4) == 0 ) write(fid_,*)

          end do

        end do

      end do

    end if

    ! ***********************
    ! active - external pairs
    ! ***********************

    t_typ='a'
    a_typ='e'

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_act_type)

      do t = first_index_(a_sym,2) , last_index_(a_sym,2)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')a_typ,'-',t_typ,&
               & ' (',a,',',t,')',vector(grad_ind)

          newline = newline + 1

          if ( mod(newline,4) == 0 ) write(fid_,*)

        end do

      end do

    end do

    if ( mod(newline,4) /= 0 ) write(fid_,*)

    return

  end subroutine print_vector

  subroutine transpose_matrix(fock)

    implicit none
    
    type(fock_info) :: fock

    integer :: num_p,p,q,p_sym

    ! subroutine to symmetrize the inactive-inactive and active-active 
    ! blocks of the Fock matrix
        
    do p_sym = 1 , nirrep_

      num_p = ndocpi_(p_sym) + nactpi_(p_sym)

      if ( num_p == 0 ) cycle

      do p = 1 , num_p

        do q = p+1 , num_p

          fock%occ(p_sym)%val(p,q) = fock%occ(p_sym)%val(q,p)

        end do

      end do

    end do

  end subroutine transpose_matrix

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

  subroutine check_max_gradient()

    implicit none

    integer  :: grad_ind,a_sym,t_sym,i,a,t,u
    real(wp) :: abs_grad_val,grad_tol

    ! ********************************
    ! external - doubly-occupied pairs
    ! ********************************

    max_grad_val_ = 0.0_wp
    max_grad_ind_ = 0
    max_grad_typ_ = 0
    max_grad_sym_ = 0

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < max_grad_val_ ) cycle

          max_grad_val_    = abs_grad_val 

          max_grad_ind_(1) = a
          max_grad_ind_(2) = i

          max_grad_typ_    = 2

          max_grad_sym_    = a_sym

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

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < max_grad_val_ ) cycle

          max_grad_val_    = abs_grad_val

          max_grad_ind_(1) = t
          max_grad_ind_(2) = i

          max_grad_typ_    = 1

          max_grad_sym_    = t_sym

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

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < max_grad_val_ ) cycle

          max_grad_val_    = abs_grad_val

          max_grad_ind_(1) = a
          max_grad_ind_(2) = t

          max_grad_typ_    = 4

          max_grad_sym_    = a_sym

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

            abs_grad_val = abs(orbital_gradient_(grad_ind))

            if ( abs_grad_val < max_grad_val_ ) cycle

            max_grad_val_    = abs_grad_val

            max_grad_ind_(1) = t
            max_grad_ind_(2) = u

            max_grad_typ_    = 3

            max_grad_sym_    = t_sym

          end do

        end do

      end do

    end if

    ! *************************************************************
    ! check how many elements are within 75% of the largest element
    ! *************************************************************

    grad_tol         = 0.75_wp * max_grad_val_

    norm_grad_large_ = 0.0_wp

    n_grad_large_    = 0

    ! ********************************
    ! external - doubly-occupied pairs
    ! ********************************

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < grad_tol ) cycle

          norm_grad_large_ = norm_grad_large_ + abs_grad_val * abs_grad_val

          n_grad_large_ = n_grad_large_ + 1

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

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < grad_tol ) cycle

          norm_grad_large_ = norm_grad_large_ + abs_grad_val * abs_grad_val

          n_grad_large_ = n_grad_large_ + 1

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

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < grad_tol ) cycle

          norm_grad_large_ = norm_grad_large_ + abs_grad_val * abs_grad_val

          n_grad_large_ = n_grad_large_ + 1

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

            abs_grad_val = abs(orbital_gradient_(grad_ind))

            if ( abs_grad_val < grad_tol ) cycle

            norm_grad_large_ = norm_grad_large_ + abs_grad_val * abs_grad_val

            n_grad_large_ = n_grad_large_ + 1

          end do

        end do

      end do

    end if

    norm_grad_large_ = sqrt(norm_grad_large_)

    return

  end subroutine check_max_gradient

  subroutine compute_orbital_gradient()

    ! subroutine to compute the orbital gradient without explicit storage
    ! of the generalized Fock matrix

    implicit none

    integer :: grad_ind
    integer :: a_sym,t_sym
    integer :: i,t,u,a
    integer :: ia,it,i_i,a_i,t_i
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

          i_i      = trans_%class_to_irrep_map(i)
          a_i      = trans_%class_to_irrep_map(a)

          ia       = ints_%gemind(i,a)

          orbital_gradient_(grad_ind) = 4.0_wp*( fock_i_%occ(a_sym)%val(a_i,i_i) &
                                       &       + fock_a_%occ(a_sym)%val(a_i,i_i) )

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

          i_i      = trans_%class_to_irrep_map(i)
          t_i      = trans_%class_to_irrep_map(t)

          it       = ints_%gemind(i,t)

          orbital_gradient_(grad_ind) = 4.0_wp*( fock_i_%occ(t_sym)%val(t_i,i_i)    & 
                                      &  +       fock_a_%occ(t_sym)%val(t_i,i_i) )  &
                                      & - 2.0_wp * ( q_(t - ndoc_tot_,i) + z_(t - ndoc_tot_,i))

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

            orbital_gradient_(grad_ind) = 2.0_wp * ( q_(u - ndoc_tot_,t) +         &
                 & z_(u - ndoc_tot_,t) - q_(t - ndoc_tot_,u) - z_(t - ndoc_tot_,u))

          end do
 
        end do

      end do

    end if

    grad_norm_ = sqrt(ddot(rot_pair_%n_tot,orbital_gradient_,1,orbital_gradient_,1))

    return 

  end subroutine compute_orbital_gradient

  subroutine compute_z(den1)
    implicit none
    ! function to compute contraction of the density with the inactive Fock matrix according to
    ! z(m,t) = sum_u { den1(tu) * f_i(mu) } t,u \in A m \in D,A,E
    real(wp), intent(in) :: den1(:)
    integer :: t,u,tu_den,mu_int,t_sym,m,m_class,m_i,u_i
    real(wp) :: val

    ! initialize output

    z_ = 0.0_wp

    ! loop over symmetries for t

    do t_sym = 1 , nirrep_

      ! loop over u indeces

      do t = first_index_(t_sym,2) , last_index_(t_sym,2)
 
        ! m_class < u_class ==> m < u

        ! loop over m orbital classes

        do m_class = 1 , 3
 
          do m = first_index_(t_sym,m_class) , last_index_(t_sym,m_class)

            m_i = trans_%class_to_irrep_map(m)

            ! initialize contraction value

            val = 0.0_wp

            do u = first_index_(t_sym,2) , last_index_(t_sym,2) 

              ! integral/density addressing

              tu_den = dens_%gemind(t,u)
              u_i    = trans_%class_to_irrep_map(u)

              val = val + den1(tu_den) * fock_i_%occ(t_sym)%val(m_i,u_i)

            end do ! end u loop

            z_( t - ndoc_tot_ , m ) = val

          end do ! end m loop

        end do ! end m_class loop 

      end do ! end t loop

    end do ! end t_sym loop

!
!    ! loop over symmetries for t
!
!    do t_sym = 1 , nirrep_
!
!      ! loop over u indeces
!
!      do t = first_index_(t_sym,2) , last_index_(t_sym,2)
!
!        ! m_class < u_class ==> m < u
!
!        m_class = 1 
!
!        ! loop over m indeces
!
!        do m = first_index_(t_sym,m_class) , last_index_(t_sym,m_class)
!
!          m_i = trans_%class_to_irrep_map(m)
!
!          ! initialize contraction value
!
!          val = 0.0_wp
!
!          do u = first_index_(t_sym,2) , last_index_(t_sym,2)
!
!            ! integral/density addressing
!
!            tu_den = dens_%gemind(t,u)
!            u_i    = trans_%class_to_irrep_map(u)
!
!            val = val + den1(tu_den) * fock_i_%occ(t_sym)%val(u_i,m_i)
!
!          end do ! end u loop
!
!          z_( t - ndoc_tot_ , m ) = val
!
!        end do ! end m loop
!
!        ! m_class = u_class ==> u <= m
!
!        m_class = 2
!
!        ! loop over m indeces
!
!        do m = first_index_(t_sym,m_class) , last_index_(t_sym,m_class)
!
!          m_i = trans_%class_to_irrep_map(m)
!
!          ! initialize contraction value
!
!          val = 0.0_wp
!
!          do u = first_index_(t_sym,2) , m - 1
!
!            ! integral/density addressing
!
!            tu_den = dens_%gemind(t,u)
!            u_i    = trans_%class_to_irrep_map(u)
!
!            val = val + den1(tu_den) * fock_i_%occ(t_sym)%val(m_i,u_i)
!
!          end do ! end u loop
!          
!          tu_den = dens_%gemind(t,t)
!          
!          val    = val + den1(tu_den) * fock_i_%occ(t_sym)%val(m_i,m_i)
!
!          do u = m + 1 , last_index_(t_sym,2)
!
!            ! integral/density addressing
!
!            tu_den = dens_%gemind(t,u)
!            u_i    = trans_%class_to_irrep_map(u)
!
!            val = val + den1(tu_den) * fock_i_%occ(t_sym)%val(u_i,m_i)
!
!          end do ! end u loop
!
!          z_( t - ndoc_tot_ , m ) = val
!
!        end do ! end m loop
!
!        ! m_class > u_class ==> u < m
!
!        m_class = 3
!
!        ! loop over m indeces
!
!        do m = first_index_(t_sym,m_class) , last_index_(t_sym,m_class)
!
!          m_i = trans_%class_to_irrep_map(m)
!
!          ! initialize contraction value
!
!          val = 0.0_wp
!
!          do u = first_index_(t_sym,2) , last_index_(t_sym,2)
!
!            ! integral/density addressing
!
!            tu_den = dens_%gemind(t,u)
!            u_i    = trans_%class_to_irrep_map(u)
!
!            val = val + den1(tu_den) * fock_i_%occ(t_sym)%val(m_i,u_i)
!
!          end do ! end u loop
!
!          z_( t - ndoc_tot_ , m ) = val
!
!        end do ! end m loop
!
!      end do ! end t loop
!
!    end do ! end t_sym loop

    return

  end subroutine compute_z

  subroutine compute_q_df(den2,int2)

    implicit none

    real(wp), intent(in) :: den2(:),int2(:)

    integer :: tu_sym,t_sym,u_sym,v_sym,w_sym,p_sym
    integer :: p_class
    integer :: t,u,v,w,p,tu_den,vw_den,tuvw
    integer :: vdf,wdf,pdf,udf,den_off,den_ind,tu_int
    integer(ip) :: vw_df,pu_df,nQ
  
    real(wp) :: val,ddot
 
    nQ = df_vars_%nQ

    ! initialize
 
    q_ = 0.0_wp
 
    do tu_sym = 1 , nirrep_

      if ( dens_%ngempi(tu_sym) == 0 ) cycle

      qint_%tuQ(tu_sym)%val = 0.0_wp

    end do

    ! assemble intermediates

    tu_sym = 1

    do t_sym = 1 , nirrep_

      do t = first_index_(t_sym,2) , last_index_(t_sym,2)

#ifdef OMP
!$omp parallel shared(t_sym,t,tu_sym,first_index_,last_index_, &
!$omp int2,den2,df_vars_,dens_,qint_) num_threads(nthread_use_)
!$omp do private(v,vdf,w,wdf,vw_df,vw_den,den_ind,u,tu_den,v_sym)
#endif

        do u = first_index_(t_sym,2) , t

          tu_den = dens_%gemind(t,u)

          do v_sym = 1 , nirrep_

            do v = first_index_(v_sym,2) , last_index_(v_sym,2)

              vdf = df_vars_%class_to_df_map(v)

              do w = first_index_(v_sym,2) , v - 1

                wdf = df_vars_%class_to_df_map(w)

                vw_df =  df_pq_index(vdf,wdf)

                vw_den = dens_%gemind(v,w)

                den_ind = pq_index(tu_den,vw_den) 
                ! update array

                call daxpy(nQ,2.0_wp*den2(den_ind),int2(vw_df+1:vw_df+nQ),1,qint_%tuQ(tu_sym)%val(:,tu_den),1)

              end do ! end w loop

              vw_df =  df_pq_index(vdf,vdf)

              vw_den = dens_%gemind(v,v)

              den_ind = pq_index(tu_den,vw_den)

              ! update array

              call daxpy(nQ,den2(den_ind),int2(vw_df+1:vw_df+nQ),1,qint_%tuQ(tu_sym)%val(:,tu_den),1) 

            end do ! end v loop

          end do ! end v_sym loop

        end do ! end u loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

      end do ! end t loop

    end do ! end t_sym loop

    do tu_sym = 2 , nirrep_
 
      den_off = dens_%offset(tu_sym)
 
      do t_sym = 1 , nirrep_

        u_sym=group_mult_tab_(tu_sym,t_sym)
         
        if ( u_sym > t_sym ) cycle

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

#ifdef OMP
!$omp parallel shared(t_sym,u_sym,tu_sym,t,den_off,      &
!$omp first_index_,last_index_,int2,den2,df_vars_,dens_, &
!$omp qint_) num_threads(nthread_use_)
!$omp do private(v,vdf,w,wdf,vw_df,vw_den,den_ind,tu_den,u,v_sym,w_sym)
#endif

          do u = first_index_(u_sym,2) , last_index_(u_sym,2)

            tu_den  = dens_%gemind(t,u)

            do v_sym = 1 , nirrep_

              w_sym = group_mult_tab_(tu_sym,v_sym)

              if ( w_sym > v_sym ) cycle

              do v = first_index_(v_sym,2) , last_index_(v_sym,2)

                vdf = df_vars_%class_to_df_map(v)

                do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                  wdf = df_vars_%class_to_df_map(w)

                  vw_df =  df_pq_index(vdf,wdf)

                  vw_den = dens_%gemind(v,w)

                  den_ind = pq_index(tu_den,vw_den) + den_off

                  ! update array

                  call daxpy(nQ,2.0_wp*den2(den_ind),int2(vw_df+1:vw_df+nQ),1,qint_%tuQ(tu_sym)%val(:,tu_den),1)

                end do ! end w loop

              end do ! end v loop

            end do ! end v_sym loop

          end do ! end u_loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

        end do ! end t_loop

      end do ! end t_sym loop
  
    end do ! end tu_sym loop

    ! compute Q-elements

    do p_sym = 1 , nirrep_

      do p_class = 1 , 3

        do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

          pdf = df_vars_%class_to_df_map(p)

          do t = first_index_(p_sym,2) , last_index_(p_sym,2)

            val = 0.0_wp

            do tu_sym = 1 , nirrep_

              u_sym = group_mult_tab_(tu_sym,p_sym)

              do u = first_index_(u_sym,2) , last_index_(u_sym,2)

                tu_den = dens_%gemind(t,u)

                udf = df_vars_%class_to_df_map(u)

                pu_df = df_pq_index(pdf,udf)

                val   = val + ddot(nQ,int2(pu_df+1:pu_df+nQ),1,qint_%tuQ(tu_sym)%val(:,tu_den),1)

              end do ! end u loop

            end do ! end tu_sym loop

            q_(t - ndoc_tot_ , p ) = val

          end do ! end t loop

        end do ! end p loop

      end do ! and p_class loop

    end do ! end p_sym_loop

    return
  
  end subroutine compute_q_df

  subroutine compute_q_df_old(den2,int2)
    implicit none
    ! this subroutine computes the "auxiliary Q" matrix according to Eq. 12.5.14 in Helgaker on page 622
    ! Q(v,m) = \SUM[w,x,y \in A] { d2(vw|xy) * g(mw|xy) }
    ! Because only those D2 elements with all four indeces \in A are nonzero, we have v \in A
    ! Also, m is a general orbital index && sym_m == sym_v
    real(wp), intent(in) :: den2(:),int2(:)
    integer :: mw_sym,x_sym,y_sym,w_sym,m_sym,m_class,m,v,w,x,y,mdf,vdf,wdf,xdf,ydf
    integer :: xy_den,vw_den
    integer(ip) :: mw,xy,nQ
    integer :: den_ind,den_sym_offset
    real(wp) :: val,int_val,ddot
    real(wp) :: v_mw(df_vars_%nQ)

    nQ = int(df_vars_%nQ,kind=ip)

    ! initialize
    q_ = 0.0_wp

    ! loop over classes for m

    do m_class = 1 , 3

      ! loop over symmetries for m

      do m_sym = 1 , nirrep_

        ! loop over orbital index m

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! m index in df order

          mdf = df_vars_%class_to_df_map(m)

          ! loop over orbital index v

#ifdef OMP
!$omp parallel shared(m_class,nQ,m_sym,m,mdf,first_index_,last_index_, &
!$omp df_vars_,den2,int2,dens_,q_) num_threads(nthread_use_)
!$omp do private(val,v,v_mw,vdf,w_sym,mw_sym,den_Sym_offset,w,wdf,mw,  &
!$omp vw_den,x_sym,y_sym,x,y,xdf,ydf,xy,xy_den,den_ind,int_val)
#endif

          do v = first_index_(m_sym,2) , last_index_(m_sym,2)

            ! v index in df order

            vdf = df_vars_%class_to_df_map(v)

            ! initialize matrix element
  
            val = 0.0_wp

            ! loop over symmetries for w

            do w_sym = 1 , nirrep_

              ! symmetry of mw geminal

              mw_sym = group_mult_tab_(m_sym,w_sym)

              ! vw density offset 
             
              den_sym_offset = dens_%offset(mw_sym)

              ! loop over w orbital index
 
              do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                ! w indef in df order
           
                wdf    = df_vars_%class_to_df_map(w)

                ! mw geminal index in df order

                mw     = df_pq_index(mdf,wdf) 

                call dcopy(df_vars_%nQ,int2(mw+1:mw+nQ),1,v_mw,1)
              
                ! vw geminal index for density

                vw_den = dens_%gemind(v,w)

                ! loop over x_sym

                do x_sym = 1 , nirrep_

                  ! figure out symmetry of y

                  y_sym = group_mult_tab_(mw_sym,x_sym)

                  if ( y_sym == x_sym ) then
                 
                    ! loop over x orbital index
 
                    do x = first_index_(x_sym,2) , last_index_(x_sym,2)

                      ! x index in df order
     
                      xdf = df_vars_%class_to_df_map(x)

                      ! ***********************************************
                      ! loop over y orbital index (x>y) --> factor of 2
                      ! ***********************************************

                      do y = first_index_(y_sym,2) , x - 1

                        ! y index in df order
                        
                        ydf     = df_vars_%class_to_df_map(y)

                        ! xy geminal index and density index

                        xy_den  = dens_%gemind(x,y)
                        den_ind = pq_index(vw_den,xy_den) + den_sym_offset
                    
                        ! xy geminal index in df order
                        xy  = df_pq_index(xdf,ydf)   

                        int_val = ddot(df_vars_%nQ,v_mw,1,int2(xy+1:xy+nQ),1)

                        ! update matrix element

                        val = val + 2.0_wp * int_val * den2(den_ind)

                      end do ! end y loop

                      ! **********************
                      ! x == y --> factor of 1
                      ! **********************

                      xy_den  = dens_%gemind(x,x)
                      den_ind = pq_index(vw_den,xy_den) + den_sym_offset

                      ! xy geminal index in df order
                      xy  = df_pq_index(xdf,xdf) 

                      int_val = ddot(df_vars_%nQ,v_mw,1,int2(xy+1:xy+nQ),1)

                      ! update matrix element

                      val = val + int_val * den2(den_ind)

                    end do ! end x loop

                  elseif ( y_sym < x_sym) then

                    ! loop over x orbital index

                    do x = first_index_(x_sym,2) , last_index_(x_sym,2)

                      ! x index in df order

                      xdf = df_vars_%class_to_df_map(x)

                      ! ***********************************************
                      ! loop over y orbital index (x>y) --> factor of 2
                      ! ***********************************************

                      do y = first_index_(y_sym,2) , last_index_(y_sym,2)

                        ! y index in df order

                        ydf     = df_vars_%class_to_df_map(y)

                        ! xy geminal index and density index

                        xy_den  = dens_%gemind(x,y)
                        den_ind = pq_index(vw_den,xy_den) + den_sym_offset

                        ! xy geminal index in df order
                        xy  = df_pq_index(xdf,ydf) 

                        int_val = ddot(df_vars_%nQ,v_mw,1,int2(xy+1:xy+nQ),1)

                        ! update matrix element

                        val = val + 2.0_wp * int_val * den2(den_ind)

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

      end do ! end m_sym loop 

    end do ! end m_class loop

    return

  end subroutine compute_q_df_old

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

  subroutine compute_f_a_df_coulomb(den1,int2)

    implicit none

    real(wp), intent(in) :: den1(:), int2(:)

    integer :: p_class,q_class,q_class_max
    integer :: tu_den,qt_int
    integer :: t_sym,u_sym,p_sym,tu_sym,q_sym,qt_sym
    integer :: t,u,p,q,q_max,q_min
    integer :: p_i,q_i
    integer(ip) :: tu_df,pq_df,qu_df,pt_df,nQ
    integer :: tdf,udf,pdf,qdf

    real(wp) :: val,ddot

    nQ = df_vars_%nQ

    ! initialize
    do p_sym = 1 , nirrep_

      if ( allocated( fock_a_%occ(p_sym)%val ) ) fock_a_%occ(p_sym)%val = 0.0_wp
      if ( allocated( fock_a_%ext(p_sym)%val ) ) fock_a_%ext(p_sym)%val = 0.0_wp

    end do

    ! *** Coulomb terms ***
    
    ! initialize
    qint_%tuQ(1)%val(:,1) = 0.0_wp 

    ! gather intermediate
 
    do t_sym = 1 , nirrep_

      do t = first_index_(t_sym,2) , last_index_(t_sym,2)

        tdf = df_vars_%class_to_df_map(t)

        do u = first_index_(t_sym,2) , t - 1

          udf = df_vars_%class_to_df_map(u)

          tu_den = dens_%gemind(t,u)

          tu_df  = df_pq_index(tdf,udf) 

          call daxpy(nQ,2.0_wp*den1(tu_den),int2(tu_df+1:tu_df+nQ),1,qint_%tuQ(1)%val(:,1),1)        

        end do ! end u loop

        tu_den = dens_%gemind(t,t)

        tu_df = df_pq_index(tdf,tdf)

        call daxpy(nQ,den1(tu_den),int2(tu_df+1:tu_df+nQ),1,qint_%tuQ(1)%val(:,1),1)

      end do ! end t loop

    end do ! end t_sym loop

    ! calculate fock matrix element

    do p_class = 1 , 3

      do q_class = 1 , 3

        if ( q_class > p_class ) cycle

        do p_sym = 1 , nirrep_

          do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)
 
            pdf = df_vars_%class_to_df_map(p)

            q_max = last_index_(p_sym,q_class)

            q_min = first_index_(p_sym,q_class)

            if ( p_class == q_class ) q_max = p
   
            if ( q_class == 3 ) q_min = p

            do q = q_min , q_max

              qdf = df_vars_%class_to_df_map(q)

              pq_df = df_pq_index(pdf,qdf)

              val = ddot(nQ,int2(pq_df+1:pq_df+nQ),1,qint_%tuQ(1)%val(:,1),1)

              if ( ( p_class == 3 ) .and. ( q_class == 3 ) ) then

                p_i = trans_%class_to_irrep_map(p)-ndocpi_(p_sym)-nactpi_(p_sym)
                fock_a_%ext(p_sym)%val(p_i) = val

              else

                p_i = trans_%class_to_irrep_map(p)
                q_i = trans_%class_to_irrep_map(q)
                fock_a_%occ(p_sym)%val(p_i,q_i) = val

              endif
              
            end do ! end q loop

          end do ! end p loop
          
        end do ! end p_sym loop     

      end do ! end q_class loop

    end do ! end p_class loop

  end subroutine compute_f_a_df_coulomb

  subroutine compute_f_a_df_exchange(den1,int2)
    implicit none
    ! this subroutine computes the "active fock matrix" according to Eq. 12.5.13 in Helgaker on page 622
    ! F_a(m|n) = h(m|n) + SUM[v,w \in A] { g(mn|vw) - 0.5 g(mw|vn) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes.
    ! v and w belong to the active orbital class
    real(wp), intent(in) :: den1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,v,w,w_sym,mdf,ndf,wdf,vdf
    integer :: den_ind,n_first,m_i,n_i
    integer(ip) :: mw,vn,wn,ww,nQ 
    real(wp) :: val,dval,ddot
    real(wp) :: v_mw(df_vars_%nQ)

    nQ = int(df_vars_%nQ,kind=ip)

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! orbital index in df order
          mdf = df_vars_%class_to_df_map(m)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          n_first = first_index_(m_sym,m_class)

          ! only compute diagonal elements for external-external block

          if ( m_class == 3 ) n_first = m

#ifdef OMP
!$omp parallel shared(ints_,nQ,dens_,m_class,m_sym,m,mdf,first_index_,  &
!$omp last_index_,df_vars_,den1,int2,n_first,fock_a_,trans_) num_threads(nthread_use_)
!$omp do private(n,ndf,val,v_mw,w_sym,w,wdf,v,den_ind,dval,vdf, &
!$omp ww,mw,wn,vn,n_i,m_i)
#endif

          do n = n_first , m

            ! orbital index in df order
            ndf = df_vars_%class_to_df_map(n)

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

                ! loop over v indeces

                do v = first_index_(w_sym,2) , last_index_(w_sym,2) !w - 1 

                  ! 1-e density element
                  den_ind = dens_%gemind(v,w)
                  dval    = den1(den_ind)

                  ! orbital index in df order
                  vdf     = df_vars_%class_to_df_map(v)

                  ! 2-e exchange contribution - g(mw|vn)
                  vn      = df_pq_index(vdf,ndf) 

                  ! contract with integral/density matrix elements
                  val     = val - dval * ddot(df_vars_%nQ,v_mw,1,int2(vn+1:vn+nQ),1)

                end do ! end v loop

              end do ! end w loop

            end do ! end w_sym loop

            ! save Fock matrix element

            if ( m_class == 3 ) then

              m_i = trans_%class_to_irrep_map(m)-ndocpi_(m_sym)-nactpi_(m_sym)

              fock_a_%ext(m_sym)%val(m_i) = fock_a_%ext(m_sym)%val(m_i) + 0.5_wp * val

            else

              m_i = trans_%class_to_irrep_map(m)
              n_i = trans_%class_to_irrep_map(n)

              fock_a_%occ(m_sym)%val(m_i,n_i) = fock_a_%occ(m_sym)%val(m_i,n_i) + 0.5_wp * val

            endif

          end do ! end n loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

#ifdef OMP
!$omp parallel shared(n_class,nQ,ints_,dens_,m_class,m_sym,m,mdf,first_index_, &
!$omp last_index_,df_vars_,den1,int2,fock_a_,trans_) num_threads(nthread_use_)
!$omp do private(n,ndf,val,v_mw,w_sym,w,wdf,v,den_ind,dval,vdf,&
!$omp mw,wn,n_i,m_i)
#endif

            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              ! orbital index in df order
              ndf = df_vars_%class_to_df_map(n)

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

                  ! loop over v indeces

                  do v = first_index_(w_sym,2) , last_index_(w_sym,2) !w -1 

                    ! 1-e density element
                    den_ind = dens_%gemind(v,w)
                    dval    = den1(den_ind)

                    ! orbital index in df order
                    vdf     = df_vars_%class_to_df_map(v)

                    ! 2-e exchange contribution - g(mw|vn)
                    vn      = df_pq_index(vdf,ndf) 

                    ! contract with integral/density matrix elements
                    val     = val - dval * ddot(df_vars_%nQ,v_mw,1,int2(vn+1:vn+nQ),1)

                  end do ! end v loop

                end do ! end w loop

              end do ! end w_sym loop

              ! save Fock matrix element

              m_i = trans_%class_to_irrep_map(m)
              n_i = trans_%class_to_irrep_map(n)

              fock_a_%occ(m_sym)%val(m_i,n_i) = fock_a_%occ(m_sym)%val(m_i,n_i) + 0.5_wp * val

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
  end subroutine compute_f_a_df_exchange

  subroutine compute_f_a(den1,int2)
    implicit none
    ! this subroutine computes the "active fock matrix" according to Eq. 12.5.13 in Helgaker on page 622
    ! F_a(m|n) = h(m|n) + SUM[v,w \in A] { g(mn|vw) - 0.5 g(mw|vn) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes.
    ! v and w belong to the active orbital class
    real(wp), intent(in) :: den1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,v,w,w_sym,mw_sym
    integer :: mn,vw,mw,vn,den_ind,n_first,m_i,n_i
    integer :: int_ind,sym_offset
    real(wp) :: val,ival,dval

    ! initialize
    do m_sym = 1 , nirrep_

      if ( allocated( fock_a_%occ(m_sym)%val ) ) fock_a_%occ(m_sym)%val = 0.0_wp
      if ( allocated( fock_a_%ext(m_sym)%val ) ) fock_a_%ext(m_sym)%val = 0.0_wp

    end do

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          n_first = first_index_(m_sym,m_class)
        
          if ( m_class == 3 ) n_first = m

          do n = n_first , m

            ! mn-geminal index
            mn  = ints_%gemind(m,n)
 
            ! initialize Fock matrix element
            val = 0.0_wp
 
            ! loop over irreps for v

            do w_sym = 1 , nirrep_
  
              mw_sym     = group_mult_tab_(m_sym,w_sym)
              sym_offset = ints_%offset(mw_sym)


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
                  int_ind    = sym_offset + pq_index(mw,vn)
                  ival       = ival - 0.5_wp * int2(int_ind)

                  ! contract with integral/density matrix elements
                  val        = val + dval * ival

                end do ! end v loop

              end do ! end w loop

            end do ! end w_sym loop

            ! save Fock matrix element

            if ( m_class == 3 ) then

              m_i                         = trans_%class_to_irrep_map(m)-&
                                          & ndocpi_(m_sym)-nactpi_(m_sym)
              fock_a_%ext(m_sym)%val(m_i) = val

            else

              m_i                             = trans_%class_to_irrep_map(m)
              n_i                             = trans_%class_to_irrep_map(n)
              fock_a_%occ(m_sym)%val(m_i,n_i) = val

            endif

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
  
                mw_sym     = group_mult_tab_(m_sym,w_sym)
                sym_offset = ints_%offset(mw_sym)

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
                    int_ind    = sym_offset + pq_index(mw,vn)
                    ival       = ival - 0.5_wp * int2(int_ind)

                    ! contract with integral/density matrix elements
                    val        = val + dval * ival

                  end do ! end v loop

                end do ! end w loop
  
              end do ! end w_sym loop

              ! save Fock matrix element

              m_i = trans_%class_to_irrep_map(m)
              n_i = trans_%class_to_irrep_map(n)

              fock_a_%occ(m_sym)%val(m_i,n_i) = val

            end do ! end n loop

          end do ! and n_class loop

        end do ! end m loop

      end do ! end m_sym loop

    end do ! end m_class loop
   
    return
  end subroutine compute_f_a

  subroutine compute_f_i_df_coulomb(int1,int2)

    implicit none

    real(wp), intent(in) :: int1(:),int2(:)

    integer :: i_sym,p_sym
    integer :: i,p,q,p_i,q_i,q_min,q_max
    integer :: pq_int
    integer :: p_class,q_class
    integer :: idf,pdf,qdf
    integer(ip) :: ii_df,pq_df,nQ
    real(wp) :: val,ddot

    nQ = int(df_vars_%nQ,kind=ip)

    ! initialize
    do p_sym = 1 , nirrep_

      if ( allocated( fock_i_%occ(p_sym)%val ) ) fock_i_%occ(p_sym)%val = 0.0_wp
      if ( allocated( fock_i_%ext(p_sym)%val ) ) fock_i_%ext(p_sym)%val = 0.0_wp

    end do

    qint_%tuQ(1)%val(:,1) = 0.0_wp

    ! calculate intermediate

    do i_sym = 1 , nirrep_

      do i = first_index_(i_sym,1) , last_index_(i_sym,1)

        idf   = df_vars_%class_to_df_map(i)

        ii_df = df_pq_index(idf,idf)

        call daxpy(nQ,2.0_wp,int2(ii_df+1:ii_df+nQ),1,qint_%tuQ(1)%val(:,1),1)

      end do ! end i loop

    end do ! end i_sym loop

    ! calculate Fock matrix element

    do p_class = 1 , 3

      do q_class = 1 , 3

        if ( q_class > p_class ) cycle

        do p_sym = 1 , nirrep_

          do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

            pdf = df_vars_%class_to_df_map(p)

            q_max = last_index_(p_sym,q_class)

            q_min = first_index_(p_sym,q_class)

            if ( p_class == q_class ) q_max = p

            if ( q_class == 3 ) q_min = p

            do q = q_min , q_max

              qdf    = df_vars_%class_to_df_map(q)

              pq_df  = df_pq_index(pdf,qdf)

              pq_int = ints_%gemind(p,q) 

              val = int1(pq_int) + & 
                  & ddot(nQ,int2(pq_df+1:pq_df+nQ),1,qint_%tuQ(1)%val(:,1),1)

              if ( ( p_class == 3 ) .and. ( q_class == 3 ) ) then

                p_i = trans_%class_to_irrep_map(p)-ndocpi_(p_sym)-nactpi_(p_sym)
                fock_i_%ext(p_sym)%val(p_i) = val

              else

                p_i = trans_%class_to_irrep_map(p)
                q_i = trans_%class_to_irrep_map(q)
                fock_i_%occ(p_sym)%val(p_i,q_i) = val

              endif

            end do ! end q loop

          end do ! end p loop

        end do ! end p_sym loop     

      end do ! end q_class loop

    end do ! end p_class loop

    return

  end subroutine compute_f_i_df_coulomb

  subroutine compute_f_i_df_exchange(int1,int2)
    ! this subroutine computes the "inactive fock matrix" according to Eq. 12.5.12 in Helgaker on page 622
    ! F_i(m|n) = h(m|n) + SUM[i \in D] { 2 g(mn|ii) - g(mi|in) }
    ! using 3-index 2-e integrals
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes. 
    implicit none
    real(wp), intent(in) :: int1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,i,i_sym,mdf,ndf,idf
    integer :: n_first,m_i,n_i
    integer(ip) :: mi,in,nQ
    real(wp) :: val,int_val,ddot

    nQ = int(df_vars_%nQ,kind=ip) 

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! orbital index in df order
          mdf = df_vars_%class_to_df_map(m)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          n_first = first_index_(m_sym,m_class)

          ! only compute diagonal elements for external-external block

          if ( m_class == 3 ) n_first = m

#ifdef OMP
!$omp parallel shared(m,mdf,first_index_,last_index_,fock_i_,int2,int1,  &
!$omp ints_,df_vars_,m_class,m_sym,n_first,trans_) num_threads(nthread_use_)
!$omp do private(n,val,ndf,idf,mi,in,m_i,n_i,i,i_sym)
#endif

          do n = n_first , m

            ! orbital index in df order
            ndf = df_vars_%class_to_df_map(n)

            val = 0.0_wp

            ! loop over irreps for i

            do i_sym = 1 , nirrep_

              ! loop over i indeces

              do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                ! orbital index in df order
                idf        = df_vars_%class_to_df_map(i) 

                ! geminal indeces in df order
                mi         = df_pq_index(mdf,idf) 
                in         = df_pq_index(idf,ndf)

                ! 2-e exchange contribution - g(mi|in)

                val        = val -  ddot(df_vars_%nQ,int2(mi+1:mi+nQ),1,int2(in+1:in+nQ),1)

              end do ! end i loop

            end do ! end i_sym loop

            ! save Fock matrix element

            if ( m_class == 3 ) then

              m_i = trans_%class_to_irrep_map(m)-ndocpi_(m_sym)-nactpi_(m_sym)
              fock_i_%ext(m_sym)%val(m_i) = fock_i_%ext(m_sym)%val(m_i) + val

            else

              m_i = trans_%class_to_irrep_map(m)
              n_i = trans_%class_to_irrep_map(n)
              fock_i_%occ(m_sym)%val(m_i,n_i) = fock_i_%occ(m_sym)%val(m_i,n_i) + val

            endif

          end do ! end n loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

#ifdef OMP
!$omp parallel shared(n_class,m,mdf,first_index_,last_index_,fock_i_,int2,int1, &
!$omp ints_,df_vars_,m_class,m_sym,trans_) num_threads(nthread_use_)
!$omp do private(n,val,ndf,idf,mi,in,m_i,n_i,i_sym,i)
#endif
            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              ! orbital index in df order
              ndf     = df_vars_%class_to_df_map(n)

              val = 0.0_wp

              ! loop over irreps for i

              do i_sym = 1 , nirrep_

                ! loop over i indeces

                do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                  ! orbital index in df order
                  idf        = df_vars_%class_to_df_map(i)

                  ! geminal indeces in df order
                  mi         = df_pq_index(mdf,idf) 
                  in         = df_pq_index(idf,ndf) 
 
                  ! 2-e exchange contribution - g(mi|in)

                  val        = val - ddot(df_vars_%nQ,int2(mi+1:mi+nQ),1,int2(in+1:in+nQ),1)

                end do ! end i loop

              end do ! end i_sym loop

              ! save Fock matrix element

              m_i              = trans_%class_to_irrep_map(m)
              n_i              = trans_%class_to_irrep_map(n)

              fock_i_%occ(m_sym)%val(m_i,n_i) = fock_i_%occ(m_sym)%val(m_i,n_i) + val

            end do ! end n_loop

#ifdef OMP
!$omp end do nowait
!$omp end parallel
#endif

          end do ! end n_class loop

        end do ! end m_loop

      end do ! end m_sym loop

    end do ! end m_class loop

    return

  end subroutine compute_f_i_df_exchange

  subroutine compute_f_i(int1,int2)
    ! this subroutine computes the "inactive fock matrix" according to Eq. 12.5.12 in Helgaker on page 622
    ! F_i(m|n) = h(m|n) + SUM[i \in D] { 2 g(mn|ii) - g(mi|in) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes. 
    implicit none
    real(wp), intent(in) :: int1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,i,i_sym,mi_sym
    integer :: mn,ii,mi,in,m_i,n_i,n_first
    integer :: int_ind,sym_offset
    real(wp) :: val

    ! initialize
    do m_sym = 1 , nirrep_

      if ( allocated( fock_i_%occ(m_sym)%val ) )     fock_i_%occ(m_sym)%val = 0.0_wp
      if ( allocated( fock_i_%ext(m_sym)%val ) ) fock_i_%ext(m_sym)%val = 0.0_wp

    end do

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)


          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          n_first = first_index_(m_sym,m_class)

          ! only compute diagonal elements for external-external block

          if ( m_class == 3 ) n_first = m

          do n = n_first , m

            ! mn-geminal index
            mn = ints_%gemind(m,n)
            ! 1-e contribution h(m|n)
            val = int1(mn)

            ! loop over irreps for i

            do i_sym = 1 , nirrep_

              mi_sym     = group_mult_tab_(m_sym,i_sym)
              sym_offset = ints_%offset(mi_sym)

              ! loop over i indeces

              do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                ! 2-e coulomb contribution 2 g(mn|ii)
                ii         = ints_%gemind(i,i)
                int_ind    = pq_index(ii,mn)
                val        = val + 2.0_wp * int2(int_ind)

                ! 2-e exchange contribution - g(mi|in)
                mi         = ints_%gemind(m,i)
                in         = ints_%gemind(i,n)
                int_ind    = sym_offset + pq_index(mi,in)
                val        = val - int2(int_ind) 

              end do ! end i loop

            end do ! end i_sym loop

            ! save Fock matrix element

            if ( m_class == 3 ) then

              m_i                             = trans_%class_to_irrep_map(m)-ndocpi_(m_sym)-nactpi_(m_sym)
              fock_i_%ext(m_sym)%val(m_i)     = val

            else

              m_i                             = trans_%class_to_irrep_map(m)
              n_i                             = trans_%class_to_irrep_map(n)
              fock_i_%occ(m_sym)%val(m_i,n_i) = val

            endif

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

                mi_sym     = group_mult_tab_(m_sym,i_sym)
                sym_offset = ints_%offset(mi_sym)

                ! loop over i indeces

                do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                  ! 2-e coulomb contribution 2 g(mn|ii)
                  ii         = ints_%gemind(i,i)
                  int_ind    = pq_index(ii,mn)
                  val        = val + 2.0_wp * int2(int_ind)

                  ! 2-e exchange contribution - g(mi|in)
                  mi         = ints_%gemind(m,i)
                  in         = ints_%gemind(i,n)
                  int_ind    = sym_offset + pq_index(mi,in)
                  val        = val - int2(int_ind)

                end do ! end i loop
 
              end do ! end i_sym loop

              ! save Fock matrix element

              m_i              = trans_%class_to_irrep_map(m)
              n_i              = trans_%class_to_irrep_map(n)

              fock_i_%occ(m_sym)%val(m_i,n_i) = val

            end do ! end n_loop

          end do ! end n_class loop

        end do ! end m_loop

      end do ! end m_sym loop
  
    end do ! end m_class loop

    return
  end subroutine compute_f_i

  subroutine print_f_matrix(mat)

    implicit none

    type(fock_info) :: mat

    integer :: m_sym,n_class,m_class,m,n,m_i,n_i
    integer :: iprint

    do m_sym = 1 , nirrep_

      if ( ndocpi_(m_sym) + nactpi_(m_sym) > 0 ) then !allocated ( mat%occ(m_sym)%val ) ) then

        write(fid_,'(a,1x,i1,2(i3,1x))')'lower-triangular elements for irrep',m_sym,&
              & ndocpi_(m_sym) + nactpi_(m_sym) + nextpi_(m_sym) , ndocpi_(m_sym) + nactpi_(m_sym)

        iprint=0
 
        do m_class = 1 , 3 
 
          do n_class = 1 , 2

            do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class) 

              do n = first_index_(m_sym,n_class) , min( m,last_index_(m_sym,n_class) )

                iprint = iprint + 1

                m_i              = trans_%class_to_irrep_map(m)
                n_i              = trans_%class_to_irrep_map(n)
 
                write(fid_,'(2(i3,1x),es20.13,1x)',advance='no')m_i,n_i,mat%occ(m_sym)%val(m_i,n_i)

                if ( iprint < 4 ) cycle

                iprint = 0
                write(fid_,*)

              end do

            end do
  
          end do

        end do

        if ( iprint > 0 ) write(fid_,*)

      end if

      if ( nextpi_(m_sym) > 0 ) then

         write(fid_,'(a,1x,i1,2(i3,1x))')'diagonal elements for irrep',m_sym,nextpi_(m_sym)

         iprint = 0

         do m = first_index_(m_sym,3) , last_index_(m_sym,3)
           
            iprint = iprint + 1
 
            m_i = trans_%class_to_irrep_map(m)-ndocpi_(m_sym)-nactpi_(m_sym)

            write(fid_,'(2(i3,1x),es20.13,1x)',advance='no')m_i,n_i,mat%ext(m_sym)%val(m_i)

            if ( iprint < 4 ) cycle

            iprint = 0
            write(fid_,*)

         end do

         if ( iprint > 0 ) write(fid_,*)

      end if

    end do

    return

  end subroutine print_f_matrix

  subroutine print_q_matrix(mat)

    implicit none

    real(wp), intent(in) :: mat(:,:)
 
    integer :: m,n,m_i,n_i,m_class,m_sym,iprint

    do m_sym = 1 , nirrep_

      if ( nactpi_(m_sym) == 0 ) cycle

      write(fid_,'(a,1x,i1,2(i3,1x))')'matrix elements for irrep',m_sym,&
         & nactpi_(m_sym) , ndocpi_(m_sym) + nactpi_(m_sym) + nextpi_(m_sym)

      iprint = 0

      do m_class = 1 , 3

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          do n = first_index_(m_sym,2) , last_index_(m_sym,2)

            iprint = iprint + 1

            m_i              = trans_%class_to_irrep_map(m)
            n_i              = trans_%class_to_irrep_map(n)

            write(fid_,'(2(i3,1x),es20.13,1x)',advance='no')m_i,n_i,mat(m-ndocpi_(m_sym),n)

            if ( iprint < 4 ) cycle

            iprint = 0
            write(fid_,*)

          end do

        end do

      end do

      if ( iprint > 0 ) write(fid_,*)

    end do

    return

  end subroutine print_q_matrix

  subroutine allocate_temporary_fock_matrices()
    implicit none
    integer :: i_sym,ntot,nocc
    ! the total number of nonzero LT elements in the active/inactive Fock matrix 
    ! is equal to the number of geminals in the totally symmetric irrep

    ! allocate Fock matrices for the external block, only the diagonal elements are needed)

    allocate(fock_i_%occ(nirrep_),fock_i_%ext(nirrep_))

    do i_sym = 1 , nirrep_

      nocc = nactpi_(i_sym)+ndocpi_(i_sym)
      ntot = nocc + nextpi_(i_sym)

      allocate(fock_i_%occ(i_sym)%val(ntot,nocc)) 

      allocate(fock_i_%ext(i_sym)%val(nextpi_(i_sym)))

    end do

    allocate(fock_a_%occ(nirrep_),fock_a_%ext(nirrep_))

    do i_sym = 1 , nirrep_

      nocc = nactpi_(i_sym)+ndocpi_(i_sym)
      ntot = nocc + nextpi_(i_sym)

      allocate(fock_a_%occ(i_sym)%val(ntot,nocc))

      allocate(fock_a_%ext(i_sym)%val(nextpi_(i_sym)))

    end do

    ! the row index here will be between 1-nact, while the column index is arbitrary

    allocate(q_(nact_tot_,nmo_tot_))

    allocate(z_(nact_tot_,nmo_tot_))

    return
  end subroutine allocate_temporary_fock_matrices

  subroutine deallocate_temporary_fock_matrices()
    implicit none

    integer :: i_sym

    if ( allocated(fock_i_%occ) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(fock_i_%occ(i_sym)%val) ) cycle

        deallocate(fock_i_%occ(i_sym)%val)

      end do

      deallocate(fock_i_%occ)

    endif

    if ( allocated(fock_i_%ext) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(fock_i_%ext(i_sym)%val) ) cycle

        deallocate(fock_i_%ext(i_sym)%val)

      end do

      deallocate(fock_i_%ext)

    end if

    if ( allocated(fock_a_%occ) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(fock_a_%occ(i_sym)%val) ) cycle

        deallocate(fock_a_%occ(i_sym)%val)

      end do

      deallocate(fock_a_%occ)

    endif

    if ( allocated(fock_a_%ext) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(fock_a_%ext(i_sym)%val) ) cycle

        deallocate(fock_a_%ext(i_sym)%val)

      end do

      deallocate(fock_a_%ext)

    end if

    if (allocated(q_))        deallocate(q_)
    if (allocated(z_))        deallocate(z_)
    return
  end subroutine deallocate_temporary_fock_matrices

  subroutine allocate_qint()
    implicit none

    integer :: i_sym

    allocate(qint_%tuQ(nirrep_))

    do i_sym = 1 , nirrep_

      allocate(qint_%tuQ(i_sym)%val(df_vars_%nQ,df_vars_%noccgempi(i_sym)))

    end do

    return
  end subroutine allocate_qint

  subroutine deallocate_qint()
    implicit none

    integer :: i_sym

    if ( allocated(qint_%tuQ) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(qint_%tuQ(i_sym)%val) ) cycle

        deallocate(qint_%tuQ(i_sym)%val)

      end do

      deallocate(qint_%tuQ)

    end if

  end subroutine deallocate_qint

end module focas_gradient
