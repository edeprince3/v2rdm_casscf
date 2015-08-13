module focas_hessian

  use focas_data

  implicit none

  contains

    subroutine diagonal_inverse_hessian_preconditioner(grad,f_i,f_a,q,z,den1)
      implicit none
      real(wp), intent(inout) :: grad(:)
      real(wp), intent(in) :: f_i(:),f_a(:),den1(:)
      real(wp), intent(in) :: q(:,:),z(:,:)
      integer :: error

      ! active-doubly occupied pairs

      if ( rot_pair_%n_ad > 0 ) error=diagonal_hessian_ad(grad,f_i,f_a,q,z,den1)

      ! active-active pairs

      if ( rot_pair_%n_aa > 0 ) error=diagonal_hessian_aa(grad)

      ! external-doubly occupied pairs

      if ( rot_pair_%n_ed > 0 ) error=diagonal_hessian_ed(grad,f_i,f_a)

      ! external-active pairs

      if ( rot_pair_%n_ea > 0 ) error=diagonal_hessian_ea(grad,f_i,f_a,q,z,den1)

      return
    end subroutine diagonal_inverse_hessian_preconditioner

    integer function diagonal_hessian_ad(grad,f_i,f_a,q,z,den1)
      implicit none
      ! subroutine to compute the diagonal Hessian matrix element H(ti|ti) according to Eq. 4.7c of
      ! Chaban, Schmidt, Gordon, Theor. Chem. Acc., 97, 88-95 (1997) 
      ! H(ti|ti) = 2 * [ 2 * { ( f_i(t,t) + f_a(t,t) ) - ( f_i(i,i) + f_a(i,i) ) }  + d(t,t) * ( f_i(i,i) + f_a(i,i) ) - q(t,t) -z(t) ]
      ! where z(t,t) = sum_{u} [ d(t,u) * f_i(t,u) ] ] && a \in ext & t,u \in act
      real(wp), intent(inout) :: grad(:) 
      real(wp), intent(in) :: f_i(:),f_a(:),q(:,:),z(:,:),den1(:)

      integer :: grad_ind,t,i,tt_int,ii_int,tt_den,t_sym
      real(wp) :: h_val,fac

      ! loop over all active - doubly-occupied pairs

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

        do i = first_index_(t_sym,1) , last_index_(t_sym,1)

          do t = first_index_(t_sym,2) , last_index_(t_sym,2)

            ! look up integral/density indeces
            ii_int  = ints_%gemind(i,i)
            tt_int  = ints_%gemind(t,t)
            tt_den  = dens_%gemind(t,t)

            ! update gradient index
            grad_ind = grad_ind + 1
 
            ! compute Hessian value
            h_val    = 2.0_wp * (                                                                 &
                       2.0_wp * ( ( f_i(tt_int) + f_a(tt_int) ) - ( f_i(ii_int) + f_a(ii_int) ) ) &
                     + den1(tt_den) * ( f_i(ii_int) + f_a(ii_int) ) - q(t-ndoc_tot_,t) - z(t-ndoc_tot_,t) )
            
            fac = 0.0_wp
            if ( h_val > 0.0_wp ) fac = 1.0_wp/h_val

            grad(grad_ind) = fac * grad(grad_ind)

          end do

        end do

      end do

      diagonal_hessian_ad = 0

      return

    end function diagonal_hessian_ad

    integer function diagonal_hessian_aa(grad)
      implicit none
      ! subroutine to compute the diagonal Hessian matrix element H(tu|tu)
      ! currently, this is not implemented, so just set to 1.0 for now
      real(wp),intent(inout) :: grad(:)
      integer :: t,u,t_sym,grad_ind
      real(wp) :: h_val,fac

      ! loop over all active - active pairs

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1

            ! compute Hessian value
            h_val    = 50.0_wp

            fac = 0.0_wp
            if ( h_val > 0.0_wp ) fac = 1.0_wp/h_val

            grad(grad_ind) = fac * grad(grad_ind)

          end do
      
        end do

      end do

      diagonal_hessian_aa = 0

      return

    end function diagonal_hessian_aa

    integer function diagonal_hessian_ed(grad,f_i,f_a)
      implicit none
      ! subroutine to compute the diagonal Hessian matrix element H(ia|ia)according to Eq. 4.7a of
      ! Chaban, Schmidt, Gordon, Theor. Chem. Acc., 97, 88-95 (1997) 
      ! H(ai|ai) = 4 * (f_i(a,a) + f_a(a,a) ) - 4 * ( f_i(i,i) + f_a(i,i)) 
      ! a \in ext & i \in act
      real(wp), intent(inout) :: grad(:)
      real(wp), intent(in) :: f_i(:),f_a(:)

      integer :: a,i,grad_ind,aa,ii,a_sym
      real(wp) :: h_val,fac

      ! loop over all external - doubly-occupied pairs

      do a_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

        do i = first_index_(a_sym,1) , last_index_(a_sym,1)

          do a = first_index_(a_sym,3) , last_index_(a_sym,3)

            ! update Hessian index

            grad_ind = grad_ind + 1

            ! inactive/active fock matrix indeces

            aa      = ints_%gemind(a,a)
            ii      = ints_%gemind(i,i)

            ! calculate diagonal Hessian element

            h_val    = 4.0 * ( f_i(aa) + f_a(aa) - f_i(ii) - f_a(ii) )

            fac = 0.0_wp
            if ( h_val > 0.0_wp ) fac = 1.0_wp/h_val

            grad(grad_ind) = fac * grad(grad_ind)

          end do

        end do
 
      end do

      diagonal_hessian_ed = 0

      return

    end function diagonal_hessian_ed

    integer function diagonal_hessian_ea(grad,f_i,f_a,q,z,den1)
      implicit none
      ! subroutine to compute the diagonal Hessian matrix element H(ai|ai) according to Eq. 4.7b of
      ! Chaban, Schmidt, Gordon, Theor. Chem. Acc., 97, 88-95 (1997) 
      ! H(at|at) = 2 * [  d(t,t) * ( f_i(a,a) + f_a(a,a) ) - q(t,t) - z(t,t) ] 
      ! where z(t,t) = sum_{u} [ d(t,u) * f_i(t,u) ]  &&  a \in ext & t,u \in act 
      real(wp), intent(inout) :: grad(:)
      real(wp), intent(in) :: f_i(:),f_a(:),z(:,:),q(:,:),den1(:)

      integer :: a,t,grad_ind,tt_den,aa_int,a_sym
      real(wp) :: h_val,fac

      ! loop over all external - doubly-occupied pairs

      do a_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_act_type)

        do t = first_index_(a_sym,2) , last_index_(a_sym,2)

          do a = first_index_(a_sym,3) , last_index_(a_sym,3)

            ! integral/density indeces
  
            tt_den  = dens_%gemind(t,t)
            aa_int  = ints_%gemind(a,a)

            ! update diagonal Hessian index
  
            grad_ind = grad_ind + 1

            ! calculate diagonal Hessian element

            h_val    = 2.0_wp * ( den1(tt_den) * ( f_i(aa_int) + f_a(aa_int) ) - q(t-ndoc_tot_,t) - z(t-ndoc_tot_,t) )

            fac = 0.0_wp
            if ( h_val > 0.0_wp ) fac = 1.0_wp/h_val

            grad(grad_ind) = fac * grad(grad_ind)

          end do

        end do

      end do

      diagonal_hessian_ea = 0

      return

    end function diagonal_hessian_ea

    subroutine allocate_hessian_data()
      implicit none
      if (allocated(diagonal_orbital_hessian_)) call deallocate_hessian_data()

      ! allocate(diagonal_orbital_hessian_(rot_pair_%n_tot))

      return
    end subroutine allocate_hessian_data

    subroutine deallocate_hessian_data
      implicit none

      if (allocated(diagonal_orbital_hessian_)) deallocate(diagonal_orbital_hessian_)
      return
    end subroutine deallocate_hessian_data

end module focas_hessian
