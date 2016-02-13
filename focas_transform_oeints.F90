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

module focas_transform_oeints

  use focas_data

  implicit none

  contains

    integer function transform_mocoeff(mo_coeff)

      implicit none
      real(wp) :: mo_coeff(:,:)
      integer :: i_sym,num_i,i_offset,max_nmopi,first_i,last_i
      real(wp), allocatable :: block_tmp(:,:)

      transform_mocoeff = 1

      if ( size(mo_coeff,dim = 1) /= nmo_tot_ ) then
        if ( log_print_ == 1 ) then
          write(fid_,'(a)')'dimension of mo_coeff(:,:) and nmo_tot_ do not match'
        end if
        return
      end if

      max_nmopi = maxval(trans_%nmopi)
   
      allocate(block_tmp(max_nmopi,max_nmopi))

      do i_sym = 1 , nirrep_

        num_i    = trans_%nmopi(i_sym)

        if ( num_i == 0 ) cycle

        i_offset = trans_%offset(i_sym)

        first_i  = i_offset + 1
        last_i   = i_offset + num_i

#ifdef BLAS

        call dgemm('N','N',num_i,num_i,num_i,1.0_wp,mo_coeff(first_i:last_i,first_i:last_i),num_i,&
                     trans_%u_irrep_block(i_sym)%val,num_i,0.0_wp,block_tmp,max_nmopi)
        call dgemm('T','N',num_i,num_i,num_i,1.0_wp,trans_%u_irrep_block(i_sym)%val,num_i, &
                     block_tmp,max_nmopi,0.0_wp,mo_coeff(first_i:last_i,first_i:last_i),num_i)

#else

        block_tmp(1:num_i,1:num_i) = matmul(mo_coeff(first_i:last_i,first_i:last_i),&
                                          & trans_%u_irrep_block(i_sym)%val)

        mo_coeff(first_i:last_i,first_i:last_i) = matmul(transpose(trans_%u_irrep_block(i_sym)%val), &
                                                       & block_tmp(1:num_i,1:num_i))

#endif

      end do

      deallocate(block_tmp)

      transform_mocoeff = 0

      return

    end function transform_mocoeff

    integer function transform_oeints(int1)
      implicit none
  
      real(wp) :: int1(:)
      real(wp), allocatable :: irrep_block_1(:,:),irrep_block_2(:,:)
 
      integer :: max_nmopi,ij_class,i_offset,num_i,i_class,j_class,i_irrep,j_irrep,i_sym

      ! initialize return value

      transform_oeints = 1

      ! determine maximum block size and allocate temporary matrix
 
      max_nmopi = maxval(trans_%nmopi)
 
      allocate(irrep_block_1(max_nmopi,max_nmopi),irrep_block_2(max_nmopi,max_nmopi))

      ! loop over irreps for i

      do i_sym = 1 , nirrep_

        ! initialize temporary matrix

        irrep_block_1 = 0.0_wp

        ! number of orbitals in this block as well as offset for indexing

        num_i    = trans_%nmopi(i_sym)

        if ( num_i == 0 ) cycle

        i_offset = trans_%offset(i_sym)

        ! **************
        ! *** GATHER ***
        ! ************** 

        ! loop over symmetry-reduced indeces for i 

        do i_irrep = 1 , num_i

          ! class-reduced index for i needed for integral addressing

          i_class = trans_%irrep_to_class_map(i_irrep+i_offset)

          ! loop over symmetry-reduced indeces for j

          do j_irrep = 1 , i_irrep

            ! class reduced index for j needed for integral addressing

            j_class  = trans_%irrep_to_class_map(j_irrep+i_offset)

            ! one-electron index

            ij_class = ints_%gemind(i_class,j_class) 

            ! save integral making sure that the current irrep block is symmetric

            irrep_block_1(j_irrep,i_irrep) = int1(ij_class)
            irrep_block_1(i_irrep,j_irrep) = irrep_block_1(j_irrep,i_irrep)

          end do ! end j_irrep loop

        end do ! end i_irrep loop

        ! *****************
        ! *** TRANSFORM ***
        ! *****************

!#ifdef BLAS
!
        call dgemm('N','N',num_i,num_i,num_i,1.0_wp,irrep_block_1,max_nmopi,&
                     trans_%u_irrep_block(i_sym)%val,num_i,0.0_wp,irrep_block_2,max_nmopi)
        call dgemm('T','N',num_i,num_i,num_i,1.0_wp,trans_%u_irrep_block(i_sym)%val,num_i, &
                     irrep_block_2,max_nmopi,0.0_wp,irrep_block_1,max_nmopi)
!
!#else
!
!        irrep_block_2(1:num_i,1:num_i) = matmul(irrep_block_1(1:num_i,1:num_i),trans_%u_irrep_block(i_sym)%val)
!
!        irrep_block_1(1:num_i,1:num_i) = matmul(transpose(trans_%u_irrep_block(i_sym)%val),irrep_block_2(1:num_i,1:num_i))
!
!#endif

        ! ***************
        ! *** SCATTER ***
        ! ***************

        do i_irrep = 1 , num_i

          ! class-reduced index for i needed for integral addressing

          i_class = trans_%irrep_to_class_map(i_irrep+i_offset)

          ! loop over symmetry-reduced indeces for j

          do j_irrep = 1, i_irrep

            ! class reduced index for j needed for integral addressing

            j_class  = trans_%irrep_to_class_map(j_irrep+i_offset)

            ! one-electron index

            ij_class = ints_%gemind(i_class,j_class)

            ! save integral

            int1(ij_class) = irrep_block_1(j_irrep,i_irrep)
 
          end do ! end j_irrep loop

        end do ! end i_irrep loop

      end do ! end i_sym loop

      ! deallocate temporary matrices

      deallocate(irrep_block_1,irrep_block_2)

      ! no errors encountered

      transform_oeints = 0

      return
  
    end function transform_oeints

end module focas_transform_oeints
