module focas_transform_oeints

  use focas_data

  implicit none

  contains

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

        irrep_block_2(1:num_i,1:num_i) = matmul(irrep_block_1(1:num_i,1:num_i),trans_%u_irrep_block(i_sym)%val)

        irrep_block_1(1:num_i,1:num_i) = matmul(transpose(trans_%u_irrep_block(i_sym)%val),irrep_block_2(1:num_i,1:num_i))

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
