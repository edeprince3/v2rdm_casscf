module focas_transform_driver

  use focas_data
  use focas_transform_teints
  use focas_transform_oeints

  implicit none

  contains

    subroutine transform_driver(int1,int2,mo_coeff)
      implicit none
      real(wp) :: int2(:),int1(:),mo_coeff(:,:)
      integer :: error

      ! 1-e integrals
      error = transform_oeints(int1)
      if ( error /= 0 ) stop   

      ! 2-e integrals
      if ( df_vars_%use_df_teints == 0 ) then
        error = transform_teints(int2)
      else
        error = transform_teints_df(int2)
      end if
      if ( error /= 0 ) stop      

      ! mo_coeff matrix
      error = transform_mocoeff(mo_coeff)
      if ( error /= 0 ) stop

      return
    end subroutine transform_driver

    subroutine allocate_transformation_matrices()
      implicit none

      integer :: nmo_i_sym,i_sym,i_class

      ! figure out the total number of mos per irrep

      allocate(trans_%npairpi(nirrep_))
      allocate(trans_%nmopi(nirrep_))

      ! figure out number of orbitals per irrep

      trans_%nmopi = 0

      ! loop over orbital symmetries

      do i_sym = 1 , nirrep_
 
        ! loop over orbital classes

        do i_class = 1 , 3

          trans_%nmopi(i_sym) = trans_%nmopi(i_sym) + ( last_index_(i_sym,i_class) - first_index_(i_sym,i_class) ) + 1

        end do ! end i_class loop
 
      end do ! end i_sym loop

      ! allocate the remaining indexing arrays

      allocate(trans_%offset(nirrep_))
      allocate(trans_%irrep_to_class_map(nmo_tot_))
      allocate(trans_%class_to_irrep_map(nmo_tot_))

      ! allocate matrices for symmetry blocks of transformation array

      allocate(trans_%u_irrep_block(nirrep_))

      ! loop over symmetry blocks

      do i_sym = 1 , nirrep_

        ! number of orbitals in this block

        nmo_i_sym = trans_%nmopi(i_sym)

        ! allocate this block

        allocate(trans_%u_irrep_block(i_sym)%val(nmo_i_sym,nmo_i_sym))

      end do

      return

    end subroutine allocate_transformation_matrices

    subroutine determine_transformation_maps()
      implicit none
      ! subroutine to figure out mapping arrays between symmetry order and class order
 
      integer :: i,i_sym,i_class,i_sym_off,i_in_irrep

      ! figure out offset array for each symmetry block

      trans_%offset(1) = 0

      ! loop over remaining symmetries for i
 
      do i_sym = 2 , nirrep_    

        trans_%offset(i_sym) = trans_%offset(i_sym-1) + trans_%nmopi(i_sym-1)

      end do ! end i_sym loop

      ! now, figure out mapping arrays

      ! loop over symmetries for i

      do i_sym = 1 , nirrep_
 
        ! initialize reduced index counter

        i_in_irrep = 0
  
        ! symmetry offset
 
        i_sym_off = trans_%offset(i_sym)
          
        ! loop over orbital types for i     

        do i_class = 1 , 3

          ! loop over i indeces

          do i = first_index_(i_sym,i_class) , last_index_(i_sym,i_class)

            ! update reduced index counter

            i_in_irrep = i_in_irrep + 1

            ! save mappings     
 
            trans_%irrep_to_class_map(i_in_irrep+i_sym_off) = i  
            trans_%class_to_irrep_map(i)                    = i_in_irrep

          end do ! end i loop

        end do ! end i_class loop

      end do ! end i_sym loop

      return
    end subroutine determine_transformation_maps

    subroutine deallocate_transformation_matrices()
      implicit none
      integer :: i_sym
      if (allocated(trans_%npairpi))            deallocate(trans_%npairpi)
      if (allocated(trans_%nmopi))              deallocate(trans_%nmopi)
      if (allocated(trans_%offset))             deallocate(trans_%offset)
      if (allocated(trans_%irrep_to_class_map)) deallocate(trans_%irrep_to_class_map)
      if (allocated(trans_%class_to_irrep_map)) deallocate(trans_%class_to_irrep_map)
      if (allocated(trans_%u_irrep_block)) then

        ! loop over irreps for i

        do i_sym = 1 , nirrep_

          if ( .not.allocated(trans_%u_irrep_block(i_sym)%val) ) cycle

          deallocate(trans_%u_irrep_block(i_sym)%val)

        end do ! end i_sym loop

        deallocate(trans_%u_irrep_block)

      endif 

      return
    end subroutine deallocate_transformation_matrices

end module focas_transform_driver
