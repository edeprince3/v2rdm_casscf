module focas_transform_teints
  
  use focas_data

  implicit none
  
  contains

    integer function transform_teints_df(int2)
      implicit none

      real(wp) :: int2(:)

      type tmp_matrix
        real(wp), allocatable :: int_tmp(:)
        real(wp), allocatable :: tmp(:,:)
        real(wp), allocatable :: mat(:,:)
      end type tmp_matrix

      type(tmp_matrix), allocatable :: aux(:)

      integer :: i_thread,Q,row,col,row_sym,col_sym,int_ind,max_nmopi
      integer :: col_min,col_max,row_min,row_max,ncol,nrow,row_off,col_off
      integer :: first_Q(nthread_use_),last_Q(nthread_use_),omp_get_thread_num

      integer :: tmpp(nthread_use_)

      max_nmopi = maxval(trans_%nmopi)

      transform_teints_df = allocate_tmp_matrices()

      transform_teints_df = setup_Q_bounds()

      ! loop over threads

#ifdef OMP
!$omp parallel shared(first_Q,last_Q,int2,df_vars_,nirrep_) num_threads(nthread_use_)
!$omp do private(i_thread,Q,col,row,int_ind,col_sym,row_sym,col_min,col_max,row_min,row_max,ncol,nrow,col_off,row_off)
#endif
      do i_thread = 1 , nthread_use_

        tmpp(i_thread)=omp_get_thread_num()

        ! loop over Q indeces to be transformed by this thread
 
        do Q = first_Q(i_thread) , last_Q(i_thread)

          ! initialize temporary matrix

          aux(i_thread)%mat = 0.0_wp

          ! *************************************************************
          ! *** GATHER ( only LT row > col elements are accessed in int2)
          ! *************************************************************

          ! loop over column indeces

          int_ind = Q

          do col = 1 , nmo_tot_

            ! loop over row indeces
 
            do row = 1 , col

              ! save off-diagonal elements
 
              aux(i_thread)%mat(col,row) = int2(int_ind)

              int_ind = int_ind + df_vars_%nQ

            end do ! end row loop 

          end do ! end col loop

          ! **************************************************************
          ! *** TRANSFORM ( only lower triangular blocks are transformed )
          ! **************************************************************

          ! loop over orbital symmetries in the columns

          do col_sym = 1 , nirrep_
 
            ncol    = trans_%nmopi(col_sym)
  
            if ( ncol == 0 ) cycle

            col_off = trans_%offset(col_sym)

            col_min = col_off + 1
            
            col_max = col_off + ncol

            ! symmetrize the diagonal block since only the LT elements were saved in the above gather

            call symmetrize_diagonal_block(aux(i_thread)%mat(col_min:col_max,col_min:col_max),ncol)

#ifdef BLAS
            call dgemm('N','N',ncol,ncol,ncol,1.0_wp,aux(i_thread)%mat(col_min:col_max,col_min:col_max),ncol,&
                     & trans_%u_irrep_block(col_sym)%val,ncol,0.0_wp,aux(i_thread)%tmp,max_nmopi)
            call dgemm('T','N',ncol,ncol,ncol,1.0_wp,trans_%u_irrep_block(col_sym)%val,ncol, &
                     & aux(i_thread)%tmp,max_nmopi,0.0_wp,aux(i_thread)%mat(col_min:col_max,col_min:col_max),ncol)
#else
            aux(i_thread)%tmp(1:ncol,1:ncol) = matmul(aux(i_thread)%mat(col_min:col_max,col_min:col_max),&
                                             & trans_%u_irrep_block(col_sym)%val)
            aux(i_thread)%mat(col_min:col_max,col_min:col_max) = matmul(transpose(trans_%u_irrep_block(col_sym)%val),&
                                                               & aux(i_thread)%tmp(1:ncol,1:ncol))
#endif

            do row_sym = col_sym + 1 , nirrep_

              nrow    = trans_%nmopi(row_sym)

              if ( nrow == 0 ) cycle

              row_off = trans_%offset(row_sym)

              row_min = row_off + 1
 
              row_max = row_off + nrow

#ifdef BLAS
              call dgemm('N','N',nrow,ncol,ncol,1.0_wp,aux(i_thread)%mat(row_min:row_max,col_min:col_max),nrow,&
                         trans_%u_irrep_block(col_sym)%val,ncol,0.0_wp,aux(i_thread)%tmp,max_nmopi)
              call dgemm('T','N',nrow,ncol,nrow,1.0_wp,trans_%u_irrep_block(row_sym)%val,nrow, &
                         aux(i_thread)%tmp,max_nmopi,0.0_wp,aux(i_thread)%mat(row_min:row_max,col_min:col_max),nrow)
#else
              aux(i_thread)%tmp(1:nrow,1:ncol) = matmul(aux(i_thread)%mat(row_min:row_max,col_min:col_max),&
                                               & trans_%u_irrep_block(col_sym)%val)
              aux(i_thread)%mat(row_min:row_max,col_min:col_max) = matmul(transpose(trans_%u_irrep_block(row_sym)%val),&
                                                                 & aux(i_thread)%tmp(1:nrow,1:ncol))
#endif


            end do ! end row_sym loop

          end do ! end col_sym loop

          ! *************************************************************
          ! *** SCATTER (only LT row > col elements are accessed in int2)
          ! *************************************************************

          int_ind = Q

          do col = 1 , nmo_tot_

            ! loop over row indeces

            do row = 1 , col

              ! save off-diagonal elements

              int2(int_ind) = aux(i_thread)%mat(col,row)

              int_ind = int_ind + df_vars_%nQ

            end do ! end row loop 

          end do ! end col loop

        end do ! end Q loop

      end do ! end i_thread loop
#ifdef OMP
!$omp end do
!$omp end parallel
#endif

      transform_teints_df = deallocate_tmp_matrices()

      return

      contains

        subroutine symmetrize_diagonal_block(diag_block,ndim)
        
         implicit none

         ! simple function to symmetrize a square matrix for which onyl the LT elements 
         ! were stored
      
         integer, intent(in)     :: ndim
         real(wp), intent(inout) :: diag_block(ndim,ndim)

         integer :: row,col

         do col = 1 , ndim

           do row = col + 1 , ndim

             diag_block(col,row) = diag_block(row,col)

           end do 
  
         end do

         return

        end subroutine symmetrize_diagonal_block

        integer function setup_Q_bounds()

          implicit none

          ! simple function to determine first and last auxiliary function index for each thread

          integer :: nQ_ave,i

          nQ_ave = df_vars_%nQ / nthread_use_

          do i = 1 , nthread_use_

            first_Q(i) = ( i - 1 ) * nQ_ave + 1
 
            last_Q(i)  = i * nQ_ave

          end do

          last_Q(nthread_use_) = df_vars_%nQ

          setup_Q_bounds = 0

          return

        end function setup_Q_bounds

        integer function allocate_tmp_matrices()

          implicit none

          ! simple function to allocate temporary matrices for 3-index integral transformation

          integer :: i
   
          allocate(aux(nthread_use_))
         
          do i = 1 , nthread_use_

            allocate(aux(i)%mat(nmo_tot_,nmo_tot_))
            allocate(aux(i)%tmp(max_nmopi,max_nmopi))

          end do

          allocate_tmp_matrices = 0

          return

        end function allocate_tmp_matrices
     
        integer function deallocate_tmp_matrices()

          implicit none

          ! simple function to deallocate temporary matrices for 3-indes integral transformation

          integer :: i

          if ( .not. allocated(aux) ) return
 
          do i = 1 , size(aux)

            if ( allocated(aux(i)%mat) ) deallocate(aux(i)%mat)
            if ( allocated(aux(i)%tmp) ) deallocate(aux(i)%tmp)

          end do

          deallocate(aux)

          deallocate_tmp_matrices = 0

          return

        end function deallocate_tmp_matrices

    end function transform_teints_df 

    integer function transform_teints(int2)
      implicit none

      ! driver to transform integrals g(ij|kl) --> g(rs|tu)
      ! ij irreps transformed sequentially

      real(wp) :: int2(:)

      integer :: ij_sym,first_ij,last_ij
      integer :: irrep_block_transformed(nirrep_)

      ! initialize error flag

      transform_teints = 1

      first_ij = ints_%offset(1) + 1
      last_ij  = first_ij + ints_%nnzpi(1) - 1

      irrep_block_transformed(1) = transform_teints_g0_block(int2(first_ij:last_ij),1)

      do ij_sym = 2 , nirrep_

        ! first/last nnz integral index in this irrep

        first_ij = ints_%offset(ij_sym) + 1
        last_ij  = first_ij + ints_%nnzpi(ij_sym) - 1

        ! transform integrals

        irrep_block_transformed(ij_sym) = transform_teints_irrep_block(int2(first_ij:last_ij),ij_sym)

        ! something went wrong, so return with error flag

        if ( irrep_block_transformed(ij_sym) /= 0 ) return

      end do

      ! if we got this far, the transformations are done

      transform_teints = 0

      return
    end function transform_teints


    integer function transform_teints_irrep_block(int_block,ij_sym)
      implicit none
      ! subroutine to transform an irrep-block (not totally symmetric) of two-electron integrals
      ! the algorithm performs the transformation in two half transfroms.
      ! 1) The first half transform consists of two steps:
      !    a) GATHER:    for a given ij pair, collect and store all integrals (ij|kl)
      !    b) TRANSFORM: for a given ij pair, perform two dense matrix-matrix multiplications to 
      !                   transform (ij|kl) --> (ij|tu)
      !
      ! To do the second half transform, half-transformed integrals need to be stored for ALL ij and tu
      ! This means, additional storage that is roughly twice the size of the LT integral block 
      !
      ! 2) The second half transform consists of three steps
      !    a) GATHER:    for a given tu pair, collect and store all integrals (ij|tu)
      !    c) TRANSFORM: for a given ij pair, perform two dense matrix-matrix multiplications to 
      !                  transform (ij|tu) --> (rs|tu)
      !    d) SCATTER:   scatter (rs|tu) into original integral array
      ! 
      ! NOTE: Currently, the GATHER and SCATTER operations do more work than the minimum required; however,
      !       after timing the total time spent in each type of operation, it is clear that the GATHER and
      !       SCATTER operations are an almost negligible fraction of the total time and thus, further
      !       optimization will be skipped. 
      integer, intent(in) :: ij_sym
      real(wp) :: int_block(:)

      type int_scr
        real(wp), allocatable :: mat(:,:)
      end type int_scr

      type int_tmp
        type(int_scr), allocatable :: irrep(:)
      end type int_tmp

      type(int_tmp), allocatable :: A(:)

      real(wp), allocatable :: A_tilde(:,:),B_tilde(:,:)

      integer :: nnz_ij,ij,max_nmopi,error

      integer :: i_class,j_class,k_class,l_class,t_class,u_class,tu_class,kl_class,ij_class
      integer :: i_irrep,j_irrep,k_irrep,l_irrep,t_irrep,u_irrep
      integer :: i_sym,j_sym,k_sym,l_sym,t_sym,u_sym 
      integer :: num_i,num_j,num_k,num_l,num_t,num_u
      integer :: i_offset,j_offset,k_offset,l_offset,t_offset,u_offset

      integer(ip) :: ijkl_class,rstu_class

      ! initialize error flag

      transform_teints_irrep_block = 1
 
      ! the number of nnz integrals is not the same as the length of the vector int2

      if ( size(int_block,dim=1) /= ints_%nnzpi(ij_sym) ) return

      ! allocate temporary matrices

      ! maximum matrix size and number of matrices

      max_nmopi = maxval(trans_%nmopi)

      nnz_ij    = ints_%ngempi(ij_sym)

      ! allocate temporary arrays

      transform_teints_irrep_block = allocate_transform_scr()

      if ( transform_teints_irrep_block /= 0 ) stop

      ! ********************************
      ! *** 1ST HALF TRANSFORM (kl-->tu)
      ! ********************************

      ! loop over ij geminals

      do ij_class = 1 , nnz_ij

        ! loop over symmetries for k

        do k_sym = 1 , nirrep_

          ! figure out symmetry of l such that ij_sym == kl_sym

          l_sym = group_mult_tab_(ij_sym,k_sym)

          ! make sure that we are only addressing unique integrals with k_sym > l_sym
 
          if ( k_sym < l_sym) cycle

          ! *******************************************
          ! *** GATHER *** A[ij][k_sym](k,l) = g(ij|kl)
          ! *******************************************

          ! here, we make use of the fact that this matrix has at most 8 nnz blocks,
          ! with the blocks indexed by k_sym

          ! offsets for figuring out full k & l indeces

          k_offset = trans_%offset(k_sym)
          l_offset = trans_%offset(l_sym)

          ! number of orbitals with symmetry k_sym and l_sym

          num_k    = trans_%nmopi(k_sym)
          num_l    = trans_%nmopi(l_sym) 

          ! loop over symmetry-reduced k indeces

          do k_irrep = 1 , num_k

            ! class k index for integral addressing

            k_class = trans_%irrep_to_class_map(k_irrep+k_offset)

            ! loop over symmetry reduced l indeces

            do l_irrep = 1 , num_l

              ! class l index for integral addressing

              l_class    = trans_%irrep_to_class_map(l_irrep+l_offset)

              ! integral address (class order)

              kl_class   = ints_%gemind(k_class,l_class)

              ijkl_class = pq_index(ij_class,kl_class) 

              ! save integral

              A(ij_class)%irrep(k_sym)%mat(k_irrep,l_irrep) = int_block(ijkl_class) 

            end do ! end l loop

          end do ! end k loop

          !*********************************************
          ! *** TRANSFORM A[ij] = C^T * A[ij][k_sym] * C
          !*********************************************

          ! A_tilde[ij][k_sym] = A[ij][k_sym] * C[l_sym]
          ! followed by
          ! A = (C[k_sym])^T * A_tilde[]ij][k_sym]

          ! only need to do work if there are orbitals with symmetry k_sym and l_sym 

          if ( ( num_k == 0 ) .or. ( num_l == 0 ) ) cycle

# ifdef BLAS
          call dgemm('N','N',num_k,num_l,num_l,1.0_wp,A(ij_class)%irrep(k_sym)%mat,num_k,&
                     trans_%u_irrep_block(l_sym)%val,num_l,0.0_wp,A_tilde,max_nmopi)
          call dgemm('T','N',num_k,num_l,num_k,1.0_wp,trans_%u_irrep_block(k_sym)%val,num_k, &
                     A_tilde,max_nmopi,0.0_wp,A(ij_class)%irrep(k_sym)%mat,num_k)
# else
          A_tilde(1:num_k,1:num_l)       = matmul(A(ij_class)%irrep(k_sym)%mat,trans_%u_irrep_block(l_sym)%val)
          A(ij_class)%irrep(k_sym)%mat   = matmul(transpose(trans_%u_irrep_block(k_sym)%val),A_tilde(1:num_k,1:num_l))
# endif
          
        end do ! end k_sym loop

      end do ! end ij_class loop      

      ! at this point, we have transformed the second index and 
      ! A[ij](t,u) = g(ij|tu)

      ! **********************
      ! *** 2ND HALF TRANSFORM
      ! **********************

      ! loop over symmetries for t


      do t_sym = 1 , nirrep_

        ! determine symmetry for u

        u_sym = group_mult_tab_(t_sym,ij_sym)

        ! make sure that we only address integrals with t_sym > u_sym

        if ( t_sym < u_sym ) cycle

        ! figure out number of orbitals with these symmetries

        num_t    = trans_%nmopi(t_sym)
        num_u    = trans_%nmopi(u_sym)

        ! this information is used in the scatter operation below

        t_offset = trans_%offset(t_sym)
        u_offset = trans_%offset(u_sym)
  
        ! loop over t indeces

        do t_irrep = 1 , num_t

          ! this information is used in the scatter operation below
  
          t_class = trans_%irrep_to_class_map(t_offset+t_irrep) 

          ! loop over u indeces

          do u_irrep = 1 , num_u

            ! this information is used in the scatter operation below

            u_class = trans_%irrep_to_class_map(u_offset+u_irrep)

            tu_class = ints_%gemind(t_class,u_class)   

            ! **********
            ! *** GATHER
            ! **********

            ! loop over symmetries for i

            do i_sym = 1 , nirrep_
 
              ! corresponding symmetry for j

              j_sym = group_mult_tab_(i_sym,ij_sym)

              ! make sure that we are only addressing unique integrals with i_sym > j_sym

              if ( i_sym < j_sym ) cycle

              ! number of i/j orbitals

              num_i = trans_%nmopi(i_sym)
              num_j = trans_%nmopi(j_sym)              

              ! offsets for indexing

              i_offset = trans_%offset(i_sym)
              j_offset = trans_%offset(j_sym)

              ! zero out temporary matrix

              B_tilde = 0.0_wp

              ! loop over i indeces

              do i_irrep = 1 , num_i

                ! i index in class order

                i_class = trans_%irrep_to_class_map(i_irrep+i_offset)

                do j_irrep = 1 , num_j
  
                  ! j index in class order

                  j_class      = trans_%irrep_to_class_map(j_irrep+j_offset)                
                  
                  ij_class     = ints_%gemind(i_class,j_class)

                  ! save the corresponding matrix element

                  B_tilde(i_irrep,j_irrep) = A(ij_class)%irrep(t_sym)%mat(t_irrep,u_irrep)

                end do ! end j_irrep loop

              end do ! end i_irrep loop

              ! *************
              ! *** TRANSFORM
              ! *************
              
              ! A_tilde[tu][i_sym] = B_tilde[tu][i_sym] * C[[j_sym]
              ! followed by
              ! B_tilde[tu][i_sym] = (C[i_sym])^T * A_tilde[tu][i_sym]

              ! only need to do work if there are orbitals with symmetry i_sym and j_sym

              if ( ( num_i == 0 ) .or. ( num_j == 0 ) ) cycle 
 
# ifdef BLAS
              call dgemm('N','N',num_i,num_j,num_j,1.0_wp,B_tilde,max_nmopi,&
                     trans_%u_irrep_block(j_sym)%val,num_j,0.0_wp,A_tilde,max_nmopi)
              call dgemm('T','N',num_i,num_j,num_i,1.0_wp,trans_%u_irrep_block(i_sym)%val,num_i, &
                     A_tilde,max_nmopi,0.0_wp,B_tilde,max_nmopi)
# else
              A_tilde(1:num_i,1:num_j) = matmul(B_tilde(1:num_i,1:num_j),trans_%u_irrep_block(j_sym)%val)
              B_tilde(1:num_i,1:num_j) = matmul(transpose(trans_%u_irrep_block(i_sym)%val),A_tilde(1:num_i,1:num_j))

# endif
              ! ***********
              ! *** SCATTER
              ! ***********

              ! at this point ij-->rs && kl-->tu so B_tilde[tu][i_sym]%mat(i,j) = g(rs|tu)

              ! loop over indeces for i

              do i_irrep = 1 , num_i

                ! i index in class order

                i_class = trans_%irrep_to_class_map(i_irrep+i_offset)

                do j_irrep = 1 , num_j

                  ! j index in class order

                  j_class               = trans_%irrep_to_class_map(j_irrep+j_offset)

                  ! ij geminal index in class order

                  ij_class              = ints_%gemind(i_class,j_class)

                  ! integral address in class order

                  rstu_class            = pq_index(ij_class,tu_class)

                  ! save integral

                  int_block(rstu_class) = B_tilde(i_irrep,j_irrep)

                end do ! end j_irrep loop

              end do ! end i_irrep loop

            end do ! end i_sym loop

          end do ! end u_irrep loop

        end do ! end t_irrep loop

      end do ! end t_sym loop

      ! ******** END ACTUAL WORK

      ! deallocate temporary arrays

      transform_teints_irrep_block = deallocate_transform_scr()

      return

      contains

        integer function allocate_transform_scr()

          integer :: ij,k_sym,l_sym,num_k,num_l
          ! function to allocate scratch arrays for integral transformation

          allocate_transform_scr = 0 

          if (allocated(A)) allocate_transform_scr = deallocate_transform_scr()

          if ( allocate_transform_scr /= 0 ) stop

          allocate(A(nnz_ij))

          do ij = 1 , nnz_ij

            allocate(A(ij)%irrep(nirrep_))

            do k_sym = 1 , nirrep_

              l_sym = group_mult_tab_(k_sym,ij_sym)

              if ( k_sym < l_sym ) cycle

              num_k = trans_%nmopi(k_sym)
              num_l = trans_%nmopi(l_sym)

              ! do not allocate if there are no orbitals with these symmetries

              if ( ( num_k == 0 ) .or. ( num_l == 0 ) ) cycle         

              allocate(A(ij)%irrep(k_sym)%mat(num_k,num_l))

              A(ij)%irrep(k_sym)%mat = 0.0_wp

            end do

          end do
          
          if (allocated(A_tilde)) deallocate(A_tilde)

          allocate(A_tilde(max_nmopi,max_nmopi)) 

          if (allocated(B_tilde)) deallocate(B_tilde)

          allocate(B_tilde(max_nmopi,max_nmopi))

          allocate_transform_scr = 0

          return
        end function allocate_transform_scr

        integer function deallocate_transform_scr()

          ! function to deallocate scratch arrays for integral transformation
          integer :: ij,k_sym

          deallocate_transform_scr = 1

          if (.not.allocated(A)) return

          do ij = 1 , nnz_ij
            
            if (.not.allocated(A(ij)%irrep)) cycle

            do k_sym = 1 , nirrep_

              if (allocated(A(ij)%irrep(k_sym)%mat)) deallocate(A(ij)%irrep(k_sym)%mat)

            end do

            deallocate(A(ij)%irrep)

          end do

          deallocate(A)

          if (allocated(A_tilde)) deallocate(A_tilde)

          if (allocated(B_tilde)) deallocate(B_tilde)

          deallocate_transform_scr = 0

          return

        end function deallocate_transform_scr

    end function transform_teints_irrep_block

    integer function transform_teints_g0_block(int_block,ij_sym)
      implicit none
      ! subroutine to transform an irrep-block (totally symmetric) of two-electron integrals
      ! the algorithm performs the transformation in two half transfroms.
      ! 1) The first half transform consists of two steps:
      !    a) GATHER:    for a given ij pair, collect and store all integrals (ij|kl)
      !    b) TRANSFORM: for a given ij pair, perform two dense matrix-matrix multiplications to 
      !                   transform (ij|kl) --> (ij|tu)
      !
      ! To do the second half transform, half-transformed integrals need to be stored for ALL ij and tu
      ! This means, additional storage that is roughly twice the size of the LT integral block 
      !
      ! 2) The second half transform consists of three steps
      !    a) GATHER:    for a given tu pair, collect and store all integrals (ij|tu)
      !    c) TRANSFORM: for a given ij pair, perform two dense matrix-matrix multiplications to 
      !                  transform (ij|tu) --> (rs|tu)
      !    d) SCATTER:   scatter (rs|tu) into original integral array
      ! 
      ! NOTE: Currently, the GATHER and SCATTER operations do more work than the minimum required; however,
      !       after timing the total time spent in each type of operation, it is clear that the GATHER and
      !       SCATTER operations are an almost negligible fraction of the total time and thus, further
      !       optimization will be skipped. 

      integer, intent(in) :: ij_sym
      real(wp) :: int_block(:)

      type int_scr
        real(wp), allocatable :: mat(:,:)
      end type int_scr

      type int_tmp
        type(int_scr), allocatable :: irrep(:)
      end type int_tmp

      type(int_tmp), allocatable :: A(:)

      real(wp), allocatable :: A_tilde(:,:),B_tilde(:,:)

      integer :: nnz_ij,ij,max_nmopi,error

      integer :: i_class,j_class,k_class,l_class,t_class,u_class,tu_class,kl_class,ij_class
      integer :: i_irrep,j_irrep,k_irrep,l_irrep,t_irrep,u_irrep
      integer :: i_sym,j_sym,k_sym,l_sym,t_sym,u_sym 
      integer :: num_i,num_j,num_k,num_l,num_t,num_u
      integer :: i_offset,j_offset,k_offset,l_offset,t_offset,u_offset

      integer(ip) :: ijkl_class,rstu_class

      integer :: first_ij(nthread_use_),last_ij_(nthread_use_)

      ! initialize error flag

      transform_teints_g0_block = 1
 
      ! the number of nnz integrals is not the same as the length of the vector int2

      if ( size(int_block,dim=1) /= ints_%nnzpi(ij_sym) ) return

      ! maximum matrix size and number of matrices

      max_nmopi = maxval(trans_%nmopi)

      nnz_ij    = ints_%ngempi(ij_sym)

      ! allocate temporary arrays

      transform_teints_g0_block = allocate_transform_scr_g0()

      if ( transform_teints_g0_block /= 0 ) stop

      ! ********************************
      ! *** 1ST HALF TRANSFORM (kl-->tu)
      ! ********************************

      ! loop over ij geminals

      do ij_class = 1 , nnz_ij

        ! loop over symmetries for k

        do k_sym = 1 , nirrep_

          ! figure out symmetry of l such that ij_sym == kl_sym

          ! *******************************************
          ! *** GATHER *** A[ij][k_sym](k,l) = g(ij|kl)
          ! *******************************************

          ! here, we make use of the fact that this matrix has at most 8 nnz blocks,
          ! with the blocks indexed by k_sym

          ! offsets for figuring out full k & l indeces

          k_offset = trans_%offset(k_sym)

          ! number of orbitals with symmetry k_sym 

          num_k    = trans_%nmopi(k_sym)

          ! loop over symmetry-reduced k indeces

          do k_irrep = 1 , num_k

            ! class k index for integral addressing

            k_class = trans_%irrep_to_class_map(k_irrep+k_offset)

            ! loop over symmetry reduced l indeces

            do l_irrep = 1 , k_irrep

              ! class l index for integral addressing

              l_class    = trans_%irrep_to_class_map(l_irrep+k_offset)

              ! integral address (class order)

              kl_class   = ints_%gemind(k_class,l_class)
              ijkl_class = pq_index(ij_class,kl_class) 

              ! save integral

              A(ij_class)%irrep(k_sym)%mat(k_irrep,l_irrep) = int_block(ijkl_class) 
              A(ij_class)%irrep(k_sym)%mat(l_irrep,k_irrep) = A(ij_class)%irrep(k_sym)%mat(k_irrep,l_irrep)

            end do ! end l loop

          end do ! end k loop

          !*********************************************
          ! *** TRANSFORM A[ij] = C^T * A[ij][k_sym] * C
          !*********************************************

          ! A_tilde[ij][k_sym] = A[ij][k_sym] * C[k_sym]
          ! followed by
          ! A = (C[k_sym])^T * A_tilde[]ij][k_sym]

          ! only need to do work if there are orbitals with symmetry k_sym  

          if ( num_k == 0 )  cycle

# ifdef BLAS
          call dgemm('N','N',num_k,num_k,num_k,1.0_wp,A(ij_class)%irrep(k_sym)%mat,num_k,&
                     trans_%u_irrep_block(k_sym)%val,num_k,0.0_wp,A_tilde,max_nmopi)
          call dgemm('T','N',num_k,num_k,num_k,1.0_wp,trans_%u_irrep_block(k_sym)%val,num_k, &
                     A_tilde,max_nmopi,0.0_wp,A(ij_class)%irrep(k_sym)%mat,num_k)
# else
          A_tilde(1:num_k,1:num_k)       = matmul(A(ij_class)%irrep(k_sym)%mat,trans_%u_irrep_block(k_sym)%val)
          A(ij_class)%irrep(k_sym)%mat   = matmul(transpose(trans_%u_irrep_block(k_sym)%val),A_tilde(1:num_k,1:num_k))
# endif

        end do ! end k_sym loop

      end do ! end ij_class loop      

      ! at this point, we have transformed the second index and 
      ! A[ij](t,u) = g(ij|tu)

      ! **********************
      ! *** 2ND HALF TRANSFORM
      ! **********************

      ! loop over symmetries for t

      do t_sym = 1 , nirrep_

        ! figure out number of orbitals with these symmetries

        num_t    = trans_%nmopi(t_sym)

        ! this information is used in the scatter operation below

        t_offset = trans_%offset(t_sym)
  
        ! loop over t indeces

        do t_irrep = 1 , num_t

          ! this information is used in the scatter operation below
  
          t_class = trans_%irrep_to_class_map(t_offset+t_irrep) 

          ! loop over u indeces

          do u_irrep = 1 , t_irrep

            ! this information is used in the scatter operation below

            u_class = trans_%irrep_to_class_map(t_offset+u_irrep)

            tu_class = ints_%gemind(t_class,u_class)   

            ! **********
            ! *** GATHER
            ! **********

            ! loop over symmetries for i

            do i_sym = 1 , nirrep_
 
              num_i = trans_%nmopi(i_sym)

              ! offsets for indexing

              i_offset = trans_%offset(i_sym)

              ! zero out temporary matrix

              B_tilde = 0.0_wp

              ! loop over i indeces

              do i_irrep = 1 , num_i

                ! i index in class order

                i_class = trans_%irrep_to_class_map(i_irrep+i_offset)

                do j_irrep = 1 , i_irrep
  
                  ! j index in class order

                  j_class      = trans_%irrep_to_class_map(j_irrep+i_offset)                
                  
                  ij_class     = ints_%gemind(i_class,j_class)

                  ! save the corresponding matrix element

                  B_tilde(i_irrep,j_irrep) = A(ij_class)%irrep(t_sym)%mat(t_irrep,u_irrep)
                  B_tilde(j_irrep,i_irrep) = B_tilde(i_irrep,j_irrep)

                end do ! end j_irrep loop

              end do ! end i_irrep loop

              ! *************
              ! *** TRANSFORM
              ! *************

              ! A_tilde[tu][i_sym] = B_tilde[tu][i_sym] * C[[i_sym]
              ! followed by
              ! B_tilde[tu][i_sym] = (C[i_sym])^T * A_tilde[tu][i_sym]

              ! only need to do work if there are orbitals with symmetry k_sym  
 
              if ( num_i == 0 )  cycle

# ifdef BLAS
              call dgemm('N','N',num_i,num_i,num_i,1.0_wp,B_tilde,max_nmopi,&
                     trans_%u_irrep_block(i_sym)%val,num_i,0.0_wp,A_tilde,max_nmopi)
              call dgemm('T','N',num_i,num_i,num_i,1.0_wp,trans_%u_irrep_block(i_sym)%val,num_i, &
                     A_tilde,max_nmopi,0.0_wp,B_tilde,max_nmopi)
# else
              A_tilde(1:num_i,1:num_i) = matmul(B_tilde(1:num_i,1:num_i),trans_%u_irrep_block(i_sym)%val)
              B_tilde(1:num_i,1:num_i) = matmul(transpose(trans_%u_irrep_block(i_sym)%val),A_tilde(1:num_i,1:num_i))
# endif

              ! ***********
              ! *** SCATTER
              ! ***********

              ! at this point ij-->rs && kl-->tu so B_tilde[tu][i_sym]%mat(i,j) = g(rs|tu)

              ! loop over indeces for i

              do i_irrep = 1 , num_i

                ! i index in class order

                i_class = trans_%irrep_to_class_map(i_irrep+i_offset)

                do j_irrep = 1 , i_irrep

                  ! j index in class order

                  j_class               = trans_%irrep_to_class_map(j_irrep+i_offset)

                  ! ij geminal index in class order

                  ij_class              = ints_%gemind(i_class,j_class)

                  ! integral address in class order

                  rstu_class            = pq_index(ij_class,tu_class)

                  ! save integral

                  int_block(rstu_class) = B_tilde(i_irrep,j_irrep)


                end do ! end j_irrep loop

              end do ! end i_irrep loop

            end do ! end i_sym loop

          end do ! end u_irrep loop

        end do ! end t_irrep loop

      end do ! end t_sym loop

      ! ******** END ACTUAL WORK

      ! deallocate temporary arrays

      transform_teints_g0_block = deallocate_transform_scr_g0()

      return

      contains

        integer function allocate_transform_scr_g0()

          integer :: ij,k_sym,l_sym,num_k,num_l
          ! function to allocate scratch arrays for integral transformation

          allocate_transform_scr_g0 = 0 

          if (allocated(A)) allocate_transform_scr_g0 = deallocate_transform_scr_g0()

          if ( allocate_transform_scr_g0 /= 0 ) stop

          allocate(A(nnz_ij))

          do ij = 1 , nnz_ij

            allocate(A(ij)%irrep(nirrep_))

            do k_sym = 1 , nirrep_

              num_k = trans_%nmopi(k_sym)
          
              ! do not allocate if there are no orbitals with this symemtry
            
              if ( num_k == 0 ) cycle
 
              allocate(A(ij)%irrep(k_sym)%mat(num_k,num_k))

              A(ij)%irrep(k_sym)%mat = 0.0_wp

            end do

          end do
          
          if (allocated(A_tilde)) deallocate(A_tilde)

          allocate(A_tilde(max_nmopi,max_nmopi)) 

          if (allocated(B_tilde)) deallocate(B_tilde)

          allocate(B_tilde(max_nmopi,max_nmopi))

          allocate_transform_scr_g0 = 0

          return
        end function allocate_transform_scr_g0

        integer function deallocate_transform_scr_g0()

          ! function to deallocate scratch arrays for integral transformation
          integer :: ij,k_sym

          deallocate_transform_scr_g0 = 1

          if (.not.allocated(A)) return

          do ij = 1 , nnz_ij
            
            if (.not.allocated(A(ij)%irrep)) cycle

            do k_sym = 1 , nirrep_

              if (allocated(A(ij)%irrep(k_sym)%mat)) deallocate(A(ij)%irrep(k_sym)%mat)

            end do

            deallocate(A(ij)%irrep)

          end do

          deallocate(A)

          if (allocated(A_tilde)) deallocate(A_tilde)

          if (allocated(B_tilde)) deallocate(B_tilde)

          deallocate_transform_scr_g0 = 0

          return

        end function deallocate_transform_scr_g0

    end function transform_teints_g0_block

end module focas_transform_teints
