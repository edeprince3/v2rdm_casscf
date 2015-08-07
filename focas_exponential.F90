module focas_exponential

  use focas_data

  implicit none

  real(wp), parameter :: max_error_tolerance = 1.0e-10_wp

  contains

    subroutine compute_exponential(kappa_in)
      implicit none
      ! subroutine to compute matrix exponential of a skew-symmetric real matrix K
      ! U = exp(K)
      ! Since the matrix K is block-diagonal, U is as well. Thus, the matrix
      ! exponential is computed for each block separately
      
      integer :: i_sym,max_nmopi,error
      real(wp), intent(in) :: kappa_in(:)
      real(wp), allocatable :: k_block(:,:)

      ! maximum size of temporary matrix

      max_nmopi = maxval(trans_%nmopi)
 
      ! allocate temporary matrix

      allocate(k_block(max_nmopi,max_nmopi))

      do i_sym=1,nirrep_

        error = gather_kappa_block(kappa_in,k_block,i_sym)
        if (error /= 0) stop

        error = compute_block_exponential(k_block,i_sym,max_nmopi)
        if (error /= 0) stop
         
      end do

      deallocate(k_block)

      return
 
    end subroutine compute_exponential

    integer function compute_block_exponential(K,block_sym,max_dim)
      implicit none
      ! function to compute the matrix exponential of a matrix according to
      ! U = exp(K) = X * cos(d) * X^T + K * X * d^(-1) * sin(d) * X^T
      ! where X and d are solutions of the eigenvalue equation K^2 * X = lambda * X
      ! and d = sqrt(-lambda)
      integer, intent(in)  :: block_sym,max_dim
      real(wp), intent(in) :: K(max_dim,max_dim)

      real(wp), allocatable :: K2(:,:),X(:,:),tmp_mat_1(:,:),tmp_mat_2(:,:)
      real(wp), allocatable :: d(:),work(:)
      integer, allocatable :: isuppz(:),iwork(:)
      integer :: iwork_tmp(1)
      real(wp) :: work_tmp(1,1),val

      integer :: block_dim,i,j, il,iu,neig_found,lwork,liwork,success
      real(wp) :: vl,vu,diag_tol,max_normalization_error,max_orthogonality_error,ddot

      ! initialize return value
      compute_block_exponential = 1

      ! dimension of actual block
      
      block_dim = trans_%nmopi(block_sym)

      ! ***************************
      ! allocate temporary matrices
      ! ***************************

      allocate(K2(block_dim,block_dim))
      allocate(d(block_dim))
      allocate(X(block_dim,block_dim))

      ! ***************
      ! compute and K^2
      ! ***************

#ifdef BLAS

      do i = 1 , block_dim
        call dcopy(block_dim,K(1:block_dim,i),1,X(:,i),1)
      end do
      call dgemm('n','n',block_dim,block_dim,block_dim,1.0_wp,K,max_dim,X,block_dim,0.0_wp,K2,block_dim)

#else

      K2 = matmul(K(1:block_dim,1:block_dim),K(1:block_dim,1:block_dim))

#endif         

      ! ***************
      ! diagonalize K^2
      ! ***************

      il=1
      iu=block_dim
      vl=-huge(1.0_wp)
      vu=0.0_wp
      diag_tol=2.0_wp*epsilon(1.0_wp)
      allocate(isuppz(2*block_dim))
      ! figure out optimal dimensions for iwork and work
      call dsyevr('v','a','u',block_dim,K2,block_dim,vl,vu,il,iu,diag_tol,neig_found,d,X,block_dim,isuppz,work_tmp,-1,iwork_tmp,-1,success)
      if ( success /= 0 ) return
      liwork=iwork_tmp(1)
      lwork=int(work_tmp(1,1))
      allocate(work(lwork),iwork(liwork))
      call dsyevr('v','a','u',block_dim,K2,block_dim,vl,vu,il,iu,diag_tol,neig_found,d,X,block_dim,isuppz,work,lwork,iwork,liwork,success)
      deallocate(isuppz,work,iwork)
      if ( success /= 0 ) return

      ! ***********************************************************************************
      ! compute exponential U = exp(K) = X * cos(d) * X^(T) + K * X * d^(-1) * sin(d) * X^T
      ! first compute the terms involving sin(d) followed by the terms involving cos(d)
      ! ***********************************************************************************

      ! scale eigenvalues
      do i = 1 , block_dim
        if (d(i) < 0.0_wp ) then 
          d(i) = sqrt(-d(i))
        else
          d(i) = sqrt(d(i))
        end if
      end do

      allocate(tmp_mat_1(block_dim,block_dim))
      allocate(tmp_mat_2(block_dim,block_dim))

      ! compute d^(-1) * sin(d) * X^T
      ! since both d and sin(d) are diagonal matrices, d^(-1) * sin(d) is also diagonal with diagonal elements sin(d(i))/d(i)
      ! since d^(-1) * sin(d) is diagonal, we can perform the matrix product C = D * M efficiently by recognizing that the ith
      ! row of C is just a scaled version of the ith row of M ... C(:,i) = D(i,i) * M(:,i) 
      ! note to self :: below, we are storing the transpose of d^(-1) * sin(d) * X^T since we directly copy rows --> rows

      do i = 1 , block_dim
        val = 1.0_wp
        if ( d(i) /= 0.0_wp ) val = sin(d(i)) / d(i) 

#ifdef BLAS

        call dcopy(block_dim,X(:,i),1,tmp_mat_1(:,i),1)
        call dscal(block_dim,val,tmp_mat_1(:,i),1)

#else

        tmp_mat_1(:,i) = val * X(:,i)

#endif      

      end do

      ! now do the rest in two steps to compute the full second term
      ! tmp_mat_2 = X * tmp_mat_1^T
      ! U = K * tmp_mat_2

#ifdef BLAS

      call dgemm('n','t',block_dim,block_dim,block_dim,1.0_wp,X,block_dim,tmp_mat_1, &
                & block_dim,0.0_wp,tmp_mat_2,block_dim)
      call dgemm('n','n',block_dim,block_dim,block_dim,1.0_wp,K,max_dim,tmp_mat_2, &
                & block_dim,0.0_wp,trans_%u_irrep_block(block_sym)%val,block_dim)

#else

      tmp_mat_2                           = matmul(X,transpose(tmp_mat_1))
      trans_%u_irrep_block(block_sym)%val = matmul(K(1:block_dim,1:block_dim),tmp_mat_2)

#endif        
     
      ! ****************************
      ! *** COMPUTE X * cos(d) * X^T
      ! ****************************

      ! tmp_mat_1 =  cos(d) * X

      do i = 1 , block_dim
        val = cos(d(i))

#ifdef BLAS

        call dcopy(block_dim,X(:,i),1,tmp_mat_1(:,i),1)
        call dscal(block_dim,val,tmp_mat_1(:,i),1)        

#else

       tmp_mat_1(:,i) = val * X(:,i)

#endif  
      end do 

      ! tmp_mat_2 = X * (tmp_mat_1)^T
      ! U = U + tmp_mat_2

#ifdef BLAS

      call dgemm('n','t',block_dim,block_dim,block_dim,1.0_wp,X,block_dim,tmp_mat_1, & 
                & block_dim,1.0_wp,trans_%u_irrep_block(block_sym)%val,block_dim)

#else

      tmp_mat_2                           = matmul(X,transpose(tmp_mat_1))
      trans_%u_irrep_block(block_sym)%val = trans_%u_irrep_block(block_sym)%val + tmp_mat_2

#endif

      ! check orthonormality error

      max_orthogonality_error = 0.0_wp
      max_normalization_error = 0.0_wp

      do i = 1 , block_dim

        ! compute norm of vector

#ifdef BLAS

        val = ddot(block_dim,trans_%u_irrep_block(block_sym)%val(:,i),1,&
              trans_%u_irrep_block(block_sym)%val(:,i),1)

#else

        val = dot_product(trans_%u_irrep_block(block_sym)%val(:,i),trans_%u_irrep_block(block_sym)%val(:,i))

#endif
        if ( abs( 1.0_wp - val ) > max_normalization_error ) max_normalization_error = abs ( 1.0_wp - val )

        do j = 1 , i - 1

#ifdef BLAS

          val = ddot(block_dim,trans_%u_irrep_block(block_sym)%val(:,i),1,&
              trans_%u_irrep_block(block_sym)%val(:,j),1)

#else

          val = dot_product(trans_%u_irrep_block(block_sym)%val(:,i),trans_%u_irrep_block(block_sym)%val(:,j))

#endif         

          if ( abs( val ) > max_orthogonality_error ) max_orthogonality_error = abs ( val )
           
        end do      
 
      end do

      if ( ( max_normalization_error > max_error_tolerance ) .or. &
         & ( max_orthogonality_error > max_error_tolerance ) ) then
        write(*,'(a,1x,i1,5x,a,1x,i3,5x,a,1x,es10.3,5x,a,1x,es10.3)')'irrep:',block_sym,'nmo:',block_dim,&
             & 'max(normalization_error):',max_normalization_error,'max(orthogonality_error):',max_orthogonality_error
      endif 

      ! ******************************
      ! deallocate tempporary matrices
      ! ******************************
     
      deallocate(tmp_mat_1,tmp_mat_2)
      deallocate(K2,d,X)

      compute_block_exponential = 0

      return
    end function compute_block_exponential

    integer function gather_kappa_block(kappa_in,block,block_sym)
      implicit none
      real(wp) :: kappa_in(:)
      real(wp) :: block(:,:)
      integer, intent(in) :: block_sym
      integer :: i,j,n_ij,ij,i_irrep,j_irrep,i_class,j_class,j_class_start,j_start

      gather_kappa_block = 1

      ! number of orbital pairs in this block
      n_ij = trans_%npairpi(block_sym)

      if ( n_ij == 0 ) then
        gather_kappa_block = 0
        return
      end if

      ! figure out first/last orbital pair index for this block 

      ij = 0
      if ( block_sym > 1 ) ij = sum(trans_%npairpi(1:block_sym-1))
      
      ! initialize matrix block
      block = 0.0_wp

      ! loop over i/j pairs in this block

      ! the loop structure below cycles through the possible rotation pairs
      ! rotation pairs are sorted according to symmetry and for each irrep,
      ! the rotation pairs are sorted according to orbital classes: ad,ed,aa,ea
      ! for each pair j>i; for more details, see subroutine setup_rotation_indeces in focas_main.F90

 
      do i_class = 1 , 3

        j_class_start = i_class + 1

        if ( ( include_aa_rot_ == 1 ) .and. ( i_class == 2 ) ) j_class_start = i_class

        do j_class = j_class_start , 3

          do i = first_index_(block_sym,i_class) , last_index_(block_sym,i_class)

            j_start = first_index_(block_sym,j_class)

            if ( i_class == j_class ) j_start = i + 1

            do j = j_start , last_index_(block_sym,j_class)

              ! figure out symmetry_reduced indeces

              i_irrep                = trans_%class_to_irrep_map(i)
              j_irrep                = trans_%class_to_irrep_map(j)

              ! update gradient index

              ij                     = ij + 1

              ! copy into matrix

              block(j_irrep,i_irrep) =  kappa_in(ij)
              block(i_irrep,j_irrep) = -kappa_in(ij)

            end do

          end do

        end do    

      end do ! end ij loop

      gather_kappa_block = 0

    end function gather_kappa_block

end module focas_exponential
