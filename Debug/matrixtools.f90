 ! files containing the routines for matrix manipulation


    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function FillMatrixCondition(i,j,kl,ku, Nsolve)
    ! returns boolean if position i,j has to be filled in the band diagonal matrix for Lapack solver DGBTRF/DGBTRS.
      implicit none
      integer(8):: i,j, 		&
		 ku,kl, 	&
		 Nsolve, 	& !size of the matrix A
		 FillMatrixCondition

      FillMatrixCondition=((max(1,j-ku) <= i) .AND. (min(Nsolve,j+kl) >= i))
      return
    end function FillMatrixCondition

    function FillMatrixRow(i,j,kl,ku, Nsolve)
    ! returns row index to fill the band diagonal matrix for Lapack solver DGBTRF/DGBTRS.
      implicit none
      integer(8):: i,j,		&!Input matrix indexes
		 ku,kl,		&!Number of upper and lower diagonals
		 Nsolve, 	&!number of rows of the full band diagonal matrix
		 FillMatrixRow	 !Output for ABsolve solve matrix
      logical FillMatrixCondition
      
      if(FillMatrixCondition(i,j,kl,ku, Nsolve)) then
	FillMatrixRow=kl+ku+1+i-j
      else
	FillMatrixRow=-1 !induce a crash
      end if
      return
    end function FillMatrixRow

    ! Returns the inverse of a sparse matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    subroutine SolveSparse(A, X, B, KL, KU)
      real(8), dimension(:,:), intent(in) :: A
      real(8), dimension(size(A,1),1):: X, B
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info, ldb, KL, KU, i, j
      real(8), dimension(2*KL+KU+1, 2*KL+KU+1) :: AB
      integer dimAx, dimAy
      dimAx=size(A,1); dimAy=size(A,2)

      ! External procedures defined in LAPACK
      ! external DGBSV !LU inversion routine for sparse matrix

      ! Store A in Ainv to prevent it from being overwritten by LAPACK

      do i=1,dimAx
	do j=1,dimAy
	  if(max(1,j-KU) <= i .AND. min(M,j+KL) >= i) then
	    AB(KU+1+i-j,j)=A(i,j)
	  end if
	end do
      end do
!                   if(mod(nbiter,iterOut).eq.0) then
! 		  write(91,'(5E30.15)') AB !BsolveP !Amatrix
! 		  write(91,*)
! 	    end if
      !
      n = size(A,1)

      ! DGBSV computes the solution to a real system of linear equations
      ! A * X = B, where A is a band matrix of order N with KL subdiagonals
      ! and KU superdiagonals, and X and B are N-by-NRHS matrices.
       CALL DGBSV( n, KL, KU, 1, AB, 2*KL+KU+1, ipiv, B, n, INFO )

      if (info /= 0) then
	stop 'Matrix is numerically singular!'
      end if

      if (info /= 0) then
	stop 'Matrix inversion failed!'
      end if
    end subroutine SolveSparse
