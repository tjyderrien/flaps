! test program for sparse matrix inversion
! TJY Derrien

program SparseInvert

implicit none

real(8) Amatrix(1:4,1:4), Bvector(1:4), Xvector(1:4)
integer i,j 
Amatrix(:,:)=0d0; Bvector(:)=0d0; Xvector(:)=0d0

do i=1, 4
  do j=1, 4
    Amatrix(i,i)=2d0; Amatrix(i,i-1)=1d0; Amatrix(i,i+1)=-1d0
  end do
end do
    
write(*,*) transpose(Amatrix)
    
Xvector(:)=1d0

Bvector(:)=2d0

call SolveSparse((Amatrix),Xvector,Bvector,1,1)

write(*,*) Bvector


  contains 

    subroutine SolveSparse(A, X, B, KL, KU)
      real(8), dimension(:,:), intent(in) :: A
      real(8), dimension(size(A,1),1):: X, B
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info, ldb, KL, KU, i, j
      real(8), dimension(2*KL+KU+1, 2*KL+KU+1) :: AB
      integer dimAx, dimAy
      dimAx=size(A,1); dimAy=size(A,2)
      n=dimAx*dimAy
      info=0
      ldb=0
      i=0; j=0
      
      ! External procedures defined in LAPACK
      ! external DGBSV !LU inversion routine for sparse matrix

      do i=1,dimAx
	do j=1,dimAy
	  if(max(1,j-KU) <= i .AND. min(N,j+KL) >= i) then
	    AB(KL+KU+1+i-j,j)=A(i,j)
	  end if
	end do
      end do
      
      n = size(A,1)

      ! DGBSV computes the solution to a real system of linear equations
      ! A * X = B, where A is a band matrix of order N with KL subdiagonals
      ! and KU superdiagonals, and X and B are N-by-NRHS matrices.
       CALL DGBSV( n, KL, KU, 1, AB, 2*KL+KU+1, ipiv, B, n, INFO )

       write(*,*) Info
       
      if (info /= 0) then
	stop 'Matrix is numerically singular!'
      end if

      if (info /= 0) then
	stop 'Matrix inversion failed!'
      end if
    end subroutine SolveSparse


end program SparseInvert