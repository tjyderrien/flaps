!------------------------------------------------------------------------------
!> @file mie.f90
!
! DESCRIPTION:
!> @brief This is for the future Poisson calculation
!
!> @author
!> Thibault J.Y. Derrien
!
!> @date
!> 15 Jun 2016 - Initial Version
!------------------------------------------------------------------------------


subroutine poisson_init_dual( M,N, x, y, xDualSW, yDualSW, xDualSE, yDualSE, &
                                         xDualNE, yDualNE, xDualNW, yDualNW, &
                                         xDual, yDual)
  implicit none
  integer,                 intent(in)    :: M, N
  real(8), dimension(M,N), intent(in)    :: x, y
  real(8), dimension(M,N), intent(inout) :: xDualSW, yDualSW, xDualSE, yDualSE, &
                                            xDualNE, yDualNE, xDualNW, yDualNW, &
                                            xDual, yDual

  integer :: i, j

  xDualSW(:,:) = 0.d0; yDualSW(:,:) = 0.d0
  xDualSE(:,:) = 0.d0; yDualSE(:,:) = 0.d0
  xDualNE(:,:) = 0.d0; yDualNE(:,:) = 0.d0
  xDualNW(:,:) = 0.d0; yDualNW(:,:) = 0.d0
  xDual(:,:) = 0.d0; yDual(:,:) = 0.d0

  !TODO: OpemMP parallelisation here?

  ! interpolation and preparation of resolution
  do j=2,N-1
    do i=2,M-1
        !Dual mesh calculation
        xDualSW(i,j) = 0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)) !x(i-1/2,j-1/2)
        yDualSW(i,j) = 0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j))
        xDualSE(i,j) = 0.25d0*(x(i,j-1)+x(i+1,j-1)+x(i+1,j)+x(i,j)) !x(i+1/2,j-1/2)
        yDualSE(i,j) = 0.25d0*(y(i,j-1)+y(i+1,j-1)+y(i+1,j)+y(i,j))
        xDualNE(i,j) = 0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)) !x(i+1/2,j+1/2)
        yDualNE(i,j) = 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))
        xDualNW(i,j) = 0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)) !x(i-1/2,j+1/2)
        yDualNW(i,j) = 0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1))

        xDual(i-1,j-1) = xDualSW(i,j); yDual(i-1,j-1)=yDualSW(i,j) ! to make for (i,j) from 2,2 to M,N
        xDual(i-1,N-1) = xDualNW(i,N-1); yDual(i-1,N-1) = yDualNW(i,N-1);
        xDual(M-1,j-1) = xDualSE(M-1,j); yDual(M-1,j-1) = yDualSE(M-1,j)
        xDual(M-1,N-1) = xDualNE(M-1,N-1); yDual(M-1,N-1) = yDualNE(M-1,N-1);
    end do
  end do

end subroutine poisson_init_dual
