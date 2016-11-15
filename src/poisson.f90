!! Copyright (C) 2012-2016 T. J.-Y. Derrien, N. Tancogne-Dejean
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>

!------------------------------------------------------------------------------
!> @file poisson.f90
!
! DESCRIPTION:
!> @brief This is for the future Poisson calculation
!------------------------------------------------------------------------------


subroutine poisson_init_dual( M,N, x, y, xDualSW, yDualSW, xDualSE, yDualSE, &
                                         xDualNE, yDualNE, xDualNW, yDualNW, &
                                         xDual, yDual)
  use Maths_m
  implicit none

  integer,                 intent(in)    :: M, N
  real(8), dimension(M,N), intent(in)    :: x, y
  real(8), dimension(M,N), intent(inout) :: xDualSW, yDualSW, xDualSE, yDualSE, &
                                            xDualNE, yDualNE, xDualNW, yDualNW, &
                                            xDual, yDual

  integer :: i, j

  xDualSW(:,:) = M_ZERO; yDualSW(:,:) = M_ZERO
  xDualSE(:,:) = M_ZERO; yDualSE(:,:) = M_ZERO
  xDualNE(:,:) = M_ZERO; yDualNE(:,:) = M_ZERO
  xDualNW(:,:) = M_ZERO; yDualNW(:,:) = M_ZERO
  xDual(:,:) = M_ZERO; yDual(:,:) = M_ZERO

  !TODO: OpemMP parallelisation here?

  !TODO: Variables should be Np and Mp, xP and yP here no?
  !      TJYD@NTD: Yes! 
  ! interpolation and preparation of resolution
  do j=2,N-1
    do i=2,M-1
        !Dual mesh calculation
        !TODO: This should be in mesh.f90, no? Why working on the dual mesh for Poisson? This is only for Crossed derivativatives (cf transport).
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




subroutine poisson_init_normal_cellarea( Mp, Np, xP, yP, NormalNxP, NormalNyP, NormalSxP, NormalSyP, NormalExP, NormalEyP, &
                                            NormalWxP, NormalWyP, CellAreaNP, CellAreaSP, CellAreaEP, CellAreaWP )
  use Maths_m
  implicit none

  integer,                 intent(in)      :: Mp, Np
  real(8), dimension(Mp,Np), intent(in)    :: xP, yP
  real(8), dimension(Mp,Np), intent(inout) :: NormalNxP, NormalNyP, NormalSxP, NormalSyP, NormalExP, NormalEyP, &
                                            NormalWxP, NormalWyP, CellAreaNP, CellAreaSP, CellAreaEP, CellAreaWP

  integer :: i, j


   do j=2,Np-1
    do i=2,Mp-1
        NormalNxP(i,j)=4d0*(0.25d0*yP(i-1,j)+0.25d0*yP(i-1,j+1)-0.25d0*yP(i+1,j+1)-0.25d0*yP(i+1,j))/(xP(i-1,j)**2 &
                      +M_TWO*xP(i-1,j)*xP(i-1,j+1)-M_TWO*xP(i-1,j)*xP(i+1,j+1)-M_TWO*xP(i-1,j)*xP(i+1,j)+xP(i-1,j+1)**2 &
                      -M_TWO*xP(i-1,j+1)*xP(i+1,j+1)-M_TWO*xP(i-1,j+1)*xP(i+1,j)+xP(i+1,j+1)**2+M_TWO*xP(i+1,j+1)*xP(i+1,j) &
                      +xP(i+1,j)**2+yP(i-1,j)**2+M_TWO*yP(i-1,j)*yP(i-1,j+1)-M_TWO*yP(i-1,j)*yP(i+1,j+1)-M_TWO*yP(i-1,j)*yP(i+1,j) &
                      +yP(i-1,j+1)**2-M_TWO*yP(i-1,j+1)*yP(i+1,j+1)-M_TWO*yP(i-1,j+1)*yP(i+1,j)+yP(i+1,j+1)**2 &
                      +M_TWO*yP(i+1,j+1)*yP(i+1,j)+yP(i+1,j)**2)**(M_HALF)
        NormalNyP(i,j)=-4d0*(0.25d0*xP(i-1,j)+0.25d0*xP(i-1,j+1)-0.25d0*xP(i+1,j+1)-0.25d0*xP(i+1,j))/(xP(i-1,j)**2 &
                +M_TWO*xP(i-1,j)*xP(i-1,j+1)-M_TWO*xP(i-1,j)*xP(i+1,j+1)-M_TWO*xP(i-1,j)*xP(i+1,j)+xP(i-1,j+1)**2 &
                -M_TWO*xP(i-1,j+1)*xP(i+1,j+1)-M_TWO*xP(i-1,j+1)*xP(i+1,j)+xP(i+1,j+1)**2+M_TWO*xP(i+1,j+1)*xP(i+1,j) &
                +xP(i+1,j)**2+yP(i-1,j)**2+M_TWO*yP(i-1,j)*yP(i-1,j+1)-M_TWO*yP(i-1,j)*yP(i+1,j+1)-M_TWO*yP(i-1,j)*yP(i+1,j) &
                +yP(i-1,j+1)**2-M_TWO*yP(i-1,j+1)*yP(i+1,j+1)-M_TWO*yP(i-1,j+1)*yP(i+1,j)+yP(i+1,j+1)**2 &
                +M_TWO*yP(i+1,j+1)*yP(i+1,j)+yP(i+1,j)**2)**(M_HALF)
        NormalSxP(i,j)=-4d0*(0.25d0*yP(i-1,j-1)+0.25d0*yP(i-1,j)-0.25d0*yP(i+1,j)-0.25d0*yP(i+1,j-1))/(xP(i-1,j-1)**2 &
                +M_TWO*xP(i-1,j-1)*xP(i-1,j)-M_TWO*xP(i-1,j-1)*xP(i+1,j)-M_TWO*xP(i-1,j-1)*xP(i+1,j-1)+xP(i-1,j)**2 &
                -M_TWO*xP(i-1,j)*xP(i+1,j)-M_TWO*xP(i-1,j)*xP(i+1,j-1)+xP(i+1,j)**2+M_TWO*xP(i+1,j)*xP(i+1,j-1) &
                +xP(i+1,j-1)**2+yP(i-1,j-1)**2+M_TWO*yP(i-1,j-1)*yP(i-1,j)-M_TWO*yP(i-1,j-1)*yP(i+1,j) &
                -M_TWO*yP(i-1,j-1)*yP(i+1,j-1)+yP(i-1,j)**2-M_TWO*yP(i-1,j)*yP(i+1,j)-M_TWO*yP(i-1,j)*yP(i+1,j-1) &
                +yP(i+1,j)**2+M_TWO*yP(i+1,j)*yP(i+1,j-1)+yP(i+1,j-1)**2)**(M_HALF)
        NormalSyP(i,j)=4d0*(0.25d0*xP(i-1,j-1)+0.25d0*xP(i-1,j)-0.25d0*xP(i+1,j)-0.25d0*xP(i+1,j-1))/(xP(i-1,j-1)**2 &
                      +M_TWO*xP(i-1,j-1)*xP(i-1,j)-M_TWO*xP(i-1,j-1)*xP(i+1,j)-M_TWO*xP(i-1,j-1)*xP(i+1,j-1)+xP(i-1,j)**2 &
                      -M_TWO*xP(i-1,j)*xP(i+1,j)-M_TWO*xP(i-1,j)*xP(i+1,j-1)+xP(i+1,j)**2+M_TWO*xP(i+1,j)*xP(i+1,j-1) &
                      +xP(i+1,j-1)**2+yP(i-1,j-1)**2+M_TWO*yP(i-1,j-1)*yP(i-1,j)-M_TWO*yP(i-1,j-1)*yP(i+1,j) &
                      -M_TWO*yP(i-1,j-1)*yP(i+1,j-1)+yP(i-1,j)**2-M_TWO*yP(i-1,j)*yP(i+1,j)-M_TWO*yP(i-1,j)*yP(i+1,j-1) &
                      +yP(i+1,j)**2+M_TWO*yP(i+1,j)*yP(i+1,j-1)+yP(i+1,j-1)**2)**(M_HALF)
        NormalExP(i,j)=-4d0*(0.25d0*yP(i+1,j-1)+0.25d0*yP(i,j-1)-0.25d0*yP(i+1,j+1)-0.25d0*yP(i,j+1))/(xP(i+1,j-1)**2 &
                +M_TWO*xP(i+1,j-1)*xP(i,j-1)-M_TWO*xP(i+1,j-1)*xP(i+1,j+1)-M_TWO*xP(i+1,j-1)*xP(i,j+1)+xP(i,j-1)**2 &
                -M_TWO*xP(i,j-1)*xP(i+1,j+1)-M_TWO*xP(i,j-1)*xP(i,j+1)+xP(i+1,j+1)**2+M_TWO*xP(i+1,j+1)*xP(i,j+1) &
                +xP(i,j+1)**2+yP(i+1,j-1)**2+M_TWO*yP(i+1,j-1)*yP(i,j-1)-M_TWO*yP(i+1,j-1)*yP(i+1,j+1) &
                -M_TWO*yP(i+1,j-1)*yP(i,j+1)+yP(i,j-1)**2-M_TWO*yP(i,j-1)*yP(i+1,j+1)-M_TWO*yP(i,j-1)*yP(i,j+1) &
                +yP(i+1,j+1)**2+M_TWO*yP(i+1,j+1)*yP(i,j+1)+yP(i,j+1)**2)**(M_HALF)
        NormalEyP(i,j)=4d0*(0.25d0*xP(i+1,j-1)+0.25d0*xP(i,j-1)-0.25d0*xP(i+1,j+1)-0.25d0*xP(i,j+1))/(xP(i+1,j-1)**2 &
                      +M_TWO*xP(i+1,j-1)*xP(i,j-1)-M_TWO*xP(i+1,j-1)*xP(i+1,j+1)-M_TWO*xP(i+1,j-1)*xP(i,j+1)+xP(i,j-1)**2 &
                      -M_TWO*xP(i,j-1)*xP(i+1,j+1)-M_TWO*xP(i,j-1)*xP(i,j+1)+xP(i+1,j+1)**2+M_TWO*xP(i+1,j+1)*xP(i,j+1) &
                      +xP(i,j+1)**2+yP(i+1,j-1)**2+M_TWO*yP(i+1,j-1)*yP(i,j-1)-M_TWO*yP(i+1,j-1)*yP(i+1,j+1) &
                      -M_TWO*yP(i+1,j-1)*yP(i,j+1)+yP(i,j-1)**2-M_TWO*yP(i,j-1)*yP(i+1,j+1)-M_TWO*yP(i,j-1)*yP(i,j+1) &
                      +yP(i+1,j+1)**2+M_TWO*yP(i+1,j+1)*yP(i,j+1)+yP(i,j+1)**2)**(M_HALF)
        NormalWxP(i,j)=4d0*(0.25d0*yP(i-1,j-1)+0.25d0*yP(i,j-1)-0.25d0*yP(i,j+1)-0.25d0*yP(i-1,j+1))/(xP(i-1,j-1)**2 &
                      +M_TWO*xP(i-1,j-1)*xP(i,j-1)-M_TWO*xP(i-1,j-1)*xP(i,j+1)-M_TWO*xP(i-1,j-1)*xP(i-1,j+1)+xP(i,j-1)**2 &
                      -M_TWO*xP(i,j-1)*xP(i,j+1)-M_TWO*xP(i,j-1)*xP(i-1,j+1)+xP(i,j+1)**2+M_TWO*xP(i,j+1)*xP(i-1,j+1) &
                      +xP(i-1,j+1)**2+yP(i-1,j-1)**2+M_TWO*yP(i-1,j-1)*yP(i,j-1)-M_TWO*yP(i-1,j-1)*yP(i,j+1) &
                      -M_TWO*yP(i-1,j-1)*yP(i-1,j+1)+yP(i,j-1)**2-M_TWO*yP(i,j-1)*yP(i,j+1)-M_TWO*yP(i,j-1)*yP(i-1,j+1) &
                      +yP(i,j+1)**2+M_TWO*yP(i,j+1)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(M_HALF)
        NormalWyP(i,j)=-4d0*(0.25d0*xP(i-1,j-1)+0.25d0*xP(i,j-1)-0.25d0*xP(i,j+1)-0.25d0*xP(i-1,j+1))/(xP(i-1,j-1)**2 &
                +M_TWO*xP(i-1,j-1)*xP(i,j-1)-M_TWO*xP(i-1,j-1)*xP(i,j+1)-M_TWO*xP(i-1,j-1)*xP(i-1,j+1)+xP(i,j-1)**2 &
                -M_TWO*xP(i,j-1)*xP(i,j+1)-M_TWO*xP(i,j-1)*xP(i-1,j+1)+xP(i,j+1)**2+M_TWO*xP(i,j+1)*xP(i-1,j+1) &
                +xP(i-1,j+1)**2+yP(i-1,j-1)**2+M_TWO*yP(i-1,j-1)*yP(i,j-1)-M_TWO*yP(i-1,j-1)*yP(i,j+1) &
                -M_TWO*yP(i-1,j-1)*yP(i-1,j+1)+yP(i,j-1)**2-M_TWO*yP(i,j-1)*yP(i,j+1)-M_TWO*yP(i,j-1)*yP(i-1,j+1) &
                +yP(i,j+1)**2+M_TWO*yP(i,j+1)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(M_HALF)

        CellAreaNP(i,j)=0.25d0*(xP(i+1,j+1)**2+M_TWO*xP(i+1,j+1)*xP(i+1,j)-M_TWO*xP(i+1,j+1)*xP(i-1,j) &
                -M_TWO*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-M_TWO*xP(i+1,j)*xP(i-1,j)-M_TWO*xP(i+1,j)*xP(i-1,j+1) &
                +xP(i-1,j)**2+M_TWO*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+M_TWO*yP(i+1,j+1)*yP(i+1,j) &
                -M_TWO*yP(i+1,j+1)*yP(i-1,j)-M_TWO*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-M_TWO*yP(i+1,j)*yP(i-1,j) &
                -M_TWO*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+M_TWO*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(M_HALF)
        CellAreaSP(i,j)=0.25d0*(xP(i+1,j)**2+M_TWO*xP(i+1,j)*xP(i+1,j-1)-M_TWO*xP(i+1,j)*xP(i-1,j-1)-M_TWO*xP(i+1,j)*xP(i-1,j) &
                +xP(i+1,j-1)**2-M_TWO*xP(i+1,j-1)*xP(i-1,j-1)-M_TWO*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2 &
                +M_TWO*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+M_TWO*yP(i+1,j)*yP(i+1,j-1) &
                -M_TWO*yP(i+1,j)*yP(i-1,j-1)-M_TWO*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-M_TWO*yP(i+1,j-1)*yP(i-1,j-1) &
                -M_TWO*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+M_TWO*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(M_HALF)
        CellAreaEP(i,j)=0.25d0*(xP(i+1,j+1)**2+M_TWO*xP(i+1,j+1)*xP(i,j+1) &
                       -M_TWO*xP(i+1,j+1)*xP(i+1,j-1)-M_TWO*xP(i+1,j+1)*xP(i,j-1) &
                        +xP(i,j+1)**2-M_TWO*xP(i,j+1)*xP(i+1,j-1) &
                       -M_TWO*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+M_TWO*xP(i+1,j-1)*xP(i,j-1) &
                +xP(i,j-1)**2+yP(i+1,j+1)**2+M_TWO*yP(i+1,j+1)*yP(i,j+1)-M_TWO*yP(i+1,j+1)*yP(i+1,j-1)-M_TWO*yP(i+1,j+1)*yP(i,j-1) &
                +yP(i,j+1)**2-M_TWO*yP(i,j+1)*yP(i+1,j-1)-M_TWO*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+M_TWO*yP(i+1,j-1)*yP(i,j-1) &
                +yP(i,j-1)**2)**(M_HALF)
        CellAreaWP(i,j)=0.25d0*(xP(i,j+1)**2+M_TWO*xP(i,j+1)*xP(i-1,j+1) &
                       -M_TWO*xP(i,j+1)*xP(i-1,j-1)-M_TWO*xP(i,j+1)*xP(i,j-1) &
                        +xP(i-1,j+1)**2-M_TWO*xP(i-1,j+1)*xP(i-1,j-1) &
                        -M_TWO*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+M_TWO*xP(i-1,j-1)*xP(i,j-1) &
                        +xP(i,j-1)**2+yP(i,j+1)**2+M_TWO*yP(i,j+1)*yP(i-1,j+1) &
                        -M_TWO*yP(i,j+1)*yP(i-1,j-1)-M_TWO*yP(i,j+1)*yP(i,j-1) &
                        +yP(i-1,j+1)**2-M_TWO*yP(i-1,j+1)*yP(i-1,j-1)-M_TWO*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2 &
                        +M_TWO*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(M_HALF)
    end do
  end do

!      do i=2,M-1
!
!         !TODO: I think that this should not be here
!         ! TJYD: I don't think it is actually used. Can be removed!  
!         GradNeX(i,N) = GradNeX(i,N-1)
!         GradNeY(i,N) = GradNeY(i,N-1)
!      end do


   write(*,*) "VESSEL CHECK"
   write(*,*) "North", NormalNxP(Mp/2,Np-1), NormalNyP(Mp/2,Np-1)
   write(*,*) "South", NormalSxP(Mp/2,2), NormalSyP(Mp/2,2)
   write(*,*) "East", NormalExP(Mp-1,Np-1), NormalEyP(Mp-1,Np-1)
   write(*,*) "West", NormalWxP(Mp-1,Np-1), NormalWyP(Mp-1,Np-1)
   write(*,*) CellAreaNP(Mp/2,Np/2), CellAreaSP(Mp/2,Np/2), CellAreaEP(Mp/2,Np/2), CellAreaWP(Mp/2,Np/2)

end subroutine poisson_init_normal_cellarea

