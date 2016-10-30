!! Copyright (C) 2012-2016 T. J.-Y. Derrien
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
!> @file grid.f90
!
! DESCRIPTION:
!> @brief Everything related to the grid
!
!> @author
!> Thibault J.Y. Derrien
!
!> @date
!> 07 Jun 2016 - Initial Version
!------------------------------------------------------------------------------

subroutine allocate_NormCurviTangent(M, N, NormalN, NormalS, NormalE, &
                                 NormalW, TangentNx, TangentNy, TangentSx, TangentSy, &
                                 TangentEx, TangentEy, TangentWx, TangentWy, CurviNx, CurviNy, &
                                 CurviSx, CurviSy, CurviEx, CurviEy, CurviWx, CurviWy )
  use Types_m
  implicit none

  integer, intent(in)    :: M, N
  type(VectorField), intent(inout) :: NormalN, NormalS, NormalW, NormalE ! normal to quadrangle elements

  real(8), dimension(:,:), allocatable,intent(inout) ::  TangentWx, TangentWy, &                ! Tangent to quadrangle elements
                             TangentEx, TangentEy, &
                             TangentNx, TangentNy, &
                             TangentSx, TangentSy, &
                             CurviWx, CurviWy, &                 ! Unit vector between cell centers
                             CurviEx, CurviEy, &
                             CurviNx, CurviNy, &
                             CurviSx, CurviSy

  allocate(NormalN%x(1:M,1:N))
  allocate(NormalN%y(1:M,1:N))
  allocate(NormalN%N(1:M,1:N))
  allocate(NormalS%x(1:M,1:N))
  allocate(NormalS%y(1:M,1:N))
  allocate(NormalS%N(1:M,1:N))
  allocate(NormalE%x(1:M,1:N))
  allocate(NormalE%y(1:M,1:N))
  allocate(NormalE%N(1:M,1:N))
  allocate(NormalW%x(1:M,1:N))
  allocate(NormalW%y(1:M,1:N))
  allocate(NormalW%N(1:M,1:N))

end subroutine allocate_NormCurviTangent

subroutine compute_distances(M, N, x, y, DistN, DistS, DistE, DistW, DistDualN, DistDualS, &
                                         DistDualE, DistDualW, CellAreaN, CellAreaS, CellAreaE, CellAreaW )
  use Maths_m
  implicit none

  integer, intent(in) :: M, N

  real(8), intent(in) :: x(1:M, 1:N), y(1:M, 1:N)                 ! needle position indexes
  real(8), intent(inout) :: CellAreaN(1:M, 1:N), CellAreaS(1:M, 1:N), &                        ! area of the finite elements
                            CellAreaE(1:M, 1:N), CellAreaW(1:M, 1:N), &
                            DistW(1:M,1:N), DistE(1:M,1:N), &                 ! distance to the center of neighboor cells
                            DistN(1:M,1:N), DistS(1:M,1:N), &
                            DistDualW(1:M,1:N), DistDualE(1:M,1:N), &                 ! distance of the element side (equal to area in 2D)
                            DistDualN(1:M,1:N), DistDualS(1:M,1:N)

  integer :: i, j

   do j=2,N-1
      do i=2,M-1

        CellAreaN(i,j)=sqrt(( 0.25d0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))-0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)))**2 &
                  +(  0.25d0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))-0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)))**2) !0.25d0*(x(i+1,j+1)**2+M_TWO*x(i+1,j+1)*x(i+1,j)-M_TWO*x(i+1,j+1)*x(i-1,j)-M_TWO*x(i+1,j+1)*x(i-1,j+1)+x(i+1,j)**2-M_TWO*x(i+1,j)*x(i-1,j)-M_TWO*x(i+1,j)*x(i-1,j+1)+x(i-1,j)**2+M_TWO*x(i-1,j)*x(i-1,j+1)+x(i-1,j+1)**2+y(i+1,j+1)**2+M_TWO*y(i+1,j+1)*y(i+1,j)-M_TWO*y(i+1,j+1)*y(i-1,j)-M_TWO*y(i+1,j+1)*y(i-1,j+1)+y(i+1,j)**2-M_TWO*y(i+1,j)*y(i-1,j)-M_TWO*y(i+1,j)*y(i-1,j+1)+y(i-1,j)**2+M_TWO*y(i-1,j)*y(i-1,j+1)+y(i-1,j+1)**2)**(M_HALF)
        CellAreaS(i,j)=sqrt( (0.25d0*(x(i+1,j-1)+x(i,j-1)+x(i+1,j)+x(i,j))-0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i-1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i+1,j-1)+y(i,j-1)+y(i+1,j)+y(i,j))-0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i-1,j)+y(i,j)))**2 ) !0.25d0*(x(i+1,j)**2+M_TWO*x(i+1,j)*x(i+1,j-1)-M_TWO*x(i+1,j)*x(i-1,j-1)-M_TWO*x(i+1,j)*x(i-1,j)+x(i+1,j-1)**2-M_TWO*x(i+1,j-1)*x(i-1,j-1)-M_TWO*x(i+1,j-1)*x(i-1,j)+x(i-1,j-1)**2+M_TWO*x(i-1,j-1)*x(i-1,j)+x(i-1,j)**2+y(i+1,j)**2+M_TWO*y(i+1,j)*y(i+1,j-1)-M_TWO*y(i+1,j)*y(i-1,j-1)-M_TWO*y(i+1,j)*y(i-1,j)+y(i+1,j-1)**2-M_TWO*y(i+1,j-1)*y(i-1,j-1)-M_TWO*y(i+1,j-1)*y(i-1,j)+y(i-1,j-1)**2+M_TWO*y(i-1,j-1)*y(i-1,j)+y(i-1,j)**2)**(M_HALF)
        CellAreaE(i,j)=sqrt( (0.25d0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))-0.25d0*(x(i+1,j-1)+x(i,j-1)+x(i+1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))-0.25d0*(y(i+1,j-1)+y(i,j-1)+y(i+1,j)+y(i,j)))**2 ) !0.25d0*(x(i+1,j+1)**2+M_TWO*x(i+1,j+1)*x(i,j+1)-M_TWO*x(i+1,j+1)*x(i+1,j-1)-M_TWO*x(i+1,j+1)*x(i,j-1)+x(i,j+1)**2-M_TWO*x(i,j+1)*x(i+1,j-1)-M_TWO*x(i,j+1)*x(i,j-1)+x(i+1,j-1)**2+M_TWO*x(i+1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i+1,j+1)**2+M_TWO*y(i+1,j+1)*y(i,j+1)-M_TWO*y(i+1,j+1)*y(i+1,j-1)-M_TWO*y(i+1,j+1)*y(i,j-1)+y(i,j+1)**2-M_TWO*y(i,j+1)*y(i+1,j-1)-M_TWO*y(i,j+1)*y(i,j-1)+y(i+1,j-1)**2+M_TWO*y(i+1,j-1)*y(i,j-1)+y(i,j-1)**2)**(M_HALF)
        CellAreaW(i,j)=sqrt( (0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j))-0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i-1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j))-0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i-1,j)+y(i,j)))**2 ) !0.25d0*(x(i,j+1)**2+M_TWO*x(i,j+1)*x(i-1,j+1)-M_TWO*x(i,j+1)*x(i-1,j-1)-M_TWO*x(i,j+1)*x(i,j-1)+x(i-1,j+1)**2-M_TWO*x(i-1,j+1)*x(i-1,j-1)-M_TWO*x(i-1,j+1)*x(i,j-1)+x(i-1,j-1)**2+M_TWO*x(i-1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i,j+1)**2+M_TWO*y(i,j+1)*y(i-1,j+1)-M_TWO*y(i,j+1)*y(i-1,j-1)-M_TWO*y(i,j+1)*y(i,j-1)+y(i-1,j+1)**2-M_TWO*y(i-1,j+1)*y(i-1,j-1)-M_TWO*y(i-1,j+1)*y(i,j-1)+y(i-1,j-1)**2+M_TWO*y(i-1,j-1)*y(i,j-1)+y(i,j-1)**2)**(M_HALF)

        DistN(i,j)=sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
        DistS(i,j)=sqrt((x(i,j)-x(i,j-1))**2+(y(i,j)-y(i,j-1))**2)
        DistE(i,j)=sqrt((x(i+1,j)-x(i,j))**2+(y(i+1,j)-y(i,j))**2)
        DistW(i,j)=sqrt((x(i,j)-x(i-1,j))**2+(y(i,j)-y(i-1,j))**2)

      end do
   end do

   do i=2,M-1
      !NORTH
         CellAreaN(i,N)=M_ZERO !M_HALF*(x(i+1,N)**2-M_TWO*x(i+1,N)*x(i-1,N)+x(i-1,N)**2+y(i+1,N)**2-M_TWO*y(i+1,N)*y(i-1,N)+y(i-1,N)**2)**M_HALF
         CellAreaS(i,N)=0.25d0*(x(i+1,N)**2+M_TWO*x(i+1,N)*x(i+1,N-1)-M_TWO*x(i+1,N)*x(i-1,N-1)-M_TWO*x(i+1,N)*x(i-1,N) &
                +x(i+1,N-1)**2-M_TWO*x(i+1,N-1)*x(i-1,N-1)-M_TWO*x(i+1,N-1)*x(i-1,N)+x(i-1,N-1)**2+M_TWO*x(i-1,N)*x(i-1,N-1) &
                +x(i-1,N)**2+y(i+1,N)**2+M_TWO*y(i+1,N)*y(i+1,N-1)-M_TWO*y(i+1,N)*y(i-1,N-1)-M_TWO*y(i+1,N)*y(i-1,N) &
                +y(i+1,N-1)**2-M_TWO*y(i+1,N-1)*y(i-1,N-1)-M_TWO*y(i+1,N-1)*y(i-1,N)+y(i-1,N-1)**2 &
                                                        +M_TWO*y(i-1,N)*y(i-1,N-1)+y(i-1,N)**2)**M_HALF
         CellAreaW(i,N)=0.25d0*(x(i,N)**2+M_TWO*x(i,N)*x(i-1,N)-M_TWO*x(i,N)*x(i-1,N-1)-M_TWO*x(i,N)*x(i,N-1)+x(i-1,N)**2 &
                -M_TWO*x(i-1,N)*x(i-1,N-1)-M_TWO*x(i-1,N)*x(i,N-1)+x(i-1,N-1)**2+M_TWO*x(i-1,N-1)*x(i,N-1)+x(i,N-1)**2 &
                +y(i,N)**2+M_TWO*y(i,N)*y(i-1,N)-M_TWO*y(i,N)*y(i-1,N-1)-M_TWO*y(i,N)*y(i,N-1)+y(i-1,N)**2 &
                -M_TWO*y(i-1,N)*y(i-1,N-1)-M_TWO*y(i-1,N)*y(i,N-1)+y(i-1,N-1)**2+M_TWO*y(i-1,N-1)*y(i,N-1) &
                +y(i,N-1)**2)**M_HALF
         CellAreaE(i,N)=0.25d0*(x(i+1,N)**2+M_TWO*x(i+1,N)*x(i,N)-M_TWO*x(i+1,N)*x(i+1,N-1)-M_TWO*x(i+1,N)*x(i,N-1)+x(i,N)**2 &
                -M_TWO*x(i,N)*x(i+1,N-1)-M_TWO*x(i,N)*x(i,N-1)+x(i+1,N-1)**2+M_TWO*x(i+1,N-1)*x(i,N-1)+x(i,N-1)**2 &
                +y(i,N)**2+M_TWO*y(i,N)*y(i+1,N)-M_TWO*y(i,N)*y(i+1,N-1)-M_TWO*y(i,N)*y(i,N-1)+y(i+1,N)**2 &
                -M_TWO*y(i+1,N)*y(i+1,N-1)-M_TWO*y(i+1,N)*y(i,N-1)+y(i+1,N-1)**2+M_TWO*y(i+1,N-1)*y(i,N-1) &
                +y(i,N-1)**2)**M_HALF

         DistN(i,N)=0d0 !sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
         DistS(i,N)=sqrt((x(i,N)-x(i,N-1))**2+(y(i,N)-y(i,N-1))**2)
         DistE(i,N)=sqrt((x(i+1,N)-x(i,N))**2+(y(i+1,N)-y(i,N))**2)
         DistW(i,N)=sqrt((x(i,N)-x(i-1,N))**2+(y(i,N)-y(i-1,N))**2)

     !SOUTH
         CellAreaN(i,1)=0.25d0*(x(i+1,2)**2+M_TWO*x(i+1,2)*x(i+1,1)-M_TWO*x(i+1,2)*x(i-1,1)-M_TWO*x(i+1,2)*x(i-1,2) &
                 +x(i+1,1)**2-M_TWO*x(i+1,1)*x(i-1,1)-M_TWO*x(i+1,1)*x(i-1,2)+x(i-1,1)**2+M_TWO*x(i-1,1)*x(i-1,2)+x(i-1,2)**2 &
                 +y(i+1,2)**2+M_TWO*y(i+1,2)*y(i+1,1)-M_TWO*y(i+1,2)*y(i-1,1)-M_TWO*y(i+1,2)*y(i-1,2)+y(i+1,1)**2 &
                 -M_TWO*y(i+1,1)*y(i-1,1)-M_TWO*y(i+1,1)*y(i-1,2)+y(i-1,1)**2+M_TWO*y(i-1,1)*y(i-1,2)+y(i-1,2)**2)**M_HALF
         CellAreaS(i,1)=0d0 !sqrt((M_HALF*(x(i+1,1)+x(i,1))-M_HALF*(x(i-1,1)+x(i,1)))**2+(M_HALF*(y(i+1,1)+y(i,1))-M_HALF*(y(i-1,1)+y(i,1)))**2)
         CellAreaW(i,1)=0.25d0*(x(i,1)**2-M_TWO*x(i,1)*x(i,2)+M_TWO*x(i,1)*x(i-1,1)-M_TWO*x(i,1)*x(i-1,2)+x(i,2)**2 &
                  -M_TWO*x(i,2)*x(i-1,1)+M_TWO*x(i,2)*x(i-1,2)+x(i-1,1)**2&
                 -M_TWO*x(i-1,1)*x(i-1,2)+x(i-1,2)**2+y(i,1)**2-M_TWO*y(i,1)*y(i,2) &
                  +M_TWO*y(i,1)*y(i-1,1)-M_TWO*y(i,1)*y(i-1,2)+y(i,2)**2 &
                  -M_TWO*y(i,2)*y(i-1,1)+M_TWO*y(i,2)*y(i-1,2)+y(i-1,1)**2 &
                  -M_TWO*y(i-1,1)*y(i-1,2)+y(i-1,2)**2)**M_HALF
         CellAreaE(i,1)=0.25d0*(x(i+1,2)**2-M_TWO*x(i+1,2)*x(i+1,1)-M_TWO*x(i+1,2)*x(i,1)&
                                +M_TWO*x(i+1,2)*x(i,2)+x(i+1,1)**2+M_TWO*x(i+1,1)*x(i,1) &
                                -M_TWO*x(i+1,1)*x(i,2)+x(i,1)**2-M_TWO*x(i,1)*x(i,2)+x(i,2)**2 &
                                +y(i+1,2)**2-M_TWO*y(i+1,2)*y(i+1,1)-M_TWO*y(i+1,2)*y(i,1) &
                                +M_TWO*y(i+1,2)*y(i,2)+y(i+1,1)**2+M_TWO*y(i+1,1)*y(i,1)&
                                -M_TWO*y(i+1,1)*y(i,2)+y(i,1)**2-M_TWO*y(i,1)*y(i,2)+y(i,2)**2)**M_HALF

         DistN(i,1)=sqrt( (x(i,2)-x(i,1))**2 + (y(i,2)-y(i,1))**2 )
         DistS(i,1)=M_ZERO
         DistE(i,1)=sqrt( (x(i+1,1)-x(i,1))**2 + (y(i+1,1)-y(i,1))**2 )
         DistW(i,1)=sqrt( (x(i,1)-x(i-1,1))**2 + (y(i,1)-y(i-1,1))**2 )
   end do

   do j=2,N-1
     !WEST
           CellAreaN(1,j)=0.25d0*(x(2,j+1)**2+M_TWO*x(2,j+1)*x(2,j)-M_TWO*x(2,j+1)*x(1,j)-M_TWO*x(2,j+1)*x(1,j+1)+x(2,j)**2 &
                  -M_TWO*x(2,j)*x(1,j)-M_TWO*x(2,j)*x(1,j+1)+x(1,j)**2+M_TWO*x(1,j)*x(1,j+1)+x(1,j+1)**2+y(2,j+1)**2 &
                  +M_TWO*y(2,j+1)*y(2,j)-M_TWO*y(2,j+1)*y(1,j)-M_TWO*y(2,j+1)*y(1,j+1)+y(2,j)**2-M_TWO*y(2,j)*y(1,j) &
                  -M_TWO*y(2,j)*y(1,j+1)+y(1,j)**2+M_TWO*y(1,j)*y(1,j+1)+y(1,j+1)**2)**M_HALF
          CellAreaS(1,j)=0.25d0*(x(2,j)**2+M_TWO*x(2,j)*x(2,j-1)-M_TWO*x(2,j)*x(1,j)-M_TWO*x(2,j)*x(1,j-1)+x(2,j-1)**2 &
                  -M_TWO*x(2,j-1)*x(1,j)-M_TWO*x(2,j-1)*x(1,j-1)+x(1,j)**2+M_TWO*x(1,j)*x(1,j-1)+x(1,j-1)**2+y(2,j)**2 &
                  +M_TWO*y(2,j)*y(2,j-1)-M_TWO*y(2,j)*y(1,j)-M_TWO*y(2,j)*y(1,j-1)+y(2,j-1)**2-M_TWO*y(2,j-1)*y(1,j) &
                  -M_TWO*y(2,j-1)*y(1,j-1)+y(1,j)**2+M_TWO*y(1,j)*y(1,j-1)+y(1,j-1)**2)**M_HALF
          CellAreaW(1,j)=sqrt( (M_HALF*(x(1,j+1)+x(1,j))-M_HALF*(x(1,j-1)+x(1,j)))**2 + (M_HALF*(y(1,j+1)+y(1,j)) &
                        -M_HALF*(y(1,j-1)+y(1,j)))**2 )
          CellAreaE(1,j)=0.25d0*(x(2,j+1)**2+M_TWO*x(2,j+1)*x(1,j+1)-M_TWO*x(2,j+1)*x(2,j-1)-M_TWO*x(2,j+1)*x(1,j-1) &
                  +x(1,j+1)**2-M_TWO*x(1,j+1)*x(2,j-1)-M_TWO*x(1,j+1)*x(1,j-1)+x(2,j-1)**2+M_TWO*x(2,j-1)*x(1,j-1) &
                  +x(1,j-1)**2+y(2,j+1)**2+M_TWO*y(2,j+1)*y(1,j+1)-M_TWO*y(2,j+1)*y(2,j-1)-M_TWO*y(2,j+1)*y(1,j-1) &
                  +y(1,j+1)**2-M_TWO*y(1,j+1)*y(2,j-1)-M_TWO*y(1,j+1)*y(1,j-1)+y(2,j-1)**2+M_TWO*y(2,j-1)*y(1,j-1) &
                  +y(1,j-1)**2)**M_HALF

          DistN(1,j)=sqrt( (x(1,j+1)-x(1,j))**2 + (y(1,j+1)-y(1,j))**2 )
          DistS(1,j)=sqrt( (x(1,j)-x(1,j-1))**2 + (y(1,j)-y(1,j-1))**2 )
          DistE(1,j)=sqrt( (x(2,j)-x(1,j))**2 + (y(2,j)-y(1,j))**2 )
          DistW(1,j)=M_ZERO

     !EAST
          CellAreaN(M,j)=0.25d0*(x(M,j+1)**2+M_TWO*x(M,j+1)*x(M,j)-M_TWO*x(M,j+1)*x(M-1,j)-M_TWO*x(M,j+1)*x(M-1,j+1)+x(M,j)**2 &
                  -M_TWO*x(M,j)*x(M-1,j)-M_TWO*x(M,j)*x(M-1,j+1)+x(M-1,j)**2+M_TWO*x(M-1,j)*x(M-1,j+1)+x(M-1,j+1)**2 &
                  +y(M,j+1)**2+M_TWO*y(M,j+1)*y(M,j)-M_TWO*y(M,j+1)*y(M-1,j)-M_TWO*y(M,j+1)*y(M-1,j+1)+y(M,j)**2 &
                  -M_TWO*y(M,j)*y(M-1,j)-M_TWO*y(M,j)*y(M-1,j+1)+y(M-1,j)**2+M_TWO*y(M-1,j)*y(M-1,j+1)+y(M-1,j+1)**2)**M_HALF
          CellAreaS(M,j)=0.25d0*(x(M,j)**2+M_TWO*x(M,j)*x(M,j-1)-M_TWO*x(M,j)*x(M-1,j-1)-M_TWO*x(M,j)*x(M-1,j)+x(M,j-1)**2 &
                  -M_TWO*x(M-1,j-1)*x(M,j-1)-M_TWO*x(M,j-1)*x(M-1,j)+x(M-1,j-1)**2+M_TWO*x(M-1,j-1)*x(M-1,j)+x(M-1,j)**2 &
                  +y(M,j)**2+M_TWO*y(M,j)*y(M,j-1)-M_TWO*y(M,j)*y(M-1,j-1)-M_TWO*y(M,j)*y(M-1,j)+y(M,j-1)**2 &
                  -M_TWO*y(M-1,j-1)*y(M,j-1)-M_TWO*y(M,j-1)*y(M-1,j)+y(M-1,j-1)**2+M_TWO*y(M-1,j-1)*y(M-1,j) &
                  +y(M-1,j)**2)**M_HALF
          CellAreaW(M,j)=0.25d0*(x(M,j+1)**2+M_TWO*x(M,j+1)*x(M-1,j+1)-M_TWO*x(M,j+1)*x(M-1,j-1)-M_TWO*x(M,j+1)*x(M,j-1) &
                  +x(M-1,j+1)**2-M_TWO*x(M-1,j+1)*x(M-1,j-1)-M_TWO*x(M-1,j+1)*x(M,j-1)+x(M-1,j-1)**2+M_TWO*x(M-1,j-1)*x(M,j-1) &
                  +x(M,j-1)**2+y(M,j+1)**2+M_TWO*y(M,j+1)*y(M-1,j+1)-M_TWO*y(M,j+1)*y(M-1,j-1)-M_TWO*y(M,j+1)*y(M,j-1) &
                  +y(M-1,j+1)**2-M_TWO*y(M-1,j+1)*y(M-1,j-1)-M_TWO*y(M-1,j+1)*y(M,j-1)+y(M-1,j-1)**2 &
                  +M_TWO*y(M-1,j-1)*y(M,j-1)+y(M,j-1)**2)**M_HALF
          CellAreaE(M,j)=sqrt( (M_HALF*(x(M,j+1)+x(M,j))-M_HALF*(x(M,j-1)+x(M,j)) )**2 + (M_HALF*(y(M,j+1)+y(M,j)) &
                        -M_HALF*(y(M,j-1)+y(M,j)) )**2)

          DistN(M,j)=sqrt( (x(M,j+1)-x(M,j))**2 + (y(M,j+1)-y(M,j))**2 )
          DistS(M,j)=sqrt( (x(M,j-1)-x(M,j))**2 + (y(M,j-1)-y(M,j))**2 )
          DistE(M,j)=M_ZERO
          DistW(M,j)=sqrt( (x(M-1,j)-x(M,j))**2 + (y(M-1,j)-y(M,j))**2 )

   end do

   !North-east
        CellAreaE(M,N)=sqrt((x(M,N)-M_HALF*(x(M,N)+x(M,N-1)))**2+(y(M,N)-M_HALF*(y(M,N)+y(M,N-1)))**2)
        CellAreaN(M,N)=sqrt((x(M,N)-M_HALF*(x(M-1,N)+x(M,N)))**2+(y(M,N)-M_HALF*(y(M-1,N)+y(M,N)))**2)
        CellAreaS(M,N)=0.25d0*(x(M,N)**2+M_TWO*x(M,N)*x(M,N-1)-M_TWO*x(M,N)*x(M-1,N-1)-M_TWO*x(M-1,N)*x(M,N)+x(M,N-1)**2 &
                -M_TWO*x(M-1,N-1)*x(M,N-1)-M_TWO*x(M-1,N)*x(M,N-1)+x(M-1,N-1)**2+M_TWO*x(M-1,N)*x(M-1,N-1)+x(M-1,N)**2 &
                +y(M,N)**2+M_TWO*y(M,N)*y(M,N-1)-M_TWO*y(M,N)*y(M-1,N-1)-M_TWO*y(M-1,N)*y(M,N)+y(M,N-1)**2 &
                -M_TWO*y(M-1,N-1)*y(M,N-1)-M_TWO*y(M-1,N)*y(M,N-1)+y(M-1,N-1)**2+M_TWO*y(M-1,N)*y(M-1,N-1) &
                +y(M-1,N)**2)**M_HALF
        CellAreaW(M,N)=0.25d0*(x(M-1,N)**2+M_TWO*x(M-1,N)*x(M,N)-M_TWO*x(M-1,N)*x(M-1,N-1)-M_TWO*x(M-1,N)*x(M,N-1)+x(M,N)**2 &
                -M_TWO*x(M,N)*x(M-1,N-1)-M_TWO*x(M,N)*x(M,N-1)+x(M-1,N-1)**2+M_TWO*x(M-1,N-1)*x(M,N-1)+x(M,N-1)**2 &
                +y(M-1,N)**2+M_TWO*y(M-1,N)*y(M,N)-M_TWO*y(M-1,N)*y(M-1,N-1)-M_TWO*y(M-1,N)*y(M,N-1)+y(M,N)**2 &
                -M_TWO*y(M,N)*y(M-1,N-1)-M_TWO*y(M,N)*y(M,N-1)+y(M-1,N-1)**2+M_TWO*y(M-1,N-1)*y(M,N-1) &
                +y(M,N-1)**2)**M_HALF

        DistN(M,N)=M_ZERO
        DistS(M,N)=sqrt((x(M,N)-x(M,N-1))**2+(y(M,N)-y(M,N-1))**2)
        DistE(M,N)=M_ZERO
        DistW(M,N)=sqrt((x(M,N)-x(M-1,N))**2+(y(M,N)-y(M-1,N))**2)

   !South-east
        CellAreaS(M,1)=sqrt(( x(M,1)-M_HALF*(x(M,1)+x(M-1,1)) )**2+( y(M,1)-M_HALF*(y(M,1)+y(M-1,1)) )**2)
        CellAreaE(M,1)=sqrt( ( x(M,1) - M_HALF*(x(M,1)+x(M,2)) )**2 + ( y(M,1) - M_HALF*(y(M,1)+y(M,2)) )**2 )
        CellAreaN(M,1)=0.25d0*(x(M,1)**2+M_TWO*x(M,1)*x(M,2)-M_TWO*x(M-1,1)*x(M,1)-M_TWO*x(M-1,2)*x(M,1)+x(M,2)**2 &
                -M_TWO*x(M-1,1)*x(M,2)-M_TWO*x(M-1,2)*x(M,2)+x(M-1,1)**2+M_TWO*x(M-1,1)*x(M-1,2)+x(M-1,2)**2 &
                +y(M,1)**2+M_TWO*y(M,1)*y(M,2)-M_TWO*y(M-1,1)*y(M,1)-M_TWO*y(M-1,2)*y(M,1)+y(M,2)**2 &
                -M_TWO*y(M-1,1)*y(M,2)-M_TWO*y(M-1,2)*y(M,2)+y(M-1,1)**2+M_TWO*y(M-1,1)*y(M-1,2) &
                +y(M-1,2)**2)**M_HALF
        CellAreaW(M,1)=0.25d0*(x(M-1,1)**2-M_TWO*x(M-1,1)*x(M-1,2)+M_TWO*x(M-1,1)*x(M,1)-M_TWO*x(M-1,1)*x(M,2) &
                +x(M-1,2)**2-M_TWO*x(M-1,2)*x(M,1)+M_TWO*x(M-1,2)*x(M,2)+x(M,1)**2-M_TWO*x(M,1)*x(M,2)+x(M,2)**2 &
                +y(M-1,1)**2-M_TWO*y(M-1,1)*y(M-1,2)+M_TWO*y(M-1,1)*y(M,1)-M_TWO*y(M-1,1)*y(M,2)+y(M-1,2)**2 &
                -M_TWO*y(M-1,2)*y(M,1)+M_TWO*y(M-1,2)*y(M,2)+y(M,1)**2-M_TWO*y(M,1)*y(M,2)+y(M,2)**2)**M_HALF

         DistN(M,1)=sqrt( (x(M,2)-x(M,1))**2 + (y(M,2)-y(M,1))**2 )
         DistS(M,1)=M_ZERO
         DistE(M,1)=M_ZERO !sqrt( (x(i+1,1)-x(i,1))**2 + (y(i+1,1)-y(i,1))**2 )
         DistW(M,1)=sqrt( (x(M,1)-x(M-1,1))**2 + (y(M,1)-y(M-1,1))**2 )


  !South-west

        CellAreaS(1,1)=sqrt((x(1,1)-M_HALF*(x(2,1)+x(1,1)))**2+(y(1,1)-M_HALF*(y(2,1)+y(1,1)))**2)
        CellAreaW(1,1)=sqrt((x(1,1)-M_HALF*(x(1,2)+x(1,1)))**2+(y(1,1)-M_HALF*(y(1,2)+y(1,1)))**2)
        CellAreaN(1,1)=sqrt((0.25d0*(x(2,2)+x(2,1)+x(1,2)+x(1,1))-M_HALF*(x(1,1)+x(1,2)))**2 &
                +(0.25d0*(y(2,2)+y(2,1)+y(1,2)+y(1,1))-M_HALF*(y(1,2)+y(1,1)))**2) !sqrt((M_HALF*(x(1,2)+x(1,1))-0.25d0*(x(2,2)+x(2,1)+x(1,1)+x(1,2)))**M_TWO+(M_HALF*(y(1,2)+y(1,1))-0.25d0*(y(2,2)+y(2,1)+y(1,1)+y(1,2)))**M_TWO) !0.25d0*(x(1,1)**2+M_TWO*x(1,1)*x(1,2)-M_TWO*x(1,1)*x(2,1)-M_TWO*x(1,1)*x(2,2)+x(1,2)**2-M_TWO*x(1,2)*x(2,1)-M_TWO*x(1,2)*x(2,2)+x(2,1)**2+M_TWO*x(2,1)*x(2,2)+x(2,2)**2+y(1,1)**2+M_TWO*y(1,1)*y(1,2)-M_TWO*y(1,1)*y(2,1)-M_TWO*y(1,1)*y(2,2)+y(1,2)**2-M_TWO*y(1,2)*y(2,1)-M_TWO*y(1,2)*y(2,2)+y(2,1)**2+M_TWO*y(2,1)*y(2,2)+y(2,2)**2)**(M_HALF);
        CellAreaE(1,1)=sqrt((0.25d0*(x(2,2)+x(2,1)+x(1,2)+x(1,1))-M_HALF*(x(2,1)+x(1,1)))**2 &
                +(0.25d0*(y(2,2)+y(2,1)+y(1,2)+y(1,1))-M_HALF*(y(2,1)+y(1,1)))**2) !0.25d0*(x(1,1)**2-M_TWO*x(1,1)*x(1,2)+M_TWO*x(1,1)*x(2,1)-M_TWO*x(1,1)*x(2,2)+x(1,2)**2-M_TWO*x(1,2)*x(2,1)+M_TWO*x(1,2)*x(2,2)+x(2,1)**2-M_TWO*x(2,1)*x(2,2)+x(2,2)**2+y(1,1)**2-M_TWO*y(1,1)*y(1,2)+M_TWO*y(1,1)*y(2,1)-M_TWO*y(1,1)*y(2,2)+y(1,2)**2-M_TWO*y(1,2)*y(2,1)+M_TWO*y(1,2)*y(2,2)+y(2,1)**2-M_TWO*y(2,1)*y(2,2)+y(2,2)**2)**(M_HALF)

        DistN(1,1)=sqrt((x(1,2)-x(1,1))**2+(y(1,2)-y(1,1))**2)
        DistS(1,1)=M_ZERO
        DistE(1,1)=sqrt((x(2,1)-x(1,1))**2+(y(2,1)-y(1,1))**2)
        DistW(1,1)=M_ZERO

  !North-West
        CellAreaN(1,N)=sqrt((x(1,N)-M_HALF*(x(1,N)+x(2,N)))**2+(y(1,N)-M_HALF*(y(1,N)+y(2,N)))**2)
        CellAreaW(1,N)=sqrt((x(1,N)-M_HALF*(x(1,N)+x(1,N-1)))**2+((y(1,N)-M_HALF*(y(1,N)+y(1,N-1))))**2)
        CellAreaS(1,N)=sqrt((0.25d0*(x(2,N-1)+x(1,N-1)+x(2,N)+x(1,N))-M_HALF*(x(1,N-1)+x(1,N)))**2 &
                      +(0.25d0*(y(2,N-1)+y(1,N-1)+y(2,N)+y(1,N))-M_HALF*(y(1,N-1)+y(1,N)))**2) !sqrt((M_HALF*(x(1,N-1)+x(1,N))-0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)))**2+(M_HALF*(y(1,N-1)+y(1,N))-0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)))**2)
        CellAreaE(1,N)=sqrt((M_HALF*(x(2,N)+x(1,N))-0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)))**2 &
                      +(M_HALF*(y(2,N)+y(1,N))-0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)))**2) !sqrt((0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1))-(M_HALF*(x(2,N)+x(1,N))))**2+(0.25d0*(y(1,N)+y(1,N-1)+y(2,N)+y(2,N-1))-M_HALF*(y(2,N)+y(1,N)))**2)

         DistN(1,N)=M_ZERO !sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
         DistS(1,N)=sqrt((x(1,N)-x(1,N-1))**2+(y(1,N)-y(1,N-1))**2)
         DistE(1,N)=sqrt((x(2,N)-x(1,N))**2+(y(2,N)-y(1,N))**2)
         DistW(1,N)=M_ZERO !sqrt((x(i,N)-x(i-1,N))**2+(y(i,N)-y(i-1,N))**2)


  ! distance for cross diffusion, equal to cell border "areas" in 2D
  DistDualN(:,:)=CellAreaN(:,:)
  DistDualS(:,:)=CellAreaS(:,:)
  DistDualE(:,:)=CellAreaE(:,:)
  DistDualW(:,:)=CellAreaW(:,:)


end subroutine compute_distances

subroutine compute_norm_tan_curv(M, N, x, y, NormalN, NormalS, NormalE, NormalW, &
                                 TangentNx, TangentNy, TangentSx, TangentSy, &
                                 TangentEx, TangentEy, TangentWx, TangentWy, CurviNx, CurviNy, &
                                 CurviSx, CurviSy, CurviEx, CurviEy, CurviWx, CurviWy )
  use Maths_m
  use Types_m
  implicit none

  integer, intent(in)    :: M, N
  real(8), intent(in)    :: x(1:M, 1:N), y(1:M, 1:N)                 ! needle position indexes
  type(VectorField), intent(inout) :: NormalN, NormalS, NormalW, NormalE ! normal to quadrangle elements

  real(8), intent(inout) ::  TangentWx(1:M,1:N), TangentWy(1:M,1:N), &                ! Tangent to quadrangle elements
                             TangentEx(1:M,1:N), TangentEy(1:M,1:N), &
                             TangentNx(1:M,1:N), TangentNy(1:M,1:N), &
                             TangentSx(1:M,1:N), TangentSy(1:M,1:N), &
                             CurviWx(1:M,1:N), CurviWy(1:M,1:N), &                 ! Unit vector between cell centers
                             CurviEx(1:M,1:N), CurviEy(1:M,1:N), &
                             CurviNx(1:M,1:N), CurviNy(1:M,1:N), &
                             CurviSx(1:M,1:N), CurviSy(1:M,1:N)

  real(8) :: Tangent, Normal
  integer :: i, j

    do j=2,N-1
      do i=2,M-1

        !TODO: Optimise by defining a stencil object

        NormalN%x(i,j)=Normal(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                     0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)),1) !4d0*(0.25d0*y(i-1,j)+0.25d0*y(i-1,j+1)-0.25d0*y(i+1,j+1)-0.25d0*y(i+1,j))/(x(i-1,j)**2+M_TWO*x(i-1,j)*x(i-1,j+1)-M_TWO*x(i-1,j)*x(i+1,j+1)-M_TWO*x(i-1,j)*x(i+1,j)+x(i-1,j+1)**2-M_TWO*x(i-1,j+1)*x(i+1,j+1)-M_TWO*x(i-1,j+1)*x(i+1,j)+x(i+1,j+1)**2+M_TWO*x(i+1,j+1)*x(i+1,j)+x(i+1,j)**2+y(i-1,j)**2+M_TWO*y(i-1,j)*y(i-1,j+1)-M_TWO*y(i-1,j)*y(i+1,j+1)-M_TWO*y(i-1,j)*y(i+1,j)+y(i-1,j+1)**2-M_TWO*y(i-1,j+1)*y(i+1,j+1)-M_TWO*y(i-1,j+1)*y(i+1,j)+y(i+1,j+1)**2+M_TWO*y(i+1,j+1)*y(i+1,j)+y(i+1,j)**2)**(M_HALF)
        NormalN%y(i,j)=Normal(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                     0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)),2) !NormalNy(i,j)=-4d0*(0.25d0*x(i-1,j)+0.25d0*x(i-1,j+1)-0.25d0*x(i+1,j+1)-0.25d0*x(i+1,j))/(x(i-1,j)**2+M_TWO*x(i-1,j)*x(i-1,j+1)-M_TWO*x(i-1,j)*x(i+1,j+1)-M_TWO*x(i-1,j)*x(i+1,j)+x(i-1,j+1)**2-M_TWO*x(i-1,j+1)*x(i+1,j+1)-M_TWO*x(i-1,j+1)*x(i+1,j)+x(i+1,j+1)**2+M_TWO*x(i+1,j+1)*x(i+1,j)+x(i+1,j)**2+y(i-1,j)**2+M_TWO*y(i-1,j)*y(i-1,j+1)-M_TWO*y(i-1,j)*y(i+1,j+1)-M_TWO*y(i-1,j)*y(i+1,j)+y(i-1,j+1)**2-M_TWO*y(i-1,j+1)*y(i+1,j+1)-M_TWO*y(i-1,j+1)*y(i+1,j)+y(i+1,j+1)**2+M_TWO*y(i+1,j+1)*y(i+1,j)+y(i+1,j)**2)**(M_HALF)
        NormalS%x(i,j)=Normal(0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)),1) ! -4d0*(0.25d0*y(i-1,j-1)+0.25d0*y(i-1,j)-0.25d0*y(i+1,j)-0.25d0*y(i+1,j-1))/(x(i-1,j-1)**2+M_TWO*x(i-1,j-1)*x(i-1,j)-M_TWO*x(i-1,j-1)*x(i+1,j)-M_TWO*x(i-1,j-1)*x(i+1,j-1)+x(i-1,j)**2-M_TWO*x(i-1,j)*x(i+1,j)-M_TWO*x(i-1,j)*x(i+1,j-1)+x(i+1,j)**2+M_TWO*x(i+1,j)*x(i+1,j-1)+x(i+1,j-1)**2+y(i-1,j-1)**2+M_TWO*y(i-1,j-1)*y(i-1,j)-M_TWO*y(i-1,j-1)*y(i+1,j)-M_TWO*y(i-1,j-1)*y(i+1,j-1)+y(i-1,j)**2-M_TWO*y(i-1,j)*y(i+1,j)-M_TWO*y(i-1,j)*y(i+1,j-1)+y(i+1,j)**2+M_TWO*y(i+1,j)*y(i+1,j-1)+y(i+1,j-1)**2)**(M_HALF)
        NormalS%y(i,j)=Normal(0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)),2) ! 4d0*(0.25d0*x(i-1,j-1)+0.25d0*x(i-1,j)-0.25d0*x(i+1,j)-0.25d0*x(i+1,j-1))/(x(i-1,j-1)**2+M_TWO*x(i-1,j-1)*x(i-1,j)-M_TWO*x(i-1,j-1)*x(i+1,j)-M_TWO*x(i-1,j-1)*x(i+1,j-1)+x(i-1,j)**2-M_TWO*x(i-1,j)*x(i+1,j)-M_TWO*x(i-1,j)*x(i+1,j-1)+x(i+1,j)**2+M_TWO*x(i+1,j)*x(i+1,j-1)+x(i+1,j-1)**2+y(i-1,j-1)**2+M_TWO*y(i-1,j-1)*y(i-1,j)-M_TWO*y(i-1,j-1)*y(i+1,j)-M_TWO*y(i-1,j-1)*y(i+1,j-1)+y(i-1,j)**2-M_TWO*y(i-1,j)*y(i+1,j)-M_TWO*y(i-1,j)*y(i+1,j-1)+y(i+1,j)**2+M_TWO*y(i+1,j)*y(i+1,j-1)+y(i+1,j-1)**2)**(M_HALF)
        NormalE%x(i,j)=Normal(0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j+1)+x(i+1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i+1,j+1)+y(i+1,j)+y(i,j+1)+y(i,j)), 1) ! -4d0*(0.25d0*y(i+1,j-1)+0.25d0*y(i,j-1)-0.25d0*y(i+1,j+1)-0.25d0*y(i,j+1))/(x(i+1,j-1)**2+M_TWO*x(i+1,j-1)*x(i,j-1)-M_TWO*x(i+1,j-1)*x(i+1,j+1)-M_TWO*x(i+1,j-1)*x(i,j+1)+x(i,j-1)**2-M_TWO*x(i,j-1)*x(i+1,j+1)-M_TWO*x(i,j-1)*x(i,j+1)+x(i+1,j+1)**2+M_TWO*x(i+1,j+1)*x(i,j+1)+x(i,j+1)**2+y(i+1,j-1)**2+M_TWO*y(i+1,j-1)*y(i,j-1)-M_TWO*y(i+1,j-1)*y(i+1,j+1)-M_TWO*y(i+1,j-1)*y(i,j+1)+y(i,j-1)**2-M_TWO*y(i,j-1)*y(i+1,j+1)-M_TWO*y(i,j-1)*y(i,j+1)+y(i+1,j+1)**2+M_TWO*y(i+1,j+1)*y(i,j+1)+y(i,j+1)**2)**(M_HALF)
        NormalE%y(i,j)=Normal(0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j+1)+x(i+1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i+1,j+1)+y(i+1,j)+y(i,j+1)+y(i,j)), 2) ! 4d0*(0.25d0*x(i+1,j-1)+0.25d0*x(i,j-1)-0.25d0*x(i+1,j+1)-0.25d0*x(i,j+1))/(x(i+1,j-1)**2+M_TWO*x(i+1,j-1)*x(i,j-1)-M_TWO*x(i+1,j-1)*x(i+1,j+1)-M_TWO*x(i+1,j-1)*x(i,j+1)+x(i,j-1)**2-M_TWO*x(i,j-1)*x(i+1,j+1)-M_TWO*x(i,j-1)*x(i,j+1)+x(i+1,j+1)**2+M_TWO*x(i+1,j+1)*x(i,j+1)+x(i,j+1)**2+y(i+1,j-1)**2+M_TWO*y(i+1,j-1)*y(i,j-1)-M_TWO*y(i+1,j-1)*y(i+1,j+1)-M_TWO*y(i+1,j-1)*y(i,j+1)+y(i,j-1)**2-M_TWO*y(i,j-1)*y(i+1,j+1)-M_TWO*y(i,j-1)*y(i,j+1)+y(i+1,j+1)**2+M_TWO*y(i+1,j+1)*y(i,j+1)+y(i,j+1)**2)**(M_HALF)
        NormalW%x(i,j)=Normal(0.25d0*(x(i-1,j+1)+x(i-1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i-1,j)+y(i,j+1)+y(i,j)), &
                     0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 1) ! 4d0*(0.25d0*y(i-1,j-1)+0.25d0*y(i,j-1)-0.25d0*y(i,j+1)-0.25d0*y(i-1,j+1))/(x(i-1,j-1)**2+M_TWO*x(i-1,j-1)*x(i,j-1)-M_TWO*x(i-1,j-1)*x(i,j+1)-M_TWO*x(i-1,j-1)*x(i-1,j+1)+x(i,j-1)**2-M_TWO*x(i,j-1)*x(i,j+1)-M_TWO*x(i,j-1)*x(i-1,j+1)+x(i,j+1)**2+M_TWO*x(i,j+1)*x(i-1,j+1)+x(i-1,j+1)**2+y(i-1,j-1)**2+M_TWO*y(i-1,j-1)*y(i,j-1)-M_TWO*y(i-1,j-1)*y(i,j+1)-M_TWO*y(i-1,j-1)*y(i-1,j+1)+y(i,j-1)**2-M_TWO*y(i,j-1)*y(i,j+1)-M_TWO*y(i,j-1)*y(i-1,j+1)+y(i,j+1)**2+M_TWO*y(i,j+1)*y(i-1,j+1)+y(i-1,j+1)**2)**(M_HALF)
        NormalW%y(i,j)=Normal(0.25d0*(x(i-1,j+1)+x(i-1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i-1,j)+y(i,j+1)+y(i,j)), &
                     0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 2)! -4d0*(0.25d0*x(i-1,j-1)+0.25d0*x(i,j-1)-0.25d0*x(i,j+1)-0.25d0*x(i-1,j+1))/(x(i-1,j-1)**2+M_TWO*x(i-1,j-1)*x(i,j-1)-M_TWO*x(i-1,j-1)*x(i,j+1)-M_TWO*x(i-1,j-1)*x(i-1,j+1)+x(i,j-1)**2-M_TWO*x(i,j-1)*x(i,j+1)-M_TWO*x(i,j-1)*x(i-1,j+1)+x(i,j+1)**2+M_TWO*x(i,j+1)*x(i-1,j+1)+x(i-1,j+1)**2+y(i-1,j-1)**2+M_TWO*y(i-1,j-1)*y(i,j-1)-M_TWO*y(i-1,j-1)*y(i,j+1)-M_TWO*y(i-1,j-1)*y(i-1,j+1)+y(i,j-1)**2-M_TWO*y(i,j-1)*y(i,j+1)-M_TWO*y(i,j-1)*y(i-1,j+1)+y(i,j+1)**2+M_TWO*y(i,j+1)*y(i-1,j+1)+y(i-1,j+1)**2)**(M_HALF)


        TangentNx(i,j)=-Tangent(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)),0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                        0.25d0*(x(i,j)+x(i,j+1)+x(i-1,j+1)+x(i-1,j)),0.25d0*(y(i,j)+y(i,j+1)+y(i-1,j+1)+y(i-1,j)), 1)
        TangentNy(i,j)=-Tangent(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)),0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                        0.25d0*(x(i,j)+x(i,j+1)+x(i-1,j+1)+x(i-1,j)),0.25d0*(y(i,j)+y(i,j+1)+y(i-1,j+1)+y(i-1,j)), 2)
        TangentSx(i,j)=Tangent(0.25d0*(x(i,j)+x(i,j-1)+x(i+1,j-1)+x(i+1,j)),0.25d0*(y(i,j)+y(i,j-1)+y(i+1,j-1)+y(i+1,j)), &
                        0.25d0*(x(i-1,j)+x(i-1,j-1)+x(i,j-1)+x(i,j)),0.25d0*(y(i-1,j)+y(i-1,j-1)+y(i,j-1)+y(i,j)), 1)
        TangentSy(i,j)=Tangent(0.25d0*(x(i,j)+x(i,j-1)+x(i+1,j-1)+x(i+1,j)),0.25d0*(y(i,j)+y(i,j-1)+y(i+1,j-1)+y(i+1,j)), &
                        0.25d0*(x(i-1,j)+x(i-1,j-1)+x(i,j-1)+x(i,j)),0.25d0*(y(i-1,j)+y(i-1,j-1)+y(i,j-1)+y(i,j)), 2)
        TangentEx(i,j)=Tangent(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                        0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j-1)+x(i,j-1)),0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j-1)+y(i,j-1)), 1)
        TangentEy(i,j)=Tangent(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                        0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j-1)+x(i,j-1)),0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j-1)+y(i,j-1)), 2)
        TangentWx(i,j)=-Tangent(0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)),0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1)), &
                        0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)),0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j)), 1)
        TangentWy(i,j)=-Tangent(0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)),0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1)), &
                        0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)),0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j)), 2)

        CurviNx(i,j)=-Tangent(x(i,j),y(i,j),x(i,j+1),y(i,j+1),1)
        CurviNy(i,j)=-Tangent(x(i,j),y(i,j),x(i,j+1),y(i,j+1),2)
        CurviSx(i,j)=-Tangent(x(i,j),y(i,j),x(i,j-1),y(i,j-1),1)
        CurviSy(i,j)=-Tangent(x(i,j),y(i,j),x(i,j-1),y(i,j-1),2)
        CurviEx(i,j)=-Tangent(x(i,j),y(i,j),x(i+1,j),y(i+1,j),1)
        CurviEy(i,j)=-Tangent(x(i,j),y(i,j),x(i+1,j),y(i+1,j),2)
        CurviWx(i,j)=-Tangent(x(i,j),y(i,j),x(i-1,j),y(i-1,j),1)
        CurviWy(i,j)=-Tangent(x(i,j),y(i,j),x(i-1,j),y(i-1,j),2)

      end do
   end do

   do i=2,M-1

        ! NORTH
         NormalE%x(i,N)=-Normal(M_HALF*(x(i,N)+x(i+1,N)), M_HALF*(y(i,N)+y(i+1,N)), &
                               0.25d0*(x(i+1,N)+x(i,N)+x(i+1,N-1)+x(i,N-1)), &
                               0.25d0*(y(i+1,N)+y(i,N)+y(i+1,N-1)+y(i,N-1)), 1)
         NormalE%y(i,N)=-Normal(M_HALF*(x(i,N)+x(i+1,N)), M_HALF*(y(i,N)+y(i+1,N)), &
                               0.25d0*(x(i+1,N)+x(i,N)+x(i+1,N-1)+x(i,N-1)), &
                               0.25d0*(y(i+1,N)+y(i,N)+y(i+1,N-1)+y(i,N-1)), 2)
         NormalW%x(i,N)=Normal(M_HALF*(x(i-1,N)+x(i,N)), M_HALF*(y(i-1,N)+y(i,N)), &
                              0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), &
                              0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),1)
         NormalW%y(i,N)=Normal(M_HALF*(x(i-1,N)+x(i,N)), M_HALF*(y(i-1,N)+y(i,N)), &
                              0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), &
                              0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),2)
         NormalN%x(i,N)=Normal(M_HALF*(x(i,N)+x(i+1,N)),M_HALF*(y(i,N)+y(i+1,N)), &
                              M_HALF*(x(i,N)+x(i-1,N)), M_HALF*(y(i,N)+y(i-1,N)),1)
         NormalN%y(i,N)=Normal(M_HALF*(x(i,N)+x(i+1,N)),M_HALF*(y(i,N)+y(i+1,N)), &
                              M_HALF*(x(i,N)+x(i-1,N)), M_HALF*(y(i,N)+y(i-1,N)),2)
         NormalS%x(i,N)=-Normal(0.25d0*(x(i+1,N-1)+x(i+1,N)+x(i,N)+x(i,N-1)), &
                               0.25d0*(y(i+1,N-1)+y(i+1,N)+y(i,N)+y(i,N-1)), &
                               0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), &
                               0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),1)
         NormalS%y(i,N)=-Normal(0.25d0*(x(i+1,N-1)+x(i+1,N)+x(i,N)+x(i,N-1)), &
                               0.25d0*(y(i+1,N-1)+y(i+1,N)+y(i,N)+y(i,N-1)), &
                               0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), &
                               0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),2)

         TangentWx(i,N)=-Tangent(M_HALF*(x(i-1,N)+x(i,N)),M_HALF*(y(i-1,N)+y(i,N)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N)+x(i,N-1)), &
                  0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N)+y(i,N-1)),1)
         TangentWy(i,N)=-Tangent(M_HALF*(x(i-1,N)+x(i,N)),M_HALF*(y(i-1,N)+y(i,N)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N)+x(i,N-1)), &
                  0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N)+y(i,N-1)),2)
         TangentEx(i,N)=-Tangent(0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)), &
                  M_HALF*(x(i+1,N)+x(i,N)), M_HALF*(y(i+1,N)+y(i,N)), 1)
         TangentEy(i,N)=-Tangent(0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)), &
                  M_HALF*(x(i+1,N)+x(i,N)), M_HALF*(y(i+1,N)+y(i,N)), 2)
         TangentNx(i,N)=-Tangent(M_HALF*(x(i,N)+x(i+1,N)),M_HALF*(y(i,N)+y(i+1,N)), &
                                 M_HALF*(x(i,N)+x(i-1,N)),M_HALF*(y(i,N)+y(i-1,N)),1)
         TangentNy(i,N)=-Tangent(M_HALF*(x(i,N)+x(i+1,N)),M_HALF*(y(i,N)+y(i+1,N)), &
                                 M_HALF*(x(i,N)+x(i-1,N)),M_HALF*(y(i,N)+y(i-1,N)),2)
         TangentSx(i,N)=-Tangent(0.25d0*(x(i-1,N-1)+x(i,N-1)+x(i-1,N)+x(i,N)), 0.25d0*(y(i-1,N-1)+y(i,N-1)+y(i-1,N)+y(i,N)),  &
                  0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)),1)
         TangentSy(i,N)=-Tangent(0.25d0*(x(i-1,N-1)+x(i,N-1)+x(i-1,N)+x(i,N)), 0.25d0*(y(i-1,N-1)+y(i,N-1)+y(i-1,N)+y(i,N)), &
                  0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)),2)

         CurviNx(i,N)=CurviNx(i,N-1) !-Tangent(x(i,j),y(i,j),x(i,j+1),y(i,j+1),1)
         CurviNy(i,N)=CurviNy(i,N-1) !-Tangent(x(i,j),y(i,j),x(i,j+1),y(i,j+1),2)
         CurviSx(i,N)=-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),1)
         CurviSy(i,N)=-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),2)
         CurviEx(i,N)=-Tangent(x(i,N),y(i,N),x(i+1,N),y(i+1,N),1)
         CurviEy(i,N)=-Tangent(x(i,N),y(i,N),x(i+1,N),y(i+1,N),2)
         CurviWx(i,N)=-Tangent(x(i,N),y(i,N),x(i-1,N),y(i-1,N),1)
         CurviWy(i,N)=-Tangent(x(i,N),y(i,N),x(i-1,N),y(i-1,N),2)


        ! SOUTH
         NormalE%x(i,1)=-Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), &
                M_HALF*(x(i,1)+x(i+1,1)),M_HALF*(y(i,1)+y(i+1,1)),1)
         NormalE%y(i,1)=-Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), &
                M_HALF*(x(i,1)+x(i+1,1)),M_HALF*(y(i,1)+y(i+1,1)),2)
         NormalW%x(i,1)=Normal(0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)), &
                M_HALF*(x(i-1,1)+x(i,1)),M_HALF*(y(i-1,1)+y(i,1)),1)
         NormalW%y(i,1)=Normal(0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)), &
                M_HALF*(x(i-1,1)+x(i,1)),M_HALF*(y(i-1,1)+y(i,1)),2)
         NormalN%x(i,1)=Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), &
                0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),1)
         NormalN%y(i,1)=Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), &
                0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),2)
         NormalS%x(i,1)=Normal(M_HALF*(x(i-1,1)+x(i,1)),M_HALF*(y(i-1,1)+y(i,1)), &
                               M_HALF*(x(i+1,1)+x(i,1)),M_HALF*(y(i+1,1)+y(i,1)),1)
         NormalS%y(i,1)=Normal(M_HALF*(x(i-1,1)+x(i,1)),M_HALF*(y(i-1,1)+y(i,1)), &
                               M_HALF*(x(i+1,1)+x(i,1)),M_HALF*(y(i+1,1)+y(i,1)),2)

         TangentNx(i,1)=Tangent(0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)), &
                  0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),1)
         TangentNy(i,1)=Tangent(0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)), &
                  0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),2)
         TangentSx(i,1)=-Tangent(M_HALF*(x(i-1,1)+x(i,1)),M_HALF*(y(i-1,1)+y(i,1)), &
                                 M_HALF*(x(i,1)+x(i+1,1)),M_HALF*(y(i,1)+y(i+1,1)),1)
         TangentSy(i,1)=-Tangent(M_HALF*(x(i-1,1)+x(i,1)),M_HALF*(y(i-1,1)+y(i,1)), &
                                 M_HALF*(x(i,1)+x(i+1,1)),M_HALF*(y(i,1)+y(i+1,1)),2)
         TangentEx(i,1)=-Tangent(M_HALF*(x(i+1,1)+x(i,1)),M_HALF*(y(i+1,1)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)), &
                  0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),1)
         TangentEy(i,1)=-Tangent(M_HALF*(x(i+1,1)+x(i,1)),M_HALF*(y(i+1,1)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)), &
                  0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),2)
         TangentWx(i,1)=Tangent(M_HALF*(x(i-1,1)+x(i,1)),M_HALF*(y(i-1,1)+y(i,1)),0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)), &
                  0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),1)
         TangentWy(i,1)=Tangent(M_HALF*(x(i-1,1)+x(i,1)),M_HALF*(y(i-1,1)+y(i,1)),0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)), &
                  0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),2)

         CurviNx(i,1)=-Tangent(x(i,1),y(i,1),x(i,2),y(i,2),1)
         CurviNy(i,1)=-Tangent(x(i,1),y(i,1),x(i,2),y(i,2),2)
         CurviSx(i,1)=CurviSx(i,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),1)
         CurviSy(i,1)=CurviSy(i,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),2)
         CurviEx(i,1)=-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),1)
         CurviEy(i,1)=-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),2)
         CurviWx(i,1)=-Tangent(x(i,1),y(i,1),x(i-1,1),y(i-1,1),1)
         CurviWy(i,1)=-Tangent(x(i,1),y(i,1),x(i-1,1),y(i-1,1),2)


      end do

        do j=2,N-1
!
          ! WEST
          NormalE%x(1,j)=-Normal(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)), &
                        0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),1)
          NormalE%y(1,j)=-Normal(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)), &
                        0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),2)
          NormalN%x(1,j)=Normal(0.25d0*(x(2,j+1)+x(1,j+1)+x(2,j)+x(1,j)),0.25d0*(y(2,j+1)+y(1,j+1)+y(2,j)+y(1,j)), &
                      M_HALF*(x(1,j+1)+x(1,j)),M_HALF*(y(1,j+1)+y(1,j)),1)
          NormalN%y(1,j)=Normal(0.25d0*(x(2,j+1)+x(1,j+1)+x(2,j)+x(1,j)),0.25d0*(y(2,j+1)+y(1,j+1)+y(2,j)+y(1,j)), &
                      M_HALF*(x(1,j+1)+x(1,j)),M_HALF*(y(1,j+1)+y(1,j)),2)
          NormalS%x(1,j)=-Normal(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)), &
                        M_HALF*(x(1,j-1)+x(1,j)),M_HALF*(y(1,j-1)+y(1,j)),1)
          NormalS%y(1,j)=-Normal(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)), &
                        M_HALF*(x(1,j-1)+x(1,j)),M_HALF*(y(1,j-1)+y(1,j)),2)
          NormalW%x(1,j)=Normal(M_HALF*(x(1,j+1)+x(1,j)),M_HALF*(y(1,j+1)+y(1,j)),M_HALF*(x(1,j-1)+x(1,j)), &
                        M_HALF*(y(1,j-1)+y(1,j)),1)
          NormalW%y(1,j)=Normal(M_HALF*(x(1,j+1)+x(1,j)),M_HALF*(y(1,j+1)+y(1,j)),M_HALF*(x(1,j-1)+x(1,j)), &
                        M_HALF*(y(1,j-1)+y(1,j)),2)

!           write(*,*) "[Debug]", NormalWx(1,j)**2+NormalWy(1,j)**2

          TangentNx(1,j)=-Tangent(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) &
                  +y(1,j)),M_HALF*(x(1,j+1)+x(1,j)),M_HALF*(y(1,j+1)+y(1,j)),1)
          TangentNy(1,j)=-Tangent(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) &
                  +y(1,j)),M_HALF*(x(1,j+1)+x(1,j)),M_HALF*(y(1,j+1)+y(1,j)),2)
          TangentSx(1,j)=-Tangent(M_HALF*(x(1,j-1)+x(1,j)),M_HALF*(y(1,j-1)+y(1,j)),0.25d0*(x(1,j-1)+x(1,j) &
                  +x(2,j-1)+x(2,j)),0.25d0*(y(1,j-1)+y(1,j)+y(2,j-1)+y(2,j)),1)
          TangentSy(1,j)=-Tangent(M_HALF*(x(1,j-1)+x(1,j)),M_HALF*(y(1,j-1)+y(1,j)),0.25d0*(x(1,j-1)+x(1,j) &
                  +x(2,j-1)+x(2,j)),0.25d0*(y(1,j-1)+y(1,j)+y(2,j-1)+y(2,j)),2)
          TangentWx(1,j)=Tangent(M_HALF*(x(1,j-1)+x(1,j)),M_HALF*(y(1,j-1)+y(1,j)),M_HALF*(x(1,j+1)+x(1,j)), &
                  M_HALF*(y(1,j+1)+y(1,j)),1)
          TangentWy(1,j)=Tangent(M_HALF*(x(1,j-1)+x(1,j)),M_HALF*(y(1,j-1)+y(1,j)),M_HALF*(x(1,j+1)+x(1,j)), &
                  M_HALF*(y(1,j+1)+y(1,j)),2)
          TangentEx(1,j)=-Tangent(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1) &
                  +y(1,j)),0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) &
                  +y(1,j)),1)
          TangentEy(1,j)=-Tangent(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1) &
                  +y(1,j)),0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) &
                  +y(1,j)),2)

          CurviNx(1,j)=-Tangent(x(1,j),y(1,j),x(1,j+1),y(1,j+1),1)
          CurviNy(1,j)=-Tangent(x(1,j),y(1,j),x(1,j+1),y(1,j+1),2)
          CurviSx(1,j)=-Tangent(x(1,j),y(1,j),x(1,j-1),y(1,j-1),1)
          CurviSy(1,j)=-Tangent(x(1,j),y(1,j),x(1,j-1),y(1,j-1),2)
          CurviEx(1,j)=-Tangent(x(1,j),y(1,j),x(2,j),y(2,j),1)
          CurviEy(1,j)=-Tangent(x(1,j),y(1,j),x(2,j),y(2,j),2)
          CurviWx(1,j)=CurviWx(2,j) !-Tangent(x(1,j),y(1,j),x(i-1,N),y(i-1,N),1)
          CurviWy(1,j)=CurviWy(2,j) !-Tangent(x(1,j),y(1,j),x(i-1,N),y(i-1,N),2)


          ! EAST
          TangentNx(M,j)=-Tangent(M_HALF*(x(M,j+1)+x(M,j)),M_HALF*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) &
                          +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1)
          TangentNy(M,j)=-Tangent(M_HALF*(x(M,j+1)+x(M,j)),M_HALF*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) &
                          +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2)
          TangentEx(M,j)=-Tangent(M_HALF*(x(M,j-1)+x(M,j)),M_HALF*(y(M,j-1)+y(M,j)),M_HALF*(x(M,j+1)+x(M,j)), &
                          M_HALF*(y(M,j+1)+y(M,j)),1)
          TangentEy(M,j)=-Tangent(M_HALF*(x(M,j-1)+x(M,j)),M_HALF*(y(M,j-1)+y(M,j)),M_HALF*(x(M,j+1)+x(M,j)), &
                          M_HALF*(y(M,j+1)+y(M,j)),2)
          TangentSx(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) &
                          +y(M,j-1)+y(M,j)),M_HALF*(x(M,j-1)+x(M,j)),M_HALF*(y(M,j-1)+y(M,j)),1)
          TangentSy(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) &
                          +y(M,j-1)+y(M,j)),M_HALF*(x(M,j-1)+x(M,j)),M_HALF*(y(M,j-1)+y(M,j)),2)
          TangentWx(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) &
                          +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)), &
                          0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1)
          TangentWy(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) &
                          +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)), &
                          0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),2)

!           if(y(1,j)<0d0) then
!             TangentNx(M,j)=-TangentNx(M,j)
!             TangentNy(M,j)=-TangentNy(M,j)
!             TangentSx(M,j)=-TangentSx(M,j)
!             TangentSy(M,j)=-TangentSy(M,j)
!             TangentWx(M,j)=-TangentWx(M,j)
!             TangentWy(M,j)=-TangentWy(M,j)
!           end if

          NormalE%x(M,j)=-Normal(M_HALF*(x(M,j+1)+x(M,j)),M_HALF*(y(M,j+1)+y(M,j)),M_HALF*(x(M,j-1)+x(M,j)), &
                        M_HALF*(y(M,j-1)+y(M,j)),1);
          NormalE%y(M,j)=-Normal(M_HALF*(x(M,j+1)+x(M,j)),M_HALF*(y(M,j+1)+y(M,j)),M_HALF*(x(M,j-1)+x(M,j)), &
                        M_HALF*(y(M,j-1)+y(M,j)),2);
          NormalW%x(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) &
                       +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)), &
                      0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1)
          NormalW%y(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) &
          +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) &
          +y(M,j-1)+y(M,j)),2)
          NormalN%x(M,j)=Normal(M_HALF*(x(M,j+1)+x(M,j)),M_HALF*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1) &
                      +x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1)
          NormalN%y(M,j)=Normal(M_HALF*(x(M,j+1)+x(M,j)),M_HALF*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) &
                      +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2)
          NormalS%x(M,j)=-Normal(M_HALF*(x(M,j-1)+x(M,j)),M_HALF*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1) &
                      +x(M,j-1)+x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),1)
          NormalS%y(M,j)=-Normal(M_HALF*(x(M,j-1)+x(M,j)),M_HALF*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1) &
                      +x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),2)

          CurviNx(M,j)=-Tangent(x(M,j),y(M,j),x(M,j+1),y(M,j+1),1)
          CurviNy(M,j)=-Tangent(x(M,j),y(M,j),x(M,j+1),y(M,j+1),2)
          CurviSx(M,j)=-Tangent(x(M,j),y(M,j),x(M,j-1),y(M,j-1),1)
          CurviSy(M,j)=-Tangent(x(M,j),y(M,j),x(M,j-1),y(M,j-1),2)
          CurviEx(M,j)=CurviEx(M-1,j) !Tangent(x(M,j),y(M,j),x(i+1,N),y(i+1,N),1)
          CurviEy(M,j)=CurviEy(M-1,j) !-Tangent(x(M,j),y(M,j),x(i+1,N),y(i+1,N),2)
          CurviWx(M,j)=-Tangent(x(M,j),y(M,j),x(M-1,j),y(M-1,j),1)
          CurviWy(M,j)=-Tangent(x(M,j),y(M,j),x(M-1,j),y(M-1,j),2)

      end do

      !North-East
        NormalE%x(M,N)=-Normal(x(M,N),y(M,N),M_HALF*(x(M,N-1)+x(M,N)),M_HALF*(y(M,N-1)+y(M,N)),1)
        NormalE%y(M,N)=-Normal(x(M,N),y(M,N),M_HALF*(x(M,N-1)+x(M,N)),M_HALF*(y(M,N-1)+y(M,N)),2)
        NormalN%x(M,N)=Normal(x(M,N),y(M,N),M_HALF*(x(M-1,N)+x(M,N)),M_HALF*(y(M-1,N)+y(M,N)),1)
        NormalN%y(M,N)=Normal(x(M,N),y(M,N),M_HALF*(x(M-1,N)+x(M,N)),M_HALF*(y(M-1,N)+y(M,N)),2)
        NormalW%x(M,N)=Normal(M_HALF*(x(M-1, N)+x(M, N)),M_HALF*(y(M-1, N)+y(M, N)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),1)
        NormalW%y(M,N)=Normal(M_HALF*(x(M-1, N)+x(M, N)),M_HALF*(y(M-1, N)+y(M, N)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),2)
        NormalS%x(M,N)=-Normal(M_HALF*(x(M, N)+x(M, N-1)),M_HALF*(y(M, N)+y(M, N-1)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),1)
        NormalS%y(M,N)=-Normal(M_HALF*(x(M, N)+x(M, N-1)),M_HALF*(y(M, N)+y(M, N-1)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),2)

        TangentNx(M,N)=-Tangent(x(M,N),y(M,N),M_HALF*(x(M-1,N)+x(M,N)),M_HALF*(y(M-1,N)+y(M,N)),1)
        TangentNy(M,N)=-Tangent(x(M,N),y(M,N),M_HALF*(x(M-1,N)+x(M,N)),M_HALF*(y(M-1,N)+y(M,N)),2)
        TangentEx(M,N)=-Tangent(M_HALF*(x(M,N-1)+x(M,N)),M_HALF*(y(M,N-1)+y(M,N)),x(M,N),y(M,N),1)
        TangentEy(M,N)=-Tangent(M_HALF*(x(M,N-1)+x(M,N)),M_HALF*(y(M,N-1)+y(M,N)),x(M,N),y(M,N),2)
        TangentWx(M,N)=-Tangent(M_HALF*(x(M-1,N)+x(M,N)),M_HALF*(y(M-1,N)+y(M,N)),0.25d0*(x(M-1,N-1) &
                        +x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),1)
        TangentWy(M,N)=-Tangent(M_HALF*(x(M-1,N)+x(M,N)),M_HALF*(y(M-1,N)+y(M,N)),0.25d0*(x(M-1,N-1)+x(M-1,N) &
                        +x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),2)
        TangentSx(M,N)=-Tangent(0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N) &
                        +y(M,N-1)+y(M,N)),M_HALF*(x(M,N-1)+x(M,N)),M_HALF*(y(M,N-1)+y(M,N)),1)
        TangentSy(M,N)=-Tangent(0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N) &
                        +y(M,N-1)+y(M,N)),M_HALF*(x(M,N-1)+x(M,N)),M_HALF*(y(M,N-1)+y(M,N)),2)

        CurviNx(M,N)=CurviNx(M,N-1)
        CurviNy(M,N)=CurviNy(M,N-1)
        CurviSx(M,N)=-Tangent(x(M,N),y(M,N),x(M,N-1),y(M,N-1),1)
        CurviSy(M,N)=-Tangent(x(M,N),y(M,N),x(M,N-1),y(M,N-1),2)
        CurviEx(M,N)=CurviEx(M-1,N)
        CurviEy(M,N)=CurviEy(M-1,N)
        CurviWx(M,N)=-Tangent(x(M,N),y(M,N),x(M-1,N),y(M-1,N),1)
        CurviWy(M,N)=-Tangent(x(M,N),y(M,N),x(M-1,N),y(M-1,N),2)



      !South-East
        NormalS%x(M,1)=Normal(M_HALF*(x(M-1,1)+x(M,1)),M_HALF*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 1)
        NormalS%y(M,1)=Normal(M_HALF*(x(M-1,1)+x(M,1)),M_HALF*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 2)
        NormalE%x(M,1)=Normal(x(M,1), y(M,1), M_HALF*(x(M,1)+x(M,2)), M_HALF*(y(M,1)+y(M,2)), 1)
        NormalE%y(M,1)=Normal(x(M,1), y(M,1), M_HALF*(x(M,1)+x(M,2)), M_HALF*(y(M,1)+y(M,2)), 2)
        NormalW%x(M,1)=Normal(0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1) &
                    +y(M, 2)),M_HALF*(x(M-1, 1)+x(M, 1)),M_HALF*(y(M-1, 1)+y(M, 1)),1);
        NormalW%y(M,1)=Normal(0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1) &
                    +y(M, 2)),M_HALF*(x(M-1, 1)+x(M, 1)),M_HALF*(y(M-1, 1)+y(M, 1)),2);
        NormalN%x(M,1)=Normal(M_HALF*(x(M, 1)+x(M, 2)),M_HALF*(y(M, 1)+y(M, 2)),0.25d0*(x(M-1, 1)+x(M-1, 2) &
                    +x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),1);
        NormalN%y(M,1)=Normal(M_HALF*(x(M, 1)+x(M, 2)),M_HALF*(y(M, 1)+y(M, 2)),0.25d0*(x(M-1, 1)+x(M-1, 2) &
                    +x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),2);

        TangentNx(M,1)=-Tangent(M_HALF*(x(M,1)+x(M,2)),M_HALF*(y(M,1)+y(M,2)),0.25d0*(x(M,1)+x(M,2)+x(M-1,1) &
                    +x(M-1,2)),0.25d0*(y(M,1)+y(M,2)+y(M-1,1)+y(M-1,2)),1)
        TangentNy(M,1)=-Tangent(M_HALF*(x(M,1)+x(M,2)),M_HALF*(y(M,1)+y(M,2)),0.25d0*(x(M,1)+x(M,2)+x(M-1,1) &
                    +x(M-1,2)),0.25d0*(y(M,1)+y(M,2)+y(M-1,1)+y(M-1,2)),2)
        TangentWx(M,1)=-Tangent(0.25d0*(x(M-1,2)+x(M-1,1)+x(M,2)+x(M,1)),0.25d0*(y(M-1,2)+y(M-1,1)+y(M,2) &
                    +y(M,1)),M_HALF*(x(M-1,1)+x(M,1)), M_HALF*(y(M-1,1)+y(M,1)),1)
        TangentWy(M,1)=-Tangent(0.25d0*(x(M-1,2)+x(M-1,1)+x(M,2)+x(M,1)),0.25d0*(y(M-1,2)+y(M-1,1)+y(M,2) &
                    +y(M,1)),M_HALF*(x(M-1,1)+x(M,1)), M_HALF*(y(M-1,1)+y(M,1)),2)
        TangentEx(M,1)=-Tangent(x(M,1),y(M,1),M_HALF*(x(M,1)+x(M,2)),M_HALF*(y(M,1)+y(M,2)),1)
        TangentEy(M,1)=-Tangent(x(M,1),y(M,1),M_HALF*(x(M,1)+x(M,2)),M_HALF*(y(M,1)+y(M,2)),2)
        TangentSx(M,1)=-Tangent(M_HALF*(x(M-1,1)+x(M,1)), M_HALF*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 1)
        TangentSy(M,1)=-Tangent(M_HALF*(x(M-1,1)+x(M,1)), M_HALF*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 2)

         CurviNx(M,1)=-Tangent(x(M,1),y(M,1),x(M,2),y(M,2),1)
         CurviNy(M,1)=-Tangent(x(M,1),y(M,1),x(M,2),y(M,2),2)
         CurviSx(M,1)=CurviSx(M,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),1)
         CurviSy(M,1)=CurviSy(M,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),2)
         CurviEx(M,1)=CurviEx(M-1,1) !-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),1)
         CurviEy(M,1)=CurviEy(M-1,1) !-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),2)
         CurviWx(M,1)=-Tangent(x(M,1),y(M,1),x(M-1,1),y(M-1,1),1)
         CurviWy(M,1)=-Tangent(x(M,1),y(M,1),x(M-1,1),y(M-1,1),2)



      !South-West
        NormalS%x(1,1)=Normal(x(1,1),y(1,1),M_HALF*(x(2,1)+x(1,1)), M_HALF*(y(2,1)+y(1,1)), 1)
        NormalS%y(1,1)=Normal(x(1,1),y(1,1),M_HALF*(x(2,1)+x(1,1)), M_HALF*(y(2,1)+y(1,1)), 2)
        NormalW%x(1,1)=-Normal(x(1,1),y(1,1),M_HALF*(x(1,2)+x(1,1)), M_HALF*(y(1,2)+y(1,1)), 1)
        NormalW%y(1,1)=-Normal(x(1,1),y(1,1),M_HALF*(x(1,2)+x(1,1)), M_HALF*(y(1,2)+y(1,1)), 2)

        NormalE%x(1,1)=-Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),M_HALF*(x(1, 1)+x(2, 1)),M_HALF*(y(1, 1)+y(2, 1)),1);
        NormalE%y(1,1)=-Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),M_HALF*(x(1, 1)+x(2, 1)),M_HALF*(y(1, 1)+y(2, 1)),2);
        NormalN%x(1,1)=Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),M_HALF*(x(1, 1)+x(1, 2)),M_HALF*(y(1, 1)+y(1, 2)),1);
        NormalN%y(1,1)=Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),M_HALF*(x(1, 1)+x(1, 2)),M_HALF*(y(1, 1)+y(1, 2)),2);

        TangentNx(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),1)
        TangentNy(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),2)
        TangentEx(1,1)=Tangent(x(1,1),y(1,1),x(2,1),y(2,1),1)
        TangentEy(1,1)=Tangent(x(1,1),y(1,1),x(2,1),y(2,1),2)
        TangentWx(1,1)=Tangent(x(1,1),y(1,1),M_HALF*(x(1,1)+x(1,2)),M_HALF*(y(1,1)+y(1,2)),1)
        TangentWy(1,1)=Tangent(x(1,1),y(1,1),M_HALF*(x(1,1)+x(1,2)),M_HALF*(y(1,1)+y(1,2)),2)
        TangentSx(1,1)=-Tangent(x(1,1),y(1,1),M_HALF*(x(1,1)+x(2,1)),M_HALF*(y(1,1)+y(2,1)),1)
        TangentSy(1,1)=-Tangent(x(1,1),y(1,1),M_HALF*(x(1,1)+x(2,1)),M_HALF*(y(1,1)+y(2,1)),2)

        CurviNx(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),1)
        CurviNy(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),2)
        CurviSx(1,1)=CurviSx(1,2)
        CurviSy(1,1)=CurviSy(1,2)
        CurviEx(1,1)=-Tangent(x(1,1),y(1,1),x(2,1),y(2,1),1)
        CurviEy(1,1)=-Tangent(x(1,1),y(1,1),x(2,1),y(2,1),2)
        CurviWx(1,1)=CurviWx(2,1)
        CurviWy(1,1)=CurviWy(2,1)


      !North-West
        NormalN%x(1,N)=-Normal(x(1,N),y(1,N),M_HALF*(x(1,N)+x(2,N)),M_HALF*(y(1,N)+y(2,N)),1)
        NormalN%y(1,N)=-Normal(x(1,N),y(1,N),M_HALF*(x(1,N)+x(2,N)),M_HALF*(y(1,N)+y(2,N)),2)
        NormalW%x(1,N)=Normal(x(1,N),y(1,N),M_HALF*(x(1,N)+x(1,N-1)),M_HALF*(y(1,N)+y(1,N-1)),1)
        NormalW%y(1,N)=Normal(x(1,N),y(1,N),M_HALF*(x(1,N)+x(1,N-1)),M_HALF*(y(1,N)+y(1,N-1)),2)
        NormalE%x(1,N)=-Normal(M_HALF*(x(1, N)+x(2, N)), &
                            M_HALF*(y(1, N)+y(2, N)), &
                            0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)) ,1);
        NormalE%y(1,N)=-Normal(M_HALF*(x(1, N)+x(2, N)), &
                            M_HALF*(y(1, N)+y(2, N)), &
                            0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), 2)
        NormalS%x(1,N)=-Normal(0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                            M_HALF*(x(1, N-1)+x(1, N)), &
                            M_HALF*(y(1, N-1)+y(1, N)),1);
        NormalS%y(1,N)=-Normal(0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                            M_HALF*(x(1, N-1)+x(1, N)), &
                            M_HALF*(y(1, N-1)+y(1, N)), 2)
        TangentNx(1,N)=-Tangent(M_HALF*(x(2,N)+x(1,N)),M_HALF*(y(2,N)+y(1,N)),x(1,N),y(1,N),1)
        TangentNy(1,N)=-Tangent(M_HALF*(x(2,N)+x(1,N)),M_HALF*(y(2,N)+y(1,N)),x(1,N),y(1,N),2)
        TangentWx(1,N)=-Tangent(x(1,N),y(1,N),M_HALF*(x(1,N-1)+x(1,N)), M_HALF*(y(1,N-1)+y(1,N)), 1)
        TangentWy(1,N)=-Tangent(x(1,N),y(1,N),M_HALF*(x(1,N-1)+x(1,N)), M_HALF*(y(1,N-1)+y(1,N)), 2)
        TangentEx(1,N)=-Tangent(0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1)),0.25d0*(y(1,N)+y(1,N-1) &
                        +y(2,N)+y(2,N-1)),M_HALF*(x(2,N)+x(1,N)),M_HALF*(y(2,N)+y(1,N)),1)
        TangentEy(1,N)=-Tangent(0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1)),0.25d0*(y(1,N)+y(1,N-1) &
                        +y(2,N)+y(2,N-1)),M_HALF*(x(2,N)+x(1,N)),M_HALF*(y(2,N)+y(1,N)),2)
        TangentSx(1,N)=-Tangent(M_HALF*(x(1,N-1)+x(1,N)),M_HALF*(y(1,N-1)+y(1,N)),0.25d0*(x(2,N-1) &
                        +x(2,N)+x(1,N-1)+x(1,N)),0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)),1)
        TangentSy(1,N)=-Tangent(M_HALF*(x(1,N-1)+x(1,N)),M_HALF*(y(1,N-1)+y(1,N)),0.25d0*(x(2,N-1) &
                        +x(2,N)+x(1,N-1)+x(1,N)),0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)),2)

        CurviNx(1,N)=CurviNx(1,N-1)
        CurviNy(1,N)=CurviNy(1,N-1)
        CurviSx(1,N)=-Tangent(x(1,N),y(1,N),x(1,N-1),y(1,N-1),1)
        CurviSy(1,N)=-Tangent(x(1,N),y(1,N),x(1,N-1),y(1,N-1),2)
        CurviEx(1,N)=-Tangent(x(1,N),y(1,N),x(2,N),y(2,N),1)
        CurviEy(1,N)=-Tangent(x(1,N),y(1,N),x(2,N),y(2,N),2)
        CurviWx(1,N)=CurviWx(2,N)
        CurviWy(1,N)=CurviWy(2,N)

end subroutine compute_norm_tan_curv

subroutine compute_cellvol(M, N, x, y, CellVol, InvCellVol )
  use Maths_m
  implicit none

  integer, intent(in)    :: M, N
  real(8), intent(in)    :: x(1:M, 1:N), y(1:M, 1:N)                 ! needle position indexes
  real(8), intent(inout) :: CellVol(1:M,1:N), InvCellVol(1:M,1:N)

  integer :: i, j

  !$OMP PARALLEL DEFAULT(NONE) SHARED(x, y, N, M, CellVol, InvCellVol )
  !$OMP DO COLLAPSE(2)
  do j=2, N-1
    do i=2, M-1
      ! define the volume of elementary cell around a point everywhere but not on boundaries
!         CellVol(i,j)=0.25d0*(AreaElement(x(i-1,j-1),y(i-1,j-1),x(i+1,j-1),y(i+1,j-1),x(i+1,j+1),y(i+1,j+1),x(i-1,j+1),y(i-1,j+1)))
      CellVol(i,j)=0.125d0* AreaElement((x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)), & !x(i-1/2,j-1/2)
                                        (y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j)), &
                                        (x(i,j-1)+x(i+1,j-1)+x(i+1,j)+x(i,j)), &        !x(i+1/2,j-1/2)
                                        (y(i,j-1)+y(i+1,j-1)+y(i+1,j)+y(i,j)), &
                                        (x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), &        !x(i+1/2,j+1/2)
                                        (y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                                        (x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)), &        !x(i-1/2,j+1/2)
                                        (y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1)))

      InvCellVol(i,j) = M_ONE/CellVol(i,j)
    end do
  end do
  !$OMP END DO

  !$OMP DO
  do i=2,M-1
    ! NORTH
    CellVol(i,N)=0.25d0*AreaElement(x(i-1,N),y(i-1,N),x(i+1,N),y(i+1,N),x(i+1,N-1),y(i+1,N-1),x(i-1,N-1),y(i-1,N-1))
    ! SOUTH
    CellVol(i,1)=0.25d0*AreaElement(x(i+1,1),y(i+1,1),x(i+1,2),y(i+1,2),x(i-1,2),y(i-1,2),x(i-1,1),y(i-1,1))
  end do
  !$OMP END DO

  !$OMP DO
  do j=2,N-1
     ! WEST
     CellVol(1,j)=0.25d0*AreaElement(x(1,j-1),y(1,j-1),x(2,j-1),y(2,j-1),x(2,j+1),y(2,j+1),x(1,j+1),y(1,j+1))
     ! EAST
     CellVol(M,j)=0.25d0*AreaElement(x(M-1,j-1),y(M-1,j-1),x(M,j-1),y(M,j-1),x(M,j+1),y(M,j+1),x(M-1,j+1),y(M-1,j+1))
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !North-East
  CellVol(M,N)=AreaElement(0.25d0*(x(M-1, N-1)+x(M, N-1)+x(M, N)+x(M-1, N)), &
                           0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N)+y(M,N-1)), &
                           M_HALF*(x(M, N)+x(M, N-1)), &
                           M_HALF*(y(M, N)+y(M, N-1)), &
                           x(M,N), &
                           y(M,N), &
                           M_HALF*(x(M-1, N)+x(M, N)), &
                           M_HALF*(y(M-1, N)+y(M, N)))
!         CellVol(M,N)=0.25d0*AreaElement(x(M-1,N-1),y(M-1,N-1),x(M,N-1),y(M,N-1), x(M,N), y(M,N), x(M-1, N), y(M-1, N))
  !South-East
  CellVol(M,1)=AreaElement(M_HALF*(x(M-1, 1)+x(M, 1)),M_HALF*(y(M-1, 1)+y(M, 1)), &
                          x(M,1), &
                          y(M,1), &
                          M_HALF*(x(M, 1)+x(M, 2)), &
                          M_HALF*(y(M, 1)+y(M, 2)), &
                          0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)), &
                          0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)))
!         CellVol(M,1)=0.25d0*AreaElement(x(M-1,1),y(M-1,1),x(M,1),y(M,1),x(M,2),y(M,2),x(M-1,2),y(M-1,2))
  !South-West
  CellVol(1,1)=AreaElement(x(1,1), &
                           y(1,1), &
                           M_HALF*(x(1, 1)+x(2, 1)), &
                           M_HALF*(y(1, 1)+y(2, 1)), &
                           0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)), &
                           0.25d0*(y(1, 1)+y(1, 2)+y(2, 1)+y(2, 2)), &
                           M_HALF*(x(1, 1)+x(1, 2)), &
                           M_HALF*(y(1, 1)+y(1, 2)))
!         CellVol(1,1)=0.25d0*AreaElement(x(1,1),y(1,1),x(2,1),y(2,1),x(2,2),y(2,2),x(1,2),y(1,2))
!         CellVol(1,1)=M_TWO*CellVol(1,1)
  !North-West
  CellVol(1,N)=AreaElement(M_HALF*(x(1, N-1)+x(1, N)), &
                                M_HALF*(y(1, N-1)+y(1, N)), &
                                0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                                0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                                M_HALF*(x(1, N)+x(2, N)), &
                                M_HALF*(y(1, N)+y(2, N)), &
                                x(1,N), &
                                y(1,N))
!         CellVol(1,N)=0.25d0*AreaElement(x(1,N-1),y(1,N-1),x(2,N-1),y(2,N-1),x(2,N),y(2,N),x(1,N),y(1,N))
!         CellVol(1,N)=M_TWO*CellVol(1,N)
end subroutine compute_cellvol

subroutine deallocate_Norm(NormalN, NormalS, NormalE, NormalW)
  use Types_m
  implicit none

  type(VectorField), intent(inout) :: NormalN, NormalS, NormalW, NormalE ! normal to quadrangle elements

  deallocate(NormalN%x, NormalN%y, NormalN%N)
  deallocate(NormalS%x, NormalS%y, NormalS%N)
  deallocate(NormalE%x, NormalE%y, NormalE%N)
  deallocate(NormalW%x, NormalW%y, NormalW%N)
end subroutine deallocate_Norm
