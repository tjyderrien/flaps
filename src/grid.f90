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

subroutine compute_distances(M, N, x, y, DistN, DistS, DistE, DistW, DistDualN, DistDualS, DistDualE, DistDualW, CellAreaN, CellAreaS, CellAreaE, CellAreaW )
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
                  +(  0.25d0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))-0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)))**2) !0.25d0*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)-2d0*x(i+1,j+1)*x(i-1,j)-2d0*x(i+1,j+1)*x(i-1,j+1)+x(i+1,j)**2-2d0*x(i+1,j)*x(i-1,j)-2d0*x(i+1,j)*x(i-1,j+1)+x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)+x(i-1,j+1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)-2d0*y(i+1,j+1)*y(i-1,j)-2d0*y(i+1,j+1)*y(i-1,j+1)+y(i+1,j)**2-2d0*y(i+1,j)*y(i-1,j)-2d0*y(i+1,j)*y(i-1,j+1)+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)
        CellAreaS(i,j)=sqrt( (0.25d0*(x(i+1,j-1)+x(i,j-1)+x(i+1,j)+x(i,j))-0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i-1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i+1,j-1)+y(i,j-1)+y(i+1,j)+y(i,j))-0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i-1,j)+y(i,j)))**2 ) !0.25d0*(x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)-2d0*x(i+1,j)*x(i-1,j-1)-2d0*x(i+1,j)*x(i-1,j)+x(i+1,j-1)**2-2d0*x(i+1,j-1)*x(i-1,j-1)-2d0*x(i+1,j-1)*x(i-1,j)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)+x(i-1,j)**2+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)-2d0*y(i+1,j)*y(i-1,j-1)-2d0*y(i+1,j)*y(i-1,j)+y(i+1,j-1)**2-2d0*y(i+1,j-1)*y(i-1,j-1)-2d0*y(i+1,j-1)*y(i-1,j)+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)+y(i-1,j)**2)**(0.5d0)
        CellAreaE(i,j)=sqrt( (0.25d0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))-0.25d0*(x(i+1,j-1)+x(i,j-1)+x(i+1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))-0.25d0*(y(i+1,j-1)+y(i,j-1)+y(i+1,j)+y(i,j)))**2 ) !0.25d0*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)-2d0*x(i+1,j+1)*x(i+1,j-1)-2d0*x(i+1,j+1)*x(i,j-1)+x(i,j+1)**2-2d0*x(i,j+1)*x(i+1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)-2d0*y(i+1,j+1)*y(i+1,j-1)-2d0*y(i+1,j+1)*y(i,j-1)+y(i,j+1)**2-2d0*y(i,j+1)*y(i+1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)
        CellAreaW(i,j)=sqrt( (0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j))-0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i-1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j))-0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i-1,j)+y(i,j)))**2 ) !0.25d0*(x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)-2d0*x(i,j+1)*x(i-1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i-1,j-1)-2d0*x(i-1,j+1)*x(i,j-1)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)-2d0*y(i,j+1)*y(i-1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i-1,j-1)-2d0*y(i-1,j+1)*y(i,j-1)+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)

        DistN(i,j)=sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
        DistS(i,j)=sqrt((x(i,j)-x(i,j-1))**2+(y(i,j)-y(i,j-1))**2)
        DistE(i,j)=sqrt((x(i+1,j)-x(i,j))**2+(y(i+1,j)-y(i,j))**2)
        DistW(i,j)=sqrt((x(i,j)-x(i-1,j))**2+(y(i,j)-y(i-1,j))**2)

      end do
   end do

   do i=2,M-1
      !NORTH
         CellAreaN(i,N)=0d0 !0.5d0*(x(i+1,N)**2-2d0*x(i+1,N)*x(i-1,N)+x(i-1,N)**2+y(i+1,N)**2-2d0*y(i+1,N)*y(i-1,N)+y(i-1,N)**2)**0.5d0
         CellAreaS(i,N)=0.25d0*(x(i+1,N)**2+2d0*x(i+1,N)*x(i+1,N-1)-2d0*x(i+1,N)*x(i-1,N-1)-2d0*x(i+1,N)*x(i-1,N) &
                +x(i+1,N-1)**2-2d0*x(i+1,N-1)*x(i-1,N-1)-2d0*x(i+1,N-1)*x(i-1,N)+x(i-1,N-1)**2+2d0*x(i-1,N)*x(i-1,N-1) &
                +x(i-1,N)**2+y(i+1,N)**2+2d0*y(i+1,N)*y(i+1,N-1)-2d0*y(i+1,N)*y(i-1,N-1)-2d0*y(i+1,N)*y(i-1,N) &
                +y(i+1,N-1)**2-2d0*y(i+1,N-1)*y(i-1,N-1)-2d0*y(i+1,N-1)*y(i-1,N)+y(i-1,N-1)**2+2d0*y(i-1,N)*y(i-1,N-1)+y(i-1,N)**2)**0.5d0
         CellAreaW(i,N)=0.25d0*(x(i,N)**2+2d0*x(i,N)*x(i-1,N)-2d0*x(i,N)*x(i-1,N-1)-2d0*x(i,N)*x(i,N-1)+x(i-1,N)**2 &
                -2d0*x(i-1,N)*x(i-1,N-1)-2d0*x(i-1,N)*x(i,N-1)+x(i-1,N-1)**2+2d0*x(i-1,N-1)*x(i,N-1)+x(i,N-1)**2 &
                +y(i,N)**2+2d0*y(i,N)*y(i-1,N)-2d0*y(i,N)*y(i-1,N-1)-2d0*y(i,N)*y(i,N-1)+y(i-1,N)**2 &
                -2d0*y(i-1,N)*y(i-1,N-1)-2d0*y(i-1,N)*y(i,N-1)+y(i-1,N-1)**2+2d0*y(i-1,N-1)*y(i,N-1) &
                +y(i,N-1)**2)**0.5d0
         CellAreaE(i,N)=0.25d0*(x(i+1,N)**2+2d0*x(i+1,N)*x(i,N)-2d0*x(i+1,N)*x(i+1,N-1)-2d0*x(i+1,N)*x(i,N-1)+x(i,N)**2 &
                -2d0*x(i,N)*x(i+1,N-1)-2d0*x(i,N)*x(i,N-1)+x(i+1,N-1)**2+2d0*x(i+1,N-1)*x(i,N-1)+x(i,N-1)**2 &
                +y(i,N)**2+2d0*y(i,N)*y(i+1,N)-2d0*y(i,N)*y(i+1,N-1)-2d0*y(i,N)*y(i,N-1)+y(i+1,N)**2 &
                -2d0*y(i+1,N)*y(i+1,N-1)-2d0*y(i+1,N)*y(i,N-1)+y(i+1,N-1)**2+2d0*y(i+1,N-1)*y(i,N-1) &
                +y(i,N-1)**2)**0.5d0

         DistN(i,N)=0d0 !sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
         DistS(i,N)=sqrt((x(i,N)-x(i,N-1))**2+(y(i,N)-y(i,N-1))**2)
         DistE(i,N)=sqrt((x(i+1,N)-x(i,N))**2+(y(i+1,N)-y(i,N))**2)
         DistW(i,N)=sqrt((x(i,N)-x(i-1,N))**2+(y(i,N)-y(i-1,N))**2)

     !SOUTH
         CellAreaN(i,1)=0.25d0*(x(i+1,2)**2+2d0*x(i+1,2)*x(i+1,1)-2d0*x(i+1,2)*x(i-1,1)-2d0*x(i+1,2)*x(i-1,2) &
                 +x(i+1,1)**2-2d0*x(i+1,1)*x(i-1,1)-2d0*x(i+1,1)*x(i-1,2)+x(i-1,1)**2+2d0*x(i-1,1)*x(i-1,2)+x(i-1,2)**2 &
                 +y(i+1,2)**2+2d0*y(i+1,2)*y(i+1,1)-2d0*y(i+1,2)*y(i-1,1)-2d0*y(i+1,2)*y(i-1,2)+y(i+1,1)**2 &
                 -2d0*y(i+1,1)*y(i-1,1)-2d0*y(i+1,1)*y(i-1,2)+y(i-1,1)**2+2d0*y(i-1,1)*y(i-1,2)+y(i-1,2)**2)**(0.5d0)
         CellAreaS(i,1)=0d0 !sqrt((0.5d0*(x(i+1,1)+x(i,1))-0.5d0*(x(i-1,1)+x(i,1)))**2+(0.5d0*(y(i+1,1)+y(i,1))-0.5d0*(y(i-1,1)+y(i,1)))**2)
         CellAreaW(i,1)=0.25d0*(x(i,1)**2-2d0*x(i,1)*x(i,2)+2d0*x(i,1)*x(i-1,1)-2d0*x(i,1)*x(i-1,2)+x(i,2)**2 &
                  -2d0*x(i,2)*x(i-1,1)+2d0*x(i,2)*x(i-1,2)+x(i-1,1)**2-2d0*x(i-1,1)*x(i-1,2)+x(i-1,2)**2+y(i,1)**2-2d0*y(i,1)*y(i,2) &
                  +2d0*y(i,1)*y(i-1,1)-2d0*y(i,1)*y(i-1,2)+y(i,2)**2-2d0*y(i,2)*y(i-1,1)+2d0*y(i,2)*y(i-1,2)+y(i-1,1)**2 &
                  -2d0*y(i-1,1)*y(i-1,2)+y(i-1,2)**2)**(0.5d0)
         CellAreaE(i,1)=0.25d0*(x(i+1,2)**2-2d0*x(i+1,2)*x(i+1,1)-2d0*x(i+1,2)*x(i,1)+2d0*x(i+1,2)*x(i,2)+x(i+1,1)**2+2d0*x(i+1,1)*x(i,1) &
                  -2d0*x(i+1,1)*x(i,2)+x(i,1)**2-2d0*x(i,1)*x(i,2)+x(i,2)**2+y(i+1,2)**2-2d0*y(i+1,2)*y(i+1,1)-2d0*y(i+1,2)*y(i,1) &
                  +2d0*y(i+1,2)*y(i,2)+y(i+1,1)**2+2d0*y(i+1,1)*y(i,1)-2d0*y(i+1,1)*y(i,2)+y(i,1)**2-2d0*y(i,1)*y(i,2)+y(i,2)**2)**(0.5d0)

         DistN(i,1)=sqrt( (x(i,2)-x(i,1))**2 + (y(i,2)-y(i,1))**2 )
         DistS(i,1)=0d0
         DistE(i,1)=sqrt( (x(i+1,1)-x(i,1))**2 + (y(i+1,1)-y(i,1))**2 )
         DistW(i,1)=sqrt( (x(i,1)-x(i-1,1))**2 + (y(i,1)-y(i-1,1))**2 )
   end do

   do j=2,N-1
     !WEST
           CellAreaN(1,j)=0.25d0*(x(2,j+1)**2+2d0*x(2,j+1)*x(2,j)-2d0*x(2,j+1)*x(1,j)-2d0*x(2,j+1)*x(1,j+1)+x(2,j)**2 &
                  -2d0*x(2,j)*x(1,j)-2d0*x(2,j)*x(1,j+1)+x(1,j)**2+2d0*x(1,j)*x(1,j+1)+x(1,j+1)**2+y(2,j+1)**2 &
                  +2d0*y(2,j+1)*y(2,j)-2d0*y(2,j+1)*y(1,j)-2d0*y(2,j+1)*y(1,j+1)+y(2,j)**2-2d0*y(2,j)*y(1,j) &
                  -2d0*y(2,j)*y(1,j+1)+y(1,j)**2+2d0*y(1,j)*y(1,j+1)+y(1,j+1)**2)**0.5d0
          CellAreaS(1,j)=0.25d0*(x(2,j)**2+2d0*x(2,j)*x(2,j-1)-2d0*x(2,j)*x(1,j)-2d0*x(2,j)*x(1,j-1)+x(2,j-1)**2 &
                  -2d0*x(2,j-1)*x(1,j)-2d0*x(2,j-1)*x(1,j-1)+x(1,j)**2+2d0*x(1,j)*x(1,j-1)+x(1,j-1)**2+y(2,j)**2 &
                  +2d0*y(2,j)*y(2,j-1)-2d0*y(2,j)*y(1,j)-2d0*y(2,j)*y(1,j-1)+y(2,j-1)**2-2d0*y(2,j-1)*y(1,j) &
                  -2d0*y(2,j-1)*y(1,j-1)+y(1,j)**2+2d0*y(1,j)*y(1,j-1)+y(1,j-1)**2)**0.5d0
          CellAreaW(1,j)=sqrt( (0.5d0*(x(1,j+1)+x(1,j))-0.5d0*(x(1,j-1)+x(1,j)))**2 + (0.5d0*(y(1,j+1)+y(1,j)) &
                        -0.5d0*(y(1,j-1)+y(1,j)))**2 )
          CellAreaE(1,j)=0.25d0*(x(2,j+1)**2+2d0*x(2,j+1)*x(1,j+1)-2d0*x(2,j+1)*x(2,j-1)-2d0*x(2,j+1)*x(1,j-1) &
                  +x(1,j+1)**2-2d0*x(1,j+1)*x(2,j-1)-2d0*x(1,j+1)*x(1,j-1)+x(2,j-1)**2+2d0*x(2,j-1)*x(1,j-1) &
                  +x(1,j-1)**2+y(2,j+1)**2+2d0*y(2,j+1)*y(1,j+1)-2d0*y(2,j+1)*y(2,j-1)-2d0*y(2,j+1)*y(1,j-1) &
                  +y(1,j+1)**2-2d0*y(1,j+1)*y(2,j-1)-2d0*y(1,j+1)*y(1,j-1)+y(2,j-1)**2+2d0*y(2,j-1)*y(1,j-1) &
                  +y(1,j-1)**2)**0.5d0

          DistN(1,j)=sqrt( (x(1,j+1)-x(1,j))**2 + (y(1,j+1)-y(1,j))**2 )
          DistS(1,j)=sqrt( (x(1,j)-x(1,j-1))**2 + (y(1,j)-y(1,j-1))**2 )
          DistE(1,j)=sqrt( (x(2,j)-x(1,j))**2 + (y(2,j)-y(1,j))**2 )
          DistW(1,j)=0d0

     !EAST
          CellAreaN(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M,j)-2d0*x(M,j+1)*x(M-1,j)-2d0*x(M,j+1)*x(M-1,j+1)+x(M,j)**2 &
                  -2d0*x(M,j)*x(M-1,j)-2d0*x(M,j)*x(M-1,j+1)+x(M-1,j)**2+2d0*x(M-1,j)*x(M-1,j+1)+x(M-1,j+1)**2 &
                  +y(M,j+1)**2+2d0*y(M,j+1)*y(M,j)-2d0*y(M,j+1)*y(M-1,j)-2d0*y(M,j+1)*y(M-1,j+1)+y(M,j)**2 &
                  -2d0*y(M,j)*y(M-1,j)-2d0*y(M,j)*y(M-1,j+1)+y(M-1,j)**2+2d0*y(M-1,j)*y(M-1,j+1)+y(M-1,j+1)**2)**0.5d0
          CellAreaS(M,j)=0.25d0*(x(M,j)**2+2d0*x(M,j)*x(M,j-1)-2d0*x(M,j)*x(M-1,j-1)-2d0*x(M,j)*x(M-1,j)+x(M,j-1)**2 &
                  -2d0*x(M-1,j-1)*x(M,j-1)-2d0*x(M,j-1)*x(M-1,j)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M-1,j)+x(M-1,j)**2 &
                  +y(M,j)**2+2d0*y(M,j)*y(M,j-1)-2d0*y(M,j)*y(M-1,j-1)-2d0*y(M,j)*y(M-1,j)+y(M,j-1)**2 &
                  -2d0*y(M-1,j-1)*y(M,j-1)-2d0*y(M,j-1)*y(M-1,j)+y(M-1,j-1)**2+2d0*y(M-1,j-1)*y(M-1,j) &
                  +y(M-1,j)**2)**0.5d0
          CellAreaW(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M-1,j+1)-2d0*x(M,j+1)*x(M-1,j-1)-2d0*x(M,j+1)*x(M,j-1) &
                  +x(M-1,j+1)**2-2d0*x(M-1,j+1)*x(M-1,j-1)-2d0*x(M-1,j+1)*x(M,j-1)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M,j-1) &
                  +x(M,j-1)**2+y(M,j+1)**2+2d0*y(M,j+1)*y(M-1,j+1)-2d0*y(M,j+1)*y(M-1,j-1)-2d0*y(M,j+1)*y(M,j-1) &
                  +y(M-1,j+1)**2-2d0*y(M-1,j+1)*y(M-1,j-1)-2d0*y(M-1,j+1)*y(M,j-1)+y(M-1,j-1)**2 &
                  +2d0*y(M-1,j-1)*y(M,j-1)+y(M,j-1)**2)**0.5d0
          CellAreaE(M,j)=sqrt( (0.5d0*(x(M,j+1)+x(M,j))-0.5d0*(x(M,j-1)+x(M,j)) )**2 + (0.5d0*(y(M,j+1)+y(M,j)) &
                        -0.5d0*(y(M,j-1)+y(M,j)) )**2)

          DistN(M,j)=sqrt( (x(M,j+1)-x(M,j))**2 + (y(M,j+1)-y(M,j))**2 )
          DistS(M,j)=sqrt( (x(M,j-1)-x(M,j))**2 + (y(M,j-1)-y(M,j))**2 )
          DistE(M,j)=0d0
          DistW(M,j)=sqrt( (x(M-1,j)-x(M,j))**2 + (y(M-1,j)-y(M,j))**2 )

   end do

   !North-east
        CellAreaE(M,N)=sqrt((x(M,N)-0.5d0*(x(M,N)+x(M,N-1)))**2+(y(M,N)-0.5d0*(y(M,N)+y(M,N-1)))**2)
        CellAreaN(M,N)=sqrt((x(M,N)-0.5d0*(x(M-1,N)+x(M,N)))**2+(y(M,N)-0.5d0*(y(M-1,N)+y(M,N)))**2)
        CellAreaS(M,N)=0.25d0*(x(M,N)**2+2d0*x(M,N)*x(M,N-1)-2d0*x(M,N)*x(M-1,N-1)-2d0*x(M-1,N)*x(M,N)+x(M,N-1)**2 &
                -2d0*x(M-1,N-1)*x(M,N-1)-2d0*x(M-1,N)*x(M,N-1)+x(M-1,N-1)**2+2d0*x(M-1,N)*x(M-1,N-1)+x(M-1,N)**2 &
                +y(M,N)**2+2d0*y(M,N)*y(M,N-1)-2d0*y(M,N)*y(M-1,N-1)-2d0*y(M-1,N)*y(M,N)+y(M,N-1)**2 &
                -2d0*y(M-1,N-1)*y(M,N-1)-2d0*y(M-1,N)*y(M,N-1)+y(M-1,N-1)**2+2d0*y(M-1,N)*y(M-1,N-1) &
                +y(M-1,N)**2)**(0.5d0)
        CellAreaW(M,N)=0.25d0*(x(M-1,N)**2+2d0*x(M-1,N)*x(M,N)-2d0*x(M-1,N)*x(M-1,N-1)-2d0*x(M-1,N)*x(M,N-1)+x(M,N)**2 &
                -2d0*x(M,N)*x(M-1,N-1)-2d0*x(M,N)*x(M,N-1)+x(M-1,N-1)**2+2d0*x(M-1,N-1)*x(M,N-1)+x(M,N-1)**2 &
                +y(M-1,N)**2+2d0*y(M-1,N)*y(M,N)-2d0*y(M-1,N)*y(M-1,N-1)-2d0*y(M-1,N)*y(M,N-1)+y(M,N)**2 &
                -2d0*y(M,N)*y(M-1,N-1)-2d0*y(M,N)*y(M,N-1)+y(M-1,N-1)**2+2d0*y(M-1,N-1)*y(M,N-1) &
                +y(M,N-1)**2)**(0.5d0)

        DistN(M,N)=0d0
        DistS(M,N)=sqrt((x(M,N)-x(M,N-1))**2+(y(M,N)-y(M,N-1))**2)
        DistE(M,N)=0d0
        DistW(M,N)=sqrt((x(M,N)-x(M-1,N))**2+(y(M,N)-y(M-1,N))**2)

   !South-east
        CellAreaS(M,1)=sqrt(( x(M,1)-0.5d0*(x(M,1)+x(M-1,1)) )**2+( y(M,1)-0.5d0*(y(M,1)+y(M-1,1)) )**2)
        CellAreaE(M,1)=sqrt( ( x(M,1) - 0.5d0*(x(M,1)+x(M,2)) )**2 + ( y(M,1) - 0.5d0*(y(M,1)+y(M,2)) )**2 )
        CellAreaN(M,1)=0.25d0*(x(M,1)**2+2d0*x(M,1)*x(M,2)-2d0*x(M-1,1)*x(M,1)-2d0*x(M-1,2)*x(M,1)+x(M,2)**2 &
                -2d0*x(M-1,1)*x(M,2)-2d0*x(M-1,2)*x(M,2)+x(M-1,1)**2+2d0*x(M-1,1)*x(M-1,2)+x(M-1,2)**2 &
                +y(M,1)**2+2d0*y(M,1)*y(M,2)-2d0*y(M-1,1)*y(M,1)-2d0*y(M-1,2)*y(M,1)+y(M,2)**2 &
                -2d0*y(M-1,1)*y(M,2)-2d0*y(M-1,2)*y(M,2)+y(M-1,1)**2+2d0*y(M-1,1)*y(M-1,2) &
                +y(M-1,2)**2)**(0.5d0)
        CellAreaW(M,1)=0.25d0*(x(M-1,1)**2-2d0*x(M-1,1)*x(M-1,2)+2d0*x(M-1,1)*x(M,1)-2d0*x(M-1,1)*x(M,2) &
                +x(M-1,2)**2-2d0*x(M-1,2)*x(M,1)+2d0*x(M-1,2)*x(M,2)+x(M,1)**2-2d0*x(M,1)*x(M,2)+x(M,2)**2 &
                +y(M-1,1)**2-2d0*y(M-1,1)*y(M-1,2)+2d0*y(M-1,1)*y(M,1)-2d0*y(M-1,1)*y(M,2)+y(M-1,2)**2 &
                -2d0*y(M-1,2)*y(M,1)+2d0*y(M-1,2)*y(M,2)+y(M,1)**2-2d0*y(M,1)*y(M,2)+y(M,2)**2)**(0.5d0)

         DistN(M,1)=sqrt( (x(M,2)-x(M,1))**2 + (y(M,2)-y(M,1))**2 )
         DistS(M,1)=0d0
         DistE(M,1)=0d0 !sqrt( (x(i+1,1)-x(i,1))**2 + (y(i+1,1)-y(i,1))**2 )
         DistW(M,1)=sqrt( (x(M,1)-x(M-1,1))**2 + (y(M,1)-y(M-1,1))**2 )


  !South-west

        CellAreaS(1,1)=sqrt((x(1,1)-0.5d0*(x(2,1)+x(1,1)))**2+(y(1,1)-0.5d0*(y(2,1)+y(1,1)))**2)
        CellAreaW(1,1)=sqrt((x(1,1)-0.5d0*(x(1,2)+x(1,1)))**2+(y(1,1)-0.5d0*(y(1,2)+y(1,1)))**2)
        CellAreaN(1,1)=sqrt((0.25d0*(x(2,2)+x(2,1)+x(1,2)+x(1,1))-0.5d0*(x(1,1)+x(1,2)))**2 &
                +(0.25d0*(y(2,2)+y(2,1)+y(1,2)+y(1,1))-0.5d0*(y(1,2)+y(1,1)))**2) !sqrt((0.5d0*(x(1,2)+x(1,1))-0.25d0*(x(2,2)+x(2,1)+x(1,1)+x(1,2)))**2d0+(0.5d0*(y(1,2)+y(1,1))-0.25d0*(y(2,2)+y(2,1)+y(1,1)+y(1,2)))**2d0) !0.25d0*(x(1,1)**2+2d0*x(1,1)*x(1,2)-2d0*x(1,1)*x(2,1)-2d0*x(1,1)*x(2,2)+x(1,2)**2-2d0*x(1,2)*x(2,1)-2d0*x(1,2)*x(2,2)+x(2,1)**2+2d0*x(2,1)*x(2,2)+x(2,2)**2+y(1,1)**2+2d0*y(1,1)*y(1,2)-2d0*y(1,1)*y(2,1)-2d0*y(1,1)*y(2,2)+y(1,2)**2-2d0*y(1,2)*y(2,1)-2d0*y(1,2)*y(2,2)+y(2,1)**2+2d0*y(2,1)*y(2,2)+y(2,2)**2)**(0.5d0);
        CellAreaE(1,1)=sqrt((0.25d0*(x(2,2)+x(2,1)+x(1,2)+x(1,1))-0.5d0*(x(2,1)+x(1,1)))**2 &
                +(0.25d0*(y(2,2)+y(2,1)+y(1,2)+y(1,1))-0.5d0*(y(2,1)+y(1,1)))**2) !0.25d0*(x(1,1)**2-2d0*x(1,1)*x(1,2)+2d0*x(1,1)*x(2,1)-2d0*x(1,1)*x(2,2)+x(1,2)**2-2d0*x(1,2)*x(2,1)+2d0*x(1,2)*x(2,2)+x(2,1)**2-2d0*x(2,1)*x(2,2)+x(2,2)**2+y(1,1)**2-2d0*y(1,1)*y(1,2)+2d0*y(1,1)*y(2,1)-2d0*y(1,1)*y(2,2)+y(1,2)**2-2d0*y(1,2)*y(2,1)+2d0*y(1,2)*y(2,2)+y(2,1)**2-2d0*y(2,1)*y(2,2)+y(2,2)**2)**(0.5d0)

        DistN(1,1)=sqrt((x(1,2)-x(1,1))**2+(y(1,2)-y(1,1))**2)
        DistS(1,1)=0d0
        DistE(1,1)=sqrt((x(2,1)-x(1,1))**2+(y(2,1)-y(1,1))**2)
        DistW(1,1)=0d0

  !North-West
        CellAreaN(1,N)=sqrt((x(1,N)-0.5d0*(x(1,N)+x(2,N)))**2+(y(1,N)-0.5d0*(y(1,N)+y(2,N)))**2)
        CellAreaW(1,N)=sqrt((x(1,N)-0.5d0*(x(1,N)+x(1,N-1)))**2+((y(1,N)-0.5d0*(y(1,N)+y(1,N-1))))**2)
        CellAreaS(1,N)=sqrt((0.25d0*(x(2,N-1)+x(1,N-1)+x(2,N)+x(1,N))-0.5d0*(x(1,N-1)+x(1,N)))**2 &
                      +(0.25d0*(y(2,N-1)+y(1,N-1)+y(2,N)+y(1,N))-0.5d0*(y(1,N-1)+y(1,N)))**2) !sqrt((0.5d0*(x(1,N-1)+x(1,N))-0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)))**2+(0.5d0*(y(1,N-1)+y(1,N))-0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)))**2)
        CellAreaE(1,N)=sqrt((0.5d0*(x(2,N)+x(1,N))-0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)))**2 &
                      +(0.5d0*(y(2,N)+y(1,N))-0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)))**2) !sqrt((0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1))-(0.5d0*(x(2,N)+x(1,N))))**2+(0.25d0*(y(1,N)+y(1,N-1)+y(2,N)+y(2,N-1))-0.5d0*(y(2,N)+y(1,N)))**2)

         DistN(1,N)=0d0 !sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
         DistS(1,N)=sqrt((x(1,N)-x(1,N-1))**2+(y(1,N)-y(1,N-1))**2)
         DistE(1,N)=sqrt((x(2,N)-x(1,N))**2+(y(2,N)-y(1,N))**2)
         DistW(1,N)=0d0 !sqrt((x(i,N)-x(i-1,N))**2+(y(i,N)-y(i-1,N))**2)


  ! distance for cross diffusion, equal to cell border "areas" in 2D
  DistDualN(:,:)=CellAreaN(:,:)
  DistDualS(:,:)=CellAreaS(:,:)
  DistDualE(:,:)=CellAreaE(:,:)
  DistDualW(:,:)=CellAreaW(:,:)


end subroutine compute_distances

subroutine compute_norm_tan_curv(M, N, x, y, NormalNx, NormalNy, NormalSx, NormalSy, NormalEx, NormalEy, NormalWx, NormalWy, TangentNx, TangentNy, TangentSx, TangentSy, &
      TangentEx, TangentEy, TangentWx, TangentWy, CurviNx, CurviNy, CurviSx, CurviSy, CurviEx, CurviEy, CurviWx, CurviWy )
  implicit none

  integer, intent(in)    :: M, N
  real(8), intent(in)    :: x(1:M, 1:N), y(1:M, 1:N)                 ! needle position indexes
  real(8), intent(inout) ::  NormalWx(1:M,1:N), NormalWy(1:M,1:N),  &
                             NormalEx(1:M,1:N), NormalEy(1:M,1:N),  &                ! normal to quadrangle elements
                             NormalNx(1:M,1:N), NormalNy(1:M,1:N),  &
                             NormalSx(1:M,1:N), NormalSy(1:M,1:N),  &
                             TangentWx(1:M,1:N), TangentWy(1:M,1:N), &                ! Tangent to quadrangle elements
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

        NormalNx(i,j)=Normal(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                     0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)),1) !4d0*(0.25d0*y(i-1,j)+0.25d0*y(i-1,j+1)-0.25d0*y(i+1,j+1)-0.25d0*y(i+1,j))/(x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)-2d0*x(i-1,j)*x(i+1,j+1)-2d0*x(i-1,j)*x(i+1,j)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i+1,j+1)-2d0*x(i-1,j+1)*x(i+1,j)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)+x(i+1,j)**2+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)-2d0*y(i-1,j)*y(i+1,j+1)-2d0*y(i-1,j)*y(i+1,j)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i+1,j+1)-2d0*y(i-1,j+1)*y(i+1,j)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)+y(i+1,j)**2)**(0.5d0)
        NormalNy(i,j)=Normal(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                     0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)),2) !NormalNy(i,j)=-4d0*(0.25d0*x(i-1,j)+0.25d0*x(i-1,j+1)-0.25d0*x(i+1,j+1)-0.25d0*x(i+1,j))/(x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)-2d0*x(i-1,j)*x(i+1,j+1)-2d0*x(i-1,j)*x(i+1,j)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i+1,j+1)-2d0*x(i-1,j+1)*x(i+1,j)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)+x(i+1,j)**2+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)-2d0*y(i-1,j)*y(i+1,j+1)-2d0*y(i-1,j)*y(i+1,j)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i+1,j+1)-2d0*y(i-1,j+1)*y(i+1,j)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)+y(i+1,j)**2)**(0.5d0)
        NormalSx(i,j)=Normal(0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)),1) ! -4d0*(0.25d0*y(i-1,j-1)+0.25d0*y(i-1,j)-0.25d0*y(i+1,j)-0.25d0*y(i+1,j-1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)-2d0*x(i-1,j-1)*x(i+1,j)-2d0*x(i-1,j-1)*x(i+1,j-1)+x(i-1,j)**2-2d0*x(i-1,j)*x(i+1,j)-2d0*x(i-1,j)*x(i+1,j-1)+x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)+x(i+1,j-1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)-2d0*y(i-1,j-1)*y(i+1,j)-2d0*y(i-1,j-1)*y(i+1,j-1)+y(i-1,j)**2-2d0*y(i-1,j)*y(i+1,j)-2d0*y(i-1,j)*y(i+1,j-1)+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)+y(i+1,j-1)**2)**(0.5d0)
        NormalSy(i,j)=Normal(0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)),2) ! 4d0*(0.25d0*x(i-1,j-1)+0.25d0*x(i-1,j)-0.25d0*x(i+1,j)-0.25d0*x(i+1,j-1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)-2d0*x(i-1,j-1)*x(i+1,j)-2d0*x(i-1,j-1)*x(i+1,j-1)+x(i-1,j)**2-2d0*x(i-1,j)*x(i+1,j)-2d0*x(i-1,j)*x(i+1,j-1)+x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)+x(i+1,j-1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)-2d0*y(i-1,j-1)*y(i+1,j)-2d0*y(i-1,j-1)*y(i+1,j-1)+y(i-1,j)**2-2d0*y(i-1,j)*y(i+1,j)-2d0*y(i-1,j)*y(i+1,j-1)+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)+y(i+1,j-1)**2)**(0.5d0)
        NormalEx(i,j)=Normal(0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j+1)+x(i+1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i+1,j+1)+y(i+1,j)+y(i,j+1)+y(i,j)), 1) ! -4d0*(0.25d0*y(i+1,j-1)+0.25d0*y(i,j-1)-0.25d0*y(i+1,j+1)-0.25d0*y(i,j+1))/(x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)-2d0*x(i+1,j-1)*x(i+1,j+1)-2d0*x(i+1,j-1)*x(i,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i+1,j+1)-2d0*x(i,j-1)*x(i,j+1)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)+x(i,j+1)**2+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)-2d0*y(i+1,j-1)*y(i+1,j+1)-2d0*y(i+1,j-1)*y(i,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i+1,j+1)-2d0*y(i,j-1)*y(i,j+1)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)+y(i,j+1)**2)**(0.5d0)
        NormalEy(i,j)=Normal(0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j+1)+x(i+1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i+1,j+1)+y(i+1,j)+y(i,j+1)+y(i,j)), 2) ! 4d0*(0.25d0*x(i+1,j-1)+0.25d0*x(i,j-1)-0.25d0*x(i+1,j+1)-0.25d0*x(i,j+1))/(x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)-2d0*x(i+1,j-1)*x(i+1,j+1)-2d0*x(i+1,j-1)*x(i,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i+1,j+1)-2d0*x(i,j-1)*x(i,j+1)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)+x(i,j+1)**2+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)-2d0*y(i+1,j-1)*y(i+1,j+1)-2d0*y(i+1,j-1)*y(i,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i+1,j+1)-2d0*y(i,j-1)*y(i,j+1)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)+y(i,j+1)**2)**(0.5d0)
        NormalWx(i,j)=Normal(0.25d0*(x(i-1,j+1)+x(i-1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i-1,j)+y(i,j+1)+y(i,j)), &
                     0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 1) ! 4d0*(0.25d0*y(i-1,j-1)+0.25d0*y(i,j-1)-0.25d0*y(i,j+1)-0.25d0*y(i-1,j+1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)-2d0*x(i-1,j-1)*x(i,j+1)-2d0*x(i-1,j-1)*x(i-1,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i,j+1)-2d0*x(i,j-1)*x(i-1,j+1)+x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)+x(i-1,j+1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)-2d0*y(i-1,j-1)*y(i,j+1)-2d0*y(i-1,j-1)*y(i-1,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i,j+1)-2d0*y(i,j-1)*y(i-1,j+1)+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)
        NormalWy(i,j)=Normal(0.25d0*(x(i-1,j+1)+x(i-1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i-1,j)+y(i,j+1)+y(i,j)), &
                     0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 2)! -4d0*(0.25d0*x(i-1,j-1)+0.25d0*x(i,j-1)-0.25d0*x(i,j+1)-0.25d0*x(i-1,j+1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)-2d0*x(i-1,j-1)*x(i,j+1)-2d0*x(i-1,j-1)*x(i-1,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i,j+1)-2d0*x(i,j-1)*x(i-1,j+1)+x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)+x(i-1,j+1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)-2d0*y(i-1,j-1)*y(i,j+1)-2d0*y(i-1,j-1)*y(i-1,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i,j+1)-2d0*y(i,j-1)*y(i-1,j+1)+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)


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
         NormalEx(i,N)=-Normal(0.5d0*(x(i,N)+x(i+1,N)), 0.5d0*(y(i,N)+y(i+1,N)), 0.25d0*(x(i+1,N)+x(i,N)+x(i+1,N-1)+x(i,N-1)), &
                0.25d0*(y(i+1,N)+y(i,N)+y(i+1,N-1)+y(i,N-1)), 1)
         NormalEy(i,N)=-Normal(0.5d0*(x(i,N)+x(i+1,N)), 0.5d0*(y(i,N)+y(i+1,N)), 0.25d0*(x(i+1,N)+x(i,N)+x(i+1,N-1)+x(i,N-1)), &
                0.25d0*(y(i+1,N)+y(i,N)+y(i+1,N-1)+y(i,N-1)), 2)
         NormalWx(i,N)=Normal(0.5d0*(x(i-1,N)+x(i,N)), 0.5d0*(y(i-1,N)+y(i,N)), 0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), &
                0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),1);
         NormalWy(i,N)=Normal(0.5d0*(x(i-1,N)+x(i,N)), 0.5d0*(y(i-1,N)+y(i,N)), 0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), &
                0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),2);
         NormalNx(i,N)=Normal(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),1)
         NormalNy(i,N)=Normal(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),2)
         NormalSx(i,N)=-Normal(0.25d0*(x(i+1,N-1)+x(i+1,N)+x(i,N)+x(i,N-1)), 0.25d0*(y(i+1,N-1)+y(i+1,N)+y(i,N)+y(i,N-1)), &
                0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)),0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),1)
         NormalSy(i,N)=-Normal(0.25d0*(x(i+1,N-1)+x(i+1,N)+x(i,N)+x(i,N-1)), 0.25d0*(y(i+1,N-1)+y(i+1,N)+y(i,N)+y(i,N-1)), &
                0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)),0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),2)

         TangentWx(i,N)=-Tangent(0.5d0*(x(i-1,N)+x(i,N)),0.5d0*(y(i-1,N)+y(i,N)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N)+x(i,N-1)), &
                  0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N)+y(i,N-1)),1)
         TangentWy(i,N)=-Tangent(0.5d0*(x(i-1,N)+x(i,N)),0.5d0*(y(i-1,N)+y(i,N)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N)+x(i,N-1)), &
                  0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N)+y(i,N-1)),2)
         TangentEx(i,N)=-Tangent(0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)), &
                  0.5d0*(x(i+1,N)+x(i,N)), 0.5d0*(y(i+1,N)+y(i,N)), 1)
         TangentEy(i,N)=-Tangent(0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)), &
                  0.5d0*(x(i+1,N)+x(i,N)), 0.5d0*(y(i+1,N)+y(i,N)), 2)
         TangentNx(i,N)=-Tangent(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),1)
         TangentNy(i,N)=-Tangent(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),2)
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
         NormalEx(i,1)=-Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), &
                0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),1)
         NormalEy(i,1)=-Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), &
                0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),2)
         NormalWx(i,1)=Normal(0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)), &
                0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),1)
         NormalWy(i,1)=Normal(0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)), &
                0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),2)
         NormalNx(i,1)=Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), &
                0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),1)
         NormalNy(i,1)=Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), &
                0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),2)
         NormalSx(i,1)=Normal(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),1)
         NormalSy(i,1)=Normal(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),2)

         TangentNx(i,1)=Tangent(0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)), &
                  0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),1)
         TangentNy(i,1)=Tangent(0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)), &
                  0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),2)
         TangentSx(i,1)=-Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),1)
         TangentSy(i,1)=-Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),2)
         TangentEx(i,1)=-Tangent(0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)), &
                  0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),1)
         TangentEy(i,1)=-Tangent(0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)), &
                  0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),2)
         TangentWx(i,1)=Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)), &
                  0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),1)
         TangentWy(i,1)=Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)), &
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
          NormalEx(1,j)=-Normal(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)), &
                        0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),1)
          NormalEy(1,j)=-Normal(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)), &
                        0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),2)
          NormalNx(1,j)=Normal(0.25d0*(x(2,j+1)+x(1,j+1)+x(2,j)+x(1,j)),0.25d0*(y(2,j+1)+y(1,j+1)+y(2,j)+y(1,j)), &
                      0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),1)
          NormalNy(1,j)=Normal(0.25d0*(x(2,j+1)+x(1,j+1)+x(2,j)+x(1,j)),0.25d0*(y(2,j+1)+y(1,j+1)+y(2,j)+y(1,j)), &
                      0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),2)
          NormalSx(1,j)=-Normal(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)), &
                        0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),1)
          NormalSy(1,j)=-Normal(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)), &
                        0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),2)
          NormalWx(1,j)=Normal(0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),0.5d0*(x(1,j-1)+x(1,j)), &
                        0.5d0*(y(1,j-1)+y(1,j)),1)
          NormalWy(1,j)=Normal(0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),0.5d0*(x(1,j-1)+x(1,j)), &
                        0.5d0*(y(1,j-1)+y(1,j)),2)

!           write(*,*) "[Debug]", NormalWx(1,j)**2+NormalWy(1,j)**2

          TangentNx(1,j)=-Tangent(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) &
                  +y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),1)
          TangentNy(1,j)=-Tangent(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) &
                  +y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),2)
          TangentSx(1,j)=-Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.25d0*(x(1,j-1)+x(1,j) &
                  +x(2,j-1)+x(2,j)),0.25d0*(y(1,j-1)+y(1,j)+y(2,j-1)+y(2,j)),1)
          TangentSy(1,j)=-Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.25d0*(x(1,j-1)+x(1,j) &
                  +x(2,j-1)+x(2,j)),0.25d0*(y(1,j-1)+y(1,j)+y(2,j-1)+y(2,j)),2)
          TangentWx(1,j)=Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)), &
                  0.5d0*(y(1,j+1)+y(1,j)),1)
          TangentWy(1,j)=Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)), &
                  0.5d0*(y(1,j+1)+y(1,j)),2)
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
          TangentNx(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) &
                          +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1)
          TangentNy(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) &
                          +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2)
          TangentEx(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)), &
                          0.5d0*(y(M,j+1)+y(M,j)),1)
          TangentEy(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)), &
                          0.5d0*(y(M,j+1)+y(M,j)),2)
          TangentSx(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) &
                          +y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),1)
          TangentSy(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) &
                          +y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),2)
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

          NormalEx(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)), &
                        0.5d0*(y(M,j-1)+y(M,j)),1);
          NormalEy(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)), &
                        0.5d0*(y(M,j-1)+y(M,j)),2);
          NormalWx(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) &
                       +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)), &
                      0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1)
          NormalWy(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) &
          +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) &
          +y(M,j-1)+y(M,j)),2)
          NormalNx(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1) &
                      +x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1)
          NormalNy(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) &
                      +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2)
          NormalSx(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1) &
                      +x(M,j-1)+x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),1)
          NormalSy(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1) &
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
        NormalEx(M,N)=-Normal(x(M,N),y(M,N),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),1)
        NormalEy(M,N)=-Normal(x(M,N),y(M,N),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),2)
        NormalNx(M,N)=Normal(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),1)
        NormalNy(M,N)=Normal(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),2)
        NormalWx(M,N)=Normal(0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),1)
        NormalWy(M,N)=Normal(0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),2)
        NormalSx(M,N)=-Normal(0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),1)
        NormalSy(M,N)=-Normal(0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),2)

        TangentNx(M,N)=-Tangent(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),1)
        TangentNy(M,N)=-Tangent(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),2)
        TangentEx(M,N)=-Tangent(0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),x(M,N),y(M,N),1)
        TangentEy(M,N)=-Tangent(0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),x(M,N),y(M,N),2)
        TangentWx(M,N)=-Tangent(0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),0.25d0*(x(M-1,N-1) &
                        +x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),1)
        TangentWy(M,N)=-Tangent(0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),0.25d0*(x(M-1,N-1)+x(M-1,N) &
                        +x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),2)
        TangentSx(M,N)=-Tangent(0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N) &
                        +y(M,N-1)+y(M,N)),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),1)
        TangentSy(M,N)=-Tangent(0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N) &
                        +y(M,N-1)+y(M,N)),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),2)

        CurviNx(M,N)=CurviNx(M,N-1)
        CurviNy(M,N)=CurviNy(M,N-1)
        CurviSx(M,N)=-Tangent(x(M,N),y(M,N),x(M,N-1),y(M,N-1),1)
        CurviSy(M,N)=-Tangent(x(M,N),y(M,N),x(M,N-1),y(M,N-1),2)
        CurviEx(M,N)=CurviEx(M-1,N)
        CurviEy(M,N)=CurviEy(M-1,N)
        CurviWx(M,N)=-Tangent(x(M,N),y(M,N),x(M-1,N),y(M-1,N),1)
        CurviWy(M,N)=-Tangent(x(M,N),y(M,N),x(M-1,N),y(M-1,N),2)



      !South-East
        NormalSx(M,1)=Normal(0.5d0*(x(M-1,1)+x(M,1)),0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 1)
        NormalSy(M,1)=Normal(0.5d0*(x(M-1,1)+x(M,1)),0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 2)
        NormalEx(M,1)=Normal(x(M,1), y(M,1), 0.5d0*(x(M,1)+x(M,2)), 0.5d0*(y(M,1)+y(M,2)), 1)
        NormalEy(M,1)=Normal(x(M,1), y(M,1), 0.5d0*(x(M,1)+x(M,2)), 0.5d0*(y(M,1)+y(M,2)), 2)
        NormalWx(M,1)=Normal(0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1) &
                    +y(M, 2)),0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)),1);
        NormalWy(M,1)=Normal(0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1) &
                    +y(M, 2)),0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)),2);
        NormalNx(M,1)=Normal(0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)),0.25d0*(x(M-1, 1)+x(M-1, 2) &
                    +x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),1);
        NormalNy(M,1)=Normal(0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)),0.25d0*(x(M-1, 1)+x(M-1, 2) &
                    +x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),2);

        TangentNx(M,1)=-Tangent(0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),0.25d0*(x(M,1)+x(M,2)+x(M-1,1) &
                    +x(M-1,2)),0.25d0*(y(M,1)+y(M,2)+y(M-1,1)+y(M-1,2)),1)
        TangentNy(M,1)=-Tangent(0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),0.25d0*(x(M,1)+x(M,2)+x(M-1,1) &
                    +x(M-1,2)),0.25d0*(y(M,1)+y(M,2)+y(M-1,1)+y(M-1,2)),2)
        TangentWx(M,1)=-Tangent(0.25d0*(x(M-1,2)+x(M-1,1)+x(M,2)+x(M,1)),0.25d0*(y(M-1,2)+y(M-1,1)+y(M,2) &
                    +y(M,1)),0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)),1)
        TangentWy(M,1)=-Tangent(0.25d0*(x(M-1,2)+x(M-1,1)+x(M,2)+x(M,1)),0.25d0*(y(M-1,2)+y(M-1,1)+y(M,2) &
                    +y(M,1)),0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)),2)
        TangentEx(M,1)=-Tangent(x(M,1),y(M,1),0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),1)
        TangentEy(M,1)=-Tangent(x(M,1),y(M,1),0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),2)
        TangentSx(M,1)=-Tangent(0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 1)
        TangentSy(M,1)=-Tangent(0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 2)

         CurviNx(M,1)=-Tangent(x(M,1),y(M,1),x(M,2),y(M,2),1)
         CurviNy(M,1)=-Tangent(x(M,1),y(M,1),x(M,2),y(M,2),2)
         CurviSx(M,1)=CurviSx(M,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),1)
         CurviSy(M,1)=CurviSy(M,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),2)
         CurviEx(M,1)=CurviEx(M-1,1) !-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),1)
         CurviEy(M,1)=CurviEy(M-1,1) !-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),2)
         CurviWx(M,1)=-Tangent(x(M,1),y(M,1),x(M-1,1),y(M-1,1),1)
         CurviWy(M,1)=-Tangent(x(M,1),y(M,1),x(M-1,1),y(M-1,1),2)



      !South-West
        NormalSx(1,1)=Normal(x(1,1),y(1,1),0.5d0*(x(2,1)+x(1,1)), 0.5d0*(y(2,1)+y(1,1)), 1)
        NormalSy(1,1)=Normal(x(1,1),y(1,1),0.5d0*(x(2,1)+x(1,1)), 0.5d0*(y(2,1)+y(1,1)), 2)
        NormalWx(1,1)=-Normal(x(1,1),y(1,1),0.5d0*(x(1,2)+x(1,1)), 0.5d0*(y(1,2)+y(1,1)), 1)
        NormalWy(1,1)=-Normal(x(1,1),y(1,1),0.5d0*(x(1,2)+x(1,1)), 0.5d0*(y(1,2)+y(1,1)), 2)

        NormalEx(1,1)=-Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)),1);
        NormalEy(1,1)=-Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)),2);
        NormalNx(1,1)=Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)),1);
        NormalNy(1,1)=Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)),2);

        TangentNx(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),1)
        TangentNy(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),2)
        TangentEx(1,1)=Tangent(x(1,1),y(1,1),x(2,1),y(2,1),1)
        TangentEy(1,1)=Tangent(x(1,1),y(1,1),x(2,1),y(2,1),2)
        TangentWx(1,1)=Tangent(x(1,1),y(1,1),0.5d0*(x(1,1)+x(1,2)),0.5d0*(y(1,1)+y(1,2)),1)
        TangentWy(1,1)=Tangent(x(1,1),y(1,1),0.5d0*(x(1,1)+x(1,2)),0.5d0*(y(1,1)+y(1,2)),2)
        TangentSx(1,1)=-Tangent(x(1,1),y(1,1),0.5d0*(x(1,1)+x(2,1)),0.5d0*(y(1,1)+y(2,1)),1)
        TangentSy(1,1)=-Tangent(x(1,1),y(1,1),0.5d0*(x(1,1)+x(2,1)),0.5d0*(y(1,1)+y(2,1)),2)

        CurviNx(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),1)
        CurviNy(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),2)
        CurviSx(1,1)=CurviSx(1,2)
        CurviSy(1,1)=CurviSy(1,2)
        CurviEx(1,1)=-Tangent(x(1,1),y(1,1),x(2,1),y(2,1),1)
        CurviEy(1,1)=-Tangent(x(1,1),y(1,1),x(2,1),y(2,1),2)
        CurviWx(1,1)=CurviWx(2,1)
        CurviWy(1,1)=CurviWy(2,1)


      !North-West
        NormalNx(1,N)=-Normal(x(1,N),y(1,N),0.5d0*(x(1,N)+x(2,N)),0.5d0*(y(1,N)+y(2,N)),1)
        NormalNy(1,N)=-Normal(x(1,N),y(1,N),0.5d0*(x(1,N)+x(2,N)),0.5d0*(y(1,N)+y(2,N)),2)
        NormalWx(1,N)=Normal(x(1,N),y(1,N),0.5d0*(x(1,N)+x(1,N-1)),0.5d0*(y(1,N)+y(1,N-1)),1)
        NormalWy(1,N)=Normal(x(1,N),y(1,N),0.5d0*(x(1,N)+x(1,N-1)),0.5d0*(y(1,N)+y(1,N-1)),2)
        NormalEx(1,N)=-Normal(0.5d0*(x(1, N)+x(2, N)), &
                            0.5d0*(y(1, N)+y(2, N)), &
                            0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)) ,1);
        NormalEy(1,N)=-Normal(0.5d0*(x(1, N)+x(2, N)), &
                            0.5d0*(y(1, N)+y(2, N)), &
                            0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), 2)
        NormalSx(1,N)=-Normal(0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                            0.5d0*(x(1, N-1)+x(1, N)), &
                            0.5d0*(y(1, N-1)+y(1, N)),1);
        NormalSy(1,N)=-Normal(0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                            0.5d0*(x(1, N-1)+x(1, N)), &
                            0.5d0*(y(1, N-1)+y(1, N)), 2)
        TangentNx(1,N)=-Tangent(0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),x(1,N),y(1,N),1)
        TangentNy(1,N)=-Tangent(0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),x(1,N),y(1,N),2)
        TangentWx(1,N)=-Tangent(x(1,N),y(1,N),0.5d0*(x(1,N-1)+x(1,N)), 0.5d0*(y(1,N-1)+y(1,N)), 1)
        TangentWy(1,N)=-Tangent(x(1,N),y(1,N),0.5d0*(x(1,N-1)+x(1,N)), 0.5d0*(y(1,N-1)+y(1,N)), 2)
        TangentEx(1,N)=-Tangent(0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1)),0.25d0*(y(1,N)+y(1,N-1) &
                        +y(2,N)+y(2,N-1)),0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),1)
        TangentEy(1,N)=-Tangent(0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1)),0.25d0*(y(1,N)+y(1,N-1) &
                        +y(2,N)+y(2,N-1)),0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),2)
        TangentSx(1,N)=-Tangent(0.5d0*(x(1,N-1)+x(1,N)),0.5d0*(y(1,N-1)+y(1,N)),0.25d0*(x(2,N-1) &
                        +x(2,N)+x(1,N-1)+x(1,N)),0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)),1)
        TangentSy(1,N)=-Tangent(0.5d0*(x(1,N-1)+x(1,N)),0.5d0*(y(1,N-1)+y(1,N)),0.25d0*(x(2,N-1) &
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
  implicit none

  integer, intent(in)    :: M, N
  real(8), intent(in)    :: x(1:M, 1:N), y(1:M, 1:N)                 ! needle position indexes
  real(8), intent(inout) :: CellVol(1:M,1:N), InvCellVol(1:M,1:N)

  real(8) :: AreaElement
  integer :: i, j

  do j=2, N-1
    do i=2, M-1
      ! define the volume of elementary cell around a point everywhere but not on boundaries
!         CellVol(i,j)=0.25d0*(AreaElement(x(i-1,j-1),y(i-1,j-1),x(i+1,j-1),y(i+1,j-1),x(i+1,j+1),y(i+1,j+1),x(i-1,j+1),y(i-1,j+1)))
      CellVol(i,j)=AreaElement(0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)), & !x(i-1/2,j-1/2)
                               0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j)), &
                               0.25d0*(x(i,j-1)+x(i+1,j-1)+x(i+1,j)+x(i,j)), &        !x(i+1/2,j-1/2)
                               0.25d0*(y(i,j-1)+y(i+1,j-1)+y(i+1,j)+y(i,j)), &
                               0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), &        !x(i+1/2,j+1/2)
                               0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                               0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)), &        !x(i-1/2,j+1/2)
                               0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1)))

      InvCellVol(i,j) = 1.0d0/CellVol(i,j)
    end do
  end do

  do i=2,M-1
    ! NORTH
    CellVol(i,N)=0.25d0*AreaElement(x(i-1,N),y(i-1,N),x(i+1,N),y(i+1,N),x(i+1,N-1),y(i+1,N-1),x(i-1,N-1),y(i-1,N-1))
    ! SOUTH
    CellVol(i,1)=0.25d0*AreaElement(x(i+1,1),y(i+1,1),x(i+1,2),y(i+1,2),x(i-1,2),y(i-1,2),x(i-1,1),y(i-1,1))
  end do

  do j=2,N-1
     ! WEST
     CellVol(1,j)=0.25d0*AreaElement(x(1,j-1),y(1,j-1),x(2,j-1),y(2,j-1),x(2,j+1),y(2,j+1),x(1,j+1),y(1,j+1))
     ! EAST
     CellVol(M,j)=0.25d0*AreaElement(x(M-1,j-1),y(M-1,j-1),x(M,j-1),y(M,j-1),x(M,j+1),y(M,j+1),x(M-1,j+1),y(M-1,j+1))
  end do

  !North-East
  CellVol(M,N)=AreaElement(0.25d0*(x(M-1, N-1)+x(M, N-1)+x(M, N)+x(M-1, N)),0.25d0*(y(M-1,N-1)+y(M-1,N) &
                      +y(M,N)+y(M,N-1)), 0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)), x(M,N),y(M,N), &
                      0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)))
!         CellVol(M,N)=0.25d0*AreaElement(x(M-1,N-1),y(M-1,N-1),x(M,N-1),y(M,N-1), x(M,N), y(M,N), x(M-1, N), y(M-1, N))
  !South-East
  CellVol(M,1)=AreaElement(0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)), x(M,1),y(M,1), &
                      0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)), 0.25d0*(x(M-1, 1)+x(M-1, 2) &
                      +x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)))
!         CellVol(M,1)=0.25d0*AreaElement(x(M-1,1),y(M-1,1),x(M,1),y(M,1),x(M,2),y(M,2),x(M-1,2),y(M-1,2))
  !South-West
  CellVol(1,1)=AreaElement(x(1,1),y(1,1),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)), &
                      0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2)+y(2, 1) &
                      +y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)))
!         CellVol(1,1)=0.25d0*AreaElement(x(1,1),y(1,1),x(2,1),y(2,1),x(2,2),y(2,2),x(1,2),y(1,2))
!         CellVol(1,1)=2d0*CellVol(1,1)
  !North-West
  CellVol(1,N)=AreaElement(0.5d0*(x(1, N-1)+x(1, N)), &
                                0.5d0*(y(1, N-1)+y(1, N)), &
                                0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                                0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                                0.5d0*(x(1, N)+x(2, N)), &
                                0.5d0*(y(1, N)+y(2, N)), &
                                x(1,N), &
                                y(1,N))
!         CellVol(1,N)=0.25d0*AreaElement(x(1,N-1),y(1,N-1),x(2,N-1),y(2,N-1),x(2,N),y(2,N),x(1,N),y(1,N))
!         CellVol(1,N)=2d0*CellVol(1,N)
end subroutine compute_cellvol
