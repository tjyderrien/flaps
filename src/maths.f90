!------------------------------------------------------------------------------
!> @file maths.f90
!
! DESCRIPTION:
!> @brief A set of mat functions
!
!> @author
!> Thibault J.Y. Derrien
!> @date
!> 01 Jun 2016 - Initial Version
!------------------------------------------------------------------------------






!    function InterpolateSelner(Z, xS, yS, xT, yT, dimXs, dimYs, dimXt, dimYt)
!    ! from Z expressed on a source mesh (xS,yS),
!    ! returns Z expressed on a target mesh (xT, yT) by bilinear interapolation for convex irregular mesh
!    ! [Selner & Wastermann, J. Comp. Phys 79, 1-11 (1988)]
!
!      integer(8) dimXs, dimYs, &                                !dimension of the source mesh
!                dimXt, dimYt                                !dimensions of the target mesh
!
!      real(8) :: InterpolateSelner(1:dimXt,1:dimYt)                 ! new value of Z
!      real(8) Z(1:dimXs,1:dimYs), &                                !source data
!              xS(1:dimXs,1:dimYs), yS(1:dimXs,1:dimYs), &
!              xT(1:dimXt, 1:dimYt), yT(1:dimXt, 1:dimYt)
!      real(8) alpha1, alpha2, AireElementVoisin, AireSum
!      real(8) CoordMatrix(1:2,1:2),&                        ! small matrix for base change
!              CoordsT(1,1:2), NewCoordsT(1,1:2),&                ! coordinate vectors of target point
!              CoordsS(1,1:2), NewCoordsS(1,1:2),&                ! coordinate vectors of source points
!              p, q                                        ! simple coefficients
!
!      integer(8) iS, jS, iT, jT                                 ! source indexes, target indexes
!      integer(8) Neighbours(1:4,1:2) !contain the index of the nearest points in source mesh
!
!      !$O M P PARALLEL DO
!      do iT=1, dimXt
!        !$O M P DO
!        do jT=1, dimYt
!          ! catch the neighbours of point (i,j) contained in source mesh
!          Neighbours=FindNeighbours(xT(iT,jT), yT(iT,jT), xS, yS, numberOfNeighbours, dimXs, dimYs)
!          if(Neighbours(1,1).eq.0) then !point is outside of the found element
!!             write(*,*) AireElementVoisin, "!=", AireSum
!            InterpolateSelner(iT,jT)=0d0
!          else !interpolate for real
!!             write(*,*) AireElementVoisin, "==", AireSum
!
!            ! the coordinate of target point is transformed into squared element
!            CoordsT(1,1)=xT(iT,jT); CoordsT(1,2)=yT(iT,jT)
!!             CoordMatrix(1,1)=Neighbours(4,1); CoordMatrix(1,2)=Neighbours(4,2);
!!             CoordMatrix(2,1)=Neighbours(2,1); CoordMatrix(2,2)=Neighbours(2,2)
!! il semble qu  il y a une erreur dans les deux lignes precedentes ! test
!            CoordMatrix(1,1)=Neighbours(2,1); CoordMatrix(1,2)=Neighbours(4,1)
!            CoordMatrix(2,1)=Neighbours(2,2); CoordMatrix(2,2)=Neighbours(4,2)
!            CoordMatrix=LaInv(CoordMatrix)
!
!            NewCoordsT=transpose(matmul(CoordMatrix, transpose(CoordsT)))
!            CoordsS(1,1)=Neighbours(3,1); CoordsS(1,2)=Neighbours(3,2)
!            NewCoordsS=transpose(matmul(CoordMatrix, transpose(CoordsS)))
!
!            if((NewCoordsS(1,1)-1d0)<=1d-15) then
!              alpha2=NewCoordsT(1,2)/(1d0+(NewCoordsT(1,1)*(NewCoordsS(1,2)-1d0)))
!            else
!              p=0.5d0*(1d0+NewCoordsT(1,1)*(NewCoordsS(1,2)-1d0)-NewCoordsT(1,2)*(NewCoordsS(1,1)-1d0))
!              q=NewCoordsT(1,2)*(NewCoordsS(1,1)-1d0)
!              alpha2=(-p+(p**2+q)**0.5d0)/(NewCoordsS(1,1)-1d0)
!            end if
!
!            alpha1=NewCoordsT(1,1)/(1d0+alpha2*(NewCoordsS(1,1)-1d0))
!
!            ! use the 4 neighbours to interpolate
!            InterpolateSelner(iT,jT)=(1d0-alpha1)*(1d0-alpha2)*Z(Neighbours(1,1),Neighbours(1,2)) &
!                                + alpha1*(1d0-alpha2)*Z(Neighbours(2,1),Neighbours(2,2)) &
!                                + alpha1*alpha2*Z(Neighbours(3,1),Neighbours(3,2)) &
!                                + (1d0-alpha1)*alpha2*Z(Neighbours(4,1),Neighbours(4,2))
!          end if
!        end do
!        !$O M P END DO
!      end do
!      !$O M P END PARALLEL DO
!      ! $ OMP BARRIER
!    end function InterpolateSelner


!     subroutine InterpolateBiCubic(Z, xS, yS, xT, yT, dimXs, dimYs, dimXt, dimYt, Interpolated, DerivativeX, DerivativeY)
!     ! use bicubic interpolation
!
!     integer(8) dimXs, dimYs, &                                !dimension of the source mesh
!                 dimXt, dimYt                                !dimensions of the target mesh
!
! !       real(8) :: InterpolateBiCubic(1:dimXt,1:dimYt)                 ! new value of Z
!       real(8) Z(1:dimXs,1:dimYs), &                                !source data
!               xS(1:dimXs,1:dimYs), yS(1:dimXs,1:dimYs), &
!               xT(1:dimXt, 1:dimYt), yT(1:dimXt, 1:dimYt), &
!               xSs(1:dimXs*dimYs), ySs(1:dimXs*dimYs), & !same in 1D matrixes
!               xTt(1:dimXt*dimYt), yTt(1:dimXt*dimYt), &
!               ShepardSource(1:dimXs*dimYs), &
!               ShepardX(1:dimXs*dimYs), ShepardY(1:dimXs*dimYs), &
!               ShepardWeight(1:dimXs*dimYs), ShepardRadius(1:dimXs*dimYs), &
!               ShepardA(1:9,dimXs*dimYs), &
!               ShepardValues(1:dimXt, 1:dimYt), &
!               ShepardNext(1:dimXs*dimYs), &
!               ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, CSTVAL, &
!               Interpolated(1:dimXt, 1:dimYt), &
!               DerivativeX(1:dimXt, 1:dimYt), &
!               DerivativeY(1:dimXt, 1:dimYt)
!
!       integer(8) iS, jS, iT, jT, &
!                   ShepardCells(1:dimXt, 1:dimYt)
!
!       integer(4) ShepardError
!
!       ! 1/ transfer everything in vectors
!       !$O M P PARALLEL DO
!       do iS=1,dimXs
!         !$O M P DO
!         do jS=1,dimYs
!           ShepardX((iS-1)*dimXs+jS)=xS(iS,jS)
!           ShepardY((iS-1)*dimXs+jS)=yS(iS,jS)
!           ShepardSource((iS-1)*dimXs+jS)=Z(iS,jS)
!         end do
!         !$O M P END DO
!       end do
!       !$O M P END PARALLEL DO
!
!       ! 2/ compute Cshep2
!       CALL Cshep2(dimXs*dimYs, ShepardX, ShepardY, ShepardSource, 11, 15, dimXt, &
!              ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, &
!              ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA, ShepardError)
!
!       ! 3/ compute the interpolation
!       if(ShepardError.eq.0) then
!         !$O M P PARALLEL DO
!         do i=1,dimXt
!           !$O M P DO
!           do j=1,dimYt
! !             Interpolated(i,j)=CSTVAL(xT(i,j), yT(i,j), dimXs*dimYs, ShepardX, ShepardY, ShepardSource, dimXt, ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA)
!             call CS2GRD(xT(i,j), yT(i,j), dimXs*dimYs, ShepardX, ShepardY, ShepardSource, dimXt, ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA, Interpolated(i,j), DerivativeX(i,j), DerivativeY(i,j), ShepardError)
!           end do
!           !$O M P END DO
!         end do
!         !$O M P END PARALLEL DO
!
!         ! 4/ extract the result
! !         InterpolateBiCubic(:,:)=Interpolated(:,:)
!       else
!         write(*,*) 'BicubicInterp. error code = ', ShepardError
!       end if
!
!     end subroutine InterpolateBiCubic


!
!     function AreaElement(x1,y1,x2,y2,x3,y3,x4,y4)
!     ! works with convex elements!
!       real(8) x1, y1, x2, y2, x3,y3,x4,y4, AreaElement
! !       AreaElement=0.5d0*abs((x3-x1)*(y4-y2)-(y3-y1)*(x4-x2))
!       AreaElement=0.5d0*abs((x3-x1)*(y2-y4)-(y3-y1)*(x2-x4))
!     end function AreaElement
!
!     function AreaTri(xA, yA, xB, yB, xP, yP)
!       real(8) xA, yA, xB, yB, xP, yP, AreaTri
!       AreaTri=0.5d0*abs((xA-xP)*(yB-yP)-(xB-xP)*(yA-yP))
!     end function AreaTri
!
!     function Distance(x1, y1, x2, y2)
!       real(8) Distance, x1, y1, x2, y2
!       Distance=sqrt((x2-x1)**2+(y2-y1)**2)
!     end function Distance
!
!     function Normal(x1,y1,x2,y2,axis)
!     !compute the normal vector to the line defined by given distances in cartesian plane
!       real(8) x1, y1, x2, y2
!       real(8) NormalX, NormalY, Dx, Dy, Normal
!       real(8) normVector
!       integer(8) axis
!       Dx=x1-x2; Dy=y1-y2
!       normVector=sqrt(Dx**2+Dy**2)
!       NormalX=-Dy/normVector
!       NormalY=Dx/normVector
!       if(axis.eq.1) then
!         Normal=NormalX
!       else
!         Normal=NormalY
!       end if
!     end function Normal
!
!     function Tangent(x1,y1,x2,y2,axis)
!     !compute the tangent vector to the line defined by given distances in cartesian plane
!       real(8) x1, y1, x2, y2
!       real(8) TangentX, TangentY, Dx, Dy, Tangent
!       real(8) normVector
!       integer(8) axis
!       Dx=x1-x2; Dy=y1-y2
!       normVector=sqrt(Dx**2+Dy**2)
!       TangentX=Dx/normVector
!       TangentY=Dy/normVector
!       if(axis.eq.1) then
!         Tangent=TangentX
!       else
!         Tangent=TangentY
!       end if
!     end function Tangent


!   write(*,*) "[TEST] Area function"
!   write(*,*) AreaElement(5.99d-7,-2.5d-7, 7.68d-7,-2.85d-7, 6.11d-7,-2.11d-7,2.71d-7, -1.66d-7)
!   write(*,*) "[TEST] Triangle area function"
!   write(*,*) "5.99d-7,-2.5d-7 ; 7.68d-7,-2.85d-7 ; 6.11d-7,-2.11d-7 ; 2.71d-7, -1.66d-7"
!   write(*,*) "P 3.26d-7, 3.82d-8"
!   write(*,*) AreaTri(5.99d-7,-2.5d-7,7.68d-7,-2.85d-7,5.74d-7,-2.07d-7)
!   write(*,*) AreaTri(5.74d-7,-2.07d-7,7.68d-7,-2.85d-7,6.11d-7,-2.11d-7)
!   write(*,*) AreaTri(5.74d-7,-2.07d-7,6.11d-7,-2.11d-7, 2.71d-7, -1.66d-7)
!   write(*,*) AreaTri(5.74d-7,-2.07d-7,2.71d-7, -1.66d-7,5.99d-7,-2.5d-7)
!   write(*,*) "Sum=", AreaTri(5.99d-7,-2.5d-7,7.68d-7,-2.85d-7,5.74d-7,-2.07d-7) + &
!   AreaTri(5.74d-7,-2.07d-7,7.68d-7,-2.85d-7,6.11d-7,-2.11d-7) + &
!   AreaTri(5.74d-7,-2.07d-7,6.11d-7,-2.11d-7, 2.71d-7, -1.66d-7) + &
!   AreaTri(5.74d-7,-2.07d-7,2.71d-7, -1.66d-7,5.99d-7,-2.5d-7)
