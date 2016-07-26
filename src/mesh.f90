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
!> @file mesh.f90
!
! DESCRIPTION:
!> @brief A set a routines related to the meshs
!
!> @author
!> Nicolas Tancogne-Dejean
!> @date
!> 01 Jun 2016 - Initial Version
!> 15 Jun 2016 - Adding the check_divergences routine
!------------------------------------------------------------------------------


!> Allocate a MeshValues type
!> @param MeshValues Structure containing Mesh values
subroutine initmesh( mesh, M, N )
    use Types_m
    implicit none
    type(MeshValues), intent(INOUT) :: mesh
    integer, intent(IN) :: M, N

    allocate(mesh%Te(M,N))
    allocate(mesh%Th(M,N))
    allocate(mesh%Ts(M,N))
    allocate(mesh%Ne(M,N))
    allocate(mesh%Nh(M,N))

    mesh%M = M
    mesh%N = N
end subroutine initmesh

!> Deallocate a MeshValues type
!> @param MeshValues Structure containing Mesh values
subroutine releasemesh( mesh )
    use Types_m
    implicit none
    type(MeshValues), intent(INOUT) :: mesh

    deallocate(mesh%Te)
    deallocate(mesh%Th)
    deallocate(mesh%Ts)
    deallocate(mesh%Ne)
    deallocate(mesh%Nh)
end subroutine releasemesh

!> Copy data from oldmesh to newmesh
subroutine copy_mesh(newmesh, oldmesh)
  use Types_m
  implicit none
  type(MeshValues), intent(IN)  :: oldmesh
  type(MeshValues), intent(INOUT) :: newmesh

 !TODO: Print something to error.dat
  if(newmesh%M /= oldmesh%M .or. newmesh%N /= oldmesh%N ) &
    stop 'Invalid meshes used in copy_mesh'

  newmesh%Te(1:newmesh%M,1:newmesh%N) = oldmesh%Te(1:newmesh%M,1:newmesh%N)
  newmesh%Th(1:newmesh%M,1:newmesh%N) = oldmesh%Th(1:newmesh%M,1:newmesh%N)
  newmesh%Ts(1:newmesh%M,1:newmesh%N) = oldmesh%Ts(1:newmesh%M,1:newmesh%N)
  newmesh%Ne(1:newmesh%M,1:newmesh%N) = oldmesh%Ne(1:newmesh%M,1:newmesh%N)
  newmesh%Nh(1:newmesh%M,1:newmesh%N) = oldmesh%Nh(1:newmesh%M,1:newmesh%N)

end subroutine copy_mesh

! interpolation bilineaire ponderee par les aires
subroutine bilinear_interpol_dual(mesh, dual, InvCellVol )
   use Types_m
   implicit none
   type(MeshValues), intent(INOUT) :: mesh, dual
   real(8), intent(IN)             :: InvCellVol(mesh%M,mesh%N)

   integer :: i,j
   real(8) :: weight

    ! interpolation on dual mesh
    ! InterpolateBiCubic(phi_source, x_s, y_s, x_t, y_t, SizeXs, SizeYs, SizeXt, SizeYt, phi_target, Grad(phi)_targetX, Grad(phi)_targetY)

!     if(UseInterpolation.eq.1) then
!     !$OMP SECTIONS
!       !$OMP SECTION
!       call InterpolateBiCubic(Ne, x, y, xDual, yDual, M, N, M-1, N-1, NeDual, DummyDual, DummyDual)
!       !$OMP SECTION
!       call InterpolateBiCubic(Nh, x, y, xDual, yDual, M, N, M-1, N-1, NhDual, DummyDual, DummyDual)
!       !$OMP SECTION
!       call InterpolateBiCubic(Te, x, y, xDual, yDual, M, N, M-1, N-1, TeDual, DummyDual, DummyDual)
!       !$OMP SECTION
!       call InterpolateBiCubic(Th, x, y, xDual, yDual, M, N, M-1, N-1, ThDual, DummyDual, DummyDual)
!       !$OMP SECTION
!       call InterpolateBiCubic(Ts, x, y, xDual, yDual, M, N, M-1, N-1, TsDual, DummyDual, DummyDual)
!       !$OMP SECTION
!   !     call InterpolateBiCubic(intensity, x, y, xDual, yDual, M, N, M-1, N-1, intensityDual, DummyDual, DummyDual)
!   !     intensityDual=InterpolateSelner(intensity, x, y, xDual, yDual, M, N, M-1, N-1)
!
!   ! !       DEBUG: test de la fonction d'interpolation
!       call InterpolateBiCubic(intensity, x, y, xDual, yDual, M, N, M-1, N-1, intensityDual, DummyDual, DummyDual)
!   !     intensityDual=InterpolateSelner(intensity, x, y, xDual, yDual, M, N, M-1, N-1)
!
!     !$OMP END SECTIONS
!     end if

   !$OMP DO  COLLAPSE(2)
      do j=1, dual%N
        do i=1, dual%M

        weight = 1.0d0/ ( InvCellVol(i,j) + InvCellVol(i+1,j) + InvCellVol(i,j+1) + InvCellVol(i+1,j+1) )

        ! interpolation bilineaire ponderee par les aires
        dual%Ne(i,j) = ( mesh%Ne(i,j)* InvCellVol(i,j) + mesh%Ne(i+1,j)*InvCellVol(i+1,j) + mesh%Ne(i,j+1)*InvCellVol(i,j+1) + mesh%Ne(i+1,j+1) &
                *InvCellVol(i+1,j+1) ) * weight
        dual%Nh(i,j) = ( mesh%Nh(i,j)*InvCellVol(i,j) + mesh%Nh(i+1,j)*InvCellVol(i+1,j) + mesh%Nh(i,j+1)*InvCellVol(i,j+1) + mesh%Nh(i+1,j+1) &
                *InvCellVol(i+1,j+1) ) * weight
        dual%Te(i,j) = ( mesh%Te(i,j)*InvCellVol(i,j) + mesh%Te(i+1,j)*InvCellVol(i+1,j) + mesh%Te(i,j+1)*InvCellVol(i,j+1) + mesh%Te(i+1,j+1) &
                *InvCellVol(i+1,j+1) ) * weight
        dual%Th(i,j) = ( mesh%Th(i,j)*InvCellVol(i,j) + mesh%Th(i+1,j)*InvCellVol(i+1,j) + mesh%Th(i,j+1)*InvCellVol(i,j+1) + mesh%Th(i+1,j+1) &
                *InvCellVol(i+1,j+1) ) * weight
        dual%Ts(i,j) = ( mesh%Ts(i,j)*InvCellVol(i,j) + mesh%Ts(i+1,j)*InvCellVol(i+1,j) + mesh%Ts(i,j+1)*InvCellVol(i,j+1) + mesh%Ts(i+1,j+1) &
                *InvCellVol(i+1,j+1) ) * weight
        end do
      end do
      !$OMP END DO

 end subroutine

 subroutine check_divergences(mesh, maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN, i, j, x, y, t)
   use Types_m
   implicit none

   type(MeshValues), intent(IN) :: mesh
   integer, intent(in)          :: i, j
   real(8), intent(in)          :: x, y, t, maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN

   logical Diverged

   Diverged = .false.

    if(isnan(mesh%Te(i,j))) then
      write(95,*) "Divergence of Te at t=", t, "x(",i,j,")=", x, "y(",i,j,")=",y
      Diverged=.true.
    end if
    if(isnan(mesh%Th(i,j))) then
      write(95,*) "Divergence of Th at t=", t, "x(",i,j,")=", x, "y(",i,j,")=",y
      Diverged=.true.
    end if
    if(isnan(mesh%Ts(i,j))) then
      write(95,*) "Divergence of Ts at t=", t, "x(",i,j,")=", x, "y(",i,j,")=",y
      Diverged=.true.
    end if
    if(isnan(mesh%Ne(i,j))) then
      write(95,*) "Divergence of Ne at t=", t, "x(",i,j,")=", x, "y(",i,j,")=",y
      Diverged=.true.
    end if
    if(isnan(mesh%Nh(i,j))) then
      write(95,*) "Divergence of Nh at t=", t, "x(",i,j,")=", x, "y(",i,j,")=",y
      Diverged=.true.
    end if

    if(maxCFLxT.gt.1d0 .OR. maxCFLyT.gt.1d0) then
      write(95,*) "Bad convergence for Te,Th. t=", t, "(CFLx,CFLy)=", maxCFLxT, maxCFLyT
      Diverged=.true.
    end if
    if(maxCFLxN.gt.1d0 .OR. maxCFLyN.gt.1d0) then
      write(95,*) "Bad convergence for Ne,Nh. t=", t, "(CFLx,CFLy)=", maxCFLxN, maxCFLyN
      Diverged=.true.
    end if

    if(Diverged .eqv. .true.) then
          write(*,*) "Divergence detected. Please check error.dat for more information."
          stop
    end if

 end subroutine check_divergences
