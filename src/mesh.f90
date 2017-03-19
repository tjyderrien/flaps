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
!> @file mesh.f90
!
! DESCRIPTION:
!> @brief A set a routines related to the meshs
!------------------------------------------------------------------------------

module Mesh_m

  implicit none

  private

  integer, public, parameter :: NORTH = 1
  integer, public, parameter :: SOUTH = 2
  integer, public, parameter :: EAST = 3
  integer, public, parameter :: WEST = 4

end module Mesh_m

!> Allocate a MeshValues type
!> @param MeshValues Structure containing Mesh values
subroutine initmesh( mesh, M, N )
    use Mesh_m
    use Types_m
    implicit none
    type(MeshValues), intent(INOUT) :: mesh
    integer, intent(IN) :: M, N

    integer :: i,j

    allocate(mesh%Te(M*N))
    allocate(mesh%Th(M*N))
    allocate(mesh%Ts(M*N))
    allocate(mesh%Ne(M*N))
    allocate(mesh%Nh(M*N))

    allocate(mesh%Map(4,M*N))

    do j = 1, N
      do i = 1,M
        mesh%Map(NORTH,i+(j-1)*M) =  i+j*M
        mesh%Map(SOUTH,i+(j-1)*M) =  i+(j-2)*M
        mesh%Map(EAST,i+(j-1)*M) =  i+1+(j-1)*M
        mesh%Map(WEST,i+(j-1)*M) =  i-1+(j-1)*M
      end do
    end do

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
    deallocate(mesh%Map)
end subroutine releasemesh

!> Copy data from oldmesh to newmesh
subroutine copy_mesh(newmesh, oldmesh)
  use Types_m
  implicit none
  type(MeshValues), intent(IN)  :: oldmesh
  type(MeshValues), intent(INOUT) :: newmesh

  integer :: i

 !TODO: Print something to error.dat
  if(newmesh%M /= oldmesh%M .or. newmesh%N /= oldmesh%N ) &
    stop 'Invalid meshes used in copy_mesh'

  !$OMP DO
  do i=1, newmesh%N*newmesh%M
      newmesh%Te(i) = oldmesh%Te(i)
      newmesh%Th(i) = oldmesh%Th(i)
      newmesh%Ts(i) = oldmesh%Ts(i)
      newmesh%Ne(i) = oldmesh%Ne(i)
      newmesh%Nh(i) = oldmesh%Nh(i)

      newmesh%Map(1:4,i) = oldmesh%Map(1:4,i)
  end do
  !$OMP END DO

end subroutine copy_mesh

! interpolation bilineaire ponderee par les aires
subroutine bilinear_interpol_dual(mesh, dual, InvCellVol )
   use Maths_m
   use Profiler_m
   use Types_m
   implicit none
   type(MeshValues), intent(INOUT) :: mesh, dual
   real(8), intent(IN)             :: InvCellVol(mesh%M,mesh%N)

   integer :: i,j
   real(8) :: weight

   type(Profiler), save :: prof

   call Profiler_start(prof, 'BILINEAR_INTERP')

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

   !$OMP PARALLEL DEFAULT(NONE) SHARED(dual, InvCellVol, mesh) &
   !$OMP PRIVATE(weight)
   !$OMP DO  COLLAPSE(2)
   do j=1, dual%N
     do i=1, dual%M

       weight = M_ONE/ ( InvCellVol(i,j) + InvCellVol(i+1,j) + InvCellVol(i,j+1) + InvCellVol(i+1,j+1) )

       ! interpolation bilineaire ponderee par les aires
       dual%Ne(i+(j-1)*dual%M) = ( mesh%Ne(i  +(j-1)*mesh%M)* InvCellVol(i,j)           &
                                 + mesh%Ne(i+1+(j-1)*mesh%M)*InvCellVol(i+1,j)          &
                                 + mesh%Ne(i  +    j*mesh%M)*InvCellVol(i,j+1)          &
                                 + mesh%Ne(i+1+    j*mesh%M)*InvCellVol(i+1,j+1) ) * weight

       dual%Nh(i+(j-1)*dual%M) = ( mesh%Nh(i  +(j-1)*mesh%M)*InvCellVol(i,j)            &
                                 + mesh%Nh(i+1+(j-1)*mesh%M)*InvCellVol(i+1,j)          &
                                 + mesh%Nh(i  +    j*mesh%M)*InvCellVol(i,j+1)          &
                                 + mesh%Nh(i+1+    j*mesh%M)*InvCellVol(i+1,j+1) ) * weight

       dual%Te(i+(j-1)*dual%M) = ( mesh%Te(i  +(j-1)*mesh%M)*InvCellVol(i,j)            &
                                 + mesh%Te(i+1+(j-1)*mesh%M)*InvCellVol(i+1,j)          &
                                 + mesh%Te(i  +    j*mesh%M)*InvCellVol(i,j+1)          &
                                 + mesh%Te(i+1+    j*mesh%M)*InvCellVol(i+1,j+1) ) * weight

       dual%Th(i+(j-1)*dual%M) = ( mesh%Th(i  +(j-1)*mesh%M)*InvCellVol(i,j)            &
                                 + mesh%Th(i+1+(j-1)*mesh%M)*InvCellVol(i+1,j)          &
                                 + mesh%Th(i  +    j*mesh%M)*InvCellVol(i,j+1)          &
                                 + mesh%Th(i+1+    j*mesh%M)*InvCellVol(i+1,j+1) ) * weight

       dual%Ts(i+(j-1)*dual%M) = ( mesh%Ts(i  +(j-1)*mesh%M)*InvCellVol(i,j)            &
                                 + mesh%Ts(i+1+(j-1)*mesh%M)*InvCellVol(i+1,j)          &
                                 + mesh%Ts(i  +    j*mesh%M)*InvCellVol(i,j+1)          &
                                 + mesh%Ts(i+1+    j*mesh%M)*InvCellVol(i+1,j+1) ) * weight
       end do
     end do
     !$OMP END DO
     !$OMP END PARALLEL

   call profiler_stop(prof)
 end subroutine

 subroutine check_divergences(mesh, maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN, x, y, t, nbiter, nmin, Params)
   use Maths_m
   use Output_m
   use Types_m
   implicit none

   type(MeshValues), intent(in) :: mesh
   real(8), intent(in)          :: t, maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN
   real(8), intent(in)          :: x(mesh%M, mesh%N), y(mesh%M, mesh%N)
   integer, intent(in)          :: nbiter, nmin
   type(InputParameters), intent(in) :: Params

   integer(8)          :: i, j
   logical Diverged

   Diverged = .false.
    call output_open(ErrorFile%unit, 'output/error.log', (nbiter/=nmin .or. Params%RestartCalc == 1))
    !$OMP PARALLEL DEFAULT(NONE) SHARED(Diverged, x, y, mesh, t, ErrorFile)
    !$OMP DO
    do i=1,mesh%M*mesh%N
      if(isnan(mesh%Te(i))) then
        write(ErrorFile%unit,*) "Divergence of Te at t=", t, " ind=", i
        Diverged=.true.
      end if
      if(isnan(mesh%Th(i))) then
        write(ErrorFile%unit,*) "Divergence of Th at t=", t, " ind=", i
        Diverged=.true.
      end if
      if(isnan(mesh%Ts(i))) then
        write(ErrorFile%unit,*) "Divergence of Ts at t=", t, " ind=", i
        Diverged=.true.
      end if
      if(isnan(mesh%Ne(i))) then
        write(ErrorFile%unit,*) "Divergence of Ne at t=", t, " ind=", i
        Diverged=.true.
      end if
      if(isnan(mesh%Nh(i))) then
        write(ErrorFile%unit,*) "Divergence of Nh at t=", t, " ind=", i
        Diverged=.true.
      end if

      if(Diverged .eqv. .true.) then
        write(*,*) "Divergence detected. Please check error.dat for more information."
        call StopProgram ()
      end if
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    if(maxCFLxT.gt.M_ONE .OR. maxCFLyT.gt.M_ONE) then
      write(ErrorFile%unit,*) "Bad convergence for Te,Th. t=", t, "(CFLx,CFLy)=", maxCFLxT, maxCFLyT
      Diverged=.true.
    end if
    if(maxCFLxN.gt.M_ONE .OR. maxCFLyN.gt.M_ONE) then
      write(ErrorFile%unit,*) "Bad convergence for Ne,Nh. t=", t, "(CFLx,CFLy)=", maxCFLxN, maxCFLyN
      Diverged=.true.
    end if

 end subroutine check_divergences
