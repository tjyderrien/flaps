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
!------------------------------------------------------------------------------


!> Allocate a MeshValues type
!> @param MeshValues Structure containing Mesh values
subroutine initmesh( mesh, M, N )
    use Types
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
    use Types
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
  use Types
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
   use Types
   implicit none
   type(MeshValues), intent(INOUT) :: mesh, dual
   real(8), intent(IN)             :: InvCellVol(mesh%M,mesh%N)

   integer :: i,j
   real(8) :: weight

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
