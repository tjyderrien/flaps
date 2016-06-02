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
subroutine allocatemesh( mesh, M, N )
    use Types
    implicit none
    type(MeshValues) :: mesh
    integer :: M, N

    allocate(mesh%Te(M,N))
    allocate(mesh%Th(M,N))
    allocate(mesh%Ts(M,N))
    allocate(mesh%Ne(M,N))
    allocate(mesh%Nh(M,N))
end subroutine allocatemesh

!> Deallocate a MeshValues type
!> @param MeshValues Structure containing Mesh values
subroutine deallocatemesh( mesh )
    use Types
    implicit none
    type(MeshValues) :: mesh

    deallocate(mesh%Te)
    deallocate(mesh%Th)
    deallocate(mesh%Ts)
    deallocate(mesh%Ne)
    deallocate(mesh%Nh)
end subroutine deallocatemesh


