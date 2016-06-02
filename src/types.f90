!------------------------------------------------------------------------------
!> @file mesh.f90
!
! MODULE: Types
!
!> @author
!> Nicolas Tancogne-Dejean
!
! DESCRIPTION:
!> This module defines various types used in the code
!
!> @date
!> 01 Jun 2016 - Initial Version
!------------------------------------------------------------------------------

module Types

  !> Parameters for defining the mesh
  type MeshValues
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :, :) :: Te !> electron temperature
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :, :) :: Th !> hole temperature
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :, :) :: Ts !> lattice temperature
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :, :) :: Ne !> electron density
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :, :) :: Nh !> hole density

    INTEGER :: M,N !> Mesh dimension
   end type MeshValues

end module Types

