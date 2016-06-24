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

module Types_m

  !> Parameters for defining the mesh
  type MeshValues
    real(8), allocatable, dimension( :, :) :: Te !> electron temperature
    real(8), allocatable, dimension( :, :) :: Th !> hole temperature
    real(8), allocatable, dimension( :, :) :: Ts !> lattice temperature
    real(8), allocatable, dimension( :, :) :: Ne !> electron density
    real(8), allocatable, dimension( :, :) :: Nh !> hole density

    integer :: M,N !> Mesh dimension
   end type MeshValues

   !> Paramters of the laser
   type LaserParams
     real(8) :: lambda    !> laser wavelength (m)
     real(8) :: fluence   !> laser fluence (J.m-2)
     real(8) :: tau       !> FWHM pulse duration (s)
     real(8) :: spotX     !> FWHM spot size in X direction (1030nm: 400nm x 50nm ; 515nm: 50um x 50 um ; 343 nm: 50um x 100nm)
     real(8) :: spotY     !> FWHM spot size in Y direction
     real(8) :: xCenter   !> X position of the max of the intensity (1030nm: 1um 0um, 515nm: idem, 343nm: 100nm x 200nm)
     real(8) :: yCenter   !> Y position of the max of the intensity
     real(8) :: omega     !> Laser frequency
     real(8) :: k
     real(8) :: inv_omega
     real(8) :: E, inv_E  !> Laser Photon energy, and its inverse
   end type LaserParams


end module Types_m

