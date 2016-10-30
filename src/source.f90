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
!> @file source.f90
!
! DESCRIPTION:
!> @brief This files contains routines for the laser source
!
!> @author
!> Nicolas Tancogne-Dejean
!
!> @date
!> 08 Jun 2016 - Initial Version
!------------------------------------------------------------------------------


module laser_m


  private

  public ::     &
    laser,      &
    init_laser


  !> Structure defining a laser
  type laser
    real(8) :: lambda    !> laser wavelength (m)
    real(8) :: fluence   !> laser fluence (J.m-2)
    real(8) :: tau       !> FWHM pulse duration (s)
    real(8) :: spotX     !> FWHM spot size in X direction (1030nm: 400nm x 50nm ; 515nm: 50um x 50 um ; 343 nm: 50um x 100nm)
    real(8) :: spotY     !> FWHM spot size in Y direction
    real(8) :: xCenter   !> X position of the max of the intensity (1030nm: 1um 0um, 515nm: idem, 343nm: 100nm x 200nm)
    real(8) :: yCenter   !> Y position of the max of the intensity
    real(8) :: omega     !> laser pulsation (s**-1)
    real(8) :: k
    real(8) :: inv_omega
    real(8) :: E, inv_E  !> Laser Photon energy, and its inverse
  end type laser

contains

  !> Initialize laser paramters
  subroutine init_laser( this )
   use Maths_m
   use Types_m
   implicit none

   type(Laser), intent(inout) :: this

    this%lambda    = 515d-9
    this%fluence   = 10d0
    this%tau       = 40d-15
    this%spotX     = 50d-6
    this%spotY     = 50d-6
    this%xCenter   = 1000d-9
    this%yCenter   = 0d0*200d-9

    this%omega     = 2d0*M_PI*c/this%lambda
    this%k         = 2d0*M_PI/this%lambda
    this%inv_omega = 1.0d0/this%omega
    this%E         = hbar*this%omega
    this%inv_E     = 1.0d0/this%E

  end subroutine init_laser

end module laser_m
