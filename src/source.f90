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

!> Initialize laser paramters
subroutine init_laser( laser )
 use Maths_m
 use Types_m
 implicit none

 type(LaserParams), intent(inout) :: laser

  laser%lambda    = 515d-9          !laser wavelength (m)
  laser%fluence   = 10d0            !laser fluence (J.m-2)
  laser%tau       = 40d-15          !FWHM pulse duration (s)
  laser%spotX     = 50d-6           !FWHM spot size in X direction (1030nm: 400nm x 50nm ; 515nm: 50um x 50 um ; 343 nm: 50um x 100nm)
  laser%spotY     = 50d-6           !FWHM spot size in Y direction
  laser%xCenter   = 1000d-9         ! X position of the max of the intensity (1030nm: 1um 0um, 515nm: idem, 343nm: 100nm x 200nm)
  laser%yCenter   = 0d0*200d-9      ! Y position of the max of the intensity

  laser%omega     = 2d0*Pi*c/laser%lambda !laser pulsation (s**-1)
  laser%k         = 2d0*pi/laser%lambda
  laser%inv_omega = 1.0d0/laser%omega
  laser%E         = hbar*laser%omega
  laser%inv_E     = 1.0d0/laser%E

end subroutine init_laser
