!------------------------------------------------------------------------------
!> @file mie.f90
!
! DESCRIPTION:
!> @brief A set of math functions for Mie scattering
!
!> @author
!> Thibault J.Y. Derrien
!
!> @date
!> 01 Jun 2016 - Initial Version
!> 07 Jun 2016 - Creating the Mie_m modules and moving routines - NTD
!------------------------------------------------------------------------------


module Mie_m
  use Maths_m
  implicit none


  !Some parameters for Bessel functions
  integer, private, parameter :: maxBesselOrder=20  !> Max of terms in series of Bessel for Mie scattering
  integer, private, parameter :: besselArray=1      !TODO: Explain what it is, if really useful


contains

    function MieScattering(r, phi, radius, dielectric, k)
      implicit none

      complex(8) :: MieScattering, dielectric
      complex(8) total

      real(8) :: r, phi, radius, k
      real(8) ireal
      integer(8) i

      total=M_ZERO
      ! Just test functions to validate
!       total=BesselJ(1d0, Unit*k*r) !test: success
!       total=BesselJ(-1d0, Unit*k*r) !test: success
!         total=Hankel1(1d0, Unit*k*r) !test: succes, undefined for z=0
!         total=Hankel1(-1d0, Unit*k*r) !test: success, undefined for z=0
!         total=BesselJprime(1d0, Unit*k*r) !test: success
!         total=BesselJprime(-1d0, Unit*k*r) !test: success
!         total=Hankel1prime(1d0, Unit*k*r) !test: success
!         total=Hankel1prime(-1d0, Unit*k*r) !test: failed, strong divergence while r->0

      do i=1, 2*maxBesselOrder+1
        ireal=real(i-maxBesselOrder-1) !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
!         write(*,*) i, ireal
        total=total+M_IM**ireal * exp(M_IM*ireal*phi) * BesselJ(ireal, &
                sqrt(dielectric)*k*r) * MieCoeff1(ireal, radius, dielectric, k)
      end do

      MieScattering=total
      return
    end function MieScattering

    function MieScatteringTE1(r, phi, radius, dielectric, k)
      implicit none

      complex(8) :: MieScatteringTE1, dielectric
      complex(8) :: total

      real(8) :: r, phi, radius, k
      real(8) :: ireal
      integer(8) :: i


      !!Careful !! This function is very sensitive to noise.

      total=M_ZERO

      do i=1, 2*maxBesselOrder+1
        ireal=real(i-maxBesselOrder-1) !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
!         write(*,*) i, ireal
        total=total + ( M_IM**ireal * exp(M_IM*ireal*phi) * BesselJ(ireal, &
              sqrt(dielectric)*k*r) * ireal * MieCoeff3(ireal, radius, dielectric, k) ) !original !!
      end do
      if(r.eq.0d0) then
        MieScatteringTE1=M_ZERO
      else
        MieScatteringTE1=-total/(dielectric*k*r)
!         MieScatteringTE1=Zero
      end if
      return
    end function MieScatteringTE1

    function MieScatteringTE2(r, phi, radius, dielectric, k)
      implicit none

      complex(8) :: MieScatteringTE2, dielectric
      complex(8) :: total

      real(8) :: r, phi, radius, k
      real(8) :: ireal
      integer(8) :: i

      total=M_ZERO

      do i=1, 2*maxBesselOrder+1
        ireal=real(i-maxBesselOrder-1) !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
!         write(*,*) i, ireal
        total=total + ( M_IM**ireal * exp(M_IM*ireal*phi) * &
              BesselJprime(ireal, sqrt(dielectric)*k*r) * MieCoeff3(ireal, radius, dielectric, k) )
      end do

      MieScatteringTE2=-total*M_IM/sqrt(dielectric)
!         MieScatteringTE2=Zero !debug
      return
    end function MieScatteringTE2

    function MieCoeff1(order, radius, dielectric, k)
      implicit none
      complex(8) MieCoeff1, dielectric
      real(8) :: k, radius
      real(8) :: order

!         MieCoeff1=Unit !debug
      MieCoeff1=(BesselJ(order, M_ONE*k*radius) - MieCoeff2(order, radius, dielectric, k)  &
          * Hankel1(order, M_ONE*k*radius)) / (BesselJ(order, k*radius*sqrt(dielectric)))
      return
    end function MieCoeff1

    function MieCoeff2(order, radius, dielectric, k)
      implicit none
      complex(8) :: MieCoeff2, dielectric
      real(8) :: k, radius
      real(8) :: order

        MieCoeff2= ( (sqrt(dielectric) * BesselJprime(order, k*radius*sqrt(dielectric)) &
                      * BesselJ(order, M_ONE*k*radius) ) - (BesselJ(order,sqrt(dielectric)*k*radius) &
                      *BesselJprime(order, M_ONE*k*radius)) ) &
                    / (( sqrt(dielectric)*BesselJprime(order,k*radius*sqrt(dielectric)) &
                        *Hankel1(order, M_ONE*k*radius) ) - ( BesselJ(order, sqrt(dielectric)*k*radius) &
                        *Hankel1prime(order, M_ONE*k*radius) )) !original

!         value1=radius*sqrt(dielectric)
!         MieCoeff2=BesselJprime(order, value1) !debug

!         MieCoeff2=BesselJprime(order, Unit*k*radius)
!          write(*,*) order, MieCoeff2
      return
    end function MieCoeff2

    function MieCoeff3(order, radius, dielectric, k)
      implicit none
      complex(8) MieCoeff3, dielectric
      real(8) :: k, radius
      real(8) :: order


!         MieCoeff3=Unit !debug
      MieCoeff3=(BesselJ(order, M_ONE*k*radius) - MieCoeff4(order, radius, dielectric, k)  &
          * Hankel1(order, M_ONE*k*radius)) / (BesselJ(order, k*radius*sqrt(dielectric)))
      return
    end function MieCoeff3

    function MieCoeff4(order, radius, dielectric, k)
      implicit none
      complex(8) :: MieCoeff4, dielectric
      real(8) :: k, radius
      real(8) :: order


        MieCoeff4= ( ( BesselJprime(order, k*radius*sqrt(dielectric))  &
                     * BesselJ(order, M_ONE*k*radius) ) - sqrt(dielectric) * (BesselJ(order,sqrt(dielectric)*k*radius) &
                     * BesselJprime(order, M_ONE*k*radius)) ) &
                    / (( BesselJprime(order,k*radius*sqrt(dielectric)) * Hankel1(order, M_ONE*k*radius) ) &
                    - sqrt(dielectric) * ( BesselJ(order, sqrt(dielectric)*k*radius) &
                    * Hankel1prime(order, M_ONE*k*radius) )) !original
!         MieCoeff4=Unit
      return
    end function MieCoeff4


  complex(8) function BesselJ(order, z)
      implicit none
      !
      complex(8), intent(in) :: z
      real(8),    intent(in) :: order
      !
      real(8) zR, zC
      integer(8) nz, ierr
      real(8) cyr(1:besselArray), cyi(1:besselArray)

      external ZBESJ
!       external ZABS

      cyr(:)=0.d0; cyi(:)=0.d0
      ierr=0; nz=0

      zR=real(z)
      zC=aimag(z)

!       write(*,*) (zR, zC)

!       CALL zbesj(1.d0, 0.d0, 0.d0, 1, besselArray, cyr, cyi, nz, ierr)

!       write(*,*) "Bessel 1", order
      CALL ZBESJ(zR, zC, abs(order), 1, besselArray, cyr, cyi, nz, ierr)

      if(ierr.ne.0) then
        write(*,*) "BesselJ is not well configured."
        write(*,*) z, cyr, cyi, ierr !, ZABS(zR, zC)
      end if
      BesselJ=M_ONE*cyr(besselArray)+M_IM*cyi(besselArray)

!       write(*,*) "Bessel", order

      if(order .lt. 0d0) then
        BesselJ=(-1d0)**(abs(order)) * BesselJ
      end if

      return
    end function BesselJ

   complex(8) function BesselJprime(order, z)
    implicit none
    !
    complex(8), intent(in) :: z
    real(8),    intent(in) :: order
    !
    external zbesj

      BesselJprime=0.5d0*(BesselJ(order-1d0,z)-BesselJ(order+1d0,z)) !Abramovitz, Eq. (9.1.27)

!! other form of the relation
!       if(z .eq. Zero) then
!         BesselJprime=Zero
!       else
!         BesselJprime=order*BesselJ(order,z)/z-BesselJ(order+1d0,z) !other form (still given in Abramovitz)
!       end if

      !debug
!       BesselJprime=Unit
!       end if
    end function BesselJprime

    complex(8) function Hankel1(order, z)
      implicit none
      !
      complex(8), intent(in) :: z
      real(8),    intent(in) :: order
      !
      real(8) zR, zC
      integer(8) nz, ierr
      real(8) cyr(1:besselArray), cyi(1:besselArray)
      external ZBESH
!       external ZABS
      cyr(:)=0.d0; cyi(:)=0.d0
      ierr=0; nz=0
      zR=real(z)
      zC=aimag(z)

!       write(*,*) zR, zC

      CALL ZBESH(zR, zC, abs(order), 1, 1, besselArray, cyr, cyi, nz, ierr)

      if(z .eq. M_ZERO) then
        Hankel1=M_ZERO
      else
        if(ierr.ne.0) then
          write(*,*) "Hankel1 is not well configured."
          write(*,*) z, order, ierr !, ZABS(zR, zC)
        end if
      end if
      Hankel1=M_ONE*cyr(besselArray)+M_IM*cyi(besselArray)

!       write(*,*) "Hankel", order, Hankel1

      if(order .lt. 0d0) then
        Hankel1=exp(M_IM*abs(order)*pi) * Hankel1
      end if
    return
    end function Hankel1

    complex(8) function Hankel1prime(order, z)
      implicit none
      complex(8), intent(in) :: z
      real(8),    intent(in) :: order
      external zbesh

      Hankel1prime=0.5d0*(Hankel1(order-1d0,z)-Hankel1(order+1d0,z))
! debug
!         Hankel1prime=Unit
    end function Hankel1prime


end module Mie_m
