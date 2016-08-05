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


   !------------------------------------------------------------------
   !This routine computes the intensity for the entire grid with one call
   !------------------------------------------------------------------
    subroutine ComputeIntensity_batch(Params, mesh, laser, intensity, Dielectric, OpticalIndex, OpticalDamping, Reflectivity, &
                                      absorptionDrudeE, absorptionDrudeH, OnePhotonIonizationRate0, TwoPhotonIonizationRate0, &
                                      t, t0, sigmaTau, I0, sigmaX, sigmaY, x, y, x0, y0, DefectThickness, BandBendingInFDTD,  &
                                      sigmaX1, sigmaY1, sigmaX2, sigmaY2, sigmaX3, sigmaY3, sigmaX4, sigmaY4, sigmaX5, sigmaY5, &
                                      sigmaX6, sigmaY6, sigmaX7, sigmaY7, sigmaX8, sigmaY8, sigmaX9, sigmaY9, x1, y1, x2, y2, &
                                      x3, y3, x4, y4, x5, EintFieldR,  &
                                      y5, x6, y6, x7, y7, x8, y8, x9, y9, I1, I2, I3, I4, I5, I6, I7, I8, I9 )
      use Maths_m
      use Types_m
      implicit none

      type(InputParameters),  intent(in) :: Params
      type(MeshValues),  intent(in)      :: mesh
      type(LaserParams), intent(in)      :: laser
      real(8),        intent(inout)      :: intensity(mesh%M, mesh%N)
      complex(8),        intent(in)      :: Dielectric(mesh%M, mesh%N)
      real(8),           intent(in)      :: OpticalIndex(mesh%M, mesh%N)
      real(8),           intent(in)      :: OpticalDamping(mesh%M, mesh%N)
      real(8),           intent(in)      :: Reflectivity(mesh%M, mesh%N)
      real(8),           intent(in)      :: absorptionDrudeE(mesh%M, mesh%N)
      real(8),           intent(in)      :: absorptionDrudeH(mesh%M, mesh%N)
      real(8),           intent(in)      :: x(mesh%M, mesh%N), y(mesh%M, mesh%N)
      real(8),           intent(in)      :: EintFieldR(mesh%M, mesh%N)
      real(8),           intent(in)      :: OnePhotonIonizationRate0, TwoPhotonIonizationRate0
      real(8),           intent(in)      :: t, t0, sigmaTau, I0, sigmaX, sigmaY, x0, y0, DefectThickness
      integer(8),        intent(in)      :: BandBendingInFDTD
      real(8),           intent(in)      :: sigmaX1, sigmaY1, sigmaX2, sigmaY2, sigmaX3, sigmaY3, sigmaX4, sigmaY4, &
                                            sigmaX5, sigmaY5, sigmaX6, sigmaY6, sigmaX7, sigmaY7, sigmaX8, sigmaY8, &
                                            x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8, &
                                            I1, I2, I3, I4, I5, I6, I7, I8, sigmaX9, sigmaY9, x9, y9, I9

      integer :: i, j
      real(8) :: ConstBLx, ConstBLy
      real(8) :: OnePhotonIonizationRate

      !TODO: It is almost impossible to read, and per se to debug such a code.

      !TODO: @TYJD: Stop doing such coding style vandalism ;)

      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N !(optimized)
        do i=1, mesh%M
            !local intensity
        if(Params%UseMieScattering .eq. -1) then
  !         ! DEBUG ZONE
  ! !         if(laser%lambda.eq.343d-9) then
  !         ! uniform distribution like in Elena's paper
          intensity(i,j)=(1d0-0e0*reflectivity(i,j))*OpticalIndex(i,j)*I0*exp(-.5d0*((t-t0)/sigmaTau)**2)
  !         intensity(i,j)=I0*exp(-.5d0*((t-t0)/sigmaTau)**2)*exp(-0.5d0*(((y(i,j)-500d-9)/sigmaY)**2+(x(i,j)/sigmaX)**2))
  !         ! with just nothing
  ! !           intensity(i,j)=(1d0-reflectivity(i,j))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2)*exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2)
  !         ! with beer-lambert
  ! !         intensity(i,j)=(-(absorptionDrude(i,j) &
  ! !                         +4d0*pi/laser%lambda*aimag(sqrt(epsilonInf)))*intensity(i,j-1) &
  ! !                         -TwoPhotonIonizationRate0*intensity(i,j-1)**2 &
  ! !                         )*x &
  ! !                         +intensity(i,j-1)
  ! !          end if

        else if(Params%UseMieScattering .eq. 0) then
          ! WITH EXTERNALLY ADJUSTED INPUTS
  !        !Lumerical mode already contains the reflectivity. Although, it doesn't consider change of optical index with ionization.
          if(laser%lambda.eq.1030d-9) then
            ConstBLx=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate0+1d0*TwoPhotonIonizationRate0) &
                      / (1d0*exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*x0) &
                      * (OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
            ConstBLy=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate0+1d0*TwoPhotonIonizationRate0) &
                      / (1d0*exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*y0) &
                      * (OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
  ! case 1030 nm distribution
            !initial field distribution
            if(BandBendingInFDTD.eq.0) then
              intensity(i,j)=(1d0-0d0*reflectivity(i,Params%N))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2) &
                      *exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2)
            else
            intensity(i,j)=(1d0-0d0*reflectivity(i,Params%N))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2) &
                            *(exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
                            + I0*1d-4*( &
        exp(-0.5d0*(((x(i,j)-x(i,Params%N))**2+(y(i,j)-y(i,Params%N))**2)/((DefectThickness)/(2d0*sqrt2ln2))**2)) &
       +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness)/(2d0*sqrt2ln2))**2)) &
       +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness)/(2d0*sqrt2ln2))**2)) &
                            ))
            end if
            ! corrections from FDTD calculations and recovering non-linear processes
            intensity(i,j)=intensity(i,j)*((OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
                          / (-TwoPhotonIonizationRate0+exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j) &
                              +absorptionDrudeH(i,j))*(x(i,j)-x0))  * ConstBLx * (OnePhotonIonizationRate0 + &
                              absorptionDrudeE(i,j) + absorptionDrudeH(i,j)))) &
                          *((OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
                          / (-TwoPhotonIonizationRate0+exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j) &
                              +absorptionDrudeH(i,j))*abs(y(i,j)-y0))  * ConstBLy * (OnePhotonIonizationRate0 + &
                              absorptionDrudeE(i,j) + absorptionDrudeH(i,j))))
          end if
  ! case 515 nm distribution
          if(laser%lambda.eq.515d-9) then
            ConstBLx=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate0+1d0*TwoPhotonIonizationRate0) &
                      / (1d0*exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*x0) &
                      * (OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
            ConstBLy=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate0+1d0*TwoPhotonIonizationRate0) &
                      / (1d0*exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*y0) &
                      * (OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
            intensity(i,j)=(1d0-0d0*reflectivity(i,Params%N))*exp(-.5d0*((t-t0)/sigmaTau)**2) & !TODO: How could you code something like that??
                          *( &
                          I1*exp(-.5d0*((x(i,j)-x1)/sigmaX1)**2)*exp(-.5d0*((y(i,j)-y1)/sigmaY1)**2) + &
                          I2*exp(-.5d0*((x(i,j)-x2)/sigmaX2)**2)*exp(-.5d0*((y(i,j)-y2)/sigmaY2)**2) + &
                          I3*exp(-.5d0*((x(i,j)-x3)/sigmaX3)**2)*exp(-.5d0*((y(i,j)-y3)/sigmaY3)**2) + &
                          I4*exp(-.5d0*((x(i,j)-x4)/sigmaX4)**2)*exp(-.5d0*((y(i,j)-y4)/sigmaY4)**2) + &
                          I5*exp(-.5d0*((x(i,j)-x5)/sigmaX5)**2)*exp(-.5d0*((y(i,j)-y5)/sigmaY5)**2) + &
                          I6*exp(-.5d0*((x(i,j)-x6)/sigmaX6)**2)*exp(-.5d0*((y(i,j)-y6)/sigmaY6)**2) + &
                          I7*exp(-.5d0*((x(i,j)-x7)/sigmaX7)**2)*exp(-.5d0*((y(i,j)-y7)/sigmaY7)**2) + &
                          I8*exp(-.5d0*((x(i,j)-x8)/sigmaX8)**2)*exp(-.5d0*((y(i,j)-y8)/sigmaY8)**2) + &
                          I9/2d0 &
  !                         *(1d0+cos(2d0*pi*x(i,j)/periodX)*sin(2d0*pi*y(i,j)/periodY)) &
                          *exp(-.5d0*((x(i,j)-x9)/sigmaX9)**2)*exp(-.5d0*((y(i,j)-y9)/sigmaY9)**2)) &
                          *((OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
                          / (-TwoPhotonIonizationRate0+exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j) &
                              +absorptionDrudeH(i,j))*(x(i,j)-x0))  * ConstBLx * (OnePhotonIonizationRate0 + &
                              absorptionDrudeE(i,j) + absorptionDrudeH(i,j)))) &
                          *((OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
                          / (-TwoPhotonIonizationRate0+exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j) &
                              +absorptionDrudeH(i,j))*abs(y(i,j)-y0))  * ConstBLy * (OnePhotonIonizationRate0 + &
                              absorptionDrudeE(i,j) + absorptionDrudeH(i,j))))
          end if
          if(laser%lambda.eq.343d-9) then
            intensity(i,j)= (1d0-0e0*reflectivity(i,Params%N))* &
                            I0*exp(-.5d0*((t-t0)/sigmaTau)**2) &
                            *( & !TODO: Use OnePhotonIonizationRate0 here
                            exp(-(OnePhotonIonizationRate()+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
                            *abs(y(i,j)-y(i,Params%N)) & !introduce discontinuity !
                            ) &
                            * exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
  !                           + exp(-(OnePhotonIonizationRate(laser%lambda, epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*abs(y(i,j)-y(i,1))) &
  !                           *exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
                            )
  !                           *exp(-(OnePhotonIonizationRate(laser%lambda, epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*abs(x(i,j)-x(i,N)))
          end if

        ! USING MIE SCATTERING ANALYTICAL FORMULAS
        else if(Params%UseMieScattering .eq. 1) then
        ! calculate electric field inside the tip
!           EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie, abs(ContourYofX(x(i,j), NeedleRadius, NeedleAngle)), Dielectric(i,j)) !*sqrt(2d0*laser%fluence/(c*epsilon0*tau))
          ! debug formula for constant cone radius
!           EintField(i,j)=MieScattering(abs(y(i,j)), phiMie, 100d-9, epsilonInf)
!           EintField(i,j)=sqrt(EintField(i,j)*conjg(EintField(i,j))) !complex to real
          intensity(i,j)=I0 *OpticalIndex(i,j)* EintFieldR(i,j)**2 * exp(-.5d0*((t-t0)/sigmaTau)**2) !laser laser%fluence and reflectivity is inside the field
        else
          write(*,*) "Input ERROR. Check the MieScattering parameter."
          stop
        end if

        if(intensity(i,j) < 1d-20) then
          intensity(i,j)=0d0
        end if
        end do
      end do
      !$OMP END DO

    end subroutine ComputeIntensity_batch


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
