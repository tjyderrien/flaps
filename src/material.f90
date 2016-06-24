!------------------------------------------------------------------------------
!> @file material.f90
!
! DESCRIPTION:
!> @brief This file contains everything related to the properties and parameters of the material.
!
!> @author
!> Thibault J.Y. Derrien
!
!> @date
!> 15 Jun 2016 - Initial Version
!------------------------------------------------------------------------------

    real(8) function DensityOfState(mDOS, T)
    use Maths_m
    implicit none
    real(8), intent(in) :: mDOS, T
      DensityOfState = 2d0*(mDOS*kb*T/(2d0*pi*hbar**2))**(1.5d0)
      return
    end function DensityOfState

    subroutine TabCreateFL( M, N, FermiTableE, FermiTableH)
      implicit none
      integer, intent(in)    :: M, N
      real(8) :: FermiTableE(1:M,1:N), FermiTableH(1:M,1:N)

      integer :: unit1, unit2
      unit1=15; unit2=16
      open (unit1,file='FermiDatasE.dat')
      open (unit2,file='FermiDatasH.dat')
      read (unit1,*) FermiTableE(:,:) !, FermiTableE(2,:) !, FermiTableE(:,3), FermiTableE(:,4), &
!             FermiTableE(:,5), FermiTableE(:,6), FermiTableE(:,7), FermiTableE(:,8), &
!             FermiTableE(:,9)
      read (unit2,*) FermiTableH(:,:) !1), FermiTableH(:,2), FermiTableH(:,3), FermiTableH(:,4), &
!              FermiTableH(:,5), FermiTableH(:,6), FermiTableH(:,7), FermiTableH(:,8), &
!             FermiTableH(:,9)
! 222        format (1F10.2, 3x, 1F10.2, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x)
      close(unit1); close(unit2)

      !    fi_min=0.012d0
      !    fistep=0.2d0
      !    fi_max=fi_num*fistep
      write(*,*) FermiTableE(3,463), FermiTableH(3,450)
! stop
      return
    end subroutine TabCreateFL



    complex(8) function DielectricConstant(lambda)
      implicit none

      real(8), intent(in) :: lambda

      if(lambda.eq.1030d-9) then
        DielectricConstant=(12.8d0,0.001414418d0)
      end if

      if(lambda.eq.800d-9) then
        DielectricConstant=(13.46d0,0.048d0)
      end if

      if(lambda.eq.515d-9) then
        DielectricConstant=(17.8254d0,0.50669d0) !refractiveindex.info
      end if

      if(lambda.eq.343d-9) then
        DielectricConstant=(18.81766303d0,31.5464d0)
      end if
      return
    end function DielectricConstant

    complex(8) function DielectricFunction(epsilonInf, ne, nuColl, me, laser)
      use Maths_m
      use Types_m

      implicit none
      complex(8), intent(in) ::  epsilonInf
      real(8), intent(in) ::  ne, nuColl, me
      type(LaserParams), intent(in) :: laser

      real(8) omegape2

      omegape2=ec*ec*ne/me/epsilon0
      DielectricFunction=epsilonInf-omegape2*laser%inv_omega*laser%inv_omega/(M_ONE+M_IM*nuColl*laser%inv_omega)
      return
    end function DielectricFunction

    complex(8) function DielectricFunctionDrude(density, Collision, mass, laser)
      use Maths_m
      use Types_m

      implicit none

      real(8), intent(in)           :: density, Collision, mass
      type(LaserParams), intent(in) :: laser

      real(8) omegape2

      omegape2=ec*ec*density/(mass*epsilon0)
      DielectricFunctionDrude=M_ONE-M_ONE*omegape2*laser%inv_omega*laser%inv_omega/(M_ONE+M_IM*Collision*laser%inv_omega)
      return
    end function DielectricFunctionDrude

    real(8) function ephCollisionFrequency(ne)
      implicit none

      real(8), intent(in) :: ne

      real(8) nth
      nth=6.02d26 !m-3 (Sjodin, PRL 1998)
      ephCollisionFrequency=((240d-15)*(1d0+(ne/nth)**2))**(-1d0)
!       CollisionFrequency=1d14 !
      ! CollisionFrequency=1d13 !
      !CollisionFrequency=5d13 !
      return
    end function ephCollisionFrequency

    real(8) function CollisionFrequency()
      implicit none

      CollisionFrequency=1d15
      return
    end function CollisionFrequency

    real(8) function ImpactIonizationRate(Te, Ne, Ts, ImpactOff)
      use Maths_m
      implicit none

      real(8), intent(in)    :: Te, Ne, Ts
      integer(8), intent(in) :: ImpactOff

      real(8) :: EgapValue

      ImpactIonizationRate = 3.6d10*exp(-1.5d0*EgapValue(Ne, Ts)/kb/Te)
      if(ImpactOff.eq.1) then !TODO: This is dirty, should be putted outside
        ImpactIonizationRate=0d0
      end if
      return
    end function ImpactIonizationRate

    real(8) function OnePhotonIonizationRate()
      implicit none

!       OnePhotonIonizationRate=4d0*pi/laser%lambda*aimag(sqrt(epsilonLinear))
      OnePhotonIonizationRate = 3.4536819356d6 !extracted from WC Dash and R Newman, Phys Rev 99, 1151 (1955)
      return
    end function OnePhotonIonizationRate

     real(8) function TwoPhotonIonizationRate(lambda)
      implicit none

      real(8), intent(in) :: lambda
      if(lambda.eq.1030d-9) then
        TwoPhotonIonizationRate=1.933288399d-11
      end if

      if(lambda.eq.800d-9) then
        TwoPhotonIonizationRate=1.857135194d-11
!         TwoPhotonIonizationRate=0d0
      end if

      if(lambda.eq.515d-9) then
        TwoPhotonIonizationRate=1.512238197d-11
!                TwoPhotonIonizationRate=0d0
      end if

      if(lambda.eq.343d-9) then
        TwoPhotonIonizationRate=0d0
      end if
      return
    end function TwoPhotonIonizationRate

    real(8) function EgapValue(Ne, Ts)
      use Maths_m

      implicit none

      real(8), intent(in) :: Ne, Ts

      EgapValue=ec*1.16d0
!        EgapValue=ec*(1.1692d0-4.9d-4*Ts**2/(Ts+655d0)-1.5d-10*Ne**(1d0/3d0)) !Korfiatis 2007

!      EgapValue=ec*(1.16d0-(7.02d-4*Ts**2)/(Ts+1108d0)-1.5d-10*Ne**(0.33333d0)) !Driel 1987
!       EgapValue=ec*(1.1692d0-4.9d-4*Ts**2/(Ts+655d0))
!         EgapValue=ec*(1.1692d0) !-4.9d-4*Ts**2/(Ts+655d0))

      if(EgapValue < 0d0) then
        EgapValue=0d0
      end if
      return
    end function EgapValue

    real(8) function LatticeHeatCapacity(T, SiDensity)
      implicit none
      real(8) T, SiDensity
!       LatticeHeatCapacity=1d3*SiDensity*0.2703d0/(exp(63.456d0/T)+0.84586d0) !bad fit ...
!         LatticeHeatCapacity=1d3*SiDensity*0.412920554599445d0/(exp(88.1830102582422d0/T)-0.676494557497076d0) !!better fit on Flubacher BUT INDUCES A SUPER BUG (+170 K with 3rd order time integration).
!         LatticeHeatCapacity=1d6*(1.978d0+3.54d-4*T-3.68d0*T**(-2)) !Driel 1987 - not very good, BUT WORKS.
!       LatticeHeatCapacity=1d3*SiDensity*(0.899d0*dexp(5.455d-05*T)-0.959d0*dexp(-0.004218d0*T))! very good exp fit on Okothin, BUT INDUCES A nonlinearity at the beginning (+50 K with 3rd order integration)
!         LatticeHeatCapacity=1D3*SiDensity*(1.239d0*sin(0.001413d0*T-0.1806d0) + 0.3168d0*sin(0.003343d0*T+0.7648d0) + 0.01947d0*sin(0.00904d0*T+0.4528d0) + 0.04943d0*sin(0.007262d0*T-0.6642d0)) !Fitted on Otokhin, but is it stable ?
!         LatticeHeatCapacity=1D3*SiDensity*(2.36d-16*T**5 -1.707d-12*T**4 + 4.619d-09*T**3 -5.912d-06*T**2 + 0.003733d0*T -0.0494d0) ! 5th order polynomial fit on Okhonin
!         LatticeHeatCapacity=1d3*SiDensity*(0.4135d0*T-0.4071d0*T**1.002d0) !Driel style (1)
        LatticeHeatCapacity=1d3*SiDensity*(-0.003592d0*T+0.01458d0*T**0.8316d0) !Driel style (2, better ?)
    end function LatticeHeatCapacity

    integer(8) function FermiIndex(NeNc, FermiMaxLines)
      implicit none
      ! Input: Value of density/DOS
      ! returns the index to take in the Fermi files
      real(8) NeNc, NeNc0, dNeNc
      integer(8), intent(in) :: FermiMaxLines
      NeNc0=1d-38
      dNeNc=1.03d0 !NeNc=NeNc0*dNeNc**n

      FermiIndex=nint(log10(NeNc/NeNc0)/log10(dNeNc)+1d0)
      if(FermiIndex < 1 .OR. FermiIndex > FermiMaxLines) then
        write(*,*) "FermiIndex problem: NeNc=", NeNc, "FermiIndex=", FermiIndex
      end if
      return
    end function FermiIndex
