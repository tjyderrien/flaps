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
!> 31 Jul 2016 - NTD : Adding the routine ComputeConductivities_batch
!------------------------------------------------------------------------------


!------------------------------------------------------------------
    subroutine TabCreateFL(FermiMaxLines, FermiTableE, FermiTableH)
      implicit none
      integer, intent(in)    :: FermiMaxLines
      real(8) :: FermiTableE(1:9,1:FermiMaxLines), FermiTableH(1:9,1:FermiMaxLines)

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



!------------------------------------------------------------------
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

!------------------------------------------------------------------
    real(8) function DensityOfState(mDOS, T)
    use Maths_m
    implicit none
    real(8), intent(in) :: mDOS, T
      DensityOfState = 2d0*(mDOS*kb*T/(2d0*pi*hbar**2))**(1.5d0)
      return
    end function DensityOfState

!------------------------------------------------------------------
  !Routine that computes both electron and hole density of states, on the full grid
  subroutine DensitiesOfState_batch(mesh, DOSe, DOSh, meDOS, mhDOS)
    use Maths_m
    use Types_m
    implicit none

    type(MeshValues),  intent(in)    :: mesh
    real(8),           intent(inout)    :: DOSe(mesh%M,mesh%N), DOSh(mesh%M,mesh%N)
    real(8), intent(in) :: meDOS, mhDOS

    real(8) :: coefE, coefH
    integer :: i, j

    coefE = meDOS*kb/(2d0*pi*hbar**2)
    coefH = mhDOS*kb/(2d0*pi*hbar**2)

    !$OMP DO COLLAPSE(2)
    do j=2, mesh%N-1 !(optimized)
      do i=2, mesh%M-1
        DOSe(i,j) = 2d0*(coefE*mesh%Te(i,j))**(1.5d0)
        DOSh(i,j) = 2d0*(coefH*mesh%Th(i,j))**(1.5d0)
      end do
    end do
    !$OMP END DO

  end subroutine DensitiesOfState_batch


!------------------------------------------------------------------
    !TODO: Do we need DielectricFunctionDrude? Is it just possible to compute it with epsilonInf=1 ?
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

!------------------------------------------------------------------
   !This routine computes the dielectric function for the entire grid with one call
    subroutine DielectricFunction_batch(mesh, N, Dielectric, epsilonInf, nuColl, me, laser)
      use Maths_m
      use Types_m
      implicit none

      type(MeshValues),  intent(in)    :: mesh
      real(8),           intent(in)    :: N(mesh%M,mesh%N)
      complex(8),        intent(inout) :: Dielectric(mesh%M, mesh%N)
      complex(8),        intent(in)    :: epsilonInf
      real(8),           intent(in)    :: nuColl, me
      type(LaserParams), intent(in)    :: laser

      complex(8) :: coef
      integer :: i, j

      coef=ec*ec/me/epsilon0*laser%inv_omega**2/(M_ONE+M_IM*nuColl*laser%inv_omega)

      !$OMP DO COLLAPSE(2)
      do j=2, mesh%N-1 !(optimized)
        do i=2, mesh%M-1
          Dielectric(i,j) =epsilonInf-N(i,j)*coef
        end do
      end do
      !$OMP END DO

    end subroutine DielectricFunction_batch

!------------------------------------------------------------------
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

!------------------------------------------------------------------
       !This routine computes the Drude dielectric function for the entire grid with one call
    subroutine DielectricFunctionDrude_batch(mesh, N, Dielectric, Collision, mass, laser)
      use Maths_m
      use Types_m
      implicit none

      type(MeshValues),  intent(in)    :: mesh
      real(8),           intent(in)    :: N(mesh%M,mesh%N)
      complex(8),        intent(inout) :: Dielectric(mesh%M,mesh%N)
      real(8),           intent(in)    :: Collision, mass
      type(LaserParams), intent(in)    :: laser

      complex(8) :: coef
      integer :: i, j

      coef= M_ONE*ec*ec/(mass*epsilon0)*laser%inv_omega**2/(M_ONE+M_IM*Collision*laser%inv_omega)

      !$OMP DO COLLAPSE(2)
      do j=2, mesh%N-1 !(optimized)
        do i=2, mesh%M-1
          Dielectric(i,j) =M_ONE-coef*N(i,j)
        end do
      end do
      !$OMP END DO

    end subroutine DielectricFunctionDrude_batch


!------------------------------------------------------------------

    !TODO: Create a batch version of this routine
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


    !TODO: Create a batch version of this routine
    real(8) function CollisionFrequency()
      implicit none

      CollisionFrequency=1d15
      return
    end function CollisionFrequency

    !TODO: Create a batch version of this routine
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


   !-------------------------------------------------------------------------------------
   !> Computes the electron and mobilities for the entire mesh
   !-------------------------------------------------------------------------------------
    subroutine ComputeMobilities_batch(mesh, mobilityE, mobilityH, &
                                       FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                       ColFermi0, ColFermiHalf, nuColl, me)
      use Maths_m
      use Types_m
      implicit none

      type(MeshValues),  intent(in)    :: mesh
      real(8),           intent(inout) :: mobilityE(mesh%M,mesh%N)
      real(8),           intent(inout) :: mobilityH(mesh%M,mesh%N)
      real(8),           intent(in)    :: FermiTableE(mesh%M,mesh%N)
      real(8),           intent(in)    :: FermiTableH(mesh%M,mesh%N)
      integer(8),        intent(in)    :: FermiIndexE(mesh%M,mesh%N)
      integer(8),        intent(in)    :: FermiIndexH(mesh%M,mesh%N)
      integer(8),        intent(in)    :: ColFermi0, ColFermiHalf
      real(8),           intent(in)    :: nuColl, me

      integer :: i, j
      real(8) :: coef

      coef = ec/(me*nuColl)

      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N !(optimized)
        do i=1, mesh%M
          mobilityE(i,j)=coef*FermiTableE(ColFermi0, FermiIndexE(i,j))/FermiTableE(ColFermiHalf, FermiIndexE(i,j))
          mobilityH(i,j)=coef*FermiTableH(ColFermi0, FermiIndexH(i,j))/FermiTableH(ColFermiHalf, FermiIndexH(i,j))
        end do
      end do
      !$OMP END DO

    end subroutine ComputeMobilities_batch

   !-------------------------------------------------------------------------------------
   !> Computes the conductivities for the entire mesh
   !-------------------------------------------------------------------------------------
    subroutine ComputeConductivities_batch(mesh, kappae, kappah, kappas, mobilityE, mobilityH, &
                                           FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                           ColFermi0, ColFermi1, ColFermi2, ConductivityFix)
      use Maths_m
      use Types_m
      implicit none

      type(MeshValues),  intent(in)    :: mesh
      real(8),           intent(inout) :: kappae(mesh%M,mesh%N)
      real(8),           intent(inout) :: kappah(mesh%M,mesh%N)
      real(8),           intent(inout) :: kappas(mesh%M,mesh%N)
      real(8),           intent(in)    :: mobilityE(mesh%M,mesh%N)
      real(8),           intent(in)    :: mobilityH(mesh%M,mesh%N)
      real(8),           intent(in)    :: FermiTableE(mesh%M,mesh%N)
      real(8),           intent(in)    :: FermiTableH(mesh%M,mesh%N)
      integer(8),        intent(in)    :: FermiIndexE(mesh%M,mesh%N)
      integer(8),        intent(in)    :: FermiIndexH(mesh%M,mesh%N)
      integer(8),        intent(in)    :: ColFermi0, ColFermi1, ColFermi2, ConductivityFix

      integer :: i, j

      !Elena Silaeva fit on: Kazan et al, Journal of Applied Physics, 2010, 107, 083503
      real(8), parameter :: aa = -8.992d0
      real(8), parameter :: bb = 68.265d0
      real(8), parameter :: cc = -.4075612391d0
      real(8), parameter :: dd = .315984470d0
      real(8), parameter :: ee = -.4756634637d0
      real(8), parameter :: ff = 2.403533689d0 !TODO: If possible, use notations of the original paper

      if(ConductivityFix.eq.-1) then
        kappae(:,:) = 0.d0
        kappah(:,:) = 0.d0
        kappas(:,:) = 0.d0
        return
      end if

      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N
        do i=1, mesh%M

          !TODO: These FermiTable etc, can we precompute them?
          ! thermal coefficients
          kappae(i,j)=kb2*inv_ec*mesh%Ne(i,j)*mobilityE(i,j)*mesh%Te(i,j)* &
             ( 6d0*  FermiTableE(ColFermi2,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) &
              -4d0*( FermiTableE(ColFermi1,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)))**2 )

          kappah(i,j)=kb2*inv_ec*mesh%Nh(i,j)*mobilityH(i,j)*mesh%Th(i,j)* &
            (  6d0*  FermiTableH(ColFermi2,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) &
              -4d0*( FermiTableH(ColFermi1,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)))**2)
!         kappas(i,j)=-.1412d0*Ts(i,j)**(1.38961d0)+0.638157d0*Ts(i,j)**(1.14013d0) !mingo till 300 K, Nano Letters, 2003, 3, 1713-1716

          !Elena Silaeva fit on: Kazan et al, Journal of Applied Physics, 2010, 107, 083503
          kappas(i,j)=max(0.d0, &
                    (aa + bb/(1d0+exp(cc-1.0d0*mesh%Ts(i,j)+dd))*(1d0-1d0/(1d0+exp(ee-2.0d0*mesh%Ts(i,j)+ff)))))
!        ! correction considering Fick diffusion in energy
!         if(ConductivityFix.eq.1) then
!             kappae(i,j)=kappae(i,j) + kb2*Te(i,j)*Ne(i,j)*mobilityE(i,j) / ec &
!                       * (etae - 2d0*FermiTableE(ColFermi1,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) )**2
!             kappah(i,j)=kappah(i,j) + kb2*Th(i,j)*Nh(i,j)*mobilityH(i,j) / ec &
!                       * (etah - 2d0*FermiTableH(ColFermi1,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) )**2
!         else if(ConductivityFix.eq.2) then
!             kappae(i,j)=kappae(i,j) + 2d0*kb2*Te(i,j)*FermiTableE(ColFermi1,FermiIndexE(i,j))*mobilityE(i,j)*FermiTableE(ColFermiHalf, FermiIndexE(i,j))*Ne(i,j) * &
!                         (2d0*FermiTableE(ColFermi1, FermiIndexE(i,j))*FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)) &
!                         /FermiTableE(ColFermiHalf,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) - 1.5d0) * &
!                         (FermiTableE(ColFermi0, FermiIndexE(i,j))*ec*FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))**(-1e0)
!
!             kappah(i,j)=kappah(i,j) + 2d0*kb2*Te(i,j)*FermiTableH(ColFermi1,FermiIndexH(i,j))*mobilityH(i,j)*FermiTableH(ColFermiHalf, FermiIndexH(i,j))*Ne(i,j) * &
!                         (2d0*FermiTableH(ColFermi1, FermiIndexH(i,j))*FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)) &
!                         /FermiTableH(ColFermiHalf,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) - 1.5d0) * &
!                         (FermiTableH(ColFermi0, FermiIndexH(i,j))*ec*FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)))**(-1.)
        end do
      end do
      !$OMP END DO

    end subroutine ComputeConductivities_batch


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

    !TODO: Create a batch version of this routine
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

    !TODO: Create a batch version of this routine
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

    !TODO: Create a batch version of this routine
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
