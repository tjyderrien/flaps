!! Copyright (C) 2012-2017 T. J.-Y. Derrien, N. Tancogne-Dejean
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
!------------------------------------------------------------------------------

module Material_m
  use Maths_m
  use Profiler_m
  use Types_m
  use laser_m

  implicit none

  private

  public ::                  &
    Material,                &
    Laser,                   &
    get_collision_frequency, &
    init_material,           &
    evaluate_bandgap,        &
    OnePhotonIonizationRate, &
    TwoPhotonIonizationRate

  !Material index
  integer, parameter ::   &
       Si  = 0,           &
       ZnO = 1

  !The different model for the band-gap value
  integer, parameter ::            &
       !Model for the band-gap of silicon
       EG_SI_CONSTANT      = 0,    &
       EG_SI_KORFIATIS07   = 1,    &
       EG_SI_VANDRIEL87    = 2,    &
       !Model for the band-gap of ZnO
       EG_ZNO_CONSTANT     = 3

  !Collision frequency
  real, parameter::            &
    COL_FREQ_SI_CONSTANT = 1d15, & !Silaeva et al, New Journal of Physics 15, 089401 (2013)
    COL_FREQ_ZNO_CONSTANT= 1d15    !Dufft and Bonse, Journal of Applied Physics 105, 034908 (2009)

  !Static dielectric permittivity. Useful in case of constant field before the laser irradiation.
  real(8), parameter::         &
    EPS_INF_SI    =  11.66570433d0, & !,0.01404457712d0)  !Palik, E. D. Handbook of Optical Constants of Solids Academic Press, 1985
    !TODO: Do you really believe al these digits? TJYD: The higher, the better. SPP spectroscopy has no limits ;)
    EPS_INF_ZNO   = (1.9d0)**2        !TODO: HASARDOUS, based on refractiveindex.info, with incomplete data. 

! === Not necessary ===
!   real(8), parameter ::        &
!     IR1P_SI       =  3.4536819356d6 
!     !TODO: This must be calculated calculated directly from the wavelenth-dependent dielectric Permittivity! 

  !Two-photon ionization rate
  !TODO: To be implemented

  type Material
    integer :: Id              !< The material ID
    integer :: Eg_model        !< The model for the band-gap
    real(8) :: AugerRateE      !< Auger rate for electrons
    real(8) :: AugerRateH      !< Auger rate for holes
    real(8) :: EpsStatic       !< dielectric constant for static field: useful when Poisson will be solved for static fields
    real(8) :: Dielectric	 !< dielectric permittivity at a given wavelength !TODO: How to put a function of wavelength and temperature here? 
  end type Material

  contains

    subroutine init_material(this, AugerOff, source)
      type(Material), intent(inout) :: this
      integer,        intent(in)    :: AugerOff
      type(Laser),    intent(in)    :: source
      
      complex(8) :: DielectricConstant

      !TODO: For the moment, this is hard-coded. In the future, this will be obtained from the input file
      this%Id = Si
      this%Eg_model = EG_SI_CONSTANT
      if(AugerOff.eq.0) then
        this%AugerRateE=2.3d-43
        this%AugerRateH=7.8d-44
      else
        this%AugerRateE=M_ZERO
        this%AugerRateH=M_ZERO
      end if

      !Lets select some hardcoded values, depending on the material
      select case(this%Id)
      case(Si)
        this%EpsStatic = DielectricConstant(source%lambda)
      case(ZnO)
        print *, 'Static value for ZnO not implemeted.'
        this%EpsStatic = EPS_INF_ZNO
        call StopProgram()
      case default
        print *, 'Bad value for material ID.'
        call StopProgram()
      end select


      !We need to check if the models are compatible with the material selected
      select case(this%Id)
      case(Si)
        if(this%Eg_model < EG_SI_CONSTANT .or. this%Eg_model > EG_SI_VANDRIEL87) then
          print *, 'Selected band-gap model is not compatible with silicon.'
          call StopProgram()
        end if
      case(ZnO)
        if(this%Eg_model < EG_ZNO_CONSTANT) then
          print *, 'Selected band-gap model is not compatible with ZnO.'
          call StopProgram()
        end if
      case default
        print *, 'Bad value for material ID.'
        call StopProgram()
      end select
    end subroutine init_material

    !------------------------------------------------------------------
    !> Evaluate the bandgap of a material from the density and the lattice temperature depending on selected model.
    !------------------------------------------------------------------
    subroutine evaluate_bandgap(this, mesh, N, Ts, Eg)
      type(Material),   intent(in)    :: this
      type(MeshValues), intent(in)    :: mesh
      real(8),          intent(in)    :: N(:,:), Ts(:,:)
      real(8),          intent(out)   :: Eg(:,:)

      integer :: i,j
      type(Profiler), save :: prof

      call Profiler_start(prof, 'BAND_GAP')

      !TODO: Is seems that these three models have a very similar parametrization.
      !This implies one implementation and coefficients outside

      select case(this%Eg_model)
      case(EG_SI_CONSTANT)!TODO: We need a proper reference for this
        Eg(1:mesh%M,1:mesh%N)=ec*1.16d0
      case(EG_SI_KORFIATIS07) !REF: Korfiatis, D. P., KA Th Thoma, and J. C. Vardaxoglou. "Conditions for femtosecond laser melting of silicon." Journal of Physics D: Applied Physics 40.21 (2007): 6803.
        !$OMP PARALLEL DO DEFAULT(NONE) SHARED (mesh, N, Ts, Eg) COLLAPSE(2)
        do j=1,mesh%N
          do i=1,mesh%M
            Eg(i,j)=ec*(1.1692d0-4.9d-4*Ts(i,j)**2/(Ts(i,j)+655d0)-1.5d-10*N(i,j)**(1d0/3d0)) !Korfiatis 2007
          end do
        end do
        !$OMP END PARALLEL DO
      case(EG_SI_VANDRIEL87) !REF: Van Driel, Henry M. "Kinetics of high-density plasmas generated in Si by 1.06-and 0.53-μm picosecond laser pulses." Physical Review B 35.15 (1987): 8166.
        !$OMP PARALLEL DO DEFAULT(NONE) SHARED (mesh, N, Ts, Eg) COLLAPSE(2)
        do j=1,mesh%N
          do i=1,mesh%M
            Eg(i,j)=ec*(1.1692d0-(7.02d-4*Ts(i,j)**2)/(Ts(i,j)+1108d0)-1.5d-10*N(i,j)**(1d0/3d0)) !Driel 1987
          end do
        end do
        !$OMP END PARALLEL DO
      case default
        Eg(1:mesh%M,1:mesh%N) = M_ZERO
      end select

      call profiler_stop(prof)
    end subroutine evaluate_bandgap

    !------------------------------------------------------------------
    !TODO: Create a batch version of this routine
    real(8) function get_collision_frequency(this) result(colfreq)
      type(Material),   intent(in)    :: this


      select case(this%Id)
      case(Si)
        colfreq=COL_FREQ_SI_CONSTANT
      case(ZnO)
        colfreq=COL_FREQ_ZNO_CONSTANT
      case default
        print *, 'Bad value for material ID.'
        call StopProgram()
      end select
    end function get_collision_frequency

    !------------------------------------------------------------------
    real(8) function OnePhotonIonizationRate(this, lambda)
      real(8),          intent(in)    :: lambda
      type(Material),   intent(in)    :: this
      
      complex(8) :: epsilonOmega, DielectricConstant
      
      epsilonOmega=DielectricConstant(lambda)
      
!       select case(this%Id)
!       case(Si)  !basically valid for any band gap material
	 OnePhotonIonizationRate = 4d0*M_PI/lambda*aimag(sqrt(epsilonOmega))
!       case(ZnO) 
!         OnePhotonIonizationRate = 4d0*M_PI/lambda*aimag(sqrt(epsilonOmega))
!       case default
!         print *, 'OnePhotonAbsorption: invalid material choice. '
!         call StopProgram()
!       end select
    end function OnePhotonIonizationRate

    !TODO: For Si, to be replaced by the model given in Bristow, Alan D., Nir Rotenberg, and Henry M. Van Driel. "Two-photon absorption and Kerr coefficients of silicon for 850–2200 nm." Appl. phys. lett 90.19 (2007): 191104.
    !We should also be able to select a tabulated Keldysh model from that.
    !------------------------------------------------------------------
    real(8) function TwoPhotonIonizationRate(this, lambda)
	 type(Material),   intent(in)    :: this
      real(8), intent(in) :: lambda

      select case (this%Id)
	   case (Si)
		if(lambda.eq.(1030d-9)) then
		  TwoPhotonIonizationRate=1.933288399d-11
		elseif (lambda.eq.(800d-9)) then
		  TwoPhotonIonizationRate=1.857135194d-11
		elseif(lambda.eq.(515d-9)) then
		  TwoPhotonIonizationRate=1.512238197d-11
		elseif(lambda.eq.(343d-9)) then
		  TwoPhotonIonizationRate=M_ZERO
		else
		  write(*,'(a)') 'Wavelength is not in database for material ID.'
		  call StopProgram()
		end if
	   case (ZnO)
		write(*,'(a)') 'Wavelength is not in database for material ID.'
		call StopProgram()
	   case default
		write(*,'(a)') 'Wavelength is not in database for material ID.'
		call StopProgram()
	 end select
    end function TwoPhotonIonizationRate

end module Material_m

!------------------------------------------------------------------
    complex(8) pure function DielectricConstant(lambda)
!> Dielectric constant for silicon material at some particular wavelengths. 
!> TODO: interface with SPP-extended-theory to obtain any value in spectrum
      implicit none
!TODO: NTD: This should not be hardcoded but should be in an external file. (Not clear how to do this properly).
!TODO: TJYD: The plan is to connect with SPP-extended-theory where Palik data [Palik, Edward D., ed. "Handbook of optical constants of solids." (1998).] are directly giving this coefficient. 
      real(8), intent(in) :: lambda

      if(lambda.eq.1030d-9) then
        DielectricConstant=(12.8d0,0.001414418d0)
        return
      end if

      if(lambda.eq.800d-9) then
        DielectricConstant=(13.64d0,0.048d0)
        return
      end if

      if(lambda.eq.515d-9) then
        DielectricConstant=(17.8254d0,0.50669d0) !refractiveindex.info
        return
      end if

      if(lambda.eq.343d-9) then
        DielectricConstant=(18.81766303d0,31.5464d0)
        return
      end if
    end function DielectricConstant

    !------------------------------------------------------------------
    !> Density of states for a 3D gas of electrons (sure?)
    pure real(8) function DensityOfState(mDOS, T)
      use Maths_m
      implicit none
      real(8), intent(in) :: mDOS, T

      DensityOfState = M_TWO*(mDOS*kb*T/(M_TWO*M_PI*hbar**2))**(1.5d0)
    end function DensityOfState

!------------------------------------------------------------------
  !Routine that computes both electron and hole density of states, on the full grid
  subroutine DensitiesOfState_batch(mesh, DOSe, DOSh, meDOS, mhDOS)
    use Maths_m
    use Profiler_m
    use Types_m
    use laser_m

    implicit none

    type(MeshValues),  intent(in)    :: mesh
    real(8),           intent(inout) :: DOSe(mesh%M,mesh%N), DOSh(mesh%M,mesh%N)
    real(8),           intent(in)    :: meDOS, mhDOS

    real(8) :: coefE, coefH
    integer :: i, j
    type(Profiler), save :: prof

    call Profiler_start(prof, 'DENSITY_OF_STATES')


    coefE = meDOS*kb/(M_TWO*M_PI*hbar**2)
    coefH = mhDOS*kb/(M_TWO*M_PI*hbar**2)

    !$OMP PARALLEL DEFAULT(NONE) SHARED (coefE, coefH, mesh, DOSe, DOSh)
    !$OMP DO COLLAPSE(2)
    do j=1, mesh%N !(optimized)
      do i=1, mesh%M
        DOSe(i,j) = M_TWO*sqrt((coefE*mesh%Te(i,j))**3)
        DOSh(i,j) = M_TWO*sqrt((coefH*mesh%Th(i,j))**3)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call profiler_stop(prof)

  end subroutine DensitiesOfState_batch


   !------------------------------------------------------------------
   !This routine computes the dielectric function for the entire grid with one call
    subroutine DielectricFunction_batch(mesh, Dielectric, OpticalIndex, OpticalDamping, Reflectivity, &
                                        epsilonInf, nuColl, me, source)
      use Maths_m
      use Profiler_m
      use Laser_m
      use Types_m
      implicit none

      type(MeshValues),  intent(in)    :: mesh
      complex(8),        intent(inout) :: Dielectric(mesh%M, mesh%N)
      real(8),        intent(inout)    :: OpticalIndex(mesh%M, mesh%N)
      real(8),        intent(inout)    :: OpticalDamping(mesh%M, mesh%N)
      real(8),        intent(inout)    :: Reflectivity(mesh%M, mesh%N)
      complex(8),        intent(in)    :: epsilonInf
      real(8),           intent(in)    :: nuColl, me
      type(Laser), intent(in)          :: source

      complex(8) :: coef, sqrtEps
      integer :: i, j
      type(Profiler), save :: prof

      call Profiler_start(prof, 'DIELECTRIC_FUNCTION')


      coef=ec*ec/me/epsilon0*source%inv_omega**2/(M_ONE+M_IM*nuColl*source%inv_omega)
      !$OMP PARALLEL DEFAULT(NONE) SHARED(mesh, Dielectric, epsilonInf, coef, &
      !$OMP OpticalIndex, OpticalDamping, Reflectivity ) &
      !$OMP PRIVATE(sqrtEps)
      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N !(optimized)
        do i=1, mesh%M
          Dielectric(i,j) =epsilonInf-mesh%Ne(i,j)*coef
          sqrtEps = sqrt(Dielectric(i,j))
          OpticalIndex(i,j)   = real(sqrtEps)
          OpticalDamping(i,j) = aimag(sqrtEps)
          Reflectivity(i,j)=  ( (OpticalIndex(i,j)-M_ONE)**2 + OpticalDamping(i,j)**2 ) &
                             /( (OpticalIndex(i,j)+M_ONE)**2 + OpticalDamping(i,j)**2 )
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      call profiler_stop(prof)
    end subroutine DielectricFunction_batch

!------------------------------------------------------------------
       !This routine computes the Drude dielectric function for the entire grid with one call
    subroutine ComputeDielectricFunctionDrude_batch(Params, mesh, N, Dielectric, absorptionDrude, Collision, mass, source)
      use Maths_m
      use Profiler_m
      use Laser_m
      use Types_m
      implicit none

      type(InputParameters), intent(in)    :: Params
      type(MeshValues),      intent(in)    :: mesh
      real(8),               intent(in)    :: N(mesh%M,mesh%N)
      complex(8),            intent(inout) :: Dielectric(mesh%M,mesh%N)
      real(8),               intent(inout) :: absorptionDrude(mesh%M,mesh%N)
      real(8),               intent(in)    :: Collision, mass
      type(Laser),           intent(in)    :: source

      complex(8) :: coef
      integer :: i, j
      type(Profiler), save :: prof

      !$OMP MASTER
      call Profiler_start(prof, 'DRUDE')
      !$OMP END MASTER

      coef= ec*ec/(mass*epsilon0)*source%inv_omega**2/(M_ONE+M_IM*Collision*source%inv_omega)

      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N !(optimized)
        do i=1, mesh%M
          Dielectric(i,j) =M_ONE-coef*N(i,j)
        end do
      end do
      !$OMP END DO

      if(Params%DrudeHeating==1) then
        !$OMP DO COLLAPSE(2)
        do j=1, mesh%N
          do i=1, mesh%M
           ! absorptionDrude(i,j)=2d0*source%k*aimag(sqrt(Dielectric(i,j)))
            absorptionDrude(i,j)=M_TWO*source%k*sqrt( M_HALF*( abs(Dielectric(i,j)) - real(Dielectric(i,j)) ) )
          end do
        end do
        !$OMP END DO
      end if

      !$OMP MASTER
      call profiler_stop(prof)
      !$OMP END MASTER

    end subroutine ComputeDielectricFunctionDrude_batch


!------------------------------------------------------------------

    !TODO: Create a batch version of this routine
    pure real(8) function ephCollisionFrequency(ne)
      use Maths_m
      implicit none

      real(8), intent(in) :: ne

      real(8), parameter :: inv_nth=M_ONE/6.02d26 !inversion of m-3

      !> [Sjodin, Theodore, Hrvoje Petek, and Hai-Lung Dai.
      !> "Ultrafast carrier dynamics in silicon: A two-color 
      !> transient reflection grating study on a (111) surface." 
      !> Physical review letters 81.25 (1998): 5664.]
      ephCollisionFrequency=M_ONE/((240d-15)*(M_ONE+(ne*inv_nth)**2))
    end function ephCollisionFrequency

!------------------------------------------------------------------

    !TODO: Create a batch version of this routine
    pure real(8) function ImpactIonizationRate(Te, Eg, ImpactOff)
      use Maths_m
      implicit none

      real(8), intent(in)    :: Te, Eg
      integer, intent(in)    :: ImpactOff

      real(8), parameter     :: inv_kb = -1.5d0/kb

      if(ImpactOff.eq.1) then
        ImpactIonizationRate=M_ZERO
        return
      end if

      !> Reference: [Driel, H. V. Kinetics of high-density plasmas generated in Si 
      !> by 1.06- and 0.53-$m picosecond laser pulses Phys. Rev. B, 1987, 35, 8166-8176
      ImpactIonizationRate = 3.6d10*exp(inv_kb*Eg/Te)

    end function ImpactIonizationRate


   !-------------------------------------------------------------------------------------
   !> Computes the electron and hole mobilities for the entire mesh
   !-------------------------------------------------------------------------------------
    subroutine ComputeMobilities_batch(mesh, mobilityE, mobilityH, &
                                       FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                       ColFermi0, ColFermiHalf, nuColl, me)
      use Maths_m
      use Profiler_m
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
      type(Profiler), save :: prof

      call Profiler_start(prof, 'MOBILITIES')

      coef = ec/(me*nuColl)

      !$OMP PARALLEL DEFAULT(NONE) SHARED (coef, mesh, mobilityE, mobilityH, FermiTableE, &
      !$OMP FermiTableH, ColFermi0, ColFermiHalf, FermiIndexE, FermiIndexH )
      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N !(optimized)
        do i=1, mesh%M
          mobilityE(i,j)=coef*FermiTableE(ColFermi0, FermiIndexE(i,j))/FermiTableE(ColFermiHalf, FermiIndexE(i,j))
          mobilityH(i,j)=coef*FermiTableH(ColFermi0, FermiIndexH(i,j))/FermiTableH(ColFermiHalf, FermiIndexH(i,j))
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      call profiler_stop(prof)

    end subroutine ComputeMobilities_batch

   !-------------------------------------------------------------------------------------
   !> Computes the diffusion terms for the entire mesh
   !-------------------------------------------------------------------------------------
    subroutine UpdateDiffusions_batch(Params, mesh, diffusionE, diffusionH, mobilityE, mobilityH, &
                                       FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                       ColFermiHalf, ColFermiMenusHalf)
      use Maths_m
      use Profiler_m
      use Types_m
      implicit none

      type(InputParameters), intent(in)    :: Params
      type(MeshValues),      intent(in)    :: mesh
      real(8),               intent(inout) :: diffusionE(mesh%M,mesh%N)
      real(8),               intent(inout) :: diffusionH(mesh%M,mesh%N)
      real(8),               intent(in)    :: mobilityE(mesh%M,mesh%N)
      real(8),               intent(in)    :: mobilityH(mesh%M,mesh%N)
      real(8),               intent(in)    :: FermiTableE(mesh%M,mesh%N)
      real(8),               intent(in)    :: FermiTableH(mesh%M,mesh%N)
      integer(8),            intent(in)    :: FermiIndexE(mesh%M,mesh%N)
      integer(8),            intent(in)    :: FermiIndexH(mesh%M,mesh%N)
      integer(8),            intent(in)    :: ColFermiMenusHalf, ColFermiHalf

      integer :: i, j
      type(Profiler), save :: prof

      if(Params%TransportModel.eq.-1) then
        return
      end if

      call Profiler_start(prof, 'DIFFUSIONS')


      !$OMP PARALLEL DEFAULT(NONE) SHARED (mesh, diffusionE, diffusionH, &
      !$OMP mobilityE, mobilityH, FermiTableE, FermiTableH, ColFermiHalf, &
      !$OMP FermiIndexE, FermiIndexH, ColFermiMenusHalf)
      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N
        do i=1, mesh%M
          diffusionE(i,j)=mobilityE(i,j)*kb*mesh%Te(i,j)*inv_ec &
                *FermiTableE(ColFermiHalf,FermiIndexE(i,j))/FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j))
          diffusionH(i,j)=mobilityH(i,j)*kb*mesh%Th(i,j)*inv_ec &
                *FermiTableH(ColFermiHalf,FermiIndexH(i,j))/FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j))
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      call profiler_stop(prof)

    end subroutine UpdateDiffusions_batch

   !-------------------------------------------------------------------------------------
   !> Computes the drift vectors for the entire mesh
   !-------------------------------------------------------------------------------------
    subroutine UpdateDriftVectors_batch(mesh, JeX, JeY, JhX, JhY, mobilityE, mobilityH, Ex, Ey, DriftOn)
      use Maths_m
      use Profiler_m
      use Types_m
      implicit none

      type(MeshValues),  intent(in)    :: mesh
      real(8),           intent(inout) :: JeX(mesh%M,mesh%N), JeY(mesh%M,mesh%N)
      real(8),           intent(inout) :: JhX(mesh%M,mesh%N), JhY(mesh%M,mesh%N)
      real(8),           intent(in)    :: mobilityE(mesh%M,mesh%N)
      real(8),           intent(in)    :: mobilityH(mesh%M,mesh%N)
      real(8),           intent(in)    :: Ex(mesh%M,mesh%N)
      real(8),           intent(in)    :: Ey(mesh%M,mesh%N)
      integer(8),        intent(in)    :: DriftOn

      integer :: i, j
      type(Profiler), save :: prof

      if(DriftOn.eq.0) then !No need to update these values
        return
      endif

      call Profiler_start(prof, 'DRIFT_VECTORS')


      !$OMP PARALLEL DEFAULT(NONE) SHARED (JeX, JeY, JhX, JhY, mobilityE, mobilityH, mesh, Ex, Ey)
      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N !(optimized)
        do i=1, mesh%M
          JeX(i,j)=-mobilityE(i,j)*mesh%Ne(i,j)*Ex(i,j)
          JeY(i,j)=-mobilityE(i,j)*mesh%Ne(i,j)*Ey(i,j)
          JhX(i,j)= mobilityH(i,j)*mesh%Nh(i,j)*Ex(i,j)
          JhY(i,j)= mobilityH(i,j)*mesh%Nh(i,j)*Ey(i,j)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      call profiler_stop(prof)
    end subroutine UpdateDriftVectors_batch

   !-------------------------------------------------------------------------------------
   !> Computes the electron, hole and lattice heat capacities for the entire mesh
   !-------------------------------------------------------------------------------------
    subroutine ComputeHeatCapacities_batch(mesh, Ce, Ch, Cs, invCe,invCh, invCs, &
                                           FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                           ColFermiThreeHalf, ColFermiHalf, ColFermiMenusHalf, ColFermiEta  )
      use Maths_m
      use Profiler_m
      use Types_m
      implicit none

      type(MeshValues),  intent(in)    :: mesh
      real(8),           intent(inout) :: Ce(mesh%M,mesh%N)
      real(8),           intent(inout) :: Ch(mesh%M,mesh%N)
      real(8),           intent(inout) :: Cs(mesh%M,mesh%N)
      real(8),           intent(inout) :: invCe(mesh%M,mesh%N)
      real(8),           intent(inout) :: invCh(mesh%M,mesh%N)
      real(8),           intent(inout) :: invCs(mesh%M,mesh%N)
      real(8),           intent(in)    :: FermiTableE(mesh%M,mesh%N)
      real(8),           intent(in)    :: FermiTableH(mesh%M,mesh%N)
      integer(8),        intent(in)    :: FermiIndexE(mesh%M,mesh%N)
      integer(8),        intent(in)    :: FermiIndexH(mesh%M,mesh%N)
      integer(8),        intent(in)    :: ColFermiThreeHalf, ColFermiHalf, ColFermiMenusHalf, ColFermiEta

      real(8) :: tmp, LatticeHeatCapacity
      integer :: i, j
      type(Profiler), save :: prof

      call Profiler_start(prof, 'HEAT_CAPACITIES')

      !$OMP PARALLEL DEFAULT(NONE) SHARED (mesh, Ce, Ch, Cs, ColFermiThreeHalf, FermiIndexE, FermiIndexH, &
      !$OMP FermiTableE, FermiTableH, ColFermiEta, ColFermiHalf, ColFermiMenusHalf, invCe, invCh, invCs) &
      !$OMP PRIVATE(tmp)
      !$OMP DO COLLAPSE(2)
      do j=1, mesh%N !(optimized)
        do i=1, mesh%M
          tmp = FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))
          Ce(i,j)=1.5d0*mesh%Ne(i,j)*kb*(tmp-FermiTableE(ColFermiEta,FermiIndexE(i,j)) &
                             *(M_ONE-(tmp/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))* &
                                     (FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))))/FermiTableE(ColFermiHalf,FermiIndexE(i,j))

          tmp = FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))
          Ch(i,j)=1.5d0*mesh%Nh(i,j)*kb*(tmp-FermiTableH(ColFermiEta,FermiIndexH(i,j)) &
                               *(M_ONE-(tmp/FermiTableH(ColFermiHalf,FermiIndexH(i,j)))* &
                                    (FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)))))/FermiTableH(ColFermiHalf,FermiIndexH(i,j))

          Cs(i,j)=LatticeHeatCapacity(mesh%Ts(i,j))

          invCe(i,j) = M_ONE/Ce(i,j)
          invCh(i,j) = M_ONE/Ch(i,j)
          invCs(i,j) = M_ONE/Cs(i,j)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      call profiler_stop(prof)
    end subroutine ComputeHeatCapacities_batch


    pure real(8) function LatticeHeatCapacity(T)
        implicit none
        real(8), intent(in) :: T

        real(8), parameter::  SiDensity=2.329d3            !Silicon rest density
!       LatticeHeatCapacity=1d3*SiDensity*0.2703d0/(exp(63.456d0/T)+0.84586d0) !bad fit ...
!         LatticeHeatCapacity=1d3*SiDensity*0.412920554599445d0/(exp(88.1830102582422d0/T)-0.676494557497076d0) !!better fit on Flubacher BUT INDUCES A SUPER BUG (+170 K with 3rd order time integration).
!         LatticeHeatCapacity=1d6*(1.978d0+3.54d-4*T-3.68d0*T**(-2)) !Driel 1987 - not very good, BUT WORKS.
!       LatticeHeatCapacity=1d3*SiDensity*(0.899d0*dexp(5.455d-05*T)-0.959d0*dexp(-0.004218d0*T))! very good exp fit on Okothin, BUT INDUCES A nonlinearity at the beginning (+50 K with 3rd order integration)
!         LatticeHeatCapacity=1D3*SiDensity*(1.239d0*sin(0.001413d0*T-0.1806d0) + 0.3168d0*sin(0.003343d0*T+0.7648d0) + 0.01947d0*sin(0.00904d0*T+0.4528d0) + 0.04943d0*sin(0.007262d0*T-0.6642d0)) !Fitted on Otokhin, but is it stable ?
!         LatticeHeatCapacity=1D3*SiDensity*(2.36d-16*T**5 -1.707d-12*T**4 + 4.619d-09*T**3 -5.912d-06*T**2 + 0.003733d0*T -0.0494d0) ! 5th order polynomial fit on Okhonin
!         LatticeHeatCapacity=1d3*SiDensity*(0.4135d0*T-0.4071d0*T**1.002d0) !Driel style (1)
        LatticeHeatCapacity=1d3*SiDensity*(-0.003592d0*T+0.01458d0*T**0.8316d0) !Driel style (2, better ?)
    end function LatticeHeatCapacity

   !-------------------------------------------------------------------------------------
   !> Updates the conductivities for the entire mesh
   !-------------------------------------------------------------------------------------
    subroutine UpdateConductivities_batch(mesh, kappae, kappah, kappas, mobilityE, mobilityH, &
                                           FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                           ColFermi0, ColFermi1, ColFermi2, TransportModel)
      use Maths_m
      use Profiler_m
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
      integer(8),        intent(in)    :: ColFermi0, ColFermi1, ColFermi2, TransportModel

      integer :: i, j
      !TODO: We have to find a nomeclature and a name for this model
      !TODO: Is there other models?
      !Elena Silaeva fit on: Kazan et al, Journal of Applied Physics, 2010, 107, 083503
      real(8), parameter :: aa = -8.992d0
      real(8), parameter :: bb = 68.265d0
      real(8), parameter :: cc = -.4075612391d0
      real(8), parameter :: dd = 2.315984470d0
      real(8), parameter :: ee = -.4756634637d0
      real(8), parameter :: ff = 2.403533689d0 !TODO: If possible, use notations of the original paper

      type(Profiler), save :: prof

      call Profiler_start(prof, 'CONDUCTIVITIES')


      if(TransportModel.eq.-1) then !No need to update the conductivity
        return
      else if (TransportModel .eq. 0 ) then
        !$OMP PARALLEL DEFAULT(NONE) SHARED (mesh, kappae, kappah, kappas, mobilityE, mobilityH, &
        !$OMP FermiTableE, FermiIndexE, FermiTableH, FermiIndexH, ColFermi1, ColFermi2, ColFermi0)
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
!           kappas(i,j)=-.1412d0*Ts(i,j)**(1.38961d0)+0.638157d0*Ts(i,j)**(1.14013d0) !mingo till 300 K, Nano Letters, 2003, 3, 1713-1716

            !Elena Silaeva fit on: Kazan et al, Journal of Applied Physics, 2010, 107, 083503
            kappas(i,j)=max(M_ZERO, &
                      (aa + bb/(1d0+exp(cc-M_ONE*mesh%Ts(i,j)+dd))*(M_ONE-M_ONE/(M_ONE+exp(ee-M_TWO*mesh%Ts(i,j)+ff)))))
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
!      else if(TransportModel.eq.1) then
!      TODO: Add reference (Tritt or Chen et al 2005 ?)
!        !$OMP DO COLLAPSE(2)
!        do j=1, mesh%N
!          do i=1, mesh%M
!             kappae(i,j)=kappae(i,j) + kb2*Te(i,j)*Ne(i,j)*mobilityE(i,j) / ec &
!                       * (etae - 2d0*FermiTableE(ColFermi1,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) )**2
!             kappah(i,j)=kappah(i,j) + kb2*Th(i,j)*Nh(i,j)*mobilityH(i,j) / ec &
!                       * (etah - 2d0*FermiTableH(ColFermi1,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) )**2
!
!           end do
!         end do
!        !$OMP END DO
!     else if(TransportModel.eq.2) then
!        !$OMP DO COLLAPSE(2)
!        do j=1, mesh%N
!          do i=1, mesh%M
!             kappae(i,j)=kappae(i,j) + 2d0*kb2*Te(i,j)*FermiTableE(ColFermi1,FermiIndexE(i,j))&
!                         *mobilityE(i,j)*FermiTableE(ColFermiHalf, FermiIndexE(i,j))*Ne(i,j) * &
!                         (2d0*FermiTableE(ColFermi1, FermiIndexE(i,j))*FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)) &
!                         /FermiTableE(ColFermiHalf,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) - 1.5d0) * &
!                         (FermiTableE(ColFermi0, FermiIndexE(i,j))*ec*FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))**(-1e0)
!
!             kappah(i,j)=kappah(i,j) + 2d0*kb2*Te(i,j)*FermiTableH(ColFermi1,FermiIndexH(i,j))&
!                         *mobilityH(i,j)*FermiTableH(ColFermiHalf, FermiIndexH(i,j))*Ne(i,j) * &
!                         (2d0*FermiTableH(ColFermi1, FermiIndexH(i,j))*FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)) &
!                         /FermiTableH(ColFermiHalf,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) - 1.5d0) * &
!                         (FermiTableH(ColFermi0, FermiIndexH(i,j))*ec*FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)))**(-1.)
!        end do
!      end do
!      !$OMP END DO
    end if

    call profiler_stop(prof)

    end subroutine UpdateConductivities_batch

   !-------------------------------------------------------------------------------------
   !> Computes the couplings for the entire mesh
   !-------------------------------------------------------------------------------------
    subroutine UpdateCouplings_batch(Params, mesh, CouplingE, CouplingH, Ce, Ch)
      use Maths_m
      use Profiler_m
      use Types_m
      implicit none

      type(InputParameters), intent(in)    :: Params
      type(MeshValues),      intent(in)    :: mesh
      real(8),               intent(inout) :: CouplingE(mesh%M,mesh%N)
      real(8),               intent(inout) :: CouplingH(mesh%M,mesh%N)
      real(8),               intent(in)    :: Ce(mesh%M,mesh%N)
      real(8),               intent(in)    :: Ch(mesh%M,mesh%N)
      !
      integer :: i, j
      real(8) :: nuColleph!        electron-phonon collision frequency
      real(8) :: ephCollisionFrequency
      !
      type(Profiler), save :: prof

      call Profiler_start(prof, 'COUPLINGS')

      if(Params%CouplingDebug.eq.1) then !No need to update the couplings, as they are zero
          return
      end if
      !TODO: why repeating the case on electrons, here? 
      !TODO: add parallel zone here
      if(Params%HolesOff.eq.0) then
        !$OMP PARALLEL DEFAULT(NONE) SHARED (CouplingE, CouplingH)
        !$OMP FIRSTPRIVATE (Ce, Ch, mesh, nuColleph)
        !$OMP DO COLLAPSE(2)
        do j=1, mesh%N !(optimized)
          do i=1, mesh%M
            ! optical coefficients
            nuColleph=ephCollisionFrequency(mesh%Ne(i,j))
            CouplingE(i,j)=Ce(i,j)*nuColleph*(mesh%Te(i,j)-mesh%Ts(i,j))
            CouplingH(i,j)=Ch(i,j)*nuColleph*(mesh%Th(i,j)-mesh%Ts(i,j))
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
      else !In this case no need to update CouplingH
       !$OMP PARALLEL DEFAULT(NONE) SHARED (CouplingE)
       !$OMP FIRSTPRIVATE (Ce, mesh)
       !$OMP DO COLLAPSE(2)
       do j=1, mesh%N !(optimized)
         do i=1, mesh%M
           ! optical coefficients
           nuColleph=ephCollisionFrequency(mesh%Ne(i,j))
           CouplingE(i,j)=Ce(i,j)*nuColleph*(mesh%Te(i,j)-mesh%Ts(i,j))
         end do
       end do
       !$OMP END DO
       !
      end if
      !
      call profiler_stop(prof)
    !
    end subroutine UpdateCouplings_batch



   !------------------------------------------------------------------
   !> Computes the Sources Gains and Losses terms for electron and holes
   !------------------------------------------------------------------
   subroutine ComputeGainsAndLosses(Params, mesh, source, matter, Egap, intensity, OnePhotonIonizationRate0, &
                              TwoPhotonIonizationRate0, absorptionDrudeE, absorptionDrudeH, &
                              Ce, Ch, CeOld, ChOld, dt, me, mh, &
                              GainsE, GainsH, SourceUe, SourceUh, SourceE, SourceH, LossesE, LossesH, ImpactOff )
     use Material_m
     use Maths_m
     use Profiler_m
     use Laser_m
     use Types_m
     implicit none

     type(InputParameters), intent(in)                 :: Params
     type(MeshValues),      intent(in)                 :: mesh
     type(Laser),           intent(in)                 :: source
     type(Material),        intent(in)                 :: matter
     real(8), dimension(mesh%M,mesh%N), intent(inout)  :: Egap, GainsE, GainsH, SourceUe, SourceUh, &
                                                          SourceE, SourceH, LossesE, LossesH
     real(8), dimension(mesh%M,mesh%N), intent(in)     :: intensity, absorptionDrudeE, absorptionDrudeH, &
                                                          Ce, Ch, CeOld, ChOld
     real(8),                           intent(in)     :: dt, me, mh, OnePhotonIonizationRate0, &
                                                          TwoPhotonIonizationRate0
     integer,                           intent(in)     :: ImpactOff

     real(8) :: Int2, ImpactIonizationRate, work
     integer :: i,j
     real(8), dimension(mesh%M,mesh%N) :: EgapH

     type(Profiler), save :: prof

     call Profiler_start(prof, 'GAINS_AND_LOSSES')

     ! free-carrier balance sources
     call evaluate_bandgap(matter, mesh, mesh%Ne, mesh%Ts, Egap)

     call evaluate_bandgap(matter, mesh, mesh%Nh, mesh%Ts, EgapH)

     ! calculation of sources
     !$OMP PARALLEL DEFAULT(NONE) SHARED (Params, Egap, EgapH, mesh, matter,       &
     !$OMP GainsE, SourceUe, SourceE, Ce, CeOld, dt, intensity, source,     &
     !$OMP OnePhotonIonizationRate0, TwoPhotonIonizationRate0, me, mh,      &
     !$OMP absorptionDrudeE, LossesE, ImpactOff,    &
     !$OMP GainsH, SourceUh, SourceH, LossesH, Ch, ChOld, absorptionDrudeH) &
     !$OMP PRIVATE(Int2, work)
     !$OMP DO  COLLAPSE(2)
     do j=1,Params%N
       do i=1,Params%M

        Int2 = intensity(i,j)**2
        work = ImpactIonizationRate(mesh%Te(i,j),Egap(i,j), ImpactOff)

        GainsE(i,j)=(OnePhotonIonizationRate0*intensity(i,j)*source%inv_E & 
                    +M_HALF*TwoPhotonIonizationRate0*Int2*source%inv_E &
                    +work*mesh%Ne(i,j))! *(4d0*SiDensity-Ne(i,j))/(4d0*SiDensity) !use Old Ne here!

        SourceUe(i,j)= ((source%E-Egap(i,j))*OnePhotonIonizationRate0*intensity(i,j) &
                     + M_HALF*(2d0*source%E - Egap(i,j))*TwoPhotonIonizationRate0*Int2 )*source%inv_E*((me)/(me+mh))&
                     - Egap(i,j)*work*mesh%Ne(i,j) &
                     + absorptionDrudeE(i,j)*intensity(i,j) &
                     + Egap(i,j)*(matter%AugerRateE*mesh%Nh(i,j) * mesh%Ne(i,j)**2)

        !SourceE(i,j) = SourceE(i,j) - diffNe(i,j)*(1.5d0*kb*Te(i,j))*(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))
        SourceE(i,j) = SourceUe(i,j) - mesh%Te(i,j) * (Ce(i,j)-CeOld(i,j))/dt

        !LossesE(i,j)=AugerRateE * (mesh%Ne(i,j))**2d0 * mesh%Nh(i,j) + AugerRateH * (mesh%Nh(i,j))**2d0 * mesh%Ne(i,j) !use Old Ne, Nh here!
        LossesE(i,j)=mesh%Ne(i,j) * mesh%Nh(i,j) * ( matter%AugerRateE * mesh%Ne(i,j) + matter%AugerRateH * mesh%Nh(i,j) ) !This is more perfomant like that

        work = ImpactIonizationRate(mesh%Th(i,j),EgapH(i,j), ImpactOff)

        GainsH(i,j)=(OnePhotonIonizationRate0*intensity(i,j)*source%inv_E &
                    +M_HALF*TwoPhotonIonizationRate0*Int2*source%inv_E &
                    +work*mesh%Nh(i,j)) !*(4d0*SiDensity-Ne(i,j))/(4d0*SiDensity) !use Old Nh here

        SourceUh(i,j)=((source%E-Egap(i,j))* OnePhotonIonizationRate0*intensity(i,j) &
                     + M_HALF*(M_TWO*source%E - Egap(i,j))*TwoPhotonIonizationRate0*Int2)*source%inv_E * ((me)/(me+mh))  &
                     - Egap(i,j)*work*mesh%Nh(i,j) &
                     + absorptionDrudeH(i,j)*intensity(i,j) &
                     + Egap(i,j)*(matter%AugerRateH*mesh%Ne(i,j) * mesh%Nh(i,j)**2)

        !SourceH(i,j) = SourceH(i,j) - diffNh(i,j)*(1.5d0*kb*Th(i,j)*(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j))))
        SourceH(i,j) = SourceUh(i,j) - mesh%Th(i,j) * (Ch(i,j)-ChOld(i,j))/dt

        LossesH(i,j)=LossesE(i,j)

        ! first order precision
!         SourceS(i,j) = 0d0 - Ts(i,j)*(Cs(i,j)-CsOld(i,j))/Cs(i,j)

        ! second order precisino
!         SourceS(i,j) = 0d0 - (3d0*Cs(i,j)-4d0*CsOld(i,j)+CsPrev(i,j))/(2d0*dt) * Ts(i,j)

        ! third order precision in already included in the scheme (Maple generated since complexity increases substancially)

        !TODO: to be implemented
        ! VeX(i,j)=0d0
        ! VeY(i,j)=0d0
        ! VhX(i,j)=0d0
        ! VhY(i,j)=0d0
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call profiler_stop(prof)
   end subroutine ComputeGainsAndLosses

!------------------------------------------------------------------
    !TODO: This is not a material related property, this should not be in this file
    subroutine TabCreateFL(FermiMaxLines, FermiTableE, FermiTableH)
      implicit none
      integer, intent(in)    :: FermiMaxLines
      real(8) :: FermiTableE(1:9,1:FermiMaxLines), FermiTableH(1:9,1:FermiMaxLines)

      integer :: unit1, unit2
      unit1=15; unit2=16
      open (unit1,file='FermiDatasE.dat')
      open (unit2,file='FermiDatasH.dat')
      read (unit1,*) FermiTableE(1:9,1:FermiMaxLines) !, FermiTableE(2,:) !, FermiTableE(:,3), FermiTableE(:,4), &
!             FermiTableE(:,5), FermiTableE(:,6), FermiTableE(:,7), FermiTableE(:,8), &
!             FermiTableE(:,9)
      read (unit2,*) FermiTableH(1:9,1:FermiMaxLines) !1), FermiTableH(:,2), FermiTableH(:,3), FermiTableH(:,4), &
!              FermiTableH(:,5), FermiTableH(:,6), FermiTableH(:,7), FermiTableH(:,8), &
!             FermiTableH(:,9)
! 222        format (1F10.2, 3x, 1F10.2, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x)
      close(unit1); close(unit2)

      !    fi_min=0.012d0
      !    fistep=0.2d0
      !    fi_max=fi_num*fistep
      write(*,*) FermiTableE(3,463), FermiTableH(3,450)
! stop
    end subroutine TabCreateFL

    !------------------------------------------------------------------
    !TODO: Create a batch version of this routine
    !TODO: This is not a material related property, this should not be in this file
    integer(8) function FermiIndex(NeNc, FermiMaxLines)
      use Maths_m
      implicit none
      ! Input: Value of density/DOS
      ! returns the index to take in the Fermi files
      real(8) NeNc, NeNc0, dNeNc
      integer(8), intent(in) :: FermiMaxLines
      NeNc0=1d-38
      dNeNc=1.03d0 !NeNc=NeNc0*dNeNc**n

      FermiIndex=nint(log10(NeNc/NeNc0)/log10(dNeNc)+M_ONE)
      if(FermiIndex < 1 .OR. FermiIndex > FermiMaxLines) then
        write(*,*) "FermiIndex problem: NeNc=", NeNc, "FermiIndex=", FermiIndex
      end if
    end function FermiIndex
