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
!> @file main_explicit_Mie.f90
!
! DESCRIPTION:
!> @brief **** calculate temperature distribution in 2D in a tip ****
!
!> @author
!> Thibault J.Y. Derrien
!> Laboratoire Hubert Curien, UMR CNRS, St-Etienne
!> ANR Ultrasonde
!
!> @date
!> Jul-Dec 2012 - Initial Version
!------------------------------------------------------------------------------

program Flaps
USE OMP_LIB
! include 'Bivariate.f'
! USE Bivariate
USE libmsh2vf !Script provided by A. Mouton, Univ Lille1, France for GMSH interfacing
!   use control_file !Script provided by Jason Blevins, Ohio State University

use Maths_m
use Mie_m
use Output_m
use Types_m

implicit none

    type(MeshValues) :: mesh, dual, newmesh
    type(LaserParams):: laser
    type(InputParameters) :: Params

    real(8), parameter::Tout=80d0 ,&  !external temperature (K)
                        potential0=7d3,&         ! potential at the bottom of the needle ; default = 7d3
                        potentialNull=0d0 !, &
 !                       phiMie0=1d0*acos(-1d0)                ! Mie scattering: plane angle in cylindrical coordinates
    
    real(8), parameter:: dt0=1d-18,& !time step (s)
                        tmax=0d-15,& !stop time
                        coeffDilaDt=2d0        ,& !diltation coeff before dt change
                        xmin=-10d-6       ,& !mesh min
                        xmax=10d-6       ,& !mesh max
                        ymin=-10d-6       ,&                
                        ymax=10d-6       ,&
                        tCenter=0d0             !time of gaussian intensity maximum
     real(8)               tmin                !max absolute time
    
                        
    integer(8), parameter::  iterOut=1000       ,& ! number of iterations between each stdout
                        iterOutMaps=1000      ,& ! number of outputs for maps between each stdout
                          VirtualPoints=3, & !number of virtual points to exclude from the GMSH file (locate them at the beginning!)
                          Mv=101       ,& !number of celles in the Vessel domain (larger) X direction
                          Nv=101        ,& !number of celles in the Vessel domain (larger) Y direction
                          MeshChoice=1       ,& !0: rectangle (xmin,xmax)(ymin,ymax). 1: Experimental cones, 2: Cone in a vessel (HS), 3: import GMSH (working)
                          MeshIterations=500000        ,&        !number of iterations to calculate meshNeedle
                          MeshIterationsVessel=100*Mv,&        !number of iterations to calculate meshVessel
                          MeshShift=1       ,&         !number of cells x N in the tip, 343 nm: 2; 515 nm: 3;
                          FermiMaxLines=3584        ,&        ! >= number of lines in Fermi file
                          AdaptativeTimeStep=1
       !                   SORiterations=1        ,&        !iteration number for over-relaxation method
       !                   InterpolateMethod=1        ,&        ! 0: linear, 1: bicubic
       !                   UseInterpolation=0        ,&
       !                   SolveImplicit=1, &        ! 0: use explicit schemes, 1: use implicit scheme (band diagonal matrixes)
       !                   numberOfNeighbours=4
                          
                          
        
    real(8), parameter::   NeedleAngleDeg=1d0        ,& ! deg
                          NeedleRadius=4d-9        ,& !m
                          NeedleLength=3d-6        , &        !m
                          SORcoeff=1.2d0        ,&        ! near 1
                          MeshDensity0=5d5        ,& ! amplitude of source for mesh refinement
                          MeshConvergenceEpsilon=1d-10, & !error tolerance on meshing convergency
                          MeshDamping0=1d8, &         ! damping coefficient (m^-1) for mesh refinement
                          ActivateInduction=0d0, &
                          CrossCoeff=-1d0, &                ! 0d0: OFF, 1d0: ON
                          maxCFL=1d-3                        ! maximum admitted on CFL condition for any time step increase
                          
    integer(8), parameter::  DrudeHeating=1, &        ! free-carrier absorption, 0: Drude heating OFF, 1: enabled (1-epsDrude)
                            ConductivityFix=-1, &        ! 2: Consider ambipolar diffusion in equations (but careful with boundary conditions)
                                                ! 1: consider Tritt particle transport (great expression), but Dumber field is needed !!! -> Poisson ! 
                                                ! 0: only fourier conductivity
                                                !-1: diffusion and conductivity OFF
                            CouplingDebug=0        ,&        !0: e/h - lattice coupling enabled, 1: disabled
                            AugerOff=0       ,&
                            ImpactOff=0       ,&
                            ConvectionEnergy=0        ,&         !0: work with Te, no convection. 1: work with Ue, convection
                            DisableCrossDiffusion=0, &
                            PoissonOn=0       ,& !0: Poisson solver is OFF. 1: Calculation of potential ON. 
                            DriftOn=0       ,& !0: Drift is disabled. 1: Enabled. 
                            CathodeZone=1        ,& !1: on the needle bottom, 0: on back vessel (not physical but stable)
                            ShiftFixedPotential=-1, &        ! while CathodeZone=1, use to adjust the number of points on which tension is applied
                            PoissonSolver=0        ,& !0: Full matrix inversion once, 1: SOR iterative for each dt
                            InterpolateOff=0,         &        !just to test speedup...
                            BandBendingInFDTD=0        ,&        !use the interpolation of FDTD 1030 nm with band-bending contribution
!                            PolarizationSource=0, &        ! 0: source TE, 1: source TM
                            !TODO: MieScattering=1 crashed!
                            NewtonIterations=1000, &
                            ExpNeedleType=0
        
    real(8), parameter::  epsilonStatic0=11.66570433d0 !,0.01404457712d0)                ! dielectric constant for static field

    real(8):: me       ,&    ! electron effective mass for conductivity !0.24 (source ?)
              mh       ,&    ! hole effective mass for conductivity !0.81 (source ?)
              meDOS       ,& ! electron effective mass for DOS
              mhDOS          ! hole effective mass for DOS

    real(8), parameter:: Ne0=1d5                                , &
                         Nh0=1d5, &                                !initial density (to calculate auto using fermi!) : at 80 K, Ne=1d5
                         Nlimit=1d0       ,&! lowest possible density
                         Nborder=1d23       ,&! density on boundaries to consider defect layer
                         DefectThickness=1d-7
!                      epsilonStatic0=(11.66570433d0,0.01404457712d0)                ! dielectric constant for static field
    
    integer(8)         nbiter, i, j, k, nmax, NeedleIndexX, NeedleIndexY, maxFermiIndexE, maxFermiIndexH, &
                Mp, Np, RunningIndex
    real(8)         t, t0, dx, dy, x0, y0, dt, dt2, dt3, dt4, h1, h2, h3
    real(8)         Te0, Th0, I0 !initial values of the problem

    real(8)     nuColl, &!        total collision frequency
                etae, etah, &        ! reduced chemical Fermi potential
                work

    real(8), allocatable, dimension(:,:) :: &
                Ue, & !electron energy
                Uh, & !hole energy
                UeNew, & !electron energy
                UhNew, & !hole energy
                TsOld, & !lattice temperature (time n-1)
                TsPrev, & !lattice temperature (time n-2)
                GradNeX,& !Grad(Ne)_x
                GradNeY,& !Grad(Ne)_y
                intensity, & !propagated intensity
                intensityDual, & !intensity dual, just for test of the function
                reflectivity, & ! surface reflectivity
                absorptionDrudeE, & ! absorption coefficient
                absorptionDrudeH, & ! absorption coefficient
                diffusionE, & !fick diffusion coefficient for electrons
                diffusionH, & !fick diffusion coefficient for holes
                GainsE, GainsH, &
                LossesE, LossesH, &
                kappae, kappah, kappas, &
                Ce, Ch, Cs, &
                invCe, invCh, invCs, &
                CeOld, ChOld, CsOld, &
                CsPrev, CsPrev2, &
                CouplingE, CouplingH, &
                mobilityE, mobilityH, &
                Egap, &                        ! local gap value
                SourceE, SourceH, & ! heating sources
                SourceUe, SourceUh, & ! free carrier thermal energy sources
                diffNe, diffNh, &         ! just for derivation in time
                x, y, &                 ! needle position indexes
                xNew, yNew, &           ! used to converge the mesh parallely
                xV, yV, &                 ! vessel position indexes
                xDualSW, yDualSW, &                 ! dual mesh position
                xDualSE, yDualSE, &
                xDualNE, yDualNE, &
                xDualNW, yDualNW, &
                xDual, yDual, &
                CFLxT, CFLyT, CFLxN, CFLyN, CFLxTs, CFLyTs, &
                ThermalEnergy, LaserEnergy, &
                TotalElectrons, TotalHoles, &
                DOSe, DOSh, &
                FermiRatioE, FermiRatioH, &
                JeX, JeY, & ! drift vectors for particles
                JhX, JhY, & ! drift vectors for particles
                VeX, VeY, &
                VhX, VhY, &
                Ex, Ey, &        ! fields in the main domain
                potentialNeedle, &        ! potential in the needle
                epsilonNeedle, &        ! dielectric static in the needle
                MaxHeating, &
                MaxHeatingTime, &
                CellAreaN, CellAreaS, &                        ! area of the finite elements
                CellAreaE, CellAreaW, &
                CellVol, InvCellVol,&          !volume of needle mesh cells
                DistW, DistE, &                 ! distance to the center of neighboor cells
                DistN, DistS, &
                DistDualW, DistDualE, &                 ! distance of the element side (equal to area in 2D)
                DistDualN, DistDualS, &
                ShapeFactorNormalE, ShapeFactorNormalW, &
                ShapeFactorNormalS, ShapeFactorNormalN, &
                ShapeFactorTangentE, ShapeFactorTangentW, &
                ShapeFactorTangentS, ShapeFactorTangentN, &
                EintFieldR, EintFieldI, &
                phiMie, &
                Radius

    type(VectorField) :: NormalN, NormalS, NormalW, NormalE ! normal to quadrangle elements
                
                 !TODO: Use dimension
    real(8), allocatable :: CurviWx(:,:), CurviWy(:,:), &                 ! Unit vector between cell centers
                            CurviEx(:,:), CurviEy(:,:), &
                            CurviNx(:,:), CurviNy(:,:), &
                            CurviSx(:,:), CurviSy(:,:),  &
                            TangentWx(:,:), TangentWy(:,:), &                ! Tangent to quadrangle elements
                            TangentEx(:,:), TangentEy(:,:), &
                            TangentNx(:,:), TangentNy(:,:), &
                            TangentSx(:,:), TangentSy(:,:)

    integer(8), allocatable, dimension(:,:) :: &
                FermiIndexE, FermiIndexH, &
                MeshVertice ! data from the GMSH file
    
    real(8) phiMie0
    real(8) Int2
    integer(8) PolarizationSource             ! Value of the Mie angle that will be distributed on various processors
    
    real(8), allocatable, target :: FermiTableE(:,:),&
                                     FermiTableH(:,:),& !reduced Fermi level for electrons and holes
                                     Bvector(:),         &
                                     Xvector(:), XvectorPrev(:),         &
                                     spectralNorm(:),                         & !objects for matrix inversion calculation
                                     xP(:,:), yP(:,:),                        &                 ! vessel position indexes
                                     ExPoisson(:, :), EyPoisson(:, :),                 & ! electric field in the vessel
                                     potential(:, :),                         &        ! electric potential in the vessel
                                     DielectricStatic(:, :),                 & ! dielectric constant for the static field
                                     DummyVessel(:,:),                        & ! option argument for interpolation
                                     NeP(:,:), NhP(:,:)        , &                  !interpolated Ne,Nh in the vessel
                                     CellVolume(:,:), &                        ! volume of the vessel/poisson mesh elements
                                     NormalWxP(:,:), NormalWyP(:,:), &
                                     NormalExP(:,:), NormalEyP(:,:), &                ! normal to quadrangle elements
                                     NormalNxP(:,:), NormalNyP(:,:), &
                                     NormalSxP(:,:), NormalSyP(:,:), &
                                     CellAreaNP(:,:), CellAreaSP(:,:), &
                                     CellAreaEP(:,:), CellAreaWP(:,:)

    integer(8), allocatable, target:: FixedPotentialIndex(:,:)         !array of points where potential has been fixed

! mesh interface with gmsh management
    CHARACTER(LEN=60)               :: namefile_msh
!
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: vertices
    INTEGER, DIMENSION(:,:), POINTER          :: points, segments, triangles, quadrangles, boundedges, edges
    INTEGER, DIMENSION(:), POINTER            :: dim_physical_entities, id_physical_entities, idvertices
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: name_physical_entities
!     INTEGER                                   :: kmouton
    INTEGER                           :: nb_vertices, nb_triangles, nb_quadrangles, nb_edges, nb_boundedges
    
    complex(8), allocatable, dimension(:,:) :: &
                Dielectric, &! solid dielectric function under laser illumination
                DielectricDrudeE, & ! Drude part of dielectric function under laser illumination
                DielectricDrudeH, &
                EintField, EintField2        !Ez internal field for Mie scattering theory

    real(8), allocatable, dimension(:,:) :: &
                OpticalIndex, &  ! solid optical index under laser illumination
                OpticalDamping   ! solid optical damping under laser illumination
                
    complex(8) epsilonInf !, SORsum !material constant

   
    real(8) AugerRateE, AugerRateH, & ! For performances
            sigmaTau, sigmaX, sigmaY, &
            maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN, maxCFLxTs, maxCFLyTs, &
            maxTe, minTe, maxTh, minTh, maxTs, minTs, maxIntensity, maxNe, minNe, maxNh, minNh, &
            maxSourceE, maxGainsE, maxSourceH, maxGainsH, maxGap, maxDiffNe, maxDiffNh, &
            TotalLaserEnergy, TotalThermalEnergy, ElectronPotentialEnergy, ElectronKineticEnergy, &
            NeedleHeight, NeedleA, NeedleB, Needlet0Limit, NeedleXParam, NeedleYParam, NeedleAngle, &
            localT, P2critic, TotalNumOfE, TotalNumOfH, &
            xmin2, xmax2, ymin2, ymax2, &
            MeshConvergence, MeshConvergenceOld
    real :: cpuefficiency, cpu_timestep_duration, ElapsedTime
            
    real(8) sigmaX1, sigmaX2, sigmaX3, sigmaX4, sigmaX5, sigmaX6, sigmaX7, sigmaX8, sigmaX9, &
            sigmaY1, sigmaY2, sigmaY3, sigmaY4, sigmaY5, sigmaY6, sigmaY7, sigmaY8, sigmaY9, &
            spotX1, spotX2, spotX3, spotX4, spotX5, spotX6, spotX7, spotX8, spotX9, &
            spotY1, spotY2, spotY3, spotY4, spotY5, spotY6, spotY7, spotY8, spotY9, &
            x1, x2, x3, x4, x5, x6, x7, x8, x9, &
            y1, y2, y3, y4, y5, y6, y7, y8, y9, &
            I1, I2, I3, I4, I5, I6, I7, I8, I9, &
            periodX, periodY, &
            NeTotal, NhTotal, &
            meshParameterTmax, meshParameterTmin, &
            meshStepDt, localTmin, localTmax, localdT, &
            IntensityEnergy, LaserIntensityEnergy, ElectronEnergy, HoleEnergy, LatticeEnergy, &
            TotalMeshVolume, &
            OnePhotonIonizationRate0, TwoPhotonIonizationRate0
            

!    integer(4) unit1, unit2
    integer(8) ColFermiNeNc, ColFermiEta, ColFermi0, ColFermi1, ColFermi2, &
               ColFermiHalf, ColFermiThreeHalf, ColFermiMenusHalf

    !! FUNCTIONS CALLS
     real(8) ConeExp1Radius, ConeExp2Radius !, Interpolate
     real(8) ConeExp1, ConeExp2, DensityOfState, EgapValue, TwoPhotonIonizationRate, OnePhotonIonizationRate, &
             CollisionFrequency, LatticeHeatCapacity, ImpactIonizationRate
     complex(8) DielectricConstant ! DielectricFunction, DielectricFunctionDrude,

            
!OPENMP declarations
    integer :: myid, nthreads
!#IFDEF (OMP_NUM_THREADS)
!    integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!#ENDIF

  !TODO: move to material.f90
  !TODO: Should be a parameter, to guaranty no modification
  me=0.5d0*me0       ! electron effective mass for conductivity !0.24 (source ?)
  mh=0.5d0*me0       ! hole effective mass for conductivity !0.81 (source ?)
  meDOS=0.36d0*me0   ! electron effective mass for DOS
  mhDOS=0.81d0*me0   ! hole effective mass for DOS


! call omp_set_num_threads(16)
    
   !  !!********OpenMP Test*************
  myid=1 !if openMP is off, then test is disabled; else: will be set to 0 by OpenMP
  nthreads=1 !idem

  !$OMP PARALLEL default(none) private(myid) &
  !$OMP shared(nthreads)
  ! Determine the number of threads and their id
        myid = OMP_GET_THREAD_NUM()
        PRINT *, 'Hello from thread =', myid !TODO: Replace the print*,... by stdout, in the perspective of MPI
        nthreads = OMP_GET_NUM_THREADS()
  !$OMP BARRIER

  !TODO: Replace the write(*,... by stdout, in the perspective of MPI
  if (myid==0) then 
    write(*,'(a)') 'OpenMP TEST'
    write(*,'(a,i1)') 'Number of Threads = ', nthreads
    write(*,'(a)') '*******************************'
  end if 
  !$OMP END PARALLEL
! !!******* END OpenMP test

  call TimerInit()

  call InitInputParameter( Params )
  call LoadInputParameters( "flaps.in", Params )
  call CheckValidityInputParameters( Params )

  allocate(Ue(1:Params%M, 1:Params%N))
  allocate(Uh(1:Params%M, 1:Params%N), & !hole energy
                 UeNew(1:Params%M, 1:Params%N), & !electron energy
                         UhNew(1:Params%M, 1:Params%N), & !hole energy
                TsOld(1:Params%M, 1:Params%N), & !lattice temperature (time n-1)
                TsPrev(1:Params%M,1:Params%N), & !lattice temperature (time n-2)
                GradNeX(1:Params%M, 1:Params%N),& !Grad(Ne)_x
                GradNeY(1:Params%M, 1:Params%N),& !Grad(Ne)_y
                intensity(1:Params%M, 1:Params%N), & !propagated intensity
                intensityDual(1:(Params%M-1),1:(Params%N-1)), & !intensity dual, just for test of the function
                reflectivity(1:Params%M,1:Params%N), & ! surface reflectivity
                absorptionDrudeE(1:Params%M, 1:Params%N), & ! absorption coefficient
                absorptionDrudeH(1:Params%M, 1:Params%N), & ! absorption coefficient
                diffusionE(1:Params%M, 1:Params%N), & !fick diffusion coefficient for electrons
                diffusionH(1:Params%M,1:Params%N), & !fick diffusion coefficient for holes
                GainsE(1:Params%M, 1:Params%N), GainsH(1:Params%M, 1:Params%N), &
                LossesE(1:Params%M, 1:Params%N), LossesH(1:Params%M, 1:Params%N), &
                kappae(1:Params%M, 1:Params%N), kappah(1:Params%M, 1:Params%N), kappas(1:Params%M,1:Params%N), &
                Ce(1:Params%M, 1:Params%N), Ch(1:Params%M, 1:Params%N), Cs(1:Params%M, 1:Params%N), &
                invCe(1:Params%M, 1:Params%N), invCh(1:Params%M, 1:Params%N), invCs(1:Params%M, 1:Params%N), &
                CeOld(1:Params%M, 1:Params%N), ChOld(1:Params%M,1:Params%N), CsOld(1:Params%M, 1:Params%N), &
                CsPrev(1:Params%M, 1:Params%N), CsPrev2(1:Params%M, 1:Params%N), &
                CouplingE(1:Params%M, 1:Params%N), CouplingH(1:Params%M, 1:Params%N), &
                mobilityE(1:Params%M, 1:Params%N), mobilityH(1:Params%M, 1:Params%N), &
                Egap(1:Params%M, 1:Params%N), &                        ! local gap value
                SourceE(1:Params%M, 1:Params%N), SourceH(1:Params%M, 1:Params%N), & ! heating sources
                SourceUe(1:Params%M, 1:Params%N), SourceUh(1:Params%M, 1:Params%N), & ! free carrier thermal energy sources
                diffNe(1:Params%M, 1:Params%N), diffNh(1:Params%M, 1:Params%N), &         ! just for derivation in time
                x(1:Params%M, 1:Params%N), y(1:Params%M, 1:Params%N), &                 ! needle position indexes
                xNew(1:Params%M, 1:Params%N), yNew(1:Params%M, 1:Params%N), &           ! used to converge the mesh parallely
                xV(1:Mv, 1:Nv), yV(1:Mv, 1:Nv), &                 ! vessel position indexes
                xDualSW(1:Params%M, 1:Params%N), yDualSW(1:Params%M, 1:Params%N), &                 ! dual mesh position
                xDualSE(1:Params%M, 1:Params%N), yDualSE(1:Params%M, 1:Params%N), &
                xDualNE(1:Params%M, 1:Params%N), yDualNE(1:Params%M, 1:Params%N), &
                xDualNW(1:Params%M, 1:Params%N), yDualNW(1:Params%M, 1:Params%N), &
                xDual(1:(Params%M-1), 1:(Params%N-1)), yDual(1:(Params%M-1), 1:(Params%N-1)), &
                CFLxT(1:Params%M, 1:Params%N), CFLyT(1:Params%M, 1:Params%N), &
                CFLxN(1:Params%M, 1:Params%N), CFLyN(1:Params%M, 1:Params%N), &
                CFLxTs(1:Params%M,1:Params%N), CFLyTs(1:Params%M,1:Params%N), &
                ThermalEnergy(1:Params%M, 1:Params%N), LaserEnergy(1:Params%M, 1:Params%N), &
                TotalElectrons(1:Params%M, 1:Params%N), TotalHoles(1:Params%M, 1:Params%N), &
                DOSe(1:Params%M, 1:Params%N), DOSh(1:Params%M, 1:Params%N), &
                FermiRatioE(1:Params%M,1:Params%N), FermiRatioH(1:Params%M,1:Params%N), &
                JeX(1:Params%M, 1:Params%N), JeY(1:Params%M, 1:Params%N), & ! drift vectors for particles
                JhX(1:Params%M, 1:Params%N), JhY(1:Params%M, 1:Params%N), & ! drift vectors for particles
                VeX(1:Params%M, 1:Params%N), VeY(1:Params%M, 1:Params%N), &
                VhX(1:Params%M, 1:Params%N), VhY(1:Params%M, 1:Params%N), &
                Ex(1:Params%M,1:Params%N), Ey(1:Params%M,1:Params%N), &        ! fields in the main domain
                potentialNeedle(1:Params%M,1:Params%N), &        ! potential in the needle
                epsilonNeedle(1:Params%M,1:Params%N), &        ! dielectric static in the needle
                MaxHeating(1:Params%M, 1:Params%N), &
                MaxHeatingTime(1:Params%M, 1:Params%N), &
                CellAreaN(1:Params%M, 1:Params%N), CellAreaS(1:Params%M, 1:Params%N), &                        ! area of the finite elements
                CellAreaE(1:Params%M, 1:Params%N), CellAreaW(1:Params%M, 1:Params%N), &
                CellVol(1:Params%M, 1:Params%N), InvCellVol(1:Params%M, 1:Params%N),&          !volume of needle mesh cells
                DistW(1:Params%M,1:Params%N), DistE(1:Params%M,1:Params%N), &                 ! distance to the center of neighboor cells
                DistN(1:Params%M,1:Params%N), DistS(1:Params%M,1:Params%N), &
                DistDualW(1:Params%M,1:Params%N), DistDualE(1:Params%M,1:Params%N), &                 ! distance of the element side (equal to area in 2D)
                DistDualN(1:Params%M,1:Params%N), DistDualS(1:Params%M,1:Params%N), &
                ShapeFactorNormalE(1:Params%M, 1:Params%N), ShapeFactorNormalW(1:Params%M, 1:Params%N), &
                ShapeFactorNormalS(1:Params%M, 1:Params%N), ShapeFactorNormalN(1:Params%M, 1:Params%N), &
                ShapeFactorTangentE(1:Params%M, 1:Params%N), ShapeFactorTangentW(1:Params%M, 1:Params%N), &
                ShapeFactorTangentS(1:Params%M, 1:Params%N), ShapeFactorTangentN(1:Params%M, 1:Params%N), &
                EintFieldR(1:Params%M,1:Params%N), EintFieldI(1:Params%M, 1:Params%N), &
                phiMie(1:Params%M, 1:Params%N), &
                Radius(1:Params%M, 1:Params%N))
 allocate(  FermiIndexE(1:Params%M,1:Params%N), FermiIndexH(1:Params%M,1:Params%N), &
                MeshVertice(1:Params%M, 1:Params%N) )
 allocate(  Dielectric(1:Params%M,1:Params%N), &! solid dielectric function under laser illumination
            OpticalIndex(1:Params%M,1:Params%N), &
            OpticalDamping(1:Params%M,1:Params%N), &
                DielectricDrudeE(1:Params%M,1:Params%N), & ! Drude part of dielectric function under laser illumination
                DielectricDrudeH(1:Params%M,1:Params%N), &
                EintField(1:Params%M,1:Params%N), EintField2(1:Params%M,1:Params%N))

!******** READ PARAMETER INPUT FILE ***********
CALL control_file(phiMie0, PolarizationSource) !read Miescattering parameters into external file
write(*,*) "Importing data on Polarization."
!******** READ GMSH MESH FILE ************
namefile_msh='libs/gmsh/mesh.msh'
RunningIndex=1 !gonna be used to mesh down
! CALL extract_parameters(namefile_msh, namefile_vf)
! We read the .msh file

WRITE(*,*) 'Loading GMSH mesh...'

CALL read_msh_file(namefile_msh, vertices, points, segments, triangles, quadrangles, dim_physical_entities, &
& id_physical_entities, name_physical_entities, idvertices)
! We verify if the triangles are sorted in trigonometric sense and we bring correction if necessary
IF (associated(triangles)) THEN
    CALL correct_orientation(vertices, triangles)
END IF
! We verify if the quadrangles are sorted in trigonometric sense and we bring correction if necessary
IF (associated(quadrangles)) THEN
    CALL correct_orientation(vertices, quadrangles)
END IF

nb_vertices = size(vertices,1)
IF (associated(triangles)) THEN
    nb_triangles = size(triangles,1)
ELSE
    nb_triangles = 0
END IF
IF (associated(quadrangles)) THEN
    nb_quadrangles = size(quadrangles,1)
ELSE
    nb_quadrangles = 0
END IF
nb_edges = size(edges,1)
nb_boundedges = size(boundedges,1)

! We build the edges
CALL compute_edges(triangles, quadrangles, edges)
! For each edge, we identify the physical zone in which it is included
CALL identify_physical_zone_for_edges(edges, segments, triangles, quadrangles)
! We extract the bound edges
CALL extract_boundedges(edges, boundedges)
! We are mainly interested in getting the point data.
    
WRITE(*,FMT=*) "MESH INFORMATIONS"
WRITE(*,FMT=*) "[GMSH] Quadrangles number"
WRITE(*,FMT=*) nb_quadrangles
WRITE(*,FMT=*) "[GMSH] Node number"
WRITE(*,FMT=*) nb_vertices

!***** Attribute positions of the nodes by recursivity

RunningIndex=VirtualPoints

! Corners
RunningIndex=RunningIndex+1; MeshVertice(1,1) = RunningIndex
RunningIndex=RunningIndex+1; MeshVertice(Params%M,1) = RunningIndex
RunningIndex=RunningIndex+1; MeshVertice(Params%M,Params%N) = RunningIndex
RunningIndex=RunningIndex+1; MeshVertice(1,Params%N) = RunningIndex

!North
do i=2, Params%M-1
  RunningIndex=RunningIndex+1; MeshVertice(i, 1)=RunningIndex
end do
!East
do j=2, Params%N-1
  RunningIndex=RunningIndex+1; MeshVertice(Params%M, j)=RunningIndex
end do
!South
do i=Params%M-1,2, -1
  RunningIndex=RunningIndex+1; MeshVertice(i,Params%N)=RunningIndex
end do
!West
do j=Params%N-1, 2, -1
  RunningIndex=RunningIndex+1; MeshVertice(1,j)=RunningIndex
end do

!rest of the domain
do i=2,Params%M-1
  do j=2, Params%N-1
    RunningIndex=RunningIndex+1; MeshVertice(i,j)=RunningIndex
  end do
end do


WRITE(*,*) 'Latest running index while remeshing', RunningIndex


!***** Compute geometrical data


!**** INITIALIZATION

  call init_laser(laser)

  tmin=tCenter-5d0*laser%tau


  dt=dt0
  dt2=dt0
  dt3=dt0
  dt4=dt0

  ColFermiNeNc=2; ColFermiEta=3; ColFermi0=4; ColFermi1=5; ColFermi2=6; ColFermiHalf=7; 
  ColFermiThreeHalf=8; ColFermiMenusHalf=9;

! test field
!TODO: move to material.f90
  if(AugerOff.eq.0) then
    AugerRateE=2.3d-43
    AugerRateH=7.8d-44
  else
    AugerRateE=0d0
    AugerRateH=0d0
  end if

    !TODO: Move to LaserParams
    sigmaTau=laser%tau/(2d0*sqrt2ln2)
    sigmaX=laser%spotX/(2d0*sqrt2ln2)
    sigmaY=laser%spotY/(2d0*sqrt2ln2)

    !TODO: Move to LaserParams
    I0=laser%fluence/laser%tau * sqrt(4d0 * log(2d0) / pi)
    
    if(laser%lambda.eq.515d-9) then
      if(PolarizationSource.eq.0) then
        x1=1.5d-7; y1=0.d0; I1=0d0*I0; spotX1=100d-9; spotY1=50d-9; !
        x2=3.3d-7; y2=0.d0; I2=0d0*9d0*I0; spotX2=50d-9; spotY2=50d-9; !3.53W
        x3=6.5d-7; y3=-3d-8; I3=0d0*17d0*I0; spotX3=70d-9; spotY3=150d-9; !28W
        x4=6.5d-7; y4=3d-8; I4=0d0*17d0*I0; spotX4=70d-9; spotY4=150d-9; !28W
        x5=9.5d-7; y5=0d0; I5=0d0*7d0*I0; spotX5=50d-9; spotY5=50d-9; !2.74W
        x6=1.2d-6; y6=0d0; I6=0d0*8d0*I0; spotX6=50d-9; spotY6=50d-9; !3.14W
        x7=1.37d-6; y7=0d0; I7=0d0*9d0*I0; spotX7=50d-9; spotY7=50d-9; !3.53W
        x8=6.2d-7; y8=-40d-9; I8=0d0*0d0*I0; spotX8=250d-9; spotY8=50d-9;
        x9=1.3d-6; y9=-20d-9; I9=0d0*0d0*I0; spotX9=500d-9; spotY9=100d-9;
      else
        periodX=200d-9; periodY=50d-9
        x1=1.0d-7; y1=0d-9; I1=13d0*I0; spotX1=120d-9; spotY1=70d-9; !17.15W
        x2=2.7d-7; y2=10d-9; I2=10d0*I0; spotX2=70d-9; spotY2=60d-9; !6.6 W
        x3=4.5d-7; y3=30d-9; I3=12d0*I0; spotX3=70d-9; spotY3=70d-9; !9.2 W unstable
        x4=6.2d-7; y4=40d-9; I4=28d0*I0; spotX4=70d-9; spotY4=120d-9; !36 W unstable
        x5=9.0d-7; y5=45d-9; I5=7d0*I0; spotX5=50d-9; spotY5=100d-9; !13.35 W unstable
        x6=1.1d-6; y6=60d-9; I6=15d0*I0; spotX6=50d-9; spotY6=70d-9; !8.24 W unstable
        x7=1.30d-6; y7=70d-9; I7=I0; spotX7=50d-9; spotY7=40d-9; !5W
        x8=6.2d-7; y8=-40d-9; I8=I0; spotX8=250d-9; spotY8=50d-9; !30W
        x9=1.3d-6; y9=-20d-9; I9=13d0*I0; spotX9=500d-9; spotY9=100d-9; !102W
      end if
    end if

    !TODO: I think that one invented arrays and loops for handling similar situations ;)
    sigmaX1=spotX1/(2e0*sqrt2ln2)
    sigmaY1=spotY1/(2e0*sqrt2ln2)
    sigmaX2=spotX2/(2e0*sqrt2ln2)
    sigmaY2=spotY2/(2e0*sqrt2ln2)
    sigmaX3=spotX3/(2e0*sqrt2ln2)
    sigmaY3=spotY3/(2e0*sqrt2ln2)
    sigmaX4=spotX4/(2e0*sqrt2ln2)
    sigmaY4=spotY4/(2e0*sqrt2ln2)
    sigmaX5=spotX5/(2e0*sqrt2ln2)
    sigmaY5=spotY5/(2e0*sqrt2ln2)
    sigmaX6=spotX6/(2e0*sqrt2ln2)
    sigmaY6=spotY6/(2e0*sqrt2ln2)
    sigmaX7=spotX7/(2e0*sqrt2ln2)
    sigmaY7=spotY7/(2e0*sqrt2ln2)
    sigmaX8=spotX8/(2e0*sqrt2ln2)
    sigmaY8=spotY8/(2e0*sqrt2ln2)
    sigmaX9=spotX9/(2e0*sqrt2ln2)
    sigmaY9=spotY9/(2e0*sqrt2ln2)
    
    t0=tCenter; x0=laser%xCenter; y0=laser%yCenter;
    
    ! We open the different files
    call InitOutputs()
    
    
  TotalLaserEnergy=0d0; 
  IntensityEnergy=0d0;
  LaserIntensityEnergy=0d0; 
  ElectronKineticEnergy=0d0; ElectronPotentialEnergy=0d0; ElectronEnergy=0d0
  HoleEnergy=0d0
  LatticeEnergy=0d0
  TotalThermalEnergy=0d0; 
  cpuefficiency=0d0; cpu_timestep_duration=0d0
  
  
!***************** MESH GENERATION *****************
  write(*,*) "[Mesh] Building..."


  !Allocate the mesh and the dual mesh
  call initmesh(mesh, Params%M, Params%N)
  call initmesh(dual, Params%M-1, Params%N-1)
  call initmesh(newmesh, Params%M, Params%N)

  ! building rectangular mesh 
  dx=(xmax-xmin)/(Params%M+1)
  dy=(ymax-ymin)/(Params%N+1)
  do i=1,Params%M
    do j=1,Params%N
      x(i,j)=0d0
      y(i,j)=0d0
    end do
  end do
  
  write(*,*)

  if(MeshChoice.eq.0) then  !rectangular mesh as main domain
    do i=1,Params%M
      do j=1,Params%N
        x(i,j)=dx*real(i)+xmin
        y(i,j)=dy*real(j)+ymin
        write(MeshInfo%unit,885, advance='yes') x(i,j), y(i,j), i, j ! writing as main mesh
      end do
      write(MeshInfo%unit,*) " "
    end do
  end if
!   else ! conical mesh as main domain
    
  if(MeshChoice==1 .OR. MeshChoice==2) then !conical mesh as main domain
    !boudary definition
    NeedleIndexX=Params%M
    NeedleIndexY=Params%N
    NeedleAngle=NeedleAngleDeg*pi/180d0
    NeedleA=NeedleRadius/(tan(NeedleAngle/2d0)**2)
    NeedleB=NeedleRadius/tan(NeedleAngle/2d0)
    Needlet0Limit=acos(NeedleRadius/(NeedleLength*tan(real(NeedleAngle)/2d0)**2+NeedleRadius))
    
!building of the conical mesh: write boundaries, then solve laplace, and iterate
    
    NeedleXParam=(Needlet0Limit-0d0)/NeedleIndexX !dt0 for X (0,t0)
    NeedleYParam=(Needlet0Limit-0d0)/NeedleIndexY !dt0 for Y (0,t0)
    
    P2critic=(real(Params%N)-1d0)**2/(2d0*real(j)-real(Params%N)-1)
    !! this limit applies if hyperbolic contour is chosen
!     if(real(MeshShift)**2 > P2critic) then
!       write(*,*) "Dilatation of the apex is too large. Reduce <= ", floor(sqrt(P2critic))
!       stop
!     end if

!! just to check which point has not been defined - every point should be replaced !
    x(:,:)=-1d0
    y(:,:)=-1d0
    
    MeshConvergence=0d0
    
    inner: do k=1,MeshIterations
    !!! HYPERBOLIC CONTOUR
!     !1/ tip apex, 2/ bottom, 3/ upper, 4/ cone base
!         do j=1, N
!         ! t = (-p*dt0,p*dt0), p integer defined by MeshShift
!             localT=real(MeshShift)*real(MeshShift)*Needlet0Limit*( (2d0*real(j-1)) / (real(N)-1d0) - 1d0 )/real(N-1)
!             x(1,j)=NeedleA*(1d0/cos( localT ) - 1d0)
!             y(1,j)=NeedleB*tan(localT)
!         end do
! !         
!         do i=1,M
!             localT=Needlet0Limit*( ( (real(i)-0d0)/(real(M)-1d0) )*(real(MeshShift)*real(MeshShift)/(real(N)-1d0)-1d0)+1d0)
!             x(i,1)=NeedleA*(1d0/cos( - localT )-1d0)
!             y(i,1)=NeedleB*tan( - localT)
!         end do
! ! ! !  
!          do i=1, M
!             !t = (-t0;-p*dt0), p integer        
!             localT=Needlet0Limit*( ( (real(i)-0d0)/(real(M)-1d0) )*(real(MeshShift)*real(MeshShift)/(real(N)-1d0)-1d0)+1d0)
!             x(i,N)=NeedleA*(1d0/cos( localT )-1d0)
!             y(i,N)=NeedleB*tan( localT)
!         end do
!         
!         do j=1, N
!             localT=Needlet0Limit*((2d0*(real(j)-1d0)/(real(N)-1d0))-1d0)
!             x(M,j)=NeedleA*(1d0/cos(Needlet0Limit)-1d0)                !base cone boundary
!             y(M,j)=NeedleB*tan(localT)         !base cone boundary
!         end do
!       !!END HYPERBOLIC CONTOURS

!!!!!!! Experimental cone of the mesh

! contour is defined by a parameter t
        if(ExpNeedleType.eq.1) then
        !  meshParameterTmax=1.4d0 !maximum t parameter
        !  meshParameterTmin=-1.4d0        !minimum t parameter
          meshParameterTmax=0.7866818869d0
          meshParameterTmin=-0.7866818869d0
        else
!           meshParameterTmax=2.471556d0 !maximum t parameter
!           meshParameterTmin=-3.35d0        !minimum t parameter
          ! reduced length cone for accelerated calculations
          meshParameterTmax=0.2361d0 !maximum t parameter
          meshParameterTmin=-0.2531506894d0        !minimum t parameter
        end if
        meshStepDt=(meshParameterTmax-meshParameterTmin)/(real(2*(Params%M-1)+Params%N)) !step of parameter t to define corners of the mesh

        !!! BOTTOM
        localTmin=meshParameterTmin !tmin in the segment (e.g. contour section)
        localTmax=-real(MeshShift)*meshStepDt !tmax in the line (e.g. contour section)
        
        if(localTmax<localTmin) then
          write(*,*) "MeshShift is too large. Reduce.", -localTmax/meshStepDt
          stop
        end if
        
        localdT=(localTmax-localTmin)/real(Params%M-1) !parameter t to distribute the nodes on the segment
        !$OMP PARALLEL DEFAULT(none) SHARED(x, y, Params, localTmin, localdT, localTmax, &
        !$OMP meshStepDt, meshParameterTmax, meshParameterTmin) &
        !$OMP PRIVATE(localT)
        !$OMP DO
        do i=1,Params%M
          localT=localTmax-real(i-1)*localdT
          if(ExpNeedleType.eq.1) then
            x(i,1)=ConeExp2(localT)
          else
            x(i,1)=ConeExp1(localT)
          end if
          y(i,1)=localT
        end do
        !$OMP END DO

        !!! TIP APEX
        localTmin=-real(MeshShift)*meshStepDt !tmin in the segment (e.g. contour section)
        localTmax=real(MeshShift)*meshStepDt !tmax in the line (e.g. contour section)
        
        !$OMP MASTER
        if(localTmax<localTmin) then
          write(*,*) "MeshShift is nul ?! MeshShift=", MeshShift
          stop
        end if
        !$OMP END MASTER
        
        localdT=(localTmax-localTmin)/real(Params%N-1) !parameter t to distribute the nodes on the segment
        !$OMP DO
        do j=2,Params%N
          localT=localTmin+real(j-1)*localdT
          if(ExpNeedleType.eq.1) then
            x(1,j)=ConeExp2(localT)
          else
            x(1,j)=ConeExp1(localT)
          end if
          y(1,j)=localT
        end do
        !$OMP END DO

        !!! UP
        localTmin=real(MeshShift)*meshStepDt !tmin in the segment (e.g. contour section)
        localTmax=meshParameterTmax !tmax in the line (e.g. contour section)
        
        if(localTmax<localTmin) then
          write(*,*) "MeshShift is too large. Reduce.", localTmax/meshStepDt
          stop
        end if
        
        localdT=(localTmax-localTmin)/real(Params%M-1) !parameter t to distribute the nodes on the segment

        if(ExpNeedleType.eq.1) then
          !$OMP DO
          do i=1,Params%M
            localT=localTmin+real(i-1)*localdT
            x(i,Params%N)=ConeExp2(localT)
            y(i,Params%N)=localT
          end do
          !$OMP END DO
        else
          !$OMP DO
          do i=1,Params%M
            localT=localTmin+real(i-1)*localdT
            x(i,Params%N)=ConeExp1(localT)
            y(i,Params%N)=localT
          end do
          !$OMP END DO
        end if


        !!! BACKSIDE
        localTmin=meshParameterTmin !tmin in the segment (e.g. contour section)
        localTmax=meshParameterTmax !tmax in the line (e.g. contour section)
        localdT=(localTmax-localTmin)/real(Params%N-1) !parameter t to distribute the nodes on the segment

        !$OMP DO
        do j=1,Params%N
          localT=localTmin+real(j-1)*localdT
          if(ExpNeedleType.eq.1) then
            x(Params%M,j)=ConeExp2(localTmax)
          else
            x(Params%M,j)=ConeExp1(localTmax)
          end if
          y(Params%M,j)=localT
        end do
       !$OMP END DO
       !$OMP END PARALLEL
        
!         ! sort points before solving
!         do j=1,M
!           do i=1,M
!             if(x(i,1) > x(j,1)) then 
!               call swap(x(i,1), x(j,1))
!               call swap(y(i,1), y(j,1))
!             end if
!             if(x(i,N) > x(i+1,N)) then 
!               call swap(x(i,N), x(j,N))
!               call swap(y(i,N), y(j,N))
!             end if
!           end do
!         end do
        
!         ! resolution of laplace (sequential)
!         do i=2, M-1
!           do j=2,N-1
!             x(i,j)=(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1))/4d0
!             y(i,j)=(y(i+1,j)+y(i-1,j)+y(i,j+1)+y(i,j-1))/4d0
!           end do
!         end do
        

        ! save the contours (can be optimized)
        xNew(:,:)=x(:,:)
        yNew(:,:)=y(:,:)


        MeshConvergenceOld=MeshConvergence
        MeshConvergence=0d0

        !$OMP PARALLEL DEFAULT(none) SHARED(x, y, xNew, yNew, Params, MeshConvergence)
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, Params%N-2)
        do i=2, Params%M-1
          do j=2,Params%N-1
            xNew(i,j)=(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1))*.25d0
            yNew(i,j)=(y(i+1,j)+y(i-1,j)+y(i,j+1)+y(i,j-1))*.25d0
          end do        
        end do
        !$OMP END DO

        !$OMP DO COLLAPSE(2) REDUCTION(+:MeshConvergence)  SCHEDULE(DYNAMIC, Params%N)
        do i=1, Params%M
          do j=1, Params%N
            MeshConvergence=MeshConvergence+(xNew(i,j)**2+yNew(i,j)**2)
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! replace with the new mesh
        x(:,:)=xNew(:,:)
        y(:,:)=yNew(:,:)

        
        if(mod(k,MeshIterations/1000).eq.0) then
          write(*,*) "[Mesh] Convergence (", k, ")=", 1d-6*abs(MeshConvergence-MeshConvergenceOld)
        end if
        
        if(abs(MeshConvergence-MeshConvergenceOld) < MeshConvergenceEpsilon*1d6 ) then
          write(*,*) "[Mesh] Convergence is reached. The show can go on."
          write(*,*) "[Mesh) ITERATIONS=", k
          exit inner
        end if
        
    end do inner !end of iterations for mesh

    x(:,:)=1d-6*x(:,:)
    y(:,:)=1d-6*y(:,:)
    
    !writing of the mesh
    do i=1,Params%M
      do j=1,Params%N
        write(MeshInfo%unit, 885, advance='yes') x(i,j), y(i,j), i, j
        885        FORMAT (1E15.8, 3x, 1E15.8, 3x, I4, 3X, I4)
        !write(99,*)
        end do
        write(MeshInfo%unit,*) " "
    end do

    
!     do i=1,M-1
!       do j=1, N-1
!         dx(i,j)=x(i+1,j)-x(i,j)
!         dy(i,j)=y(i,j+1)-y(i,j)
!       end do
!     end do
    ! now we can calculate height of the cone...
    NeedleHeight=2d0*tan(.5d0*NeedleAngle) * sqrt(1d0 - NeedleRadius**2/((tan(NeedleAngle/2d0))**4)* &
                (x(Params%M,Params%N/2)+NeedleRadius**2/((tan(0.5d0*NeedleAngle))**2)**2)) &
              * (x(Params%M,Params%N/2)+NeedleRadius/(tan(NeedleAngle/2d0))**2)
  end if

  !rectangular mesh as Vessel domain, containing the cone
  if(MeshChoice.eq.2) then
    if(CathodeZone.eq.1) then
      xmax2=NeedleLength*2d0
      xmin2=-NeedleLength*1.0d0
    else
      xmax2=NeedleLength*1.0d0
      xmin2=-NeedleLength*1.0d0
    end if
  ymin2=-NeedleHeight
  ymax2=NeedleHeight
  ! adapted grid generation
  ! loop on iterations
  
  dx=(xmax2-xmin2)/(real(Mv-1))
  dy=(ymax2-ymin2)/(real(Nv-1))
  do k=1,MeshIterationsVessel
  ! define boundaries
  do i=1,Mv
    xV(i,1)=dx*real(i-1)+xmin2; yV(i,1)=ymin2
    xV(i,Nv)=xV(i,1); yV(i,Nv)=ymax2
  end do
  do j=1,Nv
    xV(Mv,j)=xmax2; yV(Mv,j)=dy*real(j-1)+ymin2
    xV(1,j)=xmin2; yV(1,j)=yV(Mv,j)
  end do

! ! sort points before solving
!   do j=1,Mv
!     do i=1,Mv
!       if(xV(i,1) > xV(j,1)) then 
!         call swap(xV(i,1), xV(j,1))
!         call swap(yV(i,1), yV(j,1))
!       end if
!       if(xV(i,Nv) > xV(i+1,Nv)) then 
!         call swap(xV(i,Nv), xV(j,Nv))
!         call swap(yV(i,Nv), yV(j,Nv))
!       end if
!     end do
!   end do
    
!   do k=1,MeshIterationsVessel
!     solve Laplace equations
    do i=2, Mv-1
      do j=2,Nv-1
!         if(xV(i,j) < NeedleLength .AND. xV(i,j)>0d0 .AND. abs(yV(i,j)) < NeedleHeight/2d0) then !if this point is included in the needle
!           MeshDensity(i,j)=MeshDensity0* &
!                         exp(-0.5d0*( (xV(Mv/2,Nv/2)/(NeedleRadius/(2d0*sqrt(2d0*log(2.)))))**2 &
!                         +(yV(Mv/2,Nv/2)/(NeedleRadius/(2d0*sqrt(2d0*log(2.)))))**2 ) )
                    !*( & !lets define a mesh density function near from boundaries
!                     exp(-xV(i,j)*MeshDamping0) * exp(-0.5d0*(yV(i,j))**2) + exp(-MeshDamping0*abs(yV(i,j)-NeedleHeight/2d0)))
!         else 
!           MeshDensity(i,j)=0d0
!         endif
        xV(i,j)=(xV(i+1,j)+xV(i-1,j)+xV(i,j+1)+xV(i,j-1))/4d0
        yV(i,j)=(yV(i+1,j)+yV(i-1,j)+yV(i,j+1)+yV(i,j-1))/4d0
        
        end do
    end do
     
  end do !end of mesh iterations
  
  !write the new mesh to its file
  
  
  ! regular grid generation
!     dx=(xmax-xmin)/real(Mv+1)
!     dy=(ymax-ymin)/real(Nv+1)
!     do i=1,Mv
!       do j=1,Nv
!         xV(i,j)=0d0
!         yV(i,j)=0d0
!       end do
!     end do

!     do i=1,Mv
!       do j=1,Nv
!         xV(i,j)=dx*real(i)+xmin
!         yV(i,j)=dy*real(j)+ymin
!       end do
!     end do
    ! use to select the mesh in which Poisson eq is solved
  end if
  
  !Attribute the node positions
  if(MeshChoice.eq.3) then
    do i=1, Params%M
      do j=1, Params%N
    x(i,j)=vertices(idvertices(MeshVertice(i,j)),1)
    y(i,j)=vertices(idvertices(MeshVertice(i,j)),2)
      end do
    end do

    x(:,:)=1d-6*x(:,:)
    y(:,:)=1d-6*y(:,:)
    
    !writing of the mesh
    do i=1,Params%M
      do j=1,Params%N
        write(MeshInfo%unit, 8852, advance='yes') x(i,j), y(i,j), i, j, MeshVertice(i,j)
        8852        FORMAT (1E15.8, 3x, 1E15.8, 3x, I4, 3X, I4, 3X, I8)
        !write(99,*)
        end do
        write(MeshInfo%unit,*) " "
    end do
    
  end if

  if(MeshChoice.eq.2) then
     Mp=Mv
     Np=Nv
     allocate(xP(1:Mp, 1:Np))
     allocate(yP(1:Mp, 1:Np))
     xP(:,:)=xV(:,:)
     yP(:,:)=yV(:,:)
  else 
     Mp=Params%M
     Np=Params%N
     allocate(xP(1:Mp, 1:Np))
     allocate(yP(1:Mp, 1:Np))
     xP(:,:)=x(:,:)
     yP(:,:)=y(:,:)
  end if
  
  do i=1,Mp
    do j=1,Np
      write(MeshVessel%unit,885) xP(i,j), yP(i,j), i, j
    end do
    write(MeshVessel%unit,*) " "
  end do
      
  
  write (*,*) "Vessel cells:", Mp, "*", Np,"=", Mp*Np
  write (*,*) "Matter cells:", Params%M, "*", Params%N, "=", Params%M*Params%N
  
!   allocate(Amatrix(1:Mp*Np,1:Mp*Np))
  allocate(Bvector(1:Mp*Np))
  allocate(Xvector(1:Mp*Np))
  allocate(XvectorPrev(1:Mp*Np))
  allocate(spectralNorm(1:Mp*Np))
  allocate(ExPoisson(1:Mp,1:Np))
  allocate(EyPoisson(1:Mp,1:Np))
  allocate(potential(1:Mp,1:Np))
  allocate(DielectricStatic(1:Mp,1:Np))
  allocate(DummyVessel(1:Mp, 1:Np))
  allocate(NeP(1:Mp, 1:Np))
  allocate(NhP(1:Mp, 1:Np))
  allocate(FixedPotentialIndex(1:Params%N, 1:2))
  allocate(CellVolume(1:Mp, 1:Np))
  allocate(NormalNxP(1:Mp,1:Np)); allocate(NormalNyP(1:Mp,1:Np))
  allocate(NormalSxP(1:Mp,1:Np)); allocate(NormalSyP(1:Mp,1:Np))
  allocate(NormalExP(1:Mp,1:Np)); allocate(NormalEyP(1:Mp,1:Np))
  allocate(NormalWxP(1:Mp,1:Np)); allocate(NormalWyP(1:Mp,1:Np))
  allocate(CellAreaNP(1:Mp,1:Np)); allocate(CellAreaSP(1:Mp,1:Np))
  allocate(CellAreaEP(1:Mp,1:Np)); allocate(CellAreaWP(1:Mp,1:Np))
  
  call flush(MeshInfo%unit); call flush(MeshVessel%unit)
  
  ! test of neighboors ! USEFUL FOR POISSON EQUATION ONLY
!  write(*,*) "[TEST] Testing neighborhood of a internal point..."
!  write(*,'(2I10.1)') transpose(FindNeighbours(x(2,2),y(2,2),xP,yP,numberOfNeighbours,Mp,Np))
!  write(*,*) "[TEST] Testing neighborhood of a external point..."
!  write(*,'(2I10.1)') transpose(FindNeighbours(xP(34,35),yP(34,35),x,y,numberOfNeighbours,M,N))

!*********** SOURCE IMPORT **************
  ! importing Fermi functions
   allocate(FermiTableE(1:9, 1:FermiMaxLines))
   allocate(FermiTableH(1:9, 1:FermiMaxLines))
   call TabCreateFL(FermiMaxLines, FermiTableE, FermiTableH)
   FermiTableE(:,:)=1d0; FermiTableH(:,:)=1d0; ! TODO: before publishing, this must work without inducing noise!
    !uncomment if you want to disable fermi-dirac. Dont forget to lock the FermiIndexes also.
!************ INITIALIZATION ************

  write(Parameters%unit,*) "========== CONE PARAMETERS ========="
  write(Parameters%unit,*) "Cone length=", NeedleLength*1d6, "um"
  write(Parameters%unit,*) "Cone height=", NeedleHeight*1d6, "um"
  write(Parameters%unit,*) "Cone angle=", NeedleAngleDeg, "deg"
  write(Parameters%unit,*) "Cone curvature radius=", NeedleRadius*1d9, "nm"
  !TODO: Move to source
  write(Parameters%unit,*)
  write(Parameters%unit,*) "=========== LASER PARAMETERS ========"
  write(Parameters%unit,*) "Laser fluence=", laser%fluence*1d-4, "J.cm-2"
  write(Parameters%unit,*) "Laser pulse duration=", laser%tau*1d15, "fs"
  write(Parameters%unit,*) "Laser wavelength=", laser%lambda*1d9, "nm"
  write(Parameters%unit,*) "Laser spot position: (X,Y)=", x0*1d6, y0*1d6, "um"
  write(Parameters%unit,*) "Laser spot size: (Sx, Sy)=", laser%spotX*1d6, laser%spotY*1d6, "um"
  write(Parameters%unit,*) "Mie scattering:", Params%UseMieScattering
  write(Parameters%unit,*) "Laser polarization", PolarizationSource
  write(Parameters%unit,*)
  write(Parameters%unit,*) "============ MESH PARAMETERS =========="
  write(Parameters%unit,*) "Mesh size", Params%M, "x", Params%N
  write(Parameters%unit,*) "Dilatation time ratio=", coeffDilaDt
  write(Parameters%unit,*) "Mesh shift=", MeshShift
  write(Parameters%unit,*)
  write(Parameters%unit,*) "============ TIME CONTROL =========="
  write(Parameters%unit,*) "Initial timestep=", dt0
  write(Parameters%unit,*) "Maximal timestep=", tmin
  write(Parameters%unit,*) "Maximum time t=", tmax
  write(Parameters%unit,*) "Enable adaptative timestep=", AdaptativeTimeStep
  write(Parameters%unit,*) "Time output each ", iterOut, "iterations."
  write(Parameters%unit,*) "Map output each", iterOutMaps*iterOut, "iterations."
  
  call flush(Parameters%unit)
  
  call TimerStart( )

  nmax=int((tmax-tmin)/dt, 8)
  
  Ex(:,:)=0d0 !-1d10
  Ey(:,:)=0d0 !-1d9 !0d0
  DummyVessel(:,:)=0d0
  
  Te0=Tout !1400d0
  Th0=Tout
  
  ! $ OMP DO
  do j=1,Params%N
    do i=1,Params%M
        
        VeX(i,j)=0d0
        VeY(i,j)=0d0
        VhX(i,j)=0d0
        VhY(i,j)=0d0
        newmesh%Te(i,j)=Te0
        newmesh%Th(i,j)=Th0
        newmesh%Ts(i,j)=Tout
        TsOld(i,j)=Tout
        TsPrev(i,j)=Tout
        
        if(BandBendingInFDTD.eq.1) then
           newmesh%Ne(i,j)=Ne0+Nborder*(exp(-0.5d0*(((x(i,j)-x(i,Params%N))**2+(y(i,j)-y(i,Params%N))**2) &
                    /((DefectThickness)/(2d0*sqrt2ln2))**2)) &
                  +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness) & 
                    /(2d0*sqrt2ln2))**2)) &
                  +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness) & 
                    /(2d0*sqrt2ln2))**2)) &
                  )
          newmesh%Nh(i,j)=Nh0+Nborder*(exp(-0.5d0*(((x(i,j)-x(i,Params%N))**2+(y(i,j)-y(i,Params%N))**2) &
                    /((DefectThickness)/(2d0*sqrt2ln2))**2)) &
                  +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness) & 
                    /(2d0*sqrt2ln2))**2)) &
                  +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness) & 
                    /(2d0*sqrt2ln2))**2)) &
                  )
        else
          newmesh%Ne(i,j)=Ne0
          newmesh%Nh(i,j)=Nh0
!           Ne(i,j)=Ne0 !NeNew(i,j)
!           Nh(i,j)=Nh0 !NhNew(i,j)
        end if
        
    end do
  end do
  ! $ OMP END DO
        
  call copy_mesh(mesh, newmesh)

  ! $ OMP DO
  do j=1,Params%N
    do i=1,Params%M
        
        ! initialise variables to calculate Ce, Ch
        DOSe(i,j)=DensityOfState(meDOS, mesh%Te(i,j))
        DOSh(i,j)=DensityOfState(mhDOS, mesh%Th(i,j))
        FermiRatioE(i,j)=mesh%Ne(i,j)/DOSe(i,j)
        FermiRatioH(i,j)=mesh%Nh(i,j)/DOSh(i,j)
        FermiIndexE(i,j)=1! FermiIndex(FermiRatioE(i,j), FermiMaxLines) !1
        FermiIndexH(i,j)=1! FermiIndex(FermiRatioH(i,j), FermiMaxLines) !1
!         write(*,*) "iter=", nbiter, "DOS=", DOSe(i,j), DOSh(i,j)

      etae=FermiTableE(ColFermiEta,FermiIndexE(i,j))
        etah=FermiTableH(ColFermiEta,FermiIndexH(i,j))

        ! calculate semi-classical heat capacity
        CeOld(i,j)=1.5d0*mesh%Ne(i,j)*kb*(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j)) &
                              /FermiTableE(ColFermiHalf,FermiIndexE(i,j)) &
                  -etae*(1d0-(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j)) &
                              /FermiTableE(ColFermiHalf,FermiIndexE(i,j)))* &
                  (FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)) &
                              /FermiTableE(ColFermiHalf,FermiIndexE(i,j)))))
        ChOld(i,j)=1.5d0*mesh%Nh(i,j)*kb*(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j)) &
                              /FermiTableH(ColFermiHalf,FermiIndexH(i,j)) &
                  -etah*(1d0-(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j)) &
                              /FermiTableH(ColFermiHalf,FermiIndexH(i,j)))* &
                  FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)) &
                              /FermiTableH(ColFermiHalf,FermiIndexH(i,j))))
         CsOld(i,j)=LatticeHeatCapacity(Tout)
         CsPrev(i,j)=LatticeHeatCapacity(Tout)
         CsPrev2(i,j)=LatticeHeatCapacity(Tout)


        Ce(i,j)=CeOld(i,j)
        Ch(i,j)=ChOld(i,j)
        Cs(i,j)=LatticeHeatCapacity(Tout)

        
        UeNew(i,j)=newmesh%Te(i,j)*CeOld(i,j)
        UhNew(i,j)=newmesh%Th(i,j)*ChOld(i,j)
        
        intensity(i,j)=0d0
        MaxHeating(i,j)=0d0
        MaxHeatingTime(i,j)=0d0
        epsilonNeedle(i,j)=epsilonStatic0-1d0
        
    end do
   end do
   ! $ OMP END DO
!!!! THIS PART IS ONLY USEFUL IF Poisson equation is solved   
!   if(CathodeZone.eq.1) then 
!    !one need the indexes of the needle points where potential is imposed: matrix Nx2
!    write(*,*) "Indexes of needle base on the vessel mesh"
!    do j=1, N
!      SomeNeighbours=FindNeighbours(x(M,j),y(M,j),xP,yP,numberOfNeighbours,Mp,Np) !OK
!      FixedPotentialIndex(j,1:2)=SomeNeighbours(2,1:2) !select the second point for each element corresponding to needle base
! !     write(*,*) FixedPotentialIndex(j,:) !OK even if points are not very regular!
!    end do
!   end if
   
   potential(:,:)=0d0! (0d0,0d0)
   
    ! $ O M P PARALLEL DEFAULT (SHARED) 
   
!   if(InterpolateOff.eq.0) then 
!    ! from needle mesh to vessel mesh
!      if(InterpolateMethod.eq.0) then
!        DielectricStatic=InterpolateSelner(epsilonNeedle,x,y,xP,yP,M,N,Mp,Np)
!!       else 
!!         DielectricStatic=InterpolateBiCubic(epsilonNeedle, x, y, xP, yP, M, N, Mp, Np)
!!         call InterpolateBiCubic(epsilonNeedle, x, y, xP, yP, M, N, Mp, Np, DielectricStatic, DummyVessel, DummyVessel)
!      end if
!   endif
   ! $ O M P END PARALLEL
   DielectricStatic(:,:)=DielectricStatic(:,:)+1d0 !so that equals 1 outside needle, and equals epsilonStatic0 inside
   write(*,*) 'Done.'
   
   ! potential on vessel boundaries
   potential(Mp,:)=potentialNull
   potential(:,1)=potentialNull
   potential(:,Np)=potentialNull
   potential(1,:)=potentialNull
   ! potential on right boundary on the needle
   potentialNeedle(Params%M,:)=potential0
   
   ! areas of the cell boudaries e.g. distances in this 2D code 
   ! and normal vectors (Nx, Ny) the four poles of quadrangle elements
   write(*,*) "[Mesh] Calculation of normals and distances."

   allocate(NormalN%x(1:Params%M,1:Params%N));
   allocate(NormalN%y(1:Params%M,1:Params%N));
   allocate(NormalN%N(1:Params%M,1:Params%N))
   allocate(NormalS%x(1:Params%M,1:Params%N));
   allocate(NormalS%y(1:Params%M,1:Params%N));
   allocate(NormalS%N(1:Params%M,1:Params%N))
   allocate(NormalE%x(1:Params%M,1:Params%N));
   allocate(NormalE%y(1:Params%M,1:Params%N));
   allocate(NormalE%N(1:Params%M,1:Params%N))
   allocate(NormalW%x(1:Params%M,1:Params%N));
   allocate(NormalW%y(1:Params%M,1:Params%N));
   allocate(NormalW%N(1:Params%M,1:Params%N))

    allocate(CurviWx(1:Params%M,1:Params%N));
    allocate(CurviWy(1:Params%M,1:Params%N))
    allocate(CurviEx(1:Params%M,1:Params%N));
    allocate(CurviEy(1:Params%M,1:Params%N))
    allocate(CurviNx(1:Params%M,1:Params%N));
    allocate(CurviNy(1:Params%M,1:Params%N))
    allocate(CurviSx(1:Params%M,1:Params%N));
    allocate(CurviSy(1:Params%M,1:Params%N))
    allocate(TangentWx(1:Params%M,1:Params%N));
    allocate(TangentWy(1:Params%M,1:Params%N))
    allocate(TangentEx(1:Params%M,1:Params%N));
    allocate(TangentEy(1:Params%M,1:Params%N))
    allocate(TangentNx(1:Params%M,1:Params%N));
    allocate(TangentNy(1:Params%M,1:Params%N))
    allocate(TangentSx(1:Params%M,1:Params%N));
    allocate(TangentSy(1:Params%M,1:Params%N))

   call compute_distances(Params%M, Params%N, x, y, DistN, DistS, DistE, DistW, DistDualN, &
                          DistDualS, DistDualE, DistDualW, CellAreaN, CellAreaS, CellAreaE, CEllAreaW )
   !
   call compute_norm_tan_curv(Params%M, Params%N, x, y, NormalN%x, NormalN%y, NormalS%x, NormalS%y, &
                              NormalE%x, NormalE%y, NormalW%x, NormalW%y, TangentNx, TangentNy, TangentSx, TangentSy, &
                              TangentEx, TangentEy, TangentWx, TangentWy, &
                              CurviNx, CurviNy, CurviSx, CurviSy, CurviEx, CurviEy, CurviWx, CurviWy )
   !
   call compute_cellvol(Params%M, Params%N, x, y, CellVol, InvCellVol )
   

!     if(DisableCrossDiffusion.eq.1) then
!       ! disable flux parallel to element boundaries
!       TangentNx(:,:)=0d0; TangentNy(:,:)=0d0; 
!       TangentSx(:,:)=0d0; TangentSy(:,:)=0d0
!       TangentEx(:,:)=0d0; TangentEy(:,:)=0d0
!       TangentWx(:,:)=0d0; TangentWy(:,:)=0d0
!       ! unnormalize flux normal to element boundaries
!       CurviNx(:,:)=NormalN%x(:,:); CurviNy(:,:)=NormalN%y(:,:)
!       CurviSx(:,:)=NormalS%x(:,:); CurviSy(:,:)=NormalS%y(:,:)
!       CurviEx(:,:)=NormalE%x(:,:); CurviEy(:,:)=NormalE%y(:,:)
!       CurviWx(:,:)=NormalW%x(:,:); CurviWy(:,:)=NormalW%y(:,:)
!     end if
   
   write(*,*) "NEEDLE CHECK"
   write(*,*) "North", NormalN%x(Params%M/2,Params%N-1), NormalN%y(Params%M/2,Params%N-1)
   write(*,*) "South", NormalS%x(Params%M/2,2), NormalS%y(Params%M/2,2)
   write(*,*) "East", NormalE%x(Params%M-1,Params%N-1), NormalE%y(Params%M-1,Params%N-1)
   write(*,*) "West", NormalW%x(Params%M-1,Params%N-1), NormalW%y(Params%M-1,Params%N-1)
   write(*,*) CellAreaN(Params%M/2,Params%N/2), CellAreaS(Params%M/2,Params%N/2), &
              CellAreaE(Params%M/2,Params%N/2), CellAreaW(Params%M/2,Params%N/2)
   
   
   call poisson_init_normal_cellarea( Mp, Np, xP, yP, NormalNxP, NormalNyP, NormalSxP, NormalSyP, NormalExP, NormalEyP, &
                                            NormalWxP, NormalWyP, CellAreaNP, CellAreaSP, CellAreaEP, CellAreaWP )


  call poisson_init_dual( Params%M,Params%N, x, y, xDualSW, yDualSW, xDualSE, yDualSE, &
                           xDualNE, yDualNE, xDualNW, yDualNW, xDual, yDual)

   ! Drift initialization
   JeX(:,:)=0d0
   JeY(:,:)=0d0
   JhX(:,:)=0d0
   JhY(:,:)=0d0
  
   
   !! defining material index and ionization constants
   epsilonInf=DielectricConstant(laser%lambda)
   OnePhotonIonizationRate0=OnePhotonIonizationRate()
   TwoPhotonIonizationRate0=TwoPhotonIonizationRate(laser%lambda)
   
  write(*,*) 'epsilon(', 1d9*laser%lambda, 'nm)=', epsilonInf
  write(*,*) 'Re(sqrt(epsilon))=', real(sqrt(epsilonInf))
if(Params%UseMieScattering.eq.1) then
  write(*,*) 'Computing the Mie scattering field distribution...'
  write(*,*) 'Angle Mie =', phiMie0
  write(*,*) 'Polarization TM ? ', PolarizationSource

  !$OMP PARALLEL DEFAULT(NONE) SHARED(x, y, Params, &
  !$OMP phiMie, Radius, phiMie0)
  !$OMP DO COLLAPSE(2)
  do j=1,Params%N
    do i=1,Params%M
        if(y(i,j)<0d0) then
          phiMie(i,j)=phiMie0+pi
        else
          phiMie(i,j)=phiMie0
        end if
!           find the radius for the cylindrical Mie scattering model
!           Radius(i,j)=y(i,j)
        if(ExpNeedleType.eq.0) then
          Radius(i,j)=ConeExp1Radius(1d6*y(i,j), 1d6*x(i,j), 0d0)
        else
          Radius(i,j)=ConeExp2Radius(1d6*y(i,j), 1d6*x(i,j), 0d0)
        end if
      end do
    end do
  !$OMP END DO
  !$OMP END PARALLEL
    

  !$OMP PARALLEL DEFAULT(NONE) SHARED(x, y, laser, Params, epsilonInf, &
  !$OMP phiMie, Radius, EintField,EintField2, PolarizationSource)
  !$OMP DO COLLAPSE(2)
   do j=1,Params%N
    do i=1,Params%M
          if(PolarizationSource.eq.1) then !TM polarization, Bassel et al scattering on a cylinder
          ! formula for an experimental needle with interpolated radius
!             write(*,*) "TM polarization selected."
            EintField(i,j)= M_ONE * MieScattering(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf, laser%k) ! * sqrt(2d0*laser%fluence/(c*epsilon0*laser%tau))
            EintField2(i,j)=M_ZERO
          ! formula with a super mistake on radius
!           EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), 0.5d0*(y(i,N)-y(i,1)), epsilonInf) ! * sqrt(2d0*laser%fluence/(c*epsilon0*laser%tau))

          ! formulas for an hyperbolic needle
!           EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), abs(ContourYofX(x(i,j), NeedleRadius, NeedleAngle)), epsilonInf) ! * sqrt(2d0*laser%fluence/(c*epsilon0*laser%tau))

          ! formula for debug, using a constant radius
!         EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), 100d-9, epsilonInf) ! * sqrt(2d0*laser%fluence/(c*epsilon0*laser%tau)) !with a constant radius

          else !TE polarization
!             write(*,*) "TE polarization selected."
            EintField2(i,j)=M_ONE * MieScatteringTE2(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf, laser%k)
            EintField(i,j) =M_ONE * MieScatteringTE1(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf, laser%k)
          end if
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

!     EintFieldR=sqrt(EintField * conjg(EintField))
    EintFieldR=real(sqrt( EintField * conjg(EintField) + EintField2 * conjg(EintField2) ))
    

    
    do j=1,Params%N
      do i=1,Params%M
        write(Field%unit, 891, advance='yes') x(i,j), y(i,j), (EintFieldR(i,j)**2d0)**0.5d0, Radius(i,j)
891        FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
      end do
    end do


  write(*,*) 'Done.'
 end if
 
  t=tmin
  
  write(*,*) 'Elapsed time : ', ElapsedTime ( )
  write(*,*) "[Poisson eq.] Filling matrix"
!   if(PoissonOn.eq.1 .AND. PoissonSolver.eq.0) then
!     !******************* Let's make a Gauss inversion of Amatrix here !
!     ! matrix filling first
!     Amatrix(:,:)=0d0;
!     do i=1,Mp
!       do j=1, Np
!           Amatrix(OneDindex(i,j),OneDindex(i,j))=1d0 ! Diag(A)=1d0 for Dirichlet conditions !OK
!       end do
!     end do
!     
!     do i=2, Mp-1 !pour chaque point ou l on va calculer le potentiel
!         do j=2, Np-1 
!           Amatrix(OneDindex(i,j),OneDindex(i,j))=(1d0/8d0)*(-DielectricStatic(i+1,j)-DielectricStatic(i,j))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i,j)+xP(i,j)**2+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)-(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i-1,j))*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2 &
 
!           -2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i-1,j)**2-2d0*xP(i-1,j)*xP(i,j)+xP(i,j)**2+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)+(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i,j+1))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)/(xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i,j)+xP(i,j)**2+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i,j)+yP(i,j)**2)**(0.5d0)-(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i,j-1))*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)- &
!           2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5d0)/(xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i,j)+xP(i,j)**2+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
!          
!           Amatrix(OneDindex(i,j),OneDindex(i,j-1))=-(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i,j-1))*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)-2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5d0)/(xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i,j)+xP(i,j)**2+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
!           ! 0.25d0*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)-2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5) !a(i,j-1)
!           Amatrix(OneDindex(i,j),OneDindex(i,j+1))=(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i,j+1))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)/(xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i,j)+xP(i,j)**2+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
!           ! 0.25d0*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5) !a(i,j+1)
!           Amatrix(OneDindex(i,j),OneDindex(i-1,j))=-(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i-1,j))*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i-1,j)**2-2d0*xP(i-1,j)*xP(i,j)+xP(i,j)**2+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)
!           !0.25d0*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5) !a(i-1,j)
!           Amatrix(OneDindex(i,j),OneDindex(i+1,j))=(1d0/8d0)*(DielectricStatic(i+1,j)+DielectricStatic(i,j))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i,j)+xP(i,j)**2+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! !           0.25d0*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5) !a(i+1,j)
! 
!       end do
!     end do !well filled
!     
!     if(CathodeZone.eq.1) then
!       !diagonal equals to 1
!       do k=FixedPotentialIndex(1,2)-ShiftFixedPotential, FixedPotentialIndex(N,2)+1+ShiftFixedPotential ! filling the indexes of fixed potential with 1d0, not localized on boundaries
!               do j=FixedPotentialIndex(1,2)-ShiftFixedPotential, FixedPotentialIndex(N,2)+1+ShiftFixedPotential
!           !other points have to be zero
!           Amatrix(OneDindex(FixedPotentialIndex(1,1), k),&
!           OneDindex(FixedPotentialIndex(1,1), j))=0d0
!         end do
!         Amatrix(OneDindex(FixedPotentialIndex(1,1), k),&
!                 OneDindex(FixedPotentialIndex(1,1), k))=1d0
! 
!       end do
!     end if
!     
!     
!     
!     
! !     then matrix inversion
!     write(*,*) "[Poisson eq.] Preparing Matrix inversion..."
!     Amatrix=LaInv(Amatrix)
!     write(*,*) "[Poisson eq.] Matrix inverted."    
! !     do k=1,Mp*Np
! !       imax=maxval(abs(Amatrix(:,k)))
! !       if (Amatrix(imax, k).eq.0d0 ) then
! !         write(*,*) '[Gauss inversion] Singular matrix!'
! !       end if
! !       do i=1,Mp*Np
! !         call swap(Amatrix(i,k), Amatrix(i,imax))
! !       end do
! !       do i=k+1, Mp*Np
! !         do j=k+1,Mp*Np
! !           Amatrix(i,j)=Amatrix(i,j)-Amatrix(k,j)*(Amatrix(i,k)/Amatrix(k,k))
! !         end do
! !         Amatrix(i,k)=0d0
! !       end do
! !     end do
!   
!   ! now use this matrix at each time step (carefull: it increases round off errors)
!   end if
 
  ! initialisation of the displayed values
  maxIntensity=0d0; maxTe=0d0; minTe=1d10; maxTh=0d0; minTh=1d10; maxTs=0d0; minTs=1d10; 
  maxNe=0d0; maxNh=0d0; minNe=1d50; minNh=1d50
  maxCFLxT=0d0; maxCFLyT=0d0; maxCFLxN=0d0; maxCFLyN=0d0; maxCFLxTs=0d0; maxCFLyTs=0d0; 
  maxFermiIndexE=0; maxFermiIndexH=0; 
 
  !TODO: NTD: Why DistX are recomputed here? Same for CellAreaX
  !TODO: TJYD: To treat boundary conditions and treat everything with a loop on the complete mesh. This should be kept. 
  do i=1,Params%M
    CellAreaN(i,Params%N)=0d0
!     CellAreaN(i,N-1)=0d0
!     CellAreaS(i,2)=0d0
    CellAreaS(i,1)=0d0
    
    DistN(i,Params%N-1)=sqrt((0.5d0*(x(i,Params%N-1)+x(i,Params%N))-x(i,Params%N-1))**2 &
                            +(0.5d0*(y(i,Params%N-1)+y(i,Params%N))-y(i,Params%N-1))**2)
    DistS(i,2)=sqrt((x(i,2)-0.5d0*(x(i,2)+x(i,1)))**2+(y(i,2)-0.5d0*(y(i,2)+y(i,1)))**2) 
  end do
   
  do j=1,Params%N
    CellAreaW(1,j)=0d0
!     CellAreaW(2,j)=0d0
!     CellAreaE(M-1,j)=0d0
    CellAreaE(Params%M,j)=0d0
    
    DistW(2,j)=sqrt((x(2,j)-0.5d0*(x(2,j)+x(1,j)))**2+((y(2,j)-0.5d0*(y(2,j)+y(1,j))))**2)
    DistE(Params%M-1,j)=sqrt((0.5d0*(x(Params%M-1,j)+x(Params%M,j))-x(Params%M-1,j))**2&
                            +(0.5d0*(y(Params%M-1,j)+y(Params%M,j))-y(Params%M-1,j))**2)
  end do
  
     
    do j=2, Params%N-1
      do i=2, Params%M-1
        ! Calcul de Grad(Ne) sur le maillage direct
        ! Première estimation peu stable
        GradNeX(i,j) = 0.5d0 * InvCellVol(i,j) * &
                      ( (mesh%Ne(i,j) + mesh%Ne(i,j+1)) * CellAreaN(i,j) * NormalN%x(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i,j-1)) * CellAreaS(i,j) * NormalS%x(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i-1,j)) * CellAreaW(i,j) * NormalW%x(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i+1,j)) * CellAreaE(i,j) * NormalE%x(i,j) )
        GradNeY(i,j) = 0.5d0 * InvCellVol(i,j) * &
                      ( (mesh%Ne(i,j) + mesh%Ne(i,j+1)) * CellAreaN(i,j) * NormalN%y(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i,j-1)) * CellAreaS(i,j) * NormalS%y(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i-1,j)) * CellAreaW(i,j) * NormalW%y(i,j) &
                      * (mesh%Ne(i,j) + mesh%Ne(i+1,j)) * CellAreaE(i,j) * NormalE%y(i,j) )
                      
        ! interpolation lineaire des valeurs de phi sur les bords de cellules
!         phiN=Ne(i,j)+GradNeX(i,j)*(0.5d0*(xDual(i,j)+xDual(i-1,j)))+GradNeY(i,j)*(0.5d0*(yDual(i,j)+yDual(i-1,j)))
!         phiS=Ne(i,j)+GradNeX(i,j)*(0.5d0*(xDual(i,j-1)+xDual(i-1,j-1)))+GradNeY(i,j)*(0.5d0*(yDual(i,j-1)+yDual(i-1,j-1)))
        ! recalcul du gradient avec phi_bords
        
        end do
      end do

  
  do j=1,Params%N
    do i=1,Params%M
       NormalN%N(i,j) = (NormalN%x(i,j)**2+NormalN%y(i,j)**2)
       NormalE%N(i,j) = (NormalE%x(i,j)**2+NormalE%y(i,j)**2)
       NormalW%N(i,j) = (NormalW%x(i,j)**2+NormalW%y(i,j)**2)
       NormalS%N(i,j) = (NormalS%x(i,j)**2+NormalS%y(i,j)**2)
    end do
  end do
  
  do j=2, Params%N-1
    do i=2, Params%M-1
      NormalN%N(i,j) = NormalN%N(i,j)/DistN(i,j)
      NormalW%N(i,j) = NormalW%N(i,j)/DistW(i,j)
      NormalE%N(i,j) = NormalE%N(i,j)/DistE(i,j)
      NormalS%N(i,j) = NormalS%N(i,j)/DistS(i,j)
      ShapeFactorNormalE(i,j)  = CellAreaE(i,j) / (CurviEx(i,j)*NormalE%x(i,j)+CurviEy(i,j)*NormalE%y(i,j))
      ShapeFactorNormalW(i,j)  = CellAreaW(i,j) / (CurviWx(i,j)*NormalW%x(i,j)+CurviWy(i,j)*NormalW%y(i,j))
      ShapeFactorNormalN(i,j)  = CellAreaN(i,j) / (CurviNx(i,j)*NormalN%x(i,j)+CurviNy(i,j)*NormalN%y(i,j))
      ShapeFactorNormalS(i,j)  = CellAreaS(i,j) / (CurviSx(i,j)*NormalS%x(i,j)+CurviSy(i,j)*NormalS%y(i,j))
      ShapeFactorTangentE(i,j) = CrossCoeff * ( CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j) ) / DistDualE(i,j)
      ShapeFactorTangentW(i,j) = CrossCoeff * ( CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j) ) / DistDualW(i,j)
      ShapeFactorTangentN(i,j) = CrossCoeff * ( CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j) ) / DistDualN(i,j)
      ShapeFactorTangentS(i,j) = CrossCoeff * ( CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j) ) / DistDualS(i,j)
    end do
  end do
  
  !TODO: Use an OutputData here
  open(102,FILE='meshElements.dat', access='sequential',status='unknown') ! format 881
  do i=1,Params%M
        do j=1,Params%N
          write(102, 881, advance="YES") i, j, x(i,j), y(i,j), NormalN%x(i,j), & !5
                                  NormalN%y(i,j), NormalS%x(i,j), NormalS%y(i,j), NormalE%x(i,j), NormalE%y(i,j), & !10
                                  NormalW%x(i,j), NormalW%y(i,j), CellAreaN(i,j), CellAreaS(i,j), CellAreaE(i,j), & !15
                                  CellAreaW(i,j), CellVol(i,j), TangentNx(i,j), TangentNy(i,j), TangentSx(i,j), & !20
                                  TangentSy(i,j), TangentEx(i,j), TangentEy(i,j), TangentWx(i,j), TangentWy(i,j), & !25
                                  CurviNx(i,j), CurviNy(i,j), CurviSx(i,j), CurviSy(i,j), CurviEx(i,j), & !30
                                  CurviEy(i,j), CurviWx(i,j), CurviWy(i,j), DistN(i,j), DistS(i,j), & !35
                                  DistE(i,j), DistW(i,j), DistDualN(i,j), DistDualS(i,j), DistDualE(i,j),  & !40
                                  DistDualW(i,j), xDualSW(i,j), yDualSW(i,j), xDualSE(i,j), yDualSE(i,j), & !45
                                  xDualNE(i,j), yDualNE(i,j), xDualNW(i,j), yDualNW(i,j) !, xDual(i,j), & !50
!                                  yDual(i,j)

881 FORMAT (2(I5, 3x), 2(1E16.8, 3x), 44(1E12.5,3x), 1E12.5)
        end do
      end do
   close(102)

   deallocate(CurviWx,CurviWy,CurviEx,CurviEy,CurviNx,CurviNy,CurviSx,CurviSy)

   deallocate(TangentWx,TangentWy,TangentEx,TangentEy,TangentNx,TangentNy,TangentSx,TangentSy)

  ! Steps for 3rd order integration
  
  h1=dt
  h2=dt+dt2
  h3=dt+dt2+dt3
  
  write(*,*) "Starting time loop."
  call TimerStart( )

  !***************************************************************
  !***************** temporal loop *******************************
  !***************************************************************
  do nbiter=1, nmax
    
    t=t+dt; 
    if(mod(nbiter,iterOut).eq.0) then 
      write(*,*) "Time=", t, "Intensity=", maxIntensity
      write(*,*) minTe, "< Te < ", maxTe
      write(*,*) minTh, "< Th < ", maxTh
      write(*,*) minTs, "< Ts <", maxTs
      write(*,*) minNe, "< Ne <", maxNe
      write(*,*) minNh, "< Nh <", maxNh
      write(*,*) "CFL_Te=", maxCFLxT+maxCFLyT, "maxCFL_Ne=", maxCFLxN+maxCFLyN, &
                   "maxCFL_Ts=",maxCFLxTs+maxCFLyTs, "CPU=", cpuefficiency, &
                   "NumThreads=", nthreads, "Elapsed time=", cpu_timestep_duration
      write(*,*) "CFL_Limit=", maxCFL
      write(*,*) "dt_init=", dt0, "dt=", dt
    end if
    
    maxIntensity=0d0; maxTe=0d0; minTe=1d10; maxTh=0d0; minTh=1d10; maxTs=0d0; minTs=1d10; 
    maxNe=0d0; minNe=1d50; maxNh=0d0; minNh=1d50; maxCFLxT=0d0; maxCFLyT=0d0; maxCFLxN=0d0; maxCFLyN=0d0; maxCFLxTs=0d0; 
    maxCFLyTs=0d0; maxSourceE=0d0; maxSourceH=0d0; maxGainsE=0d0; maxGainsH=0d0
    
   !$OMP PARALLEL DEFAULT(NONE) SHARED (dt, dt2, dt3, dt4, UeNew, UhNew, TsOld, TsPrev, &
   !$OMP& mesh, newmesh, dual, intensityDual, laser, Params, I0, &
   !$OMP& Ue, Uh, GradNeX, GradNeY, intensity, reflectivity, FermiTableE, FermiTableH, &
   !$OMP& Dielectric, DielectricDrudeE, DielectricDrudeH, absorptionDrudeE, absorptionDrudeH, &
   !$OMP& x, y, diffusionE, diffusionH, GainsE, GainsH, LossesE, LossesH, &
   !$OMP& kappae, kappah, kappas, Ce, CeOld, Ch, ChOld, Cs, CsOld, CsPrev, CsPrev2, CouplingE, CouplingH, &
   !$OMP& mobilityE, mobilityH, Egap, me, mh, meDOS, mhDOS, DOSe, DOSh, &
   !$OMP& SourceE, SourceH, SourceUe, SourceUh, diffNe, diffNh, CFLxT, CFLyT, CFLxN, CFLyN, CFLxTs, CFLyTs, &
   !$OMP& ThermalEnergy, LaserEnergy, epsilonInf, FermiIndexE, FermiIndexH, FermiRatioE, FermiRatioH, &
   !$OMP& JeX, JeY, JhX, JhY, VeX, VeY, VhX, VhY, DielectricStatic, Xvector, XvectorPrev, Bvector, xV, yV, xP, yP, & !Amatrix
   !$OMP& spectralNorm, Ex, Ey, ExPoisson, EyPoisson, potential, potentialNeedle, NeP, NhP, FixedPotentialIndex, &
   !$OMP& NormalN, NormalS, NormalE, NormalW, &
   !$OMP& CellVolume, CellAreaN, CellAreaS, CellAreaE, CellAreaW, CellVol, InvCellVol, &
   !$OMP& DistN, DistS, DistE, DistW, DistDualN, DistDualS, DistDualE, DistDualW, &
   !$OMP& EintField, EintFieldI, EintFieldR, NeTotal, NhTotal, &
   !$OMP& ShapeFactorNormalE, ShapeFactorNormalN, ShapeFactorNormalS, ShapeFactorNormalW, &
   !$OMP& ShapeFactorTangentE, ShapeFactorTangentN, ShapeFactorTangentS, ShapeFactorTangentW, &
   !$OMP& ColFermi0, ColFermi1, ColFermi2, ColFermiEta, ColFermiHalf, &
   !$OMP& x1, x2, x3, x4, x5, x6, x7, x8, x9, &
   !$OMP& y1, y2, y3, y4, y5, y6, y7, y8, y9, &
   !$OMP& sigmaX1, sigmaX2, sigmaX3, sigmaX4, sigmaX5, sigmaX6, sigmaX7, &
   !$OMP& sigmaX8, sigmaX9, &
   !$OMP& sigmaY1, sigmaY2, sigmaY3, sigmaY4, sigmaY5, sigmaY6, sigmaY7, &
   !$OMP& sigmaY8, sigmaY9, &
   !$OMP& ColFermiMenusHalf, ColFermiNeNc, ColFermiThreeHalf, &
   !$OMP& t, t0, x0, y0, I1, I2, I3, I4, I5, I6, I7, &
   !$OMP& I8, I9, OpticalIndex, OpticalDamping, &
   !$OMP& h1, h2, h3, OnePhotonIonizationRate0, TwoPhotonIonizationRate0, &
   !$OMP& AugerRateE, AugerRateH, sigmaTau, sigmaX, sigmaY, dx, dy, &
   !$OMP& Mp, Np, invCe, invCh, invCs, nuColl, TotalElectrons, TotalHoles) &
   !$OMP& PRIVATE(Int2, work)


   NeTotal=0d0; NhTotal=0d0
   
   ! replacing old datas
   Ue(:,:)=UeNew(:,:)
   Uh(:,:)=UhNew(:,:)
   TsPrev(:,:)=TsOld(:,:)
   TsOld(:,:)=mesh%Ts(:,:)
   CeOld(:,:)=Ce(:,:)
   ChOld(:,:)=Ch(:,:)
   
   CsPrev2(:,:)=CsPrev(:,:)
   CsPrev(:,:)=CsOld(:,:)
   CsOld(:,:)=Cs(:,:)
   
   call copy_mesh(mesh, newmesh)
   !
   !TODO: Does this depends on the position? If yes, this has to be changed bak to an array
   nuColl=CollisionFrequency()
   !
   call DielectricFunction_batch(mesh, Dielectric, OpticalIndex, OpticalDamping, Reflectivity, &
                                 epsilonInf, nuColl, me, laser)
   !
   call DielectricFunctionDrude_batch(mesh, mesh%Ne, DielectricDrudeE, absorptionDrudeE, DrudeHeating, nuColl, me, laser)
   !
   call DielectricFunctionDrude_batch(mesh, mesh%Nh, DielectricDrudeH, absorptionDrudeH, DrudeHeating, nuColl, mh, laser)
   !
   call DensitiesOfState_batch(mesh, DOSe, DOSh, meDOS, mhDOS)
   !
   !This routine computes the intensity for the entire grid with one call
   call ComputeIntensity_batch(Params, mesh, laser, intensity, OpticalIndex, Reflectivity, &
                               absorptionDrudeE, absorptionDrudeH, OnePhotonIonizationRate0, TwoPhotonIonizationRate0, &
                               t, t0, sigmaTau, I0, sigmaX, sigmaY, x, y, x0, y0, DefectThickness, BandBendingInFDTD,  &
                               sigmaX1, sigmaY1, sigmaX2, sigmaY2, sigmaX3, sigmaY3, sigmaX4, sigmaY4, sigmaX5, sigmaY5, &
                               sigmaX6, sigmaY6, sigmaX7, sigmaY7, sigmaX8, sigmaY8, sigmaX9, sigmaY9, x1, y1, x2, y2, &
                               x3, y3, x4, y4, x5, EintFieldR,  &
                               y5, x6, y6, x7, y7, x8, y8, x9, y9, I1, I2, I3, I4, I5, I6, I7, I8, I9 )
   !
   !
   !
!!!! thermal calculations in the main domain
! calculation of sources
    !$OMP DO  COLLAPSE(2) 
    do j=1,Params%N
        do i=1,Params%M

        !         write(*,*) "Esprit es-tu la ?"
        FermiRatioE(i,j)=mesh%Ne(i,j)/DOSe(i,j)
        FermiRatioH(i,j)=mesh%Nh(i,j)/DOSh(i,j)
        FermiIndexE(i,j)=1! FermiIndex(FermiRatioE(i,j), FermiMaxLines) !1
        FermiIndexH(i,j)=1! FermiIndex(FermiRatioH(i,j), FermiMaxLines) !1
!         write(*,*) "iter=", nbiter, "DOS=", DOSe(i,j), DOSh(i,j)
!         write(*,*) "iter=", nbiter, "NeNc=", Ne(i,j)/DOSe(i,j), Nh(i,j)/DOSh(i,j)
!         write(*,*) "iter=", nbiter, "FermiIndex=", FermiIndex(Ne(i,j)/DOSe(i,j), FermiMaxLines), FermiIndex(Nh(i,j)/DOSh(i,j), FermiMaxLines)
!         write(*,*) "iter=", nbiter, "eta=", etae, etah
!         write(*,*) "FermiTables: etaE,etaH=", FermiTableE(3,463), FermiTableH(3,450)

      end do
   end do
   !$OMP END DO
   !
   !
   !Computes the electron and mobilities for the entire mesh
   call ComputeMobilities_batch(mesh, mobilityE, mobilityH, &
                                FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                ColFermi0, ColFermiHalf, nuColl, me)
   !
   !
   !Computes the diffusion terms for the entire mesh
   call ComputeDiffusions_batch(mesh, diffusionE, diffusionH, mobilityE, mobilityH, &
                                FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                ColFermiHalf, ColFermiMenusHalf, ConductivityFix)
   !
   !
   !Computes the drif vectors for the entire mesh
   call ComputeDriftVectors_batch(mesh, JeX, JeY, JhX, JhY, mobilityE, mobilityH, Ex, Ey, DriftOn)
   !
   !
   ! Computes the electron, hole and lattice heat capacities for the entire mesh
   call ComputeHeatCapacities_batch(mesh, Ce, Ch, Cs, invCe, invCh, invCs, &
                                    FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                    ColFermiThreeHalf, ColFermiHalf, ColFermiMenusHalf, ColFermiEta)
   !
   !Computes the couplings for the entire mesh
   call ComputeCouplings_batch(mesh, CouplingE, CouplingH, Ce, Ch, CouplingDebug, Params%HolesOff)
   !
   !
   ! calculation of sources
   !$OMP DO  COLLAPSE(2)
   do j=1,Params%N
     do i=1,Params%M

        ! free-carrier balance sources
        Egap(i,j)=EgapValue(mesh%Ne(i,j),mesh%Ts(i,j))

        Int2 = intensity(i,j)**2
        work = ImpactIonizationRate(mesh%Te(i,j),Egap(i,j), ImpactOff)

        GainsE(i,j)=(OnePhotonIonizationRate0*intensity(i,j)*laser%inv_E &
                    +0.5d0*TwoPhotonIonizationRate0*Int2*laser%inv_E &
                    +work*mesh%Ne(i,j))! *(4d0*SiDensity-Ne(i,j))/(4d0*SiDensity) !use Old Ne here!

        SourceUe(i,j)= ((laser%E-Egap(i,j))*OnePhotonIonizationRate0*intensity(i,j) &
                     + 0.5d0*(2d0*laser%E - Egap(i,j))*TwoPhotonIonizationRate0*Int2 )*laser%inv_E*((me)/(me+mh))&
                     - Egap(i,j)*work*mesh%Ne(i,j) &
                     + absorptionDrudeE(i,j)*intensity(i,j) &
                     + Egap(i,j)*(AugerRateE*mesh%Nh(i,j) * mesh%Ne(i,j)**2d0)

        !SourceE(i,j) = SourceE(i,j) - diffNe(i,j)*(1.5d0*kb*Te(i,j))*(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))
        SourceE(i,j) = SourceUe(i,j) - mesh%Te(i,j) * (Ce(i,j)-CeOld(i,j))/dt

        !LossesE(i,j)=AugerRateE * (mesh%Ne(i,j))**2d0 * mesh%Nh(i,j) + AugerRateH * (mesh%Nh(i,j))**2d0 * mesh%Ne(i,j) !use Old Ne, Nh here!
        LossesE(i,j)=mesh%Ne(i,j) * mesh%Nh(i,j) * ( AugerRateE * mesh%Ne(i,j) + AugerRateH * mesh%Nh(i,j) ) !This is more perfomant like that

        work = ImpactIonizationRate(mesh%Th(i,j),EgapValue(mesh%Nh(i,j),mesh%Ts(i,j)), ImpactOff)

        GainsH(i,j)=(OnePhotonIonizationRate0*intensity(i,j)*laser%inv_E &
                    +0.5d0*TwoPhotonIonizationRate0*Int2*laser%inv_E &
                    +work*mesh%Nh(i,j)) !*(4d0*SiDensity-Ne(i,j))/(4d0*SiDensity) !use Old Nh here
                    
        SourceUh(i,j)=((laser%E-Egap(i,j))* OnePhotonIonizationRate0*intensity(i,j) &
                     + 0.5d0*(2d0*laser%E - Egap(i,j))*TwoPhotonIonizationRate0*Int2)*laser%inv_E * ((me)/(me+mh))  &
                     - Egap(i,j)*work*mesh%Nh(i,j) &
                     + absorptionDrudeH(i,j)*intensity(i,j) &
                     + Egap(i,j)*(AugerRateH*mesh%Ne(i,j) * mesh%Nh(i,j)**2d0)
                     
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

    !Compute the new conductivites, based on the knowledge of densities and mobilities
    call ComputeConductivities_batch(mesh, kappae, kappah, kappas, mobilityE, mobilityH, &
                                     FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                     ColFermi0, ColFermi1, ColFermi2, ConductivityFix)
    !
    ! interpolation bilineaire ponderee par les aires
    call bilinear_interpol_dual(mesh, dual, InvCellVol)
    !
    !
    ! solving the 2D problem
    !
    !
    !
    if(Params%NeOff.eq.0) then
      call computeNe( newmesh, mesh, dual, dt, InvCellVol, GainsE, LossesE, diffusionE, &
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE%N, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW%N, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN%N, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS%N )
    else !TODO: This is redondant with copy_mesh operation at the begining of the temporal loop
      newmesh%Ne(:,:)=mesh%Ne(:,:)
      newmesh%Nh(:,:)=mesh%Nh(:,:) !TODO: Why this is updated  here? This should go with HolesOff
    end if
    !
    if(Params%HolesOff.eq.0 .AND. Params%NeOff.eq.0) then
      call computeNh( newmesh, mesh, dual, dt, InvCellVol, GainsH, LossesH, diffusionH, &
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE%N, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW%N, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN%N, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS%N )
    endif
    !
    !
    !
    if(Params%TeOff.ne.1) then
      !
      if(ConvectionEnergy.eq.0) then
        !
        call computeTe( newmesh, mesh, dual, dt, InvCellVol, kappae,  CouplingE, SourceE, invCe, &
                    ShapeFactorNormalE, ShapeFactorTangentE, NormalE%N,                          &
                    ShapeFactorNormalW, ShapeFactorTangentW, NormalW%N,                          &
                    ShapeFactorNormalN, ShapeFactorTangentN, NormalN%N,                          &
                    ShapeFactorNormalS, ShapeFactorTangentS, NormalS%N )
        !
      else
        !
        call computeUe( mesh, dt, InvCellVol, kappae,  CouplingE, SourceUe,                           &
                    invCe, Ue, UeNew, VeX, VeY, CellVol,                                              &
                    ShapeFactorNormalE, ShapeFactorTangentE, ShapeFactorNormalW, ShapeFactorTangentW, &
                    ShapeFactorNormalN, ShapeFactorTangentN, ShapeFactorNormalS, ShapeFactorTangentS, &
                    CellAreaE, CellAreaW, CellAreaN, CellAreaS,                                       &
                    NormalN, NormalS, NormalE, NormalW  )
        !
        !TODO: Should probably not be here
        call computeUh( mesh, dt, InvCellVol, kappah,  CouplingH, SourceUh,                           &
                    invCh, Uh, UhNew, VhX, VhY, CellVol,                                              &
                    ShapeFactorNormalE, ShapeFactorTangentE, ShapeFactorNormalW, ShapeFactorTangentW, &
                    ShapeFactorNormalN, ShapeFactorTangentN, ShapeFactorNormalS, ShapeFactorTangentS, &
                    CellAreaE,CellAreaW, CellAreaN,CellAreaS,                                         &
                    NormalN, NormalS, NormalE, NormalW  )
        !
      endif
      !
      if(Params%HolesOff.eq.0) then
        !
        if(ConvectionEnergy.eq.0) then
          !
          call computeTh( newmesh, mesh, dual, dt, InvCellVol, kappah,  CouplingH, SourceH, invCh,&
                    ShapeFactorNormalE, ShapeFactorTangentE, NormalE%N, &
                    ShapeFactorNormalW, ShapeFactorTangentW, NormalW%N, &
                    ShapeFactorNormalN, ShapeFactorTangentN, NormalN%N, &
                    ShapeFactorNormalS, ShapeFactorTangentS, NormalS%N )
          !
        else
          !
          !TODO: Check that we need that and not computeUh
          call computeUh_alt( mesh, dt, InvCellVol, kappah,  CouplingH, SourceUh, &
                    Ch, Uh, UhNew, VhX, VhY, CellVol, x, y, &
                    CellAreaE, NormalE%x, NormalE%y, CellAreaW, NormalW%x, NormalW%y, &
                    CellAreaN, NormalN%x, NormalN%y, CellAreaS, NormalS%x, NormalS%y  )
          !
        endif
        !
      endif
      !
    endif
    !
    if(Params%TsOff.ne.1) then
      !
      call computeTs( newmesh, mesh, dual, dt, InvCellVol, kappas, CouplingH, CouplingE, &
                    h1, h2, h3, invCs, TsPrev, TsOld, CellVol, &
                    ShapeFactorNormalE, ShapeFactorTangentE, NormalE%N, &
                    ShapeFactorNormalW, ShapeFactorTangentW, NormalW%N, &
                    ShapeFactorNormalN, ShapeFactorTangentN, NormalN%N, &
                    ShapeFactorNormalS, ShapeFactorTangentS, NormalS%N )
      !
    end if
    !
    !$OMP DO COLLAPSE(2) !(optimized)
    do j=2, Params%N-1 !(optimized)
      do i=2, Params%M-1 !(optimized)

      if(ConvectionEnergy.eq.1) then !define temperatures from energy
        newmesh%Te(i,j) = mesh%Te(i,j) + ((UeNew(i,j) -  Ue(i,j))-1.5d0*kb*mesh%Te(i,j)*(newmesh%Ne(i,j) - mesh%Ne(i,j)) &
            *FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)) ) * invCe(i,j)
        newmesh%Th(i,j) = mesh%Th(i,j) + ((UhNew(i,j) -  Uh(i,j))-1.5d0*kb*mesh%Th(i,j)*(newmesh%Nh(i,j) - mesh%Nh(i,j)) &
            *FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)) ) * invCh(i,j)
      end if
                      

      !TODO: Optimise
      CFLxT(i,j)=kappae(i,j)*invCe(i,j) * dt/(0.5d0*(DistW(i,j)+DistE(i,j)))**2
      CFLyT(i,j)=kappae(i,j)*invCe(i,j) * dt/(0.5d0*(DistN(i,j)+DistS(i,j)))**2
      CFLxTs(i,j)=kappas(i,j)*invCs(i,j) * dt/(0.5d0*(DistW(i,j)+DistE(i,j)))**2
      CFLyTs(i,j)=kappas(i,j)*invCs(i,j) * dt/(0.5d0*(DistN(i,j)+DistS(i,j)))**2
        
      !TODO: Optimise
      CFLxN(i,j)=diffusionE(i,j)*dt/(x(i,j)-x(i-1,j))**2 !+dt/(x(i,j)-x(i-1,j))*mobilityE(i,j)*sqrt(Ex(i,j)**2+Ey(i,j)**2)
      CFLyN(i,j)=diffusionE(i,j)*dt/(y(i,j)-y(i,j-1))**2 !+dt/(y(i,j)-y(i,j-1))*mobilityE(i,j)*sqrt(Ex(i,j)**2+Ey(i,j)**2)
        
      !TODO:This is only needed for a reduction, so lets do the reduction directly here
      TotalElectrons(i,j)=newmesh%Ne(i,j)*(0.125d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1)) &
                      -0.125d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1)))
      TotalHoles(i,j)=newmesh%Nh(i,j)*(0.125d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1)) &
                      -0.125d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1)))
        
        
      ThermalEnergy(i,j)=Ce(i,j)*mesh%Te(i,j)+Ch(i,j)*mesh%Th(i,j)+Cs(i,j)*mesh%Ts(i,j)


      work = intensity(i,j)/(1d0-reflectivity(i,j))
      LaserEnergy(i,j) = OnePhotonIonizationRate0 * work    & !energy loss by interband absorption
                       + TwoPhotonIonizationRate0 * work**2 & !energy loss by two photon absorption
                  + (absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*work !energy loss by carrrier heating
        
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of parallel section
    
    !BOUNDARY CONDITIONS
    
    do i=1, Params%M !North and South boundaries
      ! finite differences finite difference fashion
      if(DriftOn.eq.0) then
        newmesh%Ne(i,1)=newmesh%Ne(i,2)
        newmesh%Nh(i,1)=newmesh%Nh(i,2)
        newmesh%Ne(i,Params%N)=newmesh%Ne(i,Params%N-1)
        newmesh%Nh(i,Params%N)=newmesh%Nh(i,Params%N-1)
      end if

      UeNew(i,1)=UeNew(i,2)
      UhNew(i,1)=UhNew(i,2)
      newmesh%Te(i,1)=newmesh%Te(i,2)
      newmesh%Th(i,1)=newmesh%Th(i,2)
      newmesh%Ts(i,1)=newmesh%Ts(i,2)

      UeNew(i,Params%N)=UeNew(i,Params%N-1)
      UhNew(i,Params%N)=UhNew(i,Params%N-1)
      newmesh%Te(i,Params%N)=newmesh%Te(i,Params%N-1)
      newmesh%Th(i,Params%N)=newmesh%Th(i,Params%N-1)
      newmesh%Ts(i,Params%N)=newmesh%Ts(i,Params%N-1)
        ! includes also the corners... WHy are not they written?
        
!         potential(i,1)=0d0 !(0d0,0d0)
!                potential(i,N)=0d0 !(0d0,0d0)
    end do
    
    !$OMP PARALLEL DO  DEFAULT(NONE) SHARED(GradNeX, GradNeY, Params)
    do i=2,Params%M-1
      ! boundary condition v.n = 0 on boundaries.
      ! NORTH
      GradNeX(i,Params%N) = GradNeX(i,Params%N-1)
      GradNeY(i,Params%N) = GradNeY(i,Params%N-1)
      !
      ! SOUTH
      GradNeX(i,1) = GradNeX(i,2)
      GradNeY(i,1) = GradNeY(i,2)
    end do
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO  DEFAULT(NONE) SHARED(newmesh, Params, UeNew, UhNew, GradNeX, GradNeY)
    do j=2, Params%N-1 !West and East boundaries

      UeNew(1,j)=UeNew(2,j)
      UhNew(1,j)=UhNew(2,j)
      newmesh%Te(1,j)=newmesh%Te(2,j)
      newmesh%Th(1,j)=newmesh%Th(2,j)
      newmesh%Ts(1,j)=newmesh%Ts(2,j)

      newmesh%Te(Params%M,j)=newmesh%Te(Params%M-1,j) !Tout
      newmesh%Th(Params%M,j)=newmesh%Th(Params%M-1,j) !Tout
      newmesh%Ts(Params%M,j)=newmesh%Ts(Params%M-1,j) ! Tout !cooling by diffusion from outside, TsNew(M-1,j)

  !       potential(1,j)=0d0 !(0d0, 0d0)
  !       potential(M,j)=potential0 !(potential0, 0d0)
  
      if(DriftOn.eq.0) then
        newmesh%Ne(1,j)=newmesh%Ne(2,j)
        newmesh%Nh(1,j)=newmesh%Nh(2,j)
        ! conditions on the cone base - most important
        newmesh%Ne(Params%M,j)=newmesh%Ne(Params%M-1,j) !Ne0 replace by dynamic bnd condition with flux equal to the one of cell(M-1,j)
        newmesh%Nh(Params%M,j)=newmesh%Nh(Params%M-1,j) !Nh0 replace by dynamic bnd condition with flux equal to the one of cell(M-1,j)
      end if

      !TODO: Could we clean up these comments?

        !outlet condition on density and energy
!         write(*,*) CellAreaE(M,j), DistE(M-2,j)
!         NeNew(M,j) = -0.5d0*(diffusionE(M-2,j)+diffusionE(M-1,j))*(Ne(M-1,j)-Ne(M-2,j))/DistE(M-2,j)/(-0.5d0*diffusionE(M-1,j)-0.5d0*diffusionE(M,j))/DistE(M-1,j)+Ne(M-1,j)
!         NhNew(M,j) = -0.5d0*(diffusionH(M-2,j)+diffusionH(M-1,j))*(Nh(M-1,j)-Nh(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*diffusionH(M-1,j)-0.5d0*diffusionH(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Nh(M-1,j)
!         UeNew(M,j)=UeNew(M-1,j)
!         UhNew(M,j)=UhNew(M-1,j)
        

!         TeNew(M,j) = -0.5d0*(kappae(M-2,j)+kappae(M-1,j))*(Te(M-1,j)-Te(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappae(M-1,j)-0.5d0*kappae(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Te(M-1,j)
!         ThNew(M,j) = -0.5d0*(kappah(M-2,j)+kappah(M-1,j))*(Th(M-1,j)-Th(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappah(M-1,j)-0.5d0*kappah(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Th(M-1,j)
!         TsNew(M,j) = -0.5d0*(kappas(M-2,j)+kappas(M-1,j))*(Ts(M-1,j)-Ts(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappas(M-1,j)-0.5d0*kappas(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Ts(M-1,j)
      !
      ! WEST
      GradNeX(1,j) = GradNeX(2,j)
      GradNeY(1,j) = GradNeY(2,j)
      !
      ! EAST
      GradNeX(Params%M,j) = GradNeX(Params%M-1,j)
      GradNeY(Params%M,j) = GradNeY(Params%M-1,j)
    end do
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CHECKING the results

!     ElectronEnergy=0d0
!     HoleEnergy=0d0
!     LatticeEnergy=0d0
    
    maxCFLxN   = maxval(CFLxN)
    maxCFLyN   = maxval(CFLyN)
    maxCFLxT   = maxval(CFLxT)
    maxCFLyT   = maxval(CFLyT)
    maxCFLxTs  = maxval(CFLxTs)
    maxCFLyTs  = maxval(CFLyTs)
    maxTe      = maxval(newmesh%Te)
    minTe      = minval(newmesh%Te)
    maxTh      = maxval(newmesh%Th)
    minTh      = minval(newmesh%Th)
    maxTs      = maxval(newmesh%Ts)
    minTs      = minval(newmesh%Ts)
    maxNe      = maxval(newmesh%Ne)
    minNe      = minval(newmesh%Ne)
    maxNh      = maxval(newmesh%Nh)
    minNh      = minval(newmesh%Nh)
    maxIntensity = maxval(intensity)
    maxSourceE = maxval(SourceE)
    maxGainsE  = maxval(GainsE)
    maxSourceH = maxval(SourceH)
    maxGainsH  = maxval(GainsH)
    maxGap     = maxval(Egap)
    maxDiffNe  = maxval(diffNe)
    maxDiffNh  = maxval(diffNh)

    maxFermiIndexE = maxval(FermiIndexE)
    maxFermiIndexH = maxval(FermiIndexH)


    TotalMeshVolume=0d0
    !$OMP PARALLEL DEFAULT(NONE) SHARED(mesh, CellVol, NeTotal,NhTotal,TotalNumOfE,  &
    !$OMP TotalThermalEnergy, TotalLaserEnergy, TotalMeshVolume, TotalNumOfH, TotalElectrons, &
    !$OMP TotalHoles, ThermalEnergy, LaserEnergy)
    !$OMP DO COLLAPSE(2) REDUCTION(+:NeTotal,NhTotal,TotalNumOfE, TotalNumOfH,  &
    !$OMP TotalThermalEnergy, TotalLaserEnergy, TotalMeshVolume)
    do j=1,mesh%N
      do i=1,mesh%M
        NeTotal=NeTotal + mesh%Ne(i,j) * CellVol(i,j)
        NhTotal=NhTotal + mesh%Nh(i,j) * CellVol(i,j)
        TotalNumOfE=TotalNumOfE + TotalElectrons(i,j)
        TotalNumOfH=TotalNumOfH + TotalHoles(i,j)
        TotalThermalEnergy=TotalThermalEnergy+ThermalEnergy(i,j)
        TotalLaserEnergy=TotalLaserEnergy+LaserEnergy(i,j)
        TotalMeshVolume=TotalMeshVolume+CellVol(i,j)
     end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL

    work = I0*exp(-.5d0*((t-t0)/sigmaTau)**2)*dt
    !$OMP PARALLEL DEFAULT(NONE) SHARED(Params,OnePhotonIonizationRate0, absorptionDrudeE, &
    !$OMP absorptionDrudeH, OpticalIndex, CellVol, work, LaserIntensityEnergy )
    !$OMP DO COLLAPSE(2) REDUCTION(+:LaserIntensityEnergy)
    do j=1,Params%N
      do i=1,Params%M
           LaserIntensityEnergy = LaserIntensityEnergy &
               +(OnePhotonIonizationRate0 +absorptionDrudeE(i,j) +absorptionDrudeH(i,j))*OpticalIndex(i,j) &
                  *CellVol(i,j)*work
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    !TODO: Should be parallelized
    do j=1,Params%N
      do i=1,Params%M


        if(MaxHeating(i,j) < mesh%Ts(i,j) .AND. t > 100d0*laser%tau) then
           MaxHeating(i,j)=mesh%Ts(i,j)
           MaxHeatingTime(i,j)=t
        end if
        
!         if(ConductivityFix.eq.-1) then
!            TeNew(i,j)=Tout; ThNew(i,j)=Tout; 
!         end if

!         CAUTION: These definitions are erroneously including initial temperature into account. 
!         ElectronEnergy=ElectronEnergy+Ce(i,j)*Te(i,j)*CellVol(i,j)
!         HoleEnergy=HoleEnergy+Ch(i,j)*Th(i,j)*CellVol(i,j)
!         LatticeEnergy=LatticeEnergy+Cs(i,j)*Ts(i,j)*CellVol(i,j)

        ! calculation of the absorbed laser energy involved in the simulated slice !
        if((i.eq.1) .AND. (j.eq.(Params%N/2))) then
          IntensityEnergy=IntensityEnergy+(OnePhotonIonizationRate0+absorptionDrudeE(i,j) &
                  +absorptionDrudeH(i,j))*intensity(i,j)*CellVol(i,j)*dt
        end if

        work = EgapValue(newmesh%Ne(i,j),newmesh%Ts(i,j)) - Egap(i,j)
        ! calculation of the energy contained in the solid
        ElectronEnergy=ElectronEnergy &
          + (Ce(i,j)*(newmesh%Te(i,j)-mesh%Te(i,j))+(Ce(i,j)-CeOld(i,j))*mesh%Te(i,j)) * CellVol(i,j) & !kinetic energy
          + (mesh%Ne(i,j)*work  &
              + Egap(i,j)*(newmesh%Ne(i,j)-mesh%Ne(i,j))) * CellVol(i,j) !potential energy
              
        ElectronKineticEnergy=ElectronKineticEnergy &
                   +(Ce(i,j)*(newmesh%Te(i,j)-mesh%Te(i,j))+(Ce(i,j)-CeOld(i,j))*mesh%Te(i,j)) * CellVol(i,j) !kinetic energy
        ElectronPotentialEnergy=ElectronPotentialEnergy &
                   +(mesh%Ne(i,j)*work  &
                   + Egap(i,j)*(newmesh%Ne(i,j)-mesh%Ne(i,j))) * CellVol(i,j) !potential energy
!         ElectronEnergy=ElectronKineticEnergy+ElectronPotentialEnergy !already summed over time
        
        HoleEnergy=HoleEnergy+(Ch(i,j)*(newmesh%Th(i,j)-mesh%Th(i,j))+(Ch(i,j)-ChOld(i,j))*mesh%Th(i,j)) * CellVol(i,j) !kinetic energy
        
        LatticeEnergy=LatticeEnergy+((Cs(i,j)*(newmesh%Ts(i,j)-mesh%Ts(i,j)))+0d0*(Cs(i,j)-CsOld(i,j))*mesh%Ts(i,j))*CellVol(i,j) !dCs/dt=0, 20150426, TJYD.

     end do
   end do


   call check_divergences(mesh, maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN, x, y, t )


   do j=1,Params%N
     do i=1,Params%M

        if(real(FermiIndexE(i,j)) > real(FermiMaxLines) .OR. real(FermiIndexE(i,j)) < 1d0) then
          write(*,*) "t,i,j,FermiIndexE(i,j)=", t,i,j,FermiIndexE(i,j)
        end if
        if(real(FermiIndexH(i,j)) > real(FermiMaxLines) .OR. real(FermiIndexH(i,j)) < 1d0) then
          write(*,*) "t,i,j,FermiIndexH(i,j)=", t,i,j,FermiIndexH(i,j)
        end if
        
        if(real(FermiRatioE(i,j)) < 0d0 .OR. real(FermiRatioH(i,j)) < 0d0) then
          write(*,*) "Problem in DOS or Ne. DOS(i,j)=", i,j,DOSe(i,j), DOSh(i,j), "Ne,h(i,j)=", mesh%Ne(i,j), mesh%Nh(i,j)
        end if
        
     end do
   end do

   call flush(Error%unit)



   do i=1,Params%M
     do j=1,Params%N
              ! lets change dt when fast reponse is finished in order to catch the long one. 

        diffNe(i,j)=(newmesh%Ne(i,j)-mesh%Ne(i,j))/dt
        diffNh(i,j)=(newmesh%Nh(i,j)-mesh%Nh(i,j))/dt
        
                
        if(mod(nbiter,iterOut*iterOutMaps).eq.0) then 
          write(Depth%unit,887, advance="yes") t, x(i,j), y(i,j), intensity(i,j), mesh%Te(i,j), & !5
                        mesh%Th(i,j), mesh%Ts(i,j), mesh%Ne(i,j), mesh%Nh(i,j), reflectivity(i,j), & !10
                        absorptionDrudeE(i,j), absorptionDrudeH(i,j), diffNe(i,j), diffNh(i,j), TotalElectrons(i,j), & !15
                        TotalHoles(i,j), real(FermiIndexE(i,j)), REAL(FermiIndexH(i,j)), FermiRatioE(i,j), FermiRatioH(i,j), & !20
                        SourceE(i,j), SourceH(i,j), GainsE(i,j), GainsH(i,j), LossesE(i,j), & !25
                        LossesH(i,j), real(DielectricDrudeE(i,j)), aimag(DielectricDrudeE(i,j)), Egap(i,j), real(Dielectric(i,j)), & !30
                        aimag(Dielectric(i,j)), MaxHeatingTime(i,j), MaxHeating(i,j), real(potentialNeedle(i,j)), Ex(i,j), & !35
                        Ey(i,j), diffusionE(i,j), diffusionH(i,j), GradNeX(i,j), GradNeY(i,j), &!40
                        real(EintField(i,j)), aimag(EintField(i,j)), EintFieldR(i,j), EintFieldI(i,j), phiMie(i,j), & !45
                        Radius(i,j)

      !TODO: Please use short notation with prenthesis !!
  887 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
  3x, 1E12.5)        
!          write(Depth%unit,886, advance='yes')
        end if
        
      end do !on Y
      
       if(mod(nbiter,iterOut*iterOutMaps).eq.0) then 
        write(Depth%unit,886, advance="yes")
886        FORMAT (3x)
       end if
    end do !on X

    if(mod(nbiter,iterOut*iterOutMaps).eq.0) &
      call flush(Depth%unit);


    ! saving the timesteps of several previous steps (used for the high order calculation of d/dt).
    dt4=dt3;
    dt3=dt2;
    dt2=dt; 
    ! chaning the timestep based on known behavior of the system
    if(AdaptativeTimeStep.eq.1) then
      if((t>1d1*laser%tau*coeffDilaDt) .AND. (dt.eq.dt0) .AND. &
        (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
        dt=10d0*dt0
      else if ((t > 0.5d2*laser%tau*coeffDilaDt) .AND. (dt.eq.10d0*dt0) .AND. (maxCFLxN+maxCFLyN &
        + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
        dt=40d0*dt0
      else if ((t > 1d3*laser%tau*coeffDilaDt) .AND. (dt.eq.40d0*dt0) .AND. (maxCFLxN+maxCFLyN + maxCFLxT &
        + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
          dt=1d3*dt0
!           else if (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs > maxCFL) then
!             dt=dt/1d1
      end if
    end if
        
    
    ! writing result for the electrostatic calculations
    if(mod(nbiter,iterOut*iterOutMaps).eq.0) then 
      do i=1,Mp
        do j=1,Np
          ! ecriture des donnees dans un fichier different
          write(DepthVessel%unit,889, advance="yes") t, xP(i,j), yP(i,j), real(potential(i,j)), real(ExPoisson(i,j)), & !
                real(EyPoisson(i,j)), DielectricStatic(i,j), NeP(i,j), NhP(i,j)
  889 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, & !TODO: Please use short notation with prenthesis !!
  1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5) 
        end do
      end do
      call flush(DepthVessel%unit);
      
!       write(91,'(100E14.5)') potential !Xvector !Amatrix
!     write(91,*)
      
    end if
       
    ! output to files
    
    cpu_timestep_duration = ElapsedTime() / real(nbiter)
    cpuefficiency=real(nthreads)/cpu_timestep_duration
    
    !This should be moved to output.F90 file
    if(mod(nbiter,iterOut).eq.0) then 
      
      write(EnergyConservation%unit,892, advance="YES") t, IntensityEnergy, ElectronEnergy, HoleEnergy, LatticeEnergy, & !5
          TotalMeshVolume, LaserIntensityEnergy, ElectronKineticEnergy, ElectronPotentialEnergy !9
      
892 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)

      write(TimeMax%unit,888, advance="YES") t, maxTe, maxTh, maxTs, maxNe, &         !5
                    maxNh, maxIntensity, TotalLaserEnergy, TotalThermalEnergy, &        !9
                    maxSourceE, maxGainsE, maxSourceH, maxGainsH, maxGap, &        !14
                    maxDiffNe, maxDiffNh, TotalNumOfE, TotalNumOfH, real(maxFermiIndexE), &        !19
                    real(maxFermiIndexH), NeTotal, NhTotal, maxCFLxT, maxCFLyT, &        !24
                    maxCFLxN, maxCFLyN, maxCFLxTs, maxCFLyTs, IntensityEnergy, &        !29
                    TotalMeshVolume, ElectronEnergy, HoleEnergy, LatticeEnergy, LaserIntensityEnergy, &         !34
                    ElectronKineticEnergy, ElectronPotentialEnergy    !36
                    
888 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1F12.8, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E19.11, 3x, 1E19.11, 3x, 1E19.11, 3x, 1E19.11, 3x, 1E19.11, &
3x, 1E19.11, 3x, 1E19.11)
                
        write(TimeApex%unit,884, advance="YES") t, mesh%Te(1,Params%N/2), mesh%Th(1,Params%N/2), &
              mesh%Ts(1,Params%N/2), mesh%Ne(1,Params%N/2), &                        !5
              mesh%Nh(1,Params%N/2), intensity(1,Params%N/2), TotalLaserEnergy, TotalThermalEnergy, &        !9
              SourceE(1,Params%N/2), GainsE(1,Params%N/2), SourceH(1,Params%N/2), GainsH(1,Params%N/2), Egap(1,Params%N/2), &                !14
              diffNe(1,Params%N/2), diffNh(1,Params%N/2), real(FermiIndexE(1,Params%N/2)),&
               real(FermiIndexH(1,Params%N/2)), Ce(2,Params%N/2), &                !19
              CeOld(2,Params%N/2), Ch(2,Params%N/2), ChOld(2,Params%N/2), Cs(2,Params%N/2), CsOld(2,Params%N/2)                               !24
              
884 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)

        write(TimeUp%unit,883, advance="YES") t, mesh%Te(Params%M/2,Params%N), mesh%Th(Params%M/2,Params%N),&
                          mesh%Ts(Params%M/2,Params%N), mesh%Ne(Params%M/2,Params%N), &
              mesh%Nh(Params%M/2,Params%N), intensity(Params%M/2,Params%N), TotalLaserEnergy, TotalThermalEnergy, &
              SourceE(Params%M/2,Params%N), GainsE(Params%M/2,Params%N), SourceH(Params%M/2,Params%N),&
               GainsH(Params%M/2,Params%N), Egap(Params%M/2,Params%N), &
              diffNe(Params%M/2,Params%N), diffNh(Params%M/2,Params%N), &
              real(FermiIndexE(Params%M/2,Params%N)), real(FermiIndexH(Params%M/2,Params%N))
              
883 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)


        write(TimeBottom%unit,882, advance="YES") t, mesh%Te(Params%M/2,1), mesh%Th(Params%M/2,1), mesh%Ne(Params%M/2,1), &
              mesh%Nh(Params%M/2,1), intensity(Params%M/2,1), TotalLaserEnergy, TotalThermalEnergy, &
              SourceE(Params%M/2,1), GainsE(Params%M/2,1), SourceH(Params%M/2,1), GainsH(Params%M/2,1), Egap(Params%M/2,1), &
              diffNe(Params%M/2,1), diffNh(Params%M/2,1), real(FermiIndexE(Params%M/2,1)), real(FermiIndexH(Params%M/2,1))
              
882 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)

    end if
    

    
890 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
    
    ! write the functions on Dual Mesh
    if(mod(nbiter,iterOut*iterOutMaps).eq.0) then
      do i=1,Params%M-1
        do j=1,Params%N-1
          
          write(DualDepth%unit, 890, advance="YES") t, xDual(i,j), yDual(i,j), dual%Te(i,j), dual%Th(i,j), & !5
                                          dual%Ts(i,j), dual%Ne(i,j), dual%Nh(i,j), intensityDual(i,j) !9
                  
        end do
      end do
      call flush(TimeMax%unit); call flush(TimeApex%unit); call flush(TimeUp%unit); call flush(DepthVessel%unit)
    end if


  end do !end of time loop

  
  call releasemesh(mesh)
  call releasemesh(dual)
  call releasemesh(newmesh)

  deallocate(NormalN%x, NormalN%y, NormalN%N)
  deallocate(NormalS%x, NormalS%y, NormalS%N)
  deallocate(NormalE%x, NormalE%y, NormalE%N)
  deallocate(NormalW%x, NormalW%y, NormalW%N)

  !TODO: Sorry but where are the file stream closed???

  call ReleaseInputParameters( Params )

end program Flaps

