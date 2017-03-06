!! Copyright (C) 2012-2016 T. J.-Y. Derrien, N. Tancogne-Dejean
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
!------------------------------------------------------------------------------

program Flaps
  USE OMP_LIB
  ! include 'Bivariate.f'
  ! USE Bivariate
  USE libmsh2vf !Script provided by A. Mouton, Univ Lille1, France for GMSH interfacing

  use Laser_m
  use Material_m
  use Maths_m
  use Mie_m
  use Output_m
  use Profiler_m
  use Restart_m
  use Timer_m
  use Types_m


  implicit none

    type(MeshValues) :: mesh, dual, newmesh
    type(Laser)      :: source
    type(InputParameters) :: Params
    type(Material)   :: matter
    type(Timer)      :: full_timer

    real(8), parameter:: potential0=7d3,&         ! potential at the bottom of the needle ; default = 7d3
                         potentialNull=M_ZERO !, &
 !                       phiMie0=1d0*acos(-1d0)                ! Mie scattering: plane angle in cylindrical coordinates
    
    real(8), parameter:: coeffDilaDt=M_TWO        ,& !diltation coeff before dt change
                        xmin=-10d-6       ,& !mesh min
                        xmax=10d-6       ,& !mesh max
                        ymin=-10d-6       ,&                
                        ymax=10d-6       ,&
                        tCenter=M_ZERO             !time of gaussian intensity maximum

     real(8)               tmin                !max absolute time
    
                        
    integer(8), parameter:: iterOutMaps=1000      ,& ! number of outputs for maps between each stdout
                          VirtualPoints=3, & !number of virtual points to exclude from the GMSH file (locate them at the beginning!)
                          Mv=101       ,& !number of celles in the Vessel domain (larger) X direction
                          Nv=101        ,& !number of celles in the Vessel domain (larger) Y direction
                          MeshChoice=1       ,& !0: rectangle (xmin,xmax)(ymin,ymax). 1: Experimental cones, 2: Cone in a vessel (HS), 3: import GMSH (working)
                          MeshIterations=500000        ,&        !number of iterations to calculate meshNeedle
                          MeshIterationsVessel=100*Mv,&        !number of iterations to calculate meshVessel
                          MeshShift=1       ,&         !number of cells x N in the tip, 343 nm: 2; 515 nm: 3;
                          FermiMaxLines=3584            ! >= number of lines in Fermi file
       !                   SORiterations=1        ,&        !iteration number for over-relaxation method
       !                   InterpolateMethod=1        ,&        ! 0: linear, 1: bicubic
       !                   UseInterpolation=0        ,&
       !                   SolveImplicit=1, &        ! 0: use explicit schemes, 1: use implicit scheme (band diagonal matrixes)
       !                   numberOfNeighbours=4
                          
                          
        
    real(8), parameter::   NeedleAngleDeg=M_ONE        ,& ! deg
                          NeedleRadius=4d-9        ,& !m
                          NeedleLength=3d-6        , &        !m
                          SORcoeff=1.2d0        ,&        ! near 1
                          MeshDensity0=5d5        ,& ! amplitude of source for mesh refinement
                          MeshConvergenceEpsilon=1d-10, & !error tolerance on meshing convergency
                          MeshDamping0=1d8, &         ! damping coefficient (m^-1) for mesh refinement
                          ActivateInduction=M_ZERO, &
                          CrossCoeff=-M_ONE, &                ! 0d0: OFF, 1d0: ON
                          maxCFL=1d-3                        ! maximum admitted on CFL condition for any time step increase
                          
    integer(8), parameter:: PoissonOn=0       ,& !0: Poisson solver is OFF. 1: Calculation of potential ON.
                            DriftOn=0       ,& !0: Drift is disabled. 1: Enabled. 
                            CathodeZone=1        ,& !1: on the needle bottom, 0: on back vessel (not physical but stable)
                            ShiftFixedPotential=-1, &        ! while CathodeZone=1, use to adjust the number of points on which tension is applied
                            PoissonSolver=0        ,& !0: Full matrix inversion once, 1: SOR iterative for each dt
                            InterpolateOff=0,         &        !just to test speedup...
                            BandBendingInFDTD=0        ,&        !use the interpolation of FDTD 1030 nm with band-bending contribution
!                            PolarizationSource=0, &        ! 0: source TE, 1: source TM
                            !MieScattering=1, &
                            NewtonIterations=1000, &
                            ExpNeedleType=0
        
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
    
    integer(8)         nbiter, i, j, k, nmin, nmax, NeedleIndexX, NeedleIndexY, maxFermiIndexE, maxFermiIndexH, &
                Mp, Np, RunningIndex
    real(8)         t, t0, dx, dy, x0, y0, dt, dt2, dt3, dt4, h1, h2, h3
    real(8)         I0 !initial values of the problem

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
                Egap, Egap_new, &                        ! local gap value
                SourceE, SourceH, & ! heating sources
                SourceUe, SourceUh, & ! free carrier thermal energy sources
                x, y, &                 ! needle position indexes
                xNew, yNew, &           ! used to converge the mesh parallely
                xV, yV, &                 ! vessel position indexes
                xDualSW, yDualSW, &                 ! dual mesh position
                xDualSE, yDualSE, &
                xDualNE, yDualNE, &
                xDualNW, yDualNW, &
                xDual, yDual, &
                CFLxT, CFLyT, CFLxN, CFLyN, CFLxTs, CFLyTs, &
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
                
                 !TODO: Use dimension. TJYD: What do you have in mind? Example? 
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
    
    real(8), allocatable, target :: FermiTableE(:,:),&
                                     FermiTableH(:,:),& !reduced Fermi level for electrons and holes
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

   
    real(8) sigmaTau, sigmaX, sigmaY, &
            maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN, maxCFLxTs, maxCFLyTs, &
            maxTe, minTe, maxTh, minTh, maxTs, minTs, maxIntensity, maxNe, minNe, maxNh, minNh, &
            maxSourceE, maxGainsE, maxSourceH, maxGainsH, maxGap, maxDiffNe, maxDiffNh, &
            TotalLaserEnergy, TotalThermalEnergy, ElectronPotentialEnergy, ElectronKineticEnergy, &
            NeedleHeight, NeedleA, NeedleB, Needlet0Limit, NeedleXParam, NeedleYParam, NeedleAngle, &
            localT, P2critic, TotalNumOfE, TotalNumOfH, &
            xmin2, xmax2, ymin2, ymax2, &
            MeshConvergence, MeshConvergenceOld
    real :: cpuefficiency, cpu_timestep_duration
            
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
     real(8) ConeExp1, ConeExp2, DensityOfState, &
             LatticeHeatCapacity
     complex(8) DielectricConstant ! DielectricFunction, DielectricFunctionDrude,

            
!OPENMP declarations
    integer :: myid, nthreads
!#IFDEF (OMP_NUM_THREADS)
!    integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!#ENDIF

  !TODO: move to material.f90
  !TODO: Should be a parameter, to guaranty no modification
  !TODO: TJYD: this should be imposed by the model we gonna call in the continuum description library. 
  me=M_HALF*me0       ! electron effective mass for conductivity !0.24 (source ?)
  mh=M_HALF*me0       ! hole effective mass for conductivity !0.81 (source ?)
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

  call Profiler_start(prof_init, 'INIT')
  call timer_init(full_timer)

  call InitInputParameter( Params )
  call LoadInputParameters( "flaps.in", Params )
  call CheckValidityInputParameters( Params )

  call init_material( matter, Params%AugerOff )

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
                Egap(1:Params%M, 1:Params%N), Egap_new(1:Params%M, 1:Params%N), &                        ! local gap value
                SourceE(1:Params%M, 1:Params%N), SourceH(1:Params%M, 1:Params%N), & ! heating sources
                SourceUe(1:Params%M, 1:Params%N), SourceUh(1:Params%M, 1:Params%N), & ! free carrier thermal energy sources
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

 !Initialisation of data
 CouplingE(:,:)=M_ZERO
 CouplingH(:,:)=M_ZERO
 kappae(:,:) = M_ZERO
 kappah(:,:) = M_ZERO
 kappas(:,:) = M_ZERO
 JeX(:,:)=M_ZERO
 JeY(:,:)=M_ZERO
 JhX(:,:)=M_ZERO
 JhY(:,:)=M_ZERO
 diffusionE(:,:)=M_ZERO
 diffusionH(:,:)=M_ZERO
 absorptionDrudeE(:,:) = M_ZERO
 absorptionDrudeH(:,:) = M_ZERO

!******** READ GMSH MESH FILE ************
!TODO: This is very dirty
namefile_msh='external_libs/gmsh/mesh.msh'
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

  call init_laser(source, Params)

  tmin=tCenter-5d0*source%tau


  dt=Params%TimeStep
  dt2=Params%TimeStep
  dt3=Params%TimeStep
  dt4=Params%TimeStep

  ColFermiNeNc=2; ColFermiEta=3; ColFermi0=4; ColFermi1=5; ColFermi2=6; ColFermiHalf=7; 
  ColFermiThreeHalf=8; ColFermiMenusHalf=9;


    !TODO: Move to LaserParams
    sigmaTau=source%tau/(M_TWO*M_SQRT2LN2)
    sigmaX=source%spotX/(M_TWO*M_SQRT2LN2)
    sigmaY=source%spotY/(M_TWO*M_SQRT2LN2)

    !TODO: Move to LaserParams
    I0=source%fluence/source%tau * sqrt(4d0 * log(M_TWO) / M_PI)
    
    if(source%lambda.eq.515d-9) then
      if(Params%PolarizationSource.eq.0) then
        x1=1.5d-7; y1=M_ZERO; I1=M_ZERO*I0; spotX1=100d-9; spotY1=50d-9; !
        x2=3.3d-7; y2=M_ZERO; I2=M_ZERO*9d0*I0; spotX2=50d-9; spotY2=50d-9; !3.53W
        x3=6.5d-7; y3=-3d-8; I3=M_ZERO*17d0*I0; spotX3=70d-9; spotY3=150d-9; !28W
        x4=6.5d-7; y4=3d-8; I4=M_ZERO*17d0*I0; spotX4=70d-9; spotY4=150d-9; !28W
        x5=9.5d-7; y5=M_ZERO; I5=M_ZERO*7d0*I0; spotX5=50d-9; spotY5=50d-9; !2.74W
        x6=1.2d-6; y6=M_ZERO; I6=M_ZERO*8d0*I0; spotX6=50d-9; spotY6=50d-9; !3.14W
        x7=1.37d-6; y7=M_ZERO; I7=M_ZERO*9d0*I0; spotX7=50d-9; spotY7=50d-9; !3.53W
        x8=6.2d-7; y8=-40d-9; I8=M_ZERO*I0; spotX8=250d-9; spotY8=50d-9;
        x9=1.3d-6; y9=-20d-9; I9=M_ZERO*I0; spotX9=500d-9; spotY9=100d-9;
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
    sigmaX1=spotX1/(M_TWO*M_SQRT2LN2)
    sigmaY1=spotY1/(M_TWO*M_SQRT2LN2)
    sigmaX2=spotX2/(M_TWO*M_SQRT2LN2)
    sigmaY2=spotY2/(M_TWO*M_SQRT2LN2)
    sigmaX3=spotX3/(M_TWO*M_SQRT2LN2)
    sigmaY3=spotY3/(M_TWO*M_SQRT2LN2)
    sigmaX4=spotX4/(M_TWO*M_SQRT2LN2)
    sigmaY4=spotY4/(M_TWO*M_SQRT2LN2)
    sigmaX5=spotX5/(M_TWO*M_SQRT2LN2)
    sigmaY5=spotY5/(M_TWO*M_SQRT2LN2)
    sigmaX6=spotX6/(M_TWO*M_SQRT2LN2)
    sigmaY6=spotY6/(M_TWO*M_SQRT2LN2)
    sigmaX7=spotX7/(M_TWO*M_SQRT2LN2)
    sigmaY7=spotY7/(M_TWO*M_SQRT2LN2)
    sigmaX8=spotX8/(M_TWO*M_SQRT2LN2)
    sigmaY8=spotY8/(M_TWO*M_SQRT2LN2)
    sigmaX9=spotX9/(M_TWO*M_SQRT2LN2)
    sigmaY9=spotY9/(M_TWO*M_SQRT2LN2)
    
    t0=tCenter; x0=source%xCenter; y0=source%yCenter;
    
    ! We open the different files
    call InitOutputs(Params%RestartCalc)
    
    
  TotalLaserEnergy=M_ZERO;
  IntensityEnergy=M_ZERO;
  LaserIntensityEnergy=M_ZERO;
  ElectronKineticEnergy=M_ZERO; ElectronPotentialEnergy=M_ZERO; ElectronEnergy=M_ZERO
  HoleEnergy=M_ZERO
  LatticeEnergy=M_ZERO
  TotalThermalEnergy=M_ZERO;
  cpuefficiency=M_ZERO; cpu_timestep_duration=M_ZERO
  
  
!***************** MESH GENERATION *****************
  write(*,*) "[Mesh] Building..."


  !Allocate the mesh and the dual mesh
  call initmesh(mesh, Params%M, Params%N)
  call initmesh(dual, Params%M-1, Params%N-1)
  call initmesh(newmesh, Params%M, Params%N)

  ! building rectangular mesh 
  dx=(xmax-xmin)/(Params%M+1)
  dy=(ymax-ymin)/(Params%N+1)
  x(1:Params%M,1:Params%N) = M_ZERO
  y(1:Params%M,1:Params%N) = M_ZERO
  
  write(*,*)

  call output_open(MeshInfo%unit,'output/mesh.dat', .false.)               ! format 885, 8852

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
    NeedleAngle=NeedleAngleDeg*M_DEG2RAD
    NeedleA=NeedleRadius/(tan(NeedleAngle*M_HALF)**2)
    NeedleB=NeedleRadius/tan(NeedleAngle*M_HALF)
    Needlet0Limit=acos(NeedleRadius/(NeedleLength*tan(NeedleAngle*M_HALF)**2+NeedleRadius))
    
!building of the conical mesh: write boundaries, then solve laplace, and iterate
    
    NeedleXParam=(Needlet0Limit-M_ZERO)/NeedleIndexX !dt0 for X (0,t0)
    NeedleYParam=(Needlet0Limit-M_ZERO)/NeedleIndexY !dt0 for Y (0,t0)
    
    P2critic=(real(Params%N)-M_ONE)**2/(M_TWO*real(j)-real(Params%N)-1)
    !! this limit applies if hyperbolic contour is chosen
!     if(real(MeshShift)**2 > P2critic) then
!       write(*,*) "Dilatation of the apex is too large. Reduce <= ", floor(sqrt(P2critic))
!       stop
!     end if

!! just to check which point has not been defined - every point should be replaced !
    x(:,:)=-M_ONE
    y(:,:)=-M_ONE
    
    MeshConvergence=M_ONE
    
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
    NeedleHeight=2d0*tan(M_HALF*NeedleAngle) * sqrt(M_ONE - NeedleRadius**2/((tan(NeedleAngle*M_HALF))**4)* &
                (x(Params%M,Params%N/2)+NeedleRadius**2/((tan(M_HALF*NeedleAngle))**2)**2)) &
              * (x(Params%M,Params%N/2)+NeedleRadius/(tan(NeedleAngle*M_HALF))**2)
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
  close(MeshInfo%unit)

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
  
  call output_open(MeshVessel%unit,'output/meshVessel.dat', .false.)
  do i=1,Mp
    do j=1,Np
      write(MeshVessel%unit,885) xP(i,j), yP(i,j), i, j
    end do
    write(MeshVessel%unit,*) " "
  end do
  close(MeshVessel%unit)
      
  
  write (*,*) "Vessel cells:", Mp, "*", Np,"=", Mp*Np
  write (*,*) "Matter cells:", Params%M, "*", Params%N, "=", Params%M*Params%N
  
!   allocate(Amatrix(1:Mp*Np,1:Mp*Np))
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
  call output_open(Parameters%unit,'output/parameters.dat', .false.)
  write(Parameters%unit,*) "========== CONE PARAMETERS ========="
  write(Parameters%unit,*) "Cone length=", NeedleLength*1d6, "um"
  write(Parameters%unit,*) "Cone height=", NeedleHeight*1d6, "um"
  write(Parameters%unit,*) "Cone angle=", NeedleAngleDeg, "deg"
  write(Parameters%unit,*) "Cone curvature radius=", NeedleRadius*1d9, "nm"
  !TODO: Move to source
  write(Parameters%unit,*)
  write(Parameters%unit,*) "=========== LASER PARAMETERS ========"
  write(Parameters%unit,*) "Laser fluence=", source%fluence*1d-4, "J.cm-2"
  write(Parameters%unit,*) "Laser pulse duration=", source%tau*1d15, "fs"
  write(Parameters%unit,*) "Laser wavelength=", source%lambda*1d9, "nm"
  write(Parameters%unit,*) "Laser spot position: (X,Y)=", x0*1d6, y0*1d6, "um"
  write(Parameters%unit,*) "Laser spot size: (Sx, Sy)=", source%spotX*1d6, source%spotY*1d6, "um"
  write(Parameters%unit,*) "Mie scattering:", Params%UseMieScattering
  write(Parameters%unit,*) "Laser polarization", Params%PolarizationSource
  write(Parameters%unit,*)
  write(Parameters%unit,*) "============ MESH PARAMETERS =========="
  write(Parameters%unit,*) "Mesh size", Params%M, "x", Params%N
  write(Parameters%unit,*) "Dilatation time ratio=", coeffDilaDt
  write(Parameters%unit,*) "Mesh shift=", MeshShift
  write(Parameters%unit,*)
  write(Parameters%unit,*) "============ TIME CONTROL =========="
  write(Parameters%unit,*) "Initial timestep=", Params%TimeStep
  write(Parameters%unit,*) "Maximal timestep=", tmin
  write(Parameters%unit,*) "Maximum time t=", Params%TimeMax
  write(Parameters%unit,*) "Enable adaptative timestep=", Params%AdaptativeTimeStep
  write(Parameters%unit,*) "Time output each ", Params%OutputIter, "iterations."
  write(Parameters%unit,*) "Map output each", iterOutMaps*Params%OutputIter, "iterations."
  close(Parameters%unit)
  
  call timer_start(full_timer)

  nmax=int((Params%TimeMax-tmin)/dt, 8)
  
  Ex(:,:)=M_ZERO !-1d10
  Ey(:,:)=M_ZERO !-1d9 !0d0
  DummyVessel(:,:)=M_ZERO
  
  work = LatticeHeatCapacity(Params%Text)
   ! $ OMP DO
  do j=1,Params%N
    do i=1,Params%M
        newmesh%Te(i,j)=Params%Text
        newmesh%Th(i,j)=Params%Text
        newmesh%Ts(i,j)=Params%Text
        TsOld(i,j)=Params%Text
        TsPrev(i,j)=Params%Text

        CsOld(i,j)=work
        CsPrev(i,j)=work
        CsPrev2(i,j)=work
    end do
  end do
  ! $ OMP END DO

  ! $ OMP DO
  do j=1,Params%N
    do i=1,Params%M
        
        VeX(i,j)=M_ZERO
        VeY(i,j)=M_ZERO
        VhX(i,j)=M_ZERO
        VhY(i,j)=M_ZERO

        if(BandBendingInFDTD.eq.1) then
           newmesh%Ne(i,j)=Ne0+Nborder*(exp(-M_HALF*(((x(i,j)-x(i,Params%N))**2+(y(i,j)-y(i,Params%N))**2) &
                    /((DefectThickness)/(2d0*M_SQRT2LN2))**2)) &
                  +exp(-M_HALF*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness) &
                    /(M_TWO*M_SQRT2LN2))**2)) &
                  +exp(-M_HALF*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness) &
                    /(M_TWO*M_SQRT2LN2))**2)) &
                  )
          newmesh%Nh(i,j)=Nh0+Nborder*(exp(-M_HALF*(((x(i,j)-x(i,Params%N))**2+(y(i,j)-y(i,Params%N))**2) &
                    /((DefectThickness)/(2d0*M_SQRT2LN2))**2)) &
                  +exp(-M_HALF*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness) &
                    /(M_TWO*M_SQRT2LN2))**2)) &
                  +exp(-M_HALF*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness) &
                    /(M_TWO*M_SQRT2LN2))**2)) &
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
                  -etae*(M_ONE-(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j)) &
                              /FermiTableE(ColFermiHalf,FermiIndexE(i,j)))* &
                  (FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)) &
                              /FermiTableE(ColFermiHalf,FermiIndexE(i,j)))))
        ChOld(i,j)=1.5d0*mesh%Nh(i,j)*kb*(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j)) &
                              /FermiTableH(ColFermiHalf,FermiIndexH(i,j)) &
                  -etah*(M_ONE-(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j)) &
                              /FermiTableH(ColFermiHalf,FermiIndexH(i,j)))* &
                  FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)) &
                              /FermiTableH(ColFermiHalf,FermiIndexH(i,j))))

        Ce(i,j)=CeOld(i,j)
        Ch(i,j)=ChOld(i,j)
        Cs(i,j)=CsPrev(i,j)

        
        UeNew(i,j)=newmesh%Te(i,j)*CeOld(i,j)
        UhNew(i,j)=newmesh%Th(i,j)*ChOld(i,j)
        
        intensity(i,j)=M_ZERO
        MaxHeating(i,j)=M_ZERO
        MaxHeatingTime(i,j)=M_ZERO
        epsilonNeedle(i,j)=matter%EpsStatic-M_ONE
        
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
   
   potential(:,:)=M_ZERO! (0d0,0d0)
   
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
   DielectricStatic(:,:)=DielectricStatic(:,:)+M_ONE !so that equals 1 outside needle, and equals epsilonStatic0 inside
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


  allocate(CurviWx(1:Params%M,1:Params%N))
  allocate(CurviWy(1:Params%M,1:Params%N))
  allocate(CurviEx(1:Params%M,1:Params%N))
  allocate(CurviEy(1:Params%M,1:Params%N))
  allocate(CurviNx(1:Params%M,1:Params%N))
  allocate(CurviNy(1:Params%M,1:Params%N))
  allocate(CurviSx(1:Params%M,1:Params%N))
  allocate(CurviSy(1:Params%M,1:Params%N))

  allocate(TangentWx(1:Params%M,1:Params%N))
  allocate(TangentWy(1:Params%M,1:Params%N))
  allocate(TangentEx(1:Params%M,1:Params%N))
  allocate(TangentEy(1:Params%M,1:Params%N))
  allocate(TangentNx(1:Params%M,1:Params%N))
  allocate(TangentNy(1:Params%M,1:Params%N))
  allocate(TangentSx(1:Params%M,1:Params%N))
  allocate(TangentSy(1:Params%M,1:Params%N))

   call allocate_NormCurviTangent(Params%M, Params%N, NormalN, NormalS, NormalE, &
                                 NormalW )


   call compute_distances(Params%M, Params%N, x, y, DistN, DistS, DistE, DistW, DistDualN, &
                          DistDualS, DistDualE, DistDualW, CellAreaN, CellAreaS, CellAreaE, CEllAreaW )
   !
   call compute_norm_tan_curv(Params%M, Params%N, x, y, NormalN, NormalS, NormalE, NormalW, &
                              TangentNx, TangentNy, TangentSx, TangentSy, &
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


!   call poisson_init_dual( Params%M,Params%N, x, y, xDualSW, yDualSW, xDualSE, yDualSE, &
!                            xDualNE, yDualNE, xDualNW, yDualNW, xDual, yDual)

   ! Drift initialization
   JeX(:,:)=M_ZERO
   JeY(:,:)=M_ZERO
   JhX(:,:)=M_ZERO
   JhY(:,:)=M_ZERO
  
   
   !! defining material index and ionization constants
   epsilonInf=DielectricConstant(source%lambda)
   OnePhotonIonizationRate0=OnePhotonIonizationRate(matter)
   TwoPhotonIonizationRate0=TwoPhotonIonizationRate(source%lambda)
   
  write(*,*) 'epsilon(', 1d9*source%lambda, 'nm)=', epsilonInf
  write(*,*) 'Re(sqrt(epsilon))=', real(sqrt(epsilonInf))
if(Params%UseMieScattering.eq.1) then
  write(*,*) 'Computing the Mie scattering field distribution...'
  write(*,*) 'Angle Mie =', Params%phiMie0
  write(*,*) 'Polarization TM ? ', Params%PolarizationSource

  !$OMP PARALLEL DEFAULT(NONE) SHARED(x, y, Params, &
  !$OMP phiMie, Radius)
  !$OMP DO COLLAPSE(2)
  do j=1,Params%N
    do i=1,Params%M
        if(y(i,j)<0d0) then
          phiMie(i,j)=Params%phiMie0+M_PI
        else
          phiMie(i,j)=Params%phiMie0
        end if
!           find the radius for the cylindrical Mie scattering model
!           Radius(i,j)=y(i,j)
        if(ExpNeedleType.eq.0) then
          Radius(i,j)=ConeExp1Radius(1d6*y(i,j), 1d6*x(i,j), M_ZERO)
        else
          Radius(i,j)=ConeExp2Radius(1d6*y(i,j), 1d6*x(i,j), M_ZERO)
        end if
      end do
    end do
  !$OMP END DO
  !$OMP END PARALLEL
    

  !$OMP PARALLEL DEFAULT(NONE) SHARED(x, y, source, Params, epsilonInf, &
  !$OMP phiMie, Radius, EintField,EintField2)
  !$OMP DO COLLAPSE(2)
   do j=1,Params%N
    do i=1,Params%M
          if(Params%PolarizationSource.eq.1) then !TM polarization, Bassel et al scattering on a cylinder
          ! formula for an experimental needle with interpolated radius
!             write(*,*) "TM polarization selected."
            EintField(i,j)= MieScattering(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf, source%k) ! * sqrt(2d0*source%fluence/(c*epsilon0*source%tau))
            EintField2(i,j)=M_ZERO
          ! formula with a super mistake on radius
!           EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), 0.5d0*(y(i,N)-y(i,1)), epsilonInf) ! * sqrt(2d0*source%fluence/(c*epsilon0*source%tau))

          ! formulas for an hyperbolic needle
!           EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), abs(ContourYofX(x(i,j), NeedleRadius, NeedleAngle)), epsilonInf) ! * sqrt(2d0*source%fluence/(c*epsilon0*source%tau))

          ! formula for debug, using a constant radius
!         EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), 100d-9, epsilonInf) ! * sqrt(2d0*source%fluence/(c*epsilon0*source%tau)) !with a constant radius

          else !TE polarization
!             write(*,*) "TE polarization selected."
            EintField2(i,j)= MieScatteringTE2(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf, source%k)
            EintField(i,j) = MieScatteringTE1(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf, source%k)
          end if
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

!     EintFieldR=sqrt(EintField * conjg(EintField))
    EintFieldR=real(sqrt( EintField * conjg(EintField) + EintField2 * conjg(EintField2) ))
    

    call output_open(Field%unit,'output/Field.dat', .false.) ! format 891
    do j=1,Params%N
      do i=1,Params%M
        write(Field%unit, 891, advance='yes') x(i,j), y(i,j), EintFieldR(i,j), Radius(i,j)
891        FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
      end do
    end do
    close(Field%unit)


  write(*,*) 'Done.'
 end if
 
  t=tmin
  
  write(*,*) 'Elapsed time : ', timer_elapsedtime( full_timer )
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
  maxIntensity=M_ZERO; maxTe=M_ZERO; minTe=1d10; maxTh=M_ZERO; minTh=1d10; maxTs=M_ZERO; minTs=1d10;
  maxNe=M_ZERO; maxNh=M_ZERO; minNe=1d50; minNh=1d50
  maxCFLxT=M_ZERO; maxCFLyT=M_ZERO; maxCFLxN=0d0; maxCFLyN=M_ZERO; maxCFLxTs=M_ZERO; maxCFLyTs=M_ZERO;
  maxFermiIndexE=0; maxFermiIndexH=0; 
 
  !NTD: Why DistX are recomputed here? Same for CellAreaX
  !TJYD: To treat boundary conditions and treat everything with a loop on the complete mesh. This should be kept.
  do i=1,Params%M
    CellAreaN(i,Params%N)=M_ZERO
!     CellAreaN(i,N-1)=M_ZERO
!     CellAreaS(i,2)=M_ZERO
    CellAreaS(i,1)=M_ZERO
    
    DistN(i,Params%N-1)=sqrt((M_HALF*(x(i,Params%N-1)+x(i,Params%N))-x(i,Params%N-1))**2 &
                            +(M_HALF*(y(i,Params%N-1)+y(i,Params%N))-y(i,Params%N-1))**2)
    DistS(i,2)=sqrt((x(i,2)-M_HALF*(x(i,2)+x(i,1)))**2+(y(i,2)-M_HALF*(y(i,2)+y(i,1)))**2)
  end do
   
  do j=1,Params%N
    CellAreaW(1,j)=M_ZERO
!     CellAreaW(2,j)=M_ZERO
!     CellAreaE(M-1,j)=M_ZERO
    CellAreaE(Params%M,j)=M_ZERO
    
    DistW(2,j)=sqrt((x(2,j)-M_HALF*(x(2,j)+x(1,j)))**2+((y(2,j)-M_HALF*(y(2,j)+y(1,j))))**2)
    DistE(Params%M-1,j)=sqrt((M_HALF*(x(Params%M-1,j)+x(Params%M,j))-x(Params%M-1,j))**2&
                            +(M_HALF*(y(Params%M-1,j)+y(Params%M,j))-y(Params%M-1,j))**2)
  end do
  
     
    do j=2, Params%N-1
      do i=2, Params%M-1
        ! Calcul de Grad(Ne) sur le maillage direct
        ! Première estimation peu stable
        GradNeX(i,j) = M_HALF * InvCellVol(i,j) * &
                      ( (mesh%Ne(i,j) + mesh%Ne(i,j+1)) * CellAreaN(i,j) * NormalN%x(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i,j-1)) * CellAreaS(i,j) * NormalS%x(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i-1,j)) * CellAreaW(i,j) * NormalW%x(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i+1,j)) * CellAreaE(i,j) * NormalE%x(i,j) )
        GradNeY(i,j) = M_HALF * InvCellVol(i,j) * &
                      ( (mesh%Ne(i,j) + mesh%Ne(i,j+1)) * CellAreaN(i,j) * NormalN%y(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i,j-1)) * CellAreaS(i,j) * NormalS%y(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i-1,j)) * CellAreaW(i,j) * NormalW%y(i,j) &
                      * (mesh%Ne(i,j) + mesh%Ne(i+1,j)) * CellAreaE(i,j) * NormalE%y(i,j) )
                      
        ! interpolation lineaire des valeurs de phi sur les bords de cellules
!         phiN=Ne(i,j)+GradNeX(i,j)*(M_HALF*(xDual(i,j)+xDual(i-1,j)))+GradNeY(i,j)*(M_HALF*(yDual(i,j)+yDual(i-1,j)))
!         phiS=Ne(i,j)+GradNeX(i,j)*(M_HALF*(xDual(i,j-1)+xDual(i-1,j-1)))+GradNeY(i,j)*(M_HALF*(yDual(i,j-1)+yDual(i-1,j-1)))
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
  call timer_start(full_timer)

  call Profiler_stop(prof_init)

  if(Params%RestartCalc == 1) then
    call Restart_load( mesh, UeNew, UhNew, TsOld, Ce, Ch, CsPrev, CsOld, Cs, t, nmin, &
          IntensityEnergy, ElectronEnergy, ElectronKineticEnergy, ElectronPotentialEnergy, HoleEnergy, LatticeEnergy )
  else
    nmin = 0
  end if

  !***************************************************************
  !***************** temporal loop *******************************
  !***************************************************************
  do nbiter=nmin+1, nmax

   call Profiler_start(prof_timeloop, 'TIME LOOP')
    
   t=t+dt;
    
   NeTotal=M_HALF
   NhTotal=M_HALF

   ! replacing old datas
   !$OMP PARALLEL DEFAULT(NONE) SHARED(newmesh, Ue, Uh, UeNew, UhNew, &
   !$OMP TsPrev, TsOld, mesh, CeOld, Ce, ChOld, Ch, CsPrev2, CsPrev, CsOld, Cs )
   !$OMP DO COLLAPSE(2)
   do j=1, newmesh%N
     do i=1, newmesh%M
       Ue(i,j)     = UeNew(i,j)
       Uh(i,j)     = UhNew(i,j)
       TsPrev(i,j) = TsOld(i,j)
       TsOld(i,j)  = mesh%Ts(i,j)
       CeOld(i,j)  = Ce(i,j)
       ChOld(i,j)  = Ch(i,j)

       CsPrev2(i,j)= CsPrev(i,j)
       CsPrev(i,j) = CsOld(i,j)
       CsOld(i,j)  = Cs(i,j)
     end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL


   !$OMP PARALLEL DEFAULT(NONE) SHARED(mesh, newmesh)
   call copy_mesh(mesh, newmesh)
   !$OMP END PARALLEL

   !TODO: Does this depends on the position? If yes, this has to be changed bak to an array. 
   !TODO: TJYD: this can depend on position, if we apply a model for the collision frequency. This would be nice, actually. 
   nuColl= get_collision_frequency(matter)
   !
   !
   call DielectricFunction_batch(mesh, Dielectric, OpticalIndex, OpticalDamping, Reflectivity, &
                                 epsilonInf, nuColl, me, source)
   !
   !$OMP PARALLEL DEFAULT(NONE) SHARED (Params, mesh, DielectricDrudeE, DielectricDrudeH, absorptionDrudeE, &
   !$OMP absorptionDrudeH, nuColl, me, mh, source)
   call ComputeDielectricFunctionDrude_batch(Params, mesh, mesh%Ne, DielectricDrudeE, absorptionDrudeE, nuColl, me, source)
   !
   call ComputeDielectricFunctionDrude_batch(Params, mesh, mesh%Nh, DielectricDrudeH, absorptionDrudeH, nuColl, mh, source)
   !$OMP END PARALLEL
   !
   !
   call DensitiesOfState_batch(mesh, DOSe, DOSh, meDOS, mhDOS)
   !
   !$OMP PARALLEL DEFAULT(NONE) SHARED (Params, mesh, source, intensity, OpticalIndex, Reflectivity, &
   !$OMP absorptionDrudeE, absorptionDrudeH, OnePhotonIonizationRate0, TwoPhotonIonizationRate0, &
   !$OMP t, t0, sigmaTau, I0, sigmaX, sigmaY, x, y, x0, y0,  &
   !$OMP sigmaX1, sigmaY1, sigmaX2, sigmaY2, sigmaX3, sigmaY3, sigmaX4, sigmaY4, sigmaX5, sigmaY5, &
   !$OMP sigmaX6, sigmaY6, sigmaX7, sigmaY7, sigmaX8, sigmaY8, sigmaX9, sigmaY9, x1, y1, x2, y2, &
   !$OMP x3, y3, x4, y4, x5, EintFieldR,  &
   !$OMP y5, x6, y6, x7, y7, x8, y8, x9, y9, I1, I2, I3, I4, I5, I6, I7, I8, I9)
   !This routine computes the intensity for the entire grid with one call
   call ComputeIntensity_batch(Params, mesh, source, intensity, OpticalIndex, Reflectivity, &
                               absorptionDrudeE, absorptionDrudeH, OnePhotonIonizationRate0, TwoPhotonIonizationRate0, &
                               t, t0, sigmaTau, I0, sigmaX, sigmaY, x, y, x0, y0, DefectThickness, BandBendingInFDTD,  &
                               sigmaX1, sigmaY1, sigmaX2, sigmaY2, sigmaX3, sigmaY3, sigmaX4, sigmaY4, sigmaX5, sigmaY5, &
                               sigmaX6, sigmaY6, sigmaX7, sigmaY7, sigmaX8, sigmaY8, sigmaX9, sigmaY9, x1, y1, x2, y2, &
                               x3, y3, x4, y4, x5, EintFieldR,  &
                               y5, x6, y6, x7, y7, x8, y8, x9, y9, I1, I2, I3, I4, I5, I6, I7, I8, I9 )
   !$OMP END PARALLEL
   !
   !
!!!! thermal calculations in the main domain
! calculation of sources
   !$OMP PARALLEL DEFAULT(NONE) SHARED (Params, mesh, DOSe, DOSh, FermiRatioH, FermiRatioE, &
   !$OMP FermiIndexE, FermiIndexH )
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
   !$OMP END PARALLEL
   !
   !
   !
   !Computes the electron and mobilities for the entire mesh
   call ComputeMobilities_batch(mesh, mobilityE, mobilityH, &
                                FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                ColFermi0, ColFermiHalf, nuColl, me)
   !
   !
   !Computes the diffusion terms for the entire mesh
   call UpdateDiffusions_batch(Params, mesh, diffusionE, diffusionH, mobilityE, mobilityH, &
                                FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                ColFermiHalf, ColFermiMenusHalf)
   !
   !
   !Computes the drif vectors for the entire mesh
   call UpdateDriftVectors_batch(mesh, JeX, JeY, JhX, JhY, mobilityE, mobilityH, Ex, Ey, DriftOn)
   !
   !
   ! Computes the electron, hole and lattice heat capacities for the entire mesh
   call ComputeHeatCapacities_batch(mesh, Ce, Ch, Cs, invCe, invCh, invCs, &
                                    FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                    ColFermiThreeHalf, ColFermiHalf, ColFermiMenusHalf, ColFermiEta)
   !
   !Updates the couplings for the entire mesh
   call UpdateCouplings_batch(Params, mesh, CouplingE, CouplingH, Ce, Ch)
   !
   !
   !Computes the Sources Gains and Losses terms for electron and holes
   call ComputeGainsAndLosses(Params, mesh, source, matter, Egap, intensity, OnePhotonIonizationRate0, &
                              TwoPhotonIonizationRate0, absorptionDrudeE, absorptionDrudeH, &
                              Ce, Ch, CeOld, ChOld, dt, me, mh, &
                              GainsE, GainsH, SourceUe, SourceUh, SourceE, SourceH, LossesE, LossesH, Params%ImpactOff )
   !
   !Compute the new conductivites, based on the knowledge of densities and mobilities
   call UpdateConductivities_batch(mesh, kappae, kappah, kappas, mobilityE, mobilityH, &
                                           FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                                           ColFermi0, ColFermi1, ColFermi2, Params%TransportModel)
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
   else !TODO: This is redondant with copy_mesh operation at the begining of the temporal loop !TJYD: True...
     newmesh%Ne(:,:)=mesh%Ne(:,:)
   end if
   !
   !
   if(Params%HolesOff.eq.0 .AND. Params%NeOff.eq.0) then
     call computeNh( newmesh, mesh, dual, dt, InvCellVol, GainsH, LossesH, diffusionH, &
                     ShapeFactorNormalE, ShapeFactorTangentE, NormalE%N, &
                     ShapeFactorNormalW, ShapeFactorTangentW, NormalW%N, &
                     ShapeFactorNormalN, ShapeFactorTangentN, NormalN%N, &
                     ShapeFactorNormalS, ShapeFactorTangentS, NormalS%N )
   else !So for this one? 
     newmesh%Nh(:,:)=mesh%Nh(:,:)
   endif
   !
   !
   !$OMP PARALLEL DEFAULT(NONE) SHARED (dt, UeNew, UhNew, &
   !$OMP mesh, newmesh, dual, source, Params,  &
   !$OMP Ue, Uh, FermiTableE, FermiTableH, &
   !$OMP x, y, diffusionE, diffusionH, GainsE, GainsH, LossesE, LossesH, &
   !$OMP kappae, kappah, kappas, Ce, CeOld, Ch, ChOld, Cs, CsOld, CsPrev, CouplingE, CouplingH, &
   !$OMP SourceE, SourceH, SourceUe, SourceUh, &
   !$OMP VeX, VeY, VhX, VhY, &
   !$OMP NormalN, NormalS, NormalE, NormalW, &
   !$OMP CellVolume, CellAreaN, CellAreaS, CellAreaE, CellAreaW, CellVol, InvCellVol, &
   !$OMP DistN, DistS, DistE, DistW, DistDualN, DistDualS, DistDualE, DistDualW, &
   !$OMP ShapeFactorNormalE, ShapeFactorNormalN, ShapeFactorNormalS, ShapeFactorNormalW, &
   !$OMP ShapeFactorTangentE, ShapeFactorTangentN, ShapeFactorTangentS, ShapeFactorTangentW, &
   !$OMP h1, h2, h3, invCe, invCh, invCs)
   !
   if(Params%TeOff.ne.1) then
     !
     if(Params%ConvectionEnergy.eq.0) then
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
       if(Params%ConvectionEnergy.eq.0) then
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
   !$OMP END PARALLEL
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
   if(Params%ConvectionEnergy.eq.1) then
     call computeConvection( mesh, newmesh, UeNew, UhNew, Ue, Uh, invCe, invCh, &
                             FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                             ColFermiThreeHalf, ColFermiHalf, ColFermiMenusHalf, ColFermiEta )
   end if
   !
   !$OMP PARALLEL DEFAULT(NONE) SHARED (dt, x, y, Params, diffusionE, kappae, kappas, &
   !$OMP CFLxT, CFLyT, CFLxN, CFLyN, CFLxTs, CFLyTs, DistN, DistS, DistE, DistW, invCe, invCh, invCs) &
   !$OMP PRIVATE(work)
   !$OMP DO COLLAPSE(2)
   do j=2, Params%N-1
     do i=2, Params%M-1
       !
       work = dt/(M_HALF*(DistW(i,j)+DistE(i,j)))**2
       CFLxT(i,j)  = kappae(i,j)*invCe(i,j) * work
       CFLxTs(i,j) = kappas(i,j)*invCs(i,j) * work
       !
       work = dt/(M_HALF*(DistN(i,j)+DistS(i,j)))**2
       CFLyT(i,j)  = kappae(i,j)*invCe(i,j) * work
       CFLyTs(i,j) = kappas(i,j)*invCs(i,j) * work
       !
       !TODO: Optimise
       CFLxN(i,j)=diffusionE(i,j)*dt/(x(i,j)-x(i-1,j))**2 !+dt/(x(i,j)-x(i-1,j))*mobilityE(i,j)*sqrt(Ex(i,j)**2+Ey(i,j)**2)
       CFLyN(i,j)=diffusionE(i,j)*dt/(y(i,j)-y(i,j-1))**2 !+dt/(y(i,j)-y(i,j-1))*mobilityE(i,j)*sqrt(Ex(i,j)**2+Ey(i,j)**2)
       !
     end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of parallel section
   !
   call applyBoundaryConditions( newmesh, UeNew, UhNew, GradNeX, GradNeY, DriftOn )
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! CHECKING the results
   !
   maxCFLxN   = maxval(CFLxN)
   maxCFLyN   = maxval(CFLyN)
   maxCFLxT   = maxval(CFLxT)
   maxCFLyT   = maxval(CFLyT)
   maxCFLxTs  = maxval(CFLxTs)
   maxCFLyTs  = maxval(CFLyTs)
   !
   !Check and control quantity are computed only when needed
   if(mod(nbiter,Params%OutputIter).eq.0) then
     !
     maxIntensity = maxval(intensity)
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
     maxSourceE = maxval(SourceE)
     maxGainsE  = maxval(GainsE)
     maxSourceH = maxval(SourceH)
     maxGainsH  = maxval(GainsH)
     maxGap     = maxval(Egap)
     maxFermiIndexE = maxval(FermiIndexE)
     maxFermiIndexH = maxval(FermiIndexH)
     !
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
     write(*,*) "dt_init=", Params%TimeStep, "dt=", dt

     TotalThermalEnergy = 0.0d0
     TotalLaserEnergy = 0.0d0

     !$OMP PARALLEL DEFAULT(NONE) SHARED(Params, mesh, newmesh, x, y, TotalElectrons, TotalHoles, &
     !$OMP Ce, Ch, Cs, reflectivity, intensity, absorptionDrudeE, TotalThermalEnergy, &
     !$OMP absorptionDrudeH, OnePhotonIonizationRate0, TwoPhotonIonizationRate0, TotalLaserEnergy ) &
     !$OMP PRIVATE(work)
     !$OMP DO COLLAPSE(2) REDUCTION(+:TotalLaserEnergy, TotalThermalEnergy)
     do j=2, Params%N-1 !(optimized)
       do i=2, Params%M-1 !(optimized)
         !
         !TODO:This is only needed for a reduction, so lets do the reduction directly here
         TotalElectrons(i,j)= newmesh%Ne(i,j)*(0.125d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1)) &
                            -0.125d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1)))
         TotalHoles(i,j)    =newmesh%Nh(i,j)*(0.125d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1)) &
                            -0.125d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1)))
         !
         !
         TotalThermalEnergy= TotalThermalEnergy + Ce(i,j)*mesh%Te(i,j)+Ch(i,j)*mesh%Th(i,j)+Cs(i,j)*mesh%Ts(i,j)
         !
         !
         work = intensity(i,j)*(1d0-reflectivity(i,j))
         TotalLaserEnergy = TotalLaserEnergy + & !TODO: is this expression valid ?!
                         OnePhotonIonizationRate0 * work    & !energy loss by one photon absorption
                       + TwoPhotonIonizationRate0 * work**2 & !energy loss by two photon absorption
                  + (absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*work !energy loss by carrrier heating
       end do
     end do
     !$OMP END DO
     !$OMP END PARALLEL

     TotalMeshVolume=0d0
     !$OMP PARALLEL DEFAULT(NONE) SHARED(mesh, CellVol, NeTotal,NhTotal,TotalNumOfE,  &
     !$OMP TotalMeshVolume, TotalNumOfH, TotalElectrons, &
     !$OMP TotalHoles)
     !$OMP DO COLLAPSE(2) REDUCTION(+:NeTotal,NhTotal,TotalNumOfE, TotalNumOfH,  &
     !$OMP TotalMeshVolume)
     do j=1,mesh%N
       do i=1,mesh%M
         NeTotal=NeTotal + mesh%Ne(i,j) * CellVol(i,j)
         NhTotal=NhTotal + mesh%Nh(i,j) * CellVol(i,j)
         TotalNumOfE=TotalNumOfE + TotalElectrons(i,j)
         TotalNumOfH=TotalNumOfH + TotalHoles(i,j)
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

    if(t > 100d0*source%tau) then
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(mesh,MaxHeating, MaxHeatingTime, t, Params ) &
      !$OMP COLLAPSE(2)
      do j=1,Params%N
        do i=1,Params%M
          if(MaxHeating(i,j) < mesh%Ts(i,j)) then
            MaxHeating(i,j)=mesh%Ts(i,j)
            MaxHeatingTime(i,j)=t
          end if
        end do
      end do
      !$OMP END PARALLEL DO
    endif

    !TODO: Add a type Energy

    ! calculation of the absorbed laser energy involved in the simulated slice !
    IntensityEnergy=IntensityEnergy+(OnePhotonIonizationRate0+absorptionDrudeE(1,Params%N/2) &
                 +absorptionDrudeH(1,Params%N/2))*intensity(1,Params%N/2)*CellVol(1,Params%N/2)*dt


    call evaluate_bandgap(matter, newmesh, newmesh%Ne,newmesh%Ts, Egap_new)

    !$OMP PARALLEL DEFAULT(NONE) SHARED(Params, mesh, newmesh, matter, Egap, Egap_new, &
    !$OMP Ce, CeOld, Ch, ChOld, Cs, CsOld, OpticalIndex, CellVol, &
    !$OMP ElectronEnergy, ElectronKineticEnergy, ElectronPotentialEnergy,HoleEnergy,LatticeEnergy) &
    !$OMP PRIVATE(work)
    !$OMP DO COLLAPSE(2) REDUCTION(+:ElectronEnergy, ElectronKineticEnergy, ElectronPotentialEnergy,HoleEnergy,LatticeEnergy)
    do j=1,Params%N
      do i=1,Params%M

!       if(TransportModel.eq.-1) then
!          TeNew(i,j)=Tout; ThNew(i,j)=Tout;
!       end if
!       CAUTION: These definitions are erroneously including initial temperature into account.
!       ElectronEnergy=ElectronEnergy+Ce(i,j)*Te(i,j)*CellVol(i,j)
!       HoleEnergy=HoleEnergy+Ch(i,j)*Th(i,j)*CellVol(i,j)
!       LatticeEnergy=LatticeEnergy+Cs(i,j)*Ts(i,j)*CellVol(i,j)

        work = Egap_new(i,j) - Egap(i,j)
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
!       ElectronEnergy=ElectronKineticEnergy+ElectronPotentialEnergy !already summed over time

        HoleEnergy=HoleEnergy+(Ch(i,j)*(newmesh%Th(i,j)-mesh%Th(i,j))+(Ch(i,j)-ChOld(i,j))*mesh%Th(i,j)) * CellVol(i,j) !kinetic energy
        
        LatticeEnergy=LatticeEnergy+((Cs(i,j)*(newmesh%Ts(i,j)-mesh%Ts(i,j)))+0d0*(Cs(i,j)-CsOld(i,j))*mesh%Ts(i,j))*CellVol(i,j) !dCs/dt=0, 20150426, TJYD.

     end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL

   call check_divergences(mesh, maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN, x, y, t, nbiter, nmin, Params )

   !TODO: Move this to check_divergences
   !$OMP PARALLEL DEFAULT(NONE) SHARED(Params, FermiIndexE, FermiIndexH, &
   !$OMP FermiRatioE, FermiRatioH, mesh, DOSe, DOSh, t)
   !$OMP DO COLLAPSE(2)
   do j=1,Params%N
     do i=1,Params%M

       if(real(FermiIndexE(i,j)) > real(FermiMaxLines) .OR. real(FermiIndexE(i,j)) < M_ONE) then
         write(*,*) "t,i,j,FermiIndexE(i,j)=", t,i,j,FermiIndexE(i,j)
       end if
       if(real(FermiIndexH(i,j)) > real(FermiMaxLines) .OR. real(FermiIndexH(i,j)) < M_ONE) then
         write(*,*) "t,i,j,FermiIndexH(i,j)=", t,i,j,FermiIndexH(i,j)
       end if
       if(real(FermiRatioE(i,j)) < M_ZERO .OR. real(FermiRatioH(i,j)) < M_ZERO) then
         write(*,*) "Problem in DOS or Ne. DOS(i,j)=", i,j,DOSe(i,j), DOSh(i,j), "Ne,h(i,j)=", mesh%Ne(i,j), mesh%Nh(i,j)
       end if
     end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL
  end if

  !TODO: There should be a module managing the adaptative time step
  ! lets change dt when fast reponse is finished in order to catch the long one.

  ! saving the timesteps of several previous steps (used for the high order calculation of d/dt).
  dt4=dt3;
  dt3=dt2;
  dt2=dt;
  ! chaning the timestep based on known behavior of the system
  if(Params%AdaptativeTimeStep.eq.1) then
    if((t>1d1*source%tau*coeffDilaDt) .AND. (dt.eq.Params%TimeStep) .AND. &
      (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
      dt=10d0*Params%TimeStep
    else if ((t > M_HALF*source%tau*coeffDilaDt) .AND. (dt.eq.10d0*Params%TimeStep) .AND. (maxCFLxN+maxCFLyN &
      + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
      dt=40d0*Params%TimeStep
    else if ((t > 1d3*source%tau*coeffDilaDt) .AND. (dt.eq.40d0*Params%TimeStep) .AND. (maxCFLxN+maxCFLyN + maxCFLxT &
      + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
        dt=1d3*Params%TimeStep
!           else if (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs > maxCFL) then
!             dt=dt/1d1
    end if
  end if

   !
   if(mod(nbiter,Params%OutputIter*iterOutMaps).eq.0  .or. nbiter == 1) then
     !
     !
     call output_open(Depth%unit,'output/Depth.dat', (nbiter/=1 .or. Params%RestartCalc == 1))               ! format 887
     do i=1,Params%M
       do j=1,Params%N
         !
         write(Depth%unit,'(46(1E12.5,3x))', advance="yes") t, x(i,j), y(i,j), intensity(i,j), mesh%Te(i,j), & !5
                        mesh%Th(i,j), mesh%Ts(i,j), mesh%Ne(i,j), mesh%Nh(i,j), reflectivity(i,j), & !10
                        absorptionDrudeE(i,j), absorptionDrudeH(i,j), TotalElectrons(i,j), & !13
                        TotalHoles(i,j), real(FermiIndexE(i,j)), REAL(FermiIndexH(i,j)), FermiRatioE(i,j), FermiRatioH(i,j), & !18
                        SourceE(i,j), SourceH(i,j), GainsE(i,j), GainsH(i,j), LossesE(i,j), & !23
                        LossesH(i,j), real(DielectricDrudeE(i,j)), aimag(DielectricDrudeE(i,j)), Egap(i,j), real(Dielectric(i,j)), & !28
                        aimag(Dielectric(i,j)), MaxHeatingTime(i,j), MaxHeating(i,j), real(potentialNeedle(i,j)), Ex(i,j), & !33
                        Ey(i,j), diffusionE(i,j), diffusionH(i,j), GradNeX(i,j), GradNeY(i,j), &!38
                        real(EintField(i,j)), aimag(EintField(i,j)), EintFieldR(i,j), EintFieldI(i,j), phiMie(i,j), & !43
                        Radius(i,j)
         !
       end do !on Y
       write(Depth%unit,'(3x)', advance="yes")
     end do !on X
     close(Depth%unit);
     !
     ! writing result for the electrostatic calculations
     call output_open(DepthVessel%unit,'output/DepthVessel.dat', (nbiter/=1 .or. Params%RestartCalc == 1))       ! format 889
     do i=1,Mp
       do j=1,Np
         ! ecriture des donnees dans un fichier different
         write(DepthVessel%unit,889, advance="yes") t, xP(i,j), yP(i,j), real(potential(i,j)), real(ExPoisson(i,j)), & !
                real(EyPoisson(i,j)), DielectricStatic(i,j), NeP(i,j), NhP(i,j)
  889 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, & !TODO: Please use short notation with prenthesis !!
  1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5) 
        end do
      end do
      close(DepthVessel%unit);
      !
      ! write the functions on Dual Mesh
      !
      call output_open(DualDepth%unit,'output/DualDepth.dat', (nbiter/=1 .or. Params%RestartCalc == 1))! format 890
      do i=1,Params%M-1
        do j=1,Params%N-1
          !
          write(DualDepth%unit, '(9(1E12.5, 3x))', advance="YES") t, xDual(i,j), yDual(i,j), dual%Te(i,j), dual%Th(i,j), & !5
                                          dual%Ts(i,j), dual%Ne(i,j), dual%Nh(i,j), intensityDual(i,j) !9
          !
        end do
      end do
      close(DualDepth%unit)
      !
    end if
    !
    !
    ! output to files
    if(mod(nbiter,Params%OutputIter).eq.0 .or. nbiter == 1) then
      !
      cpu_timestep_duration = timer_elapsedtime(full_timer) / real(nbiter)
      cpuefficiency=real(nthreads)/cpu_timestep_duration
      !
      !This should be moved to output.F90 file
      !
      call output_open(EnergyConservation%unit, 'output/EnergyConservation.dat', (nbiter/=nmin .or. Params%RestartCalc == 1)) !format 892
      write(EnergyConservation%unit,892, advance="YES") t, IntensityEnergy, ElectronEnergy, HoleEnergy, LatticeEnergy, & !5
          TotalMeshVolume, LaserIntensityEnergy, ElectronKineticEnergy, ElectronPotentialEnergy !9
892 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
      close(EnergyConservation%unit)

      call output_open(TimeMax%unit,'output/TimeMax.dat', (nbiter/=1 .or. Params%RestartCalc == 1))                ! format 888
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
       close(TimeMax%unit)

       call output_open(TimeApex%unit, 'output/TimeApex.dat', (nbiter/=1 .or. Params%RestartCalc == 1))         ! format 884
       write(TimeApex%unit,884, advance="YES") t, mesh%Te(1,Params%N/2), mesh%Th(1,Params%N/2), &
              mesh%Ts(1,Params%N/2), mesh%Ne(1,Params%N/2), &                        !5
              mesh%Nh(1,Params%N/2), intensity(1,Params%N/2), TotalLaserEnergy, TotalThermalEnergy, &        !9
              SourceE(1,Params%N/2), GainsE(1,Params%N/2), SourceH(1,Params%N/2), GainsH(1,Params%N/2), Egap(1,Params%N/2), &                !14
              real(FermiIndexE(1,Params%N/2)),&
               real(FermiIndexH(1,Params%N/2)), Ce(2,Params%N/2), &                !17
              CeOld(2,Params%N/2), Ch(2,Params%N/2), ChOld(2,Params%N/2), Cs(2,Params%N/2), CsOld(2,Params%N/2)                               !22
              
884 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
       close(TimeApex%unit)

       call output_open(TimeUp%unit, 'output/TimeUp.dat', (nbiter/=1 .or. Params%RestartCalc == 1))          ! format 883
       write(TimeUp%unit,883, advance="YES") t, mesh%Te(Params%M/2,Params%N), mesh%Th(Params%M/2,Params%N),&
                          mesh%Ts(Params%M/2,Params%N), mesh%Ne(Params%M/2,Params%N), &
              mesh%Nh(Params%M/2,Params%N), intensity(Params%M/2,Params%N), TotalLaserEnergy, TotalThermalEnergy, &
              SourceE(Params%M/2,Params%N), GainsE(Params%M/2,Params%N), SourceH(Params%M/2,Params%N),&
               GainsH(Params%M/2,Params%N), Egap(Params%M/2,Params%N), &
              real(FermiIndexE(Params%M/2,Params%N)), real(FermiIndexH(Params%M/2,Params%N))
              
883 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
      close(TimeUp%unit)


      call output_open(TimeBottom%unit, 'output/TimeBottom.dat', (nbiter/=1 .or. Params%RestartCalc == 1))         ! format 882
      write(TimeBottom%unit,882, advance="YES") t, mesh%Te(Params%M/2,1), mesh%Th(Params%M/2,1), mesh%Ne(Params%M/2,1), &
              mesh%Nh(Params%M/2,1), intensity(Params%M/2,1), TotalLaserEnergy, TotalThermalEnergy, &
              SourceE(Params%M/2,1), GainsE(Params%M/2,1), SourceH(Params%M/2,1), GainsH(Params%M/2,1), Egap(Params%M/2,1), &
              real(FermiIndexE(Params%M/2,1)), real(FermiIndexH(Params%M/2,1))
              
882 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, & !TODO: Please use short notation with prenthesis !!
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
      close(TimeBottom%unit)

    end if
    !
    if(mod(nbiter,Params%DumpInterval).eq.0) then
      call Restart_dump(mesh, UeNew, UhNew, TsOld, Ce, Ch, CsPrev, CsOld, Cs, t, nbiter, &
            IntensityEnergy, ElectronEnergy, ElectronKineticEnergy, ElectronPotentialEnergy, HoleEnergy, LatticeEnergy )
    end if
    !
  call Profiler_stop(prof_timeloop)
    !
  end do !end of time loop

  
  call releasemesh(mesh)
  call releasemesh(dual)
  call releasemesh(newmesh)

  call deallocate_Norm(NormalN, NormalS, NormalE, NormalW )

  call CloseOutputs()

  call ReleaseInputParameters( Params )

  call Profiler_write_report( )

end program Flaps

