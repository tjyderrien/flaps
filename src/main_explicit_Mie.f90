!------------------------------------------------------------------------------
!> @file main_explicit_Mie.f90
!
! DESCRIPTION:
!> @brief **** calculate temperature distribution in 2D in a tip ****
!
!> @author
!> Thibault J.Y. Derrien
!  Laboratoire Hubert Curien, UMR CNRS, St-Etienne
!  ANR Ultrasonde
!
!> @date
!> Jul-Dec 2012 - Initial Version
!------------------------------------------------------------------------------

program Flaps

! include 'Bivariate.f'
! USE Bivariate
USE libmsh2vf !Script provided by A. Mouton, Univ Lille1, France for GMSH interfacing
!   use control_file !Script provided by Jason Blevins, Ohio State University

use Types

implicit none

    type(MeshValues) :: mesh, dual, newmesh

    real(8), parameter:: lambda=515d-9         , & !laser wavelength (m)
                        fluence=0d0                , & !laser fluence (J.m-2)
                        tau=40d-15                , & !FWHM pulse duration (s)
                        spotX=50d-6                , & !FWHM spot size in X direction (1030nm: 400nm x 50nm ; 515nm: 50um x 50 um ; 343 nm: 50um x 100nm)
                        spotY=50d-6                , & !FWHM spot size in Y direction
                        xCenter=1000d-9      ,& ! X position of the max of the intensity (1030nm: 1um 0um, 515nm: idem, 343nm: 100nm x 200nm)
                        yCenter=0d0*200d-9      ,&! Y position of the max of the intensity
                        Tout=80d0 ,&  !external temperature (K)
                        potential0=7d3,&         ! potential at the bottom of the needle ; default = 7d3
                        potentialNull=0d0 !, &
 !                       phiMie0=1d0*acos(-1d0)                ! Mie scattering: plane angle in cylindrical coordinates
    
    real(8), parameter:: dt0=1d-18,& !time step (s)
                        tmax=-190d-15,& !stop time
                        coeffDilaDt=2d0        ,& !diltation coeff before dt change
                        xmin=-10d-6       ,& !mesh min
                        xmax=10d-6       ,& !mesh max
                        ymin=-10d-6       ,&                
                        ymax=10d-6       ,&
                        tCenter=0d0       ,&         !time of gaussian intensity maximum
                        tmin=tCenter-5d0*tau                     !max absolute time
    
                        
    integer(8), parameter::  iterOut=1       ,& ! number of iterations between each stdout
                        iterOutMaps=1000      ,& ! number of outputs for maps between each stdout
                          M=11   ,& !number of cells main domain X direction
                          N=11     ,& !number of cells main domain Y direection
                          VirtualPoints=3, & !number of virtual points to exclude from the GMSH file (locate them at the beginning!)
                          Mv=101       ,& !number of celles in the Vessel domain (larger) X direction
                          Nv=101        ,& !number of celles in the Vessel domain (larger) Y direction
                          MeshChoice=1       ,& !0: rectangle (xmin,xmax)(ymin,ymax). 1: Experimental cones, 2: Cone in a vessel (HS), 3: import GMSH (working)
                          MeshIterations=500000        ,&        !number of iterations to calculate meshNeedle
                          MeshIterationsVessel=100*Mv,&        !number of iterations to calculate meshVessel
                          MeshShift=1       ,&         !number of cells x N in the tip, 343 nm: 2; 515 nm: 3;
                          FermiMaxLines=3584        ,&        ! >= number of lines in Fermi file
                          SORiterations=1        ,&        !iteration number for over-relaxation method
                          InterpolateMethod=1        ,&        ! 0: linear, 1: bicubic
                          UseInterpolation=0        ,&
                          AdaptativeTimeStep=1, &
                          SolveImplicit=1, &        ! 0: use explicit schemes, 1: use implicit scheme (band diagonal matrixes)
                          numberOfNeighbours=4
                          
                          
        
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
                            NeOff=0       ,&
                            TeOff=0       ,&                   !0: Disable temperature calculations
                            HolesOff=0       ,&
                            TsOff=0       ,&
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
                            UseMieScattering=-1,&                 ! 1: Enable Mie scattering analytic formula, 0: badly fitted FDTD input, -1: constant intensity
                            maxBesselOrder=20,&                ! Max of terms in series of Bessel for Mie scattering
                            besselArray=1, &
                            NewtonIterations=1000, &
                            ExpNeedleType=0
        
    real(8), parameter:: pi=acos(-1d0),&         !pi number
                          hbar=1.05457d-34  ,&         !planck constant
                          epsilon0=8.85418781762d-12 ,& !vacuum dielectric permittivity
                          mu0=4d0*pi*1d-7       ,& ! vacuum magnetic permeability
                          ec=1.60217646d-19        ,&         !elementary charge
                          me0=9.10938188d-31        ,&         !electron mass
                          c=2.99792458d8  ,&         !light speed
                          kb=1.3806488d-23,&          !Boltzmann constant, 
                          SiDensity=2.329d3        ,&           !Silicon rest density
                          epsilonStatic0=11.66570433d0 !,0.01404457712d0)                ! dielectric constant for static field

    real(8), parameter:: omegaLaser=2d0*pi*c/lambda        ,& !laser pulsation (s**-1)
                         me=0.5d0*me0       ,& ! electron effective mass for conductivity !0.24 (source ?)
                         mh=0.5d0*me0       ,&   ! hole effective mass for conductivity !0.81 (source ?)
                         meDOS=0.36d0*me0        ,& ! electron effective mass for DOS
                         mhDOS=0.81d0*me0       ,& ! hole effective mass for DOS
                         Ne0=1d5                                , &
                         Nh0=1d5, &                                !initial density (to calculate auto using fermi!) : at 80 K, Ne=1d5
                         Nlimit=1d0       ,&! lowest possible density
                         Nborder=1d23       ,&! density on boundaries to consider defect layer
                         DefectThickness=1d-7                        
    
    complex(8), parameter:: Imaginary=(0d0,1d0), Unit=(1d0,0d0), Zero=(0d0,0d0)                ! complex unity
!                      epsilonStatic0=(11.66570433d0,0.01404457712d0)                ! dielectric constant for static field
    
    integer(8)         nbiter, i, itwo, jtwo, j, k, l, nmax, NeedleIndexX, NeedleIndexY, maxFermiIndexE, maxFermiIndexH, &
                Mp, Np, imax, NewtonIteration, RunningIndex
    logical        Diverged
    real(8)         t, t0, dx, dy, x0, y0, dt, dt1, dt2, dt3, dt4, h1, h2, h3, h4
    real(8)         Te0, Th0, I0 !initial values of the problem
    real(8)        Ue(1:M, 1:N), & !electron energy
                Uh(1:M, 1:N), & !hole energy
                UeNew(1:M, 1:N), & !electron energy
                UhNew(1:M, 1:N), & !hole energy
                TsOld(1:M, 1:N), & !lattice temperature (time n-1)
                TsPrev(1:M,1:N), & !lattice temperature (time n-2)
                GradNeX(1:M, 1:N),& !Grad(Ne)_x
                GradNeY(1:M, 1:N),& !Grad(Ne)_y
                intensity(1:M, 1:N), & !propagated intensity
                intensity2(1:M-1, 1:N-1), & !test
                intensityDual(1:M-1,1:N-1), & !intensity dual, just for test of the function
                reflectivity(1:M,1:N), & ! surface reflectivity
                absorptionDrudeE(1:M, 1:N), & ! absorption coefficient
                absorptionDrudeH(1:M, 1:N), & ! absorption coefficient
                diffusionE(1:M, 1:N), & !fick diffusion coefficient for electrons
                diffusionH(1:M,1:N), & !fick diffusion coefficient for holes
                GainsE(1:M, 1:N), GainsH(1:M, 1:N), &
                LossesE(1:M, 1:N), LossesH(1:M, 1:N), &
                kappae(1:M, 1:N), kappah(1:M, 1:N), kappas(1:M,1:N), &
                Ce(1:M, 1:N), Ch(1:M, 1:N), Cs(1:M, 1:N), &
                CeOld(1:M, 1:N), ChOld(1:M,1:N), CsOld(1:M, 1:N), &
                CsPrev(1:M, 1:N), CsPrev2(1:M, 1:N), &
                CouplingE(1:M, 1:N), CouplingH(1:M, 1:N), &
                nuColl(1:M, 1:N) , &!        total collision frequency
                nuColleph(1:M, 1:N) , &!        electron-phonon collision frequency
                mobilityE(1:M, 1:N), mobilityH(1:M, 1:N), & 
                etae(1:M,1:N), etah(1:M,1:N), &        ! reduced chemical Fermi potential
                Egap(1:M, 1:N), &                        ! local gap value
                SourceE(1:M, 1:N), SourceH(1:M, 1:N), & ! heating sources
                SourceUe(1:M, 1:N), SourceUh(1:M, 1:N), & ! free carrier thermal energy sources
                diffNe(1:M, 1:N), diffNh(1:M, 1:N), &         ! just for derivation in time
                x(1:M, 1:N), y(1:M, 1:N), &                 ! needle position indexes
                xNew(1:M, 1:N), yNew(1:M, 1:N), &                 ! to make the mesh converge
                xV(1:Mv, 1:Nv), yV(1:Mv, 1:Nv), &                 ! vessel position indexes
                xDualSW(1:M, 1:N), yDualSW(1:M, 1:N), &                 ! dual mesh position 
                xDualSE(1:M, 1:N), yDualSE(1:M, 1:N), &                 
                xDualNE(1:M, 1:N), yDualNE(1:M, 1:N), &                 
                xDualNW(1:M, 1:N), yDualNW(1:M, 1:N), &
                xDual(1:M-1,1:N-1), yDual(1:M-1, 1:N-1), &
                CFLxT(1:M, 1:N), CFLyT(1:M, 1:N), CFLxN(1:M, 1:N), CFLyN(1:M, 1:N), CFLxTs(1:M,1:N), CFLyTs(1:M,1:N), &
                ThermalEnergy(1:M, 1:N), LaserEnergy(1:M, 1:N), &
                TotalElectrons(1:M, 1:N), TotalHoles(1:M, 1:N), &
                DOSe(1:M, 1:N), DOSh(1:M, 1:N), &
                FermiRatioE(1:M,1:N), FermiRatioH(1:M,1:N), &
                OmegaX(1:M, 1:N), OmegaY(1:M, 1:N), & ! drift vectors for energy
                JeX(1:M, 1:N), JeY(1:M, 1:N), & ! drift vectors for particles
                JhX(1:M, 1:N), JhY(1:M, 1:N), & ! drift vectors for particles
                VeX(1:M, 1:N), VeY(1:M, 1:N), &
                VhX(1:M, 1:N), VhY(1:M, 1:N), &
                Ex(1:M,1:N), Ey(1:M,1:N), &        ! fields in the main domain
                potentialNeedle(1:M,1:N), &        ! potential in the needle
                DummyNeedle(1:M, 1:N), &                ! optional arguments for interpolation
                DummyDual(1:M-1, 1:N-1), &                ! optional arguments for interpolation
                epsilonNeedle(1:M,1:N), &        ! dielectric static in the needle
                MaxHeating(1:M, 1:N), &
                MaxHeatingTime(1:M, 1:N), &
                MeshDensity(1:Mv, 1:Nv)        ,&                ! function to adapt Poisson mesh on the needle mesh
                CellAreaN(1:M, 1:N), CellAreaS(1:M, 1:N), &                        ! area of the finite elements 
                CellAreaE(1:M, 1:N), CellAreaW(1:M, 1:N), &
                CellVol(1:M, 1:N), InvCellVol(1:M, 1:N),&          !volume of needle mesh cells
                NormalWx(1:M,1:N), NormalWy(1:M,1:N), NormalW2(1:M,1:N), &
                NormalEx(1:M,1:N), NormalEy(1:M,1:N), NormalE2(1:M,1:N),  &                ! normal to quadrangle elements
                NormalNx(1:M,1:N), NormalNy(1:M,1:N), NormalN2(1:M,1:N),  &
                NormalSx(1:M,1:N), NormalSy(1:M,1:N), NormalS2(1:M,1:N),  &
                TangentWx(1:M,1:N), TangentWy(1:M,1:N), &                ! Tangent to quadrangle elements
                TangentEx(1:M,1:N), TangentEy(1:M,1:N), &                
                TangentNx(1:M,1:N), TangentNy(1:M,1:N), &
                TangentSx(1:M,1:N), TangentSy(1:M,1:N), &
                CurviWx(1:M,1:N), CurviWy(1:M,1:N), &                 ! Unit vector between cell centers
                CurviEx(1:M,1:N), CurviEy(1:M,1:N), &                
                CurviNx(1:M,1:N), CurviNy(1:M,1:N), &
                CurviSx(1:M,1:N), CurviSy(1:M,1:N), &
                DistW(1:M,1:N), DistE(1:M,1:N), &                 ! distance to the center of neighboor cells
                DistN(1:M,1:N), DistS(1:M,1:N), &
                DistDualW(1:M,1:N), DistDualE(1:M,1:N), &                 ! distance of the element side (equal to area in 2D)
                DistDualN(1:M,1:N), DistDualS(1:M,1:N), &
                ShapeFactorNormalE(1:M, 1:N), ShapeFactorNormalW(1:M, 1:N), &
                ShapeFactorNormalS(1:M, 1:N), ShapeFactorNormalN(1:M, 1:N), &
                ShapeFactorTangentE(1:M, 1:N), ShapeFactorTangentW(1:M, 1:N), &
                ShapeFactorTangentS(1:M, 1:N), ShapeFactorTangentN(1:M, 1:N), &
                EintFieldR(1:M,1:N), EintFieldI(1:M, 1:N), &
                EintFieldDual(1:M-1, N-1), &
                phiMie(1:M, 1:N), &
                Radius(1:M, 1:N)
                
    integer(8)  FermiIndexE(1:M,1:N), FermiIndexH(1:M,1:N), &
                MeshVertice(1:M, 1:N), & ! data from the GMSH file
                SomeNeighbours(1:4,1:2)   !Neighbours for the fixed potential
    
    real(8) phiMie0
    real(8) Int2
    integer(8) PolarizationSource             ! Value of the Mie angle that will be distributed on various processors
    
    real(8), allocatable, target :: FermiTableE(:,:),&
                                     FermiTableH(:,:),& !reduced Fermi level for electrons and holes
                                     Amatrix(:,:), Bvector(:),         &
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
    CHARACTER(LEN=60)               :: namefile_msh, namefile_vf
!
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: vertices
    INTEGER, DIMENSION(:,:), POINTER          :: points, segments, triangles, quadrangles, boundedges, edges
    INTEGER, DIMENSION(:), POINTER            :: dim_physical_entities, id_physical_entities, idvertices
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: name_physical_entities
!     INTEGER                                   :: kmouton
    CHARACTER(LEN=1)                          :: choice
    INTEGER                           :: nb_vertices, nb_triangles, nb_quadrangles, nb_edges, nb_boundedges
    
    complex(8)         Dielectric(1:M,1:N), &! solid dielectric function under laser illumination
                DielectricDrudeE(1:M,1:N), & ! Drude part of dielectric function under laser illumination
                DielectricDrudeH(1:M,1:N), &
                EintField(1:M,1:N), EintField2(1:M,1:N)        !Ez internal field for Mie scattering theory

                
    complex(8) epsilonInf !, SORsum !material constant
   
    real(8) AugerRateE, AugerRateH, &
            sigmaTau, sigmaX, sigmaY, &
            maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN, maxCFLxTs, maxCFLyTs, &
            maxTe, minTe, maxTh, minTh, maxTs, minTs, maxIntensity, maxNe, minNe, maxNh, minNh, &
            maxEnergy, maxSourceE, maxGainsE, maxSourceH, maxGainsH, maxGap, maxDiffNe, maxDiffNh, &
            TotalLaserEnergy, TotalThermalEnergy, ElectronPotentialEnergy, ElectronKineticEnergy, &
            cpuefficiency, cpu_timestep_duration, calc_time_begin, calc_time_1, calc_time_2, calc_time_3, &
            NeedleHeight, NeedleA, NeedleB, Needlet0Limit, NeedleXParam, NeedleYParam, NeedleAngle, &
            localT, P2critic, ConstBLx, ConstBLy, TotalNumOfE, TotalNumOfH, &
            distX, distY, SORsum, xmin2, xmax2, ymin2, ymax2, &
            ErrorSum, MeshConvergence, MeshConvergenceOld
            
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
     integer(8) ConeExp1Radius, ConeExp2Radius !, Interpolate
     real(8) ConeExp1, ConeExp2, &
                 Tangent, Normal, AreaElement, AreaTri
            
    character(len=50)::format
!OPENMP declarations
    integer :: myid, nthreads
    integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

! call omp_set_num_threads(16)
    
   !  !!********OpenMP Test*************
  myid=1 !if openMP is off, then test is disabled; else: will be set to 0 by OpenMP
  nthreads=1 !idem

  !$OMP PARALLEL default(none) private(myid) &
  !$OMP shared(nthreads)
  ! Determine the number of threads and their id
        myid = OMP_GET_THREAD_NUM()
        PRINT *, 'Hello from thread =', myid
        nthreads = OMP_GET_NUM_THREADS()
  !$OMP BARRIER
  
  if (myid==0) then 
    write(*,'(a)') 'OpenMP TEST'
    write(*,'(a,i1)') 'Number of Threads = ', nthreads
    write(*,'(a)') '*******************************'
  end if 
  !$OMP END PARALLEL
! !!******* END OpenMP test


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
RunningIndex=RunningIndex+1; MeshVertice(M,1) = RunningIndex
RunningIndex=RunningIndex+1; MeshVertice(M,N) = RunningIndex
RunningIndex=RunningIndex+1; MeshVertice(1,N) = RunningIndex

!North
do i=2, M-1
  RunningIndex=RunningIndex+1; MeshVertice(i, 1)=RunningIndex
end do
!East
do j=2, N-1
  RunningIndex=RunningIndex+1; MeshVertice(M, j)=RunningIndex
end do
!South
do i=M-1,2, -1
  RunningIndex=RunningIndex+1; MeshVertice(i,N)=RunningIndex
end do
!West
do j=N-1, 2, -1
  RunningIndex=RunningIndex+1; MeshVertice(1,j)=RunningIndex
end do

!rest of the domain
do i=2,M-1
  do j=2, N-1
    RunningIndex=RunningIndex+1; MeshVertice(i,j)=RunningIndex
  end do
end do


WRITE(*,*) 'Latest running index while remeshing', RunningIndex


!***** Compute geometrical data


!**** INITIALIZATION

  Diverged=.false.

  dt=dt0
  dt2=dt0
  dt3=dt0
  dt4=dt0

  ColFermiNeNc=2; ColFermiEta=3; ColFermi0=4; ColFermi1=5; ColFermi2=6; ColFermiHalf=7; 
  ColFermiThreeHalf=8; ColFermiMenusHalf=9;

! test field

  if(AugerOff.eq.0) then
    AugerRateE=2.3d-43
    AugerRateH=7.8d-44
  else
    AugerRateE=0d0
    AugerRateH=0d0
  end if

    sigmaTau=tau/(2d0*sqrt(2d0*log(2e0)))
    sigmaX=spotX/(2d0*sqrt(2d0*log(2e0)))
    sigmaY=spotY/(2d0*sqrt(2d0*log(2e0)))

    I0=fluence/tau * sqrt(4d0 * log(2d0) / pi)
    
    if(lambda.eq.515d-9) then
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

    sigmaX1=spotX1/(2e0*sqrt(2e0*log(2e0)))
    sigmaY1=spotY1/(2e0*sqrt(2e0*log(2e0)))
    sigmaX2=spotX2/(2e0*sqrt(2e0*log(2e0)))
    sigmaY2=spotY2/(2e0*sqrt(2e0*log(2e0)))
    sigmaX3=spotX3/(2e0*sqrt(2e0*log(2e0)))
    sigmaY3=spotY3/(2e0*sqrt(2e0*log(2e0)))
    sigmaX4=spotX4/(2e0*sqrt(2e0*log(2e0)))
    sigmaY4=spotY4/(2e0*sqrt(2e0*log(2e0)))
    sigmaX5=spotX5/(2e0*sqrt(2e0*log(2e0)))
    sigmaY5=spotY5/(2e0*sqrt(2e0*log(2e0)))
    sigmaX6=spotX6/(2e0*sqrt(2e0*log(2e0)))
    sigmaY6=spotY6/(2e0*sqrt(2e0*log(2e0)))
    sigmaX7=spotX7/(2e0*sqrt(2e0*log(2e0)))
    sigmaY7=spotY7/(2e0*sqrt(2e0*log(2e0)))
    sigmaX8=spotX8/(2e0*sqrt(2e0*log(2e0)))
    sigmaY8=spotY8/(2e0*sqrt(2e0*log(2e0)))
    sigmaX9=spotX9/(2e0*sqrt(2e0*log(2e0)))
    sigmaY9=spotY9/(2e0*sqrt(2e0*log(2e0)))
    
    t0=tCenter; x0=xCenter; y0=yCenter;
    
    
  ! opening files
    open(90,FILE='LaplaceConvergence.dat', access='sequential',status='unknown')
    open(91,FILE='LaplaceMatrix.dat',access='sequential',status='unknown')
    open(92,FILE='TimeBottom.dat', access='sequential', status='unknown')         ! format 882
    open(93,FILE='TimeUp.dat', access='sequential', status='unknown')         ! format 883
    open(94,FILE='TimeApex.dat', access='sequential', status='unknown')         ! format 884
    open(95,FILE='error.dat', access='sequential', status='unknown')
    open(96,FILE='parameters.dat', access='sequential', status='unknown')
    open(97,FILE='Depth.dat',access='sequential',status='unknown')                ! format 887
    open(98,FILE='TimeMax.dat',access='sequential',status='unknown')                 ! format 888
    open(99,FILE='mesh.dat',access='sequential',status='unknown')                ! format 885, 8852
    open(100,FILE='meshVessel.dat',access='sequential',status='unknown')
    open(101,FILE='DepthVessel.dat',access='sequential',status='unknown')         ! format 889
    open(102,FILE='meshElements.dat', access='sequential',status='unknown') ! format 881
    open(103,FILE='DualDepth.dat', access='sequential', status='unknown') ! format 890
    open(104,FILE='Field.dat', access='sequential', status='unknown') ! format 891
    open(105,FILE='EnergyConservation.dat', access='sequential', status='unknown') !format 892
    ! FORMAT numbers already used for writing: 887, 886, 888, 885, 882, 883
    
  
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
  call initmesh(mesh, M, N)
  call initmesh(dual, M-1, N-1)
  call initmesh(newmesh, M, N)

  ! building rectangular mesh 
  dx=(xmax-xmin)/(M+1)
  dy=(ymax-ymin)/(N+1)
  do i=1,M
    do j=1,N
      x(i,j)=0d0
      y(i,j)=0d0
    end do
  end do
  
  write(*,*)

  if(MeshChoice.eq.0) then  !rectangular mesh as main domain
    do i=1,M
      do j=1,N
        x(i,j)=dx*real(i)+xmin
        y(i,j)=dy*real(j)+ymin
        write(99,885, advance='yes') x(i,j), y(i,j), i, j ! writing as main mesh
      end do
      write(99,*) " "
    end do
  end if
!   else ! conical mesh as main domain
    
  if(MeshChoice==1 .OR. MeshChoice==2) then !conical mesh as main domain
    !boudary definition
    NeedleIndexX=M
    NeedleIndexY=N
    NeedleAngle=NeedleAngleDeg*pi/180d0
    NeedleA=NeedleRadius/(tan(NeedleAngle/2d0)**2)
    NeedleB=NeedleRadius/tan(NeedleAngle/2d0)
    Needlet0Limit=acos(NeedleRadius/(NeedleLength*tan(real(NeedleAngle)/2d0)**2+NeedleRadius))
    
!building of the conical mesh: write boundaries, then solve laplace, and iterate
    
    NeedleXParam=(Needlet0Limit-0d0)/NeedleIndexX !dt0 for X (0,t0)
    NeedleYParam=(Needlet0Limit-0d0)/NeedleIndexY !dt0 for Y (0,t0)
    
    P2critic=(real(N)-1d0)**2/(2d0*real(j)-real(N)-1)
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
        meshStepDt=(meshParameterTmax-meshParameterTmin)/(real(2*(M-1)+N)) !step of parameter t to define corners of the mesh

        !!! BOTTOM
        localTmin=meshParameterTmin !tmin in the segment (e.g. contour section)
        localTmax=-real(MeshShift)*meshStepDt !tmax in the line (e.g. contour section)
        
        if(localTmax<localTmin) then
          write(*,*) "MeshShift is too large. Reduce.", -localTmax/meshStepDt
          stop
        end if
        
        localdT=(localTmax-localTmin)/real(M-1) !parameter t to distribute the nodes on the segment
        do i=1,M
          localT=localTmax-real(i-1)*localdT
          if(ExpNeedleType.eq.1) then
            x(i,1)=ConeExp2(localT)
          else
            x(i,1)=ConeExp1(localT)
          end if
          y(i,1)=localT
        end do

        !!! TIP APEX
        localTmin=-real(MeshShift)*meshStepDt !tmin in the segment (e.g. contour section)
        localTmax=real(MeshShift)*meshStepDt !tmax in the line (e.g. contour section)
        
        if(localTmax<localTmin) then
          write(*,*) "MeshShift is nul ?! MeshShift=", MeshShift
          stop
        end if
        
        localdT=(localTmax-localTmin)/real(N-1) !parameter t to distribute the nodes on the segment
        do j=2,N
          localT=localTmin+real(j-1)*localdT
          if(ExpNeedleType.eq.1) then
            x(1,j)=ConeExp2(localT)
          else
            x(1,j)=ConeExp1(localT)
          end if
          y(1,j)=localT
        end do

        !!! UP
        localTmin=real(MeshShift)*meshStepDt !tmin in the segment (e.g. contour section)
        localTmax=meshParameterTmax !tmax in the line (e.g. contour section)
        
        if(localTmax<localTmin) then
          write(*,*) "MeshShift is too large. Reduce.", localTmax/meshStepDt
          stop
        end if
        
        localdT=(localTmax-localTmin)/real(M-1) !parameter t to distribute the nodes on the segment
        do i=1,M
          localT=localTmin+real(i-1)*localdT
          if(ExpNeedleType.eq.1) then
            x(i,N)=ConeExp2(localT)
          else
            x(i,N)=ConeExp1(localT)
          end if
          y(i,N)=localT
        end do

        !!! BACKSIDE
        localTmin=meshParameterTmin !tmin in the segment (e.g. contour section)
        localTmax=meshParameterTmax !tmax in the line (e.g. contour section)
        localdT=(localTmax-localTmin)/real(N-1) !parameter t to distribute the nodes on the segment
        do j=1,N
          localT=localTmin+real(j-1)*localdT
          if(ExpNeedleType.eq.1) then
            x(M,j)=ConeExp2(localTmax)
          else
            x(M,j)=ConeExp1(localTmax)
          end if
          y(M,j)=localT
        end do
!         
        
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
        
        ! resolution of laplace (sequential)
        do i=2, M-1
          do j=2,N-1
            x(i,j)=(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1))/4d0
            y(i,j)=(y(i+1,j)+y(i-1,j)+y(i,j+1)+y(i,j-1))/4d0
          end do
        end do
        
        
!         ! attempt of optimization using OpenMP
!         !$OMP PARALLEL DEFAULT(private) SHARED(x, y, xNew, yNew)
!         
!         ! save the contours (can be optimized)
!         
!         xNew(:,:)=x(:,:)
!         yNew(:,:)=y(:,:)
!         
!         !$OMP DO !(optimized)
!         do i=2, M-1
!         !$OMP PARALLEL DO !(optimized)
!           do j=2,N-1
!             xNew(i,j)=(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1))/4d0
!             yNew(i,j)=(y(i+1,j)+y(i-1,j)+y(i,j+1)+y(i,j-1))/4d0
!           end do        
!         !$OMP END PARALLEL DO
!         end do
!         !$OMP END DO
!         
!         ! replace with the new mesh
!         x(:,:)=xNew(:,:)
!         y(:,:)=yNew(:,:)
!         
!         !$OMP END PARALLEL

        
        
        MeshConvergenceOld=MeshConvergence
        MeshConvergence=0d0
        do i=1, M
          do j=1, N
            MeshConvergence=MeshConvergence+(x(i,j)**2+y(i,j)**2)
          end do
        end do
        
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
    do i=1,M
      do j=1,N
        write(99, 885, advance='yes') x(i,j), y(i,j), i, j
        885        FORMAT (1E15.8, 3x, 1E15.8, 3x, I4, 3X, I4)
        !write(99,*)
        end do
        write(99,*) " "
    end do

    
!     do i=1,M-1
!       do j=1, N-1
!         dx(i,j)=x(i+1,j)-x(i,j)
!         dy(i,j)=y(i,j+1)-y(i,j)
!       end do
!     end do
    ! now we can calculate height of the cone...
    NeedleHeight=2d0*tan(.5d0*NeedleAngle) * sqrt(1d0 - NeedleRadius**2/((tan(NeedleAngle/2d0))**4)* &
                (x(M,N/2)+NeedleRadius**2/((tan(0.5d0*NeedleAngle))**2)**2)) * (x(M,N/2)+NeedleRadius/(tan(NeedleAngle/2d0))**2)
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
    do i=1, M
      do j=1, N
    x(i,j)=vertices(idvertices(MeshVertice(i,j)),1)
    y(i,j)=vertices(idvertices(MeshVertice(i,j)),2)
      end do
    end do

    x(:,:)=1d-6*x(:,:)
    y(:,:)=1d-6*y(:,:)
    
    !writing of the mesh
    do i=1,M
      do j=1,N
        write(99, 8852, advance='yes') x(i,j), y(i,j), i, j, MeshVertice(i,j)
        8852        FORMAT (1E15.8, 3x, 1E15.8, 3x, I4, 3X, I4, 3X, I8)
        !write(99,*)
        end do
        write(99,*) " "
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
     Mp=M
     Np=N
     allocate(xP(1:Mp, 1:Np))
     allocate(yP(1:Mp, 1:Np))
     xP(:,:)=x(:,:)
     yP(:,:)=y(:,:)
  end if
  
  do i=1,Mp
    do j=1,Np
      write(100,885) xP(i,j), yP(i,j), i, j
    end do
    write(100,*) " "
  end do
      
  
  write (*,*) "Vessel cells:", Mp, "*", Np,"=", Mp*Np
  write (*,*) "Matter cells:", M, "*", N, "=", M*N
  
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
  allocate(FixedPotentialIndex(1:N, 1:2))
  allocate(CellVolume(1:Mp, 1:Np))
  allocate(NormalNxP(1:Mp,1:Np)); allocate(NormalNyP(1:Mp,1:Np))
  allocate(NormalSxP(1:Mp,1:Np)); allocate(NormalSyP(1:Mp,1:Np))
  allocate(NormalExP(1:Mp,1:Np)); allocate(NormalEyP(1:Mp,1:Np))
  allocate(NormalWxP(1:Mp,1:Np)); allocate(NormalWyP(1:Mp,1:Np))
  allocate(CellAreaNP(1:Mp,1:Np)); allocate(CellAreaSP(1:Mp,1:Np))
  allocate(CellAreaEP(1:Mp,1:Np)); allocate(CellAreaWP(1:Mp,1:Np))
  
  call flush(99); call flush(100)
  
  ! test of neighboors ! USEFUL FOR POISSON EQUATION ONLY
!  write(*,*) "[TEST] Testing neighborhood of a internal point..."
!  write(*,'(2I10.1)') transpose(FindNeighbours(x(2,2),y(2,2),xP,yP,numberOfNeighbours,Mp,Np))
!  write(*,*) "[TEST] Testing neighborhood of a external point..."
!  write(*,'(2I10.1)') transpose(FindNeighbours(xP(34,35),yP(34,35),x,y,numberOfNeighbours,M,N))

!*********** SOURCE IMPORT **************
  ! importing Fermi functions
   allocate(FermiTableE(1:9, 1:FermiMaxLines))
   allocate(FermiTableH(1:9, 1:FermiMaxLines))
   call TabCreateFL !(FermiTableE, FermiTableH)
    FermiTableE(:,:)=1d0; FermiTableH(:,:)=1d0; ! uncomment if you want to disable fermi-dirac. Dont forget to lock the FermiIndexes also.
!************ INITIALIZATION ************

  write(96,*) "========== CONE PARAMETERS ========="
  write(96,*) "Cone length=", NeedleLength*1d6, "um"
  write(96,*) "Cone height=", NeedleHeight*1d6, "um"
  write(96,*) "Cone angle=", NeedleAngleDeg, "deg"
  write(96,*) "Cone curvature radius=", NeedleRadius*1d9, "nm"
  write(96,*)
  write(96,*) "=========== LASER PARAMETERS ========"
  write(96,*) "Laser fluence=", fluence*1d-4, "J.cm-2"
  write(96,*) "Laser pulse duration=", tau*1d15, "fs"
  write(96,*) "Laser wavelength=", lambda*1d9, "nm"
  write(96,*) "Laser spot position: (X,Y)=", x0*1d6, y0*1d6, "um"
  write(96,*) "Laser spot size: (Sx, Sy)=", spotX*1d6, spotY*1d6, "um"
  write(96,*) "Mie scattering:", UseMieScattering
  write(96,*) "Laser polarization", PolarizationSource
  write(96,*)
  write(96,*) "============ MESH PARAMETERS =========="
  write(96,*) "Mesh size", M, "x", N
  write(96,*) "Dilatation time ratio=", coeffDilaDt
  write(96,*) "Mesh shift=", MeshShift
  write(96,*)
  write(96,*) "============ TIME CONTROL =========="
  write(96,*) "Initial timestep=", dt0
  write(96,*) "Maximal timestep=", tmin
  write(96,*) "Maximum time t=", tmax
  write(96,*) "Enable adaptative timestep=", AdaptativeTimeStep
  write(96,*) "Time output each ", iterOut, "iterations."
  write(96,*) "Map output each", iterOutMaps*iterOut, "iterations."
  
  call flush(96)
  
  ! initialisation
  call cpu_time(calc_time_begin) !initialisation time
  
  nmax=int((tmax-tmin)/dt, 8)
  
  Ex(:,:)=0d0 !-1d10
  Ey(:,:)=0d0 !-1d9 !0d0
  DummyVessel(:,:)=0d0
  
  Te0=Tout !1400d0
  Th0=Tout
  
  ! $ OMP DO
  do j=1,N
    do i=1,M
        
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
                  newmesh%Ne(i,j)=Ne0+Nborder*(exp(-0.5d0*(((x(i,j)-x(i,N))**2+(y(i,j)-y(i,N))**2) &
                    /((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
                  +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness) & 
                    /(2d0*sqrt(2d0*log(2d0))))**2)) &
                  +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness) & 
                    /(2d0*sqrt(2d0*log(2d0))))**2)) &
                  )
                  newmesh%Nh(i,j)=Nh0+Nborder*(exp(-0.5d0*(((x(i,j)-x(i,N))**2+(y(i,j)-y(i,N))**2) &
                    /((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
                  +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness) & 
                    /(2d0*sqrt(2d0*log(2d0))))**2)) &
                  +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness) & 
                    /(2d0*sqrt(2d0*log(2d0))))**2)) &
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
  do j=1,N
    do i=1,M
        
        ! initialise variables to calculate Ce, Ch
        DOSe(i,j)=DensityOfStateE(mesh%Te(i,j))
        DOSh(i,j)=DensityOfStateH(mesh%Th(i,j))
        FermiRatioE(i,j)=mesh%Ne(i,j)/DOSe(i,j)
        FermiRatioH(i,j)=mesh%Nh(i,j)/DOSh(i,j)
        FermiIndexE(i,j)=1! FermiIndex(FermiRatioE(i,j)) !1
        FermiIndexH(i,j)=1! FermiIndex(FermiRatioH(i,j)) !1
!         write(*,*) "iter=", nbiter, "DOS=", DOSe(i,j), DOSh(i,j)
        etae(i,j)=FermiTableE(ColFermiEta,FermiIndexE(i,j))
        etah(i,j)=FermiTableH(ColFermiEta,FermiIndexH(i,j))

        ! calculate semi-classical heat capacity
        CeOld(i,j)=1.5d0*mesh%Ne(i,j)*kb*(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j)) &
                              /FermiTableE(ColFermiHalf,FermiIndexE(i,j)) &
                  -etae(i,j)*(1d0-(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j)) & 
                              /FermiTableE(ColFermiHalf,FermiIndexE(i,j)))* &
                  (FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)) & 
                              /FermiTableE(ColFermiHalf,FermiIndexE(i,j)))))
        ChOld(i,j)=1.5d0*mesh%Nh(i,j)*kb*(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j)) &
                              /FermiTableH(ColFermiHalf,FermiIndexH(i,j)) &
                  -etah(i,j)*(1d0-(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j)) & 
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
   potentialNeedle(M,:)=potential0
   
   ! areas of the cell boudaries e.g. distances in this 2D code 
   ! and normal vectors (Nx, Ny) the four poles of quadrangle elements
   write(*,*) "[Mesh] Calculation of normales and distances."
   do j=2,N-1
      do i=2,M-1
      
        DistN(i,j)=sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
        DistS(i,j)=sqrt((x(i,j)-x(i,j-1))**2+(y(i,j)-y(i,j-1))**2)
        DistE(i,j)=sqrt((x(i+1,j)-x(i,j))**2+(y(i+1,j)-y(i,j))**2)
        DistW(i,j)=sqrt((x(i,j)-x(i-1,j))**2+(y(i,j)-y(i-1,j))**2)
        
        !TODO: Optimise by defining a stencil object

        NormalNx(i,j)=Normal(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                     0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)),1) !4d0*(0.25d0*y(i-1,j)+0.25d0*y(i-1,j+1)-0.25d0*y(i+1,j+1)-0.25d0*y(i+1,j))/(x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)-2d0*x(i-1,j)*x(i+1,j+1)-2d0*x(i-1,j)*x(i+1,j)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i+1,j+1)-2d0*x(i-1,j+1)*x(i+1,j)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)+x(i+1,j)**2+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)-2d0*y(i-1,j)*y(i+1,j+1)-2d0*y(i-1,j)*y(i+1,j)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i+1,j+1)-2d0*y(i-1,j+1)*y(i+1,j)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)+y(i+1,j)**2)**(0.5d0)
        NormalNy(i,j)=Normal(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                     0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)),2) !NormalNy(i,j)=-4d0*(0.25d0*x(i-1,j)+0.25d0*x(i-1,j+1)-0.25d0*x(i+1,j+1)-0.25d0*x(i+1,j))/(x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)-2d0*x(i-1,j)*x(i+1,j+1)-2d0*x(i-1,j)*x(i+1,j)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i+1,j+1)-2d0*x(i-1,j+1)*x(i+1,j)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)+x(i+1,j)**2+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)-2d0*y(i-1,j)*y(i+1,j+1)-2d0*y(i-1,j)*y(i+1,j)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i+1,j+1)-2d0*y(i-1,j+1)*y(i+1,j)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)+y(i+1,j)**2)**(0.5d0)
        NormalSx(i,j)=Normal(0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)),1) ! -4d0*(0.25d0*y(i-1,j-1)+0.25d0*y(i-1,j)-0.25d0*y(i+1,j)-0.25d0*y(i+1,j-1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)-2d0*x(i-1,j-1)*x(i+1,j)-2d0*x(i-1,j-1)*x(i+1,j-1)+x(i-1,j)**2-2d0*x(i-1,j)*x(i+1,j)-2d0*x(i-1,j)*x(i+1,j-1)+x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)+x(i+1,j-1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)-2d0*y(i-1,j-1)*y(i+1,j)-2d0*y(i-1,j-1)*y(i+1,j-1)+y(i-1,j)**2-2d0*y(i-1,j)*y(i+1,j)-2d0*y(i-1,j)*y(i+1,j-1)+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)+y(i+1,j-1)**2)**(0.5d0)
        NormalSy(i,j)=Normal(0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)),2) ! 4d0*(0.25d0*x(i-1,j-1)+0.25d0*x(i-1,j)-0.25d0*x(i+1,j)-0.25d0*x(i+1,j-1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)-2d0*x(i-1,j-1)*x(i+1,j)-2d0*x(i-1,j-1)*x(i+1,j-1)+x(i-1,j)**2-2d0*x(i-1,j)*x(i+1,j)-2d0*x(i-1,j)*x(i+1,j-1)+x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)+x(i+1,j-1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)-2d0*y(i-1,j-1)*y(i+1,j)-2d0*y(i-1,j-1)*y(i+1,j-1)+y(i-1,j)**2-2d0*y(i-1,j)*y(i+1,j)-2d0*y(i-1,j)*y(i+1,j-1)+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)+y(i+1,j-1)**2)**(0.5d0)
        NormalEx(i,j)=Normal(0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j+1)+x(i+1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i+1,j+1)+y(i+1,j)+y(i,j+1)+y(i,j)), 1) ! -4d0*(0.25d0*y(i+1,j-1)+0.25d0*y(i,j-1)-0.25d0*y(i+1,j+1)-0.25d0*y(i,j+1))/(x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)-2d0*x(i+1,j-1)*x(i+1,j+1)-2d0*x(i+1,j-1)*x(i,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i+1,j+1)-2d0*x(i,j-1)*x(i,j+1)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)+x(i,j+1)**2+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)-2d0*y(i+1,j-1)*y(i+1,j+1)-2d0*y(i+1,j-1)*y(i,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i+1,j+1)-2d0*y(i,j-1)*y(i,j+1)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)+y(i,j+1)**2)**(0.5d0)
        NormalEy(i,j)=Normal(0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)), &
                     0.25d0*(x(i+1,j+1)+x(i+1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i+1,j+1)+y(i+1,j)+y(i,j+1)+y(i,j)), 2) ! 4d0*(0.25d0*x(i+1,j-1)+0.25d0*x(i,j-1)-0.25d0*x(i+1,j+1)-0.25d0*x(i,j+1))/(x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)-2d0*x(i+1,j-1)*x(i+1,j+1)-2d0*x(i+1,j-1)*x(i,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i+1,j+1)-2d0*x(i,j-1)*x(i,j+1)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)+x(i,j+1)**2+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)-2d0*y(i+1,j-1)*y(i+1,j+1)-2d0*y(i+1,j-1)*y(i,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i+1,j+1)-2d0*y(i,j-1)*y(i,j+1)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)+y(i,j+1)**2)**(0.5d0)
        NormalWx(i,j)=Normal(0.25d0*(x(i-1,j+1)+x(i-1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i-1,j)+y(i,j+1)+y(i,j)), &
                     0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 1) ! 4d0*(0.25d0*y(i-1,j-1)+0.25d0*y(i,j-1)-0.25d0*y(i,j+1)-0.25d0*y(i-1,j+1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)-2d0*x(i-1,j-1)*x(i,j+1)-2d0*x(i-1,j-1)*x(i-1,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i,j+1)-2d0*x(i,j-1)*x(i-1,j+1)+x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)+x(i-1,j+1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)-2d0*y(i-1,j-1)*y(i,j+1)-2d0*y(i-1,j-1)*y(i-1,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i,j+1)-2d0*y(i,j-1)*y(i-1,j+1)+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)
        NormalWy(i,j)=Normal(0.25d0*(x(i-1,j+1)+x(i-1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i-1,j)+y(i,j+1)+y(i,j)), &
                     0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 2)! -4d0*(0.25d0*x(i-1,j-1)+0.25d0*x(i,j-1)-0.25d0*x(i,j+1)-0.25d0*x(i-1,j+1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)-2d0*x(i-1,j-1)*x(i,j+1)-2d0*x(i-1,j-1)*x(i-1,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i,j+1)-2d0*x(i,j-1)*x(i-1,j+1)+x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)+x(i-1,j+1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)-2d0*y(i-1,j-1)*y(i,j+1)-2d0*y(i-1,j-1)*y(i-1,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i,j+1)-2d0*y(i,j-1)*y(i-1,j+1)+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)
        
        
        TangentNx(i,j)=-Tangent(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)),0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                        0.25d0*(x(i,j)+x(i,j+1)+x(i-1,j+1)+x(i-1,j)),0.25d0*(y(i,j)+y(i,j+1)+y(i-1,j+1)+y(i-1,j)), 1)
        TangentNy(i,j)=-Tangent(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)),0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                        0.25d0*(x(i,j)+x(i,j+1)+x(i-1,j+1)+x(i-1,j)),0.25d0*(y(i,j)+y(i,j+1)+y(i-1,j+1)+y(i-1,j)), 2)
        TangentSx(i,j)=Tangent(0.25d0*(x(i,j)+x(i,j-1)+x(i+1,j-1)+x(i+1,j)),0.25d0*(y(i,j)+y(i,j-1)+y(i+1,j-1)+y(i+1,j)), &
                        0.25d0*(x(i-1,j)+x(i-1,j-1)+x(i,j-1)+x(i,j)),0.25d0*(y(i-1,j)+y(i-1,j-1)+y(i,j-1)+y(i,j)), 1)
        TangentSy(i,j)=Tangent(0.25d0*(x(i,j)+x(i,j-1)+x(i+1,j-1)+x(i+1,j)),0.25d0*(y(i,j)+y(i,j-1)+y(i+1,j-1)+y(i+1,j)), &
                        0.25d0*(x(i-1,j)+x(i-1,j-1)+x(i,j-1)+x(i,j)),0.25d0*(y(i-1,j)+y(i-1,j-1)+y(i,j-1)+y(i,j)), 2)
        TangentEx(i,j)=Tangent(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                        0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j-1)+x(i,j-1)),0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j-1)+y(i,j-1)), 1)
        TangentEy(i,j)=Tangent(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                        0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j-1)+x(i,j-1)),0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j-1)+y(i,j-1)), 2)
        TangentWx(i,j)=-Tangent(0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)),0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1)), &
                        0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)),0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j)), 1)
        TangentWy(i,j)=-Tangent(0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)),0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1)), &
                        0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)),0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j)), 2)
        
        CurviNx(i,j)=-Tangent(x(i,j),y(i,j),x(i,j+1),y(i,j+1),1)
        CurviNy(i,j)=-Tangent(x(i,j),y(i,j),x(i,j+1),y(i,j+1),2)
        CurviSx(i,j)=-Tangent(x(i,j),y(i,j),x(i,j-1),y(i,j-1),1)
        CurviSy(i,j)=-Tangent(x(i,j),y(i,j),x(i,j-1),y(i,j-1),2)
        CurviEx(i,j)=-Tangent(x(i,j),y(i,j),x(i+1,j),y(i+1,j),1)
        CurviEy(i,j)=-Tangent(x(i,j),y(i,j),x(i+1,j),y(i+1,j),2)
        CurviWx(i,j)=-Tangent(x(i,j),y(i,j),x(i-1,j),y(i-1,j),1)
        CurviWy(i,j)=-Tangent(x(i,j),y(i,j),x(i-1,j),y(i-1,j),2)
        
        CellAreaN(i,j)=sqrt(( 0.25d0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))-0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)))**2 & 
                  +(  0.25d0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))-0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)))**2) !0.25d0*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)-2d0*x(i+1,j+1)*x(i-1,j)-2d0*x(i+1,j+1)*x(i-1,j+1)+x(i+1,j)**2-2d0*x(i+1,j)*x(i-1,j)-2d0*x(i+1,j)*x(i-1,j+1)+x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)+x(i-1,j+1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)-2d0*y(i+1,j+1)*y(i-1,j)-2d0*y(i+1,j+1)*y(i-1,j+1)+y(i+1,j)**2-2d0*y(i+1,j)*y(i-1,j)-2d0*y(i+1,j)*y(i-1,j+1)+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)
        CellAreaS(i,j)=sqrt( (0.25d0*(x(i+1,j-1)+x(i,j-1)+x(i+1,j)+x(i,j))-0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i-1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i+1,j-1)+y(i,j-1)+y(i+1,j)+y(i,j))-0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i-1,j)+y(i,j)))**2 ) !0.25d0*(x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)-2d0*x(i+1,j)*x(i-1,j-1)-2d0*x(i+1,j)*x(i-1,j)+x(i+1,j-1)**2-2d0*x(i+1,j-1)*x(i-1,j-1)-2d0*x(i+1,j-1)*x(i-1,j)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)+x(i-1,j)**2+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)-2d0*y(i+1,j)*y(i-1,j-1)-2d0*y(i+1,j)*y(i-1,j)+y(i+1,j-1)**2-2d0*y(i+1,j-1)*y(i-1,j-1)-2d0*y(i+1,j-1)*y(i-1,j)+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)+y(i-1,j)**2)**(0.5d0)
        CellAreaE(i,j)=sqrt( (0.25d0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))-0.25d0*(x(i+1,j-1)+x(i,j-1)+x(i+1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))-0.25d0*(y(i+1,j-1)+y(i,j-1)+y(i+1,j)+y(i,j)))**2 ) !0.25d0*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)-2d0*x(i+1,j+1)*x(i+1,j-1)-2d0*x(i+1,j+1)*x(i,j-1)+x(i,j+1)**2-2d0*x(i,j+1)*x(i+1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)-2d0*y(i+1,j+1)*y(i+1,j-1)-2d0*y(i+1,j+1)*y(i,j-1)+y(i,j+1)**2-2d0*y(i,j+1)*y(i+1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)
        CellAreaW(i,j)=sqrt( (0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j))-0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i-1,j)+x(i,j)))**2 &
                   + (0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j))-0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i-1,j)+y(i,j)))**2 ) !0.25d0*(x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)-2d0*x(i,j+1)*x(i-1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i-1,j-1)-2d0*x(i-1,j+1)*x(i,j-1)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)-2d0*y(i,j+1)*y(i-1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i-1,j-1)-2d0*y(i-1,j+1)*y(i,j-1)+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)
        
        ! distance for cross diffusion, equal to cell border "areas" in 2D
        DistDualN(i,j)=CellAreaN(i,j)
        DistDualS(i,j)=CellAreaS(i,j)
        DistDualE(i,j)=CellAreaE(i,j)
        DistDualW(i,j)=CellAreaW(i,j)
        
      end do
   end do
   
!     if(DisableCrossDiffusion.eq.1) then
!       ! disable flux parallel to element boundaries
!       TangentNx(:,:)=0d0; TangentNy(:,:)=0d0; 
!       TangentSx(:,:)=0d0; TangentSy(:,:)=0d0
!       TangentEx(:,:)=0d0; TangentEy(:,:)=0d0
!       TangentWx(:,:)=0d0; TangentWy(:,:)=0d0
!       ! unnormalize flux normal to element boundaries
!       CurviNx(:,:)=NormalNx(:,:); CurviNy(:,:)=NormalNy(:,:)
!       CurviSx(:,:)=NormalSx(:,:); CurviSy(:,:)=NormalSy(:,:)
!       CurviEx(:,:)=NormalEx(:,:); CurviEy(:,:)=NormalEy(:,:)
!       CurviWx(:,:)=NormalWx(:,:); CurviWy(:,:)=NormalWy(:,:)
!     end if
   
   write(*,*) "NEEDLE CHECK"
   write(*,*) "North", NormalNx(M/2,N-1), NormalNy(M/2,N-1)
   write(*,*) "South", NormalSx(M/2,2), NormalSy(M/2,2)
   write(*,*) "East", NormalEx(M-1,N-1), NormalEy(M-1,N-1)
   write(*,*) "West", NormalWx(M-1,N-1), NormalWy(M-1,N-1)
   write(*,*) CellAreaN(M/2,N/2), CellAreaS(M/2,N/2), CellAreaE(M/2,N/2), CellAreaW(M/2,N/2)
   
   do j=2,Np-1
    do i=2,Mp-1
        NormalNxP(i,j)=4d0*(0.25d0*yP(i-1,j)+0.25d0*yP(i-1,j+1)-0.25d0*yP(i+1,j+1)-0.25d0*yP(i+1,j))/(xP(i-1,j)**2 & 
                      +2d0*xP(i-1,j)*xP(i-1,j+1)-2d0*xP(i-1,j)*xP(i+1,j+1)-2d0*xP(i-1,j)*xP(i+1,j)+xP(i-1,j+1)**2 & 
                      -2d0*xP(i-1,j+1)*xP(i+1,j+1)-2d0*xP(i-1,j+1)*xP(i+1,j)+xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j) & 
                      +xP(i+1,j)**2+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)-2d0*yP(i-1,j)*yP(i+1,j+1)-2d0*yP(i-1,j)*yP(i+1,j) & 
                      +yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i+1,j+1)-2d0*yP(i-1,j+1)*yP(i+1,j)+yP(i+1,j+1)**2 & 
                      +2d0*yP(i+1,j+1)*yP(i+1,j)+yP(i+1,j)**2)**(0.5d0)
        NormalNyP(i,j)=-4d0*(0.25d0*xP(i-1,j)+0.25d0*xP(i-1,j+1)-0.25d0*xP(i+1,j+1)-0.25d0*xP(i+1,j))/(xP(i-1,j)**2 & 
                +2d0*xP(i-1,j)*xP(i-1,j+1)-2d0*xP(i-1,j)*xP(i+1,j+1)-2d0*xP(i-1,j)*xP(i+1,j)+xP(i-1,j+1)**2 & 
                -2d0*xP(i-1,j+1)*xP(i+1,j+1)-2d0*xP(i-1,j+1)*xP(i+1,j)+xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j) & 
                +xP(i+1,j)**2+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)-2d0*yP(i-1,j)*yP(i+1,j+1)-2d0*yP(i-1,j)*yP(i+1,j) & 
                +yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i+1,j+1)-2d0*yP(i-1,j+1)*yP(i+1,j)+yP(i+1,j+1)**2 & 
                +2d0*yP(i+1,j+1)*yP(i+1,j)+yP(i+1,j)**2)**(0.5d0)
        NormalSxP(i,j)=-4d0*(0.25d0*yP(i-1,j-1)+0.25d0*yP(i-1,j)-0.25d0*yP(i+1,j)-0.25d0*yP(i+1,j-1))/(xP(i-1,j-1)**2 & 
                +2d0*xP(i-1,j-1)*xP(i-1,j)-2d0*xP(i-1,j-1)*xP(i+1,j)-2d0*xP(i-1,j-1)*xP(i+1,j-1)+xP(i-1,j)**2 & 
                -2d0*xP(i-1,j)*xP(i+1,j)-2d0*xP(i-1,j)*xP(i+1,j-1)+xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1) & 
                +xP(i+1,j-1)**2+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)-2d0*yP(i-1,j-1)*yP(i+1,j) & 
                -2d0*yP(i-1,j-1)*yP(i+1,j-1)+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i+1,j)-2d0*yP(i-1,j)*yP(i+1,j-1) & 
                +yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)+yP(i+1,j-1)**2)**(0.5d0)
        NormalSyP(i,j)=4d0*(0.25d0*xP(i-1,j-1)+0.25d0*xP(i-1,j)-0.25d0*xP(i+1,j)-0.25d0*xP(i+1,j-1))/(xP(i-1,j-1)**2 &
                      +2d0*xP(i-1,j-1)*xP(i-1,j)-2d0*xP(i-1,j-1)*xP(i+1,j)-2d0*xP(i-1,j-1)*xP(i+1,j-1)+xP(i-1,j)**2 & 
                      -2d0*xP(i-1,j)*xP(i+1,j)-2d0*xP(i-1,j)*xP(i+1,j-1)+xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1) & 
                      +xP(i+1,j-1)**2+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)-2d0*yP(i-1,j-1)*yP(i+1,j) & 
                      -2d0*yP(i-1,j-1)*yP(i+1,j-1)+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i+1,j)-2d0*yP(i-1,j)*yP(i+1,j-1) & 
                      +yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)+yP(i+1,j-1)**2)**(0.5d0)
        NormalExP(i,j)=-4d0*(0.25d0*yP(i+1,j-1)+0.25d0*yP(i,j-1)-0.25d0*yP(i+1,j+1)-0.25d0*yP(i,j+1))/(xP(i+1,j-1)**2 & 
                +2d0*xP(i+1,j-1)*xP(i,j-1)-2d0*xP(i+1,j-1)*xP(i+1,j+1)-2d0*xP(i+1,j-1)*xP(i,j+1)+xP(i,j-1)**2 & 
                -2d0*xP(i,j-1)*xP(i+1,j+1)-2d0*xP(i,j-1)*xP(i,j+1)+xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1) & 
                +xP(i,j+1)**2+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)-2d0*yP(i+1,j-1)*yP(i+1,j+1) & 
                -2d0*yP(i+1,j-1)*yP(i,j+1)+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i+1,j+1)-2d0*yP(i,j-1)*yP(i,j+1) & 
                +yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)+yP(i,j+1)**2)**(0.5d0)
        NormalEyP(i,j)=4d0*(0.25d0*xP(i+1,j-1)+0.25d0*xP(i,j-1)-0.25d0*xP(i+1,j+1)-0.25d0*xP(i,j+1))/(xP(i+1,j-1)**2 & 
                      +2d0*xP(i+1,j-1)*xP(i,j-1)-2d0*xP(i+1,j-1)*xP(i+1,j+1)-2d0*xP(i+1,j-1)*xP(i,j+1)+xP(i,j-1)**2 & 
                      -2d0*xP(i,j-1)*xP(i+1,j+1)-2d0*xP(i,j-1)*xP(i,j+1)+xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1) & 
                      +xP(i,j+1)**2+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)-2d0*yP(i+1,j-1)*yP(i+1,j+1) & 
                      -2d0*yP(i+1,j-1)*yP(i,j+1)+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i+1,j+1)-2d0*yP(i,j-1)*yP(i,j+1) & 
                      +yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)+yP(i,j+1)**2)**(0.5d0)
        NormalWxP(i,j)=4d0*(0.25d0*yP(i-1,j-1)+0.25d0*yP(i,j-1)-0.25d0*yP(i,j+1)-0.25d0*yP(i-1,j+1))/(xP(i-1,j-1)**2 & 
                      +2d0*xP(i-1,j-1)*xP(i,j-1)-2d0*xP(i-1,j-1)*xP(i,j+1)-2d0*xP(i-1,j-1)*xP(i-1,j+1)+xP(i,j-1)**2 & 
                      -2d0*xP(i,j-1)*xP(i,j+1)-2d0*xP(i,j-1)*xP(i-1,j+1)+xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1) & 
                      +xP(i-1,j+1)**2+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)-2d0*yP(i-1,j-1)*yP(i,j+1) & 
                      -2d0*yP(i-1,j-1)*yP(i-1,j+1)+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j+1)-2d0*yP(i,j-1)*yP(i-1,j+1) & 
                      +yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)
        NormalWyP(i,j)=-4d0*(0.25d0*xP(i-1,j-1)+0.25d0*xP(i,j-1)-0.25d0*xP(i,j+1)-0.25d0*xP(i-1,j+1))/(xP(i-1,j-1)**2 & 
                +2d0*xP(i-1,j-1)*xP(i,j-1)-2d0*xP(i-1,j-1)*xP(i,j+1)-2d0*xP(i-1,j-1)*xP(i-1,j+1)+xP(i,j-1)**2 & 
                -2d0*xP(i,j-1)*xP(i,j+1)-2d0*xP(i,j-1)*xP(i-1,j+1)+xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1) & 
                +xP(i-1,j+1)**2+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)-2d0*yP(i-1,j-1)*yP(i,j+1) & 
                -2d0*yP(i-1,j-1)*yP(i-1,j+1)+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j+1)-2d0*yP(i,j-1)*yP(i-1,j+1) & 
                +yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)
        
        CellAreaNP(i,j)=0.25d0*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j) & 
                -2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1) & 
                +xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j) & 
                -2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j) & 
                -2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)
        CellAreaSP(i,j)=0.25d0*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j) & 
                +xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)-2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2 & 
                +2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1) & 
                -2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1) & 
                -2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5d0)
        CellAreaEP(i,j)=0.25d0*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1) & 
                +xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1) & 
                +xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1) & 
                +yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1) & 
                +yP(i,j-1)**2)**(0.5d0)
        CellAreaWP(i,j)=0.25d0*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1) & 
                +xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1) & 
                +xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1) & 
                +yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2 & 
                +2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)
    end do
  end do
  
      do i=2,M-1
        
         
        ! NORTH
        
         CellVol(i,N)=0.25d0*AreaElement(x(i-1,N),y(i-1,N),x(i+1,N),y(i+1,N),x(i+1,N-1),y(i+1,N-1),x(i-1,N-1),y(i-1,N-1))
         CellAreaN(i,N)=0d0 !0.5d0*(x(i+1,N)**2-2d0*x(i+1,N)*x(i-1,N)+x(i-1,N)**2+y(i+1,N)**2-2d0*y(i+1,N)*y(i-1,N)+y(i-1,N)**2)**0.5d0
         CellAreaS(i,N)=0.25d0*(x(i+1,N)**2+2d0*x(i+1,N)*x(i+1,N-1)-2d0*x(i+1,N)*x(i-1,N-1)-2d0*x(i+1,N)*x(i-1,N) & 
                +x(i+1,N-1)**2-2d0*x(i+1,N-1)*x(i-1,N-1)-2d0*x(i+1,N-1)*x(i-1,N)+x(i-1,N-1)**2+2d0*x(i-1,N)*x(i-1,N-1) & 
                +x(i-1,N)**2+y(i+1,N)**2+2d0*y(i+1,N)*y(i+1,N-1)-2d0*y(i+1,N)*y(i-1,N-1)-2d0*y(i+1,N)*y(i-1,N) & 
                +y(i+1,N-1)**2-2d0*y(i+1,N-1)*y(i-1,N-1)-2d0*y(i+1,N-1)*y(i-1,N)+y(i-1,N-1)**2+2d0*y(i-1,N)*y(i-1,N-1)+y(i-1,N)**2)**0.5d0
         CellAreaW(i,N)=0.25d0*(x(i,N)**2+2d0*x(i,N)*x(i-1,N)-2d0*x(i,N)*x(i-1,N-1)-2d0*x(i,N)*x(i,N-1)+x(i-1,N)**2 & 
                -2d0*x(i-1,N)*x(i-1,N-1)-2d0*x(i-1,N)*x(i,N-1)+x(i-1,N-1)**2+2d0*x(i-1,N-1)*x(i,N-1)+x(i,N-1)**2 & 
                +y(i,N)**2+2d0*y(i,N)*y(i-1,N)-2d0*y(i,N)*y(i-1,N-1)-2d0*y(i,N)*y(i,N-1)+y(i-1,N)**2 & 
                -2d0*y(i-1,N)*y(i-1,N-1)-2d0*y(i-1,N)*y(i,N-1)+y(i-1,N-1)**2+2d0*y(i-1,N-1)*y(i,N-1) & 
                +y(i,N-1)**2)**0.5d0
         CellAreaE(i,N)=0.25d0*(x(i+1,N)**2+2d0*x(i+1,N)*x(i,N)-2d0*x(i+1,N)*x(i+1,N-1)-2d0*x(i+1,N)*x(i,N-1)+x(i,N)**2 & 
                -2d0*x(i,N)*x(i+1,N-1)-2d0*x(i,N)*x(i,N-1)+x(i+1,N-1)**2+2d0*x(i+1,N-1)*x(i,N-1)+x(i,N-1)**2 & 
                +y(i,N)**2+2d0*y(i,N)*y(i+1,N)-2d0*y(i,N)*y(i+1,N-1)-2d0*y(i,N)*y(i,N-1)+y(i+1,N)**2 & 
                -2d0*y(i+1,N)*y(i+1,N-1)-2d0*y(i+1,N)*y(i,N-1)+y(i+1,N-1)**2+2d0*y(i+1,N-1)*y(i,N-1) & 
                +y(i,N-1)**2)**0.5d0
         
         NormalEx(i,N)=-Normal(0.5d0*(x(i,N)+x(i+1,N)), 0.5d0*(y(i,N)+y(i+1,N)), 0.25d0*(x(i+1,N)+x(i,N)+x(i+1,N-1)+x(i,N-1)), & 
                0.25d0*(y(i+1,N)+y(i,N)+y(i+1,N-1)+y(i,N-1)), 1)
         NormalEy(i,N)=-Normal(0.5d0*(x(i,N)+x(i+1,N)), 0.5d0*(y(i,N)+y(i+1,N)), 0.25d0*(x(i+1,N)+x(i,N)+x(i+1,N-1)+x(i,N-1)), &
                0.25d0*(y(i+1,N)+y(i,N)+y(i+1,N-1)+y(i,N-1)), 2)
         NormalWx(i,N)=Normal(0.5d0*(x(i-1,N)+x(i,N)), 0.5d0*(y(i-1,N)+y(i,N)), 0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), & 
                0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),1); 
         NormalWy(i,N)=Normal(0.5d0*(x(i-1,N)+x(i,N)), 0.5d0*(y(i-1,N)+y(i,N)), 0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), & 
                0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),2); 
         NormalNx(i,N)=Normal(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),1)
         NormalNy(i,N)=Normal(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),2)
         NormalSx(i,N)=-Normal(0.25d0*(x(i+1,N-1)+x(i+1,N)+x(i,N)+x(i,N-1)), 0.25d0*(y(i+1,N-1)+y(i+1,N)+y(i,N)+y(i,N-1)), &
                0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)),0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),1) 
         NormalSy(i,N)=-Normal(0.25d0*(x(i+1,N-1)+x(i+1,N)+x(i,N)+x(i,N-1)), 0.25d0*(y(i+1,N-1)+y(i+1,N)+y(i,N)+y(i,N-1)), & 
                0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)),0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),2) 
         
         TangentWx(i,N)=-Tangent(0.5d0*(x(i-1,N)+x(i,N)),0.5d0*(y(i-1,N)+y(i,N)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N)+x(i,N-1)), & 
                  0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N)+y(i,N-1)),1)
         TangentWy(i,N)=-Tangent(0.5d0*(x(i-1,N)+x(i,N)),0.5d0*(y(i-1,N)+y(i,N)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N)+x(i,N-1)), &
                  0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N)+y(i,N-1)),2)
         TangentEx(i,N)=-Tangent(0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)), & 
                  0.5d0*(x(i+1,N)+x(i,N)), 0.5d0*(y(i+1,N)+y(i,N)), 1)
         TangentEy(i,N)=-Tangent(0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)), &
                  0.5d0*(x(i+1,N)+x(i,N)), 0.5d0*(y(i+1,N)+y(i,N)), 2)
         TangentNx(i,N)=-Tangent(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),1)
         TangentNy(i,N)=-Tangent(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),2)
         TangentSx(i,N)=-Tangent(0.25d0*(x(i-1,N-1)+x(i,N-1)+x(i-1,N)+x(i,N)), 0.25d0*(y(i-1,N-1)+y(i,N-1)+y(i-1,N)+y(i,N)),  &
                  0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)),1)
         TangentSy(i,N)=-Tangent(0.25d0*(x(i-1,N-1)+x(i,N-1)+x(i-1,N)+x(i,N)), 0.25d0*(y(i-1,N-1)+y(i,N-1)+y(i-1,N)+y(i,N)), & 
                  0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)),2)
         
         CurviNx(i,N)=CurviNx(i,N-1) !-Tangent(x(i,j),y(i,j),x(i,j+1),y(i,j+1),1)
         CurviNy(i,N)=CurviNy(i,N-1) !-Tangent(x(i,j),y(i,j),x(i,j+1),y(i,j+1),2)
         CurviSx(i,N)=-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),1)
         CurviSy(i,N)=-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),2)
         CurviEx(i,N)=-Tangent(x(i,N),y(i,N),x(i+1,N),y(i+1,N),1)
         CurviEy(i,N)=-Tangent(x(i,N),y(i,N),x(i+1,N),y(i+1,N),2)
         CurviWx(i,N)=-Tangent(x(i,N),y(i,N),x(i-1,N),y(i-1,N),1)
         CurviWy(i,N)=-Tangent(x(i,N),y(i,N),x(i-1,N),y(i-1,N),2)
         
         DistN(i,N)=0d0 !sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
         DistS(i,N)=sqrt((x(i,N)-x(i,N-1))**2+(y(i,N)-y(i,N-1))**2)
         DistE(i,N)=sqrt((x(i+1,N)-x(i,N))**2+(y(i+1,N)-y(i,N))**2)
         DistW(i,N)=sqrt((x(i,N)-x(i-1,N))**2+(y(i,N)-y(i-1,N))**2)
         
         DistDualN(i,N)=CellAreaN(i,N)
         DistDualS(i,N)=CellAreaS(i,N)
         DistDualE(i,N)=CellAreaE(i,N)
         DistDualW(i,N)=CellAreaW(i,N)
         
         GradNeX(i,N) = GradNeX(i,N-1)
         GradNeY(i,N) = GradNeY(i,N-1)
      
        ! SOUTH
         CellVol(i,1)=0.25d0*AreaElement(x(i+1,1),y(i+1,1),x(i+1,2),y(i+1,2),x(i-1,2),y(i-1,2),x(i-1,1),y(i-1,1))
         CellAreaN(i,1)=0.25d0*(x(i+1,2)**2+2d0*x(i+1,2)*x(i+1,1)-2d0*x(i+1,2)*x(i-1,1)-2d0*x(i+1,2)*x(i-1,2) & 
                 +x(i+1,1)**2-2d0*x(i+1,1)*x(i-1,1)-2d0*x(i+1,1)*x(i-1,2)+x(i-1,1)**2+2d0*x(i-1,1)*x(i-1,2)+x(i-1,2)**2 & 
                 +y(i+1,2)**2+2d0*y(i+1,2)*y(i+1,1)-2d0*y(i+1,2)*y(i-1,1)-2d0*y(i+1,2)*y(i-1,2)+y(i+1,1)**2 & 
                 -2d0*y(i+1,1)*y(i-1,1)-2d0*y(i+1,1)*y(i-1,2)+y(i-1,1)**2+2d0*y(i-1,1)*y(i-1,2)+y(i-1,2)**2)**(0.5d0)
         CellAreaS(i,1)=0d0 !sqrt((0.5d0*(x(i+1,1)+x(i,1))-0.5d0*(x(i-1,1)+x(i,1)))**2+(0.5d0*(y(i+1,1)+y(i,1))-0.5d0*(y(i-1,1)+y(i,1)))**2)
         CellAreaW(i,1)=0.25d0*(x(i,1)**2-2d0*x(i,1)*x(i,2)+2d0*x(i,1)*x(i-1,1)-2d0*x(i,1)*x(i-1,2)+x(i,2)**2 & 
                  -2d0*x(i,2)*x(i-1,1)+2d0*x(i,2)*x(i-1,2)+x(i-1,1)**2-2d0*x(i-1,1)*x(i-1,2)+x(i-1,2)**2+y(i,1)**2-2d0*y(i,1)*y(i,2) & 
                  +2d0*y(i,1)*y(i-1,1)-2d0*y(i,1)*y(i-1,2)+y(i,2)**2-2d0*y(i,2)*y(i-1,1)+2d0*y(i,2)*y(i-1,2)+y(i-1,1)**2 & 
                  -2d0*y(i-1,1)*y(i-1,2)+y(i-1,2)**2)**(0.5d0)
         CellAreaE(i,1)=0.25d0*(x(i+1,2)**2-2d0*x(i+1,2)*x(i+1,1)-2d0*x(i+1,2)*x(i,1)+2d0*x(i+1,2)*x(i,2)+x(i+1,1)**2+2d0*x(i+1,1)*x(i,1) & 
                  -2d0*x(i+1,1)*x(i,2)+x(i,1)**2-2d0*x(i,1)*x(i,2)+x(i,2)**2+y(i+1,2)**2-2d0*y(i+1,2)*y(i+1,1)-2d0*y(i+1,2)*y(i,1) & 
                  +2d0*y(i+1,2)*y(i,2)+y(i+1,1)**2+2d0*y(i+1,1)*y(i,1)-2d0*y(i+1,1)*y(i,2)+y(i,1)**2-2d0*y(i,1)*y(i,2)+y(i,2)**2)**(0.5d0)
         
         NormalEx(i,1)=-Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), & 
                0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),1) 
         NormalEy(i,1)=-Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), & 
                0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),2)
         NormalWx(i,1)=Normal(0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)), & 
                0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),1)
         NormalWy(i,1)=Normal(0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)), & 
                0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),2)
         NormalNx(i,1)=Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), & 
                0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),1)
         NormalNy(i,1)=Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)), & 
                0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),2)
         NormalSx(i,1)=Normal(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),1) 
         NormalSy(i,1)=Normal(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),2) 
         
         TangentNx(i,1)=Tangent(0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)), & 
                  0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),1)
         TangentNy(i,1)=Tangent(0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)), & 
                  0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),2)
         TangentSx(i,1)=-Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),1)
         TangentSy(i,1)=-Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),2)
         TangentEx(i,1)=-Tangent(0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)), & 
                  0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),1)
         TangentEy(i,1)=-Tangent(0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)), & 
                  0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),2)
         TangentWx(i,1)=Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)), &
                  0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),1)
         TangentWy(i,1)=Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)), & 
                  0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),2)
         
         CurviNx(i,1)=-Tangent(x(i,1),y(i,1),x(i,2),y(i,2),1)
         CurviNy(i,1)=-Tangent(x(i,1),y(i,1),x(i,2),y(i,2),2)
         CurviSx(i,1)=CurviSx(i,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),1)
         CurviSy(i,1)=CurviSy(i,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),2)
         CurviEx(i,1)=-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),1)
         CurviEy(i,1)=-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),2)
         CurviWx(i,1)=-Tangent(x(i,1),y(i,1),x(i-1,1),y(i-1,1),1)
         CurviWy(i,1)=-Tangent(x(i,1),y(i,1),x(i-1,1),y(i-1,1),2)
         
         DistN(i,1)=sqrt( (x(i,2)-x(i,1))**2 + (y(i,2)-y(i,1))**2 )
         DistS(i,1)=0d0
         DistE(i,1)=sqrt( (x(i+1,1)-x(i,1))**2 + (y(i+1,1)-y(i,1))**2 )
         DistW(i,1)=sqrt( (x(i,1)-x(i-1,1))**2 + (y(i,1)-y(i-1,1))**2 )
         
         DistDualN(i,1)=CellAreaN(i,1)
         DistDualS(i,1)=CellAreaS(i,1)
         DistDualE(i,1)=CellAreaE(i,1)
         DistDualW(i,1)=CellAreaW(i,1)
         
      end do
      
        do j=2,N-1
!         
          ! WEST
          CellVol(1,j)=0.25d0*AreaElement(x(1,j-1),y(1,j-1),x(2,j-1),y(2,j-1),x(2,j+1),y(2,j+1),x(1,j+1),y(1,j+1))
          CellAreaN(1,j)=0.25d0*(x(2,j+1)**2+2d0*x(2,j+1)*x(2,j)-2d0*x(2,j+1)*x(1,j)-2d0*x(2,j+1)*x(1,j+1)+x(2,j)**2 & 
                  -2d0*x(2,j)*x(1,j)-2d0*x(2,j)*x(1,j+1)+x(1,j)**2+2d0*x(1,j)*x(1,j+1)+x(1,j+1)**2+y(2,j+1)**2 & 
                  +2d0*y(2,j+1)*y(2,j)-2d0*y(2,j+1)*y(1,j)-2d0*y(2,j+1)*y(1,j+1)+y(2,j)**2-2d0*y(2,j)*y(1,j) & 
                  -2d0*y(2,j)*y(1,j+1)+y(1,j)**2+2d0*y(1,j)*y(1,j+1)+y(1,j+1)**2)**0.5d0
          CellAreaS(1,j)=0.25d0*(x(2,j)**2+2d0*x(2,j)*x(2,j-1)-2d0*x(2,j)*x(1,j)-2d0*x(2,j)*x(1,j-1)+x(2,j-1)**2 & 
                  -2d0*x(2,j-1)*x(1,j)-2d0*x(2,j-1)*x(1,j-1)+x(1,j)**2+2d0*x(1,j)*x(1,j-1)+x(1,j-1)**2+y(2,j)**2 &
                  +2d0*y(2,j)*y(2,j-1)-2d0*y(2,j)*y(1,j)-2d0*y(2,j)*y(1,j-1)+y(2,j-1)**2-2d0*y(2,j-1)*y(1,j) & 
                  -2d0*y(2,j-1)*y(1,j-1)+y(1,j)**2+2d0*y(1,j)*y(1,j-1)+y(1,j-1)**2)**0.5d0
          CellAreaW(1,j)=sqrt( (0.5d0*(x(1,j+1)+x(1,j))-0.5d0*(x(1,j-1)+x(1,j)))**2 + (0.5d0*(y(1,j+1)+y(1,j)) & 
                        -0.5d0*(y(1,j-1)+y(1,j)))**2 )
          CellAreaE(1,j)=0.25d0*(x(2,j+1)**2+2d0*x(2,j+1)*x(1,j+1)-2d0*x(2,j+1)*x(2,j-1)-2d0*x(2,j+1)*x(1,j-1) & 
                  +x(1,j+1)**2-2d0*x(1,j+1)*x(2,j-1)-2d0*x(1,j+1)*x(1,j-1)+x(2,j-1)**2+2d0*x(2,j-1)*x(1,j-1) & 
                  +x(1,j-1)**2+y(2,j+1)**2+2d0*y(2,j+1)*y(1,j+1)-2d0*y(2,j+1)*y(2,j-1)-2d0*y(2,j+1)*y(1,j-1) & 
                  +y(1,j+1)**2-2d0*y(1,j+1)*y(2,j-1)-2d0*y(1,j+1)*y(1,j-1)+y(2,j-1)**2+2d0*y(2,j-1)*y(1,j-1) & 
                  +y(1,j-1)**2)**0.5d0
          
          NormalEx(1,j)=-Normal(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)), & 
                        0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),1)
          NormalEy(1,j)=-Normal(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)), & 
                        0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),2)
          NormalNx(1,j)=Normal(0.25d0*(x(2,j+1)+x(1,j+1)+x(2,j)+x(1,j)),0.25d0*(y(2,j+1)+y(1,j+1)+y(2,j)+y(1,j)), & 
                      0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),1)
          NormalNy(1,j)=Normal(0.25d0*(x(2,j+1)+x(1,j+1)+x(2,j)+x(1,j)),0.25d0*(y(2,j+1)+y(1,j+1)+y(2,j)+y(1,j)), & 
                      0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),2)
          NormalSx(1,j)=-Normal(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)), & 
                        0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),1) 
          NormalSy(1,j)=-Normal(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)), & 
                        0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),2) 
          NormalWx(1,j)=Normal(0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),0.5d0*(x(1,j-1)+x(1,j)), & 
                        0.5d0*(y(1,j-1)+y(1,j)),1)
          NormalWy(1,j)=Normal(0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),0.5d0*(x(1,j-1)+x(1,j)), & 
                        0.5d0*(y(1,j-1)+y(1,j)),2)
          
!           write(*,*) "[Debug]", NormalWx(1,j)**2+NormalWy(1,j)**2
          
          TangentNx(1,j)=-Tangent(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) & 
                  +y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),1)
          TangentNy(1,j)=-Tangent(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) & 
                  +y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),2)
          TangentSx(1,j)=-Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.25d0*(x(1,j-1)+x(1,j) & 
                  +x(2,j-1)+x(2,j)),0.25d0*(y(1,j-1)+y(1,j)+y(2,j-1)+y(2,j)),1)
          TangentSy(1,j)=-Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.25d0*(x(1,j-1)+x(1,j) & 
                  +x(2,j-1)+x(2,j)),0.25d0*(y(1,j-1)+y(1,j)+y(2,j-1)+y(2,j)),2)
          TangentWx(1,j)=Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)), & 
                  0.5d0*(y(1,j+1)+y(1,j)),1)
          TangentWy(1,j)=Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)), & 
                  0.5d0*(y(1,j+1)+y(1,j)),2)
          TangentEx(1,j)=-Tangent(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1) & 
                  +y(1,j)),0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) & 
                  +y(1,j)),1)
          TangentEy(1,j)=-Tangent(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1) & 
                  +y(1,j)),0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1) &
                  +y(1,j)),2)
          
          CurviNx(1,j)=-Tangent(x(1,j),y(1,j),x(1,j+1),y(1,j+1),1)
          CurviNy(1,j)=-Tangent(x(1,j),y(1,j),x(1,j+1),y(1,j+1),2)
          CurviSx(1,j)=-Tangent(x(1,j),y(1,j),x(1,j-1),y(1,j-1),1)
          CurviSy(1,j)=-Tangent(x(1,j),y(1,j),x(1,j-1),y(1,j-1),2)
          CurviEx(1,j)=-Tangent(x(1,j),y(1,j),x(2,j),y(2,j),1)
          CurviEy(1,j)=-Tangent(x(1,j),y(1,j),x(2,j),y(2,j),2)
          CurviWx(1,j)=CurviWx(2,j) !-Tangent(x(1,j),y(1,j),x(i-1,N),y(i-1,N),1)
          CurviWy(1,j)=CurviWy(2,j) !-Tangent(x(1,j),y(1,j),x(i-1,N),y(i-1,N),2)
          
          DistN(1,j)=sqrt( (x(1,j+1)-x(1,j))**2 + (y(1,j+1)-y(1,j))**2 )
          DistS(1,j)=sqrt( (x(1,j)-x(1,j-1))**2 + (y(1,j)-y(1,j-1))**2 )
          DistE(1,j)=sqrt( (x(2,j)-x(1,j))**2 + (y(2,j)-y(1,j))**2 )
          DistW(1,j)=0d0
          
          DistDualN(1,j)=CellAreaN(1,j)
          DistDualS(1,j)=CellAreaS(1,j)
          DistDualE(1,j)=CellAreaE(1,j)
          DistDualW(1,j)=CellAreaW(1,j)
          
          CellVol(M,j)=0.25d0*AreaElement(x(M-1,j-1),y(M-1,j-1),x(M,j-1),y(M,j-1),x(M,j+1),y(M,j+1), & 
                      x(M-1,j+1),y(M-1,j+1))
          CellAreaN(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M,j)-2d0*x(M,j+1)*x(M-1,j)-2d0*x(M,j+1)*x(M-1,j+1) & 
                      +x(M,j)**2-2d0*x(M,j)*x(M-1,j)-2d0*x(M,j)*x(M-1,j+1)+x(M-1,j)**2+2d0*x(M-1,j)*x(M-1,j+1) & 
                      +x(M-1,j+1)**2+y(M,j+1)**2+2d0*y(M,j+1)*y(M,j)-2d0*y(M,j+1)*y(M-1,j)-2d0*y(M,j+1)*y(M-1,j+1) & 
                      +y(M,j)**2-2d0*y(M,j)*y(M-1,j)-2d0*y(M,j)*y(M-1,j+1)+y(M-1,j)**2+2d0*y(M-1,j)*y(M-1,j+1) & 
                      +y(M-1,j+1)**2)**0.5d0
          CellAreaS(M,j)=0.25d0*(x(M,j)**2+2d0*x(M,j)*x(M,j-1)-2d0*x(M,j)*x(M-1,j-1)-2d0*x(M,j)*x(M-1,j) & 
                  +x(M,j-1)**2-2d0*x(M-1,j-1)*x(M,j-1)-2d0*x(M,j-1)*x(M-1,j)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M-1,j) & 
                  +x(M-1,j)**2+y(M,j)**2+2d0*y(M,j)*y(M,j-1)-2d0*y(M,j)*y(M-1,j-1)-2d0*y(M,j)*y(M-1,j)+y(M,j-1)**2 &
                  -2d0*y(M-1,j-1)*y(M,j-1)-2d0*y(M,j-1)*y(M-1,j)+y(M-1,j-1)**2+2d0*y(M-1,j-1)*y(M-1,j) & 
                  +y(M-1,j)**2)**0.5d0
          CellAreaW(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M-1,j+1)-2d0*x(M,j+1)*x(M-1,j-1)-2d0*x(M,j+1)*x(M,j-1) & 
                  +x(M-1,j+1)**2-2d0*x(M-1,j+1)*x(M-1,j-1)-2d0*x(M-1,j+1)*x(M,j-1)+x(M-1,j-1)**2 & 
                  +2d0*x(M-1,j-1)*x(M,j-1)+x(M,j-1)**2+y(M,j+1)**2+2d0*y(M,j+1)*y(M-1,j+1) & 
                  -2d0*y(M,j+1)*y(M-1,j-1)-2d0*y(M,j+1)*y(M,j-1)+y(M-1,j+1)**2-2d0*y(M-1,j+1)*y(M-1,j-1) & 
                  -2d0*y(M-1,j+1)*y(M,j-1)+y(M-1,j-1)**2+2d0*y(M-1,j-1)*y(M,j-1)+y(M,j-1)**2)**0.5d0
                  
          CellAreaE(M,j)=sqrt( (0.5d0*(x(M,j+1)+x(M,j))-0.5d0*(x(M,j-1)+x(M,j)) )**2 + (0.5d0*(y(M,j+1) & 
                    +y(M,j))-0.5d0*(y(M,j-1)+y(M,j)) )**2)
          
          TangentNx(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1) & 
                  +x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1)
          TangentNy(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1) & 
                  +x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2)
          TangentEx(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1) & 
                  +x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),1)
          TangentEy(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)), & 
                  0.5d0*(y(M,j+1)+y(M,j)),2)
          TangentSx(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) & 
                  +y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),1)
          TangentSy(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) & 
                  +y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),2)
          TangentWx(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) & 
                  +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1) &
                  +y(M-1,j)+y(M,j-1)+y(M,j)),1)
          TangentWy(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) & 
                  +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1) & 
                  +y(M-1,j)+y(M,j-1)+y(M,j)),2)
          
!           if(y(1,j)<0d0) then
!             TangentNx(M,j)=-TangentNx(M,j)
!             TangentNy(M,j)=-TangentNy(M,j)
!             TangentSx(M,j)=-TangentSx(M,j)
!             TangentSy(M,j)=-TangentSy(M,j)
!             TangentWx(M,j)=-TangentWx(M,j)
!             TangentWy(M,j)=-TangentWy(M,j)
!           end if
          
          NormalEx(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)), & 
                        0.5d0*(y(M,j-1)+y(M,j)),1); 
          NormalEy(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)), & 
                        0.5d0*(y(M,j-1)+y(M,j)),2); 
          NormalWx(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) & 
                        +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)), & 
                        0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1) 
          NormalWy(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)), & 
                       0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),2)
          NormalNx(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)), & 
                       0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1) 
          NormalNy(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)), & 
                       0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2) 
          NormalSx(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1)+x(M-1,j)+x(M,j)), &
                        0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),1) 
          NormalSy(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1)+x(M-1,j)+x(M,j)), & 
                        0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),2) 
          
          CurviNx(M,j)=-Tangent(x(M,j),y(M,j),x(M,j+1),y(M,j+1),1)
          CurviNy(M,j)=-Tangent(x(M,j),y(M,j),x(M,j+1),y(M,j+1),2)
          CurviSx(M,j)=-Tangent(x(M,j),y(M,j),x(M,j-1),y(M,j-1),1)
          CurviSy(M,j)=-Tangent(x(M,j),y(M,j),x(M,j-1),y(M,j-1),2)
          CurviEx(M,j)=CurviEx(M-1,j) !Tangent(x(M,j),y(M,j),x(i+1,N),y(i+1,N),1)
          CurviEy(M,j)=CurviEy(M-1,j) !-Tangent(x(M,j),y(M,j),x(i+1,N),y(i+1,N),2)
          CurviWx(M,j)=-Tangent(x(M,j),y(M,j),x(M-1,j),y(M-1,j),1)
          CurviWy(M,j)=-Tangent(x(M,j),y(M,j),x(M-1,j),y(M-1,j),2)
          
          DistN(M,j)=sqrt( (x(M,j+1)-x(M,j))**2 + (y(M,j+1)-y(M,j))**2 ) 
          DistS(M,j)=sqrt( (x(M,j-1)-x(M,j))**2 + (y(M,j-1)-y(M,j))**2 )
          DistE(M,j)=0d0
          DistW(M,j)=sqrt( (x(M-1,j)-x(M,j))**2 + (y(M-1,j)-y(M,j))**2 )
          
          DistDualN(M,j)=CellAreaN(M,j)
          DistDualS(M,j)=CellAreaS(M,j)
          DistDualE(M,j)=CellAreaE(M,j)
          DistDualW(M,j)=CellAreaW(M,j)
          
          ! EAST
          CellVol(M,j)=0.25d0*AreaElement(x(M-1,j-1),y(M-1,j-1),x(M,j-1),y(M,j-1),x(M,j+1),y(M,j+1),x(M-1,j+1),y(M-1,j+1))
          CellAreaN(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M,j)-2d0*x(M,j+1)*x(M-1,j)-2d0*x(M,j+1)*x(M-1,j+1)+x(M,j)**2 &
                  -2d0*x(M,j)*x(M-1,j)-2d0*x(M,j)*x(M-1,j+1)+x(M-1,j)**2+2d0*x(M-1,j)*x(M-1,j+1)+x(M-1,j+1)**2 & 
                  +y(M,j+1)**2+2d0*y(M,j+1)*y(M,j)-2d0*y(M,j+1)*y(M-1,j)-2d0*y(M,j+1)*y(M-1,j+1)+y(M,j)**2 & 
                  -2d0*y(M,j)*y(M-1,j)-2d0*y(M,j)*y(M-1,j+1)+y(M-1,j)**2+2d0*y(M-1,j)*y(M-1,j+1)+y(M-1,j+1)**2)**0.5d0
          CellAreaS(M,j)=0.25d0*(x(M,j)**2+2d0*x(M,j)*x(M,j-1)-2d0*x(M,j)*x(M-1,j-1)-2d0*x(M,j)*x(M-1,j)+x(M,j-1)**2 & 
                  -2d0*x(M-1,j-1)*x(M,j-1)-2d0*x(M,j-1)*x(M-1,j)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M-1,j)+x(M-1,j)**2 & 
                  +y(M,j)**2+2d0*y(M,j)*y(M,j-1)-2d0*y(M,j)*y(M-1,j-1)-2d0*y(M,j)*y(M-1,j)+y(M,j-1)**2 & 
                  -2d0*y(M-1,j-1)*y(M,j-1)-2d0*y(M,j-1)*y(M-1,j)+y(M-1,j-1)**2+2d0*y(M-1,j-1)*y(M-1,j) & 
                  +y(M-1,j)**2)**0.5d0
          CellAreaW(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M-1,j+1)-2d0*x(M,j+1)*x(M-1,j-1)-2d0*x(M,j+1)*x(M,j-1) & 
                  +x(M-1,j+1)**2-2d0*x(M-1,j+1)*x(M-1,j-1)-2d0*x(M-1,j+1)*x(M,j-1)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M,j-1) & 
                  +x(M,j-1)**2+y(M,j+1)**2+2d0*y(M,j+1)*y(M-1,j+1)-2d0*y(M,j+1)*y(M-1,j-1)-2d0*y(M,j+1)*y(M,j-1) & 
                  +y(M-1,j+1)**2-2d0*y(M-1,j+1)*y(M-1,j-1)-2d0*y(M-1,j+1)*y(M,j-1)+y(M-1,j-1)**2 & 
                  +2d0*y(M-1,j-1)*y(M,j-1)+y(M,j-1)**2)**0.5d0
          CellAreaE(M,j)=sqrt( (0.5d0*(x(M,j+1)+x(M,j))-0.5d0*(x(M,j-1)+x(M,j)) )**2 + (0.5d0*(y(M,j+1)+y(M,j)) & 
                        -0.5d0*(y(M,j-1)+y(M,j)) )**2)
          
          TangentNx(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) & 
                          +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1)
          TangentNy(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) & 
                          +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2)
          TangentEx(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)), & 
                          0.5d0*(y(M,j+1)+y(M,j)),1)
          TangentEy(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)), & 
                          0.5d0*(y(M,j+1)+y(M,j)),2)
          TangentSx(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) & 
                          +y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),1)
          TangentSy(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) & 
                          +y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),2)
          TangentWx(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) & 
                          +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)), & 
                          0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1)
          TangentWy(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) &
                          +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)), &
                          0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),2)
          
!           if(y(1,j)<0d0) then
!             TangentNx(M,j)=-TangentNx(M,j)
!             TangentNy(M,j)=-TangentNy(M,j)
!             TangentSx(M,j)=-TangentSx(M,j)
!             TangentSy(M,j)=-TangentSy(M,j)
!             TangentWx(M,j)=-TangentWx(M,j)
!             TangentWy(M,j)=-TangentWy(M,j)
!           end if
          
          NormalEx(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)), & 
                        0.5d0*(y(M,j-1)+y(M,j)),1); 
          NormalEy(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)), & 
                        0.5d0*(y(M,j-1)+y(M,j)),2); 
          NormalWx(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) & 
                       +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)), & 
                      0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1) 
          NormalWy(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j) & 
          +y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j) & 
          +y(M,j-1)+y(M,j)),2)
          NormalNx(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1) & 
                      +x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1) 
          NormalNy(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j) & 
                      +x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2) 
          NormalSx(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1) &
                      +x(M,j-1)+x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),1) 
          NormalSy(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1) & 
                      +x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),2) 
          
          CurviNx(M,j)=-Tangent(x(M,j),y(M,j),x(M,j+1),y(M,j+1),1)
          CurviNy(M,j)=-Tangent(x(M,j),y(M,j),x(M,j+1),y(M,j+1),2)
          CurviSx(M,j)=-Tangent(x(M,j),y(M,j),x(M,j-1),y(M,j-1),1)
          CurviSy(M,j)=-Tangent(x(M,j),y(M,j),x(M,j-1),y(M,j-1),2)
          CurviEx(M,j)=CurviEx(M-1,j) !Tangent(x(M,j),y(M,j),x(i+1,N),y(i+1,N),1)
          CurviEy(M,j)=CurviEy(M-1,j) !-Tangent(x(M,j),y(M,j),x(i+1,N),y(i+1,N),2)
          CurviWx(M,j)=-Tangent(x(M,j),y(M,j),x(M-1,j),y(M-1,j),1)
          CurviWy(M,j)=-Tangent(x(M,j),y(M,j),x(M-1,j),y(M-1,j),2)
          
          DistN(M,j)=sqrt( (x(M,j+1)-x(M,j))**2 + (y(M,j+1)-y(M,j))**2 ) 
          DistS(M,j)=sqrt( (x(M,j-1)-x(M,j))**2 + (y(M,j-1)-y(M,j))**2 )
          DistE(M,j)=0d0
          DistW(M,j)=sqrt( (x(M-1,j)-x(M,j))**2 + (y(M-1,j)-y(M,j))**2 )
          
          DistDualN(M,j)=CellAreaN(M,j)
          DistDualS(M,j)=CellAreaS(M,j)
          DistDualE(M,j)=CellAreaE(M,j)
          DistDualW(M,j)=CellAreaW(M,j)
      
      end do
      
      !North-East
        CellVol(M,N)=AreaElement(0.25d0*(x(M-1, N-1)+x(M, N-1)+x(M, N)+x(M-1, N)),0.25d0*(y(M-1,N-1)+y(M-1,N) & 
                      +y(M,N)+y(M,N-1)), 0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)), x(M,N),y(M,N), &
                      0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)))
!         CellVol(M,N)=0.25d0*AreaElement(x(M-1,N-1),y(M-1,N-1),x(M,N-1),y(M,N-1), x(M,N), y(M,N), x(M-1, N), y(M-1, N))
        CellAreaE(M,N)=sqrt((x(M,N)-0.5d0*(x(M,N)+x(M,N-1)))**2+(y(M,N)-0.5d0*(y(M,N)+y(M,N-1)))**2)
        CellAreaN(M,N)=sqrt((x(M,N)-0.5d0*(x(M-1,N)+x(M,N)))**2+(y(M,N)-0.5d0*(y(M-1,N)+y(M,N)))**2)
        CellAreaS(M,N)=0.25d0*(x(M,N)**2+2d0*x(M,N)*x(M,N-1)-2d0*x(M,N)*x(M-1,N-1)-2d0*x(M-1,N)*x(M,N)+x(M,N-1)**2 &
                -2d0*x(M-1,N-1)*x(M,N-1)-2d0*x(M-1,N)*x(M,N-1)+x(M-1,N-1)**2+2d0*x(M-1,N)*x(M-1,N-1)+x(M-1,N)**2 &
                +y(M,N)**2+2d0*y(M,N)*y(M,N-1)-2d0*y(M,N)*y(M-1,N-1)-2d0*y(M-1,N)*y(M,N)+y(M,N-1)**2 & 
                -2d0*y(M-1,N-1)*y(M,N-1)-2d0*y(M-1,N)*y(M,N-1)+y(M-1,N-1)**2+2d0*y(M-1,N)*y(M-1,N-1) & 
                +y(M-1,N)**2)**(0.5d0)
        CellAreaW(M,N)=0.25d0*(x(M-1,N)**2+2d0*x(M-1,N)*x(M,N)-2d0*x(M-1,N)*x(M-1,N-1)-2d0*x(M-1,N)*x(M,N-1)+x(M,N)**2 &
                -2d0*x(M,N)*x(M-1,N-1)-2d0*x(M,N)*x(M,N-1)+x(M-1,N-1)**2+2d0*x(M-1,N-1)*x(M,N-1)+x(M,N-1)**2 &
                +y(M-1,N)**2+2d0*y(M-1,N)*y(M,N)-2d0*y(M-1,N)*y(M-1,N-1)-2d0*y(M-1,N)*y(M,N-1)+y(M,N)**2 &
                -2d0*y(M,N)*y(M-1,N-1)-2d0*y(M,N)*y(M,N-1)+y(M-1,N-1)**2+2d0*y(M-1,N-1)*y(M,N-1) &
                +y(M,N-1)**2)**(0.5d0)
        NormalEx(M,N)=-Normal(x(M,N),y(M,N),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),1)
        NormalEy(M,N)=-Normal(x(M,N),y(M,N),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),2)
        NormalNx(M,N)=Normal(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),1)
        NormalNy(M,N)=Normal(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),2)
        NormalWx(M,N)=Normal(0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)),0.25d0*(x(M-1, N-1)+x(M, N-1) & 
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),1) 
        NormalWy(M,N)=Normal(0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),2)
        NormalSx(M,N)=-Normal(0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),1)
        NormalSy(M,N)=-Normal(0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)),0.25d0*(x(M-1, N-1)+x(M, N-1) &
                           +x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),2)
        
        TangentNx(M,N)=-Tangent(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),1)
        TangentNy(M,N)=-Tangent(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),2)
        TangentEx(M,N)=-Tangent(0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),x(M,N),y(M,N),1)
        TangentEy(M,N)=-Tangent(0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),x(M,N),y(M,N),2)
        TangentWx(M,N)=-Tangent(0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),0.25d0*(x(M-1,N-1) & 
                        +x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),1)
        TangentWy(M,N)=-Tangent(0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),0.25d0*(x(M-1,N-1)+x(M-1,N) &
                        +x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),2)
        TangentSx(M,N)=-Tangent(0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N) &
                        +y(M,N-1)+y(M,N)),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),1)
        TangentSy(M,N)=-Tangent(0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N) &
                        +y(M,N-1)+y(M,N)),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),2)
        
        CurviNx(M,N)=CurviNx(M,N-1)
        CurviNy(M,N)=CurviNy(M,N-1)
        CurviSx(M,N)=-Tangent(x(M,N),y(M,N),x(M,N-1),y(M,N-1),1)
        CurviSy(M,N)=-Tangent(x(M,N),y(M,N),x(M,N-1),y(M,N-1),2)
        CurviEx(M,N)=CurviEx(M-1,N)
        CurviEy(M,N)=CurviEy(M-1,N)
        CurviWx(M,N)=-Tangent(x(M,N),y(M,N),x(M-1,N),y(M-1,N),1)
        CurviWy(M,N)=-Tangent(x(M,N),y(M,N),x(M-1,N),y(M-1,N),2)
        
        DistN(M,N)=0d0
        DistS(M,N)=sqrt((x(M,N)-x(M,N-1))**2+(y(M,N)-y(M,N-1))**2)
        DistE(M,N)=0d0
        DistW(M,N)=sqrt((x(M,N)-x(M-1,N))**2+(y(M,N)-y(M-1,N))**2)
        
        DistDualN(M,N)=CellAreaN(M,N)
        DistDualS(M,N)=CellAreaS(M,N)
        DistDualE(M,N)=CellAreaE(M,N)
        DistDualW(M,N)=CellAreaW(M,N)
      
      !South-East
        CellVol(M,1)=AreaElement(0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)), x(M,1),y(M,1), &
                      0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)), 0.25d0*(x(M-1, 1)+x(M-1, 2) &
                      +x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)))
!         CellVol(M,1)=0.25d0*AreaElement(x(M-1,1),y(M-1,1),x(M,1),y(M,1),x(M,2),y(M,2),x(M-1,2),y(M-1,2))
        CellAreaS(M,1)=sqrt(( x(M,1)-0.5d0*(x(M,1)+x(M-1,1)) )**2+( y(M,1)-0.5d0*(y(M,1)+y(M-1,1)) )**2)
        CellAreaE(M,1)=sqrt( ( x(M,1) - 0.5d0*(x(M,1)+x(M,2)) )**2 + ( y(M,1) - 0.5d0*(y(M,1)+y(M,2)) )**2 )
        CellAreaN(M,1)=0.25d0*(x(M,1)**2+2d0*x(M,1)*x(M,2)-2d0*x(M-1,1)*x(M,1)-2d0*x(M-1,2)*x(M,1)+x(M,2)**2 &
                -2d0*x(M-1,1)*x(M,2)-2d0*x(M-1,2)*x(M,2)+x(M-1,1)**2+2d0*x(M-1,1)*x(M-1,2)+x(M-1,2)**2 &
                +y(M,1)**2+2d0*y(M,1)*y(M,2)-2d0*y(M-1,1)*y(M,1)-2d0*y(M-1,2)*y(M,1)+y(M,2)**2 &
                -2d0*y(M-1,1)*y(M,2)-2d0*y(M-1,2)*y(M,2)+y(M-1,1)**2+2d0*y(M-1,1)*y(M-1,2) &
                +y(M-1,2)**2)**(0.5d0)
        CellAreaW(M,1)=0.25d0*(x(M-1,1)**2-2d0*x(M-1,1)*x(M-1,2)+2d0*x(M-1,1)*x(M,1)-2d0*x(M-1,1)*x(M,2) &
                +x(M-1,2)**2-2d0*x(M-1,2)*x(M,1)+2d0*x(M-1,2)*x(M,2)+x(M,1)**2-2d0*x(M,1)*x(M,2)+x(M,2)**2 &
                +y(M-1,1)**2-2d0*y(M-1,1)*y(M-1,2)+2d0*y(M-1,1)*y(M,1)-2d0*y(M-1,1)*y(M,2)+y(M-1,2)**2 &
                -2d0*y(M-1,2)*y(M,1)+2d0*y(M-1,2)*y(M,2)+y(M,1)**2-2d0*y(M,1)*y(M,2)+y(M,2)**2)**(0.5d0)
        NormalSx(M,1)=Normal(0.5d0*(x(M-1,1)+x(M,1)),0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 1)
        NormalSy(M,1)=Normal(0.5d0*(x(M-1,1)+x(M,1)),0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 2)
        NormalEx(M,1)=Normal(x(M,1), y(M,1), 0.5d0*(x(M,1)+x(M,2)), 0.5d0*(y(M,1)+y(M,2)), 1)
        NormalEy(M,1)=Normal(x(M,1), y(M,1), 0.5d0*(x(M,1)+x(M,2)), 0.5d0*(y(M,1)+y(M,2)), 2)
        NormalWx(M,1)=Normal(0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1) &
                    +y(M, 2)),0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)),1); 
        NormalWy(M,1)=Normal(0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1) &
                    +y(M, 2)),0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)),2);
        NormalNx(M,1)=Normal(0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)),0.25d0*(x(M-1, 1)+x(M-1, 2) &
                    +x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),1); 
        NormalNy(M,1)=Normal(0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)),0.25d0*(x(M-1, 1)+x(M-1, 2) &
                    +x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),2);
        
        TangentNx(M,1)=-Tangent(0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),0.25d0*(x(M,1)+x(M,2)+x(M-1,1) &
                    +x(M-1,2)),0.25d0*(y(M,1)+y(M,2)+y(M-1,1)+y(M-1,2)),1)
        TangentNy(M,1)=-Tangent(0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),0.25d0*(x(M,1)+x(M,2)+x(M-1,1) &
                    +x(M-1,2)),0.25d0*(y(M,1)+y(M,2)+y(M-1,1)+y(M-1,2)),2)
        TangentWx(M,1)=-Tangent(0.25d0*(x(M-1,2)+x(M-1,1)+x(M,2)+x(M,1)),0.25d0*(y(M-1,2)+y(M-1,1)+y(M,2) &
                    +y(M,1)),0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)),1)
        TangentWy(M,1)=-Tangent(0.25d0*(x(M-1,2)+x(M-1,1)+x(M,2)+x(M,1)),0.25d0*(y(M-1,2)+y(M-1,1)+y(M,2) &
                    +y(M,1)),0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)),2)
        TangentEx(M,1)=-Tangent(x(M,1),y(M,1),0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),1)
        TangentEy(M,1)=-Tangent(x(M,1),y(M,1),0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),2)
        TangentSx(M,1)=-Tangent(0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 1)
        TangentSy(M,1)=-Tangent(0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 2)

         CurviNx(M,1)=-Tangent(x(M,1),y(M,1),x(M,2),y(M,2),1)
         CurviNy(M,1)=-Tangent(x(M,1),y(M,1),x(M,2),y(M,2),2)
         CurviSx(M,1)=CurviSx(M,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),1)
         CurviSy(M,1)=CurviSy(M,2) !-Tangent(x(i,N),y(i,N),x(i,N-1),y(i,N-1),2)
         CurviEx(M,1)=CurviEx(M-1,1) !-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),1)
         CurviEy(M,1)=CurviEy(M-1,1) !-Tangent(x(i,1),y(i,1),x(i+1,1),y(i+1,1),2)
         CurviWx(M,1)=-Tangent(x(M,1),y(M,1),x(M-1,1),y(M-1,1),1)
         CurviWy(M,1)=-Tangent(x(M,1),y(M,1),x(M-1,1),y(M-1,1),2)
         
         DistN(M,1)=sqrt( (x(M,2)-x(M,1))**2 + (y(M,2)-y(M,1))**2 )
         DistS(M,1)=0d0
         DistE(M,1)=0d0 !sqrt( (x(i+1,1)-x(i,1))**2 + (y(i+1,1)-y(i,1))**2 )
         DistW(M,1)=sqrt( (x(M,1)-x(M-1,1))**2 + (y(M,1)-y(M-1,1))**2 )
         
         DistDualN(M,1)=CellAreaN(M,1)
         DistDualS(M,1)=CellAreaS(M,1)
         DistDualE(M,1)=CellAreaE(M,1)
         DistDualW(M,1)=CellAreaW(M,1)
      
      !South-West
        CellVol(1,1)=AreaElement(x(1,1),y(1,1),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)), &
                      0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2)+y(2, 1) &
                      +y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)))
!         CellVol(1,1)=0.25d0*AreaElement(x(1,1),y(1,1),x(2,1),y(2,1),x(2,2),y(2,2),x(1,2),y(1,2))
        
!         CellVol(1,1)=2d0*CellVol(1,1)

        CellAreaS(1,1)=sqrt((x(1,1)-0.5d0*(x(2,1)+x(1,1)))**2+(y(1,1)-0.5d0*(y(2,1)+y(1,1)))**2) 
        CellAreaW(1,1)=sqrt((x(1,1)-0.5d0*(x(1,2)+x(1,1)))**2+(y(1,1)-0.5d0*(y(1,2)+y(1,1)))**2)
        CellAreaN(1,1)=sqrt((0.25d0*(x(2,2)+x(2,1)+x(1,2)+x(1,1))-0.5d0*(x(1,1)+x(1,2)))**2 &
                +(0.25d0*(y(2,2)+y(2,1)+y(1,2)+y(1,1))-0.5d0*(y(1,2)+y(1,1)))**2) !sqrt((0.5d0*(x(1,2)+x(1,1))-0.25d0*(x(2,2)+x(2,1)+x(1,1)+x(1,2)))**2d0+(0.5d0*(y(1,2)+y(1,1))-0.25d0*(y(2,2)+y(2,1)+y(1,1)+y(1,2)))**2d0) !0.25d0*(x(1,1)**2+2d0*x(1,1)*x(1,2)-2d0*x(1,1)*x(2,1)-2d0*x(1,1)*x(2,2)+x(1,2)**2-2d0*x(1,2)*x(2,1)-2d0*x(1,2)*x(2,2)+x(2,1)**2+2d0*x(2,1)*x(2,2)+x(2,2)**2+y(1,1)**2+2d0*y(1,1)*y(1,2)-2d0*y(1,1)*y(2,1)-2d0*y(1,1)*y(2,2)+y(1,2)**2-2d0*y(1,2)*y(2,1)-2d0*y(1,2)*y(2,2)+y(2,1)**2+2d0*y(2,1)*y(2,2)+y(2,2)**2)**(0.5d0);
        CellAreaE(1,1)=sqrt((0.25d0*(x(2,2)+x(2,1)+x(1,2)+x(1,1))-0.5d0*(x(2,1)+x(1,1)))**2 &
                +(0.25d0*(y(2,2)+y(2,1)+y(1,2)+y(1,1))-0.5d0*(y(2,1)+y(1,1)))**2) !0.25d0*(x(1,1)**2-2d0*x(1,1)*x(1,2)+2d0*x(1,1)*x(2,1)-2d0*x(1,1)*x(2,2)+x(1,2)**2-2d0*x(1,2)*x(2,1)+2d0*x(1,2)*x(2,2)+x(2,1)**2-2d0*x(2,1)*x(2,2)+x(2,2)**2+y(1,1)**2-2d0*y(1,1)*y(1,2)+2d0*y(1,1)*y(2,1)-2d0*y(1,1)*y(2,2)+y(1,2)**2-2d0*y(1,2)*y(2,1)+2d0*y(1,2)*y(2,2)+y(2,1)**2-2d0*y(2,1)*y(2,2)+y(2,2)**2)**(0.5d0)
        
        NormalSx(1,1)=Normal(x(1,1),y(1,1),0.5d0*(x(2,1)+x(1,1)), 0.5d0*(y(2,1)+y(1,1)), 1)
        NormalSy(1,1)=Normal(x(1,1),y(1,1),0.5d0*(x(2,1)+x(1,1)), 0.5d0*(y(2,1)+y(1,1)), 2)
        NormalWx(1,1)=-Normal(x(1,1),y(1,1),0.5d0*(x(1,2)+x(1,1)), 0.5d0*(y(1,2)+y(1,1)), 1)
        NormalWy(1,1)=-Normal(x(1,1),y(1,1),0.5d0*(x(1,2)+x(1,1)), 0.5d0*(y(1,2)+y(1,1)), 2)
        
        NormalEx(1,1)=-Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) & 
                +y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)),1); 
        NormalEy(1,1)=-Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)),2); 
        NormalNx(1,1)=Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) &
                +y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)),1); 
        NormalNy(1,1)=Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2) & 
                +y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)),2);
        
        TangentNx(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),1)
        TangentNy(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),2)
        TangentEx(1,1)=Tangent(x(1,1),y(1,1),x(2,1),y(2,1),1)
        TangentEy(1,1)=Tangent(x(1,1),y(1,1),x(2,1),y(2,1),2)
        TangentWx(1,1)=Tangent(x(1,1),y(1,1),0.5d0*(x(1,1)+x(1,2)),0.5d0*(y(1,1)+y(1,2)),1)
        TangentWy(1,1)=Tangent(x(1,1),y(1,1),0.5d0*(x(1,1)+x(1,2)),0.5d0*(y(1,1)+y(1,2)),2)
        TangentSx(1,1)=-Tangent(x(1,1),y(1,1),0.5d0*(x(1,1)+x(2,1)),0.5d0*(y(1,1)+y(2,1)),1)
        TangentSy(1,1)=-Tangent(x(1,1),y(1,1),0.5d0*(x(1,1)+x(2,1)),0.5d0*(y(1,1)+y(2,1)),2)
        
        CurviNx(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),1)
        CurviNy(1,1)=-Tangent(x(1,1),y(1,1),x(1,2),y(1,2),2)
        CurviSx(1,1)=CurviSx(1,2)
        CurviSy(1,1)=CurviSy(1,2)
        CurviEx(1,1)=-Tangent(x(1,1),y(1,1),x(2,1),y(2,1),1)
        CurviEy(1,1)=-Tangent(x(1,1),y(1,1),x(2,1),y(2,1),2)
        CurviWx(1,1)=CurviWx(2,1)
        CurviWy(1,1)=CurviWy(2,1)
        
        DistN(1,1)=sqrt((x(1,2)-x(1,1))**2+(y(1,2)-y(1,1))**2)
        DistS(1,1)=0d0
        DistE(1,1)=sqrt((x(2,1)-x(1,1))**2+(y(2,1)-y(1,1))**2)
        DistW(1,1)=0d0
        
        DistDualN(1,1)=CellAreaN(1,1)
        DistDualS(1,1)=CellAreaS(1,1)
        DistDualE(1,1)=CellAreaE(1,1)
        DistDualW(1,1)=CellAreaW(1,1)
      
      !North-West
        CellVol(1,N)=AreaElement(0.5d0*(x(1, N-1)+x(1, N)), & 
                                0.5d0*(y(1, N-1)+y(1, N)), &
                                0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                                0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                                0.5d0*(x(1, N)+x(2, N)), &
                                0.5d0*(y(1, N)+y(2, N)), &
                                x(1,N), &
                                y(1,N))
!         CellVol(1,N)=0.25d0*AreaElement(x(1,N-1),y(1,N-1),x(2,N-1),y(2,N-1),x(2,N),y(2,N),x(1,N),y(1,N))
        
!         CellVol(1,N)=2d0*CellVol(1,N)
        
        CellAreaN(1,N)=sqrt((x(1,N)-0.5d0*(x(1,N)+x(2,N)))**2+(y(1,N)-0.5d0*(y(1,N)+y(2,N)))**2)
        CellAreaW(1,N)=sqrt((x(1,N)-0.5d0*(x(1,N)+x(1,N-1)))**2+((y(1,N)-0.5d0*(y(1,N)+y(1,N-1))))**2)
        CellAreaS(1,N)=sqrt((0.25d0*(x(2,N-1)+x(1,N-1)+x(2,N)+x(1,N))-0.5d0*(x(1,N-1)+x(1,N)))**2 & 
                      +(0.25d0*(y(2,N-1)+y(1,N-1)+y(2,N)+y(1,N))-0.5d0*(y(1,N-1)+y(1,N)))**2) !sqrt((0.5d0*(x(1,N-1)+x(1,N))-0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)))**2+(0.5d0*(y(1,N-1)+y(1,N))-0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)))**2)
        CellAreaE(1,N)=sqrt((0.5d0*(x(2,N)+x(1,N))-0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)))**2 & 
                      +(0.5d0*(y(2,N)+y(1,N))-0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)))**2) !sqrt((0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1))-(0.5d0*(x(2,N)+x(1,N))))**2+(0.25d0*(y(1,N)+y(1,N-1)+y(2,N)+y(2,N-1))-0.5d0*(y(2,N)+y(1,N)))**2)
        
        NormalNx(1,N)=-Normal(x(1,N),y(1,N),0.5d0*(x(1,N)+x(2,N)),0.5d0*(y(1,N)+y(2,N)),1) 
        NormalNy(1,N)=-Normal(x(1,N),y(1,N),0.5d0*(x(1,N)+x(2,N)),0.5d0*(y(1,N)+y(2,N)),2)
        NormalWx(1,N)=Normal(x(1,N),y(1,N),0.5d0*(x(1,N)+x(1,N-1)),0.5d0*(y(1,N)+y(1,N-1)),1)
        NormalWy(1,N)=Normal(x(1,N),y(1,N),0.5d0*(x(1,N)+x(1,N-1)),0.5d0*(y(1,N)+y(1,N-1)),2)
        NormalEx(1,N)=-Normal(0.5d0*(x(1, N)+x(2, N)), &
                            0.5d0*(y(1, N)+y(2, N)), &
                            0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)) ,1); 
        NormalEy(1,N)=-Normal(0.5d0*(x(1, N)+x(2, N)), &
                            0.5d0*(y(1, N)+y(2, N)), &
                            0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), 2)
        NormalSx(1,N)=-Normal(0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                            0.5d0*(x(1, N-1)+x(1, N)), &
                            0.5d0*(y(1, N-1)+y(1, N)),1); 
        NormalSy(1,N)=-Normal(0.25d0*(x(1, N-1)+x(1, N)+x(2, N)+x(2, N-1)), &
                            0.25d0*(y(1, N-1)+y(1, N)+y(2, N)+y(2, N-1)), &
                            0.5d0*(x(1, N-1)+x(1, N)), &
                            0.5d0*(y(1, N-1)+y(1, N)), 2)
        TangentNx(1,N)=-Tangent(0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),x(1,N),y(1,N),1)
        TangentNy(1,N)=-Tangent(0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),x(1,N),y(1,N),2)
        TangentWx(1,N)=-Tangent(x(1,N),y(1,N),0.5d0*(x(1,N-1)+x(1,N)), 0.5d0*(y(1,N-1)+y(1,N)), 1)
        TangentWy(1,N)=-Tangent(x(1,N),y(1,N),0.5d0*(x(1,N-1)+x(1,N)), 0.5d0*(y(1,N-1)+y(1,N)), 2)
        TangentEx(1,N)=-Tangent(0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1)),0.25d0*(y(1,N)+y(1,N-1) & 
                        +y(2,N)+y(2,N-1)),0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),1)
        TangentEy(1,N)=-Tangent(0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1)),0.25d0*(y(1,N)+y(1,N-1) &
                        +y(2,N)+y(2,N-1)),0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),2)
        TangentSx(1,N)=-Tangent(0.5d0*(x(1,N-1)+x(1,N)),0.5d0*(y(1,N-1)+y(1,N)),0.25d0*(x(2,N-1) &
                        +x(2,N)+x(1,N-1)+x(1,N)),0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)),1)
        TangentSy(1,N)=-Tangent(0.5d0*(x(1,N-1)+x(1,N)),0.5d0*(y(1,N-1)+y(1,N)),0.25d0*(x(2,N-1) &
                        +x(2,N)+x(1,N-1)+x(1,N)),0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)),2)
        
        CurviNx(1,N)=CurviNx(1,N-1)
        CurviNy(1,N)=CurviNy(1,N-1)
        CurviSx(1,N)=-Tangent(x(1,N),y(1,N),x(1,N-1),y(1,N-1),1)
        CurviSy(1,N)=-Tangent(x(1,N),y(1,N),x(1,N-1),y(1,N-1),2)
        CurviEx(1,N)=-Tangent(x(1,N),y(1,N),x(2,N),y(2,N),1)
        CurviEy(1,N)=-Tangent(x(1,N),y(1,N),x(2,N),y(2,N),2)
        CurviWx(1,N)=CurviWx(2,N)
        CurviWy(1,N)=CurviWy(2,N)
        
         DistN(1,N)=0d0 !sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
         DistS(1,N)=sqrt((x(1,N)-x(1,N-1))**2+(y(1,N)-y(1,N-1))**2)
         DistE(1,N)=sqrt((x(2,N)-x(1,N))**2+(y(2,N)-y(1,N))**2)
         DistW(1,N)=0d0 !sqrt((x(i,N)-x(i-1,N))**2+(y(i,N)-y(i-1,N))**2)
         
         DistDualN(1,N)=CellAreaN(1,N)
         DistDualS(1,N)=CellAreaS(1,N)
         DistDualE(1,N)=CellAreaE(1,N)
         DistDualW(1,N)=CellAreaW(1,N)
      
      
   write(*,*) "VESSEL CHECK"
   write(*,*) "North", NormalNxP(Mp/2,Np-1), NormalNyP(Mp/2,Np-1)
   write(*,*) "South", NormalSxP(Mp/2,2), NormalSyP(Mp/2,2)
   write(*,*) "East", NormalExP(Mp-1,Np-1), NormalEyP(Mp-1,Np-1)
   write(*,*) "West", NormalWxP(Mp-1,Np-1), NormalWyP(Mp-1,Np-1)
   write(*,*) CellAreaNP(Mp/2,Np/2), CellAreaSP(Mp/2,Np/2), CellAreaEP(Mp/2,Np/2), CellAreaWP(Mp/2,Np/2)
   
   ! Drift initialization
   JeX(:,:)=0d0
   JeY(:,:)=0d0
   JhX(:,:)=0d0
   JhY(:,:)=0d0
  
   
   !! defining material index and ionization constants
   epsilonInf=DielectricConstant(lambda)
   OnePhotonIonizationRate0=OnePhotonIonizationRate(lambda,epsilonInf)
   TwoPhotonIonizationRate0=TwoPhotonIonizationRate(lambda)
   
  write(*,*) 'epsilon(', 1d9*lambda, 'nm)=', epsilonInf
  write(*,*) 'Re(sqrt(epsilon))=', real(sqrt(epsilonInf))
if(UseMieScattering.eq.1) then
  write(*,*) 'Computing the Mie scattering field distribution...'
  write(*,*) 'Angle Mie =', phiMie0
  write(*,*) 'Polarization TM ? ', PolarizationSource
  ! $ O M P DO
  do j=1,N
    do i=1,M
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
    ! $ OMP END DO

    

    ! $ OMP DO
   do j=1,N
    do i=1,M
          if(PolarizationSource.eq.1) then !TM polarization, Bassel et al scattering on a cylinder
          ! formula for an experimental needle with interpolated radius
!             write(*,*) "TM polarization selected."
            EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf) ! * sqrt(2d0*fluence/(c*epsilon0*tau))
            EintField2(i,j)=Zero
          ! formula with a super mistake on radius
!           EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), 0.5d0*(y(i,N)-y(i,1)), epsilonInf) ! * sqrt(2d0*fluence/(c*epsilon0*tau))

          ! formulas for an hyperbolic needle
!           EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), abs(ContourYofX(x(i,j), NeedleRadius, NeedleAngle)), epsilonInf) ! * sqrt(2d0*fluence/(c*epsilon0*tau))

          ! formula for debug, using a constant radius
!         EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie(i,j), 100d-9, epsilonInf) ! * sqrt(2d0*fluence/(c*epsilon0*tau)) !with a constant radius

          else !TE polarization
!             write(*,*) "TE polarization selected."
            EintField2(i,j)=Unit * MieScatteringTE2(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf, besselArray)
            EintField(i,j)=Unit * MieScatteringTE1(abs(y(i,j)), phiMie(i,j), 1d-6*Radius(i,j), epsilonInf, besselArray)
          end if
      end do
    end do
    ! $ OMP END DO

!     EintFieldR=sqrt(EintField * conjg(EintField))
    EintFieldR=real(sqrt( EintField * conjg(EintField) + EintField2 * conjg(EintField2) ))
    

    
    do j=1,N
      do i=1,M
        write(104, 891, advance='yes') x(i,j), y(i,j), (EintFieldR(i,j)**2d0)**0.5d0, Radius(i,j)
891        FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
      end do
    end do


  write(*,*) 'Done.'
 end if
 
  t=tmin
  
  call cpu_time(calc_time_2)
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
 
  do i=1,M
    CellAreaN(i,N)=0d0 
!     CellAreaN(i,N-1)=0d0
!     CellAreaS(i,2)=0d0
    CellAreaS(i,1)=0d0
    
    DistN(i,N-1)=sqrt((0.5d0*(x(i,N-1)+x(i,N))-x(i,N-1))**2+(0.5d0*(y(i,N-1)+y(i,N))-y(i,N-1))**2)
    DistS(i,2)=sqrt((x(i,2)-0.5d0*(x(i,2)+x(i,1)))**2+(y(i,2)-0.5d0*(y(i,2)+y(i,1)))**2) 
  end do
   
  do j=1,N
    CellAreaW(1,j)=0d0
!     CellAreaW(2,j)=0d0
!     CellAreaE(M-1,j)=0d0
    CellAreaE(M,j)=0d0
    
    DistW(2,j)=sqrt((x(2,j)-0.5d0*(x(2,j)+x(1,j)))**2+((y(2,j)-0.5d0*(y(2,j)+y(1,j))))**2)
    DistE(M-1,j)=sqrt((0.5d0*(x(M-1,j)+x(M,j))-x(M-1,j))**2+(0.5d0*(y(M-1,j)+y(M,j))-y(M-1,j))**2)
  end do
  
     
    do j=2, N-1
      do i=2, M-1
      ! define the volume of elementary cell around a point everywhere but not on boundaries
!         CellVol(i,j)=0.25d0*(AreaElement(x(i-1,j-1),y(i-1,j-1),x(i+1,j-1),y(i+1,j-1),x(i+1,j+1),y(i+1,j+1),x(i-1,j+1),y(i-1,j+1)))
        CellVol(i,j)=AreaElement(0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)), & !x(i-1/2,j-1/2)
                                 0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j)), &
                                 0.25d0*(x(i,j-1)+x(i+1,j-1)+x(i+1,j)+x(i,j)), &        !x(i+1/2,j-1/2)
                                 0.25d0*(y(i,j-1)+y(i+1,j-1)+y(i+1,j)+y(i,j)), &
                                 0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), &        !x(i+1/2,j+1/2)
                                 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
                                 0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)), &        !x(i-1/2,j+1/2)
                                 0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1)))
                                 
        InvCellVol(i,j) = 1.0d0/CellVol(i,j)
                                 
        ! Calcul de Grad(Ne) sur le maillage direct
        ! Première estimation peu stable
        GradNeX(i,j) = 0.5d0 * InvCellVol(i,j) * &
                      ( (mesh%Ne(i,j) + mesh%Ne(i,j+1)) * CellAreaN(i,j) * NormalNx(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i,j-1)) * CellAreaS(i,j) * NormalSx(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i-1,j)) * CellAreaW(i,j) * NormalWx(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i+1,j)) * CellAreaE(i,j) * NormalEx(i,j) )
        GradNeY(i,j) = 0.5d0 * InvCellVol(i,j) * &
                      ( (mesh%Ne(i,j) + mesh%Ne(i,j+1)) * CellAreaN(i,j) * NormalNy(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i,j-1)) * CellAreaS(i,j) * NormalSy(i,j) &
                      + (mesh%Ne(i,j) + mesh%Ne(i-1,j)) * CellAreaW(i,j) * NormalWy(i,j) &
                      * (mesh%Ne(i,j) + mesh%Ne(i+1,j)) * CellAreaE(i,j) * NormalEy(i,j) )
                      
        ! interpolation lineaire des valeurs de phi sur les bords de cellules
!         phiN=Ne(i,j)+GradNeX(i,j)*(0.5d0*(xDual(i,j)+xDual(i-1,j)))+GradNeY(i,j)*(0.5d0*(yDual(i,j)+yDual(i-1,j)))
!         phiS=Ne(i,j)+GradNeX(i,j)*(0.5d0*(xDual(i,j-1)+xDual(i-1,j-1)))+GradNeY(i,j)*(0.5d0*(yDual(i,j-1)+yDual(i-1,j-1)))
        ! recalcul du gradient avec phi_bords
        
        end do
      end do
  
  ! interpolation and preparation of resolution
  do j=2,N-1
    do i=2,M-1
        !Dual mesh calculation
        xDualSW(i,j) = 0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)) !x(i-1/2,j-1/2)
        yDualSW(i,j) = 0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j))
        xDualSE(i,j) = 0.25d0*(x(i,j-1)+x(i+1,j-1)+x(i+1,j)+x(i,j)) !x(i+1/2,j-1/2)
        yDualSE(i,j) = 0.25d0*(y(i,j-1)+y(i+1,j-1)+y(i+1,j)+y(i,j))
        xDualNE(i,j) = 0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)) !x(i+1/2,j+1/2)
        yDualNE(i,j) = 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))
        xDualNW(i,j) = 0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)) !x(i-1/2,j+1/2)
        yDualNW(i,j) = 0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1))
                                
        xDual(i-1,j-1) = xDualSW(i,j); yDual(i-1,j-1)=yDualSW(i,j) ! to make for (i,j) from 2,2 to M,N
        xDual(i-1,N-1) = xDualNW(i,N-1); yDual(i-1,N-1) = yDualNW(i,N-1); 
        xDual(M-1,j-1) = xDualSE(M-1,j); yDual(M-1,j-1) = yDualSE(M-1,j)
        xDual(M-1,N-1) = xDualNE(M-1,N-1); yDual(M-1,N-1) = yDualNE(M-1,N-1); 
    end do
  end do
  
  do j=1,N
    do i=1,M
       NormalN2(i,j) = (NormalNx(i,j)**2+NormalNy(i,j)**2)
       NormalE2(i,j) = (NormalEx(i,j)**2+NormalEy(i,j)**2)
       NormalW2(i,j) = (NormalWx(i,j)**2+NormalWy(i,j)**2)
       NormalS2(i,j) = (NormalSx(i,j)**2+NormalSy(i,j)**2)
    end do
  end do
  
  do j=2, N-1
    do i=2, M-1
      NormalN2(i,j) = NormalN2(i,j)/DistN(i,j)
      NormalW2(i,j) = NormalW2(i,j)/DistW(i,j)
      NormalE2(i,j) = NormalE2(i,j)/DistE(i,j)
      NormalS2(i,j) = NormalS2(i,j)/DistS(i,j)
      ShapeFactorNormalE(i,j)  = CellAreaE(i,j) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
      ShapeFactorNormalW(i,j)  = CellAreaW(i,j) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
      ShapeFactorNormalN(i,j)  = CellAreaN(i,j) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
      ShapeFactorNormalS(i,j)  = CellAreaS(i,j) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
      ShapeFactorTangentE(i,j) = CrossCoeff * ( CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j) ) / DistDualE(i,j)
      ShapeFactorTangentW(i,j) = CrossCoeff * ( CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j) ) / DistDualW(i,j)
      ShapeFactorTangentN(i,j) = CrossCoeff * ( CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j) ) / DistDualN(i,j)
      ShapeFactorTangentS(i,j) = CrossCoeff * ( CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j) ) / DistDualS(i,j)
    end do
  end do
  
  ! Steps for 3rd order integration
  
  h1=dt
  h2=dt+dt2
  h3=dt+dt2+dt3
  
  write(*,*) "Starting time loop."

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
                   "maxCFL_Ts=",maxCFLxTs+maxCFLyTs, "CPU=", cpuefficiency, "NumThreads=", nthreads, "Elapsed time=", cpu_timestep_duration
      write(*,*) "CFL_Limit=", maxCFL
      write(*,*) "dt_init=", dt0, "dt=", dt
    end if
    
    maxIntensity=0d0; maxTe=0d0; minTe=1d10; maxTh=0d0; minTh=1d10; maxTs=0d0; minTs=1d10; 
    maxNe=0d0; minNe=1d50; maxNh=0d0; minNh=1d50; maxCFLxT=0d0; maxCFLyT=0d0; maxCFLxN=0d0; maxCFLyN=0d0; maxCFLxTs=0d0; 
    maxCFLyTs=0d0; maxSourceE=0d0; maxSourceH=0d0; maxGainsE=0d0; maxGainsH=0d0
    
   !$OMP PARALLEL DEFAULT (PRIVATE) SHARED (dt, dt1, dt2, dt3, dt4, UeNew, UhNew, TsOld, TsPrev, &
   !$OMP& mesh, dual, intensityDual, &
   !$OMP& Ue, Uh, GradNeX, GradNeY, intensity, intensity2, reflectivity, FermiTableE, FermiTableH, &
   !$OMP& Dielectric, DielectricDrudeE, DielectricDrudeH, absorptionDrudeE, absorptionDrudeH, &
   !$OMP& x, y, xDual, yDual, xDualSW, yDualSW, xDualSE, yDualSE, xDualNE, yDualNE, xDualNW, yDualNW, &
   !$OMP& diffusionE, diffusionH, GainsE, GainsH, LossesE, LossesH, &
   !$OMP& kappae, kappah, kappas, Ce, CeOld, Ch, ChOld, Cs, CsOld, CsPrev, CsPrev2, CouplingE, CouplingH, &
   !$OMP& nuColl, nuColleph, mobilityE, mobilityH, etae, etah, Egap, & 
   !$OMP& SourceE, SourceH, SourceUe, SourceUh, diffNe, diffNh, CFLxT, CFLyT, CFLxN, CFLyN, CFLxTs, CFLyTs, &
   !$OMP& ThermalEnergy, LaserEnergy, epsilonInf, FermiIndexE, FermiIndexH, FermiRatioE, FermiRatioH, &
   !$OMP& OmegaX, OmegaY, JeX, JeY, JhX, JhY, VeX, VeY, VhX, VhY, DielectricStatic, Amatrix, Xvector, XvectorPrev, Bvector, xV, yV, xP, yP, &
   !$OMP& spectralNorm, Ex, Ey, ExPoisson, EyPoisson, potential, potentialNeedle, NeP, NhP, FixedPotentialIndex, &
   !$OMP& NormalN2,NormalNx, NormalNy, NormalS2, NormalSx, NormalSy, NormalE2, NormalEx, NormalEy, NormalW2, NormalWx, NormalWy, &
   !$OMP& CellVolume, CellAreaN, CellAreaS, CellAreaE, CellAreaW, CellVol, InvCellVol, CellAreaNP, CellAreaSP, CellAreaEP, CellAreaWP, &
   !$OMP& NormalNxP, NormalNyP, NormalSxP, NormalSyP, NormalExP, NormalEyP, NormalWxP, NormalWyP, TangentWx, TangentWy, TangentNx, &
   !$OMP& TangentNy, TangentSx, TangentSy, TangentEx, TangentEy, CurviNx, CurviNy, CurviSx, CurviSy, CurviEx, CurviEy, CurviWx, &
   !$OMP& CurviWy, ConstBLx, ConstBLy, DistN, DistS, DistE, DistW, DistDualN, DistDualS, DistDualE, DistDualW, &
   !$OMP& EintField, EintFieldDual, EintFieldI, EintFieldR, NeTotal, NhTotal, &
   !$OMP& ShapeFactorNormalE, ShapeFactorNormalN, ShapeFactorNormalS, ShapeFactorNormalW, &
   !$OMP& ShapeFactorTangentE, ShapeFactorTangentN, ShapeFactorTangentS, ShapeFactorTangentW ), &
   !$OMP& FIRSTPRIVATE (t, t0, x0, y0, I0, I1, I2, I3, I4, I5, I6, I7, &
   !$OMP& I8, I9, &
   !$OMP& x1, x2, x3, x4, x5, x6, x7, x8, x9, &
   !$OMP& y1, y2, y3, y4, y5, y6, y7, y8, y9, &
   !$OMP& h1, h2, h3, h4, OnePhotonIonizationRate0, TwoPhotonIonizationRate0, &
   !$OMP& sigmaX1, sigmaX2, sigmaX3, sigmaX4, sigmaX5, sigmaX6, sigmaX7, &
   !$OMP& sigmaX8, sigmaX9, &
   !$OMP& sigmaY1, sigmaY2, sigmaY3, sigmaY4, sigmaY5, sigmaY6, sigmaY7, &
   !$OMP& sigmaY8, sigmaY9, &
   !$OMP& AugerRateE, AugerRateH, sigmaTau, sigmaX, sigmaY, dx, dy, nbiter, &
   !$OMP& cpuefficiency, ColFermi0, ColFermi1, ColFermi2, ColFermiEta, ColFermiHalf, &
   !$OMP& ColFermiMenusHalf, ColFermiNeNc, ColFermiThreeHalf, SORsum, Mp, Np)

   nthreads = OMP_GET_NUM_THREADS()

   
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

    
!!!! thermal calculations in the main domain
! calculation of sources
    !$OMP DO  COLLAPSE(2) 
    do j=1,N
        do i=1,M
!      intensity(i,1)=(1d0-reflectivity(i,1))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2.-.5d0*((x(i,1)-x0)/sigmaX)**2.-.5d0*((y(i,1)-y0)/sigmaY)**2.)

        ! optical coefficients
        
!         ! DEBUG


!         Ex(i,j)=0d0 !-1d9
!         Ey(i,j)=0d0*        -1d9
        
        nuColl(i,j)=CollisionFrequency()
        nuColleph(i,j)=ephCollisionFrequency(mesh%Ne(i,j))
        Dielectric(i,j)=DielectricFunction(lambda, epsilonInf, mesh%Ne(i,j), nuColl(i,j))
        DielectricDrudeE(i,j)=DielectricFunctionDrude(mesh%Ne(i,j), nuColl(i,j),me)
        DielectricDrudeH(i,j)=DielectricFunctionDrude(mesh%Nh(i,j), nuColl(i,j),mh)

        !         write(*,*) "Esprit es-tu la ?"
        
!         write(*,*) i,j,TeNew(i,j)
        DOSe(i,j)=DensityOfStateE(mesh%Te(i,j))
        DOSh(i,j)=DensityOfStateH(mesh%Th(i,j))
        FermiRatioE(i,j)=mesh%Ne(i,j)/DOSe(i,j)
        FermiRatioH(i,j)=mesh%Nh(i,j)/DOSh(i,j)
        FermiIndexE(i,j)=1! FermiIndex(FermiRatioE(i,j)) !1
        FermiIndexH(i,j)=1! FermiIndex(FermiRatioH(i,j)) !1
!         write(*,*) "iter=", nbiter, "DOS=", DOSe(i,j), DOSh(i,j)
        etae(i,j)=FermiTableE(ColFermiEta,FermiIndexE(i,j)) 
        etah(i,j)=FermiTableH(ColFermiEta,FermiIndexH(i,j))
!         write(*,*) "iter=", nbiter, "NeNc=", Ne(i,j)/DOSe(i,j), Nh(i,j)/DOSh(i,j)
!         write(*,*) "iter=", nbiter, "FermiIndex=", FermiIndex(Ne(i,j)/DOSe(i,j)), FermiIndex(Nh(i,j)/DOSh(i,j))
!         write(*,*) "iter=", nbiter, "eta=", etae(i,j), etah(i,j)
!         write(*,*) "FermiTables: etaE,etaH=", FermiTableE(3,463), FermiTableH(3,450)

        mobilityE(i,j)=(ec/(me*nuColl(i,j)))*FermiTableE(ColFermi0, FermiIndexE(i,j))/FermiTableE(ColFermiHalf, FermiIndexE(i,j))
        mobilityH(i,j)=(ec/(mh*nuColl(i,j)))*FermiTableH(ColFermi0, FermiIndexH(i,j))/FermiTableH(ColFermiHalf, FermiIndexH(i,j))
        
!         write(*,*) "iter=", nbiter, "mobility=", mobilityE(i,j), mobilityH(i,j)
        if(DrudeHeating==1) then
          absorptionDrudeE(i,j)=4d0*pi/lambda*aimag(sqrt(DielectricDrudeE(i,j)))
          absorptionDrudeH(i,j)=4d0*pi/lambda*aimag(sqrt(DielectricDrudeH(i,j)))
        else
          absorptionDrudeE(i,j)=0d0
          absorptionDrudeH(i,j)=0d0
!           absorptionDrudeE(i,j)=sqrt(2d0)*sqrt(mu0*omegaLaser*ec*mobilityE(i,j)*Ne(i,j))
!           absorptionDrudeH(i,j)=sqrt(2d0)*sqrt(mu0*omegaLaser*ec*mobilityH(i,j)*Nh(i,j))
        end if
        
!         write(*,*) "iter=", nbiter, "absorption=", absorptionDrudeE(i,j), absorptionDrudeH(i,j)
        reflectivity(i,j)=(real(sqrt(Dielectric(i,j)))**2+aimag(sqrt(Dielectric(i,j)))**2-2d0*real(sqrt(Dielectric(i,j)))+1d0) &
                            /(real(sqrt(Dielectric(i,j)))**2+aimag(sqrt(Dielectric(i,j)))**2+2d0*real(sqrt(Dielectric(i,j))+1d0))
        !local intensity
        if(UseMieScattering .eq. -1) then
  !         ! DEBUG ZONE
  ! !         if(lambda.eq.343d-9) then
  !         ! uniform distribution like in Elena's paper
          intensity(i,j)=(1d0-0e0*reflectivity(i,j))*real(sqrt(Dielectric(i,j)))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2)
  !         intensity(i,j)=I0*exp(-.5d0*((t-t0)/sigmaTau)**2)*exp(-0.5d0*(((y(i,j)-500d-9)/sigmaY)**2+(x(i,j)/sigmaX)**2))
  !         ! with just nothing
  ! !           intensity(i,j)=(1d0-reflectivity(i,j))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2)*exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2)
  !         ! with beer-lambert
  ! !         intensity(i,j)=(-(absorptionDrude(i,j) &
  ! !                         +4d0*pi/lambda*aimag(sqrt(epsilonInf)))*intensity(i,j-1) &
  ! !                         -TwoPhotonIonizationRate0*intensity(i,j-1)**2 &
  ! !                         )*x &
  ! !                         +intensity(i,j-1)
  ! !          end if
        
        else if(UseMieScattering .eq. 0) then
          ! WITH EXTERNALLY ADJUSTED INPUTS
  !        !Lumerical mode already contains the reflectivity. Although, it doesn't consider change of optical index with ionization. 
          if(lambda.eq.1030d-9) then 
            ConstBLx=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate0+1d0*TwoPhotonIonizationRate0) &
                      / (1d0*exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*x0) &
                      * (OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
            ConstBLy=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate0+1d0*TwoPhotonIonizationRate0) &
                      / (1d0*exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*y0) &
                      * (OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
  ! case 1030 nm distribution
            !initial field distribution  
            if(BandBendingInFDTD.eq.0) then
              intensity(i,j)=(1d0-0d0*reflectivity(i,N))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2) & 
                      *exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2)
            else 
            intensity(i,j)=(1d0-0d0*reflectivity(i,N))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2) & 
                            *(exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
                            + I0*1d-4*( & 
                            exp(-0.5d0*(((x(i,j)-x(i,N))**2+(y(i,j)-y(i,N))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
                            +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
                            +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
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
          if(lambda.eq.515d-9) then 
            ConstBLx=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate0+1d0*TwoPhotonIonizationRate0) &
                      / (1d0*exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*x0) &
                      * (OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
            ConstBLy=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate0+1d0*TwoPhotonIonizationRate0) &
                      / (1d0*exp(-(OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*y0) &
                      * (OnePhotonIonizationRate0+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
            intensity(i,j)=(1d0-0d0*reflectivity(i,N))*exp(-.5d0*((t-t0)/sigmaTau)**2) &
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
          if(lambda.eq.343d-9) then
            intensity(i,j)= (1d0-0e0*reflectivity(i,N))* & 
                            I0*exp(-.5d0*((t-t0)/sigmaTau)**2) & 
                            *( &
                            exp(-(OnePhotonIonizationRate(lambda, epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) & 
                            *abs(y(i,j)-y(i,N)) & !introduce discontinuity !
                            ) & 
                            * exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
  !                           + exp(-(OnePhotonIonizationRate(lambda, epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*abs(y(i,j)-y(i,1))) & 
  !                           *exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
                            )
  !                           *exp(-(OnePhotonIonizationRate(lambda, epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*abs(x(i,j)-x(i,N)))
          end if

        ! USING MIE SCATTERING ANALYTICAL FORMULAS
        else if(UseMieScattering .eq. 1) then
        ! calculate electric field inside the tip
!           EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie, abs(ContourYofX(x(i,j), NeedleRadius, NeedleAngle)), Dielectric(i,j)) !*sqrt(2d0*fluence/(c*epsilon0*tau))
          ! debug formula for constant cone radius
!           EintField(i,j)=MieScattering(abs(y(i,j)), phiMie, 100d-9, epsilonInf)
!           EintField(i,j)=sqrt(EintField(i,j)*conjg(EintField(i,j))) !complex to real
          intensity(i,j)=I0 * real(sqrt(Dielectric(i,j))) * EintFieldR(i,j)**2 * exp(-.5d0*((t-t0)/sigmaTau)**2) !laser fluence and reflectivity is inside the field
        else 
          write(*,*) "Input ERROR. Check the MieScattering parameter."
          stop
        end if
        
        if(intensity(i,j) < 1d-20) then 
          intensity(i,j)=0d0
        end if
        

        ! free-carrier balance sources
        Egap(i,j)=EgapValue(mesh%Ne(i,j),mesh%Ts(i,j))

        diffusionE(i,j)=mobilityE(i,j)*kb*mesh%Te(i,j)/ec &
                *FermiTableE(ColFermiHalf,FermiIndexE(i,j))/FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j))
        diffusionH(i,j)=mobilityH(i,j)*kb*mesh%Th(i,j)/ec &
                *FermiTableH(ColFermiHalf,FermiIndexH(i,j))/FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j))
        if(ConductivityFix.eq.-1) then
           diffusionE(i,j)=0d0
           diffusionH(i,j)=0d0
        end if        

        JeX(i,j)=-mobilityE(i,j)*mesh%Ne(i,j)*Ex(i,j)
        JeY(i,j)=-mobilityE(i,j)*mesh%Ne(i,j)*Ey(i,j)
        JhX(i,j)=mobilityH(i,j)*mesh%Nh(i,j)*Ex(i,j)
        JhY(i,j)=mobilityH(i,j)*mesh%Nh(i,j)*Ey(i,j)

        if(DriftOn.eq.0) then 
          JeX(i,j)=0d0; JeY(i,j)=0d0; 
          JhX(i,j)=0d0; JhY(i,j)=0d0;
        endif

        Int2 = intensity(i,j)**2
        GainsE(i,j)=(OnePhotonIonizationRate0*intensity(i,j)/hbar/omegaLaser &
                    +TwoPhotonIonizationRate0*Int2/(2d0*hbar*omegaLaser) &
                    +ImpactIonizationRate(mesh%Te(i,j),mesh%Ne(i,j),mesh%Ts(i,j))*mesh%Ne(i,j))! *(4d0*SiDensity-Ne(i,j))/(4d0*SiDensity) !use Old Ne here!
                    
        GainsH(i,j)=(OnePhotonIonizationRate0*intensity(i,j)/hbar/omegaLaser &
                    +TwoPhotonIonizationRate0*Int2/(2d0*hbar*omegaLaser) &
                    +ImpactIonizationRate(mesh%Te(i,j),mesh%Ne(i,j),mesh%Ts(i,j))*mesh%Nh(i,j)) !*(4d0*SiDensity-Ne(i,j))/(4d0*SiDensity) !use Old Nh here
                    
        LossesE(i,j)=AugerRateE * (mesh%Ne(i,j))**2d0 * mesh%Nh(i,j) + AugerRateH * (mesh%Nh(i,j))**2d0 * mesh%Ne(i,j) !use Old Ne, Nh here!
        LossesH(i,j)=LossesE(i,j)

        
        ! thermal coefficients
        kappae(i,j)=kb**2*mesh%Ne(i,j)*mobilityE(i,j)*mesh%Te(i,j)/ec*(6d0*FermiTableE(ColFermi2,FermiIndexE(i,j)) &
                /FermiTableE(ColFermi0,FermiIndexE(i,j)) &
                -4d0*(FermiTableE(ColFermi1,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)))**2)
        kappah(i,j)=kb**2*mesh%Nh(i,j)*mobilityH(i,j)*mesh%Th(i,j)/ec*(6d0*FermiTableH(ColFermi2,FermiIndexH(i,j)) &
                /FermiTableH(ColFermi0,FermiIndexH(i,j)) &
                -4d0*(FermiTableH(ColFermi1,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)))**2)
!         kappas(i,j)=-.1412d0*Ts(i,j)**(1.38961d0)+0.638157d0*Ts(i,j)**(1.14013d0) !mingo till 300 K
        kappas(i,j)=(-8.992d0+68.265d0/(1d0+exp(-.4075612391d-1*mesh%Ts(i,j)+2.315984470d0))*(1d0-1d0/ &
                    (1d0+exp(-.4756634637d-2*mesh%Ts(i,j)+2.403533689d0)))) !Elena fit sur Kazan (2010)
        if(kappas(i,j).lt.0d0) then
          kappas(i,j)=0d0
        end if
        
!        ! correction considering Fick diffusion in energy
!         if(ConductivityFix.eq.1) then 
!             kappae(i,j)=kappae(i,j) + kb**2*Te(i,j)*Ne(i,j)*mobilityE(i,j) / ec &
!                       * (etae(i,j) - 2d0*FermiTableE(ColFermi1,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) )**2 
!             kappah(i,j)=kappah(i,j) + kb**2*Th(i,j)*Nh(i,j)*mobilityH(i,j) / ec &
!                       * (etah(i,j) - 2d0*FermiTableH(ColFermi1,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) )**2 
!         else if(ConductivityFix.eq.2) then 
!             kappae(i,j)=kappae(i,j) + 2d0*kb**2*Te(i,j)*FermiTableE(ColFermi1,FermiIndexE(i,j))*mobilityE(i,j)*FermiTableE(ColFermiHalf, FermiIndexE(i,j))*Ne(i,j) * &
!                         (2d0*FermiTableE(ColFermi1, FermiIndexE(i,j))*FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)) & 
!                         /FermiTableE(ColFermiHalf,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) - 1.5d0) * & 
!                         (FermiTableE(ColFermi0, FermiIndexE(i,j))*ec*FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))**(-1e0)
!                         
!             kappah(i,j)=kappah(i,j) + 2d0*kb**2*Te(i,j)*FermiTableH(ColFermi1,FermiIndexH(i,j))*mobilityH(i,j)*FermiTableH(ColFermiHalf, FermiIndexH(i,j))*Ne(i,j) * &
!                         (2d0*FermiTableH(ColFermi1, FermiIndexH(i,j))*FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)) & 
!                         /FermiTableH(ColFermiHalf,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) - 1.5d0) * & 
!                         (FermiTableH(ColFermi0, FermiIndexH(i,j))*ec*FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)))**(-1.)
!                         
!         else if(ConductivityFix.eq.-1) then

        if(ConductivityFix.eq.-1) then
           kappae(i,j)=0d0
           kappah(i,j)=0d0
           kappas(i,j)=0d0
        end if
        
        Ce(i,j)=1.5d0*mesh%Ne(i,j)*kb*(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)) &
                  -etae(i,j)*(1d0-(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))* &
                  (FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))))
        Ch(i,j)=1.5d0*mesh%Nh(i,j)*kb*(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)) &
                  -etah(i,j)*(1d0-(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)))* &
                  (FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)))))
        Cs(i,j)=LatticeHeatCapacity(mesh%Ts(i,j))

!         Ce(i,j)=Ch(i,j) ! DEBUG test
        
        CouplingE(i,j)=Ce(i,j)*nuColleph(i,j)*(mesh%Te(i,j)-mesh%Ts(i,j))
        CouplingH(i,j)=Ch(i,j)*nuColleph(i,j)*(mesh%Th(i,j)-mesh%Ts(i,j))

        if(CouplingDebug.eq.1) then
          CouplingE(i,j)=0d0
          CouplingH(i,j)=0d0
        end if

        if(HolesOff.eq.1) then
          CouplingH(i,j)=0d0
        end if
! 
!         diffNe(i,j)=0d0 !just for debug !
!         diffNh(i,j)=0d0 !just for debug !!
!         
        SourceE(i,j)= (hbar*omegaLaser-Egap(i,j))/hbar/omegaLaser*((me)/(me+mh))*OnePhotonIonizationRate0*intensity(i,j) &
                     + (2d0*hbar*omegaLaser - Egap(i,j))/(2d0*hbar*omegaLaser)* ((me)/(me+mh)) * TwoPhotonIonizationRate0*Int2 &
                     - Egap(i,j)*ImpactIonizationRate(mesh%Te(i,j),mesh%Ne(i,j),mesh%Ts(i,j))*mesh%Ne(i,j) &
                     + absorptionDrudeE(i,j)*intensity(i,j) &
                     + Egap(i,j)*(AugerRateE*mesh%Nh(i,j) * mesh%Ne(i,j)**2d0)

                     
        SourceUe(i,j) = SourceE(i,j)
!         SourceE(i,j) = SourceE(i,j) - diffNe(i,j)*(1.5d0*kb*Te(i,j))*(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))
        SourceE(i,j) = SourceE(i,j) - mesh%Te(i,j) * (Ce(i,j)-CeOld(i,j))/dt

        SourceH(i,j)=(hbar*omegaLaser-Egap(i,j))/hbar/omegaLaser * ((me)/(me+mh)) * OnePhotonIonizationRate0*intensity(i,j) &
                     + (2d0*hbar*omegaLaser - Egap(i,j))/(2d0*hbar*omegaLaser)* ((me)/(me+mh)) *TwoPhotonIonizationRate0*Int2 &
                     - Egap(i,j)*ImpactIonizationRate(mesh%Th(i,j),mesh%Nh(i,j),mesh%Ts(i,j))*mesh%Nh(i,j) &
                     + absorptionDrudeH(i,j)*intensity(i,j) &
                     + Egap(i,j)*(AugerRateH*mesh%Ne(i,j) * mesh%Nh(i,j)**2d0)
                     
        SourceUh(i,j) = SourceH(i,j)
!         SourceH(i,j) = SourceH(i,j) - diffNh(i,j)*(1.5d0*kb*Th(i,j)*(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j))))
        SourceH(i,j) = SourceH(i,j) - mesh%Th(i,j) * (Ch(i,j)-ChOld(i,j))/dt

        ! first order precision
!         SourceS(i,j) = 0d0 - Ts(i,j)*(Cs(i,j)-CsOld(i,j))/Cs(i,j)

        ! second order precisino
!         SourceS(i,j) = 0d0 - (3d0*Cs(i,j)-4d0*CsOld(i,j)+CsPrev(i,j))/(2d0*dt) * Ts(i,j)

        ! third order precision in already included in the scheme (Maple generated since complexity increases substancially)
        
        VeX(i,j)=0d0
        VeY(i,j)=0d0
        VhX(i,j)=0d0
        VhY(i,j)=0d0
        
      end do
    end do
    !$OMP END DO

    
    ! interpolation on dual mesh
    ! InterpolateBiCubic(phi_source, x_s, y_s, x_t, y_t, SizeXs, SizeYs, SizeXt, SizeYt, phi_target, Grad(phi)_targetX, Grad(phi)_targetY)
    
!     if(UseInterpolation.eq.1) then
!     !$OMP SECTIONS
!       !$OMP SECTION
!       call InterpolateBiCubic(Ne, x, y, xDual, yDual, M, N, M-1, N-1, NeDual, DummyDual, DummyDual)
!       !$OMP SECTION
!       call InterpolateBiCubic(Nh, x, y, xDual, yDual, M, N, M-1, N-1, NhDual, DummyDual, DummyDual)
!       !$OMP SECTION
!       call InterpolateBiCubic(Te, x, y, xDual, yDual, M, N, M-1, N-1, TeDual, DummyDual, DummyDual)
!       !$OMP SECTION
!       call InterpolateBiCubic(Th, x, y, xDual, yDual, M, N, M-1, N-1, ThDual, DummyDual, DummyDual)
!       !$OMP SECTION
!       call InterpolateBiCubic(Ts, x, y, xDual, yDual, M, N, M-1, N-1, TsDual, DummyDual, DummyDual)
!       !$OMP SECTION
!   !     call InterpolateBiCubic(intensity, x, y, xDual, yDual, M, N, M-1, N-1, intensityDual, DummyDual, DummyDual)
!   !     intensityDual=InterpolateSelner(intensity, x, y, xDual, yDual, M, N, M-1, N-1)
!         
!   ! !       DEBUG: test de la fonction d'interpolation
!       call InterpolateBiCubic(intensity, x, y, xDual, yDual, M, N, M-1, N-1, intensityDual, DummyDual, DummyDual)
!   !     intensityDual=InterpolateSelner(intensity, x, y, xDual, yDual, M, N, M-1, N-1)
! 
!     !$OMP END SECTIONS
!     end if

      ! interpolation bilineaire ponderee par les aires
      call bilinear_interpol_dual(mesh, dual, InvCellVol)


      ! solving the 2D problem
      
      !$OMP DO COLLAPSE(2) !(optimized)
!       do i=2, M-1
        do j=2, N-1 !(optimized)
!         do j=2, N-1
        do i=2, M-1 !(optimized)
        
     
!          ! diffusion is separated from drift
         newmesh%Ne(i,j) = mesh%Ne(i,j) + dt*InvCellVol(i,j)*(&
              ( GainsE(i,j)-LossesE(i,j) )*CellVol(i,j)& 
                  !         
                  +0.5d0*( &
                  + ShapeFactorNormalE(i,j) &
                  * ( &
                     ! direct diffusion operator over irregular mesh
                     ( diffusionE(i,j)+diffusionE(i+1,j) )*( mesh%Ne(i+1,j)-mesh%Ne(i,j) )* NormalE2(i,j) &
                     ! Cross-diffusion from [Mathur and Murthy (1997)]
                     !TODO: put proper ref here
                     + ShapeFactorTangentE(i,j) &
                     *(diffusionE(i,j)+diffusionE(i+1,j))*(dual%Ne(i,j) - dual%Ne(i,j-1))  &
                                                                                                        ) &
                  + ShapeFactorNormalW(i,j) &
                    *( &
                       ! direct diffusion operator over irregular mesh
                      -( diffusionE(i-1,j)+diffusionE(i,j) )*( mesh%Ne(i,j)-mesh%Ne(i-1,j) )* NormalW2(i,j)  &
                       ! Cross-diffusion from [Mathur and Murthy (1997)]
                      +ShapeFactorTangentW(i,j) &
                      *(diffusionE(i-1,j)+diffusionE(i,j))*( dual%Ne(i-1,j-1) - dual%Ne(i-1,j) ) &
                                                                                                             ) &     
                  + ShapeFactorNormalN(i,j) &
                    *( &
                       ! direct diffusion operator over irregular mesh
                      +( diffusionE(i,j+1)+diffusionE(i,j) )*( mesh%Ne(i,j+1)-mesh%Ne(i,j) )* NormalN2(i,j)  &
                       ! Cross-diffusion from [Mathur and Murthy (1997)]
                      +ShapeFactorTangentN(i,j) &
                      *(diffusionE(i,j+1)+diffusionE(i,j))*( dual%Ne(i-1,j) - dual%Ne(i,j) ) &
                                                                                                         ) &     
                  + ShapeFactorNormalS(i,j) &
                    *( &
                       ! direct diffusion operator over irregular mesh
                      -( diffusionE(i,j-1)+diffusionE(i,j) )*( mesh%Ne(i,j)-mesh%Ne(i,j-1) )* NormalS2(i,j)  &
                       ! Cross-diffusion from [Mathur and Murthy (1997)]
                      +ShapeFactorTangentS(i,j) &
                      *(diffusionE(i,j-1)+diffusionE(i,j))*( dual%Ne(i,j-1) - dual%Ne(i-1,j-1) ) &
                                                                                                         ) ))
            !      +  CrossCoeff*( CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j) ) &
            !        *( diffusionE(i,j)+diffusionE(i+1,j) )*( NeDual(i,j) - NeDual(i,j-1) )/DistDualE(i,j)&
            !        * CellAreaE(i,j)/( CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j) ) & 
                      ! Cross-diffusion from [Mathur and Murthy (1997)]
            !      +  CrossCoeff*( CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j) )*CellAreaW(i,j) &
            !        *(diffusionE(i-1,j)+diffusionE(i,j))*( NeDual(i-1,j-1) - NeDual(i-1,j) ) &
            !        /( (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j) )*DistDualW(i,j)) &
                      ! Cross-diffusion from [Mathur and Murthy (1997)]
            !      +  CrossCoeff*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j) &
            !        *(diffusionE(i,j+1)+diffusionE(i,j))*( NeDual(i-1,j) - NeDual(i,j) ) &
            !        /( (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))*DistDualN(i,j) ) &
                     ! Cross-diffusion from [Mathur and Murthy (1997)]
            !      +  CrossCoeff*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j) &
            !        *(diffusionE(i,j-1)+diffusionE(i,j))*( NeDual(i,j-1) - NeDual(i-1,j-1) ) &
            !        /( (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))*DistDualS(i,j) ) & 
            !       &
            !
            !
                ! convection
!               -((0.5d0*(JeX(i+1,j)+JeX(i,j))*NormalEx(i,j)+0.5d0*(JeY(i+1,j) & 
!               +JeY(i,j))*NormalEy(i,j))*CellAreaE(i,j)+(0.5d0*(JeX(i,j)+JeX(i-1,j))*NormalWx(i,j)+0.5d0*(JeY(i,j) &
!               +JeY(i-1,j))*NormalWy(i,j))*CellAreaW(i,j)+(0.5d0*(JeX(i,j)+JeX(i,j+1))*NormalNx(i,j)+0.5d0*(JeY(i,j) & 
!               +JeY(i,j+1))*NormalNy(i,j))*CellAreaN(i,j)+(0.5d0*(JeX(i,j)+JeX(i,j-1))*NormalSx(i,j)+0.5d0*(JeY(i,j) & 
!               +JeY(i,j-1))*NormalSy(i,j))*CellAreaS(i,j)) &
!
!                   ! Cross-diffusion from Mathur and Murthy 2007 without interpolation
!                   + 0d0*(&
!                     0.5d0*CellAreaE(i,j)*(diffusionE(i,j)+diffusionE(i+1,j)) * ( 0.5d0* (GradNeX(i,j)+GradNeX(i+1,j) )*NormalEx(i,j) + 0.5d0* (GradNeY(i,j)+GradNeY(i+1,j) )*NormalEy(i,j) - 0.5d0*( (GradNeX(i,j)+GradNeX(i+1,j) )*CurviEx(i,j) + ( GradNeY(i,j)+GradNeY(i+1,j) )*CurviEy(i,j) )/(NormalEx(i,j)*CurviEx(i,j)+NormalEy(i,j)*CurviEy(i,j))) &
!                   + 0.5d0*CellAreaW(i,j)*(diffusionE(i,j)+diffusionE(i-1,j)) * ( 0.5d0* (GradNeX(i,j)+GradNeX(i-1,j) )*NormalWx(i,j) + 0.5d0* (GradNeY(i,j)+GradNeY(i-1,j) )*NormalWy(i,j) - 0.5d0*( (GradNeX(i,j)+GradNeX(i-1,j) )*CurviWx(i,j) + ( GradNeY(i,j)+GradNeY(i-1,j) )*CurviWy(i,j) )/(NormalWx(i,j)*CurviWx(i,j)+NormalWy(i,j)*CurviWy(i,j))) &
!                   + 0.5d0*CellAreaN(i,j)*(diffusionE(i,j)+diffusionE(i,j+1)) * ( 0.5d0* (GradNeX(i,j)+GradNeX(i,j+1) )*NormalNx(i,j) + 0.5d0* (GradNeY(i,j)+GradNeY(i,j+1) )*NormalNy(i,j) - 0.5d0*( (GradNeX(i,j)+GradNeX(i,j+1) )*CurviNx(i,j) + ( GradNeY(i,j)+GradNeY(i,j+1) )*CurviNy(i,j) )/(NormalNx(i,j)*CurviNx(i,j)+NormalNy(i,j)*CurviNy(i,j))) &
!                   + 0.5d0*CellAreaS(i,j)*(diffusionE(i,j)+diffusionE(i,j-1)) * ( 0.5d0* (GradNeX(i,j)+GradNeX(i,j-1) )*NormalSx(i,j) + 0.5d0* (GradNeY(i,j)+GradNeY(i,j-1) )*NormalSy(i,j) - 0.5d0*( (GradNeX(i,j)+GradNeX(i,j-1) )*CurviSx(i,j) + ( GradNeY(i,j)+GradNeY(i,j-1) )*CurviSy(i,j) )/(NormalSx(i,j)*CurviSx(i,j)+NormalSy(i,j)*CurviSy(i,j))) &
!                   ) &
            
        !TODO: This is redondant with copy_mesh operation at the begining of the temporal loop
        if(NeOff.eq.1) then
          newmesh%Ne(i,j)=mesh%Ne(i,j); newmesh%Nh(i,j)=mesh%Nh(i,j)
        end if
        
        if(HolesOff.eq.0 .AND. NeOff.eq.0) then
                     
        !TODO: This can be frther optimise
        newmesh%Nh(i,j) = ( &
              (GainsH(i,j)-LossesH(i,j))*CellVol(i,j) & 
              ! drift
!               -((0.5d0*(JhX(i+1,j)+JhX(i,j))*NormalEx(i,j)+0.5d0*(JhY(i+1,j) & 
!               +JhY(i,j))*NormalEy(i,j))*CellAreaE(i,j)+(0.5d0*(JhX(i,j)+JhX(i-1,j))*NormalWx(i,j)+0.5d0*(JhY(i,j) &
!               +JhY(i-1,j))*NormalWy(i,j))*CellAreaW(i,j)+(0.5d0*(JhX(i,j)+JhX(i,j+1))*NormalNx(i,j)+0.5d0*(JhY(i,j) & 
!               +JhY(i,j+1))*NormalNy(i,j))*CellAreaN(i,j)+(0.5d0*(JhX(i,j)+JhX(i,j-1))*NormalSx(i,j)+0.5d0*(JhY(i,j) & 
!               +JhY(i,j-1))*NormalSy(i,j))*CellAreaS(i,j)) &
!               ! diffusion operator on regular mesh
!               +(0.5d0*(Nh(i+1,j)-Nh(i,j))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 & 
!               +y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0)*(diffusionH(i+1,j)+diffusionH(i,j))*CellAreaE(i,j) & 
!               -0.5d0/(x(i,j)**2 & 
!               -2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j)+y(i-1,j)**2)**(0.5d0)*(diffusionH(i-1,j) & 
!               +diffusionH(i,j))*(Nh(i,j)-Nh(i-1,j))*CellAreaW(i,j) &
!               +0.5d0*(diffusionH(i,j+1)+diffusionH(i,j))*(Nh(i,j+1) & 
!               -Nh(i,j))*CellAreaN(i,j)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j)+x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j) & 
!               +y(i,j)**2)**(0.5d0) & 
!               -0.5d0*(diffusionH(i,j)+diffusionH(i,j-1))*(Nh(i,j)-Nh(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 &
!               -2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) & 
              ! diffusion on irregular mesh
              +0.5d0*(& 
              + NormalE2(i,j)*ShapeFactorNormalE(i,j)*(diffusionH(i,j)+diffusionH(i+1,j))*(mesh%Nh(i+1,j)-mesh%Nh(i,j)) &
!               - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(diffusionH(i,j)+diffusionH(i+1,j))*(0.25d0*Nh(i+1,j+1)+0.25d0*Nh(i,j+1)-0.25d0*Nh(i+1,j-1)-0.25d0*Nh(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
              - NormalW2(i,j)*ShapeFactorNormalW(i,j)*(diffusionH(i-1,j)+diffusionH(i,j))*(mesh%Nh(i,j)-mesh%Nh(i-1,j)) &
!               - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(diffusionH(i-1,j)+diffusionH(i,j))*(0.25d0*Nh(i,j+1)+0.25d0*Nh(i-1,j+1)-0.25d0*Nh(i-1,j-1)-0.25d0*Nh(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
              + NormalN2(i,j)*ShapeFactorNormalN(i,j)*(diffusionH(i,j+1)+diffusionH(i,j))*(mesh%Nh(i,j+1)-mesh%Nh(i,j)) &
!               - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(diffusionH(i,j+1)+diffusionH(i,j))*(0.25d0*Nh(i+1,j+1)+0.25d0*Nh(i+1,j)-0.25d0*Nh(i-1,j)-0.25d0*Nh(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
              - NormalS2(i,j)*ShapeFactorNormalS(i,j)*(diffusionH(i,j-1)+diffusionH(i,j))*(mesh%Nh(i,j)-mesh%Nh(i,j-1)) &
!               - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(diffusionH(i,j-1)+diffusionH(i,j))*(0.25d0*Nh(i+1,j)+0.25d0*Nh(i+1,j-1)-0.25d0*Nh(i-1,j-1)-0.25d0*Nh(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
!                   ! Cross-diffusion from [Mathur and Murthy (1997)]
              + ShapeFactorTangentE(i,j)*ShapeFactorNormalE(i,j)*(diffusionH(i,j)+diffusionH(i+1,j)) &
                    *( dual%Nh(i,j) - dual%Nh(i,j-1) ) &
              + ShapeFactorTangentW(i,j)*ShapeFactorNormalW(i,j)*(diffusionH(i-1,j)+diffusionH(i,j)) &
                    *( dual%Nh(i-1,j-1) - dual%Nh(i-1,j) ) &
              + ShapeFactorTangentN(i,j)*ShapeFactorNormalN(i,j)*(diffusionH(i,j+1)+diffusionH(i,j)) &
                    *( dual%Nh(i-1,j) - dual%Nh(i,j) ) &
              + ShapeFactorTangentS(i,j)*ShapeFactorNormalS(i,j)*(diffusionH(i,j-1)+diffusionH(i,j)) &
                    *( dual%Nh(i,j-1) - dual%Nh(i-1,j-1) ) &
              ))*dt*InvCellVol(i,j)+mesh%Nh(i,j)
!         else
!           NhNew(i,j)=NeNew(i,j)
        endif
! ! !  
                     
! form with bug corrected in derivatives and (OmegaX, OmegaY) drift transport included in finite volumes
!     if(ConductivityFix < 2) then
        if(TeOff.ne.1) then 
          if(ConvectionEnergy.eq.0) then

          !TODO: This can be further optimise
          newmesh%Te(i,j) = &
                  0.5d0*(&
                  + NormalE2(i,j)*ShapeFactorNormalE(i,j)*(kappae(i,j)+kappae(i+1,j))*(mesh%Te(i+1,j)-mesh%Te(i,j)) &
                  - NormalW2(i,j)*ShapeFactorNormalW(i,j)*(kappae(i-1,j)+kappae(i,j))*(mesh%Te(i,j)-mesh%Te(i-1,j)) &
                  + NormalN2(i,j)*ShapeFactorNormalN(i,j)*(kappae(i,j+1)+kappae(i,j))*(mesh%Te(i,j+1)-mesh%Te(i,j)) &
                  - NormalS2(i,j)*ShapeFactorNormalS(i,j)*(kappae(i,j-1)+kappae(i,j))*(mesh%Te(i,j)-mesh%Te(i,j-1)) &
                  + ShapeFactorTangentE(i,j)*ShapeFactorNormalE(i,j)*(kappae(i,j)+kappae(i+1,j)) &
                        *( dual%Te(i,j) - dual%Te(i,j-1) ) &
                  + ShapeFactorTangentW(i,j)*ShapeFactorNormalW(i,j)*(kappae(i-1,j)+kappae(i,j)) &
                        *( dual%Te(i-1,j-1) - dual%Te(i-1,j) ) &
                  + ShapeFactorTangentN(i,j)*ShapeFactorNormalN(i,j)*(kappae(i,j+1)+kappae(i,j)) &
                        *( dual%Te(i-1,j) - dual%Te(i,j) ) &
                  + ShapeFactorTangentS(i,j)*ShapeFactorNormalS(i,j)*(kappae(i,j-1)+kappae(i,j)) &
                        *( dual%Te(i,j-1) - dual%Te(i-1,j-1) ) &
                  +2.0d0*(-CouplingE(i,j)+SourceE(i,j))*CellVol(i,j)) & !source
                /Ce(i,j) * dt * InvCellVol(i,j) &
                +mesh%Te(i,j)

            else !convective term included!
            
                UeNew(i,j) = ((SourceUe(i,j)-CouplingE(i,j))*CellVol(i,j) & 
                ! convective term for transport of the energy by the field
                -0.5d0*(((VeX(i+1,j)+VeX(i,j))*NormalEx(i,j)                   &
                        +(VeY(i+1,j)+VeY(i,j))*NormalEy(i,j)) * CellAreaE(i,j) & 
                +       ((VeX(i,j)+VeX(i-1,j))*NormalWx(i,j)                   &
                        +(VeY(i,j)+VeY(i-1,j))*NormalWy(i,j)) * CellAreaW(i,j) & 
                +       ((VeX(i,j)+VeX(i,j+1))*NormalNx(i,j)                   &
                        +(VeY(i,j)+VeY(i,j+1))*NormalNy(i,j)) * CellAreaN(i,j) &
                +       ((VeX(i,j)+VeX(i,j-1))*NormalSx(i,j)                   &
                        +(VeY(i,j)+VeY(i,j-1))*NormalSy(i,j)) * CellAreaS(i,j))&
                ! diffusive term for energy - rewrite correctly
!                 + ( (Ue(i+1,j)-Ue(i,j))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 & 
!                 + y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0)*( (kappae(i+1,j)+kappae(i,j))/(Ce(i+1,j)+Ce(i,j)) )*CellAreaE(i,j) &
!                 - 1d0/(x(i,j)**2-2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j) + & 
!                 y(i-1,j)**2)**(0.5d0)*((kappae(i-1,j)+kappae(i,j))/(Ce(i-1,j)+Ce(i,j)))*(Ue(i,j)-Ue(i-1,j))*CellAreaW(i,j) & 
!                 + 1d0*((kappae(i,j+1)+kappae(i,j))/(Ce(i,j+1)+Ce(i,j)))*(Ue(i,j+1)-Ue(i,j))*CellAreaN(i,j)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j) & 
!                 +x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j)+y(i,j)**2)**(0.5d0) & 
!                 -1d0*((kappae(i,j)+kappae(i,j-1))/(Ce(i,j)+Ce(i,j-1)))*(Ue(i,j)-Ue(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 & 
!                 -2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) &
                )*dt*InvCellVol(i,j) & 
                !TODO: It seems that there is a bug here, as CrossCoef is missing
                +0.5d0*5d0/3d0*( &
                    NormalE2(i,j)*ShapeFactorNormalE(i,j)*(kappae(i,j)/Ce(i,j)+kappae(i+1,j)/Ce(i+1,j))*(Ue(i+1,j)-Ue(i,j))  &
                    !
                  - (CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*ShapeFactorNormalE(i,j)*(kappae(i,j) &
                          /Ce(i,j)+kappae(i+1,j)/Ce(i+1,j))*0.25d0*(Ue(i+1,j+1)+Ue(i,j+1)-Ue(i+1,j-1)-Ue(i,j-1))/DistDualE(i,j) &
                    !
                  + NormalW2(i,j)*ShapeFactorNormalW(i,j)*(kappae(i-1,j)/Ce(i-1,j)+kappae(i,j)/Ce(i,j))*(Ue(i,j)-Ue(i-1,j)) &
                    !
                  - (CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*ShapeFactorNormalW(i,j)*(kappae(i-1,j) &
                          /Ce(i-1,j)+kappae(i,j)/Ce(i,j))*0.25d0*(Ue(i,j+1)+Ue(i-1,j+1)-Ue(i-1,j-1)-Ue(i,j-1))/DistDualW(i,j) &
                    !
                  + NormalN2(i,j)*ShapeFactorNormalN(i,j)*(kappae(i,j+1)/Ce(i,j+1)+kappae(i,j)/Ce(i,j))*(Ue(i,j+1)-Ue(i,j)) &
                    !
                  - (CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*ShapeFactorNormalN(i,j)*(kappae(i,j+1) &
                          /Ce(i,j+1)+kappae(i,j)/Ce(i,j))*0.25d0*(Ue(i+1,j+1)+Ue(i+1,j)-Ue(i-1,j)-Ue(i-1,j+1))/DistDualN(i,j) &
                    !
                  + NormalS2(i,j)*ShapeFactorNormalS(i,j)*(kappae(i,j-1)/Ce(i,j-1)+kappae(i,j)/Ce(i,j))*(Ue(i,j)-Ue(i,j-1)) &
                    !
                  - (CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*ShapeFactorNormalS(i,j)*(kappae(i,j-1) &
                          /Ce(i,j-1)+kappae(i,j)/Ce(i,j))*0.25d0*(Ue(i+1,j)+Ue(i+1,j-1)-Ue(i-1,j-1)-Ue(i-1,j))/DistDualS(i,j) &
!                     Korfiatis 2007 equation
!                     0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappae(i,j)*diffusionE(i,j)/Ne(i,j)+kappae(i+1,j)*diffusionE(i+1,j)/Ne(i+1,j))*(Ne(i+1,j)-Ne(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
!                   - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappae(i,j)*diffusionE(i,j)/Ne(i,j)+kappae(i+1,j)*diffusionE(i+1,j)/Ne(i+1,j))*(0.25d0*Ne(i+1,j+1)+0.25d0*Ne(i,j+1)-0.25d0*Ne(i+1,j-1)-0.25d0*Ne(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
!                   + 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappae(i-1,j)*diffusionE(i-1,j)/Ne(i-1,j)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j)-Ne(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
!                   - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappae(i-1,j)*diffusionE(i-1,j)/Ne(i-1,j)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i,j+1)+0.25d0*Ne(i-1,j+1)-0.25d0*Ne(i-1,j-1)-0.25d0*Ne(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                   + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappae(i,j+1)*diffusionE(i,j+1)/Ne(i,j+1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j+1)-Ne(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
!                   - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappae(i,j+1)*diffusionE(i,j+1)/Ne(i,j+1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i+1,j+1)+0.25d0*Ne(i+1,j)-0.25d0*Ne(i-1,j)-0.25d0*Ne(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
!                   + 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappae(i,j-1)*diffusionE(i,j-1)/Ne(i,j-1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j)-Ne(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
!                   - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappae(i,j-1)*diffusionE(i,j-1)/Ne(i,j-1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i+1,j)+0.25d0*Ne(i+1,j-1)-0.25d0*Ne(i-1,j-1)-0.25d0*Ne(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
                  ) & 
                  *dt*InvCellVol(i,j) &
                +Ue(i,j)
                
                UhNew(i,j) = ((SourceUh(i,j)-CouplingH(i,j))*CellVol(i,j) & 
                ! convective term for transport of the energy by the field
                -0.5d0*(((VhX(i+1,j)+VhX(i,j))*NormalEx(i,j)+(VhY(i+1,j)+VhY(i,j))*NormalEy(i,j)) * CellAreaE(i,j) & 
                +( (VhX(i,j)+VhX(i-1,j))*NormalWx(i,j)+(VhY(i,j)+VhY(i-1,j))*NormalWy(i,j)) * CellAreaW(i,j) & 
                +( (VhX(i,j)+VhX(i,j+1))*NormalNx(i,j)+(VhY(i,j)+VhY(i,j+1))*NormalNy(i,j)) * CellAreaN(i,j) &
                +( (VhX(i,j)+VhX(i,j-1))*NormalSx(i,j)+(VhY(i,j)+VhY(i,j-1))*NormalSy(i,j)) * CellAreaS(i,j)) &
                ! diffusive term for energy - rewrite correctly
!                 + ( (Uh(i+1,j)-Uh(i,j))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 & 
!                 + y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0)*( (kappah(i+1,j)+kappah(i,j))/(Ch(i+1,j)+Ch(i,j)) )*CellAreaE(i,j) &
!                 - 1d0/(x(i,j)**2-2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j) + & 
!                 y(i-1,j)**2)**(0.5d0)*((kappah(i-1,j)+kappah(i,j))/(Ch(i-1,j)+Ch(i,j)))*(Uh(i,j)-Uh(i-1,j))*CellAreaW(i,j) & 
!                 + 1d0*((kappah(i,j+1)+kappah(i,j))/(Ch(i,j+1)+Ch(i,j)))*(Uh(i,j+1)-Uh(i,j))*CellAreaN(i,j)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j) & 
!                 +x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j)+y(i,j)**2)**(0.5d0) & 
!                 -1d0*((kappah(i,j)+kappah(i,j-1))/(Ch(i,j)+Ch(i,j-1)))*(Uh(i,j)-Uh(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 & 
!                 -2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) &
                )*dt*InvCellVol(i,j) & 
                !TODO: It seems that CrossCoef is also missing here...
                +0.5d0*( & 
                    NormalE2(i,j)*ShapeFactorNormalS(i,j)*(kappah(i,j)/Ch(i,j)+kappah(i+1,j)/Ch(i+1,j))*(Uh(i+1,j)-Uh(i,j))    &
                    !
                  - (CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*ShapeFactorNormalE(i,j)*(kappah(i,j)/Ch(i,j)     &
                          +kappah(i+1,j)/Ch(i+1,j))*0.25d0*(Uh(i+1,j+1)+Uh(i,j+1)-Uh(i+1,j-1)-Uh(i,j-1))/DistDualE(i,j)        &
                    !
                  + NormalW2(i,j)*ShapeFactorNormalW(i,j)*(kappah(i-1,j)/Ch(i-1,j)+kappah(i,j)/Ch(i,j))*(Uh(i,j)-Uh(i-1,j))    &
                    !
                  - (CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*ShapeFactorNormalS(i,j)*(kappah(i-1,j)/Ch(i-1,j) &
                          +kappah(i,j)/Ch(i,j))*0.25d0*(Uh(i,j+1)+Uh(i-1,j+1)-Uh(i-1,j-1)-Uh(i,j-1))/DistDualW(i,j) &
                    !
                  + NormalN2(i,j)*ShapeFactorNormalN(i,j)*(kappah(i,j+1)/Ch(i,j+1)+kappah(i,j)/Ch(i,j))*(Uh(i,j+1)-Uh(i,j))    &
                    !
                  - (CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*ShapeFactorNormalN(i,j)*(kappah(i,j+1)/Ch(i,j+1) &
                          +kappah(i,j)/Ch(i,j))*0.25d0*(Uh(i+1,j+1)+Uh(i+1,j)-Uh(i-1,j)-Uh(i-1,j+1))/DistDualN(i,j)            &
                    !
                  + NormalS2(i,j)*ShapeFactorNormalS(i,j)*(kappah(i,j-1)/Ch(i,j-1)+kappah(i,j)/Ch(i,j))*(Uh(i,j)-Uh(i,j-1))    &
                    !
                  - (CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*ShapeFactorNormalS(i,j)*(kappah(i,j-1)/Ch(i,j-1) &
                        +kappah(i,j)/Ch(i,j))*0.25d0*(Uh(i+1,j)+Uh(i+1,j-1)-Uh(i-1,j-1)-Uh(i-1,j))/DistDualS(i,j) &
! ! Korfiatis 2007 equation
!                     0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappah(i,j)*diffusionH(i,j)/Nh(i,j)+kappah(i+1,j)*diffusionH(i+1,j)/Nh(i+1,j))*(Nh(i+1,j)-Nh(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
!                   - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappah(i,j)*diffusionH(i,j)/Nh(i,j)+kappah(i+1,j)*diffusionH(i+1,j)/Nh(i+1,j))*(0.25d0*Nh(i+1,j+1)+0.25d0*Nh(i,j+1)-0.25d0*Nh(i+1,j-1)-0.25d0*Nh(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
!                   + 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappah(i-1,j)*diffusionH(i-1,j)/Nh(i-1,j)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(Nh(i,j)-Nh(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
!                   - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappah(i-1,j)*diffusionH(i-1,j)/Nh(i-1,j)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(0.25d0*Nh(i,j+1)+0.25d0*Nh(i-1,j+1)-0.25d0*Nh(i-1,j-1)-0.25d0*Nh(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                   + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappah(i,j+1)*diffusionH(i,j+1)/Nh(i,j+1)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(Nh(i,j+1)-Nh(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
!                   - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappah(i,j+1)*diffusionH(i,j+1)/Nh(i,j+1)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(0.25d0*Nh(i+1,j+1)+0.25d0*Nh(i+1,j)-0.25d0*Nh(i-1,j)-0.25d0*Nh(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
!                   + 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappah(i,j-1)*diffusionH(i,j-1)/Nh(i,j-1)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(Nh(i,j)-Nh(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
!                   - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappah(i,j-1)*diffusionH(i,j-1)/Nh(i,j-1)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(0.25d0*Nh(i+1,j)+0.25d0*Nh(i+1,j-1)-0.25d0*Nh(i-1,j-1)-0.25d0*Nh(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
                  ) & 
                  *dt*InvCellVol(i,j) &
                +Uh(i,j)
                
            end if


  !  scheme ready for enhanced conductivity and drift of vector (OmegaX, OmegaY). 
  !         if(ConductivityFix < 2) then
        if(HolesOff.eq.0) then
            if(ConvectionEnergy.eq.0) then 

            !TODO: This can be further optimised
            newmesh%Th(i,j) = (0.5d0*( &
                  + NormalE2(i,j)*ShapeFactorNormalE(i,j)*(kappah(i,j)+kappah(i+1,j))*(mesh%Th(i+1,j)-mesh%Th(i,j))  &
!                   - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappah(i,j)+kappah(i+1,j))*(0.25d0*Th(i+1,j+1)+0.25d0*Th(i,j+1)-0.25d0*Th(i+1,j-1)-0.25d0*Th(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
                  - NormalW2(i,j)*ShapeFactorNormalW(i,j)*(kappah(i-1,j)+kappah(i,j))*(mesh%Th(i,j)-mesh%Th(i-1,j))  &
!                   - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappah(i-1,j)+kappah(i,j))*(0.25d0*Th(i,j+1)+0.25d0*Th(i-1,j+1)-0.25d0*Th(i-1,j-1)-0.25d0*Th(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
                  + NormalN2(i,j)*ShapeFactorNormalN(i,j)*(kappah(i,j+1)+kappah(i,j))*(mesh%Th(i,j+1)-mesh%Th(i,j))  &
!                   - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappah(i,j+1)+kappah(i,j))*(0.25d0*Th(i+1,j+1)+0.25d0*Th(i+1,j)-0.25d0*Th(i-1,j)-0.25d0*Th(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
                  - NormalS2(i,j)*ShapeFactorNormalS(i,j)*(kappah(i,j-1)+kappah(i,j))*(mesh%Th(i,j)-mesh%Th(i,j-1))  &
!                   - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappah(i,j-1)+kappah(i,j))*(0.25d0*Th(i+1,j)+0.25d0*Th(i+1,j-1)-0.25d0*Th(i-1,j-1)-0.25d0*Th(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
                  + ShapeFactorTangentE(i,j)*ShapeFactorNormalE(i,j)*(kappah(i,j)+kappah(i+1,j))*( dual%Th(i,j) - dual%Th(i,j-1) ) &
                  + ShapeFactorTangentW(i,j)*ShapeFactorNormalW(i,j)*(kappah(i-1,j)+kappah(i,j))*( dual%Th(i-1,j-1) - dual%Th(i-1,j) ) &
                  + ShapeFactorTangentN(i,j)*ShapeFactorNormalN(i,j)*(kappah(i,j+1)+kappah(i,j))*( dual%Th(i-1,j) - dual%Th(i,j) ) &
                  + ShapeFactorTangentS(i,j)*ShapeFactorNormalS(i,j)*(kappah(i,j-1)+kappah(i,j))*( dual%Th(i,j-1) - dual%Th(i-1,j-1) ) &
                  ) &
                  +(-CouplingH(i,j)+SourceH(i,j))*CellVol(i,j)) & 
                  /Ch(i,j)*dt*InvCellVol(i,j)+mesh%Th(i,j)
            

         else !convection scheme

            !TODO: This must be optmised !
            UhNew(i,j) = ((SourceUh(i,j)-CouplingH(i,j))*CellVol(i,j)-0.5d0*( & 
                   ((VhX(i+1,j)+VhX(i,j))*NormalEx(i,j)+(VhY(i+1,j)+VhY(i,j))*NormalEy(i,j)) * CellAreaE(i,j) & 
                  +((VhX(i,j)+VhX(i-1,j))*NormalWx(i,j)+(VhY(i,j)+VhY(i-1,j))*NormalWy(i,j)) * CellAreaW(i,j) & 
                  +((VhX(i,j)+VhX(i,j+1))*NormalNx(i,j)+(VhY(i,j)+VhY(i,j+1))*NormalNy(i,j)) * CellAreaN(i,j) &
                  +((VhX(i,j)+VhX(i,j-1))*NormalSx(i,j)+(VhY(i,j)+VhY(i,j-1))*NormalSy(i,j)) * CellAreaS(i,j)) &
                  + ( ((kappah(i+1,j)+kappah(i,j))/(Ch(i+1,j)+Ch(i,j))) * (Uh(i+1,j)-Uh(i,j)) & 
                  /(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 + y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0) &
                  *CellAreaE(i,j) &
                  - ((kappah(i-1,j)+kappah(i,j))/(Ch(i-1,j)+Ch(i,j))) * (Uh(i,j)-Uh(i-1,j)) &
                  /(x(i,j)**2-2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j) + y(i-1,j)**2)**(0.5d0) &
                  *CellAreaW(i,j) & 
                  + ((kappah(i,j+1)+kappah(i,j))/(Ch(i,j+1)+Ch(i,j)))*(Uh(i,j+1)-Uh(i,j)) & 
                  /(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j) + x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j)+y(i,j)**2)**(0.5d0) & 
                  *CellAreaN(i,j) &
                  -1d0*((kappah(i,j)+kappah(i,j-1))/(Ch(i,j)+Ch(i,j-1)))*(Uh(i,j)-Uh(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 & 
                  -2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) &
                  )*dt*InvCellVol(i,j)+Uh(i,j)
        end if
        
!       else
!         ThNew(i,j)=TeNew(i,j)
      end if !hole control
    end if !Te/Th control

    if(TsOff.ne.1) then
      ! version with cross-diffusion

      ! first order precision in time
!       TsNew(i,j) =  ( (&
!                     + 0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*(Ts(i+1,j)-Ts(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
!   !                   - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*(0.25d0*Ts(i+1,j+1)+0.25d0*Ts(i,j+1)-0.25d0*Ts(i+1,j-1)-0.25d0*Ts(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
!                     - 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*(Ts(i,j)-Ts(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) &
!   !                   - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*(0.25d0*Ts(i,j+1)+0.25d0*Ts(i-1,j+1)-0.25d0*Ts(i-1,j-1)-0.25d0*Ts(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                     + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*(Ts(i,j+1)-Ts(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) &
!   !                   - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*(0.25d0*Ts(i+1,j+1)+0.25d0*Ts(i+1,j)-0.25d0*Ts(i-1,j)-0.25d0*Ts(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
!                     - 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*(Ts(i,j)-Ts(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j)) &
!   !                   - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*(0.25d0*Ts(i+1,j)+0.25d0*Ts(i+1,j-1)-0.25d0*Ts(i-1,j-1)-0.25d0*Ts(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
!                     + CrossCoeff*( &
!                     + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*( TsDual(i,j) - TsDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
!                     + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*( TsDual(i-1,j-1) - TsDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                     + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*( TsDual(i-1,j) - TsDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
!                     + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*( TsDual(i,j-1) - TsDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
!                     ) &
!                     +(CouplingE(i,j)+0d0*CouplingH(i,j)) * CellVol(i,j)) &
!                     /Cs(i,j) * dt / CellVol(i,j) + SourceS(i,j) &
!                     +Ts(i,j)

        ! second order precision in time
!         TsNew(i,j) =  &
!                     ( (&
!                     + 0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*(Ts(i+1,j)-Ts(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
!   !                   - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*(0.25d0*Ts(i+1,j+1)+0.25d0*Ts(i,j+1)-0.25d0*Ts(i+1,j-1)-0.25d0*Ts(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
!                     - 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*(Ts(i,j)-Ts(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) &
!   !                   - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*(0.25d0*Ts(i,j+1)+0.25d0*Ts(i-1,j+1)-0.25d0*Ts(i-1,j-1)-0.25d0*Ts(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                     + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*(Ts(i,j+1)-Ts(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) &
!   !                   - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*(0.25d0*Ts(i+1,j+1)+0.25d0*Ts(i+1,j)-0.25d0*Ts(i-1,j)-0.25d0*Ts(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
!                     - 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*(Ts(i,j)-Ts(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j)) &
!   !                   - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*(0.25d0*Ts(i+1,j)+0.25d0*Ts(i+1,j-1)-0.25d0*Ts(i-1,j-1)-0.25d0*Ts(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
!                     + CrossCoeff*( &
!                     + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*( TsDual(i,j) - TsDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
!                     + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*( TsDual(i-1,j-1) - TsDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                     + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*( TsDual(i-1,j) - TsDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
!                     + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*( TsDual(i,j-1) - TsDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
!                     )) * 2d0/3d0 * dt/(Cs(i,j)*CellVol(i,j)) +&
!                     (CouplingE(i,j)+0d0*CouplingH(i,j)+SourceS(i,j)) * 2d0/3d0 * dt/Cs(i,j) &
!                     +4d0/3d0*Ts(i,j)-TsOld(i,j)/3d0

        ! third order precision in time, with a constant timestep dt (WORKS)
!         TsNew(i,j)=6d0/11d0*dt/Cs(i,j)/CellVol(i,j)*( (&
!                     + 0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*(Ts(i+1,j)-Ts(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
!   !                   - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*(0.25d0*Ts(i+1,j+1)+0.25d0*Ts(i,j+1)-0.25d0*Ts(i+1,j-1)-0.25d0*Ts(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
!                     - 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*(Ts(i,j)-Ts(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) &
!   !                   - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*(0.25d0*Ts(i,j+1)+0.25d0*Ts(i-1,j+1)-0.25d0*Ts(i-1,j-1)-0.25d0*Ts(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                     + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*(Ts(i,j+1)-Ts(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) &
!   !                   - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*(0.25d0*Ts(i+1,j+1)+0.25d0*Ts(i+1,j)-0.25d0*Ts(i-1,j)-0.25d0*Ts(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
!                     - 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*(Ts(i,j)-Ts(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j)) &
!   !                   - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*(0.25d0*Ts(i+1,j)+0.25d0*Ts(i+1,j-1)-0.25d0*Ts(i-1,j-1)-0.25d0*Ts(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
!                     + CrossCoeff*( &
!                     + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*( TsDual(i,j) - TsDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
!                     + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*( TsDual(i-1,j-1) - TsDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                     + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*( TsDual(i-1,j) - TsDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
!                     + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*( TsDual(i,j-1) - TsDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
!                     )) +6d0/11d0*dt/Cs(i,j)*(CouplingE(i,j)+CouplingH(i,j))+7d0/11d0*Ts(i,j)+18d0/11d0/Cs(i,j)*Ts(i,j)*CsOld(i,j)-9d0/11d0/Cs(i,j)*Ts(i,j)*CsPrev(i,j)+2d0/11d0/Cs(i,j)*Ts(i,j)*CsPrev2(i,j)-9d0/11d0*TsOld(i,j)+2d0/11d0*TsPrev(i,j)

!         ! third order precision in time, with variable timesteps

        newmesh%Ts(i,j) = ((0.5d0*( (&
              + NormalE2(i,j)*ShapeFactorNormalE(i,j)*( kappas(i,j)+kappas(i+1,j))*(mesh%Ts(i+1,j)-mesh%Ts(i,j)) &
                !
              - NormalW2(i,j)*ShapeFactorNormalW(i,j)*( kappas(i-1,j)+kappas(i,j))*(mesh%Ts(i,j)-mesh%Ts(i-1,j)) &
                !
              + NormalN2(i,j)*ShapeFactorNormalN(i,j)*( kappas(i,j+1)+kappas(i,j))*(mesh%Ts(i,j+1)-mesh%Ts(i,j)) &
                !
              - NormalS2(i,j)*ShapeFactorNormalS(i,j)*( kappas(i,j-1)+kappas(i,j))*(mesh%Ts(i,j)-mesh%Ts(i,j-1)) &
                !
              + ShapeFactorTangentE(i,j)*ShapeFactorNormalE(i,j) &
                *( kappas(i,j) + kappas(i+1,j))*( dual%Ts(i,j) - dual%Ts(i,j-1) ) &
	            !
              + ShapeFactorTangentW(i,j)*ShapeFactorNormalW(i,j) &
                *( kappas(i-1,j) + kappas(i,j))*( dual%Ts(i-1,j-1) - dual%Ts(i-1,j) ) &
                !
              + ShapeFactorTangentN(i,j)*ShapeFactorNormalN(i,j) &
                *( kappas(i,j+1) + kappas(i,j))*( dual%Ts(i-1,j) - dual%Ts(i,j) ) &
                !
              + ShapeFactorTangentS(i,j)*ShapeFactorNormalS(i,j) &
               *( kappas(i,j-1) + kappas(i,j))*( dual%Ts(i,j-1) - dual%Ts(i-1,j-1) ) &
                ) &
                    + 2d0*(CouplingE(i,j)+CouplingH(i,j)) * CellVol(i,j) &
                ) * InvCellVol(i,j) &
!                     - ((h1 * h2 + h1 * h3 &
!                     + h2 * h3) / h2 / h1 / h3 * Cs(i,j) - h2 * h3 / h1 &
!                     /(-h3 + h1) / (-h2 + h1) * CsOld(i,j) + h1 * h3 / (-h2 &
!                     + h1) / h2 / (-h3 + h2) * CsPrev(i,j) - h1 * h2 / h3 &
!                     / (h3 ** 2 - h1 * h3 - h2 * h3 + h1 * h2) * CsPrev2(i,j))*0d0 & !20150426-Temporal variation of Cs is killed here.
!                     * Ts(i,j)           
                    ) / Cs(i,j) &
                    + h2 * h3 / h1 / (-h3 + h1)  &
                    / (-h2 + h1) * mesh%Ts(i,j) - h1 * h3 / (-h2 + h1) / h2 / (-h3 + h2) &
                    * TsOld(i,j) + h1 * h2 / h3 / (h3 ** 2 - h1 * h3 - h2 &
                    * h3 + h1 * h2) * TsPrev(i,j)) / (h1 * h2 + h1 * h3 + &
                    h2 * h3) * h2 * h1 * h3
!                    ) / Cs(i,j) - h2 * h3 / h1 / (-h3 + h1) / dt2 * Ts(i,j) &
!                    - h1 * h3 / h2 / dt2 / dt3 * TsOld(i,j) &
!                    + h1 * h2 / (h3 * ( h3 - h1 - h2 ) + h1 * h2) * TsPrev(i,j)) &
!                    / (h1 * h2 + h1 * h3 + h2 * h3) * h2 * h1
                    
    end if

    if(ConvectionEnergy.eq.1) then !define temperatures from energy
      newmesh%Te(i,j) = mesh%Te(i,j) + ((UeNew(i,j) -  Ue(i,j))-1.5d0*kb*mesh%Te(i,j)*(newmesh%Ne(i,j) - mesh%Ne(i,j)) &
          *FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)) ) / Ce(i,j)
      newmesh%Th(i,j) = mesh%Th(i,j) + ((UhNew(i,j) -  Uh(i,j))-1.5d0*kb*mesh%Th(i,j)*(newmesh%Nh(i,j) - mesh%Nh(i,j)) &
          *FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)) ) / Ch(i,j)
    end if
                      

    !TODO: Optimise
    CFLxT(i,j)=kappae(i,j)/Ce(i,j) * dt/(0.5d0*(DistW(i,j)+DistE(i,j)))**2
    CFLyT(i,j)=kappae(i,j)/Ce(i,j) * dt/(0.5d0*(DistN(i,j)+DistS(i,j)))**2
    CFLxTs(i,j)=kappas(i,j)/Cs(i,j) * dt/(0.5d0*(DistW(i,j)+DistE(i,j)))**2
    CFLyTs(i,j)=kappas(i,j)/Cs(i,j) * dt/(0.5d0*(DistN(i,j)+DistS(i,j)))**2
          
!         end if
        
    !TODO: Optimise
    CFLxN(i,j)=diffusionE(i,j)*dt/(x(i,j)-x(i-1,j))**2 !+dt/(x(i,j)-x(i-1,j))*mobilityE(i,j)*sqrt(Ex(i,j)**2+Ey(i,j)**2)
    CFLyN(i,j)=diffusionE(i,j)*dt/(y(i,j)-y(i,j-1))**2 !+dt/(y(i,j)-y(i,j-1))*mobilityE(i,j)*sqrt(Ex(i,j)**2+Ey(i,j)**2)
        
    !TODO:Optimise
    TotalElectrons(i,j)=newmesh%Ne(i,j)*(0.125d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1)) &
                    -0.125d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1)))
    TotalHoles(i,j)=newmesh%Nh(i,j)*(0.125d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1)) &
                    -0.125d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1)))
        
        
    ThermalEnergy(i,j)=Ce(i,j)*mesh%Te(i,j)+Ch(i,j)*mesh%Th(i,j)+Cs(i,j)*mesh%Ts(i,j)
    LaserEnergy(i,j)=OnePhotonIonizationRate0*intensity(i,j)/(1d0-reflectivity(i,j))+ & !energy loss by interband absorption
               TwoPhotonIonizationRate0*intensity(i,j)**2/(1d0-reflectivity(i,j))**2+ & !energy loss by two photon absorption
              (absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*intensity(i,j)/(1d0-reflectivity(i,j)) !energy loss by carrrier heating
        
      end do
    end do
    !$OMP END DO

  !$OMP END PARALLEL
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of parallel section
    
    !BOUNDARY CONDITIONS
    
    do i=1, M !North and South boundaries
      ! finite differences finite difference fashion
      if(DriftOn.eq.0) then
        newmesh%Ne(i,1)=newmesh%Ne(i,2)
        newmesh%Nh(i,1)=newmesh%Nh(i,2)
        newmesh%Ne(i,N)=newmesh%Ne(i,N-1)
        newmesh%Nh(i,N)=newmesh%Nh(i,N-1)
      end if

        UeNew(i,1)=UeNew(i,2)
        UhNew(i,1)=UhNew(i,2)
        newmesh%Te(i,1)=newmesh%Te(i,2)
        newmesh%Th(i,1)=newmesh%Th(i,2)
        newmesh%Ts(i,1)=newmesh%Ts(i,2)

        UeNew(i,N)=UeNew(i,N-1)
        UhNew(i,N)=UhNew(i,N-1)
        newmesh%Te(i,N)=newmesh%Te(i,N-1)
        newmesh%Th(i,N)=newmesh%Th(i,N-1)
        newmesh%Ts(i,N)=newmesh%Ts(i,N-1)

        ! includes also the corners... WHy are not they written?
        
!         potential(i,1)=0d0 !(0d0,0d0)
!                potential(i,N)=0d0 !(0d0,0d0)
    end do
    
    do i=2,M-1
        
         ! boundary condition v.n = 0 on boundaries. 
        ! NORTH
        
         GradNeX(i,N) = GradNeX(i,N-1)
         GradNeY(i,N) = GradNeY(i,N-1)

! 
        ! SOUTH
         GradNeX(i,1) = GradNeX(i,2)
         GradNeY(i,1) = GradNeY(i,2)

     end do
    
    
      do j=2, N-1 !West and East boundaries
      ! finite differences bad fashion
        if(DriftOn.eq.0) then
                newmesh%Ne(1,j)=newmesh%Ne(2,j)
                newmesh%Nh(1,j)=newmesh%Nh(2,j)
        end if
        
        UeNew(1,j)=UeNew(2,j)
        UhNew(1,j)=UhNew(2,j)
        newmesh%Te(1,j)=newmesh%Te(2,j)
        newmesh%Th(1,j)=newmesh%Th(2,j)
        newmesh%Ts(1,j)=newmesh%Ts(2,j)
  !       potential(1,j)=0d0 !(0d0, 0d0)
  !       potential(M,j)=potential0 !(potential0, 0d0)
  
        ! conditions on the cone base - most important        
        if(DriftOn.eq.0) then
                newmesh%Ne(M,j)=newmesh%Ne(M-1,j) !Ne0 replace by dynamic bnd condition with flux equal to the one of cell(M-1,j)
                newmesh%Nh(M,j)=newmesh%Nh(M-1,j) !Nh0 replace by dynamic bnd condition with flux equal to the one of cell(M-1,j)
        end if
        
        !outlet condition on density and energy
!         write(*,*) CellAreaE(M,j), DistE(M-2,j)
        

!         NeNew(M,j) = -0.5d0*(diffusionE(M-2,j)+diffusionE(M-1,j))*(Ne(M-1,j)-Ne(M-2,j))/DistE(M-2,j)/(-0.5d0*diffusionE(M-1,j)-0.5d0*diffusionE(M,j))/DistE(M-1,j)+Ne(M-1,j)
!         NhNew(M,j) = -0.5d0*(diffusionH(M-2,j)+diffusionH(M-1,j))*(Nh(M-1,j)-Nh(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*diffusionH(M-1,j)-0.5d0*diffusionH(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Nh(M-1,j)
        
!         UeNew(M,j)=UeNew(M-1,j)
!         UhNew(M,j)=UhNew(M-1,j)
        
        newmesh%Te(M,j)=newmesh%Te(M-1,j) !Tout
        newmesh%Th(M,j)=newmesh%Th(M-1,j) !Tout
        newmesh%Ts(M,j)=newmesh%Ts(M-1,j) ! Tout !cooling by diffusion from outside, TsNew(M-1,j)
        
!         TeNew(M,j) = -0.5d0*(kappae(M-2,j)+kappae(M-1,j))*(Te(M-1,j)-Te(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappae(M-1,j)-0.5d0*kappae(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Te(M-1,j)
!         ThNew(M,j) = -0.5d0*(kappah(M-2,j)+kappah(M-1,j))*(Th(M-1,j)-Th(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappah(M-1,j)-0.5d0*kappah(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Th(M-1,j)
!         TsNew(M,j) = -0.5d0*(kappas(M-2,j)+kappas(M-1,j))*(Ts(M-1,j)-Ts(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappas(M-1,j)-0.5d0*kappas(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Ts(M-1,j)
        
      end do
      

!       if(DriftOn.eq.1) then
          
          do j=2,N-1
!         
          ! WEST
          
          GradNeX(1,j) = GradNeX(2,j)
          GradNeY(1,j) = GradNeY(2,j)

         ! EAST
          
          GradNeX(M,j) = GradNeX(M-1,j)
          GradNeY(M,j) = GradNeY(M-1,j)

        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CHECKING the results

    TotalMeshVolume=0d0
!     ElectronEnergy=0d0
!     HoleEnergy=0d0
!     LatticeEnergy=0d0
    
    !TODO: Replace with Fortran native min and max functions

    do i=1,M

      do j=1,N

        NeTotal=NeTotal + mesh%Ne(i,j) * CellVol(i,j)
        NhTotal=NhTotal + mesh%Nh(i,j) * CellVol(i,j)
      
        if(MaxHeating(i,j) < mesh%Ts(i,j) .AND. t > 100d0*tau) then
           MaxHeating(i,j)=mesh%Ts(i,j)
           MaxHeatingTime(i,j)=t
        end if

              if(maxCFLxN < CFLxN(i,j)) then 
          maxCFLxN=CFLxN(i,j)
        end if
        
        if(maxCFLyN < CFLyN(i,j)) then 
          maxCFLyN=CFLyN(i,j)
        end if
        
        if(maxCFLxT < CFLxT(i,j)) then 
          maxCFLxT=CFLxT(i,j)
        end if
        
        if(maxCFLyT < CFLyT(i,j)) then 
          maxCFLyT=CFLyT(i,j)
        end if
        
        if(maxCFLxTs < CFLxTs(i,j)) then 
          maxCFLxTs=CFLxTs(i,j)
        end if
        
        if(maxCFLyTs < CFLyTs(i,j)) then 
          maxCFLyTs=CFLyTs(i,j)
        end if
        
        if(maxTe < newmesh%Te(i,j)) then
          maxTe=newmesh%Te(i,j)
        end if
        
        if(minTe > newmesh%Te(i,j)) then
          minTe=newmesh%Te(i,j)
        end if
        
        if(maxTh < newmesh%Th(i,j)) then
          maxTh=newmesh%Th(i,j)
        end if        
        
        if(minTh > newmesh%Th(i,j)) then
          minTh=newmesh%Th(i,j)
        end if

        if(maxTs < newmesh%Ts(i,j)) then
          maxTs=newmesh%Ts(i,j)
        end if
        
        if(minTs > newmesh%Ts(i,j)) then
          minTs=newmesh%Ts(i,j)
        end if
        
        if(maxNe < newmesh%Ne(i,j)) then
          maxNe=newmesh%Ne(i,j)
        end if
        
        if(minNe > newmesh%Ne(i,j)) then
          minNe=newmesh%Ne(i,j)
        end if
        
        if(maxNh < newmesh%Nh(i,j)) then
          maxNh=newmesh%Nh(i,j)
        end if
        
        if(minNh > newmesh%Nh(i,j)) then
          minNh=newmesh%Nh(i,j)
        end if
        
        if(maxIntensity < intensity(i,j)) then 
          maxIntensity=intensity(i,j)
        end if
        
        if(maxSourceE < SourceE(i,j)) then 
          maxSourceE=SourceE(i,j)
        end if
        
        if(maxGainsE < GainsE(i,j)) then 
          maxGainsE=GainsE(i,j)
        end if
                
        if(maxSourceH < SourceH(i,j)) then 
          maxSourceH=SourceH(i,j)
        end if
        
        if(maxGainsH < GainsH(i,j)) then 
          maxGainsH=GainsH(i,j)
        end if
        
        if(maxGap < Egap(i,j)) then 
          maxGap=Egap(i,j)
        end if
        
        if(maxDiffNe < diffNe(i,j)) then 
          maxDiffNe=diffNe(i,j)
        end if
        
        if(maxDiffNh < diffNh(i,j)) then 
          maxDiffNh=diffNh(i,j)
        end if
        
        if(real(maxFermiIndexE) < real(FermiIndexE(i,j))) then
          maxFermiIndexE=FermiIndexE(i,j)
        end if
        
        if(real(maxFermiIndexH) < real(FermiIndexH(i,j))) then
          maxFermiIndexH=FermiIndexH(i,j)
        end if
        
!         if(ConductivityFix.eq.-1) then
!            TeNew(i,j)=Tout; ThNew(i,j)=Tout; 
!         end if
        
        TotalNumOfE=TotalNumOfE + TotalElectrons(i,j)
        TotalNumOfH=TotalNumOfH + TotalHoles(i,j)
        TotalThermalEnergy=TotalThermalEnergy+ThermalEnergy(i,j)
        TotalLaserEnergy=TotalLaserEnergy+LaserEnergy(i,j)
        TotalMeshVolume=TotalMeshVolume+CellVol(i,j)

!         CAUTION: These definitions are erroneously including initial temperature into account. 
!         ElectronEnergy=ElectronEnergy+Ce(i,j)*Te(i,j)*CellVol(i,j)
!         HoleEnergy=HoleEnergy+Ch(i,j)*Th(i,j)*CellVol(i,j)
!         LatticeEnergy=LatticeEnergy+Cs(i,j)*Ts(i,j)*CellVol(i,j)

        ! calculation of the absorbed laser energy involved in the simulated slice !
        if((i.eq.1) .AND. (j.eq.(N/2))) then 
	IntensityEnergy=IntensityEnergy+(OnePhotonIonizationRate0+absorptionDrudeE(i,j) & 
                  +absorptionDrudeH(i,j))*intensity(i,j)*CellVol(i,j)*dt
        end if
        LaserIntensityEnergy=LaserIntensityEnergy+(OnePhotonIonizationRate0+absorptionDrudeE(i,j) & 
                  +absorptionDrudeH(i,j))*I0*real(sqrt(Dielectric(i,j))) & 
                  *exp(-.5d0*((t-t0)/sigmaTau)**2)*CellVol(i,j)*dt
        
        ! calculation of the energy contained in the solid
        ElectronEnergy=ElectronEnergy &
          + (Ce(i,j)*(newmesh%Te(i,j)-mesh%Te(i,j))+(Ce(i,j)-CeOld(i,j))*mesh%Te(i,j)) * CellVol(i,j) & !kinetic energy
          + (mesh%Ne(i,j)*(EgapValue(newmesh%Ne(i,j), newmesh%Ts(i,j))-EgapValue(mesh%Ne(i,j),mesh%Ts(i,j)))  &
              + Egap(i,j)*(newmesh%Ne(i,j)-mesh%Ne(i,j))) * CellVol(i,j) !potential energy
              
        ElectronKineticEnergy=ElectronKineticEnergy+(Ce(i,j)*(newmesh%Te(i,j)-mesh%Te(i,j))+(Ce(i,j)-CeOld(i,j))*mesh%Te(i,j)) * CellVol(i,j) !kinetic energy
        ElectronPotentialEnergy=ElectronPotentialEnergy+(mesh%Ne(i,j)*(EgapValue(newmesh%Ne(i,j), newmesh%Ts(i,j))-EgapValue(mesh%Ne(i,j),mesh%Ts(i,j)))  &
              + Egap(i,j)*(newmesh%Ne(i,j)-mesh%Ne(i,j))) * CellVol(i,j) !potential energy
!         ElectronEnergy=ElectronKineticEnergy+ElectronPotentialEnergy !already summed over time
        
        HoleEnergy=HoleEnergy+(Ch(i,j)*(newmesh%Th(i,j)-mesh%Th(i,j))+(Ch(i,j)-ChOld(i,j))*mesh%Th(i,j)) * CellVol(i,j) !kinetic energy
        
        LatticeEnergy=LatticeEnergy+((Cs(i,j)*(newmesh%Ts(i,j)-mesh%Ts(i,j)))+0d0*(Cs(i,j)-CsOld(i,j))*mesh%Ts(i,j))*CellVol(i,j) !dCs/dt=0, 20150426, TJYD.
        
        !TODO: What the f...? I do not understant
        if(mesh%Te(i,j).ne.mesh%Te(i,j)) then
          write(95,*) "Divergence of Te at t=", t, "x(",i,j,")=", x(i,j), "y(",i,j,")=",y(i,j)
          Diverged=.true.
        end if
        if((mesh%Th(i,j)).ne.mesh%Th(i,j)) then
          write(95,*) "Divergence of Th at t=", t, "x(",i,j,")=", x(i,j), "y(",i,j,")=",y(i,j)
          Diverged=.true.
        end if
        if((mesh%Ts(i,j)).ne.mesh%Ts(i,j)) then
          write(95,*) "Divergence of Ts at t=", t, "x(",i,j,")=", x(i,j), "y(",i,j,")=",y(i,j)
          Diverged=.true.
        end if
        if((mesh%Ne(i,j)).ne.mesh%Ne(i,j)) then
          write(95,*) "Divergence of Ne at t=", t, "x(",i,j,")=", x(i,j), "y(",i,j,")=",y(i,j)
          Diverged=.true.
        end if
        if((mesh%Nh(i,j)).ne.mesh%Nh(i,j)) then
          write(95,*) "Divergence of Nh at t=", t, "x(",i,j,")=", x(i,j), "y(",i,j,")=",y(i,j)
          Diverged=.true.
        end if

        if(maxCFLxT.gt.1d0 .OR. maxCFLyT.gt.1d0) then 
          write(95,*) "Bad convergence for Te,Th. t=", t, "(CFLx,CFLy)=", maxCFLxT, maxCFLyT
          Diverged=.true.
        end if
        if(maxCFLxN.gt.1d0 .OR. maxCFLyN.gt.1d0) then
          write(95,*) "Bad convergence for Ne,Nh. t=", t, "(CFLx,CFLy)=", maxCFLxN, maxCFLyN
          Diverged=.true.
        end if
        
        if(real(FermiIndexE(i,j)) > real(FermiMaxLines) .OR. real(FermiIndexE(i,j)) < 1d0) then
          write(*,*) "t,i,j,FermiIndexE(i,j)=", t,i,j,FermiIndexE(i,j)
        end if
        if(real(FermiIndexH(i,j)) > real(FermiMaxLines) .OR. real(FermiIndexH(i,j)) < 1d0) then
          write(*,*) "t,i,j,FermiIndexH(i,j)=", t,i,j,FermiIndexH(i,j)
        end if
        
        if(real(FermiRatioE(i,j)) < 0d0 .OR. real(FermiRatioH(i,j)) < 0d0) then
          write(*,*) "Problem in DOS or Ne. DOS(i,j)=", i,j,DOSe(i,j), DOSh(i,j), "Ne,h(i,j)=", mesh%Ne(i,j), mesh%Nh(i,j)
        end if
        
        call flush(95)
        
        if(Diverged .eqv. .true.) then
          write(*,*) "Divergence detected. Please check error.dat for more information."
          stop
        end if

              ! lets change dt when fast reponse is finished in order to catch the long one. 

        diffNe(i,j)=(newmesh%Ne(i,j)-mesh%Ne(i,j))/dt
        diffNh(i,j)=(newmesh%Nh(i,j)-mesh%Nh(i,j))/dt
        
                
        if(mod(nbiter,iterOut*iterOutMaps).eq.0) then 
          write(97,887, advance="yes") t, x(i,j), y(i,j), intensity(i,j), mesh%Te(i,j), & !5
                        mesh%Th(i,j), mesh%Ts(i,j), mesh%Ne(i,j), mesh%Nh(i,j), reflectivity(i,j), & !10
                        absorptionDrudeE(i,j), absorptionDrudeH(i,j), diffNe(i,j), diffNh(i,j), TotalElectrons(i,j), & !15
                        TotalHoles(i,j), real(FermiIndexE(i,j)), REAL(FermiIndexH(i,j)), FermiRatioE(i,j), FermiRatioH(i,j), & !20
                        SourceE(i,j), SourceH(i,j), GainsE(i,j), GainsH(i,j), LossesE(i,j), & !25
                        LossesH(i,j), real(DielectricDrudeE(i,j)), aimag(DielectricDrudeE(i,j)), Egap(i,j), real(Dielectric(i,j)), & !30
                        aimag(Dielectric(i,j)), MaxHeatingTime(i,j), MaxHeating(i,j), real(potentialNeedle(i,j)), Ex(i,j), & !35
                        Ey(i,j), diffusionE(i,j), diffusionH(i,j), GradNeX(i,j), GradNeY(i,j), &!40
                        real(EintField(i,j)), aimag(EintField(i,j)), EintFieldR(i,j), EintFieldI(i,j), phiMie(i,j), & !45
                        Radius(i,j)
                        
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
!          write(97,886, advance='yes')
          call flush(97); 
        end if
        
      end do !on Y
      
       if(mod(nbiter,iterOut*iterOutMaps).eq.0) then 
        write(97,886, advance="yes")
886        FORMAT (3x)
       end if
    end do !on X


    ! saving the timesteps of several previous steps (used for the high order calculation of d/dt).
    dt4=dt3;
    dt3=dt2;
    dt2=dt; 
    ! chaning the timestep based on known behavior of the system
    if(AdaptativeTimeStep.eq.1) then
      if((t>1d1*tau*coeffDilaDt) .AND. (dt.eq.dt0) .AND. &
        (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
        dt=10d0*dt0
      else if ((t > 0.5d2*tau*coeffDilaDt) .AND. (dt.eq.10d0*dt0) .AND. (maxCFLxN+maxCFLyN &
        + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
        dt=40d0*dt0
      else if ((t > 1d3*tau*coeffDilaDt) .AND. (dt.eq.40d0*dt0) .AND. (maxCFLxN+maxCFLyN + maxCFLxT &
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
          write(101,889, advance="yes") t, xP(i,j), yP(i,j), real(potential(i,j)), real(ExPoisson(i,j)), & !
                real(EyPoisson(i,j)), DielectricStatic(i,j), NeP(i,j), NhP(i,j)
  889 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, &
  1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5) 
          call flush(101); 
        end do
      end do
      
!       write(91,'(100E14.5)') potential !Xvector !Amatrix
!     write(91,*)
      
    end if
       
    ! output to files
    call cpu_time(calc_time_3)
    
    cpu_timestep_duration = (calc_time_3-calc_time_begin) / real(nbiter)
    cpuefficiency=real(nbiter)/(calc_time_3-calc_time_begin)*real(nthreads)
    
    if(mod(nbiter,iterOut).eq.0) then 
      
      write(105,892, advance="YES") t, IntensityEnergy, ElectronEnergy, HoleEnergy, LatticeEnergy, & !5
          TotalMeshVolume, LaserIntensityEnergy, ElectronKineticEnergy, ElectronPotentialEnergy !9
      
892 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
      
      write(98,888, advance="YES") t, maxTe, maxTh, maxTs, maxNe, &         !5
                    maxNh, maxIntensity, TotalLaserEnergy, TotalThermalEnergy, cpuefficiency, &        !10
                    maxSourceE, maxGainsE, maxSourceH, maxGainsH, maxGap, &        !15
                    maxDiffNe, maxDiffNh, TotalNumOfE, TotalNumOfH, real(maxFermiIndexE), &        !20
                    real(maxFermiIndexH), NeTotal, NhTotal, maxCFLxT, maxCFLyT, &        !25
                    maxCFLxN, maxCFLyN, maxCFLxTs, maxCFLyTs, IntensityEnergy, &        !30
                    TotalMeshVolume, ElectronEnergy, HoleEnergy, LatticeEnergy, LaserIntensityEnergy, &         !35
                    ElectronKineticEnergy, ElectronPotentialEnergy, cpu_timestep_duration    !38
                    
888 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1F12.8, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E19.11, 3x, 1E19.11, 3x, 1E19.11, 3x, 1E19.11, 3x, 1E19.11, &
3x, 1E19.11, 3x, 1E19.11, 3x, 1E12.5)
                
        write(94,884, advance="YES") t, mesh%Te(1,N/2), mesh%Th(1,N/2), mesh%Ts(1,N/2), mesh%Ne(1,N/2), &                        !5
              mesh%Nh(1,N/2), intensity(1,N/2), TotalLaserEnergy, TotalThermalEnergy, cpuefficiency, &        !10
              SourceE(1,N/2), GainsE(1,N/2), SourceH(1,N/2), GainsH(1,N/2), Egap(1,N/2), &                !15
              diffNe(1,N/2), diffNh(1,N/2), real(FermiIndexE(1,N/2)), real(FermiIndexH(1,N/2)), Ce(2,N/2), &                !20
              CeOld(2,N/2), Ch(2,N/2), ChOld(2,N/2), Cs(2,N/2), CsOld(2,N/2), &                                !25
              cpu_timestep_duration !26
              
884 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5)

        write(93,883, advance="YES") t, mesh%Te(M/2,N), mesh%Th(M/2,N), mesh%Ts(M/2,N), mesh%Ne(M/2,N), &
              mesh%Nh(M/2,N), intensity(M/2,N), TotalLaserEnergy, TotalThermalEnergy, cpuefficiency, &
              SourceE(M/2,N), GainsE(M/2,N), SourceH(M/2,N), GainsH(M/2,N), Egap(M/2,N), &
              diffNe(M/2,N), diffNh(M/2,N), real(FermiIndexE(M/2,N)), real(FermiIndexH(M/2,N))
              
883 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)

        write(92,882, advance="YES") t, mesh%Te(M/2,1), mesh%Th(M/2,1), mesh%Ne(M/2,1), &
              mesh%Nh(M/2,1), intensity(M/2,1), TotalLaserEnergy, TotalThermalEnergy, cpuefficiency, &
              SourceE(M/2,1), GainsE(M/2,1), SourceH(M/2,1), GainsH(M/2,1), Egap(M/2,1), &
              diffNe(M/2,1), diffNh(M/2,1), real(FermiIndexE(M/2,1)), real(FermiIndexH(M/2,1))
              
882 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)

    end if
    
    if(nbiter < 2) then
      do i=1,M
        do j=1,N
          write(102, 881, advance="YES") i, j, x(i,j), y(i,j), NormalNx(i,j), & !5
                                  NormalNy(i,j), NormalSx(i,j), NormalSy(i,j), NormalEx(i,j), NormalEy(i,j), & !10
                                  NormalWx(i,j), NormalWy(i,j), CellAreaN(i,j), CellAreaS(i,j), CellAreaE(i,j), & !15
                                  CellAreaW(i,j), CellVol(i,j), TangentNx(i,j), TangentNy(i,j), TangentSx(i,j), & !20
                                  TangentSy(i,j), TangentEx(i,j), TangentEy(i,j), TangentWx(i,j), TangentWy(i,j), & !25
                                  CurviNx(i,j), CurviNy(i,j), CurviSx(i,j), CurviSy(i,j), CurviEx(i,j), & !30
                                  CurviEy(i,j), CurviWx(i,j), CurviWy(i,j), DistN(i,j), DistS(i,j), & !35
                                  DistE(i,j), DistW(i,j), DistDualN(i,j), DistDualS(i,j), DistDualE(i,j),  & !40
                                  DistDualW(i,j), xDualSW(i,j), yDualSW(i,j), xDualSE(i,j), yDualSE(i,j), & !45
                                  xDualNE(i,j), yDualNE(i,j), xDualNW(i,j), yDualNW(i,j) !, xDual(i,j), & !50
!                                  yDual(i,j)
                                  
881 FORMAT (I3, 3x, I3, 3x, 1E16.8, 3x, 1E16.8, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
        end do
      end do
    end if
    
890 FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
    
    ! write the functions on Dual Mesh
    if(mod(nbiter,iterOut*iterOutMaps).eq.0) then
      do i=1,M-1
        do j=1,N-1
          
          write(103, 890, advance="YES") t, xDual(i,j), yDual(i,j), dual%Te(i,j), dual%Th(i,j), & !5
                                          dual%Ts(i,j), dual%Ne(i,j), dual%Nh(i,j), intensityDual(i,j) !9
                  
        end do
      end do        
    end if
    call flush(98); call flush(94); call flush(93); call flush(101)
  end do
  
  call releasemesh(mesh)
  call releasemesh(dual)
  call releasemesh(newmesh)

  contains 
  
    function DensityOfStateE(Te)
    implicit none
    real(8) DensityOfStateE, Te
      DensityOfStateE=2d0*(meDOS*kb*Te/(2d0*pi*hbar**2))**(1.5d0)
      return
    end function DensityOfStateE
    
    function DensityOfStateH(Th)
    implicit none
    real(8) DensityOfStateH, Th
      DensityOfStateH=2d0*(mhDOS*kb*Th/(2d0*pi*hbar**2))**(1.5d0)
      return
    end function DensityOfStateH
  
    subroutine TabCreateFL !(FermiTableE, FermiTableH)
      implicit none
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

    function DielectricConstant(lambda)
      complex(8) :: DielectricConstant
      real(8) lambda
      
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
    
    function DielectricFunction(lambda, epsilonInf, ne, nuColl)
      complex(8) :: DielectricFunction
      complex(8) epsilonInf
      real(8) lambda, ne, nuColl, omegape
!      real(8) ne, nuColl, omegape
      
      omegape=sqrt(ne*ec**2/me/epsilon0)
      DielectricFunction=epsilonInf-(omegape/omegaLaser)**2*Unit/(Unit+Imaginary*nuColl/omegaLaser)
      return
    end function DielectricFunction
    
    function DielectricFunctionDrude(density, Collision, mass)
      complex(8) :: DielectricFunctionDrude
      real(8) density, Collision, omegape, mass
      
      omegape=sqrt(density*ec**2/(mass*epsilon0))
      DielectricFunctionDrude=Unit-Unit*(omegape/omegaLaser)**2*Unit/(Unit+Imaginary*Collision/omegaLaser)
      return
    end function DielectricFunctionDrude
    
    function ephCollisionFrequency(ne)
      real(8) :: ephCollisionFrequency, ne
      real(8) nth
      nth=6.02d26 !m-3 (Sjodin, PRL 1998)
      ephCollisionFrequency=((240d-15)*(1d0+(ne/nth)**2))**(-1d0)
!       CollisionFrequency=1d14 ! 
      ! CollisionFrequency=1d13 !
      !CollisionFrequency=5d13 ! 
      return
    end function ephCollisionFrequency
  
    function CollisionFrequency()
      real(8) :: CollisionFrequency    
      CollisionFrequency=1d15
      return
    end function CollisionFrequency
  
    function ImpactIonizationRate(Te, Ne, Ts)
      real(8) :: ImpactIonizationRate
      real(8) Te, Ne, Ts
      ImpactIonizationRate=3.6d10*exp(-1.5d0*EgapValue(Ne, Ts)/kb/Te)
      if(ImpactOff.eq.1) then
        ImpactIonizationRate=0d0
      end if
      return
    end function ImpactIonizationRate
    
    function OnePhotonIonizationRate(lambda,epsilonLinear)
      real(8) :: OnePhotonIonizationRate
      real(8) lambda
      complex(8) epsilonLinear
!       OnePhotonIonizationRate=4d0*pi/lambda*aimag(sqrt(epsilonLinear))
      OnePhotonIonizationRate = 3.4536819356d6 !extracted from WC Dash and R Newman, Phys Rev 99, 1151 (1955)
      return
    end function OnePhotonIonizationRate
    
    function TwoPhotonIonizationRate(lambda)
      real(8) :: TwoPhotonIonizationRate
      real(8) lambda
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
    
    function EgapValue(Ne, Ts)
      real(8) :: EgapValue
      real(8) Ne, Ts

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

    function LatticeHeatCapacity(T)
      implicit none
      real(8) :: LatticeHeatCapacity
      real(8) T
!       LatticeHeatCapacity=1d3*SiDensity*0.2703d0/(exp(63.456d0/T)+0.84586d0) !bad fit ...
!         LatticeHeatCapacity=1d3*SiDensity*0.412920554599445d0/(exp(88.1830102582422d0/T)-0.676494557497076d0) !!better fit on Flubacher BUT INDUCES A SUPER BUG (+170 K with 3rd order time integration).
!         LatticeHeatCapacity=1d6*(1.978d0+3.54d-4*T-3.68d0*T**(-2)) !Driel 1987 - not very good, BUT WORKS.
!       LatticeHeatCapacity=1d3*SiDensity*(0.899d0*dexp(5.455d-05*T)-0.959d0*dexp(-0.004218d0*T))! very good exp fit on Okothin, BUT INDUCES A nonlinearity at the beginning (+50 K with 3rd order integration)
!         LatticeHeatCapacity=1D3*SiDensity*(1.239d0*sin(0.001413d0*T-0.1806d0) + 0.3168d0*sin(0.003343d0*T+0.7648d0) + 0.01947d0*sin(0.00904d0*T+0.4528d0) + 0.04943d0*sin(0.007262d0*T-0.6642d0)) !Fitted on Otokhin, but is it stable ?
!         LatticeHeatCapacity=1D3*SiDensity*(2.36d-16*T**5 -1.707d-12*T**4 + 4.619d-09*T**3 -5.912d-06*T**2 + 0.003733d0*T -0.0494d0) ! 5th order polynomial fit on Okhonin
!         LatticeHeatCapacity=1d3*SiDensity*(0.4135d0*T-0.4071d0*T**1.002d0) !Driel style (1)
        LatticeHeatCapacity=1d3*SiDensity*(-0.003592d0*T+0.01458d0*T**0.8316d0) !Driel style (2, better ?)
    end function LatticeHeatCapacity

     
    function FermiIndex(NeNc)
      implicit none
      ! Input: Value of density/DOS
      ! returns the index to take in the Fermi files
      real(8) NeNc, NeNc0, dNeNc
      integer(8) FermiIndex
      NeNc0=1d-38
      dNeNc=1.03d0 !NeNc=NeNc0*dNeNc**n
      
      FermiIndex=nint(log10(NeNc/NeNc0)/log10(dNeNc)+1d0)
      if(FermiIndex < 1 .OR. FermiIndex > FermiMaxLines) then
        write(*,*) "FermiIndex problem: NeNc=", NeNc, "FermiIndex=", FermiIndex
      end if
      return
    end function FermiIndex



    function MieScattering(r, phi, radius, dielectric)
      implicit none

      complex(8) :: MieScattering, dielectric, BesselJ
      complex(8) total

      real(8) :: r, phi, radius, k=2d0*pi/lambda
      real(8) ireal
      integer(8) i, j


      total=Zero
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
        total=total+Imaginary**ireal * exp(Imaginary*ireal*phi) * BesselJ(ireal, &
                sqrt(dielectric)*k*r, besselArray) * MieCoeff1(ireal, radius, dielectric, besselArray)
      end do

      MieScattering=total
      return
    end function MieScattering

    function MieScatteringTE1(r, phi, radius, dielectric, besselArray)
      implicit none

      complex(8) :: MieScatteringTE1, dielectric, BesselJ
      complex(8) :: total

      real(8) :: r, phi, radius, k=2d0*pi/lambda
      real(8) :: ireal
      integer(8) :: i, j, besselArray

      !!Careful !! This function is very sensitive to noise.

      total=Zero

      do i=1, 2*maxBesselOrder+1
        ireal=real(i-maxBesselOrder-1) !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
!         write(*,*) i, ireal
        total=total + ( Imaginary**ireal * exp(Imaginary*ireal*phi) * BesselJ(ireal, &
              sqrt(dielectric)*k*r, besselArray) * ireal * MieCoeff3(ireal, radius, dielectric, besselArray) ) !original !!
      end do
      if(r.eq.0d0) then
        MieScatteringTE1=Zero
      else
        MieScatteringTE1=-total/(dielectric*k*r)
!         MieScatteringTE1=Zero
      end if
      return
    end function MieScatteringTE1

    function MieScatteringTE2(r, phi, radius, dielectric, besselArray)
      implicit none

      complex(8) :: MieScatteringTE2, dielectric, BesselJprime
      complex(8) :: total

      real(8) :: r, phi, radius, k=2d0*pi/lambda
      real(8) :: ireal
      integer(8) :: i, j

      integer(8) :: besselArray

      total=Zero

      do i=1, 2*maxBesselOrder+1
        ireal=real(i-maxBesselOrder-1) !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
!         write(*,*) i, ireal
        total=total + ( Imaginary**ireal * exp(Imaginary*ireal*phi) * &
              BesselJprime(ireal, sqrt(dielectric)*k*r, besselArray) * MieCoeff3(ireal, radius, dielectric, besselArray) )
      end do

      MieScatteringTE2=-total*Imaginary/sqrt(dielectric)
!         MieScatteringTE2=Zero !debug
      return
    end function MieScatteringTE2

    function MieCoeff1(order, radius, dielectric, besselArray)
      implicit none
      complex(8) MieCoeff1, dielectric, Hankel1, BesselJ
      real(8) :: k=2d0*pi/lambda, radius
      real(8) :: order
      integer(8) :: besselArray
!         MieCoeff1=Unit !debug
      MieCoeff1=(BesselJ(order, Unit*k*radius, besselArray) - MieCoeff2(order, radius, dielectric, besselArray)  &
          * Hankel1(order, Unit*k*radius, besselArray)) / (BesselJ(order, k*radius*sqrt(dielectric), besselArray))
      return
    end function MieCoeff1

    function MieCoeff2(order, radius, dielectric, besselArray)
      implicit none
      complex(8) :: MieCoeff2, dielectric, Hankel1prime, Hankel1, BesselJ, BesselJprime
      real(8) :: k=2d0*pi/lambda, radius
      real(8) :: order
      complex(8) value1
      integer(8) :: besselArray
        MieCoeff2= ( (sqrt(dielectric) * BesselJprime(order, k*radius*sqrt(dielectric), besselArray) &
                      * BesselJ(order, Unit*k*radius, besselArray) ) - (BesselJ(order,sqrt(dielectric)*k*radius, besselArray) &
                      *BesselJprime(order, Unit*k*radius, besselArray)) ) &
                    / (( sqrt(dielectric)*BesselJprime(order,k*radius*sqrt(dielectric), besselArray) &
                        *Hankel1(order, Unit*k*radius,besselArray) ) - ( BesselJ(order, sqrt(dielectric)*k*radius, besselArray) &
                        *Hankel1prime(order, Unit*k*radius, besselArray) )) !original

!         value1=radius*sqrt(dielectric)
!         MieCoeff2=BesselJprime(order, value1) !debug

!         MieCoeff2=BesselJprime(order, Unit*k*radius)
!          write(*,*) order, MieCoeff2
      return
    end function MieCoeff2

    function MieCoeff3(order, radius, dielectric, besselArray)
      implicit none
      complex(8) MieCoeff3, dielectric, Hankel1, BesselJ
      real(8) :: k=2d0*pi/lambda, radius
      real(8) :: order
      integer(8) :: besselArray
!         MieCoeff3=Unit !debug
      MieCoeff3=(BesselJ(order, Unit*k*radius, besselArray) - MieCoeff4(order, radius, dielectric, besselArray)  &
          * Hankel1(order, Unit*k*radius, besselArray)) / (BesselJ(order, k*radius*sqrt(dielectric), besselArray))
      return
    end function MieCoeff3

    function MieCoeff4(order, radius, dielectric, besselArray)
      implicit none
      complex(8) :: MieCoeff4, dielectric, Hankel1, Hankel1prime, BesselJ, BesselJprime
      real(8) :: k=2d0*pi/lambda, radius
      real(8) :: order
      integer(8) :: besselArray
      complex(8) value1
        MieCoeff4= ( ( BesselJprime(order, k*radius*sqrt(dielectric), besselArray)  &
                     * BesselJ(order, Unit*k*radius, besselArray) ) - sqrt(dielectric) * (BesselJ(order,sqrt(dielectric)*k*radius, besselArray) &
                     * BesselJprime(order, Unit*k*radius, besselArray)) ) &
                    / (( BesselJprime(order,k*radius*sqrt(dielectric), besselArray) * Hankel1(order, Unit*k*radius, besselArray) ) &
                    - sqrt(dielectric) * ( BesselJ(order, sqrt(dielectric)*k*radius, besselArray) &
                    * Hankel1prime(order, Unit*k*radius, besselArray) )) !original
!         MieCoeff4=Unit
      return
    end function MieCoeff4

end program Flaps

