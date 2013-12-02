! Thibault J.Y. Derrien
! Laboratoire Hubert Curien, UMR CNRS, St-Etienne
! ANR Ultrasonde
! Juillet - Decembre 2012

  !**** calculate temperature distribution in 2D in a tip ****


program Flaps

! include 'Bivariate.f'
! USE Bivariate

implicit none

    real(8), parameter:: lambda=515d-9 	, & !laser wavelength (m)
			fluence=10d0		, & !laser fluence (J.m-2)
			tau=40d-15		, & !FWHM pulse duration (s)
			spotX=50d-6		, & !FWHM spot size in X direction (1030nm: 400nm x 50nm ; 515nm: 50um x 50 um ; 343 nm: 50um x 100nm)
			spotY=50d-6		, & !FWHM spot size in Y direction
			xCenter=1000d-9			,& ! X position of the max of the intensity (1030nm: 1um 0um, 515nm: idem, 343nm: 100nm x 200nm)
			yCenter=0d0*200d-9			,&! Y position of the max of the intensity
			Tout=80d0	 		 ,&  !external temperature (K)
			potential0=7d3,& 	! potential at the bottom of the needle ; default = 7d3
			potentialNull=0d0, &
			phiMie=acos(-1d0)		! Mie scattering: polarization angle
    
    real(8), parameter:: dt0=10d-18          	,& !time step (s)
			tmax=12d-9	   	,& !stop time	
			xmin=-10d-6		,& !mesh min
			xmax=10d-6		,& !mesh max
			ymin=-10d-6		,&		
			ymax=10d-6		,&
			tCenter=0d0		,&	 !time of gaussian intensity maximum
			tmin=tCenter-5d0*tau 	   	 !max absolute time
    
			
    integer(8), parameter::  iterOut=10	,& ! number of iterations between each stdout
			  M=201	     	,& !number of cells main domain X direction
			  N=201		,& !number of cells main domain Y direection
			  Mv=101		,& !number of celles in the Vessel domain (larger) X direction
			  Nv=101	,& !number of celles in the Vessel domain (larger) Y direction
			  MeshChoice=1		,& !0: rectangle (xmin,xmax)(ymin,ymax). 1: cone, 2: cone in a vessel
			  MeshIterations=5000	,&	!number of iterations to calculate meshNeedle
			  MeshIterationsVessel=100*Mv,&	!number of iterations to calculate meshVessel
			  MeshShift=5		,& 	!number of cells x N in the tip, 343 nm: 2; 515 nm: 3;
			  FermiMaxLines=1112	,&	! >= number of lines in Fermi file
			  SORiterations=1, &	!iteration number for over-relaxation method
			  InterpolateMethod=1, &	! 0: linear, 1: bicubic
			  UseInterpolation=0, &
			  AdaptativeTimeStep=1, &
			  SolveImplicit=1	! 0: use explicit schemes, 1: use implicit scheme (band diagonal matrixes)
			  
			  
	
    real(8), parameter::   NeedleAngleDeg=10d0	,& ! deg
			  NeedleRadius=50d-9	,& !m
			  NeedleLength=3d-6	, &	!m
			  SORcoeff=1.2d0	,&	! near 1
			  MeshDensity0=5d5	,& ! amplitude of source for mesh refinement
			  MeshDamping0=1d8, &	 ! damping coefficient (m^-1) for mesh refinement
			  ActivateInduction=0d0, &
			  CrossCoeff=-1d0, &		! 0d0: OFF, 1d0: ON
			  maxCFL=1d-6			! maximum admitted on CFL condition for any time step increase
			  
    integer(8), parameter::  DrudeMode=0, &	! free-carrier absorption model, 0: classical, 1: using dispersion relation with plasma currents
			    ConductivityFix=0, &	! 2: Consider ambipolar diffusion in equations (but careful with boundary conditions)
						! 1: consider Tritt particle transport (great expression), but Dumber field is needed !!! -> Poisson ! 
						! 0: only fourier conductivity
						!-1: diffusion and conductivity OFF 
			    TeOff=0		,&		   !0: Disable temperature calculations
			    HolesOff=0		,&
			    TsOff=0		,&
			    ConvectionEnergy=0	,& 	!0: work with Te, no convection. 1: work with Ue, convection
			    DisableCrossDiffusion=0, &
			    PoissonOn=0		,& !0: Poisson solver is OFF. 1: Calculation of potential ON. 
			    DriftOn=0		,& !0: Drift is disabled. 1: Enabled. 
			    CathodeZone=1	,& !1: on the needle bottom, 0: on back vessel (not physical but stable)
			    ShiftFixedPotential=-1, &	! while CathodeZone=1, use to adjust the number of points on which tension is applied
			    PoissonSolver=0	,& !0: Full matrix inversion once, 1: SOR iterative for each dt
			    InterpolateOff=0, 	&	!just to test speedup...
			    BandBendingInFDTD=0	,&	!use the interpolation of FDTD 1030 nm with band-bending contribution
			    PolarizationSource=1, &	! 0: source Transervse electric, 1: transverse magnetic
			    UseMieScattering=0,& 		! 1: Enable Mie scattering analytic formula
			    maxBesselOrder=20,&		! Max of terms in series of Bessel for Mie scattering
			    besselArray=1
	
    real(8), parameter:: pi=3.14159265358d0 	,& 	!pi number
			  hbar=1.05457d-34   	,& 	!planck constant
			  epsilon0=8.85418781762d-12 ,& !vacuum dielectric permittivity
			  mu0=4d0*pi*1d-7		,& ! vacuum magnetic permeability
			  ec=1.60217646d-19        ,& 	!elementary charge
			  me0=9.10938188d-31        ,& 	!electron mass
			  c=2.99792458d8           ,& 	!light speed
			  kb=1.3806488d-23         ,&  	!Boltzmann constant, 
			  SiDensity=2.57d3	,&   	!Silicon rest density
			  epsilonStatic0=11.66570433d0 !,0.01404457712d0)		! dielectric constant for static field

    real(8), parameter:: omegaLaser=2d0*pi*c/lambda	,& !laser pulsation (s**-1)
			 me=0.5d0*me0				,& ! electron effective mass for conductivity !0.24 (source ?)
			 mh=0.5d0*me0				,&   ! hole effective mass for conductivity !0.81 (source ?)
			 meDOS=0.36d0*me0				,& ! electron effective mass for DOS
			 mhDOS=0.81d0*me0				,& ! hole effective mass for DOS
			 Ne0=1d5				, & 
			 Nh0=1d5, &				!initial density (to calculate auto using fermi
			 Nlimit=1d0				,&! lowest possible density
			 Nborder=1d23				,&! density on boundaries to consider defect layer
			 DefectThickness=1d-7			
    
    complex(8), parameter:: Imaginary=(0d0,1d0), Unit=(1d0,0d0), Zero=(0d0,0d0)		! complex unity
! 	             epsilonStatic0=(11.66570433d0,0.01404457712d0)		! dielectric constant for static field
    
    integer(8) 	nbiter, i, itwo, jtwo, j, k, l, nmax, NeedleIndexX, NeedleIndexY, maxFermiIndexE, maxFermiIndexH, &
		Mp, Np, imax, Nsolve, KLsolve, KUsolve, &
		Nsolve1, KLsolve1, KUsolve1, &
		Nsolve2, KLsolve2, KUsolve2, &
		Nsolve3, KLsolve3, KUsolve3, &
		Nsolve4, KLsolve4, KUsolve4, &
		Nsolve5, KLsolve5, KUsolve5
    
    logical	Diverged
    real(8) 	t, t0, dx, dy, x0, y0, dt
    real(8) 	Te0, Th0, Ts0, I0 !initial values of the problem
    real(8)	Ue(1:M, 1:N), & !electron energy
		Uh(1:M, 1:N), & !hole energy
		Te(1:M, 1:N), & !electron temperature
		Th(1:M, 1:N), & !hole temperature
		Ts(1:M, 1:N), & !lattice temperature
		Ne(1:M, 1:N), & !electron density
		Nh(1:M, 1:N), & !hole density
		UeNew(1:M, 1:N), & !electron energy
		UhNew(1:M, 1:N), & !hole energy
		TeNew(1:M, 1:N), & !electron temperature
		ThNew(1:M, 1:N), & !hole temperature
		TsNew(1:M, 1:N), & !lattice temperature
		NeNew(1:M, 1:N), & !electron density
		NhNew(1:M, 1:N), & !hole density
		GradNeX(1:M, 1:N),& !Grad(Ne)_x
		GradNeY(1:M, 1:N),& !Grad(Ne)_y
		TeDual(1:M-1,1:N-1), & !Te dual
		ThDual(1:M-1,1:N-1), & !Th dual
		TsDual(1:M-1,1:N-1), & !Ts dual
		NeDual(1:M-1,1:N-1), & !Ne dual
		NhDual(1:M-1,1:N-1), & !Nh dual
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
		CeOld(1:M, 1:N), ChOld(1:M,1:N), &
		CouplingE(1:M, 1:N), CouplingH(1:M, 1:N), &
		nuColl(1:M, 1:N) , &!	total collision frequency
		nuColleph(1:M, 1:N) , &!	electron-phonon collision frequency
		mobilityE(1:M, 1:N), mobilityH(1:M, 1:N), & 
		etae(1:M,1:N), etah(1:M,1:N), &	! reduced chemical Fermi potential
		Egap(1:M, 1:N), &			! local gap value
		SourceE(1:M, 1:N), SourceH(1:M, 1:N), & ! heating sources
		SourceUe(1:M, 1:N), SourceUh(1:M, 1:N), & ! free carrier thermal energy sources
		diffNe(1:M, 1:N), diffNh(1:M, 1:N), & 	! just for derivation in time
		x(1:M, 1:N), y(1:M, 1:N), & 		! needle position indexes
		xV(1:Mv, 1:Nv), yV(1:Mv, 1:Nv), & 		! vessel position indexes
		xDualSW(1:M, 1:N), yDualSW(1:M, 1:N), & 		! dual mesh position 
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
		Ex(1:M,1:N), Ey(1:M,1:N), &	! fields in the main domain
		potentialNeedle(1:M,1:N), &	! potential in the needle
		DummyNeedle(1:M, 1:N), &		! optional arguments for interpolation
		DummyDual(1:M-1, 1:N-1), &		! optional arguments for interpolation
		epsilonNeedle(1:M,1:N), &	! dielectric static in the needle
		MaxHeating(1:M, 1:N), &
		MaxHeatingTime(1:M, 1:N), &
		MeshDensity(1:Mv, 1:Nv)	,&		! function to adapt Poisson mesh on the needle mesh
		CellAreaN(1:M, 1:N), CellAreaS(1:M, 1:N), &			! area of the finite elements 
		CellAreaE(1:M, 1:N), CellAreaW(1:M, 1:N), &
		CellVol(1:M, 1:N), &				!volume of needle mesh cells
		NormalWx(1:M,1:N), NormalWy(1:M,1:N), &
		NormalEx(1:M,1:N), NormalEy(1:M,1:N), &		! normal to quadrangle elements
		NormalNx(1:M,1:N), NormalNy(1:M,1:N), &
		NormalSx(1:M,1:N), NormalSy(1:M,1:N), &
		TangentWx(1:M,1:N), TangentWy(1:M,1:N), &		! Tangent to quadrangle elements
		TangentEx(1:M,1:N), TangentEy(1:M,1:N), &		
		TangentNx(1:M,1:N), TangentNy(1:M,1:N), &
		TangentSx(1:M,1:N), TangentSy(1:M,1:N), &
		CurviWx(1:M,1:N), CurviWy(1:M,1:N), & 		! Unit vector between cell centers
		CurviEx(1:M,1:N), CurviEy(1:M,1:N), &		
		CurviNx(1:M,1:N), CurviNy(1:M,1:N), &
		CurviSx(1:M,1:N), CurviSy(1:M,1:N), &
		DistW(1:M,1:N), DistE(1:M,1:N), & 		! distance to the center of neighboor cells
		DistN(1:M,1:N), DistS(1:M,1:N), &
		DistDualW(1:M,1:N), DistDualE(1:M,1:N), & 		! distance of the element side (equal to area in 2D)
		DistDualN(1:M,1:N), DistDualS(1:M,1:N), &
		EintFieldR(1:M,1:N), EintFieldI(1:M, 1:N), &
		EintFieldDual(1:M-1, N-1)
		
    integer(8)	FermiIndexE(1:M,1:N), FermiIndexH(1:M,1:N), &
		SomeNeighbours(1:4,1:2) 				!Neighbours for the fixed potential
    
    
    real(8), allocatable, target :: FermiTableE(:,:) 			,&
				     FermiTableH(:,:) 			,& !reduced Fermi level for electrons and holes
				     Amatrix(:,:), Bvector(:), 	&
				     Xvector(:), XvectorPrev(:), 	&
! 				     AsolveP(:,:),
				     BsolveP(:), ABsolveP(:,:), XsolveP(:), &
				     BsolveP1(:), BsolveP2(:), BsolveP3(:), BsolveP4(:), BsolveP5(:), & 	!matrixes for implicit scheme
				     ABsolveP1(:,:), ABsolveP2(:,:), ABsolveP3(:,:), ABsolveP4(:,:), ABsolveP5(:,:), &
				     XsolvePrev(:), &
				     ErrorVec(:), &
				     spectralNorm(:), 			& !objects for matrix inversion calculation
				     xP(:,:), yP(:,:),			& 		! vessel position indexes
				     ExPoisson(:, :), EyPoisson(:, :), 		& ! electric field in the vessel
				     potential(:, :), 			&	! electric potential in the vessel
				     DielectricStatic(:, :), 		& ! dielectric constant for the static field
				     DummyVessel(:,:),			& ! option argument for interpolation
				     NeP(:,:), NhP(:,:)	, &		  !interpolated Ne,Nh in the vessel
				     CellVolume(:,:), &			! volume of the vessel/poisson mesh elements
				     NormalWxP(:,:), NormalWyP(:,:), &
				     NormalExP(:,:), NormalEyP(:,:), &		! normal to quadrangle elements
				     NormalNxP(:,:), NormalNyP(:,:), &
				     NormalSxP(:,:), NormalSyP(:,:), &
				     CellAreaNP(:,:), CellAreaSP(:,:), &
				     CellAreaEP(:,:), CellAreaWP(:,:)

    integer(8), allocatable, target:: FixedPotentialIndex(:,:), & 	!array of points where potential has been fixed
				      ipiv(:), ipiv1(:), ipiv2(:), ipiv3(:), ipiv4(:), ipiv5(:)
    
    complex(8) 	Dielectric(1:M,1:N), &! solid dielectric function under laser illumination
		DielectricDrudeE(1:M,1:N), & ! Drude part of dielectric function under laser illumination
		DielectricDrudeH(1:M,1:N), &
		EintField(1:M,1:N)	!Ez internal field for Mie scattering theory

		
    complex(8) epsilonInf !, SORsum !material constant
   
    real(8) AugerRateE, AugerRateH, &
	    sigmaTau, sigmaX, sigmaY, &
	    maxCFLxT, maxCFLyT, maxCFLxN, maxCFLyN, maxCFLxTs, maxCFLyTs, &
	    maxTe, minTe, maxTh, minTh, maxTs, minTs, maxIntensity, maxNe, minNe, maxNh, minNh, &
	    maxEnergy, maxSourceE, maxGainsE, maxSourceH, maxGainsH, maxGap, maxDiffNe, maxDiffNh, &
	    TotalLaserEnergy, TotalThermalEnergy, &
	    cpuefficiency, calc_time_begin, calc_time_1, calc_time_2, calc_time_3, &
	    NeedleHeight, NeedleA, NeedleB, Needlet0Limit, NeedleXParam, NeedleYParam, NeedleAngle, &
	    localT, P2critic, ConstBLx, ConstBLy, TotalNumOfE, TotalNumOfH, &
	    distX, distY, SORsum, xmin2, xmax2, ymin2, ymax2, &
	    ErrorSum
	    
    real(8) sigmaX1, sigmaX2, sigmaX3, sigmaX4, sigmaX5, sigmaX6, sigmaX7, sigmaX8, sigmaX9, &
	    sigmaY1, sigmaY2, sigmaY3, sigmaY4, sigmaY5, sigmaY6, sigmaY7, sigmaY8, sigmaY9, &
	    spotX1, spotX2, spotX3, spotX4, spotX5, spotX6, spotX7, spotX8, spotX9, &
	    spotY1, spotY2, spotY3, spotY4, spotY5, spotY6, spotY7, spotY8, spotY9, &
	    x1, x2, x3, x4, x5, x6, x7, x8, x9, &
	    y1, y2, y3, y4, y5, y6, y7, y8, y9, &
	    I1, I2, I3, I4, I5, I6, I7, I8, I9, &
	    periodX, periodY, &
	    NeTotal, NhTotal
	    
	    
    integer(4) unit1, unit2
    integer(8) ColFermiNeNc, ColFermiEta, ColFermi0, ColFermi1, ColFermi2, &
	       ColFermiHalf, ColFermiThreeHalf, ColFermiMenusHalf, info, &
	       info1, info2, info3, info4, info5
	    
    character(len=50)::format
!     integer(8) nthreads
!  !!********OpenMP Test*************
    integer :: myid, nthreads
    integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

    external dgbtrf, dgbtrs

!     CALL OMP_SET_NUM_THREADS(4);

  !$OMP PARALLEL default(none) private(myid) &
  !$OMP shared(nthreads)
  ! Determine the number of threads and their id
      myid = OMP_GET_THREAD_NUM()
      nthreads = OMP_GET_NUM_THREADS()
  !$OMP BARRIER
  
  if (myid==0) then 
    write(*,'(a)') 'OpenMP TEST'
    write(*,'(a,i1)') 'Number of Threads = ', nthreads
    write(*,'(a)') '*******************************'
  end if 
  !$OMP END PARALLEL
! !!******* END OpenMP test

  dt=dt0

  ColFermiNeNc=2; ColFermiEta=3; ColFermi0=4; ColFermi1=5; ColFermi2=6; ColFermiHalf=7; 
  ColFermiThreeHalf=8; ColFermiMenusHalf=9;

! test field


    AugerRateE=2.3d-43
    AugerRateH=7.8d-44

    sigmaTau=tau/(2d0*sqrt(2d0*log(2e0)))
    sigmaX=spotX/(2d0*sqrt(2d0*log(2e0)))
    sigmaY=spotY/(2d0*sqrt(2d0*log(2e0)))

    I0=fluence/tau
    
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
	x3=4.5d-7; y3=30d-9; I3=0d0*12d0*I0; spotX3=70d-9; spotY3=70d-9; !9.2 W unstable
	x4=6.2d-7; y4=40d-9; I4=0d0*28d0*I0; spotX4=70d-9; spotY4=120d-9; !36 W unstable
	x5=9.0d-7; y5=45d-9; I5=0d0*7d0*I0; spotX5=50d-9; spotY5=100d-9; !13.35 W unstable
	x6=1.1d-6; y6=60d-9; I6=0d0*15d0*I0; spotX6=50d-9; spotY6=70d-9; !8.24 W unstable
! 	x4=1.3d-6; y4=40d-9; I4=1d0*I0; spotX4=500d-9; spotY4=500d-9; !custom
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
    open(92,FILE='TimeBottom.dat', access='sequential', status='unknown') 	! format 882
    open(93,FILE='TimeUp.dat', access='sequential', status='unknown') 	! format 883
    open(94,FILE='TimeApex.dat', access='sequential', status='unknown') 	! format 884
    open(95,FILE='error.dat', access='sequential', status='unknown')
    open(96,FILE='parameters.dat', access='sequential', status='unknown')
    open(97,FILE='Depth.dat',access='sequential',status='unknown')		! format 887
        open(98,FILE='TimeMax.dat',access='sequential',status='unknown') 		! format 888
    open(99,FILE='mesh.dat',access='sequential',status='unknown')		! format 885
    open(100,FILE='meshVessel.dat',access='sequential',status='unknown')
    open(101,FILE='DepthVessel.dat',access='sequential',status='unknown') 	! format 889
    open(102,FILE='meshElements.dat', access='sequential',status='unknown') ! format 881
    open(103,FILE='DualDepth.dat', access='sequential', status='unknown') ! format 890
    ! FORMAT numbers already used for writing: 887, 886, 888, 885, 882, 883
    
  
  TotalLaserEnergy=0d0; TotalThermalEnergy=0d0
  cpuefficiency=0d0

!***************** MESH GENERATION *****************
  write(*,*) "[Mesh] Building..."
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
    
  if(MeshChoice>=1) then !conical mesh as main domain
    !boudary definition
    NeedleIndexX=M
    NeedleIndexY=N
    NeedleAngle=NeedleAngleDeg*pi/180d0
    NeedleA=NeedleRadius/(tan(NeedleAngle/2d0)**2)
    NeedleB=NeedleRadius/tan(NeedleAngle/2d0)
    Needlet0Limit=acos(NeedleRadius/(NeedleLength*tan(real(NeedleAngle)/2d0)**2+NeedleRadius))
    
!built of the conical mesh: write boundaries, then solve laplace, and iterate
    
    NeedleXParam=(Needlet0Limit-0d0)/NeedleIndexX !dt0 for X (0,t0)
    NeedleYParam=(Needlet0Limit-0d0)/NeedleIndexY !dt0 for Y (0,t0)
    
    P2critic=(real(N)-1d0)**2/(2d0*real(j)-real(N)-1)
    if(real(MeshShift)**2 > P2critic) then
      write(*,*) "Dilatation of the apex is too large. Reduce <= ", floor(sqrt(P2critic))
      stop
    end if
    
    do k=1,MeshIterations
    
!     !1/ tip apex, 2/ upper, 3/ bottom, 4/ cone base
	do j=1, N
	! t = (-p*dt0,p*dt0), p integer defined by MeshShift
	    localT=real(MeshShift)*real(MeshShift)*Needlet0Limit*( (2d0*real(j-1)) / (real(N)-1d0) - 1d0 )/real(N-1)
	    x(1,j)=NeedleA*(1d0/cos( localT ) - 1d0)
	    y(1,j)=NeedleB*tan(localT)
	    !x(1,j)=NeedleA*(1d0/cos(NeedleYParam*real( real(MeshShift)*((2d0*real(j-1))/(real(N)-1d0)-1d0)/(real(N)-1d0) )))	!tip cone boundary, 
	    !y(1,j)=NeedleB*tan(NeedleYParam*real( real(MeshShift)*((2d0*real(j-1))/(real(N)-1d0)-1d0)/(real(N)-1d0) ))		!tip cone boundary
	end do
! 	
	do i=1,M
	    localT=Needlet0Limit*( ( (real(i)-0d0)/(real(M)-1d0) )*(real(MeshShift)*real(MeshShift)/(real(N)-1d0)-1d0)+1d0)
	    x(i,1)=NeedleA*(1d0/cos( - localT )-1d0)
	    y(i,1)=NeedleB*tan( - localT)
	    ! x(i,1)=NeedleA*(1d0/cos(-NeedleXParam*real(i))) 	!upper boundary
	    ! y(i,1)=NeedleB*tan(-NeedleXParam*real(i)) 		!upper boundary
	end do
! ! !  
 	do i=1, M
	    !t = (-t0;-p*dt0), p integer	
	    localT=Needlet0Limit*( ( (real(i)-0d0)/(real(M)-1d0) ) * (real(MeshShift)*real(MeshShift)/(real(N)-1d0)-1d0)+1d0)
	    x(i,N)=NeedleA*(1d0/cos( localT )-1d0)
	    y(i,N)=NeedleB*tan( localT)
! 	    x(i,N)=NeedleA*(1d0/cos(localT)-1d0)		!bottom bounday,
! 	    y(i,N)=NeedleB*tan(NeedleXParam*real(i))			!bottom boundary
	end do

! ! ! 	NeedleXParam=(Needlet0Limit-0d0)/NeedleIndexX !dt0 for X (0,t0)
! ! ! 	NeedleYParam=(Needlet0Limit-0d0)/NeedleIndexY !dt0 for Y (0,t0)
	do j=1, N
	    localT=Needlet0Limit*((2d0*(real(j)-1d0)/(real(N)-1d0))-1d0)
	    x(M,j)=NeedleA*(1d0/cos(Needlet0Limit)-1d0)		!base cone boundary
	    y(M,j)=NeedleB*tan(localT) 	!base cone boundary
	end do
	
	! sort points before solving
	do j=1,M
	  do i=1,M
	    if(x(i,1) > x(j,1)) then 
	      call swap(x(i,1), x(j,1))
	      call swap(y(i,1), y(j,1))
	    end if
	    if(x(i,N) > x(i+1,N)) then 
	      call swap(x(i,N), x(j,N))
	      call swap(y(i,N), y(j,N))
	    end if
	  end do
	end do
	
	! resolution of laplace
	do i=2, M-1
	  do j=2,N-1
	    x(i,j)=(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1))/4d0
	    y(i,j)=(y(i+1,j)+y(i-1,j)+y(i,j+1)+y(i,j-1))/4d0
	  end do
	end do
	
    end do !end of iterations for mesh
    
    !writing of the mesh
    do i=1,M
      do j=1,N
	write(99, 885, advance='yes') x(i,j), y(i,j), i, j
	885	FORMAT (1E15.8, 3x, 1E15.8, 3x, I3, 3X, I3)
	!write(99,*)
	end do
	write(99,*) " "
    end do

    
!     do i=1,M-1
!       do j=1, N-1
! 	dx(i,j)=x(i+1,j)-x(i,j)
! 	dy(i,j)=y(i,j+1)-y(i,j)
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
! 	call swap(xV(i,1), xV(j,1))
! 	call swap(yV(i,1), yV(j,1))
!       end if
!       if(xV(i,Nv) > xV(i+1,Nv)) then 
! 	call swap(xV(i,Nv), xV(j,Nv))
! 	call swap(yV(i,Nv), yV(j,Nv))
!       end if
!     end do
!   end do
    
!   do k=1,MeshIterationsVessel
!     solve Laplace equations
    do i=2, Mv-1
      do j=2,Nv-1
! 	if(xV(i,j) < NeedleLength .AND. xV(i,j)>0d0 .AND. abs(yV(i,j)) < NeedleHeight/2d0) then !if this point is included in the needle
! 	  MeshDensity(i,j)=MeshDensity0* &
! 			exp(-0.5d0*( (xV(Mv/2,Nv/2)/(NeedleRadius/(2d0*sqrt(2d0*log(2.)))))**2 &
! 			+(yV(Mv/2,Nv/2)/(NeedleRadius/(2d0*sqrt(2d0*log(2.)))))**2 ) )
		    !*( & !lets define a mesh density function near from boundaries
! 		    exp(-xV(i,j)*MeshDamping0) * exp(-0.5d0*(yV(i,j))**2) + exp(-MeshDamping0*abs(yV(i,j)-NeedleHeight/2d0)))
! 	else 
! 	  MeshDensity(i,j)=0d0
! 	endif
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
! 	xV(i,j)=0d0
! 	yV(i,j)=0d0
!       end do
!     end do

!     do i=1,Mv
!       do j=1,Nv
! 	xV(i,j)=dx*real(i)+xmin
! 	yV(i,j)=dy*real(j)+ymin
!       end do
!     end do
    ! use to select the mesh in which Poisson eq is solved
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
!   allocate(AsolveP(1:Mp*Np,1:Mp*Np))
  allocate(BsolveP(1:Mp*Np));
  allocate(BsolveP1(1:Mp*Np)); allocate(BsolveP2(1:Mp*Np)); allocate(BsolveP3(1:Mp*Np)); allocate(BsolveP4(1:Mp*Np)); allocate(BsolveP5(1:Mp*Np))
  allocate(XsolveP(1:Mp*Np));
!   allocate(XsolveP1(1:Mp*Np)); allocate(XsolveP2(1:Mp*Np)); allocate(XsolveP3(1:Mp*Np)); allocate(XsolveP4(1:Mp*Np)); allocate(XsolveP5(1:Mp*Np))
  allocate(XsolvePrev(1:Mp*Np))
  allocate(ErrorVec(1:Mp*Np))

  allocate(ABsolveP(1:2*Np+Np+1,1:Mp*Np));
  allocate(ABsolveP1(1:2*Np+Np+1,1:Mp*Np)); allocate(ABsolveP2(1:2*Np+Np+1,1:Mp*Np)); allocate(ABsolveP3(1:2*Np+Np+1,1:Mp*Np)); allocate(ABsolveP4(1:2*Np+Np+1,1:Mp*Np)); allocate(ABsolveP5(1:2*Np+Np+1,1:Mp*Np))
  allocate(ipiv(1:Mp*Np));
  allocate(ipiv1(1:Mp*Np)); allocate(ipiv2(1:Mp*Np)); allocate(ipiv3(1:Mp*Np)); allocate(ipiv4(1:Mp*Np));allocate(ipiv5(1:Mp*Np))
  
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
  allocate(NormalNxP(1:Mp,1:Np), NormalNyP(1:Mp,:Np))
  allocate(NormalSxP(1:Mp,1:Np), NormalSyP(1:Mp,:Np))
  allocate(NormalExP(1:Mp,1:Np), NormalEyP(1:Mp,:Np))
  allocate(NormalWxP(1:Mp,1:Np), NormalWyP(1:Mp,:Np))
  allocate(CellAreaNP(1:Mp,1:Np), CellAreaSP(1:Mp,1:Np))
  allocate(CellAreaEP(1:Mp,1:Np), CellAreaWP(1:Mp,1:Np))
  
  
  call flush(99); call flush(100)
  
  ! test of neighboors
  write(*,*) "[TEST] Testing neighborhood of a internal point..."
  write(*,'(2I10.1)') transpose(FindNeighbours(x(2,2),y(2,2),xP,yP,4,Mp,Np))
  write(*,*) "[TEST] Testing neighborhood of a external point..."
  write(*,'(2I10.1)') transpose(FindNeighbours(xP(34,35),yP(34,35),x,y,4,M,N))
  
!   write(*,*) "[TEST] Area function"
!   write(*,*) AreaElement(5.99d-7,-2.5d-7, 7.68d-7,-2.85d-7, 6.11d-7,-2.11d-7,2.71d-7, -1.66d-7)
!   write(*,*) "[TEST] Triangle area function"
!   write(*,*) "5.99d-7,-2.5d-7 ; 7.68d-7,-2.85d-7 ; 6.11d-7,-2.11d-7 ; 2.71d-7, -1.66d-7"
!   write(*,*) "P 3.26d-7, 3.82d-8"
!   write(*,*) AreaTri(5.99d-7,-2.5d-7,7.68d-7,-2.85d-7,5.74d-7,-2.07d-7)
!   write(*,*) AreaTri(5.74d-7,-2.07d-7,7.68d-7,-2.85d-7,6.11d-7,-2.11d-7)
!   write(*,*) AreaTri(5.74d-7,-2.07d-7,6.11d-7,-2.11d-7, 2.71d-7, -1.66d-7)
!   write(*,*) AreaTri(5.74d-7,-2.07d-7,2.71d-7, -1.66d-7,5.99d-7,-2.5d-7)
!   write(*,*) "Sum=", AreaTri(5.99d-7,-2.5d-7,7.68d-7,-2.85d-7,5.74d-7,-2.07d-7) + &
!   AreaTri(5.74d-7,-2.07d-7,7.68d-7,-2.85d-7,6.11d-7,-2.11d-7) + &
!   AreaTri(5.74d-7,-2.07d-7,6.11d-7,-2.11d-7, 2.71d-7, -1.66d-7) + &
!   AreaTri(5.74d-7,-2.07d-7,2.71d-7, -1.66d-7,5.99d-7,-2.5d-7)
  
!*********** SOURCE IMPORT **************
  ! importing Fermi functions
   allocate(FermiTableE(1:9, 1:FermiMaxLines))
   allocate(FermiTableH(1:9, 1:FermiMaxLines))
   call TabCreateFL !(FermiTableE, FermiTableH)
!    FermiTableE(:,:)=1d0; FermiTableH(:,:)=1d0; ! uncomment if you want to disable fermi-dirac. Dont forget to lock the FermiIndexes also. 
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
  
  call flush(96)
  
  ! initialisation
  call cpu_time(calc_time_begin) !initialisation time
  
  nmax=int((tmax-tmin)/dt, 8)
  Ex(:,:)=0d0 !-1d10
  Ey(:,:)=0d0 !-1d9 !0d0
  DummyVessel(:,:)=0d0
!   AsolveP(:,:)=0d0;
  BsolveP(:)=0d0; XsolveP(:)=0d0
  Te0=Tout
  Th0=Tout
  
  !$OMP DO
  do i=1,M
    !$OMP DO
    do j=1,N
        
        VeX(i,j)=0d0
        VeY(i,j)=0d0
        VhX(i,j)=0d0
        VhY(i,j)=0d0
        TeNew(i,j)=Tout
        ThNew(i,j)=Tout
        TsNew(i,j)=Tout
        if(BandBendingInFDTD.eq.1) then
		  NeNew(i,j)=Ne0+Nborder*(exp(-0.5d0*(((x(i,j)-x(i,N))**2+(y(i,j)-y(i,N))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
		  +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
		  +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
		  )
		  NhNew(i,j)=Nh0+Nborder*(exp(-0.5d0*(((x(i,j)-x(i,N))**2+(y(i,j)-y(i,N))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
		  +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
		  +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
		  )
        else
	  NeNew(i,j)=Ne0
	  NhNew(i,j)=Nh0
	  Ne(i,j)=NeNew(i,j)
	  Nh(i,j)=NhNew(i,j)
        end if
        Ne(i,j)=NeNew(i,j)
	Nh(i,j)=NhNew(i,j)
        CeOld(i,j)=1.5d0*kb*NeNew(i,j) !approximation
        ChOld(i,j)=1.5d0*kb*NhNew(i,j) !approximation
        Ce(i,j)=CeOld(i,j)
        Ch(i,j)=ChOld(i,j)
        UeNew(i,j)=TeNew(i,j)*1.5d0*kb*NeNew(i,j) !approximation
        UhNew(i,j)=ThNew(i,j)*1.5d0*kb*NhNew(i,j) !approximation
        intensity(i,j)=0d0
        MaxHeating(i,j)=0d0
        MaxHeatingTime(i,j)=0d0
        epsilonNeedle(i,j)=epsilonStatic0-1d0
    end do
    !$OMP END DO
   end do
   !$OMP END DO
   if(CathodeZone.eq.1) then 
    !one need the indexes of the needle points where potential is imposed: matrix Nx2
    write(*,*) "Indexes of needle base on the vessel mesh"
    do j=1, N
      SomeNeighbours=FindNeighbours(x(M,j),y(M,j),xP,yP,4,Mp,Np) !OK
      FixedPotentialIndex(j,1:2)=SomeNeighbours(2,1:2) !select the second point for each element corresponding to needle base
      write(*,*) FixedPotentialIndex(j,:) !OK even if points are not very regular!
    end do
   end if
   
   potential(:,:)=0d0! (0d0,0d0)
   
   ! $ OMP PARALLEL DEFAULT (SHARED)
   write(*,*) 'Dielectric function interpolation...'
   if(InterpolateOff.eq.0) then 
    ! from needle mesh to vessel mesh
      if(InterpolateMethod.eq.0) then
	DielectricStatic=Interpolate(epsilonNeedle,x,y,xP,yP,M,N,Mp,Np)
      else 
! 	DielectricStatic=InterpolateBiCubic(epsilonNeedle, x, y, xP, yP, M, N, Mp, Np)
	call InterpolateBiCubic(epsilonNeedle, x, y, xP, yP, M, N, Mp, Np, DielectricStatic, DummyVessel, DummyVessel)
      end if
   endif
   ! $ OMP END PARALLEL
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
   do i=2,M-1
      do j=2,N-1
      
	DistN(i,j)=sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)
	DistS(i,j)=sqrt((x(i,j)-x(i,j-1))**2+(y(i,j)-y(i,j-1))**2)
	DistE(i,j)=sqrt((x(i+1,j)-x(i,j))**2+(y(i+1,j)-y(i,j))**2)
	DistW(i,j)=sqrt((x(i,j)-x(i-1,j))**2+(y(i,j)-y(i-1,j))**2)
	
	NormalNx(i,j)=Normal(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), 0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)),1) !4d0*(0.25d0*y(i-1,j)+0.25d0*y(i-1,j+1)-0.25d0*y(i+1,j+1)-0.25d0*y(i+1,j))/(x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)-2d0*x(i-1,j)*x(i+1,j+1)-2d0*x(i-1,j)*x(i+1,j)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i+1,j+1)-2d0*x(i-1,j+1)*x(i+1,j)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)+x(i+1,j)**2+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)-2d0*y(i-1,j)*y(i+1,j+1)-2d0*y(i-1,j)*y(i+1,j)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i+1,j+1)-2d0*y(i-1,j+1)*y(i+1,j)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)+y(i+1,j)**2)**(0.5d0)
	NormalNy(i,j)=Normal(0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), 0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)),2) !NormalNy(i,j)=-4d0*(0.25d0*x(i-1,j)+0.25d0*x(i-1,j+1)-0.25d0*x(i+1,j+1)-0.25d0*x(i+1,j))/(x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)-2d0*x(i-1,j)*x(i+1,j+1)-2d0*x(i-1,j)*x(i+1,j)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i+1,j+1)-2d0*x(i-1,j+1)*x(i+1,j)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)+x(i+1,j)**2+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)-2d0*y(i-1,j)*y(i+1,j+1)-2d0*y(i-1,j)*y(i+1,j)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i+1,j+1)-2d0*y(i-1,j+1)*y(i+1,j)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)+y(i+1,j)**2)**(0.5d0)
	NormalSx(i,j)=Normal(0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)),1) ! -4d0*(0.25d0*y(i-1,j-1)+0.25d0*y(i-1,j)-0.25d0*y(i+1,j)-0.25d0*y(i+1,j-1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)-2d0*x(i-1,j-1)*x(i+1,j)-2d0*x(i-1,j-1)*x(i+1,j-1)+x(i-1,j)**2-2d0*x(i-1,j)*x(i+1,j)-2d0*x(i-1,j)*x(i+1,j-1)+x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)+x(i+1,j-1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)-2d0*y(i-1,j-1)*y(i+1,j)-2d0*y(i-1,j-1)*y(i+1,j-1)+y(i-1,j)**2-2d0*y(i-1,j)*y(i+1,j)-2d0*y(i-1,j)*y(i+1,j-1)+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)+y(i+1,j-1)**2)**(0.5d0)
	NormalSy(i,j)=Normal(0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)),2) ! 4d0*(0.25d0*x(i-1,j-1)+0.25d0*x(i-1,j)-0.25d0*x(i+1,j)-0.25d0*x(i+1,j-1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)-2d0*x(i-1,j-1)*x(i+1,j)-2d0*x(i-1,j-1)*x(i+1,j-1)+x(i-1,j)**2-2d0*x(i-1,j)*x(i+1,j)-2d0*x(i-1,j)*x(i+1,j-1)+x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)+x(i+1,j-1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)-2d0*y(i-1,j-1)*y(i+1,j)-2d0*y(i-1,j-1)*y(i+1,j-1)+y(i-1,j)**2-2d0*y(i-1,j)*y(i+1,j)-2d0*y(i-1,j)*y(i+1,j-1)+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)+y(i+1,j-1)**2)**(0.5d0)
	NormalEx(i,j)=Normal(0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)), 0.25d0*(x(i+1,j+1)+x(i+1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i+1,j+1)+y(i+1,j)+y(i,j+1)+y(i,j)), 1) ! -4d0*(0.25d0*y(i+1,j-1)+0.25d0*y(i,j-1)-0.25d0*y(i+1,j+1)-0.25d0*y(i,j+1))/(x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)-2d0*x(i+1,j-1)*x(i+1,j+1)-2d0*x(i+1,j-1)*x(i,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i+1,j+1)-2d0*x(i,j-1)*x(i,j+1)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)+x(i,j+1)**2+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)-2d0*y(i+1,j-1)*y(i+1,j+1)-2d0*y(i+1,j-1)*y(i,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i+1,j+1)-2d0*y(i,j-1)*y(i,j+1)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)+y(i,j+1)**2)**(0.5d0)
	NormalEy(i,j)=Normal(0.25d0*(x(i+1,j-1)+x(i+1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i+1,j-1)+y(i+1,j)+y(i,j-1)+y(i,j)), 0.25d0*(x(i+1,j+1)+x(i+1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i+1,j+1)+y(i+1,j)+y(i,j+1)+y(i,j)), 2) ! 4d0*(0.25d0*x(i+1,j-1)+0.25d0*x(i,j-1)-0.25d0*x(i+1,j+1)-0.25d0*x(i,j+1))/(x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)-2d0*x(i+1,j-1)*x(i+1,j+1)-2d0*x(i+1,j-1)*x(i,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i+1,j+1)-2d0*x(i,j-1)*x(i,j+1)+x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)+x(i,j+1)**2+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)-2d0*y(i+1,j-1)*y(i+1,j+1)-2d0*y(i+1,j-1)*y(i,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i+1,j+1)-2d0*y(i,j-1)*y(i,j+1)+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)+y(i,j+1)**2)**(0.5d0)
	NormalWx(i,j)=Normal(0.25d0*(x(i-1,j+1)+x(i-1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i-1,j)+y(i,j+1)+y(i,j)), 0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 1) ! 4d0*(0.25d0*y(i-1,j-1)+0.25d0*y(i,j-1)-0.25d0*y(i,j+1)-0.25d0*y(i-1,j+1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)-2d0*x(i-1,j-1)*x(i,j+1)-2d0*x(i-1,j-1)*x(i-1,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i,j+1)-2d0*x(i,j-1)*x(i-1,j+1)+x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)+x(i-1,j+1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)-2d0*y(i-1,j-1)*y(i,j+1)-2d0*y(i-1,j-1)*y(i-1,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i,j+1)-2d0*y(i,j-1)*y(i-1,j+1)+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)
	NormalWy(i,j)=Normal(0.25d0*(x(i-1,j+1)+x(i-1,j)+x(i,j+1)+x(i,j)), 0.25d0*(y(i-1,j+1)+y(i-1,j)+y(i,j+1)+y(i,j)), 0.25d0*(x(i-1,j-1)+x(i-1,j)+x(i,j-1)+x(i,j)), 0.25d0*(y(i-1,j-1)+y(i-1,j)+y(i,j-1)+y(i,j)), 2)! -4d0*(0.25d0*x(i-1,j-1)+0.25d0*x(i,j-1)-0.25d0*x(i,j+1)-0.25d0*x(i-1,j+1))/(x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)-2d0*x(i-1,j-1)*x(i,j+1)-2d0*x(i-1,j-1)*x(i-1,j+1)+x(i,j-1)**2-2d0*x(i,j-1)*x(i,j+1)-2d0*x(i,j-1)*x(i-1,j+1)+x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)+x(i-1,j+1)**2+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)-2d0*y(i-1,j-1)*y(i,j+1)-2d0*y(i-1,j-1)*y(i-1,j+1)+y(i,j-1)**2-2d0*y(i,j-1)*y(i,j+1)-2d0*y(i,j-1)*y(i-1,j+1)+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)
	
	
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
	
	CellAreaN(i,j)=sqrt(( 0.25d0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))-0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j)))**2+(  0.25d0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))-0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j)))**2) !0.25d0*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)-2d0*x(i+1,j+1)*x(i-1,j)-2d0*x(i+1,j+1)*x(i-1,j+1)+x(i+1,j)**2-2d0*x(i+1,j)*x(i-1,j)-2d0*x(i+1,j)*x(i-1,j+1)+x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)+x(i-1,j+1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)-2d0*y(i+1,j+1)*y(i-1,j)-2d0*y(i+1,j+1)*y(i-1,j+1)+y(i+1,j)**2-2d0*y(i+1,j)*y(i-1,j)-2d0*y(i+1,j)*y(i-1,j+1)+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5d0)
	CellAreaS(i,j)=sqrt( (0.25d0*(x(i+1,j-1)+x(i,j-1)+x(i+1,j)+x(i,j))-0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i-1,j)+x(i,j)))**2 + (0.25d0*(y(i+1,j-1)+y(i,j-1)+y(i+1,j)+y(i,j))-0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i-1,j)+y(i,j)))**2 ) !0.25d0*(x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)-2d0*x(i+1,j)*x(i-1,j-1)-2d0*x(i+1,j)*x(i-1,j)+x(i+1,j-1)**2-2d0*x(i+1,j-1)*x(i-1,j-1)-2d0*x(i+1,j-1)*x(i-1,j)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)+x(i-1,j)**2+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)-2d0*y(i+1,j)*y(i-1,j-1)-2d0*y(i+1,j)*y(i-1,j)+y(i+1,j-1)**2-2d0*y(i+1,j-1)*y(i-1,j-1)-2d0*y(i+1,j-1)*y(i-1,j)+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i-1,j)+y(i-1,j)**2)**(0.5d0)
	CellAreaE(i,j)=sqrt( (0.25d0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))-0.25d0*(x(i+1,j-1)+x(i,j-1)+x(i+1,j)+x(i,j)))**2 + (0.25d0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))-0.25d0*(y(i+1,j-1)+y(i,j-1)+y(i+1,j)+y(i,j)))**2 ) !0.25d0*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)-2d0*x(i+1,j+1)*x(i+1,j-1)-2d0*x(i+1,j+1)*x(i,j-1)+x(i,j+1)**2-2d0*x(i,j+1)*x(i+1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)-2d0*y(i+1,j+1)*y(i+1,j-1)-2d0*y(i+1,j+1)*y(i,j-1)+y(i,j+1)**2-2d0*y(i,j+1)*y(i+1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)
	CellAreaW(i,j)=sqrt( (0.25d0*(x(i-1,j+1)+x(i,j+1)+x(i-1,j)+x(i,j))-0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i-1,j)+x(i,j)))**2 + (0.25d0*(y(i-1,j+1)+y(i,j+1)+y(i-1,j)+y(i,j))-0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i-1,j)+y(i,j)))**2 ) !0.25d0*(x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)-2d0*x(i,j+1)*x(i-1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i-1,j-1)-2d0*x(i-1,j+1)*x(i,j-1)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)-2d0*y(i,j+1)*y(i-1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i-1,j-1)-2d0*y(i-1,j+1)*y(i,j-1)+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)
	
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
   
   do i=2,Mp-1
    do j=2,Np-1
	NormalNxP(i,j)=4d0*(0.25d0*yP(i-1,j)+0.25d0*yP(i-1,j+1)-0.25d0*yP(i+1,j+1)-0.25d0*yP(i+1,j))/(xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)-2d0*xP(i-1,j)*xP(i+1,j+1)-2d0*xP(i-1,j)*xP(i+1,j)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i+1,j+1)-2d0*xP(i-1,j+1)*xP(i+1,j)+xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)+xP(i+1,j)**2+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)-2d0*yP(i-1,j)*yP(i+1,j+1)-2d0*yP(i-1,j)*yP(i+1,j)+yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i+1,j+1)-2d0*yP(i-1,j+1)*yP(i+1,j)+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)+yP(i+1,j)**2)**(0.5d0)
	NormalNyP(i,j)=-4d0*(0.25d0*xP(i-1,j)+0.25d0*xP(i-1,j+1)-0.25d0*xP(i+1,j+1)-0.25d0*xP(i+1,j))/(xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)-2d0*xP(i-1,j)*xP(i+1,j+1)-2d0*xP(i-1,j)*xP(i+1,j)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i+1,j+1)-2d0*xP(i-1,j+1)*xP(i+1,j)+xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)+xP(i+1,j)**2+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)-2d0*yP(i-1,j)*yP(i+1,j+1)-2d0*yP(i-1,j)*yP(i+1,j)+yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i+1,j+1)-2d0*yP(i-1,j+1)*yP(i+1,j)+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)+yP(i+1,j)**2)**(0.5d0)
	NormalSxP(i,j)=-4d0*(0.25d0*yP(i-1,j-1)+0.25d0*yP(i-1,j)-0.25d0*yP(i+1,j)-0.25d0*yP(i+1,j-1))/(xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)-2d0*xP(i-1,j-1)*xP(i+1,j)-2d0*xP(i-1,j-1)*xP(i+1,j-1)+xP(i-1,j)**2-2d0*xP(i-1,j)*xP(i+1,j)-2d0*xP(i-1,j)*xP(i+1,j-1)+xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)+xP(i+1,j-1)**2+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)-2d0*yP(i-1,j-1)*yP(i+1,j)-2d0*yP(i-1,j-1)*yP(i+1,j-1)+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i+1,j)-2d0*yP(i-1,j)*yP(i+1,j-1)+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)+yP(i+1,j-1)**2)**(0.5d0)
	NormalSyP(i,j)=4d0*(0.25d0*xP(i-1,j-1)+0.25d0*xP(i-1,j)-0.25d0*xP(i+1,j)-0.25d0*xP(i+1,j-1))/(xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)-2d0*xP(i-1,j-1)*xP(i+1,j)-2d0*xP(i-1,j-1)*xP(i+1,j-1)+xP(i-1,j)**2-2d0*xP(i-1,j)*xP(i+1,j)-2d0*xP(i-1,j)*xP(i+1,j-1)+xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)+xP(i+1,j-1)**2+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)-2d0*yP(i-1,j-1)*yP(i+1,j)-2d0*yP(i-1,j-1)*yP(i+1,j-1)+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i+1,j)-2d0*yP(i-1,j)*yP(i+1,j-1)+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)+yP(i+1,j-1)**2)**(0.5d0)
	NormalExP(i,j)=-4d0*(0.25d0*yP(i+1,j-1)+0.25d0*yP(i,j-1)-0.25d0*yP(i+1,j+1)-0.25d0*yP(i,j+1))/(xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)-2d0*xP(i+1,j-1)*xP(i+1,j+1)-2d0*xP(i+1,j-1)*xP(i,j+1)+xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i+1,j+1)-2d0*xP(i,j-1)*xP(i,j+1)+xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)+xP(i,j+1)**2+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)-2d0*yP(i+1,j-1)*yP(i+1,j+1)-2d0*yP(i+1,j-1)*yP(i,j+1)+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i+1,j+1)-2d0*yP(i,j-1)*yP(i,j+1)+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)+yP(i,j+1)**2)**(0.5d0)
	NormalEyP(i,j)=4d0*(0.25d0*xP(i+1,j-1)+0.25d0*xP(i,j-1)-0.25d0*xP(i+1,j+1)-0.25d0*xP(i,j+1))/(xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)-2d0*xP(i+1,j-1)*xP(i+1,j+1)-2d0*xP(i+1,j-1)*xP(i,j+1)+xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i+1,j+1)-2d0*xP(i,j-1)*xP(i,j+1)+xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)+xP(i,j+1)**2+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)-2d0*yP(i+1,j-1)*yP(i+1,j+1)-2d0*yP(i+1,j-1)*yP(i,j+1)+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i+1,j+1)-2d0*yP(i,j-1)*yP(i,j+1)+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)+yP(i,j+1)**2)**(0.5d0)
	NormalWxP(i,j)=4d0*(0.25d0*yP(i-1,j-1)+0.25d0*yP(i,j-1)-0.25d0*yP(i,j+1)-0.25d0*yP(i-1,j+1))/(xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)-2d0*xP(i-1,j-1)*xP(i,j+1)-2d0*xP(i-1,j-1)*xP(i-1,j+1)+xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i,j+1)-2d0*xP(i,j-1)*xP(i-1,j+1)+xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)-2d0*yP(i-1,j-1)*yP(i,j+1)-2d0*yP(i-1,j-1)*yP(i-1,j+1)+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j+1)-2d0*yP(i,j-1)*yP(i-1,j+1)+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)
	NormalWyP(i,j)=-4d0*(0.25d0*xP(i-1,j-1)+0.25d0*xP(i,j-1)-0.25d0*xP(i,j+1)-0.25d0*xP(i-1,j+1))/(xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)-2d0*xP(i-1,j-1)*xP(i,j+1)-2d0*xP(i-1,j-1)*xP(i-1,j+1)+xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i,j+1)-2d0*xP(i,j-1)*xP(i-1,j+1)+xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)-2d0*yP(i-1,j-1)*yP(i,j+1)-2d0*yP(i-1,j-1)*yP(i-1,j+1)+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j+1)-2d0*yP(i,j-1)*yP(i-1,j+1)+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)
	
	CellAreaNP(i,j)=0.25d0*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)
	CellAreaSP(i,j)=0.25d0*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)-2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5d0)
	CellAreaEP(i,j)=0.25d0*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)
	CellAreaWP(i,j)=0.25d0*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)
    end do
  end do
  
      do i=2,M-1
	
	 
	! NORTH
	
	 CellVol(i,N)=0.25d0*AreaElement(x(i-1,N),y(i-1,N),x(i+1,N),y(i+1,N),x(i+1,N-1),y(i+1,N-1),x(i-1,N-1),y(i-1,N-1))
	 CellAreaN(i,N)=0d0 !0.5d0*(x(i+1,N)**2-2d0*x(i+1,N)*x(i-1,N)+x(i-1,N)**2+y(i+1,N)**2-2d0*y(i+1,N)*y(i-1,N)+y(i-1,N)**2)**0.5d0
	 CellAreaS(i,N)=0.25d0*(x(i+1,N)**2+2d0*x(i+1,N)*x(i+1,N-1)-2d0*x(i+1,N)*x(i-1,N-1)-2d0*x(i+1,N)*x(i-1,N)+x(i+1,N-1)**2-2d0*x(i+1,N-1)*x(i-1,N-1)-2d0*x(i+1,N-1)*x(i-1,N)+x(i-1,N-1)**2+2d0*x(i-1,N)*x(i-1,N-1)+x(i-1,N)**2+y(i+1,N)**2+2d0*y(i+1,N)*y(i+1,N-1)-2d0*y(i+1,N)*y(i-1,N-1)-2d0*y(i+1,N)*y(i-1,N)+y(i+1,N-1)**2-2d0*y(i+1,N-1)*y(i-1,N-1)-2d0*y(i+1,N-1)*y(i-1,N)+y(i-1,N-1)**2+2d0*y(i-1,N)*y(i-1,N-1)+y(i-1,N)**2)**0.5d0
	 CellAreaW(i,N)=0.25d0*(x(i,N)**2+2d0*x(i,N)*x(i-1,N)-2d0*x(i,N)*x(i-1,N-1)-2d0*x(i,N)*x(i,N-1)+x(i-1,N)**2-2d0*x(i-1,N)*x(i-1,N-1)-2d0*x(i-1,N)*x(i,N-1)+x(i-1,N-1)**2+2d0*x(i-1,N-1)*x(i,N-1)+x(i,N-1)**2+y(i,N)**2+2d0*y(i,N)*y(i-1,N)-2d0*y(i,N)*y(i-1,N-1)-2d0*y(i,N)*y(i,N-1)+y(i-1,N)**2-2d0*y(i-1,N)*y(i-1,N-1)-2d0*y(i-1,N)*y(i,N-1)+y(i-1,N-1)**2+2d0*y(i-1,N-1)*y(i,N-1)+y(i,N-1)**2)**0.5d0
	 CellAreaE(i,N)=0.25d0*(x(i+1,N)**2+2d0*x(i+1,N)*x(i,N)-2d0*x(i+1,N)*x(i+1,N-1)-2d0*x(i+1,N)*x(i,N-1)+x(i,N)**2-2d0*x(i,N)*x(i+1,N-1)-2d0*x(i,N)*x(i,N-1)+x(i+1,N-1)**2+2d0*x(i+1,N-1)*x(i,N-1)+x(i,N-1)**2+y(i,N)**2+2d0*y(i,N)*y(i+1,N)-2d0*y(i,N)*y(i+1,N-1)-2d0*y(i,N)*y(i,N-1)+y(i+1,N)**2-2d0*y(i+1,N)*y(i+1,N-1)-2d0*y(i+1,N)*y(i,N-1)+y(i+1,N-1)**2+2d0*y(i+1,N-1)*y(i,N-1)+y(i,N-1)**2)**0.5d0
	 
	 NormalEx(i,N)=-Normal(0.5d0*(x(i,N)+x(i+1,N)), 0.5d0*(y(i,N)+y(i+1,N)), 0.25d0*(x(i+1,N)+x(i,N)+x(i+1,N-1)+x(i,N-1)), 0.25d0*(y(i+1,N)+y(i,N)+y(i+1,N-1)+y(i,N-1)), 1)
	 NormalEy(i,N)=-Normal(0.5d0*(x(i,N)+x(i+1,N)), 0.5d0*(y(i,N)+y(i+1,N)), 0.25d0*(x(i+1,N)+x(i,N)+x(i+1,N-1)+x(i,N-1)), 0.25d0*(y(i+1,N)+y(i,N)+y(i+1,N-1)+y(i,N-1)), 2)
	 NormalWx(i,N)=Normal(0.5d0*(x(i-1,N)+x(i,N)), 0.5d0*(y(i-1,N)+y(i,N)), 0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), 0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),1); 
	 NormalWy(i,N)=Normal(0.5d0*(x(i-1,N)+x(i,N)), 0.5d0*(y(i-1,N)+y(i,N)), 0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)), 0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),2); 
	 NormalNx(i,N)=Normal(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),1)
	 NormalNy(i,N)=Normal(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),2)
	 NormalSx(i,N)=-Normal(0.25d0*(x(i+1,N-1)+x(i+1,N)+x(i,N)+x(i,N-1)), 0.25d0*(y(i+1,N-1)+y(i+1,N)+y(i,N)+y(i,N-1)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)),0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),1) 
	 NormalSy(i,N)=-Normal(0.25d0*(x(i+1,N-1)+x(i+1,N)+x(i,N)+x(i,N-1)), 0.25d0*(y(i+1,N-1)+y(i+1,N)+y(i,N)+y(i,N-1)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N-1)+x(i,N)),0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N-1)+y(i,N)),2) 
	 
	 TangentWx(i,N)=-Tangent(0.5d0*(x(i-1,N)+x(i,N)),0.5d0*(y(i-1,N)+y(i,N)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N)+x(i,N-1)),0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N)+y(i,N-1)),1)
	 TangentWy(i,N)=-Tangent(0.5d0*(x(i-1,N)+x(i,N)),0.5d0*(y(i-1,N)+y(i,N)),0.25d0*(x(i-1,N-1)+x(i-1,N)+x(i,N)+x(i,N-1)),0.25d0*(y(i-1,N-1)+y(i-1,N)+y(i,N)+y(i,N-1)),2)
	 TangentEx(i,N)=-Tangent(0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)), 0.5d0*(x(i+1,N)+x(i,N)), 0.5d0*(y(i+1,N)+y(i,N)), 1)
	 TangentEy(i,N)=-Tangent(0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)), 0.5d0*(x(i+1,N)+x(i,N)), 0.5d0*(y(i+1,N)+y(i,N)), 2)
	 TangentNx(i,N)=-Tangent(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),1)
	 TangentNy(i,N)=-Tangent(0.5d0*(x(i,N)+x(i+1,N)),0.5d0*(y(i,N)+y(i+1,N)), 0.5d0*(x(i,N)+x(i-1,N)), 0.5d0*(y(i,N)+y(i-1,N)),2)
	 TangentSx(i,N)=-Tangent(0.25d0*(x(i-1,N-1)+x(i,N-1)+x(i-1,N)+x(i,N)), 0.25d0*(y(i-1,N-1)+y(i,N-1)+y(i-1,N)+y(i,N)), 0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)),1)
	 TangentSy(i,N)=-Tangent(0.25d0*(x(i-1,N-1)+x(i,N-1)+x(i-1,N)+x(i,N)), 0.25d0*(y(i-1,N-1)+y(i,N-1)+y(i-1,N)+y(i,N)), 0.25d0*(x(i+1,N-1)+x(i,N-1)+x(i+1,N)+x(i,N)), 0.25d0*(y(i+1,N-1)+y(i,N-1)+y(i+1,N)+y(i,N)),2)
	 
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
	 CellAreaN(i,1)=0.25d0*(x(i+1,2)**2+2d0*x(i+1,2)*x(i+1,1)-2d0*x(i+1,2)*x(i-1,1)-2d0*x(i+1,2)*x(i-1,2)+x(i+1,1)**2-2d0*x(i+1,1)*x(i-1,1)-2d0*x(i+1,1)*x(i-1,2)+x(i-1,1)**2+2d0*x(i-1,1)*x(i-1,2)+x(i-1,2)**2+y(i+1,2)**2+2d0*y(i+1,2)*y(i+1,1)-2d0*y(i+1,2)*y(i-1,1)-2d0*y(i+1,2)*y(i-1,2)+y(i+1,1)**2-2d0*y(i+1,1)*y(i-1,1)-2d0*y(i+1,1)*y(i-1,2)+y(i-1,1)**2+2d0*y(i-1,1)*y(i-1,2)+y(i-1,2)**2)**(0.5d0)
	 CellAreaS(i,1)=0d0 !sqrt((0.5d0*(x(i+1,1)+x(i,1))-0.5d0*(x(i-1,1)+x(i,1)))**2+(0.5d0*(y(i+1,1)+y(i,1))-0.5d0*(y(i-1,1)+y(i,1)))**2)
	 CellAreaW(i,1)=0.25d0*(x(i,1)**2-2d0*x(i,1)*x(i,2)+2d0*x(i,1)*x(i-1,1)-2d0*x(i,1)*x(i-1,2)+x(i,2)**2-2d0*x(i,2)*x(i-1,1)+2d0*x(i,2)*x(i-1,2)+x(i-1,1)**2-2d0*x(i-1,1)*x(i-1,2)+x(i-1,2)**2+y(i,1)**2-2d0*y(i,1)*y(i,2)+2d0*y(i,1)*y(i-1,1)-2d0*y(i,1)*y(i-1,2)+y(i,2)**2-2d0*y(i,2)*y(i-1,1)+2d0*y(i,2)*y(i-1,2)+y(i-1,1)**2-2d0*y(i-1,1)*y(i-1,2)+y(i-1,2)**2)**(0.5d0)
	 CellAreaE(i,1)=0.25d0*(x(i+1,2)**2-2d0*x(i+1,2)*x(i+1,1)-2d0*x(i+1,2)*x(i,1)+2d0*x(i+1,2)*x(i,2)+x(i+1,1)**2+2d0*x(i+1,1)*x(i,1)-2d0*x(i+1,1)*x(i,2)+x(i,1)**2-2d0*x(i,1)*x(i,2)+x(i,2)**2+y(i+1,2)**2-2d0*y(i+1,2)*y(i+1,1)-2d0*y(i+1,2)*y(i,1)+2d0*y(i+1,2)*y(i,2)+y(i+1,1)**2+2d0*y(i+1,1)*y(i,1)-2d0*y(i+1,1)*y(i,2)+y(i,1)**2-2d0*y(i,1)*y(i,2)+y(i,2)**2)**(0.5d0)
	 
	 NormalEx(i,1)=-Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)),0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),1) 
	 NormalEy(i,1)=-Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)),0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),2)
	 NormalWx(i,1)=Normal(0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),1)
	 NormalWy(i,1)=Normal(0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),2)
	 NormalNx(i,1)=Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)),0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),1)
	 NormalNy(i,1)=Normal(0.25d0*(x(i,1)+x(i,2)+x(i+1,1)+x(i+1,2)),0.25d0*(y(i,1)+y(i,2)+y(i+1,1)+y(i+1,2)),0.25d0*(x(i-1,1)+x(i-1,2)+x(i,1)+x(i,2)),0.25d0*(y(i-1,1)+y(i-1,2)+y(i,1)+y(i,2)),2)
	 NormalSx(i,1)=Normal(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),1) 
	 NormalSy(i,1)=Normal(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),2) 
	 
	 TangentNx(i,1)=Tangent(0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),1)
	 TangentNy(i,1)=Tangent(0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),2)
	 TangentSx(i,1)=-Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),1)
	 TangentSy(i,1)=-Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.5d0*(x(i,1)+x(i+1,1)),0.5d0*(y(i,1)+y(i+1,1)),2)
	 TangentEx(i,1)=-Tangent(0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),1)
	 TangentEy(i,1)=-Tangent(0.5d0*(x(i+1,1)+x(i,1)),0.5d0*(y(i+1,1)+y(i,1)),0.25d0*(x(i+1,2)+x(i+1,1)+x(i,2)+x(i,1)),0.25d0*(y(i+1,2)+y(i+1,1)+y(i,2)+y(i,1)),2)
	 TangentWx(i,1)=Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),1)
	 TangentWy(i,1)=Tangent(0.5d0*(x(i-1,1)+x(i,1)),0.5d0*(y(i-1,1)+y(i,1)),0.25d0*(x(i-1,2)+x(i-1,1)+x(i,2)+x(i,1)),0.25d0*(y(i-1,2)+y(i-1,1)+y(i,2)+y(i,1)),2)
	 
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
	  CellAreaN(1,j)=0.25d0*(x(2,j+1)**2+2d0*x(2,j+1)*x(2,j)-2d0*x(2,j+1)*x(1,j)-2d0*x(2,j+1)*x(1,j+1)+x(2,j)**2-2d0*x(2,j)*x(1,j)-2d0*x(2,j)*x(1,j+1)+x(1,j)**2+2d0*x(1,j)*x(1,j+1)+x(1,j+1)**2+y(2,j+1)**2+2d0*y(2,j+1)*y(2,j)-2d0*y(2,j+1)*y(1,j)-2d0*y(2,j+1)*y(1,j+1)+y(2,j)**2-2d0*y(2,j)*y(1,j)-2d0*y(2,j)*y(1,j+1)+y(1,j)**2+2d0*y(1,j)*y(1,j+1)+y(1,j+1)**2)**0.5d0
	  CellAreaS(1,j)=0.25d0*(x(2,j)**2+2d0*x(2,j)*x(2,j-1)-2d0*x(2,j)*x(1,j)-2d0*x(2,j)*x(1,j-1)+x(2,j-1)**2-2d0*x(2,j-1)*x(1,j)-2d0*x(2,j-1)*x(1,j-1)+x(1,j)**2+2d0*x(1,j)*x(1,j-1)+x(1,j-1)**2+y(2,j)**2+2d0*y(2,j)*y(2,j-1)-2d0*y(2,j)*y(1,j)-2d0*y(2,j)*y(1,j-1)+y(2,j-1)**2-2d0*y(2,j-1)*y(1,j)-2d0*y(2,j-1)*y(1,j-1)+y(1,j)**2+2d0*y(1,j)*y(1,j-1)+y(1,j-1)**2)**0.5d0
	  CellAreaW(1,j)=sqrt( (0.5d0*(x(1,j+1)+x(1,j))-0.5d0*(x(1,j-1)+x(1,j)))**2 + (0.5d0*(y(1,j+1)+y(1,j))-0.5d0*(y(1,j-1)+y(1,j)))**2 )
	  CellAreaE(1,j)=0.25d0*(x(2,j+1)**2+2d0*x(2,j+1)*x(1,j+1)-2d0*x(2,j+1)*x(2,j-1)-2d0*x(2,j+1)*x(1,j-1)+x(1,j+1)**2-2d0*x(1,j+1)*x(2,j-1)-2d0*x(1,j+1)*x(1,j-1)+x(2,j-1)**2+2d0*x(2,j-1)*x(1,j-1)+x(1,j-1)**2+y(2,j+1)**2+2d0*y(2,j+1)*y(1,j+1)-2d0*y(2,j+1)*y(2,j-1)-2d0*y(2,j+1)*y(1,j-1)+y(1,j+1)**2-2d0*y(1,j+1)*y(2,j-1)-2d0*y(1,j+1)*y(1,j-1)+y(2,j-1)**2+2d0*y(2,j-1)*y(1,j-1)+y(1,j-1)**2)**0.5d0
	  
	  NormalEx(1,j)=-Normal(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)),0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),1)
	  NormalEy(1,j)=-Normal(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)),0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),2)
	  NormalNx(1,j)=Normal(0.25d0*(x(2,j+1)+x(1,j+1)+x(2,j)+x(1,j)),0.25d0*(y(2,j+1)+y(1,j+1)+y(2,j)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),1)
	  NormalNy(1,j)=Normal(0.25d0*(x(2,j+1)+x(1,j+1)+x(2,j)+x(1,j)),0.25d0*(y(2,j+1)+y(1,j+1)+y(2,j)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),2)
	  NormalSx(1,j)=-Normal(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),1) 
	  NormalSy(1,j)=-Normal(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),2) 
	  NormalWx(1,j)=Normal(0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),1)
	  NormalWy(1,j)=Normal(0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),2)
	  
! 	  write(*,*) "[Debug]", NormalWx(1,j)**2+NormalWy(1,j)**2
	  
	  TangentNx(1,j)=-Tangent(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),1)
	  TangentNy(1,j)=-Tangent(0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),2)
	  TangentSx(1,j)=-Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.25d0*(x(1,j-1)+x(1,j)+x(2,j-1)+x(2,j)),0.25d0*(y(1,j-1)+y(1,j)+y(2,j-1)+y(2,j)),1)
	  TangentSy(1,j)=-Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.25d0*(x(1,j-1)+x(1,j)+x(2,j-1)+x(2,j)),0.25d0*(y(1,j-1)+y(1,j)+y(2,j-1)+y(2,j)),2)
	  TangentWx(1,j)=Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),1)
	  TangentWy(1,j)=Tangent(0.5d0*(x(1,j-1)+x(1,j)),0.5d0*(y(1,j-1)+y(1,j)),0.5d0*(x(1,j+1)+x(1,j)),0.5d0*(y(1,j+1)+y(1,j)),2)
	  TangentEx(1,j)=-Tangent(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)),1)
	  TangentEy(1,j)=-Tangent(0.25d0*(x(2,j-1)+x(2,j)+x(1,j-1)+x(1,j)),0.25d0*(y(2,j-1)+y(2,j)+y(1,j-1)+y(1,j)),0.25d0*(x(2,j+1)+x(2,j)+x(1,j+1)+x(1,j)),0.25d0*(y(2,j+1)+y(2,j)+y(1,j+1)+y(1,j)),2)
	  
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
	  
	  CellVol(M,j)=0.25d0*AreaElement(x(M-1,j-1),y(M-1,j-1),x(M,j-1),y(M,j-1),x(M,j+1),y(M,j+1),x(M-1,j+1),y(M-1,j+1))
	  CellAreaN(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M,j)-2d0*x(M,j+1)*x(M-1,j)-2d0*x(M,j+1)*x(M-1,j+1)+x(M,j)**2-2d0*x(M,j)*x(M-1,j)-2d0*x(M,j)*x(M-1,j+1)+x(M-1,j)**2+2d0*x(M-1,j)*x(M-1,j+1)+x(M-1,j+1)**2+y(M,j+1)**2+2d0*y(M,j+1)*y(M,j)-2d0*y(M,j+1)*y(M-1,j)-2d0*y(M,j+1)*y(M-1,j+1)+y(M,j)**2-2d0*y(M,j)*y(M-1,j)-2d0*y(M,j)*y(M-1,j+1)+y(M-1,j)**2+2d0*y(M-1,j)*y(M-1,j+1)+y(M-1,j+1)**2)**0.5d0
	  CellAreaS(M,j)=0.25d0*(x(M,j)**2+2d0*x(M,j)*x(M,j-1)-2d0*x(M,j)*x(M-1,j-1)-2d0*x(M,j)*x(M-1,j)+x(M,j-1)**2-2d0*x(M-1,j-1)*x(M,j-1)-2d0*x(M,j-1)*x(M-1,j)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M-1,j)+x(M-1,j)**2+y(M,j)**2+2d0*y(M,j)*y(M,j-1)-2d0*y(M,j)*y(M-1,j-1)-2d0*y(M,j)*y(M-1,j)+y(M,j-1)**2-2d0*y(M-1,j-1)*y(M,j-1)-2d0*y(M,j-1)*y(M-1,j)+y(M-1,j-1)**2+2d0*y(M-1,j-1)*y(M-1,j)+y(M-1,j)**2)**0.5d0
	  CellAreaW(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M-1,j+1)-2d0*x(M,j+1)*x(M-1,j-1)-2d0*x(M,j+1)*x(M,j-1)+x(M-1,j+1)**2-2d0*x(M-1,j+1)*x(M-1,j-1)-2d0*x(M-1,j+1)*x(M,j-1)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M,j-1)+x(M,j-1)**2+y(M,j+1)**2+2d0*y(M,j+1)*y(M-1,j+1)-2d0*y(M,j+1)*y(M-1,j-1)-2d0*y(M,j+1)*y(M,j-1)+y(M-1,j+1)**2-2d0*y(M-1,j+1)*y(M-1,j-1)-2d0*y(M-1,j+1)*y(M,j-1)+y(M-1,j-1)**2+2d0*y(M-1,j-1)*y(M,j-1)+y(M,j-1)**2)**0.5d0
	  CellAreaE(M,j)=sqrt( (0.5d0*(x(M,j+1)+x(M,j))-0.5d0*(x(M,j-1)+x(M,j)) )**2 + (0.5d0*(y(M,j+1)+y(M,j))-0.5d0*(y(M,j-1)+y(M,j)) )**2)
	  
	  TangentNx(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1)
	  TangentNy(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2)
	  TangentEx(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),1)
	  TangentEy(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),2)
	  TangentSx(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),1)
	  TangentSy(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),2)
	  TangentWx(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1)
	  TangentWy(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),2)
	  
! 	  if(y(1,j)<0d0) then
! 	    TangentNx(M,j)=-TangentNx(M,j)
! 	    TangentNy(M,j)=-TangentNy(M,j)
! 	    TangentSx(M,j)=-TangentSx(M,j)
! 	    TangentSy(M,j)=-TangentSy(M,j)
! 	    TangentWx(M,j)=-TangentWx(M,j)
! 	    TangentWy(M,j)=-TangentWy(M,j)
! 	  end if
	  
	  NormalEx(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),1); 
	  NormalEy(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),2); 
	  NormalWx(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1) 
	  NormalWy(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),2)
	  NormalNx(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1) 
	  NormalNy(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2) 
	  NormalSx(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1)+x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),1) 
	  NormalSy(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1)+x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),2) 
	  
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
	  CellAreaN(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M,j)-2d0*x(M,j+1)*x(M-1,j)-2d0*x(M,j+1)*x(M-1,j+1)+x(M,j)**2-2d0*x(M,j)*x(M-1,j)-2d0*x(M,j)*x(M-1,j+1)+x(M-1,j)**2+2d0*x(M-1,j)*x(M-1,j+1)+x(M-1,j+1)**2+y(M,j+1)**2+2d0*y(M,j+1)*y(M,j)-2d0*y(M,j+1)*y(M-1,j)-2d0*y(M,j+1)*y(M-1,j+1)+y(M,j)**2-2d0*y(M,j)*y(M-1,j)-2d0*y(M,j)*y(M-1,j+1)+y(M-1,j)**2+2d0*y(M-1,j)*y(M-1,j+1)+y(M-1,j+1)**2)**0.5d0
	  CellAreaS(M,j)=0.25d0*(x(M,j)**2+2d0*x(M,j)*x(M,j-1)-2d0*x(M,j)*x(M-1,j-1)-2d0*x(M,j)*x(M-1,j)+x(M,j-1)**2-2d0*x(M-1,j-1)*x(M,j-1)-2d0*x(M,j-1)*x(M-1,j)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M-1,j)+x(M-1,j)**2+y(M,j)**2+2d0*y(M,j)*y(M,j-1)-2d0*y(M,j)*y(M-1,j-1)-2d0*y(M,j)*y(M-1,j)+y(M,j-1)**2-2d0*y(M-1,j-1)*y(M,j-1)-2d0*y(M,j-1)*y(M-1,j)+y(M-1,j-1)**2+2d0*y(M-1,j-1)*y(M-1,j)+y(M-1,j)**2)**0.5d0
	  CellAreaW(M,j)=0.25d0*(x(M,j+1)**2+2d0*x(M,j+1)*x(M-1,j+1)-2d0*x(M,j+1)*x(M-1,j-1)-2d0*x(M,j+1)*x(M,j-1)+x(M-1,j+1)**2-2d0*x(M-1,j+1)*x(M-1,j-1)-2d0*x(M-1,j+1)*x(M,j-1)+x(M-1,j-1)**2+2d0*x(M-1,j-1)*x(M,j-1)+x(M,j-1)**2+y(M,j+1)**2+2d0*y(M,j+1)*y(M-1,j+1)-2d0*y(M,j+1)*y(M-1,j-1)-2d0*y(M,j+1)*y(M,j-1)+y(M-1,j+1)**2-2d0*y(M-1,j+1)*y(M-1,j-1)-2d0*y(M-1,j+1)*y(M,j-1)+y(M-1,j-1)**2+2d0*y(M-1,j-1)*y(M,j-1)+y(M,j-1)**2)**0.5d0
	  CellAreaE(M,j)=sqrt( (0.5d0*(x(M,j+1)+x(M,j))-0.5d0*(x(M,j-1)+x(M,j)) )**2 + (0.5d0*(y(M,j+1)+y(M,j))-0.5d0*(y(M,j-1)+y(M,j)) )**2)
	  
	  TangentNx(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1)
	  TangentNy(M,j)=-Tangent(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2)
	  TangentEx(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),1)
	  TangentEy(M,j)=-Tangent(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),2)
	  TangentSx(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),1)
	  TangentSy(M,j)=-Tangent(0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),2)
	  TangentWx(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1)
	  TangentWy(M,j)=-Tangent(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),2)
	  
! 	  if(y(1,j)<0d0) then
! 	    TangentNx(M,j)=-TangentNx(M,j)
! 	    TangentNy(M,j)=-TangentNy(M,j)
! 	    TangentSx(M,j)=-TangentSx(M,j)
! 	    TangentSy(M,j)=-TangentSy(M,j)
! 	    TangentWx(M,j)=-TangentWx(M,j)
! 	    TangentWy(M,j)=-TangentWy(M,j)
! 	  end if
	  
	  NormalEx(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),1); 
	  NormalEy(M,j)=-Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),2); 
	  NormalWx(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),1) 
	  NormalWy(M,j)=Normal(0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M-1,j)+x(M,j-1)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M-1,j)+y(M,j-1)+y(M,j)),2)
	  NormalNx(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),1) 
	  NormalNy(M,j)=Normal(0.5d0*(x(M,j+1)+x(M,j)),0.5d0*(y(M,j+1)+y(M,j)),0.25d0*(x(M-1,j+1)+x(M-1,j)+x(M,j+1)+x(M,j)),0.25d0*(y(M-1,j+1)+y(M-1,j)+y(M,j+1)+y(M,j)),2) 
	  NormalSx(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1)+x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),1) 
	  NormalSy(M,j)=-Normal(0.5d0*(x(M,j-1)+x(M,j)),0.5d0*(y(M,j-1)+y(M,j)),0.25d0*(x(M-1,j-1)+x(M,j-1)+x(M-1,j)+x(M,j)),0.25d0*(y(M-1,j-1)+y(M,j-1)+y(M-1,j)+y(M,j)),2) 
	  
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
	CellVol(M,N)=AreaElement(0.25d0*(x(M-1, N-1)+x(M, N-1)+x(M, N)+x(M-1, N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N)+y(M,N-1)), 0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)), x(M,N),y(M,N), 0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)))
! 	CellVol(M,N)=0.25d0*AreaElement(x(M-1,N-1),y(M-1,N-1),x(M,N-1),y(M,N-1), x(M,N), y(M,N), x(M-1, N), y(M-1, N))
	CellAreaE(M,N)=sqrt((x(M,N)-0.5d0*(x(M,N)+x(M,N-1)))**2+(y(M,N)-0.5d0*(y(M,N)+y(M,N-1)))**2)
	CellAreaN(M,N)=sqrt((x(M,N)-0.5d0*(x(M-1,N)+x(M,N)))**2+(y(M,N)-0.5d0*(y(M-1,N)+y(M,N)))**2)
	CellAreaS(M,N)=0.25d0*(x(M,N)**2+2d0*x(M,N)*x(M,N-1)-2d0*x(M,N)*x(M-1,N-1)-2d0*x(M-1,N)*x(M,N)+x(M,N-1)**2-2d0*x(M-1,N-1)*x(M,N-1)-2d0*x(M-1,N)*x(M,N-1)+x(M-1,N-1)**2+2d0*x(M-1,N)*x(M-1,N-1)+x(M-1,N)**2+y(M,N)**2+2d0*y(M,N)*y(M,N-1)-2d0*y(M,N)*y(M-1,N-1)-2d0*y(M-1,N)*y(M,N)+y(M,N-1)**2-2d0*y(M-1,N-1)*y(M,N-1)-2d0*y(M-1,N)*y(M,N-1)+y(M-1,N-1)**2+2d0*y(M-1,N)*y(M-1,N-1)+y(M-1,N)**2)**(0.5d0)
	CellAreaW(M,N)=0.25d0*(x(M-1,N)**2+2d0*x(M-1,N)*x(M,N)-2d0*x(M-1,N)*x(M-1,N-1)-2d0*x(M-1,N)*x(M,N-1)+x(M,N)**2-2d0*x(M,N)*x(M-1,N-1)-2d0*x(M,N)*x(M,N-1)+x(M-1,N-1)**2+2d0*x(M-1,N-1)*x(M,N-1)+x(M,N-1)**2+y(M-1,N)**2+2d0*y(M-1,N)*y(M,N)-2d0*y(M-1,N)*y(M-1,N-1)-2d0*y(M-1,N)*y(M,N-1)+y(M,N)**2-2d0*y(M,N)*y(M-1,N-1)-2d0*y(M,N)*y(M,N-1)+y(M-1,N-1)**2+2d0*y(M-1,N-1)*y(M,N-1)+y(M,N-1)**2)**(0.5d0)
	NormalEx(M,N)=-Normal(x(M,N),y(M,N),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),1)
	NormalEy(M,N)=-Normal(x(M,N),y(M,N),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),2)
	NormalNx(M,N)=Normal(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),1)
	NormalNy(M,N)=Normal(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),2)
	NormalWx(M,N)=Normal(0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)),0.25d0*(x(M-1, N-1)+x(M, N-1)+x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),1) 
	NormalWy(M,N)=Normal(0.5d0*(x(M-1, N)+x(M, N)),0.5d0*(y(M-1, N)+y(M, N)),0.25d0*(x(M-1, N-1)+x(M, N-1)+x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),2)
	NormalSx(M,N)=-Normal(0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)),0.25d0*(x(M-1, N-1)+x(M, N-1)+x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),1)
	NormalSy(M,N)=-Normal(0.5d0*(x(M, N)+x(M, N-1)),0.5d0*(y(M, N)+y(M, N-1)),0.25d0*(x(M-1, N-1)+x(M, N-1)+x(M, N)+x(M-1, N)),0.25d0*(y(M-1, N-1)+y(M, N-1)+y(M, N)+y(M-1, N)),2)
	
	TangentNx(M,N)=-Tangent(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),1)
	TangentNy(M,N)=-Tangent(x(M,N),y(M,N),0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),2)
	TangentEx(M,N)=-Tangent(0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),x(M,N),y(M,N),1)
	TangentEy(M,N)=-Tangent(0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),x(M,N),y(M,N),2)
	TangentWx(M,N)=-Tangent(0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),1)
	TangentWy(M,N)=-Tangent(0.5d0*(x(M-1,N)+x(M,N)),0.5d0*(y(M-1,N)+y(M,N)),0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),2)
	TangentSx(M,N)=-Tangent(0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),1)
	TangentSy(M,N)=-Tangent(0.25d0*(x(M-1,N-1)+x(M-1,N)+x(M,N-1)+x(M,N)),0.25d0*(y(M-1,N-1)+y(M-1,N)+y(M,N-1)+y(M,N)),0.5d0*(x(M,N-1)+x(M,N)),0.5d0*(y(M,N-1)+y(M,N)),2)
	
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
	CellVol(M,1)=AreaElement(0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)), x(M,1),y(M,1), 0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)), 0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)))
! 	CellVol(M,1)=0.25d0*AreaElement(x(M-1,1),y(M-1,1),x(M,1),y(M,1),x(M,2),y(M,2),x(M-1,2),y(M-1,2))
	CellAreaS(M,1)=sqrt(( x(M,1)-0.5d0*(x(M,1)+x(M-1,1)) )**2+( y(M,1)-0.5d0*(y(M,1)+y(M-1,1)) )**2)
	CellAreaE(M,1)=sqrt( ( x(M,1) - 0.5d0*(x(M,1)+x(M,2)) )**2 + ( y(M,1) - 0.5d0*(y(M,1)+y(M,2)) )**2 )
	CellAreaN(M,1)=0.25d0*(x(M,1)**2+2d0*x(M,1)*x(M,2)-2d0*x(M-1,1)*x(M,1)-2d0*x(M-1,2)*x(M,1)+x(M,2)**2-2d0*x(M-1,1)*x(M,2)-2d0*x(M-1,2)*x(M,2)+x(M-1,1)**2+2d0*x(M-1,1)*x(M-1,2)+x(M-1,2)**2+y(M,1)**2+2d0*y(M,1)*y(M,2)-2d0*y(M-1,1)*y(M,1)-2d0*y(M-1,2)*y(M,1)+y(M,2)**2-2d0*y(M-1,1)*y(M,2)-2d0*y(M-1,2)*y(M,2)+y(M-1,1)**2+2d0*y(M-1,1)*y(M-1,2)+y(M-1,2)**2)**(0.5d0)
	CellAreaW(M,1)=0.25d0*(x(M-1,1)**2-2d0*x(M-1,1)*x(M-1,2)+2d0*x(M-1,1)*x(M,1)-2d0*x(M-1,1)*x(M,2)+x(M-1,2)**2-2d0*x(M-1,2)*x(M,1)+2d0*x(M-1,2)*x(M,2)+x(M,1)**2-2d0*x(M,1)*x(M,2)+x(M,2)**2+y(M-1,1)**2-2d0*y(M-1,1)*y(M-1,2)+2d0*y(M-1,1)*y(M,1)-2d0*y(M-1,1)*y(M,2)+y(M-1,2)**2-2d0*y(M-1,2)*y(M,1)+2d0*y(M-1,2)*y(M,2)+y(M,1)**2-2d0*y(M,1)*y(M,2)+y(M,2)**2)**(0.5d0)
	NormalSx(M,1)=Normal(0.5d0*(x(M-1,1)+x(M,1)),0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 1)
	NormalSy(M,1)=Normal(0.5d0*(x(M-1,1)+x(M,1)),0.5d0*(y(M-1,1)+y(M,1)), x(M,1), y(M,1), 2)
	NormalEx(M,1)=Normal(x(M,1), y(M,1), 0.5d0*(x(M,1)+x(M,2)), 0.5d0*(y(M,1)+y(M,2)), 1)
	NormalEy(M,1)=Normal(x(M,1), y(M,1), 0.5d0*(x(M,1)+x(M,2)), 0.5d0*(y(M,1)+y(M,2)), 2)
	NormalWx(M,1)=Normal(0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)),1); 
	NormalWy(M,1)=Normal(0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),0.5d0*(x(M-1, 1)+x(M, 1)),0.5d0*(y(M-1, 1)+y(M, 1)),2);
	NormalNx(M,1)=Normal(0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)),0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),1); 
	NormalNy(M,1)=Normal(0.5d0*(x(M, 1)+x(M, 2)),0.5d0*(y(M, 1)+y(M, 2)),0.25d0*(x(M-1, 1)+x(M-1, 2)+x(M, 1)+x(M, 2)),0.25d0*(y(M-1, 1)+y(M-1, 2)+y(M, 1)+y(M, 2)),2);
	
	TangentNx(M,1)=-Tangent(0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),0.25d0*(x(M,1)+x(M,2)+x(M-1,1)+x(M-1,2)),0.25d0*(y(M,1)+y(M,2)+y(M-1,1)+y(M-1,2)),1)
	TangentNy(M,1)=-Tangent(0.5d0*(x(M,1)+x(M,2)),0.5d0*(y(M,1)+y(M,2)),0.25d0*(x(M,1)+x(M,2)+x(M-1,1)+x(M-1,2)),0.25d0*(y(M,1)+y(M,2)+y(M-1,1)+y(M-1,2)),2)
	TangentWx(M,1)=-Tangent(0.25d0*(x(M-1,2)+x(M-1,1)+x(M,2)+x(M,1)),0.25d0*(y(M-1,2)+y(M-1,1)+y(M,2)+y(M,1)),0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)),1)
	TangentWy(M,1)=-Tangent(0.25d0*(x(M-1,2)+x(M-1,1)+x(M,2)+x(M,1)),0.25d0*(y(M-1,2)+y(M-1,1)+y(M,2)+y(M,1)),0.5d0*(x(M-1,1)+x(M,1)), 0.5d0*(y(M-1,1)+y(M,1)),2)
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
	CellVol(1,1)=AreaElement(x(1,1),y(1,1),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)),0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2)+y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)))
! 	CellVol(1,1)=0.25d0*AreaElement(x(1,1),y(1,1),x(2,1),y(2,1),x(2,2),y(2,2),x(1,2),y(1,2))
	
! 	CellVol(1,1)=2d0*CellVol(1,1)

	CellAreaS(1,1)=sqrt((x(1,1)-0.5d0*(x(2,1)+x(1,1)))**2+(y(1,1)-0.5d0*(y(2,1)+y(1,1)))**2) 
	CellAreaW(1,1)=sqrt((x(1,1)-0.5d0*(x(1,2)+x(1,1)))**2+(y(1,1)-0.5d0*(y(1,2)+y(1,1)))**2)
	CellAreaN(1,1)=sqrt((0.25d0*(x(2,2)+x(2,1)+x(1,2)+x(1,1))-0.5d0*(x(1,1)+x(1,2)))**2+(0.25d0*(y(2,2)+y(2,1)+y(1,2)+y(1,1))-0.5d0*(y(1,2)+y(1,1)))**2) !sqrt((0.5d0*(x(1,2)+x(1,1))-0.25d0*(x(2,2)+x(2,1)+x(1,1)+x(1,2)))**2d0+(0.5d0*(y(1,2)+y(1,1))-0.25d0*(y(2,2)+y(2,1)+y(1,1)+y(1,2)))**2d0) !0.25d0*(x(1,1)**2+2d0*x(1,1)*x(1,2)-2d0*x(1,1)*x(2,1)-2d0*x(1,1)*x(2,2)+x(1,2)**2-2d0*x(1,2)*x(2,1)-2d0*x(1,2)*x(2,2)+x(2,1)**2+2d0*x(2,1)*x(2,2)+x(2,2)**2+y(1,1)**2+2d0*y(1,1)*y(1,2)-2d0*y(1,1)*y(2,1)-2d0*y(1,1)*y(2,2)+y(1,2)**2-2d0*y(1,2)*y(2,1)-2d0*y(1,2)*y(2,2)+y(2,1)**2+2d0*y(2,1)*y(2,2)+y(2,2)**2)**(0.5d0);
	CellAreaE(1,1)=sqrt((0.25d0*(x(2,2)+x(2,1)+x(1,2)+x(1,1))-0.5d0*(x(2,1)+x(1,1)))**2+(0.25d0*(y(2,2)+y(2,1)+y(1,2)+y(1,1))-0.5d0*(y(2,1)+y(1,1)))**2) !0.25d0*(x(1,1)**2-2d0*x(1,1)*x(1,2)+2d0*x(1,1)*x(2,1)-2d0*x(1,1)*x(2,2)+x(1,2)**2-2d0*x(1,2)*x(2,1)+2d0*x(1,2)*x(2,2)+x(2,1)**2-2d0*x(2,1)*x(2,2)+x(2,2)**2+y(1,1)**2-2d0*y(1,1)*y(1,2)+2d0*y(1,1)*y(2,1)-2d0*y(1,1)*y(2,2)+y(1,2)**2-2d0*y(1,2)*y(2,1)+2d0*y(1,2)*y(2,2)+y(2,1)**2-2d0*y(2,1)*y(2,2)+y(2,2)**2)**(0.5d0)
	
	NormalSx(1,1)=Normal(x(1,1),y(1,1),0.5d0*(x(2,1)+x(1,1)), 0.5d0*(y(2,1)+y(1,1)), 1)
	NormalSy(1,1)=Normal(x(1,1),y(1,1),0.5d0*(x(2,1)+x(1,1)), 0.5d0*(y(2,1)+y(1,1)), 2)
	NormalWx(1,1)=-Normal(x(1,1),y(1,1),0.5d0*(x(1,2)+x(1,1)), 0.5d0*(y(1,2)+y(1,1)), 1)
	NormalWy(1,1)=-Normal(x(1,1),y(1,1),0.5d0*(x(1,2)+x(1,1)), 0.5d0*(y(1,2)+y(1,1)), 2)
	
	NormalEx(1,1)=-Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2)+y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)),1); 
	NormalEy(1,1)=-Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2)+y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(2, 1)),0.5d0*(y(1, 1)+y(2, 1)),2); 
	NormalNx(1,1)=Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2)+y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)),1); 
	NormalNy(1,1)=Normal(0.25d0*(x(1, 1)+x(1, 2)+x(2, 1)+x(2, 2)),0.25d0*(y(1, 1)+y(1, 2)+y(2, 1)+y(2, 2)),0.5d0*(x(1, 1)+x(1, 2)),0.5d0*(y(1, 1)+y(1, 2)),2);
	
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
! 	CellVol(1,N)=0.25d0*AreaElement(x(1,N-1),y(1,N-1),x(2,N-1),y(2,N-1),x(2,N),y(2,N),x(1,N),y(1,N))
	
! 	CellVol(1,N)=2d0*CellVol(1,N)
	
	CellAreaN(1,N)=sqrt((x(1,N)-0.5d0*(x(1,N)+x(2,N)))**2+(y(1,N)-0.5d0*(y(1,N)+y(2,N)))**2)
	CellAreaW(1,N)=sqrt((x(1,N)-0.5d0*(x(1,N)+x(1,N-1)))**2+((y(1,N)-0.5d0*(y(1,N)+y(1,N-1))))**2)
	CellAreaS(1,N)=sqrt((0.25d0*(x(2,N-1)+x(1,N-1)+x(2,N)+x(1,N))-0.5d0*(x(1,N-1)+x(1,N)))**2+(0.25d0*(y(2,N-1)+y(1,N-1)+y(2,N)+y(1,N))-0.5d0*(y(1,N-1)+y(1,N)))**2) !sqrt((0.5d0*(x(1,N-1)+x(1,N))-0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)))**2+(0.5d0*(y(1,N-1)+y(1,N))-0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)))**2)
	CellAreaE(1,N)=sqrt((0.5d0*(x(2,N)+x(1,N))-0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)))**2+(0.5d0*(y(2,N)+y(1,N))-0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)))**2) !sqrt((0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1))-(0.5d0*(x(2,N)+x(1,N))))**2+(0.25d0*(y(1,N)+y(1,N-1)+y(2,N)+y(2,N-1))-0.5d0*(y(2,N)+y(1,N)))**2)
	
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
	TangentEx(1,N)=-Tangent(0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1)),0.25d0*(y(1,N)+y(1,N-1)+y(2,N)+y(2,N-1)),0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),1)
	TangentEy(1,N)=-Tangent(0.25d0*(x(1,N)+x(1,N-1)+x(2,N)+x(2,N-1)),0.25d0*(y(1,N)+y(1,N-1)+y(2,N)+y(2,N-1)),0.5d0*(x(2,N)+x(1,N)),0.5d0*(y(2,N)+y(1,N)),2)
	TangentSx(1,N)=-Tangent(0.5d0*(x(1,N-1)+x(1,N)),0.5d0*(y(1,N-1)+y(1,N)),0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)),0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)),1)
	TangentSy(1,N)=-Tangent(0.5d0*(x(1,N-1)+x(1,N)),0.5d0*(y(1,N-1)+y(1,N)),0.25d0*(x(2,N-1)+x(2,N)+x(1,N-1)+x(1,N)),0.25d0*(y(2,N-1)+y(2,N)+y(1,N-1)+y(1,N)),2)
	
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
  
   
   !! defining material index
   epsilonInf=DielectricConstant(lambda)
  write(*,*) 'epsilon(', 1d9*lambda, 'nm)=', epsilonInf
if(UseMieScattering.eq.1) then
  write(*,*) 'Computing the Mie scattering field distribution...'
  !$OMP DO
  do i=1,M
    !$OMP DO
    do j=1,N
      EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie, abs(ContourYofX(x(i,j), NeedleRadius, NeedleAngle)), epsilonInf) * sqrt(2d0*fluence/(c*epsilon0*tau))
    end do
    !$OMP END DO
  end do
  !$OMP END DO
  EintField=EintField*conjg(EintField) !change to modulus of the Ez field
  write(*,*) 'Smoothing the obtained intensity profile to reduce mesh size...'
  EintFieldR=EintField
!   EintFieldI=aimag(EintField) !should be equal to 0
!! SMOOTHING
!   !$OMP DO
!   do i=1,M-1
!     !$OMP DO
!     do j=1,N-1
!       0d0=( &
! 		  + 0.5d0*CellAreaE(i,j)*(Ne(i+1,j)-Ne(i,j))/DistE(i,j) * (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j)) &
! 		  - 0.5d0*CellAreaW(i,j)*(Ne(i,j)-Ne(i-1,j))/DistW(i,j) * (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j)) &
! 		  + 0.5d0*CellAreaN(i,j)*(Ne(i,j+1)-Ne(i,j))/DistN(i,j) * (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j)) &
! 		  - 0.5d0*CellAreaS(i,j)*(Ne(i,j)-Ne(i,j-1))/DistS(i,j) * (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j)) &
! ! 		  ! Cross-diffusion from [Mathur and Murthy (1997)]
! 		  + CrossCoeff*( &
! 		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*( NeDual(i,j) - NeDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
! 		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*( NeDual(i-1,j-1) - NeDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*( NeDual(i-1,j) - NeDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
! 		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*( NeDual(i,j-1) - NeDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
! 		  )
!     end do
!     !$OMP END DO
!   end do
!   !$OMP END DO
  
!   call InterpolateBiCubic(EintFieldR, x, y, xDual, yDual, M, N, M-1, N-1, EintFieldDual, DummyNeedle, DummyNeedle)
!   call InterpolateBiCubic(EintFieldDual, xDual, yDual, x, y, M-1, N-1, M, N, EintFieldI, DummyNeedle, DummyNeedle)
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
! 	  Amatrix(OneDindex(i,j),OneDindex(i,j))=1d0 ! Diag(A)=1d0 for Dirichlet conditions !OK
!       end do
!     end do
!     
!     do i=2, Mp-1 !pour chaque point ou l on va calculer le potentiel
! 	do j=2, Np-1 
! 	  Amatrix(OneDindex(i,j),OneDindex(i,j))=(1d0/8d0)*(-DielectricStatic(i+1,j)-DielectricStatic(i,j))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i,j)+xP(i,j)**2+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)-(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i-1,j))*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2 & 
! 	  -2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i-1,j)**2-2d0*xP(i-1,j)*xP(i,j)+xP(i,j)**2+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)+(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i,j+1))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)/(xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i,j)+xP(i,j)**2+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i,j)+yP(i,j)**2)**(0.5d0)-(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i,j-1))*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)- &
! 	  2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5d0)/(xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i,j)+xP(i,j)**2+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! 	 
! 	  Amatrix(OneDindex(i,j),OneDindex(i,j-1))=-(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i,j-1))*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)-2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5d0)/(xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i,j)+xP(i,j)**2+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! 	  ! 0.25d0*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)-2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5) !a(i,j-1)
! 	  Amatrix(OneDindex(i,j),OneDindex(i,j+1))=(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i,j+1))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)/(xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i,j)+xP(i,j)**2+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! 	  ! 0.25d0*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5) !a(i,j+1)
! 	  Amatrix(OneDindex(i,j),OneDindex(i-1,j))=-(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i-1,j))*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i-1,j)**2-2d0*xP(i-1,j)*xP(i,j)+xP(i,j)**2+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! 	  !0.25d0*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5) !a(i-1,j)
! 	  Amatrix(OneDindex(i,j),OneDindex(i+1,j))=(1d0/8d0)*(DielectricStatic(i+1,j)+DielectricStatic(i,j))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i,j)+xP(i,j)**2+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! ! 	  0.25d0*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5) !a(i+1,j)
! 
!       end do
!     end do !well filled
!     
!     if(CathodeZone.eq.1) then
!       !diagonal equals to 1
!       do k=FixedPotentialIndex(1,2)-ShiftFixedPotential, FixedPotentialIndex(N,2)+1+ShiftFixedPotential ! filling the indexes of fixed potential with 1d0, not localized on boundaries
!       	do j=FixedPotentialIndex(1,2)-ShiftFixedPotential, FixedPotentialIndex(N,2)+1+ShiftFixedPotential
! 	  !other points have to be zero
! 	  Amatrix(OneDindex(FixedPotentialIndex(1,1), k),&
! 	  OneDindex(FixedPotentialIndex(1,1), j))=0d0
! 	end do
! 	Amatrix(OneDindex(FixedPotentialIndex(1,1), k),&
! 		OneDindex(FixedPotentialIndex(1,1), k))=1d0
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
! ! 	write(*,*) '[Gauss inversion] Singular matrix!'
! !       end if
! !       do i=1,Mp*Np
! ! 	call swap(Amatrix(i,k), Amatrix(i,imax))
! !       end do
! !       do i=k+1, Mp*Np
! ! 	do j=k+1,Mp*Np
! ! 	  Amatrix(i,j)=Amatrix(i,j)-Amatrix(k,j)*(Amatrix(i,k)/Amatrix(k,k))
! ! 	end do
! ! 	Amatrix(i,k)=0d0
! !       end do
! !     end do
!   
!   ! now use this matrix at each time step (carefull: it increases round off errors)
!   end if
  
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
		   "maxCFL_Ts=",maxCFLxTs+maxCFLyTs, "CPU=", cpuefficiency, "NumThreads=", nthreads
      write(*,*) "NeTot=", NeTotal, " NhTotal=", NhTotal
    end if
    
    maxIntensity=0d0; maxTe=0d0; minTe=1d10; maxTh=0d0; minTh=1d10; maxTs=0d0; minTs=1d10; 
    maxNe=0d0; minNe=1d50; maxNh=0d0; minNh=1d50; maxCFLxT=0d0; maxCFLyT=0d0; maxCFLxN=0d0; maxCFLyN=0d0; maxCFLxTs=0d0; maxCFLyTs=0d0;
    maxSourceE=0d0; maxSourceH=0d0; maxGainsE=0d0; maxGainsH=0d0
    
   !$OMP PARALLEL DEFAULT (PRIVATE) SHARED (dt, UeNew, UhNew, TeNew, ThNew, TsNew, NeNew, NhNew, &
   !$OMP& TeDual, ThDual, TsDual, NeDual, NhDual, intensityDual, &
   !$OMP& Ue, Uh, Te, Th, Ts, Ne, Nh, GradNeX, GradNeY, intensity, intensity2, reflectivity, FermiTableE, FermiTableH, &
   !$OMP& Dielectric, DielectricDrudeE, DielectricDrudeH, absorptionDrudeE, absorptionDrudeH, &
   !$OMP& x, y, xDual, yDual, xDualSW, yDualSW, xDualSE, yDualSE, xDualNE, yDualNE, xDualNW, yDualNW, &
   !$OMP& diffusionE, diffusionH, GainsE, GainsH, LossesE, LossesH, &
   !$OMP& kappae, kappah, kappas, Ce, CeOld, Ch, ChOld, Cs, CouplingE, CouplingH, &
   !$OMP& nuColl, nuColleph, mobilityE, mobilityH, etae, etah, Egap, &
   !$OMP& SourceE, SourceH, SourceUe, SourceUh, diffNe, diffNh, CFLxT, CFLyT, CFLxN, CFLyN, CFLxTs, CFLyTs, &
   !$OMP& ThermalEnergy, LaserEnergy, epsilonInf, FermiIndexE, FermiIndexH, FermiRatioE, FermiRatioH, &
   !$OMP& OmegaX, OmegaY, JeX, JeY, JhX, JhY, VeX, VeY, VhX, VhY, DielectricStatic, Amatrix, Xvector, XvectorPrev, Bvector, xV, yV, xP, yP, &
   !$OMP& BsolveP, XsolveP, XsolvePrev, ABsolveP, KLsolve, KUsolve, Nsolve, ErrorVec, ipiv, info, &
   !$OMP& spectralNorm, Ex, Ey, ExPoisson, EyPoisson, potential, potentialNeedle, NeP, NhP, FixedPotentialIndex, &
   !$OMP& NormalNx, NormalNy, NormalSx, NormalSy, NormalEx, NormalEy, NormalWx, NormalWy, &
   !$OMP& CellVolume, CellAreaN, CellAreaS, CellAreaE, CellAreaW, CellVol, CellAreaNP, CellAreaSP, CellAreaEP, CellAreaWP, &
   !$OMP& NormalNxP, NormalNyP, NormalSxP, NormalSyP, NormalExP, NormalEyP, NormalWxP, NormalWyP, TangentWx, TangentWy, TangentNx, &
   !$OMP& TangentNy, TangentSx, TangentSy, TangentEx, TangentEy, CurviNx, CurviNy, CurviSx, CurviSy, CurviEx, CurviEy, CurviWx, &
   !$OMP& CurviWy, ConstBLx, ConstBLy, DistN, DistS, DistE, DistW, DistDualN, DistDualS, DistDualE, DistDualW, &
   !$OMP& EintField, EintFieldDual, EintFieldI, EintFieldR) &
   !$OMP& FIRSTPRIVATE (t, t0, x0, y0, I0, I1, I2, I3, I4, I5, I6, I7, &
   !$OMP& x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7, &
   !$OMP& sigmaX1, sigmaX2, sigmaX3, sigmaX4, sigmaX5, sigmaX6, sigmaX7, &
   !$OMP& sigmaY1, sigmaY2, sigmaY3, sigmaY4, sigmaY5, sigmaY6, sigmaY7, &
   !$OMP& AugerRateE, AugerRateH, sigmaTau, sigmaX, sigmaY, dx, dy, nbiter, &
   !$OMP& cpuefficiency, ColFermi0, ColFermi1, ColFermi2, ColFermiEta, ColFermiHalf, &
   !$OMP& ColFermiMenusHalf, ColFermiNeNc, ColFermiThreeHalf, SORsum, Mp, Np, NeTotal, NhTotal, ErrorSum)

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
   
   NeTotal=0d0; NhTotal=0d0
   
   ! replacing old datas
   
      !$OMP DO
      do i=1,M
        !$OMP PARALLEL DO
        do j=1,N
           Ue(i,j)=UeNew(i,j)
           Uh(i,j)=UhNew(i,j)
           Te(i,j)=TeNew(i,j)
	   Th(i,j)=ThNew(i,j)
	   Ts(i,j)=TsNew(i,j)
	   Ne(i,j)=NeNew(i,j)
	   Nh(i,j)=NhNew(i,j)
	   CeOld(i,j)=Ce(i,j)
	   ChOld(i,j)=Ch(i,j)
        end do
        !$OMP END PARALLEL DO
      end do
      !$OMP END DO
   

      
! if(PoissonOn==1) then
! 
!     !!!! POISSON interpolation before calculation *****
!    if(MeshChoice.eq.2) then
! !    write(*,*) "[Interpolation] TTM to Poisson"
! 
!     if(InterpolateOff.eq.0) then
!       
!     !$O M P SECTIONS
!       !$O M P SECTION
!       if(InterpolateMethod.eq.0) then
! 	NeP(:,:)=Interpolate(Ne,x,y,xP,yP,M,N,Mp,Np)
!       else 
! ! 	NeP(:,:)=InterpolateBiCubic(Ne,x,y,xP,yP,M,N,Mp,Np)
! 	call InterpolateBiCubic(Ne,x,y,xP,yP,M,N,Mp,Np, NeP, DummyVessel, DummyVessel)
!       end if
!       !$O M P SECTION
!       if(InterpolateMethod.eq.0) then
! 	NhP(:,:)=Interpolate(Nh,x,y,xP,yP,M,N,Mp,Np)
!       else 
! ! 	NhP(:,:)=InterpolateBiCubic(Nh,x,y,xP,yP,M,N,Mp,Np)
! 	call InterpolateBiCubic(Nh,x,y,xP,yP,M,N,Mp,Np, NhP, DummyVessel, DummyVessel)
!       end if
! !     !$O M P SECTION
!       if(InterpolateMethod.eq.0) then
! 	potential(:,:)=Interpolate(potentialNeedle,x,y,xP,yP,M,N,Mp,Np)
!       else 
! ! 	potential(:,:)=InterpolateBiCubic(potentialNeedle,x,y,xP,yP,M,N,Mp,Np)
! 	call InterpolateBiCubic(potentialNeedle,x,y,xP,yP,M,N,Mp,Np, potential, ExPoisson, EyPoisson)
! 	ExPoisson=-ExPoisson
! 	EyPoisson=-EyPoisson
!       end if
!       !$O M P END SECTIONS
!     endif
! 
!    else
!     NeP(:,:)=Ne(:,:)
!     NhP(:,:)=Nh(:,:)
!    end if
! 
!   !!!!! resolution of Poisson equation
!   ! write(*,*) "Test 1D index: i=10, j=10. Index=", OneDindex(10,10)
!   !     ! filling matrix & vectors
!   
!   Xvector(:)=0d0; Bvector(:)=0d0; 
! !   if(PoissonSolver.eq.1) then
! !     Amatrix(:,:)=0d0;
! !   end if
!   
!   ! initializing
!     do i=1,Mp
!       do j=1,Np
! 	Xvector(OneDindex(i,j))=potential(i,j) !OK
! 	Bvector(OneDindex(i,j))=0d0		!OK
! 	if(PoissonSolver.eq.1) then
! 		Amatrix(OneDindex(i,j),OneDindex(i,j))=1d0 ! Diag(A)=1d0 for Dirichlet boundary conditions !OK
! 	end if
!       end do
!     end do
! 
!   !!! filling out from domain boundaries
!     do i=2, Mp-1 !pour chaque point ou l on va calculer le potentiel
!       do j=2, Np-1 
! 	CellVolume(i,j)=(1d0/8d0) * ((xP(i+1,j+1)-xP(i-1,j-1)) * (yP(i-1,j+1)-yP(i+1,j-1))-(xP(i-1,j+1)-xP(i+1,j-1)) * (yP(i+1,j+1)-yP(i-1,j-1)) )
! 	if(MeshChoice<=1) then 
! 	  Bvector(OneDindex(i,j))=0d0
! 	else !this case will need interpolation between the meshes
! 	  Bvector(OneDindex(i,j))=-ActivateInduction*ec*(NhP(i,j)-NeP(i,j))/epsilon0*CellVolume(i,j) !but we should also write the values of boundaries which are very different
!         end if
!         if(PoissonSolver.eq.1) then
! 	  Amatrix(OneDindex(i,j),OneDindex(i,j))=(1d0/8d0)*(-DielectricStatic(i+1,j)-DielectricStatic(i,j))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i,j)+xP(i,j)**2+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)-(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i-1,j))*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2 & 
! 	  -2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i-1,j)**2-2d0*xP(i-1,j)*xP(i,j)+xP(i,j)**2+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)+(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i,j+1))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)/(xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i,j)+xP(i,j)**2+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i,j)+yP(i,j)**2)**(0.5d0)-(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i,j-1))*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)- &
! 	  2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5d0)/(xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i,j)+xP(i,j)**2+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! 	  Amatrix(OneDindex(i,j),OneDindex(i,j-1))=-(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i,j-1))*(xP(i+1,j)**2+2d0*xP(i+1,j)*xP(i+1,j-1)-2d0*xP(i+1,j)*xP(i-1,j-1)-2d0*xP(i+1,j)*xP(i-1,j)+xP(i+1,j-1)**2-2d0*xP(i+1,j-1)*xP(i-1,j-1)-2d0*xP(i+1,j-1)*xP(i-1,j)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i-1,j)+xP(i-1,j)**2+yP(i+1,j)**2+2d0*yP(i+1,j)*yP(i+1,j-1)-2d0*yP(i+1,j)*yP(i-1,j-1)-2d0*yP(i+1,j)*yP(i-1,j)+yP(i+1,j-1)**2-2d0*yP(i+1,j-1)*yP(i-1,j-1)-2d0*yP(i+1,j-1)*yP(i-1,j)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i-1,j)+yP(i-1,j)**2)**(0.5d0)/(xP(i,j-1)**2-2d0*xP(i,j-1)*xP(i,j)+xP(i,j)**2+yP(i,j-1)**2-2d0*yP(i,j-1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! 	  Amatrix(OneDindex(i,j),OneDindex(i,j+1))=(1d0/8d0)*(DielectricStatic(i,j)+DielectricStatic(i,j+1))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i+1,j)-2d0*xP(i+1,j+1)*xP(i-1,j)-2d0*xP(i+1,j+1)*xP(i-1,j+1)+xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i-1,j)-2d0*xP(i+1,j)*xP(i-1,j+1)+xP(i-1,j)**2+2d0*xP(i-1,j)*xP(i-1,j+1)+xP(i-1,j+1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i+1,j)-2d0*yP(i+1,j+1)*yP(i-1,j)-2d0*yP(i+1,j+1)*yP(i-1,j+1)+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i-1,j)-2d0*yP(i+1,j)*yP(i-1,j+1)+yP(i-1,j)**2+2d0*yP(i-1,j)*yP(i-1,j+1)+yP(i-1,j+1)**2)**(0.5d0)/(xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i,j)+xP(i,j)**2+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! 	  Amatrix(OneDindex(i,j),OneDindex(i-1,j))=-(1d0/8d0)*(-DielectricStatic(i,j)-DielectricStatic(i-1,j))*(xP(i,j+1)**2+2d0*xP(i,j+1)*xP(i-1,j+1)-2d0*xP(i,j+1)*xP(i-1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i-1,j+1)**2-2d0*xP(i-1,j+1)*xP(i-1,j-1)-2d0*xP(i-1,j+1)*xP(i,j-1)+xP(i-1,j-1)**2+2d0*xP(i-1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i,j+1)**2+2d0*yP(i,j+1)*yP(i-1,j+1)-2d0*yP(i,j+1)*yP(i-1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i-1,j+1)**2-2d0*yP(i-1,j+1)*yP(i-1,j-1)-2d0*yP(i-1,j+1)*yP(i,j-1)+yP(i-1,j-1)**2+2d0*yP(i-1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i-1,j)**2-2d0*xP(i-1,j)*xP(i,j)+xP(i,j)**2+yP(i-1,j)**2-2d0*yP(i-1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! 	  Amatrix(OneDindex(i,j),OneDindex(i+1,j))=(1d0/8d0)*(DielectricStatic(i+1,j)+DielectricStatic(i,j))*(xP(i+1,j+1)**2+2d0*xP(i+1,j+1)*xP(i,j+1)-2d0*xP(i+1,j+1)*xP(i+1,j-1)-2d0*xP(i+1,j+1)*xP(i,j-1)+xP(i,j+1)**2-2d0*xP(i,j+1)*xP(i+1,j-1)-2d0*xP(i,j+1)*xP(i,j-1)+xP(i+1,j-1)**2+2d0*xP(i+1,j-1)*xP(i,j-1)+xP(i,j-1)**2+yP(i+1,j+1)**2+2d0*yP(i+1,j+1)*yP(i,j+1)-2d0*yP(i+1,j+1)*yP(i+1,j-1)-2d0*yP(i+1,j+1)*yP(i,j-1)+yP(i,j+1)**2-2d0*yP(i,j+1)*yP(i+1,j-1)-2d0*yP(i,j+1)*yP(i,j-1)+yP(i+1,j-1)**2+2d0*yP(i+1,j-1)*yP(i,j-1)+yP(i,j-1)**2)**(0.5d0)/(xP(i+1,j)**2-2d0*xP(i+1,j)*xP(i,j)+xP(i,j)**2+yP(i+1,j)**2-2d0*yP(i+1,j)*yP(i,j)+yP(i,j)**2)**(0.5d0)
! ! 	  
! 	  if(CathodeZone.eq.1) then
! 	    do k=1,N ! filling the indexes of fixed potential with 1d0, not localized on boundaries
! 	      Amatrix(OneDindex(FixedPotentialIndex(k,1), FixedPotentialIndex(k,2)),&
! 		      OneDindex(FixedPotentialIndex(k,1), FixedPotentialIndex(k,2)))=1d0
! 	    end do
! 	  end if
! 	end if
!       end do
!     end do
!       
! 
!   
!   ! filling vessel boundaries
!   do j=1,Np
!     if(CathodeZone.eq.1) then
!       Bvector(OneDindex(Mp,j))=potentialNull !east: electrode. Others points are null and already defined.
!     else 
!       Bvector(OneDindex(Mp,j))=potential0 !east: electrode. Others points are null and already defined.
!     end if
!     Bvector(OneDindex(1,j))=potentialNull !apex: fictive mass here.
!   end do
!   do i=1,Mp
!     Bvector(OneDindex(i,1))=potentialNull !south: fictive mass
!     Bvector(OneDindex(i,Np))=potentialNull !north: fictive mass
!   end do
!   ! filling needle base boundary
! 
!   if(CathodeZone.eq.1) then
!     do i=FixedPotentialIndex(1,2)-ShiftFixedPotential, FixedPotentialIndex(N,2)+1+ShiftFixedPotential
!       Bvector(OneDindex(FixedPotentialIndex(1,1),i))=potential0
! !       Bvector(OneDindex(FixedPotentialIndex(1,1)+1,i))=potential0! potential0 !potential is imposed on each index in a continous way
!     end do
! !     do i=1, N
! ! !       write(*,*) FixedPotentialIndex(i,:)
! !       Bvector(OneDindex(FixedPotentialIndex(i,1),FixedPotentialIndex(i,2)))=potential0 ! every point is not written by this method
! !     end do
!   end if
! !   write(91,'(400E14.5)') Bvector
! !   call flush(91)
! ! 
! 
! !!!!********* overrelaxation algorithm till convergence, but not parallel ************
!     if(PoissonSolver.eq.1) then
! 	do k=1,SORiterations
! 	  do i=1,Mp*Np
! 	  SORsum=0d0 !(0d0,0d0)
! 	    do j=1,Mp*Np
! 	      if(j.ne.i) then 
! 		SORsum=SORsum+Amatrix(i,j)*Xvector(j)
! 	      end if
! 	    end do
!     ! 	  write(*,*) Amatrix(i,j)*Xvector(j)
!     !         write(91,'(9E14.5)') Xvector
! 	    XvectorPrev(i)=Xvector(i)
! 	    Xvector(i)=Xvector(i)+SORcoeff*((Bvector(i)-SORsum)/Amatrix(i,i)-Xvector(i))
! 	  end do
!     !       write(91,'(9E14.5)') Amatrix
! ! 	  call flush(91)
! 	  
! 	  ! convergence test
! 	  spectralNorm(k)=0d0
! 	  do j=1,Mp*Np
! 	    spectralNorm(k)=spectralNorm(k)+(Xvector(j)-XvectorPrev(j))**2
! 	  end do
! 	  spectralNorm(k)=sqrt(spectralNorm(k))
! 	  write(90,*) nbiter, k, spectralNorm(k)
! 	  
! 	end do
! 
!     !     write(*,*) Xvector, "\n"
!     else ! otherwise lets just apply the gauss inverted Amatrix to X vector
!       Xvector=matmul(Amatrix,Bvector)
!     end if
! 
! !     extract & check the solutions
!     do i=1,Mp    
!       do j=1,Np
!         potential(i,j)=Xvector(OneDindex(i,j))
! !         write(*,*) "SOR: Loop number", k, "potential(center)=", Xvector(OneDindex(3,3)), potential(3,3)	! controle du potentiel verifie juste apres le calcul !
! 	if(potential(i,j).eq."NaN") then 
! 	   write(95,*) "Divergence of potential at t=", t, "x,y=", x, y
! 	   write(*,*) "Divergence of potential at t=", t, "x,y=", x, y
! 	   Diverged=.true.
! 	end if
!       end do
!     end do
!     
!     ! calculate the fields !
!     !$O M P PARALLEL DO
!     do i=2,Mp
!       !$O M P DO
!       do j=2,Np
! 	distX=sqrt((xP(i+1,j)-xP(i-1,j))**2+(yP(i+1,j)-yP(i-1,j))**2)
! 	distY=sqrt((xP(i,j+1)-xP(i,j-1))**2+(yP(i,j+1)-yP(i,j-1))**2)
! 	! direct calculation: finite difference like
! ! 	ExPoisson(i,j)=(potential(i+1,j)-potential(i-1,j))/distX
! ! 	EyPoisson(i,j)=(potential(i,j+1)-potential(i,j-1))/distY
! 	! vector between nodes, scalar on node like
! ! 	ExPoisson(i,j)=-(-2d0*((potential(i,j)-potential(i-1,j))/distX)-ExPoisson(i-1,j))
! ! 	EyPoisson(i,j)=-(-2d0*((potential(i,j)-potential(i,j-1))/distY)-ExPoisson(i,j-1))
! 	! finite volume formulation
! 	ExPoisson(i,j) = -(0.5d0*(potential(i+1,j)+potential(i,j))*NormalExP(i,j)*CellAreaEP(i,j)+0.5d0*(potential(i-1,j)+potential(i,j))*NormalWxP(i,j)*CellAreaWP(i,j)+0.5d0*(potential(i,j+1)+potential(i,j))*NormalNxP(i,j)*CellAreaNP(i,j)+0.5d0*(potential(i,j)+potential(i,j-1))*NormalSxP(i,j)*CellAreaSP(i,j))/CellVolume(i,j)
! 	EyPoisson(i,j) = -(0.5d0*(potential(i+1,j)+potential(i,j))*NormalEyP(i,j)*CellAreaEP(i,j)+0.5d0*(potential(i-1,j)+potential(i,j))*NormalWyP(i,j)*CellAreaWP(i,j)+0.5d0*(potential(i,j+1)+potential(i,j))*NormalNyP(i,j)*CellAreaNP(i,j)+0.5d0*(potential(i,j)+potential(i,j-1))*NormalSyP(i,j)*CellAreaSP(i,j))/CellVolume(i,j)
! 	
! 	if(ExPoisson(i,j).eq."NaN") then
! 	  write(*,*) CellVolume(i,j)
! 	end if
! 	
!       end do
!       !$O M P END DO
!     end do
!     !$O M P END PARALLEL DO
! 
!     !!!! POISSON interpolation after calculation ***** dummy form working if MeshChoice<2 !!! TO WRITE
!     if(MeshChoice.eq.2) then
!     
! !       write(*,*) "[Interpolation] Poisson --> TTM"
!       if(InterpolateOff.eq.0) then
!       !$O M P SECTIONS
! 	!$O M P SECTION
! 	if(InterpolateMethod.eq.0) then 
! 	  Ex=Interpolate(ExPoisson, xP, yP, x, y, Mp, Np, M, N) !from Poisson to TTM
! 	else 
! ! 	  Ex=InterpolateBiCubic(ExPoisson, xP, yP, x, y, Mp, Np, M, N) !from Poisson to TTM
! 	  call InterpolateBiCubic(potential, xP, yP, x, y, Mp, Np, M, N, potentialNeedle, Ex,Ey) !from Poisson to TTM
! 	  Ex=-Ex
! 	  Ey=-Ey
! 	end if
! 	!$O M P SECTION
! 	if(InterpolateMethod.eq.0) then 
! 	  Ey=Interpolate(EyPoisson, xP, yP, x, y, Mp, Np, M, N)
! ! 	else 
! ! 	  Ey=InterpolateBiCubic(EyPoisson, xP, yP, x, y, Mp, Np, M, N)
! 	end if
!   
! 	!$O M P SECTION
! 	if(InterpolateMethod.eq.0) then
! 	  potentialNeedle(:,:)=Interpolate(potential, xP, yP, x, y, Mp, Np, M, N)
! ! 	else 
! ! 	  potentialNeedle(:,:)=InterpolateBiCubic(potential, xP, yP, x, y, Mp, Np, M, N)
! 	end if
! ! 	if(CathodeZone==1) then
! ! 	  potentialNeedle(M,:)=potential0
! ! 	endif
!       !$O M P END SECTIONS
!       end if
!     else 
!       Ex(:,:)=ExPoisson(:,:)
!       Ey(:,:)=EyPoisson(:,:)
!       potentialNeedle(:,:)=potential(:,:)
!     end if
! !     write(*,*) potentialNeedle(M-1,N-1)
! end if !Poisson ON
!     
    do i=1,M
      do j=1,N
	   NeTotal=NeTotal + Ne(i,j) * CellVol(i,j) 
	   NhTotal=NhTotal + Nh(i,j) * CellVol(i,j)
      end do
    end do
    
!!!! thermal calculations in the main domain
! calculation of sources
    !$OMP DO
    do i=1,M
!      intensity(i,1)=(1d0-reflectivity(i,1))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2.-.5d0*((x(i,1)-x0)/sigmaX)**2.-.5d0*((y(i,1)-y0)/sigmaY)**2.)
      !$OMP PARALLEL DO
      do j=1,N
	! optical coefficients
	
! 	! DEBUG


! 	Ex(i,j)=0d0 !-1d9
! 	Ey(i,j)=0d0*	-1d9
	
	nuColl(i,j)=CollisionFrequency()
	nuColleph(i,j)=ephCollisionFrequency()
	Dielectric(i,j)=DielectricFunction(lambda, epsilonInf, Ne(i,j), nuColl(i,j))
	DielectricDrudeE(i,j)=DielectricFunctionDrude(Ne(i,j), nuColl(i,j),me) 
	DielectricDrudeH(i,j)=DielectricFunctionDrude(Nh(i,j), nuColl(i,j),mh) 

	! 	write(*,*) "Esprit es-tu la ?"
	
! 	write(*,*) i,j,TeNew(i,j)
	DOSe(i,j)=DensityOfStateE(Te(i,j))
	DOSh(i,j)=DensityOfStateH(Th(i,j))
	FermiRatioE(i,j)=Ne(i,j)/DOSe(i,j)
	FermiRatioH(i,j)=Nh(i,j)/DOSh(i,j)
	FermiIndexE(i,j)=FermiIndex(FermiRatioE(i,j)) !1
	FermiIndexH(i,j)=FermiIndex(FermiRatioH(i,j)) !1
! 	write(*,*) "iter=", nbiter, "DOS=", DOSe(i,j), DOSh(i,j)
	etae(i,j)=FermiTableE(ColFermiEta,FermiIndexE(i,j)) 
	etah(i,j)=FermiTableH(ColFermiEta,FermiIndexH(i,j))
! 	write(*,*) "iter=", nbiter, "NeNc=", Ne(i,j)/DOSe(i,j), Nh(i,j)/DOSh(i,j)
! 	write(*,*) "iter=", nbiter, "FermiIndex=", FermiIndex(Ne(i,j)/DOSe(i,j)), FermiIndex(Nh(i,j)/DOSh(i,j))
! 	write(*,*) "iter=", nbiter, "eta=", etae(i,j), etah(i,j)
! 	write(*,*) "FermiTables: etaE,etaH=", FermiTableE(3,463), FermiTableH(3,450)

	mobilityE(i,j)=(ec/(me*nuColl(i,j)))*FermiTableE(ColFermi0, FermiIndexE(i,j))/FermiTableE(ColFermiHalf, FermiIndexE(i,j))
	mobilityH(i,j)=(ec/(mh*nuColl(i,j)))*FermiTableH(ColFermi0, FermiIndexH(i,j))/FermiTableH(ColFermiHalf, FermiIndexH(i,j))
	
! 	write(*,*) "iter=", nbiter, "mobility=", mobilityE(i,j), mobilityH(i,j)
	if(DrudeMode==0) then
	  absorptionDrudeE(i,j)=4d0*pi/lambda*aimag(sqrt(1d0+DielectricDrudeE(i,j))) ! overestimates the free-carrier absorption
	  absorptionDrudeH(i,j)=4d0*pi/lambda*aimag(sqrt(1d0+DielectricDrudeH(i,j)))
	else 
	  absorptionDrudeE(i,j)=sqrt(2d0)*sqrt(mu0*omegaLaser*ec*mobilityE(i,j)*Ne(i,j))
	  absorptionDrudeH(i,j)=sqrt(2d0)*sqrt(mu0*omegaLaser*ec*mobilityH(i,j)*Nh(i,j))
	end if
	
! 	write(*,*) "iter=", nbiter, "absorption=", absorptionDrudeE(i,j), absorptionDrudeH(i,j)
	reflectivity(i,j)=(real(sqrt(Dielectric(i,j)))**2+aimag(sqrt(Dielectric(i,j)))**2-2d0*real(sqrt(Dielectric(i,j)))+1d0) &
			    /(real(sqrt(Dielectric(i,j)))**2+aimag(sqrt(Dielectric(i,j)))**2+2d0*real(sqrt(Dielectric(i,j))+1d0))
	!local intensity
	
! 	! DEBUG ZONE
! 	if(lambda.eq.343d-9) then
! 	! uniform distribution such as Elena
! ! 	intensity(i,j)=(1d0-reflectivity(i,j))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2)
! 	intensity(i,j)=I0*exp(-.5d0*((t-t0)/sigmaTau)**2)*exp(-0.5d0*(((y(i,j)-500d-9)/sigmaY)**2+(x(i,j)/sigmaX)**2))
! 	! with just nothing
! !   	intensity(i,j)=(1d0-reflectivity(i,j))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2)*exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2)
! 	! with beer-lambert
! ! 	intensity(i,j)=(-(absorptionDrude(i,j) &
! ! 			+4d0*pi/lambda*aimag(sqrt(epsilonInf)))*intensity(i,j-1) &
! ! 			-TwoPhotonIonizationRate(lambda)*intensity(i,j-1)**2 &
! ! 			)*x &
! ! 			+intensity(i,j-1)
! !  	end if
	

	! WITH EXTERNALLY ADJUSTED INPUTS
!	!Lumerical mode already contains the reflectivity. Although, it doesn't consider change of optical index with ionization. 
	if(lambda.eq.1030d-9) then 
	  ConstBLx=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate(lambda,epsilonInf)+1d0*TwoPhotonIonizationRate(lambda)) &
		    / (1d0*exp(-(OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*x0) &
		    * (OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
	  ConstBLy=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate(lambda,epsilonInf)+1d0*TwoPhotonIonizationRate(lambda)) &
		    / (1d0*exp(-(OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*y0) &
		    * (OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
! case 1030 nm distribution
	  !initial field distribution  
	  if(BandBendingInFDTD.eq.0) then
	    intensity(i,j)=(1d0-reflectivity(i,N))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2)*exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2)
	  else 
	  intensity(i,j)=(1d0-reflectivity(i,N))*I0*exp(-.5d0*((t-t0)/sigmaTau)**2) & 
			  *(exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
			  + I0*1d-4*( & 
			   exp(-0.5d0*(((x(i,j)-x(i,N))**2+(y(i,j)-y(i,N))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
			  +exp(-0.5d0*(((x(i,j)-x(1,j))**2+(y(i,j)-y(1,j))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
			  +exp(-0.5d0*(((x(i,j)-x(i,1))**2+(y(i,j)-y(i,1))**2)/((DefectThickness)/(2d0*sqrt(2d0*log(2d0))))**2)) &
			  ))
	  end if
	  ! corrections from FDTD calculations and recovering non-linear processes
	  intensity(i,j)=intensity(i,j)*((OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
			/ (-TwoPhotonIonizationRate(lambda)+exp(-(OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j) & 
			    +absorptionDrudeH(i,j))*(x(i,j)-x0))  * ConstBLx * (OnePhotonIonizationRate(lambda,epsilonInf) + &
			    absorptionDrudeE(i,j) + absorptionDrudeH(i,j)))) &
			*((OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
			/ (-TwoPhotonIonizationRate(lambda)+exp(-(OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j) & 
			    +absorptionDrudeH(i,j))*abs(y(i,j)-y0))  * ConstBLy * (OnePhotonIonizationRate(lambda,epsilonInf) + &
			    absorptionDrudeE(i,j) + absorptionDrudeH(i,j))))
	end if
! case 515 nm distribution
	if(lambda.eq.515d-9) then 
	  ConstBLx=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate(lambda,epsilonInf)+1d0*TwoPhotonIonizationRate(lambda)) &
		    / (1d0*exp(-(OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*x0) &
		    * (OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
	  ConstBLy=(absorptionDrudeE(i,j)+absorptionDrudeH(i,j)+OnePhotonIonizationRate(lambda,epsilonInf)+1d0*TwoPhotonIonizationRate(lambda)) &
		    / (1d0*exp(-(OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*y0) &
		    * (OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)))
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
! 			*(1d0+cos(2d0*pi*x(i,j)/periodX)*sin(2d0*pi*y(i,j)/periodY)) &
			*exp(-.5d0*((x(i,j)-x9)/sigmaX9)**2)*exp(-.5d0*((y(i,j)-y9)/sigmaY9)**2)) &
			*((OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
			/ (-TwoPhotonIonizationRate(lambda)+exp(-(OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j) & 
			    +absorptionDrudeH(i,j))*(x(i,j)-x0))  * ConstBLx * (OnePhotonIonizationRate(lambda,epsilonInf) + &
			    absorptionDrudeE(i,j) + absorptionDrudeH(i,j)))) &
			*((OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) &
			/ (-TwoPhotonIonizationRate(lambda)+exp(-(OnePhotonIonizationRate(lambda,epsilonInf)+absorptionDrudeE(i,j) & 
			    +absorptionDrudeH(i,j))*abs(y(i,j)-y0))  * ConstBLy * (OnePhotonIonizationRate(lambda,epsilonInf) + &
			    absorptionDrudeE(i,j) + absorptionDrudeH(i,j))))
	end if
	if(lambda.eq.343d-9) then
	  intensity(i,j)= (1d0-reflectivity(i,N))* & 
			  I0*exp(-.5d0*((t-t0)/sigmaTau)**2) & 
			  *( &
			  exp(-(OnePhotonIonizationRate(lambda, epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j)) & 
			  *abs(y(i,j)-y(i,N)) & 
			  ) & 
			  * exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
! 			  + exp(-(OnePhotonIonizationRate(lambda, epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*abs(y(i,j)-y(i,1))) & 
! 			  *exp(-.5d0*((x(i,j)-x0)/sigmaX)**2)*exp(-.5d0*((y(i,j)-y0)/sigmaY)**2) &
			  )
! 			  *exp(-(OnePhotonIonizationRate(lambda, epsilonInf)+absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*abs(x(i,j)-x(i,N)))
	end if

	! USING MIE SCATTERING ANALYTICAL FORMULAS
	if(UseMieScattering .eq. 1) then
	! calculate electric field inside the tip
! 	  EintField(i,j)=Unit * MieScattering(abs(y(i,j)), phiMie, abs(ContourYofX(x(i,j), NeedleRadius, NeedleAngle)), Dielectric(i,j)) !*sqrt(2d0*fluence/(c*epsilon0*tau))
	  ! debug formula for constant cone radius
! 	  EintField(i,j)=MieScattering(abs(y(i,j)), phiMie, 100d-9, epsilonInf)
! 	  EintField(i,j)=sqrt(EintField(i,j)*conjg(EintField(i,j))) !complex to real
	  intensity(i,j)=0.5d0*c*epsilon0*EintField(i,j)*exp(-.5d0*((t-t0)/sigmaTau)**2) !laser fluence and reflectivity is inside the field
	end if

	
	if(intensity(i,j) < 1d-20) then 
	  intensity(i,j)=0d0
	end if
	! free-carrier balance sources
	Egap(i,j)=EgapValue(Ne(i,j),Ts(i,j))

	diffusionE(i,j)=mobilityE(i,j)*kb*Te(i,j)/ec*FermiTableE(ColFermiHalf,FermiIndexE(i,j))/FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j))
	diffusionH(i,j)=mobilityH(i,j)*kb*Th(i,j)/ec*FermiTableH(ColFermiHalf,FermiIndexH(i,j))/FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j))
	if(ConductivityFix.eq.-1) then
	   diffusionE(i,j)=0d0
	   diffusionH(i,j)=0d0
	end if	

	JeX(i,j)=-mobilityE(i,j)*Ne(i,j)*Ex(i,j)
	JeY(i,j)=-mobilityE(i,j)*Ne(i,j)*Ey(i,j)
	JhX(i,j)=mobilityH(i,j)*Nh(i,j)*Ex(i,j)
	JhY(i,j)=mobilityH(i,j)*Nh(i,j)*Ey(i,j)

	if(DriftOn.eq.0) then 
	  JeX(i,j)=0d0; JeY(i,j)=0d0; 
	  JhX(i,j)=0d0; JhY(i,j)=0d0;
	endif

	GainsE(i,j)=(OnePhotonIonizationRate(lambda,epsilonInf)*intensity(i,j)/hbar/omegaLaser &
		    +TwoPhotonIonizationRate(lambda)*intensity(i,j)**2/(2d0*hbar*omegaLaser) &
		    +ImpactIonizationRate(Te(i,j),Ne(i,j),Ts(i,j))*Ne(i,j))! *(4d0*SiDensity-Ne(i,j))/(4d0*SiDensity) !use Old Ne here! 
		    
	GainsH(i,j)=(OnePhotonIonizationRate(lambda,epsilonInf)*intensity(i,j)/hbar/omegaLaser &
		    +TwoPhotonIonizationRate(lambda)*intensity(i,j)**2/(2d0*hbar*omegaLaser) &
		    +ImpactIonizationRate(Te(i,j),Ne(i,j),Ts(i,j))*Nh(i,j)) !*(4d0*SiDensity-Ne(i,j))/(4d0*SiDensity) !use Old Nh here
		    
	LossesE(i,j)=AugerRateE*Ne(i,j)**2*Nh(i,j)+AugerRateH*Nh(i,j)**2*Ne(i,j) !use Old Ne, Nh here! 
	LossesH(i,j)=LossesE(i,j)
	
	! thermal coefficients
	kappae(i,j)=kb**2*Ne(i,j)*mobilityE(i,j)*Te(i,j)/ec*(6d0*FermiTableE(ColFermi2,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) &
							  -4d0*(FermiTableE(ColFermi1,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)))**2)
	kappah(i,j)=kb**2*Nh(i,j)*mobilityH(i,j)*Th(i,j)/ec*(6d0*FermiTableH(ColFermi2,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) &
							  -4d0*(FermiTableH(ColFermi1,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)))**2)
! 	kappas(i,j)=-.1412d0*Ts(i,j)**(1.38961d0)+0.638157d0*Ts(i,j)**(1.14013d0) !mingo till 300 K
	kappas(i,j)=(-8.992d0+68.265d0/(1d0+exp(-.4075612391d-1*Ts(i,j)+2.315984470d0))*(1d0-1d0/(1d0+exp(-.4756634637d-2*Ts(i,j)+2.403533689d0)))) !Elena fit sur Kazan (2010)
	if(kappas(i,j).lt.0d0) then
	  kappas(i,j)=0d0
	end if
	
	! correction considering Fick diffusion in energy
	if(ConductivityFix.eq.1) then 
	    kappae(i,j)=kappae(i,j) + kb**2*Te(i,j)*Ne(i,j)*mobilityE(i,j) / ec &
		      * (etae(i,j) - 2d0*FermiTableE(ColFermi1,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) )**2 
	    kappah(i,j)=kappah(i,j) + kb**2*Th(i,j)*Nh(i,j)*mobilityH(i,j) / ec &
		      * (etah(i,j) - 2d0*FermiTableH(ColFermi1,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) )**2 
	else if(ConductivityFix.eq.2) then 
	    kappae(i,j)=kappae(i,j) + 2d0*kb**2*Te(i,j)*FermiTableE(ColFermi1,FermiIndexE(i,j))*mobilityE(i,j)*FermiTableE(ColFermiHalf, FermiIndexE(i,j))*Ne(i,j) * &
			(2d0*FermiTableE(ColFermi1, FermiIndexE(i,j))*FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)) & 
			/FermiTableE(ColFermiHalf,FermiIndexE(i,j))/FermiTableE(ColFermi0,FermiIndexE(i,j)) - 1.5d0) * & 
			(FermiTableE(ColFermi0, FermiIndexE(i,j))*ec*FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))**(-1e0)
			
	    kappah(i,j)=kappah(i,j) + 2d0*kb**2*Te(i,j)*FermiTableH(ColFermi1,FermiIndexH(i,j))*mobilityH(i,j)*FermiTableH(ColFermiHalf, FermiIndexH(i,j))*Ne(i,j) * &
			(2d0*FermiTableH(ColFermi1, FermiIndexH(i,j))*FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)) & 
			/FermiTableH(ColFermiHalf,FermiIndexH(i,j))/FermiTableH(ColFermi0,FermiIndexH(i,j)) - 1.5d0) * & 
			(FermiTableH(ColFermi0, FermiIndexH(i,j))*ec*FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j)))**(-1.)
			
	else if(ConductivityFix.eq.-1) then
	   kappae(i,j)=0d0
	   kappah(i,j)=0d0
	   kappas(i,j)=0d0
	end if
	Ce(i,j)=1.5d0*Ne(i,j)*kb*(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)) &
		  -etae(i,j)*(1d0-FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))* &
		  FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))
	Ch(i,j)=1.5d0*Nh(i,j)*kb*(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)) &
		  -etah(i,j)*(1d0-(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)))* &
		  FermiTableH(ColFermiMenusHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j))))
	! Cs(i,j)=1d3*SiDensity*0.2703d0/(exp(63.456d0/Ts(i,j))+0.84586d0) !bad fit ...
	Cs(i,j)=1d3*SiDensity*0.4129d0/(exp(88.183d0/Ts(i,j))-0.676d0) !!better fit on Flubacher

	CouplingE(i,j)=Ce(i,j)*nuColleph(i,j)*(Te(i,j)-Ts(i,j))
	CouplingH(i,j)=Ch(i,j)*nuColleph(i,j)*(Th(i,j)-Ts(i,j))
	
! 	diffNe(i,j)=0d0 !just for debug !
! 	diffNh(i,j)=0d0 !just for debug !!
	
	SourceE(i,j)= (hbar*omegaLaser-Egap(i,j))/hbar/omegaLaser*OnePhotonIonizationRate(lambda,epsilonInf)*intensity(i,j) &
		     + (2d0*hbar*omegaLaser - Egap(i,j))/(2d0*hbar*omegaLaser)*TwoPhotonIonizationRate(lambda)*intensity(i,j)**2 &
		     - Egap(i,j)*ImpactIonizationRate(Te(i,j),Ne(i,j),Ts(i,j))*Ne(i,j) &
		     + absorptionDrudeE(i,j)*intensity(i,j) &
		     + Egap(i,j)*AugerRateE*Ne(i,j)**2*Nh(i,j)

		     
	SourceUe(i,j) = SourceE(i,j)
	SourceE(i,j) = SourceE(i,j) - diffNe(i,j)*(1.5d0*kb*Te(i,j))*(FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)))

	SourceH(i,j)=(hbar*omegaLaser-Egap(i,j))/hbar/omegaLaser * OnePhotonIonizationRate(lambda,epsilonInf)*intensity(i,j) &
		     + (2d0*hbar*omegaLaser - Egap(i,j))/(2d0*hbar*omegaLaser) *TwoPhotonIonizationRate(lambda)*intensity(i,j)**2 &
		     - Egap(i,j)*ImpactIonizationRate(Th(i,j),Nh(i,j),Ts(i,j))*Nh(i,j) &
		     + absorptionDrudeH(i,j)*intensity(i,j) & 
		     + Egap(i,j)*AugerRateH*Nh(i,j)**2*Ne(i,j)
		     
	SourceUh(i,j) = SourceH(i,j)
	SourceH(i,j) = SourceH(i,j) - diffNh(i,j)*(1.5d0*kb*Th(i,j)*(FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j))))
	
	VeX(i,j)=0d0
	VeY(i,j)=0d0
	VhX(i,j)=0d0
	VhY(i,j)=0d0
	
      end do
      !$OMP END PARALLEL DO
      
    end do
    !$OMP END DO
    
    ! interpolation and preparation of resolution
    
    !$OMP DO
    do i=2,M-1
      !$OMP PARALLEL DO
      do j=2,N-1
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
      !$OMP END PARALLEL DO
    end do
    !$OMP END DO
    
    ! interpolation on dual mesh
    ! InterpolateBiCubic(phi_source, x_s, y_s, x_t, y_t, SizeXs, SizeYs, SizeXt, SizeYt, phi_target, Grad(phi)_targetX, Grad(phi)_targetY)
    
    if(UseInterpolation.eq.1) then
    !$O M P SECTIONS
      !$O M P SECTION
      call InterpolateBiCubic(Ne, x, y, xDual, yDual, M, N, M-1, N-1, NeDual, DummyDual, DummyDual)
      !$O M P SECTION
      call InterpolateBiCubic(Nh, x, y, xDual, yDual, M, N, M-1, N-1, NhDual, DummyDual, DummyDual)
      !$O M P SECTION
      call InterpolateBiCubic(Te, x, y, xDual, yDual, M, N, M-1, N-1, TeDual, DummyDual, DummyDual)
      !$O M P SECTION
      call InterpolateBiCubic(Th, x, y, xDual, yDual, M, N, M-1, N-1, ThDual, DummyDual, DummyDual)
      !$O M P SECTION
      call InterpolateBiCubic(Ts, x, y, xDual, yDual, M, N, M-1, N-1, TsDual, DummyDual, DummyDual)
      !$O M P SECTION
  !     call InterpolateBiCubic(intensity, x, y, xDual, yDual, M, N, M-1, N-1, intensityDual, DummyDual, DummyDual)
  !     intensityDual=Interpolate(intensity, x, y, xDual, yDual, M, N, M-1, N-1)
	
  ! !       DEBUG: test de la fonction d'interpolation
      call InterpolateBiCubic(intensity, x, y, xDual, yDual, M, N, M-1, N-1, intensityDual, DummyDual, DummyDual)
  !     intensityDual=Interpolate(intensity, x, y, xDual, yDual, M, N, M-1, N-1)

    !$O M P END SECTIONS
    end if
    
    
    
    !$OMP DO
    do i=2, M-1
      !$OMP PARALLEL DO
      do j=2, N-1
      ! define the volume of elementary cell around a point everywhere but not on boundaries
! 	CellVol(i,j)=0.25d0*(AreaElement(x(i-1,j-1),y(i-1,j-1),x(i+1,j-1),y(i+1,j-1),x(i+1,j+1),y(i+1,j+1),x(i-1,j+1),y(i-1,j+1)))
	CellVol(i,j)=AreaElement(0.25d0*(x(i-1,j-1)+x(i,j-1)+x(i,j)+x(i-1,j)), & !x(i-1/2,j-1/2)
				 0.25d0*(y(i-1,j-1)+y(i,j-1)+y(i,j)+y(i-1,j)), &
				 0.25d0*(x(i,j-1)+x(i+1,j-1)+x(i+1,j)+x(i,j)), &	!x(i+1/2,j-1/2)
				 0.25d0*(y(i,j-1)+y(i+1,j-1)+y(i+1,j)+y(i,j)), &
				 0.25d0*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1)), &	!x(i+1/2,j+1/2)
				 0.25d0*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1)), &
				 0.25d0*(x(i-1,j)+x(i,j)+x(i,j+1)+x(i-1,j+1)), &	!x(i-1/2,j+1/2)
				 0.25d0*(y(i-1,j)+y(i,j)+y(i,j+1)+y(i-1,j+1)))
				 
	! Calcul de Grad(Ne) sur le maillage direct
	! Première estimation peu stable
	GradNeX(i,j) = 1d0/CellVol(i,j) * &
		      ( 0.5d0 * (Ne(i,j) + Ne(i,j+1)) * CellAreaN(i,j) * NormalNx(i,j) & 
		      + 0.5d0 * (Ne(i,j) + Ne(i,j-1)) * CellAreaS(i,j) * NormalSx(i,j) & 
		      + 0.5d0 * (Ne(i,j) + Ne(i-1,j)) * CellAreaW(i,j) * NormalWx(i,j) & 
		      + 0.5d0 * (Ne(i,j) + Ne(i+1,j)) * CellAreaE(i,j) * NormalEx(i,j) )
	GradNeY(i,j) = 1d0/CellVol(i,j) * &
		      ( 0.5d0 * (Ne(i,j) + Ne(i,j+1)) * CellAreaN(i,j) * NormalNy(i,j) & 
		      + 0.5d0 * (Ne(i,j) + Ne(i,j-1)) * CellAreaS(i,j) * NormalSy(i,j) & 
		      + 0.5d0 * (Ne(i,j) + Ne(i-1,j)) * CellAreaW(i,j) * NormalWy(i,j) & 
		      * 0.5d0 * (Ne(i,j) + Ne(i+1,j)) * CellAreaE(i,j) * NormalEy(i,j) )
		      
	! interpolation lineaire des valeurs de phi sur les bords de cellules
! 	phiN=Ne(i,j)+GradNeX(i,j)*(0.5d0*(xDual(i,j)+xDual(i-1,j)))+GradNeY(i,j)*(0.5d0*(yDual(i,j)+yDual(i-1,j)))
! 	phiS=Ne(i,j)+GradNeX(i,j)*(0.5d0*(xDual(i,j-1)+xDual(i-1,j-1)))+GradNeY(i,j)*(0.5d0*(yDual(i,j-1)+yDual(i-1,j-1)))
	! recalcul du gradient avec phi_bords
	
	end do
	!$OMP END PARALLEL DO
      end do
      !$OMP END DO
      
      !$OMP DO
      do i=1, M-1
	!$OMP PARALLEL DO
	do j=1, N-1
	
      ! interpolation bilineaire ponderee par les aires
	NeDual(i,j) = ( Ne(i,j)/CellVol(i,j) + Ne(i+1,j)/CellVol(i+1,j) + Ne(i,j+1)/CellVol(i,j+1) + Ne(i+1,j+1) / CellVol(i+1,j+1) ) / ( 1d0 / CellVol(i,j) + 1d0 / CellVol(i+1,j) + 1d0/CellVol(i,j+1) + 1d0/CellVol(i+1,j+1) )
	NhDual(i,j) = ( Nh(i,j)/CellVol(i,j) + Nh(i+1,j)/CellVol(i+1,j) + Nh(i,j+1)/CellVol(i,j+1) + Nh(i+1,j+1) / CellVol(i+1,j+1) ) / ( 1d0 / CellVol(i,j) + 1d0 / CellVol(i+1,j) + 1d0/CellVol(i,j+1) + 1d0/CellVol(i+1,j+1) )
	TeDual(i,j) = ( Te(i,j)/CellVol(i,j) + Te(i+1,j)/CellVol(i+1,j) + Te(i,j+1)/CellVol(i,j+1) + Te(i+1,j+1) / CellVol(i+1,j+1) ) / ( 1d0 / CellVol(i,j) + 1d0 / CellVol(i+1,j) + 1d0/CellVol(i,j+1) + 1d0/CellVol(i+1,j+1) )
	ThDual(i,j) = ( Th(i,j)/CellVol(i,j) + Th(i+1,j)/CellVol(i+1,j) + Th(i,j+1)/CellVol(i,j+1) + Th(i+1,j+1) / CellVol(i+1,j+1) ) / ( 1d0 / CellVol(i,j) + 1d0 / CellVol(i+1,j) + 1d0/CellVol(i,j+1) + 1d0/CellVol(i+1,j+1) )
	TsDual(i,j) = ( Ts(i,j)/CellVol(i,j) + Ts(i+1,j)/CellVol(i+1,j) + Ts(i,j+1)/CellVol(i,j+1) + Ts(i+1,j+1) / CellVol(i+1,j+1) ) / ( 1d0 / CellVol(i,j) + 1d0 / CellVol(i+1,j) + 1d0/CellVol(i,j+1) + 1d0/CellVol(i+1,j+1) )
	end do
	!$OMP END PARALLEL DO
      end do
      !$OMP END DO
      
      !$OMP END PARALLEL !!end of parallel section
      !=======================================================

      
      ! solving the 2D problem
      
      !================================== IMPLICIT SOLVERS ===================================
      ! A try with implicit solvers
      ! 1/ Setting the boundary conditions
      ! 2/ Filling matrix directly using band data storage (cf Lapack)
      ! 3/ Find the solution with Lapack without using a large matrix
      
      ! Analyze of the matrix shape with boundary conditions


      
      
      ! ============ Solve Ne  ===============
! if(HolesOff.eq.0) then      
!       AsolveP(:,:)=0d0;
      BsolveP1(:)=0d0; ABsolveP1(:,:)=0d0; KLsolve1=Np; KUsolve1=Np; Nsolve1=Mp*Np
      BsolveP2(:)=0d0; ABsolveP2(:,:)=0d0; KLsolve2=Np; KUsolve2=Np; Nsolve2=Mp*Np
      BsolveP3(:)=0d0; ABsolveP3(:,:)=0d0; KLsolve3=Np; KUsolve3=Np; Nsolve3=Mp*Np
      BsolveP4(:)=0d0; ABsolveP4(:,:)=0d0; KLsolve4=Np; KUsolve4=Np; Nsolve4=Mp*Np
      BsolveP5(:)=0d0; ABsolveP5(:,:)=0d0; KLsolve5=Np; KUsolve5=Np; Nsolve5=Mp*Np

      !$OMP PARALLEL SHARED(KLsolve1, KUsolve1, Nsolve1, BsolveP1, ABsolveP1, &
      !$OMP& KLsolve2, KUsolve2, Nsolve2, BsolveP2, ABsolveP2, &
      !$OMP& KLsolve3, KUsolve3, Nsolve3, BsolveP3, ABsolveP3, &
      !$OMP& KLsolve4, KUsolve4, Nsolve4, BsolveP4, ABsolveP4, &
      !$OMP& KLsolve5, KUsolve5, Nsolve5, BsolveP5, ABsolveP5, &
      !$OMP& CellAreaE, CellAreaN, CellAreaS, CellAreaW, DistN, DistS, DistE, DistW, &
      !$OMP& NormalEx, NormalEy, NormalNx, NormalNy, NormalSx, NormalSy, NormalWx, NormalWy, &
      !$OMP& CurviEx, CurviEy, CurviNx, CurviNy, CurviSx, CurviSy, CurviWx, CurviWy, CellVol, &
      !$OMP& DistDualN, DistDualS, DistDualE, DistDualW, &
      !$OMP& Ne, NeDual, GainsE, LossesE, diffusionE, &
      !$OMP& Nh, NhDual, GainsH, LossesH, diffusionH, &
      !$OMP& Te, TeDual, SourceE, CouplingE, kappaE, Ce, &
      !$OMP& Th, ThDual, SourceH, CouplingH, kappaH, Ch, &
      !$OMP& Ts, TsDual, kappaS, Cs) &
      !$OMP& FIRSTPRIVATE(Mp, Np, dt)
      
      !$OMP DO COLLAPSE(2)
      do i=2,Mp-1
	do j=2,Np-1
	  !!Ne
	  if(FillMatrixCondition(OneDindex(i,j),OneDindex(i,j),KLsolve1, KUsolve1, Nsolve1)) then
	    ABsolveP1(FillMatrixRow(OneDindex(i,j), OneDindex(i,j),KLsolve1,KUsolve1,Nsolve1),OneDindex(i,j)) = CellVol(i,j)/dt + ( &
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j)) = CellVol(i,j)/dt + ( &
			+ 0.5d0*CellAreaE(i,j)*(diffusionE(i,j ) +diffusionE(i+1,j))*(1d0)/DistE(i,j) * (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j)) &
			+ 0.5d0*CellAreaW(i,j)*(diffusionE(i-1,j)+diffusionE(i,j) ) *(1d0)/DistW(i,j) * (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j)) &
			+ 0.5d0*CellAreaN(i,j)*(diffusionE(i,j+1)+diffusionE(i,j) ) *(1d0)/DistN(i,j) * (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j)) &
			+ 0.5d0*CellAreaS(i,j)*(diffusionE(i,j-1)+diffusionE(i,j) ) *(1d0)/DistS(i,j) * (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j)) &
			)
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j-1),KLsolve1,KUsolve1,Nsolve1)) then
	    ABsolveP1(FillMatrixRow(OneDindex(i,j), OneDindex(i,j-1),KLsolve1, KUsolve1, Nsolve1), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(diffusionE(i,j-1 ) + diffusionE(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(diffusionE(i,j-1 ) + diffusionE(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j+1),KLsolve1,KUsolve1,Nsolve1)) then
	    ABsolveP1(FillMatrixRow(OneDindex(i,j), OneDindex(i,j+1),KLsolve1, KUsolve1, Nsolve1), OneDindex(i,j+1))=0.5d0*CellAreaN(i,j)*(diffusionE(i,j+1 ) + diffusionE(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j+1)) =0.5d0*CellAreaN(i,j)*(diffusionE(i,j+1 ) + diffusionE(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i-1,j),KLsolve1,KUsolve1,Nsolve1)) then
	    ABsolveP1(FillMatrixRow(OneDindex(i,j), OneDindex(i-1,j),KLsolve1, KUsolve1, Nsolve1), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(diffusionE(i-1,j ) + diffusionE(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(diffusionE(i-1,j ) + diffusionE(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i+1,j),KLsolve1,KUsolve1,Nsolve1)) then
	    ABsolveP1(FillMatrixRow(OneDindex(i,j), OneDindex(i+1,j),KLsolve1, KUsolve1, Nsolve1), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(diffusionE(i,j )   + diffusionE(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(diffusionE(i,j )   + diffusionE(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
	  end if

  ! 	  ! dirichlet condnitions on boundaries
  ! 	  AsolveP(OneDindex(i,N), OneDindex(i,N)) = 1d0; BsolveP(OneDindex(i,N)) = Ne0; AsolveP(OneDindex(i,N), OneDindex(i,N-1)) = 0d0;
  ! 	  AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0; BsolveP(OneDindex(i,1)) = Ne0; AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = 0d0;
  ! 	  AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0; BsolveP(OneDindex(M,j)) = Ne0; AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = 0d0;
  ! 	  AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0; BsolveP(OneDindex(1,j)) = Ne0; AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = 0d0;

	  !!Nh

	  if(FillMatrixCondition(OneDindex(i,j),OneDindex(i,j),KLsolve2, KUsolve2, Nsolve2)) then
	    ABsolveP2(FillMatrixRow(OneDindex(i,j), OneDindex(i,j),KLsolve2,KUsolve2,Nsolve2),OneDindex(i,j)) = CellVol(i,j)/dt + ( &
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j)) = CellVol(i,j)/dt + ( &
			+ 0.5d0*CellAreaE(i,j)*(diffusionH(i,j ) +diffusionH(i+1,j))*(1d0)/DistE(i,j) * (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j)) &
			+ 0.5d0*CellAreaW(i,j)*(diffusionH(i-1,j)+diffusionH(i,j) ) *(1d0)/DistW(i,j) * (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j)) &
			+ 0.5d0*CellAreaN(i,j)*(diffusionH(i,j+1)+diffusionH(i,j) ) *(1d0)/DistN(i,j) * (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j)) &
			+ 0.5d0*CellAreaS(i,j)*(diffusionH(i,j-1)+diffusionH(i,j) ) *(1d0)/DistS(i,j) * (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j)) &
			)
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j-1),KLsolve2,KUsolve2,Nsolve2)) then
	    ABsolveP2(FillMatrixRow(OneDindex(i,j), OneDindex(i,j-1),KLsolve2, KUsolve2, Nsolve2), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(diffusionH(i,j-1 ) + diffusionH(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(diffusionH(i,j-1 ) + diffusionH(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j+1),KLsolve2,KUsolve2,Nsolve2)) then
	    ABsolveP2(FillMatrixRow(OneDindex(i,j), OneDindex(i,j+1),KLsolve2, KUsolve2, Nsolve2), OneDindex(i,j+1))=0.5d0*CellAreaN(i,j)*(diffusionH(i,j+1 ) + diffusionH(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j+1)) =0.5d0*CellAreaN(i,j)*(diffusionH(i,j+1 ) + diffusionH(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i-1,j),KLsolve2,KUsolve2,Nsolve2)) then
	    ABsolveP2(FillMatrixRow(OneDindex(i,j), OneDindex(i-1,j),KLsolve2, KUsolve2, Nsolve2), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(diffusionH(i-1,j ) + diffusionH(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(diffusionH(i-1,j ) + diffusionH(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i+1,j),KLsolve2,KUsolve2,Nsolve2)) then
	    ABsolveP2(FillMatrixRow(OneDindex(i,j), OneDindex(i+1,j),KLsolve2, KUsolve2, Nsolve2), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(diffusionH(i,j )   + diffusionH(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(diffusionH(i,j )   + diffusionH(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
	  end if

  ! 	  ! dirichlet condnitions on boundaries
  ! 	  AsolveP(OneDindex(i,N), OneDindex(i,N)) = 1d0; BsolveP(OneDindex(i,N)) = Ne0; AsolveP(OneDindex(i,N), OneDindex(i,N-1)) = 0d0;
  ! 	  AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0; BsolveP(OneDindex(i,1)) = Ne0; AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = 0d0;
  ! 	  AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0; BsolveP(OneDindex(M,j)) = Ne0; AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = 0d0;
  ! 	  AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0; BsolveP(OneDindex(1,j)) = Ne0; AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = 0d0;

	!!Te

	  if(FillMatrixCondition(OneDindex(i,j),OneDindex(i,j),KLsolve3, KUsolve3, Nsolve3)) then
	    ABsolveP3(FillMatrixRow(OneDindex(i,j), OneDindex(i,j),KLsolve3,KUsolve3,Nsolve3),OneDindex(i,j)) = Ce(i,j)*CellVol(i,j)/dt + ( &
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j)) = CellVol(i,j)/dt + ( &
			+ 0.5d0*CellAreaE(i,j)*(kappae(i,j ) +kappae(i+1,j))*(1d0)/DistE(i,j) * (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j)) &
			+ 0.5d0*CellAreaW(i,j)*(kappae(i-1,j)+kappae(i,j) ) *(1d0)/DistW(i,j) * (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j)) &
			+ 0.5d0*CellAreaN(i,j)*(kappae(i,j+1)+kappae(i,j) ) *(1d0)/DistN(i,j) * (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j)) &
			+ 0.5d0*CellAreaS(i,j)*(kappae(i,j-1)+kappae(i,j) ) *(1d0)/DistS(i,j) * (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j)) &
			)
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j-1),KLsolve3,KUsolve3,Nsolve3)) then
	    ABsolveP3(FillMatrixRow(OneDindex(i,j), OneDindex(i,j-1),KLsolve3, KUsolve3, Nsolve3), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(kappae(i,j-1 ) + kappae(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(kappae(i,j-1 ) + kappae(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j+1),KLsolve3,KUsolve3,Nsolve3)) then
	    ABsolveP3(FillMatrixRow(OneDindex(i,j), OneDindex(i,j+1),KLsolve3, KUsolve3, Nsolve3), OneDindex(i,j+1))=0.5d0*CellAreaN(i,j)*(kappae(i,j+1 ) + kappae(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j+1)) =0.5d0*CellAreaN(i,j)*(kappae(i,j+1 ) + kappae(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i-1,j),KLsolve3,KUsolve3,Nsolve3)) then
	    ABsolveP3(FillMatrixRow(OneDindex(i,j), OneDindex(i-1,j),KLsolve3, KUsolve3, Nsolve3), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(kappae(i-1,j ) + kappae(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(kappae(i-1,j ) + kappae(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i+1,j),KLsolve3,KUsolve3,Nsolve3)) then
	    ABsolveP3(FillMatrixRow(OneDindex(i,j), OneDindex(i+1,j),KLsolve3, KUsolve3, Nsolve3), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(kappae(i,j )   + kappae(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(kappae(i,j )   + kappae(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
	  end if

  ! 	  ! dirichlet condnitions on boundaries
  ! 	  AsolveP(OneDindex(i,N), OneDindex(i,N)) = 1d0; BsolveP3(OneDindex(i,N)) = Ne0; AsolveP(OneDindex(i,N), OneDindex(i,N-1)) = 0d0;
  ! 	  AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0; BsolveP3(OneDindex(i,1)) = Ne0; AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = 0d0;
  ! 	  AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0; BsolveP3(OneDindex(M,j)) = Ne0; AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = 0d0;
  ! 	  AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0; BsolveP3(OneDindex(1,j)) = Ne0; AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = 0d0;

      !!Th
	  if(FillMatrixCondition(OneDindex(i,j),OneDindex(i,j),KLsolve4, KUsolve4, Nsolve4)) then
	    ABsolveP4(FillMatrixRow(OneDindex(i,j), OneDindex(i,j),KLsolve4,KUsolve4,Nsolve4),OneDindex(i,j)) = Ch(i,j)*CellVol(i,j)/dt + ( &
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j)) = CellVol(i,j)/dt + ( &
			+ 0.5d0*CellAreaE(i,j)*(kappah(i,j ) +kappah(i+1,j))*(1d0)/DistE(i,j) * (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j)) &
			+ 0.5d0*CellAreaW(i,j)*(kappah(i-1,j)+kappah(i,j) ) *(1d0)/DistW(i,j) * (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j)) &
			+ 0.5d0*CellAreaN(i,j)*(kappah(i,j+1)+kappah(i,j) ) *(1d0)/DistN(i,j) * (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j)) &
			+ 0.5d0*CellAreaS(i,j)*(kappah(i,j-1)+kappah(i,j) ) *(1d0)/DistS(i,j) * (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j)) &
			)
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j-1),KLsolve4,KUsolve4,Nsolve4)) then
	    ABsolveP4(FillMatrixRow(OneDindex(i,j), OneDindex(i,j-1),KLsolve4, KUsolve4, Nsolve4), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(kappah(i,j-1 ) + kappah(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(kappah(i,j-1 ) + kappah(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j+1),KLsolve4,KUsolve4,Nsolve4)) then
	    ABsolveP4(FillMatrixRow(OneDindex(i,j), OneDindex(i,j+1),KLsolve4, KUsolve4, Nsolve4), OneDindex(i,j+1))=0.5d0*CellAreaN(i,j)*(kappah(i,j+1 ) + kappah(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j+1)) =0.5d0*CellAreaN(i,j)*(kappah(i,j+1 ) + kappah(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i-1,j),KLsolve4,KUsolve4,Nsolve4)) then
	    ABsolveP4(FillMatrixRow(OneDindex(i,j), OneDindex(i-1,j),KLsolve4, KUsolve4, Nsolve4), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(kappah(i-1,j ) + kappah(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(kappah(i-1,j ) + kappah(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i+1,j),KLsolve4,KUsolve4,Nsolve4)) then
	    ABsolveP4(FillMatrixRow(OneDindex(i,j), OneDindex(i+1,j),KLsolve4, KUsolve4, Nsolve4), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(kappah(i,j )   + kappah(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(kappah(i,j )   + kappah(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
	  end if

  ! 	  ! dirichlet condnitions on boundaries
  ! 	  AsolveP(OneDindex(i,N), OneDindex(i,N)) = 1d0; BsolveP4(OneDindex(i,N)) = Ne0; AsolveP(OneDindex(i,N), OneDindex(i,N-1)) = 0d0;
  ! 	  AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0; BsolveP4(OneDindex(i,1)) = Ne0; AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = 0d0;
  ! 	  AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0; BsolveP4(OneDindex(M,j)) = Ne0; AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = 0d0;
  ! 	  AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0; BsolveP4(OneDindex(1,j)) = Ne0; AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = 0d0;


      !!Ts
	  if(FillMatrixCondition(OneDindex(i,j),OneDindex(i,j),KLsolve5, KUsolve5, Nsolve5)) then
	    ABsolveP5(FillMatrixRow(OneDindex(i,j), OneDindex(i,j),KLsolve5,KUsolve5,Nsolve5),OneDindex(i,j)) = Cs(i,j)*CellVol(i,j)/dt + ( &
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j)) = CellVol(i,j)/dt + ( &
			+ 0.5d0*CellAreaE(i,j)*(kappas(i,j ) +kappas(i+1,j))*(1d0)/DistE(i,j) * (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j)) &
			+ 0.5d0*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j) ) *(1d0)/DistW(i,j) * (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j)) &
			+ 0.5d0*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j) ) *(1d0)/DistN(i,j) * (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j)) &
			+ 0.5d0*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j) ) *(1d0)/DistS(i,j) * (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j)) &
			)
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j-1),KLsolve5,KUsolve5,Nsolve5)) then
	    ABsolveP5(FillMatrixRow(OneDindex(i,j), OneDindex(i,j-1),KLsolve5, KUsolve5, Nsolve5), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(kappas(i,j-1 ) + kappas(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j-1)) =0.5d0*CellAreaS(i,j)*(kappas(i,j-1 ) + kappas(i,j)   )*( -1d0 ) / DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i,j+1),KLsolve5,KUsolve5,Nsolve5)) then
	    ABsolveP5(FillMatrixRow(OneDindex(i,j), OneDindex(i,j+1),KLsolve5, KUsolve5, Nsolve5), OneDindex(i,j+1))=0.5d0*CellAreaN(i,j)*(kappas(i,j+1 ) + kappas(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i,j+1)) =0.5d0*CellAreaN(i,j)*(kappas(i,j+1 ) + kappas(i,j)   )*( -1d0 ) / DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i-1,j),KLsolve5,KUsolve5,Nsolve5)) then
	    ABsolveP5(FillMatrixRow(OneDindex(i,j), OneDindex(i-1,j),KLsolve5, KUsolve5, Nsolve5), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(kappas(i-1,j ) + kappas(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i-1,j)) =0.5d0*CellAreaW(i,j)*(kappas(i-1,j ) + kappas(i,j)   )*( -1d0 ) / DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))
	  end if
	  if(FillMatrixCondition(OneDindex(i,j), OneDindex(i+1,j),KLsolve5,KUsolve5,Nsolve5)) then
	    ABsolveP5(FillMatrixRow(OneDindex(i,j), OneDindex(i+1,j),KLsolve5, KUsolve5, Nsolve5), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(kappas(i,j )   + kappas(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
! 	    AsolveP(OneDindex(i,j), OneDindex(i+1,j)) =0.5d0*CellAreaE(i,j)*(kappas(i,j )   + kappas(i+1,j) )*( -1d0 ) / DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))
	  end if

  ! 	  ! dirichlet condnitions on boundaries
  ! 	  AsolveP(OneDindex(i,N), OneDindex(i,N)) = 1d0; BsolveP5(OneDindex(i,N)) = Ne0; AsolveP(OneDindex(i,N), OneDindex(i,N-1)) = 0d0;
  ! 	  AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0; BsolveP5(OneDindex(i,1)) = Ne0; AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = 0d0;
  ! 	  AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0; BsolveP5(OneDindex(M,j)) = Ne0; AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = 0d0;
  ! 	  AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0; BsolveP5(OneDindex(1,j)) = Ne0; AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = 0d0;

      
	end do !j
      end do !i
      !$OMP END DO
      
      !$OMP DO
      do i=2, Mp-1
! 	write(*,*) "Complete the boundaries in i..."
    !! Ne
	  ! neumann conditions on boundaries
	if(FillMatrixCondition(OneDindex(i,Np), OneDindex(i,Np),KLsolve1,KUsolve1,Nsolve1)) then
! 	    AsolveP(OneDindex(i,Np), OneDindex(i,Np)) = 1d0;
	  ABsolveP1(FillMatrixRow(OneDindex(i,Np), OneDindex(i,Np), KLsolve1, KUsolve1, Nsolve1), OneDindex(i,Np)) = 1d0;
	end if
	if(FillMatrixCondition(OneDindex(i,Np),OneDindex(i,Np-1),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(i,Np), OneDindex(i,Np-1)) = -1d0;
	  ABsolveP1(FillMatrixRow(OneDindex(i,Np), OneDindex(i,Np-1),KLsolve1,KUsolve1,Nsolve1),OneDindex(i,Np-1)) = -1d0;
	end if
	if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,1),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0;
	  ABsolveP1(FillMatrixRow(OneDindex(i,1), OneDindex(i,1), KLsolve1, KUsolve1, Nsolve1), OneDindex(i,1)) = 1d0;
	end if
	if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,2),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = -1d0;
	  ABsolveP1(FillMatrixRow(OneDindex(i,1), OneDindex(i,2), KLsolve1, KUsolve1, Nsolve1), OneDindex(i,2)) = -1d0;
	end if

	BsolveP1(OneDindex(i,Np)) = 0d0; BsolveP1(OneDindex(i,1)) = 0d0;

    !! Nh

	! neumann conditions on boundaries
	if(FillMatrixCondition(OneDindex(i,Np), OneDindex(i,Np), KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(i,Np), OneDindex(i,Np)) = 1d0;
	  ABsolveP2(FillMatrixRow(OneDindex(i,Np), OneDindex(i,Np), KLsolve2, KUsolve2, Nsolve2), OneDindex(i,Np)) = 1d0;
	end if
	if(FillMatrixCondition(OneDindex(i,Np),OneDindex(i,N-1),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(i,Np), OneDindex(i,N-1)) = -1d0;
	  ABsolveP2(FillMatrixRow(OneDindex(i,Np), OneDindex(i,Np-1), KLsolve2, KUsolve2, Nsolve2),OneDindex(i,Np-1)) = -1d0;
	end if
	if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,1),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0;
	  ABsolveP2(FillMatrixRow(OneDindex(i,1), OneDindex(i,1), KLsolve2, KUsolve2, Nsolve2), OneDindex(i,1)) = 1d0;
	end if
	if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,2),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = -1d0;
	  ABsolveP2(FillMatrixRow(OneDindex(i,1), OneDindex(i,2), KLsolve2, KUsolve2, Nsolve2), OneDindex(i,2)) = -1d0;
	end if

	BsolveP2(OneDindex(i,Np)) = 0d0; BsolveP2(OneDindex(i,1)) = 0d0;


    !!Te

	  ! neumann conditions on boundaries
	  if(FillMatrixCondition(OneDindex(i,N), OneDindex(i,N),KLsolve3,KUsolve3,Nsolve3)) then
! 	    AsolveP(OneDindex(i,N), OneDindex(i,N)) = 1d0;
	    ABsolveP3(FillMatrixRow(OneDindex(i,N), OneDindex(i,N), KLsolve3, KUsolve3, Nsolve3), OneDindex(i,N)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,N),OneDindex(i,N-1),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(i,N), OneDindex(i,N-1)) = -1d0;
	    ABsolveP3(FillMatrixRow(OneDindex(i,N), OneDindex(i,N-1),KLsolve3,KUsolve3,Nsolve3),OneDindex(i,N-1)) = -1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,1),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0;
	    ABsolveP3(FillMatrixRow(OneDindex(i,1), OneDindex(i,1), KLsolve3, KUsolve3, Nsolve3), OneDindex(i,1)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,2),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = -1d0;
	    ABsolveP3(FillMatrixRow(OneDindex(i,1), OneDindex(i,2), KLsolve3, KUsolve3, Nsolve3), OneDindex(i,2)) = -1d0;
	  end if

	  BsolveP3(OneDindex(i,Np)) = 0d0; BsolveP3(OneDindex(i,1)) = 0d0;
	  
      !!Th
	  ! neumann conditions on boundaries
	  if(FillMatrixCondition(OneDindex(i,N), OneDindex(i,N),KLsolve4,KUsolve4,Nsolve4)) then
! 	    AsolveP(OneDindex(i,N), OneDindex(i,N)) = 1d0;
	    ABsolveP4(FillMatrixRow(OneDindex(i,N), OneDindex(i,N), KLsolve4, KUsolve4, Nsolve4), OneDindex(i,N)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,N),OneDindex(i,N-1),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(i,N), OneDindex(i,N-1)) = -1d0;
	    ABsolveP4(FillMatrixRow(OneDindex(i,N), OneDindex(i,N-1),KLsolve4,KUsolve4,Nsolve4),OneDindex(i,N-1)) = -1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,1),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0;
	    ABsolveP4(FillMatrixRow(OneDindex(i,1), OneDindex(i,1), KLsolve4, KUsolve4, Nsolve4), OneDindex(i,1)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,2),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = -1d0;
	    ABsolveP4(FillMatrixRow(OneDindex(i,1), OneDindex(i,2), KLsolve4, KUsolve4, Nsolve4), OneDindex(i,2)) = -1d0;
	  end if

	  BsolveP4(OneDindex(i,N)) = 0d0; BsolveP4(OneDindex(i,1)) = 0d0;

      !! Ts
	    ! neumann conditions on boundaries
	  if(FillMatrixCondition(OneDindex(i,N), OneDindex(i,N),KLsolve5,KUsolve5,Nsolve5)) then
! 	    AsolveP(OneDindex(i,N), OneDindex(i,N)) = 1d0;
	    ABsolveP5(FillMatrixRow(OneDindex(i,N), OneDindex(i,N), KLsolve5, KUsolve5, Nsolve5), OneDindex(i,N)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,N),OneDindex(i,N-1),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(i,N), OneDindex(i,N-1)) = -1d0;
	    ABsolveP5(FillMatrixRow(OneDindex(i,N), OneDindex(i,N-1),KLsolve5,KUsolve5,Nsolve5),OneDindex(i,N-1)) = -1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,1),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,1)) = 1d0;
	    ABsolveP5(FillMatrixRow(OneDindex(i,1), OneDindex(i,1), KLsolve5, KUsolve5, Nsolve5), OneDindex(i,1)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(i,1),OneDindex(i,2),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(i,1), OneDindex(i,2 ) ) = -1d0;
	    ABsolveP5(FillMatrixRow(OneDindex(i,1), OneDindex(i,2), KLsolve5, KUsolve5, Nsolve5), OneDindex(i,2)) = -1d0;
	  end if

	  BsolveP5(OneDindex(i,N)) = 0d0; BsolveP5(OneDindex(i,1)) = 0d0;
	
      end do !i
      !$OMP END DO

      
!       write(*,*) "Complete the boundaries in j..."
      
      !$OMP DO
      do j=2, Np-1

    !!Ne
	if(FillMatrixCondition(OneDindex(Mp,j),OneDindex(Mp,j),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0;
	  ABsolveP1(FillMatrixRow(OneDindex(Mp,j), OneDindex(Mp,j), KLsolve1, KUsolve1, Nsolve1), OneDindex(Mp,j)) = 1d0;
	end if
	if(FillMatrixCondition(OneDindex(Mp,j),OneDindex(Mp-1,j),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = -1d0;
	  ABsolveP1(FillMatrixRow(OneDindex(Mp,j), OneDindex(Mp-1,j), KLsolve1, KUsolve1, Nsolve1), OneDindex(Mp-1,j)) = -1d0;
	end if
	if(FillMatrixCondition(OneDindex(1,j),OneDindex(1,j),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0;
	  ABsolveP1(FillMatrixRow(OneDindex(1,j), OneDindex(1,j), KLsolve1, KUsolve1, Nsolve1), OneDindex(1,j)) = 1d0;
	end if
	if(FillMatrixCondition(OneDindex(1,j),OneDindex(2,j),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = -1d0;
	  ABsolveP1(FillMatrixRow(OneDindex(1,j), OneDindex(2,j), KLsolve1, KUsolve1, Nsolve1), OneDindex(2,j)) = -1d0;
	end if

	BsolveP1(OneDindex(Mp,j)) = 0d0; BsolveP1(OneDindex(1,j)) = 0d0;

    !! Nh

	if(FillMatrixCondition(OneDindex(Mp,j),OneDindex(Mp,j),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(Mp,j), OneDindex(Mp,j)) = 1d0;
	  ABsolveP2(FillMatrixRow(OneDindex(Mp,j), OneDindex(Mp,j), KLsolve2, KUsolve2, Nsolve2), OneDindex(Mp,j)) = 1d0;
	end if
	if(FillMatrixCondition(OneDindex(Mp,j),OneDindex(Mp-1,j),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(Mp,j), OneDindex(Mp-1,j)) = -1d0;
	  ABsolveP2(FillMatrixRow(OneDindex(Mp,j), OneDindex(Mp-1,j), KLsolve2, KUsolve2, Nsolve2), OneDindex(Mp-1,j)) = -1d0;
	end if
	if(FillMatrixCondition(OneDindex(1,j),OneDindex(1,j),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0;
	  ABsolveP2(FillMatrixRow(OneDindex(1,j), OneDindex(1,j), KLsolve2, KUsolve2, Nsolve2), OneDindex(1,j)) = 1d0;
	end if
	if(FillMatrixCondition(OneDindex(1,j),OneDindex(2,j),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = -1d0;
	  ABsolveP2(FillMatrixRow(OneDindex(1,j), OneDindex(2,j), KLsolve2, KUsolve2, Nsolve2), OneDindex(2,j)) = -1d0;
	end if

	BsolveP2(OneDindex(Mp,j)) = 0d0; BsolveP2(OneDindex(1,j)) = 0d0;


    !! Te

	if(FillMatrixCondition(OneDindex(M,j),OneDindex(M,j),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0;
	    ABsolveP3(FillMatrixRow(OneDindex(M,j), OneDindex(M,j), KLsolve3, KUsolve3, Nsolve3), OneDindex(M,j)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(M,j),OneDindex(M-1,j),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = -1d0;
	    ABsolveP3(FillMatrixRow(OneDindex(M,j), OneDindex(M-1,j), KLsolve3, KUsolve3, Nsolve3), OneDindex(M-1,j)) = -1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(1,j),OneDindex(1,j),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0;
	    ABsolveP3(FillMatrixRow(OneDindex(1,j), OneDindex(1,j), KLsolve3, KUsolve3, Nsolve3), OneDindex(1,j)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(1,j),OneDindex(2,j),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = -1d0;
	    ABsolveP3(FillMatrixRow(OneDindex(1,j), OneDindex(2,j), KLsolve3, KUsolve3, Nsolve3), OneDindex(2,j)) = -1d0;
	  end if

	  BsolveP3(OneDindex(M,j)) = 0d0; BsolveP3(OneDindex(1,j)) = 0d0;
	  
      !!Th
	  if(FillMatrixCondition(OneDindex(M,j),OneDindex(M,j),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0;
	    ABsolveP4(FillMatrixRow(OneDindex(M,j), OneDindex(M,j), KLsolve4, KUsolve4, Nsolve4), OneDindex(M,j)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(M,j),OneDindex(M-1,j),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = -1d0;
	    ABsolveP4(FillMatrixRow(OneDindex(M,j), OneDindex(M-1,j), KLsolve4, KUsolve4, Nsolve4), OneDindex(M-1,j)) = -1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(1,j),OneDindex(1,j),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0;
	    ABsolveP4(FillMatrixRow(OneDindex(1,j), OneDindex(1,j), KLsolve4, KUsolve4, Nsolve4), OneDindex(1,j)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(1,j),OneDindex(2,j),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = -1d0;
	    ABsolveP4(FillMatrixRow(OneDindex(1,j), OneDindex(2,j), KLsolve4, KUsolve4, Nsolve4), OneDindex(2,j)) = -1d0;
	  end if

	  BsolveP4(OneDindex(M,j)) = 0d0; BsolveP4(OneDindex(1,j)) = 0d0;

      !! Ts
	  if(FillMatrixCondition(OneDindex(M,j),OneDindex(M,j),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(M,j), OneDindex(M,j)) = 1d0;
	    ABsolveP5(FillMatrixRow(OneDindex(M,j), OneDindex(M,j), KLsolve5, KUsolve5, Nsolve5), OneDindex(M,j)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(M,j),OneDindex(M-1,j),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(M,j), OneDindex(M-1,j)) = -1d0;
	    ABsolveP5(FillMatrixRow(OneDindex(M,j), OneDindex(M-1,j), KLsolve5, KUsolve5, Nsolve5), OneDindex(M-1,j)) = -1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(1,j),OneDindex(1,j),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(1,j)) = 1d0;
	    ABsolveP5(FillMatrixRow(OneDindex(1,j), OneDindex(1,j), KLsolve5, KUsolve5, Nsolve5), OneDindex(1,j)) = 1d0;
	  end if
	  if(FillMatrixCondition(OneDindex(1,j),OneDindex(2,j),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(1,j), OneDindex(2,j ) ) = -1d0;
	    ABsolveP5(FillMatrixRow(OneDindex(1,j), OneDindex(2,j), KLsolve5, KUsolve5, Nsolve5), OneDindex(2,j)) = -1d0;
	  end if

	  BsolveP5(OneDindex(M,j)) = 0d0; BsolveP5(OneDindex(1,j)) = 0d0;
	  
	  
      end do !j
      !$OMP END DO


      
      
      !$OMP DO COLLAPSE(2)
      do i=1,Mp
	do j=1,Np
      !! Ne
	  BsolveP1(OneDindex(i,j)) = CellVol(i,j)*Ne(i,j)/dt  + (GainsE(i,j)-LossesE(i,j))*CellVol(i,j) &
		+ CrossCoeff*( &
		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(diffusionE(i,j)+diffusionE(i+1,j))*( NeDual(i,j) - NeDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(diffusionE(i-1,j)+diffusionE(i,j))*( NeDual(i-1,j-1) - NeDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(diffusionE(i,j+1)+diffusionE(i,j))*( NeDual(i-1,j) - NeDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(diffusionE(i,j-1)+diffusionE(i,j))*( NeDual(i,j-1) - NeDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
		  )

      !! Nh
      	  BsolveP2(OneDindex(i,j)) = CellVol(i,j)*Nh(i,j)/dt  + (GainsH(i,j)-LossesH(i,j))*CellVol(i,j) &
		+ CrossCoeff*( &
		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(diffusionH(i,j)+diffusionH(i+1,j))*( NhDual(i,j) - NhDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(diffusionH(i-1,j)+diffusionH(i,j))*( NhDual(i-1,j-1) - NhDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(diffusionH(i,j+1)+diffusionH(i,j))*( NhDual(i-1,j) - NhDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(diffusionH(i,j-1)+diffusionH(i,j))*( NhDual(i,j-1) - NhDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
		  )

      !! Te
	  XsolveP(OneDindex(i,j)) = Ne(i,j)
	  BsolveP3(OneDindex(i,j)) = CellVol(i,j)*Ce(i,j)*Te(i,j)/dt  + (SourceE(i,j)-CouplingE(i,j))*CellVol(i,j) &
		+ CrossCoeff*( &
		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappae(i,j)+kappae(i+1,j))*( TeDual(i,j) - TeDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappae(i-1,j)+kappae(i,j))*( TeDual(i-1,j-1) - TeDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappae(i,j+1)+kappae(i,j))*( TeDual(i-1,j) - TeDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappae(i,j-1)+kappae(i,j))*( TeDual(i,j-1) - TeDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
		  )

      !! Th
! 	  XsolveP(OneDindex(i,j)) = Ne(i,j)
	  BsolveP4(OneDindex(i,j)) = CellVol(i,j)*Ch(i,j)*Th(i,j)/dt  + (SourceH(i,j)-CouplingH(i,j))*CellVol(i,j) &
		+ CrossCoeff*( &
		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappah(i,j)+kappah(i+1,j))*( ThDual(i,j) - ThDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappah(i-1,j)+kappah(i,j))*( ThDual(i-1,j-1) - ThDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappah(i,j+1)+kappah(i,j))*( ThDual(i-1,j) - ThDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappah(i,j-1)+kappah(i,j))*( ThDual(i,j-1) - ThDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
		  )
      !! Ts

! 	  XsolveP(OneDindex(i,j)) = Ne(i,j)
	  BsolveP5(OneDindex(i,j)) = CellVol(i,j)*Cs(i,j)*Ts(i,j)/dt  + (CouplingE(i,j)+CouplingH(i,j))*CellVol(i,j) &
		+ CrossCoeff*( &
		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*( TsDual(i,j) - TsDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*( TsDual(i-1,j-1) - TsDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*( TsDual(i-1,j) - TsDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*( TsDual(i,j-1) - TsDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
		  )
	end do
      end do
      !$OMP END DO

      !$OMP DO
      do i=1, Mp

      !!Ne
	  ! for Dirichlet conditions
! 	  BsolveP(OneDindex(i,1)) = Ne0; BsolveP(OneDindex(i,N)) = Ne0;
	  ! for Neumann conditions
	  BsolveP1(OneDindex(i,1)) = 0d0; BsolveP1(OneDindex(i,Np)) = 0d0

      !!Nh
	  ! for Dirichlet conditions
! 	  BsolveP2(OneDindex(i,1)) = Ne0; BsolveP2(OneDindex(i,N)) = Ne0;
 	  ! for Neumann conditions
	  BsolveP2(OneDindex(i,1)) = 0d0; BsolveP2(OneDindex(i,Np)) = 0d0

      !! Te
	  ! for Dirichlet conditions
! 	  BsolveP3(OneDindex(i,1)) = Ne0; BsolveP3(OneDindex(i,Np)) = Ne0;
 	  ! for Neumann conditions
 	  BsolveP3(OneDindex(i,1)) = 0d0; BsolveP3(OneDindex(i,Np)) = 0d0

      !! Th
	  ! for Dirichlet conditions
! 	  BsolveP4(OneDindex(i,1)) = Ne0; BsolveP4(OneDindex(i,Np)) = Ne0;
	  ! for Neumann conditions
	  BsolveP4(OneDindex(i,1)) = 0d0; BsolveP4(OneDindex(i,Np)) = 0d0

      !! Ts
   	  ! for Dirichlet conditions
! 	  BsolveP5(OneDindex(i,1)) = Ne0; BsolveP5(OneDindex(i,N)) = Ne0;
 	  ! for Neumann conditions
 	  BsolveP5(OneDindex(i,1)) = 0d0; BsolveP5(OneDindex(i,Np)) = 0d0

      end do
      !$OMP END DO

      !$OMP DO
      do j=1, Np
    !! Ne
! 	! for Dirichlet conditions
! 	BsolveP1(OneDindex(1,j)) = Ne0; BsolveP1(OneDindex(M,j)) = Ne0;
	! for Neumann conditions
	BsolveP1(OneDindex(1,j)) = 0d0; BsolveP1(OneDindex(Mp,j)) = 0d0;
    !! Nh
	! for Dirichlet conditions
! 	BsolveP2(OneDindex(1,j)) = Ne0; BsolveP2(OneDindex(M,j)) = Ne0;
	! for Neumann conditions
	BsolveP2(OneDindex(1,j)) = 0d0; BsolveP2(OneDindex(Mp,j)) = 0d0

    !!Te
	! for Dirichlet conditions
! 	BsolveP3(OneDindex(1,j)) = Ne0; BsolveP3(OneDindex(Mp,j)) = Ne0;
	! for Neumann conditions
	BsolveP3(OneDindex(1,j)) = 0d0; BsolveP3(OneDindex(Mp,j)) = 0d0;

    !!Th
	! for Dirichlet conditions
!  	BsolveP4(OneDindex(1,j)) = Ne0; BsolveP4(OneDindex(M,j)) = Ne0;
 	! for Neumann conditions
 	BsolveP4(OneDindex(1,j)) = 0d0; BsolveP4(OneDindex(Mp,j)) = 0d0;

    !! Ts
	! for Dirichlet conditions
	BsolveP5(OneDindex(1,j)) = Ne0; BsolveP5(OneDindex(Mp,j)) = Ne0;
 	! for Neumann conditions
 	BsolveP5(OneDindex(1,j)) = 0d0; BsolveP5(OneDindex(Mp,j)) = 0d0
 	
      end do
      !$OMP END DO

      !$OMP DO COLLAPSE(2)
      do i=1,Mp
	do j=1,Np

      !!Ne
	 ! corners dont play any role then must be deleted
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(i,j),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(1,1), OneDindex(i,j))=0d0;
	    ABsolveP1(FillMatrixRow(OneDindex(1,1), OneDindex(i,j),KLsolve1,KUsolve1,Nsolve1),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,Np),OneDindex(i,j),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(1,Np), OneDindex(i,j))=0d0
	    ABsolveP1(FillMatrixRow(OneDindex(1,Np), OneDindex(i,j),KLsolve1,KUsolve1,Nsolve1),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,1),OneDindex(i,j),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(Mp,1), OneDindex(i,j))=0d0
	    ABsolveP1(FillMatrixRow(OneDindex(Mp,1), OneDindex(i,j),KLsolve1,KUsolve1,Nsolve1),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,Np),OneDindex(i,j),KLsolve1, KUsolve1, Nsolve1)) then
! 	    AsolveP(OneDindex(Mp,Np), OneDindex(i,j))=0d0
	    ABsolveP1(FillMatrixRow(OneDindex(Mp,Np), OneDindex(i,j),KLsolve1,KUsolve1,Nsolve1),OneDindex(i,j))=0d0
	  end if
	  
      !! Nh
	  ! corners dont play any role then must be deleted
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(i,j),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(1,1), OneDindex(i,j))=0d0;
	    ABsolveP2(FillMatrixRow(OneDindex(1,1), OneDindex(i,j),KLsolve2,KUsolve2,Nsolve2),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,Np),OneDindex(i,j),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(1,Np), OneDindex(i,j))=0d0
	    ABsolveP2(FillMatrixRow(OneDindex(1,Np), OneDindex(i,j),KLsolve2,KUsolve2,Nsolve2),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,1),OneDindex(i,j),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(Mp,1), OneDindex(i,j))=0d0
	    ABsolveP2(FillMatrixRow(OneDindex(Mp,1), OneDindex(i,j),KLsolve2,KUsolve2,Nsolve2),OneDindex(i,j))=0d0
	  end if

	  if(FillMatrixCondition(OneDindex(Mp,Np),OneDindex(i,j),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(Mp,Np), OneDindex(i,j))=0d0
	    ABsolveP2(FillMatrixRow(OneDindex(Mp,Np), OneDindex(i,j),KLsolve2,KUsolve2,Nsolve2),OneDindex(i,j))=0d0
	  end if

      !! Te
	  ! corners dont play any role then must be deleted
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(i,j),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(1,1), OneDindex(i,j))=0d0;
	    ABsolveP3(FillMatrixRow(OneDindex(1,1), OneDindex(i,j),KLsolve3,KUsolve3,Nsolve3),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,N),OneDindex(i,j),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(1,N), OneDindex(i,j))=0d0
	    ABsolveP3(FillMatrixRow(OneDindex(1,N), OneDindex(i,j),KLsolve3,KUsolve3,Nsolve3),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,1),OneDindex(i,j),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(M,1), OneDindex(i,j))=0d0
	    ABsolveP3(FillMatrixRow(OneDindex(M,1), OneDindex(i,j),KLsolve3,KUsolve3,Nsolve3),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,N),OneDindex(i,j),KLsolve3, KUsolve3, Nsolve3)) then
! 	    AsolveP(OneDindex(M,N), OneDindex(i,j))=0d0
	    ABsolveP3(FillMatrixRow(OneDindex(M,N), OneDindex(i,j),KLsolve3,KUsolve3,Nsolve3),OneDindex(i,j))=0d0
	  end if

      !!Th
	  ! corners dont play any role then must be deleted
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(i,j),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(1,1), OneDindex(i,j))=0d0;
	    ABsolveP4(FillMatrixRow(OneDindex(1,1), OneDindex(i,j),KLsolve4,KUsolve4,Nsolve4),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,N),OneDindex(i,j),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(1,N), OneDindex(i,j))=0d0
	    ABsolveP4(FillMatrixRow(OneDindex(1,N), OneDindex(i,j),KLsolve4,KUsolve4,Nsolve4),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,1),OneDindex(i,j),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(M,1), OneDindex(i,j))=0d0
	    ABsolveP4(FillMatrixRow(OneDindex(M,1), OneDindex(i,j),KLsolve4,KUsolve4,Nsolve4),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,N),OneDindex(i,j),KLsolve4, KUsolve4, Nsolve4)) then
! 	    AsolveP(OneDindex(M,N), OneDindex(i,j))=0d0
	    ABsolveP4(FillMatrixRow(OneDindex(M,N), OneDindex(i,j),KLsolve4,KUsolve4,Nsolve4),OneDindex(i,j))=0d0
	  end if


      !! Ts
	  ! corners dont play any role then must be deleted
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(i,j),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(1,1), OneDindex(i,j))=0d0;
	    ABsolveP5(FillMatrixRow(OneDindex(1,1), OneDindex(i,j),KLsolve5,KUsolve5,Nsolve5),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,N),OneDindex(i,j),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(1,N), OneDindex(i,j))=0d0
	    ABsolveP5(FillMatrixRow(OneDindex(1,N), OneDindex(i,j),KLsolve5,KUsolve5,Nsolve5),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,1),OneDindex(i,j),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(M,1), OneDindex(i,j))=0d0
	    ABsolveP5(FillMatrixRow(OneDindex(M,1), OneDindex(i,j),KLsolve5,KUsolve5,Nsolve5),OneDindex(i,j))=0d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,N),OneDindex(i,j),KLsolve5, KUsolve5, Nsolve5)) then
! 	    AsolveP(OneDindex(M,N), OneDindex(i,j))=0d0
	    ABsolveP5(FillMatrixRow(OneDindex(M,N), OneDindex(i,j),KLsolve5,KUsolve5,Nsolve5),OneDindex(i,j))=0d0
	  end if
	end do
      end do
      !$OMP END DO
	  
	  ! other direction, no ?
! 	  AsolveP(OneDindex(i,j), OneDindex(1,1))=0d0; AsolveP(OneDindex(i,j), OneDindex(1,Np))=0d0
! 	  AsolveP(OneDindex(i,j), OneDindex(Mp,1))=0d0; AsolveP(OneDindex(i,j), OneDindex(M,Np))=0d0
    !! Ne
	  ! however, diagonal cannot be null
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(1,1),KLsolve1, KUsolve1, Nsolve1)) then
  ! 	    AsolveP(OneDindex(1,1), OneDindex(1,1))=1d0
	    ABsolveP1(FillMatrixRow(OneDindex(1,1), OneDindex(1,1),KLsolve1,KUsolve1,Nsolve1),OneDindex(1,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,Np),OneDindex(1,Np),KLsolve1, KUsolve1, Nsolve1)) then
  ! 	    AsolveP(OneDindex(1,Np), OneDindex(1,Np))=1d0
	    ABsolveP1(FillMatrixRow(OneDindex(1,Np), OneDindex(1,Np),KLsolve1,KUsolve1,Nsolve1),OneDindex(1,Np))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,1),OneDindex(Mp,1),KLsolve1, KUsolve1, Nsolve1)) then
  ! 	    AsolveP(OneDindex(Mp,1), OneDindex(Mp,1))=1d0
	    ABsolveP1(FillMatrixRow(OneDindex(Mp,1), OneDindex(Mp,1),KLsolve1,KUsolve1,Nsolve1),OneDindex(Mp,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,Np),OneDindex(Mp,Np),KLsolve1, KUsolve1, Nsolve1)) then
  ! 	    AsolveP(OneDindex(Mp,Np), OneDindex(Mp,Np))=1d0
	    ABsolveP1(FillMatrixRow(OneDindex(Mp,Np), OneDindex(Mp,Np),KLsolve1,KUsolve1,Nsolve1),OneDindex(Mp,Np))=1d0
	  end if

	  ! let give a value to the corners, since they are used for NeDual calculation
  ! 	  BsolveP1(OneDindex(1,1))=Ne0; BsolveP1(OneDindex(1,Np))=Ne0
  ! 	  BsolveP1(OneDindex(Mp,1))=Ne0; BsolveP1(OneDindex(Mp,Np))=Ne0
	  BsolveP1(OneDindex(1,1))=0.5d0*(Ne(1,2)+Ne(2,1)); BsolveP1(OneDindex(1,Np))=0.5d0*(Ne(1,Np-1)+Ne(2,Np))
	  BsolveP1(OneDindex(Mp,1))=0.5d0*(Ne(Mp,2)+Ne(Mp-1,1)); BsolveP1(OneDindex(Mp,Np))=0.5d0*(Ne(Mp,Np-1)+Ne(Mp-1,Np))

    !! Nh
	  ! however, diagonal cannot be null
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(1,1),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(1,1), OneDindex(1,1))=1d0
	    ABsolveP2(FillMatrixRow(OneDindex(1,1), OneDindex(1,1),KLsolve2,KUsolve2,Nsolve2),OneDindex(1,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,Np),OneDindex(1,Np),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(1,Np), OneDindex(1,Np))=1d0
	    ABsolveP2(FillMatrixRow(OneDindex(1,Np), OneDindex(1,Np),KLsolve2,KUsolve2,Nsolve2),OneDindex(1,Np))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,1),OneDindex(Mp,1),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(Mp,1), OneDindex(Mp,1))=1d0
	    ABsolveP2(FillMatrixRow(OneDindex(Mp,1), OneDindex(Mp,1),KLsolve2,KUsolve2,Nsolve2),OneDindex(Mp,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,Np),OneDindex(Mp,Np),KLsolve2, KUsolve2, Nsolve2)) then
! 	    AsolveP(OneDindex(Mp,Np), OneDindex(Mp,Np))=1d0
	    ABsolveP2(FillMatrixRow(OneDindex(Mp,Np), OneDindex(Mp,Np),KLsolve2,KUsolve2,Nsolve2),OneDindex(Mp,Np))=1d0
	  end if

	 ! let give a value to the corners, since they are used for NeDual calculation
	  BsolveP2(OneDindex(1,1))=0.5d0*(Nh(1,2)+Nh(2,1)); BsolveP2(OneDindex(1,Np))=0.5d0*(Nh(1,Np-1)+Nh(2,Np))
	  BsolveP2(OneDindex(Mp,1))=0.5d0*(Nh(Mp,2)+Nh(Mp-1,1)); BsolveP2(OneDindex(Mp,Np))=0.5d0*(Nh(Mp,Np-1)+Nh(Mp-1,Np))

    !!Te
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(1,1),KLsolve3, KUsolve3, Nsolve3)) then
    ! 	    AsolveP(OneDindex(1,1), OneDindex(1,1))=1d0
	    ABsolveP3(FillMatrixRow(OneDindex(1,1), OneDindex(1,1),KLsolve3,KUsolve3,Nsolve3),OneDindex(1,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,N),OneDindex(1,N),KLsolve3, KUsolve3, Nsolve3)) then
    ! 	    AsolveP(OneDindex(1,N), OneDindex(1,N))=1d0
	    ABsolveP3(FillMatrixRow(OneDindex(1,N), OneDindex(1,N),KLsolve3,KUsolve3,Nsolve3),OneDindex(1,N))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,1),OneDindex(M,1),KLsolve3, KUsolve3, Nsolve3)) then
    ! 	    AsolveP(OneDindex(M,1), OneDindex(M,1))=1d0
	    ABsolveP3(FillMatrixRow(OneDindex(M,1), OneDindex(M,1),KLsolve3,KUsolve3,Nsolve3),OneDindex(M,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,N),OneDindex(M,N),KLsolve3, KUsolve3, Nsolve3)) then
    ! 	    AsolveP(OneDindex(M,N), OneDindex(M,N))=1d0
	    ABsolveP3(FillMatrixRow(OneDindex(M,N), OneDindex(M,N),KLsolve3,KUsolve3,Nsolve3),OneDindex(M,N))=1d0
	  end if

	  ! let give a value to the corners, since they are used for NeDual calculation
	  BsolveP3(OneDindex(1,1))=0.5d0*(Te(1,2)+Te(2,1)); BsolveP3(OneDindex(1,N))=0.5d0*(Te(1,N-1)+Te(2,N))
	  BsolveP3(OneDindex(M,1))=0.5d0*(Te(M,2)+Te(M-1,1)); BsolveP3(OneDindex(M,N))=0.5d0*(Te(M,N-1)+Te(M-1,N))
    !! Th

	  ! however, diagonal cannot be null
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(1,1),KLsolve4, KUsolve4, Nsolve4)) then
    ! 	    AsolveP(OneDindex(1,1), OneDindex(1,1))=1d0
	    ABsolveP4(FillMatrixRow(OneDindex(1,1), OneDindex(1,1),KLsolve4,KUsolve4,Nsolve4),OneDindex(1,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,Np),OneDindex(1,Np),KLsolve4, KUsolve4, Nsolve4)) then
    ! 	    AsolveP(OneDindex(1,Np), OneDindex(1,Np))=1d0
	    ABsolveP4(FillMatrixRow(OneDindex(1,Np), OneDindex(1,Np),KLsolve4,KUsolve4,Nsolve4),OneDindex(1,Np))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,1),OneDindex(M,1),KLsolve4, KUsolve4, Nsolve4)) then
    ! 	    AsolveP(OneDindex(M,1), OneDindex(M,1))=1d0
	    ABsolveP4(FillMatrixRow(OneDindex(M,1), OneDindex(M,1),KLsolve4,KUsolve4,Nsolve4),OneDindex(M,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(M,Np),OneDindex(M,Np),KLsolve4, KUsolve4, Nsolve4)) then
    ! 	    AsolveP(OneDindex(M,Np), OneDindex(M,Np))=1d0
	    ABsolveP4(FillMatrixRow(OneDindex(M,Np), OneDindex(M,Np),KLsolve4,KUsolve4,Nsolve4),OneDindex(M,Np))=1d0
	  end if

	    ! let give a value to the corners, since they are used for NeDual calculation
	  BsolveP4(OneDindex(1,1))=0.5d0*(Th(1,2)+Th(2,1)); BsolveP4(OneDindex(1,Np))=0.5d0*(Th(1,Np-1)+Th(2,Np))
	  BsolveP4(OneDindex(Mp,1))=0.5d0*(Th(Mp,2)+Th(Mp-1,1)); BsolveP4(OneDindex(Mp,Np))=0.5d0*(Th(Mp,Np-1)+Th(Mp-1,Np))

      !! Ts
	  ! however, diagonal cannot be null
	  if(FillMatrixCondition(OneDindex(1,1),OneDindex(1,1),KLsolve5, KUsolve5, Nsolve5)) then
    ! 	    AsolveP(OneDindex(1,1), OneDindex(1,1))=1d0
	    ABsolveP5(FillMatrixRow(OneDindex(1,1), OneDindex(1,1),KLsolve5,KUsolve5,Nsolve5),OneDindex(1,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(1,Np),OneDindex(1,Np),KLsolve5, KUsolve5, Nsolve5)) then
    ! 	    AsolveP(OneDindex(1,Np), OneDindex(1,Np))=1d0
	    ABsolveP5(FillMatrixRow(OneDindex(1,Np), OneDindex(1,Np),KLsolve5,KUsolve5,Nsolve5),OneDindex(1,Np))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,1),OneDindex(Mp,1),KLsolve5, KUsolve5, Nsolve5)) then
    ! 	    AsolveP(OneDindex(Mp,1), OneDindex(Mp,1))=1d0
	    ABsolveP5(FillMatrixRow(OneDindex(Mp,1), OneDindex(Mp,1),KLsolve5,KUsolve5,Nsolve5),OneDindex(Mp,1))=1d0
	  end if
	  if(FillMatrixCondition(OneDindex(Mp,Np),OneDindex(Mp,Np),KLsolve5, KUsolve5, Nsolve5)) then
    ! 	    AsolveP(OneDindex(Mp,Np), OneDindex(Mp,Np))=1d0
	    ABsolveP5(FillMatrixRow(OneDindex(Mp,Np), OneDindex(Mp,Np),KLsolve5,KUsolve5,Nsolve5),OneDindex(Mp,Np))=1d0
	  end if

	    ! let give a value to the corners, since they are used for NeDual calculation
	  BsolveP5(OneDindex(1,1))=0.5d0*(Ts(1,2)+Ts(2,1)); BsolveP5(OneDindex(1,Np))=0.5d0*(Ts(1,Np-1)+Ts(2,Np))
	  BsolveP5(OneDindex(Mp,1))=0.5d0*(Ts(Mp,2)+Ts(Mp-1,1)); BsolveP5(OneDindex(Mp,Np))=0.5d0*(Ts(Mp,Np-1)+Ts(Mp-1,Np))

	  
  !$OMP END PARALLEL
	

! 	!!!Test of the optimised solutions (idiot, but good to check validity of the solver)
! 	call SolveSparse(AsolveP, XsolveP, BsolveP, Nsolve, Nsolve)
!!! SOLVE THE SYSTEMS

      
    !! Ne
	CALL dgbtrf(Mp*Np, Mp*Np, KLsolve1, KUsolve1, ABsolveP1, 2*KLsolve1+KUsolve1+1, ipiv1, info1)
	CALL dgbtrs('N', Nsolve1, KLsolve1, KUsolve1, 1, ABsolveP1, 2*KLsolve1+KUsolve1+1, ipiv1, BsolveP1, Nsolve1, info1)
    !! Nh
	CALL dgbtrf(Mp*Np, Mp*Np, KLsolve2, KUsolve2, ABsolveP2, 2*KLsolve2+KUsolve2+1, ipiv2, info2)
	CALL dgbtrs('N', Nsolve2, KLsolve2, KUsolve2, 1, ABsolveP2, 2*KLsolve2+KUsolve2+1, ipiv2, BsolveP2, Nsolve2, info2)
    !! Te
	CALL dgbtrf(Mp*Np, Mp*Np, KLsolve3, KUsolve3, ABsolveP3, 2*KLsolve3+KUsolve3+1, ipiv3, info3)
	CALL dgbtrs('N', Nsolve3, KLsolve3, KUsolve3, 1, ABsolveP3, 2*KLsolve3+KUsolve3+1, ipiv3, BsolveP3, Nsolve3, info3)
    !! Th
	CALL dgbtrf(Mp*Np, Mp*Np, KLsolve4, KUsolve4, ABsolveP4, 2*KLsolve4+KUsolve4+1, ipiv4, info4)
	CALL dgbtrs('N', Nsolve4, KLsolve4, KUsolve4, 1, ABsolveP4, 2*KLsolve4+KUsolve4+1, ipiv4, BsolveP4, Nsolve4, info4)
    !! Ts
	CALL dgbtrf(Mp*Np, Mp*Np, KLsolve5, KUsolve5, ABsolveP5, 2*KLsolve5+KUsolve5+1, ipiv5, info5)
	CALL dgbtrs('N', Nsolve5, KLsolve5, KUsolve5, 1, ABsolveP5, 2*KLsolve5+KUsolve5+1, ipiv5, BsolveP5, Nsolve5, info5)
	
	if ((info1 .ne. 0) .OR. (info2 .ne. 0) .OR. (info3 .ne. 0) .OR. (info4 .ne. 0) .OR. (info5 .ne. 0)) then
	  stop 'Matrix is numerically singular!'
	end if
! 
!!! Affecting the solutions

      !$OMP PARALLEL SHARED(BsolveP1, BsolveP2, BsolveP3, BsolveP4, BsolveP5, &
      !$OMP& NeNew, NhNew, TeNew, ThNew, TsNew) FIRSTPRIVATE(Mp, Np)
      
      !$OMP DO COLLAPSE(2)
      do i=1,Mp
	do j=1,Np
	  NeNew(i,j) = BsolveP1(OneDindex(i,j))
	  NhNew(i,j) = BsolveP2(OneDindex(i,j))
	  TeNew(i,j) = BsolveP3(OneDindex(i,j))
	  ThNew(i,j) = BsolveP4(OneDindex(i,j))
	  TsNew(i,j) = BsolveP5(OneDindex(i,j))
	end do
      end do
      
      !$OMP END DO
      !$OMP END PARALLEL

      !$O M P END SECTIONS
      
!       	if(mod(nbiter,iterOut).eq.0) then
! 		  write(91,'(5E18.9)') BsolveP !Amatrix
! 		  write(91,*)
! 	    end if
! 	    
      ! Checking errors
!       ErrorVec=matmul(AsolveP,XsolveP)-BsolveP
!       ErrorSum=0d0
!       do i=1,Mp*Np
! 	ErrorSum=ErrorSum+ErrorVec(i)**2
!       end do
!       write(*,*) "Norm", ErrorSum
! end if !temporary !!!      

!============ SOLVE Ts ==================


      if(TeOff.eq.0 .AND. TsOff.eq.0) then



	if (info5 .ne. 0) then
	  stop 'Matrix is numerically singular!'
	end if
      end if

      !=======================================================================================
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EXPLICIT SOLVERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! validated, works, but limited to low precision, and short max time on 8 cores x 24 h. 
      !=======================================================================================
      !$OMP DO
      do i=2, Mp-1
	!$OMP PARALLEL DO
	do j=2, Np-1
! 		 
! 	 NeNew(i,j) = (&
! 	      (GainsE(i,j)-LossesE(i,j))*CellVol(i,j) & 
! ! 	      ! convection
! ! 	      -((0.5d0*(JeX(i+1,j)+JeX(i,j))*NormalEx(i,j)+0.5d0*(JeY(i+1,j) & 
! ! 	      +JeY(i,j))*NormalEy(i,j))*CellAreaE(i,j)+(0.5d0*(JeX(i,j)+JeX(i-1,j))*NormalWx(i,j)+0.5d0*(JeY(i,j) &
! ! 	      +JeY(i-1,j))*NormalWy(i,j))*CellAreaW(i,j)+(0.5d0*(JeX(i,j)+JeX(i,j+1))*NormalNx(i,j)+0.5d0*(JeY(i,j) & 
! ! 	      +JeY(i,j+1))*NormalNy(i,j))*CellAreaN(i,j)+(0.5d0*(JeX(i,j)+JeX(i,j-1))*NormalSx(i,j)+0.5d0*(JeY(i,j) & 
! ! 	      +JeY(i,j-1))*NormalSy(i,j))*CellAreaS(i,j)) &
! 	      ! direct diffusion operator over irregular mesh
! 		  +( &
! 		  + 0.5d0*CellAreaE(i,j)*(diffusionE(i,j)+diffusionE(i+1,j))*(Ne(i+1,j)-Ne(i,j))/DistE(i,j)* (NormalEx(i,j)**2+NormalEy(i,j)**2) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j)) &
! 		  - 0.5d0*CellAreaW(i,j)*(diffusionE(i-1,j)+diffusionE(i,j))*(Ne(i,j)-Ne(i-1,j))/DistW(i,j)* (NormalWx(i,j)**2+NormalWy(i,j)**2) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j)) &
! 		  + 0.5d0*CellAreaN(i,j)*(diffusionE(i,j+1)+diffusionE(i,j))*(Ne(i,j+1)-Ne(i,j))/DistN(i,j)* (NormalNx(i,j)**2+NormalNy(i,j)**2) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j)) & 
! 		  - 0.5d0*CellAreaS(i,j)*(diffusionE(i,j-1)+diffusionE(i,j))*(Ne(i,j)-Ne(i,j-1))/DistS(i,j)* (NormalSx(i,j)**2+NormalSy(i,j)**2) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j)) &
! ! 		  ! Cross-diffusion from [Mathur and Murthy (1997)]
! 		  + CrossCoeff*( &
! 		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(diffusionE(i,j)+diffusionE(i+1,j))*( NeDual(i,j) - NeDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(diffusionE(i-1,j)+diffusionE(i,j))*( NeDual(i-1,j-1) - NeDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(diffusionE(i,j+1)+diffusionE(i,j))*( NeDual(i-1,j) - NeDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
! 		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(diffusionE(i,j-1)+diffusionE(i,j))*( NeDual(i,j-1) - NeDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! 		  ) &
! ! 		  ! Cross-diffusion from Mathur and Murthy 2007 without interpolation
! ! 		  + 0d0*(&
! ! 		    0.5d0*CellAreaE(i,j)*(diffusionE(i,j)+diffusionE(i+1,j)) * ( 0.5d0* (GradNeX(i,j)+GradNeX(i+1,j) )*NormalEx(i,j) + 0.5d0* (GradNeY(i,j)+GradNeY(i+1,j) )*NormalEy(i,j) - 0.5d0*( (GradNeX(i,j)+GradNeX(i+1,j) )*CurviEx(i,j) + ( GradNeY(i,j)+GradNeY(i+1,j) )*CurviEy(i,j) )/(NormalEx(i,j)*CurviEx(i,j)+NormalEy(i,j)*CurviEy(i,j))) &
! ! 		  + 0.5d0*CellAreaW(i,j)*(diffusionE(i,j)+diffusionE(i-1,j)) * ( 0.5d0* (GradNeX(i,j)+GradNeX(i-1,j) )*NormalWx(i,j) + 0.5d0* (GradNeY(i,j)+GradNeY(i-1,j) )*NormalWy(i,j) - 0.5d0*( (GradNeX(i,j)+GradNeX(i-1,j) )*CurviWx(i,j) + ( GradNeY(i,j)+GradNeY(i-1,j) )*CurviWy(i,j) )/(NormalWx(i,j)*CurviWx(i,j)+NormalWy(i,j)*CurviWy(i,j))) &
! ! 		  + 0.5d0*CellAreaN(i,j)*(diffusionE(i,j)+diffusionE(i,j+1)) * ( 0.5d0* (GradNeX(i,j)+GradNeX(i,j+1) )*NormalNx(i,j) + 0.5d0* (GradNeY(i,j)+GradNeY(i,j+1) )*NormalNy(i,j) - 0.5d0*( (GradNeX(i,j)+GradNeX(i,j+1) )*CurviNx(i,j) + ( GradNeY(i,j)+GradNeY(i,j+1) )*CurviNy(i,j) )/(NormalNx(i,j)*CurviNx(i,j)+NormalNy(i,j)*CurviNy(i,j))) &
! ! 		  + 0.5d0*CellAreaS(i,j)*(diffusionE(i,j)+diffusionE(i,j-1)) * ( 0.5d0* (GradNeX(i,j)+GradNeX(i,j-1) )*NormalSx(i,j) + 0.5d0* (GradNeY(i,j)+GradNeY(i,j-1) )*NormalSy(i,j) - 0.5d0*( (GradNeX(i,j)+GradNeX(i,j-1) )*CurviSx(i,j) + ( GradNeY(i,j)+GradNeY(i,j-1) )*CurviSy(i,j) )/(NormalSx(i,j)*CurviSx(i,j)+NormalSy(i,j)*CurviSy(i,j))) &
! ! 		  ) &
! 	      ))*dt/CellVol(i,j) & 
! 	      +Ne(i,j)
! 	      
! 	if(HolesOff.eq.0) then             
! 	             
! 	NhNew(i,j) = ( & 
! 	      (GainsH(i,j)-LossesH(i,j))*CellVol(i,j) & 
! 	      ! drift
! ! 	      -((0.5d0*(JhX(i+1,j)+JhX(i,j))*NormalEx(i,j)+0.5d0*(JhY(i+1,j) & 
! ! 	      +JhY(i,j))*NormalEy(i,j))*CellAreaE(i,j)+(0.5d0*(JhX(i,j)+JhX(i-1,j))*NormalWx(i,j)+0.5d0*(JhY(i,j) &
! ! 	      +JhY(i-1,j))*NormalWy(i,j))*CellAreaW(i,j)+(0.5d0*(JhX(i,j)+JhX(i,j+1))*NormalNx(i,j)+0.5d0*(JhY(i,j) & 
! ! 	      +JhY(i,j+1))*NormalNy(i,j))*CellAreaN(i,j)+(0.5d0*(JhX(i,j)+JhX(i,j-1))*NormalSx(i,j)+0.5d0*(JhY(i,j) & 
! ! 	      +JhY(i,j-1))*NormalSy(i,j))*CellAreaS(i,j)) &
! ! 	      ! diffusion operator on regular mesh
! ! 	      +(0.5d0*(Nh(i+1,j)-Nh(i,j))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 & 
! ! 	      +y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0)*(diffusionH(i+1,j)+diffusionH(i,j))*CellAreaE(i,j) & 
! ! 	      -0.5d0/(x(i,j)**2 & 
! ! 	      -2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j)+y(i-1,j)**2)**(0.5d0)*(diffusionH(i-1,j) & 
! ! 	      +diffusionH(i,j))*(Nh(i,j)-Nh(i-1,j))*CellAreaW(i,j) &
! ! 	      +0.5d0*(diffusionH(i,j+1)+diffusionH(i,j))*(Nh(i,j+1) & 
! ! 	      -Nh(i,j))*CellAreaN(i,j)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j)+x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j) & 
! ! 	      +y(i,j)**2)**(0.5d0) & 
! ! 	      -0.5d0*(diffusionH(i,j)+diffusionH(i,j-1))*(Nh(i,j)-Nh(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 &
! ! 	      -2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) & 
! 	      ! diffusion on irregular mesh
! 	      +(& 
! 	      + 0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(diffusionH(i,j)+diffusionH(i+1,j))*(Nh(i+1,j)-Nh(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
! ! 	      - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(diffusionH(i,j)+diffusionH(i+1,j))*(0.25d0*Nh(i+1,j+1)+0.25d0*Nh(i,j+1)-0.25d0*Nh(i+1,j-1)-0.25d0*Nh(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 	      - 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(diffusionH(i-1,j)+diffusionH(i,j))*(Nh(i,j)-Nh(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
! ! 	      - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(diffusionH(i-1,j)+diffusionH(i,j))*(0.25d0*Nh(i,j+1)+0.25d0*Nh(i-1,j+1)-0.25d0*Nh(i-1,j-1)-0.25d0*Nh(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 	      + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(diffusionH(i,j+1)+diffusionH(i,j))*(Nh(i,j+1)-Nh(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
! ! 	      - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(diffusionH(i,j+1)+diffusionH(i,j))*(0.25d0*Nh(i+1,j+1)+0.25d0*Nh(i+1,j)-0.25d0*Nh(i-1,j)-0.25d0*Nh(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
! 	      - 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(diffusionH(i,j-1)+diffusionH(i,j))*(Nh(i,j)-Nh(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
! ! 	      - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(diffusionH(i,j-1)+diffusionH(i,j))*(0.25d0*Nh(i+1,j)+0.25d0*Nh(i+1,j-1)-0.25d0*Nh(i-1,j-1)-0.25d0*Nh(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! ! 		  ! Cross-diffusion from [Mathur and Murthy (1997)]
! 	      + CrossCoeff*( &
! 	      + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(diffusionH(i,j)+diffusionH(i+1,j))*( NhDual(i,j) - NhDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 	      + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(diffusionH(i-1,j)+diffusionH(i,j))*( NhDual(i-1,j-1) - NhDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 	      + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(diffusionH(i,j+1)+diffusionH(i,j))*( NhDual(i-1,j) - NhDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
! 	      + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(diffusionH(i,j-1)+diffusionH(i,j))*( NhDual(i,j-1) - NhDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! 	      ) &
! 	      ))*dt/CellVol(i,j)+Nh(i,j)
! 	endif
! ! ! !  
 		    
! form with bug corrected in derivatives and (OmegaX, OmegaY) drift transport included in finite volumes
!     if(ConductivityFix < 2) then
	if(TeOff.ne.1) then 
	  if(ConvectionEnergy.eq.0) then
! 	  TeNew(i,j) = (&
! 		  + 0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappae(i,j)+kappae(i+1,j))*(Te(i+1,j)-Te(i,j)) / DistE(i,j) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j)) &
! 		  - 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappae(i-1,j)+kappae(i,j))*(Te(i,j)-Te(i-1,j)) / DistW(i,j) / (CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j)) & 
! 		  + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappae(i,j+1)+kappae(i,j))*(Te(i,j+1)-Te(i,j)) / DistN(i,j) / (CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j)) & 
! 		  - 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappae(i,j-1)+kappae(i,j))*(Te(i,j)-Te(i,j-1)) / DistS(i,j) / (CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j)) & 
! ! 		  + CrossCoeff*( &
! ! 		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappae(i,j)+kappae(i+1,j))*( TeDual(i,j) - TeDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! ! 		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappae(i-1,j)+kappae(i,j))*( TeDual(i-1,j-1) - TeDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! ! 		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappae(i,j+1)+kappae(i,j))*( TeDual(i-1,j) - TeDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
! ! 		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappae(i,j-1)+kappae(i,j))*( TeDual(i,j-1) - TeDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! ! 		  ) &
! 		  +(-CouplingE(i,j)+SourceE(i,j))*CellVol(i,j)) & 
! 		  /Ce(i,j) * dt/CellVol(i,j) &
! 		  +Te(i,j)
		      
	    else !convective term included!
	    
		UeNew(i,j) = ((SourceUe(i,j)-CouplingE(i,j))*CellVol(i,j) & 
		! convective term for transport of the energy by the field
		-((0.5d0*(VeX(i+1,j)+VeX(i,j))*NormalEx(i,j)+0.5d0*(VeY(i+1,j)+VeY(i,j))*NormalEy(i,j)) * CellAreaE(i,j) & 
		+( 0.5d0*(VeX(i,j)+VeX(i-1,j))*NormalWx(i,j)+0.5d0*(VeY(i,j)+VeY(i-1,j))*NormalWy(i,j)) * CellAreaW(i,j) & 
		+( 0.5d0*(VeX(i,j)+VeX(i,j+1))*NormalNx(i,j)+0.5d0*(VeY(i,j)+VeY(i,j+1))*NormalNy(i,j)) * CellAreaN(i,j) &
		+( 0.5d0*(VeX(i,j)+VeX(i,j-1))*NormalSx(i,j)+0.5d0*(VeY(i,j)+VeY(i,j-1))*NormalSy(i,j)) * CellAreaS(i,j)) &
		! diffusive term for energy - rewrite correctly
! 		+ ( (Ue(i+1,j)-Ue(i,j))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 & 
! 		+ y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0)*( (kappae(i+1,j)+kappae(i,j))/(Ce(i+1,j)+Ce(i,j)) )*CellAreaE(i,j) &
! 		- 1d0/(x(i,j)**2-2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j) + & 
! 		y(i-1,j)**2)**(0.5d0)*((kappae(i-1,j)+kappae(i,j))/(Ce(i-1,j)+Ce(i,j)))*(Ue(i,j)-Ue(i-1,j))*CellAreaW(i,j) & 
! 		+ 1d0*((kappae(i,j+1)+kappae(i,j))/(Ce(i,j+1)+Ce(i,j)))*(Ue(i,j+1)-Ue(i,j))*CellAreaN(i,j)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j) & 
! 		+x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j)+y(i,j)**2)**(0.5d0) & 
! 		-1d0*((kappae(i,j)+kappae(i,j-1))/(Ce(i,j)+Ce(i,j-1)))*(Ue(i,j)-Ue(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 & 
! 		-2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) &
		)*dt/CellVol(i,j) & 
		+( &
		    0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*5d0/3d0*(kappae(i,j)/Ce(i,j)+kappae(i+1,j)/Ce(i+1,j))*(Ue(i+1,j)-Ue(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
		  - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*5d0/3d0*(kappae(i,j)/Ce(i,j)+kappae(i+1,j)/Ce(i+1,j))*(0.25d0*Ue(i+1,j+1)+0.25d0*Ue(i,j+1)-0.25d0*Ue(i+1,j-1)-0.25d0*Ue(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
		  + 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*5d0/3d0*(kappae(i-1,j)/Ce(i-1,j)+kappae(i,j)/Ce(i,j))*(Ue(i,j)-Ue(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
		  - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*5d0/3d0*(kappae(i-1,j)/Ce(i-1,j)+kappae(i,j)/Ce(i,j))*(0.25d0*Ue(i,j+1)+0.25d0*Ue(i-1,j+1)-0.25d0*Ue(i-1,j-1)-0.25d0*Ue(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
		  + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*5d0/3d0*(kappae(i,j+1)/Ce(i,j+1)+kappae(i,j)/Ce(i,j))*(Ue(i,j+1)-Ue(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
		  - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*5d0/3d0*(kappae(i,j+1)/Ce(i,j+1)+kappae(i,j)/Ce(i,j))*(0.25d0*Ue(i+1,j+1)+0.25d0*Ue(i+1,j)-0.25d0*Ue(i-1,j)-0.25d0*Ue(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
		  + 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*5d0/3d0*(kappae(i,j-1)/Ce(i,j-1)+kappae(i,j)/Ce(i,j))*(Ue(i,j)-Ue(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
		  - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*5d0/3d0*(kappae(i,j-1)/Ce(i,j-1)+kappae(i,j)/Ce(i,j))*(0.25d0*Ue(i+1,j)+0.25d0*Ue(i+1,j-1)-0.25d0*Ue(i-1,j-1)-0.25d0*Ue(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! 		    Korfiatis 2007 equation
! 		    0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappae(i,j)*diffusionE(i,j)/Ne(i,j)+kappae(i+1,j)*diffusionE(i+1,j)/Ne(i+1,j))*(Ne(i+1,j)-Ne(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
! 		  - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappae(i,j)*diffusionE(i,j)/Ne(i,j)+kappae(i+1,j)*diffusionE(i+1,j)/Ne(i+1,j))*(0.25d0*Ne(i+1,j+1)+0.25d0*Ne(i,j+1)-0.25d0*Ne(i+1,j-1)-0.25d0*Ne(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 		  + 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappae(i-1,j)*diffusionE(i-1,j)/Ne(i-1,j)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j)-Ne(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
! 		  - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappae(i-1,j)*diffusionE(i-1,j)/Ne(i-1,j)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i,j+1)+0.25d0*Ne(i-1,j+1)-0.25d0*Ne(i-1,j-1)-0.25d0*Ne(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 		  + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappae(i,j+1)*diffusionE(i,j+1)/Ne(i,j+1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j+1)-Ne(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
! 		  - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappae(i,j+1)*diffusionE(i,j+1)/Ne(i,j+1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i+1,j+1)+0.25d0*Ne(i+1,j)-0.25d0*Ne(i-1,j)-0.25d0*Ne(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
! 		  + 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappae(i,j-1)*diffusionE(i,j-1)/Ne(i,j-1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j)-Ne(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
! 		  - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappae(i,j-1)*diffusionE(i,j-1)/Ne(i,j-1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i+1,j)+0.25d0*Ne(i+1,j-1)-0.25d0*Ne(i-1,j-1)-0.25d0*Ne(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
		  ) & 
		  *dt/CellVol(i,j) &
		+Ue(i,j)
		
		UhNew(i,j) = ((SourceUh(i,j)-CouplingH(i,j))*CellVol(i,j) & 
		! convective term for transport of the energy by the field
		-((0.5d0*(VhX(i+1,j)+VhX(i,j))*NormalEx(i,j)+0.5d0*(VhY(i+1,j)+VhY(i,j))*NormalEy(i,j)) * CellAreaE(i,j) & 
		+( 0.5d0*(VhX(i,j)+VhX(i-1,j))*NormalWx(i,j)+0.5d0*(VhY(i,j)+VhY(i-1,j))*NormalWy(i,j)) * CellAreaW(i,j) & 
		+( 0.5d0*(VhX(i,j)+VhX(i,j+1))*NormalNx(i,j)+0.5d0*(VhY(i,j)+VhY(i,j+1))*NormalNy(i,j)) * CellAreaN(i,j) &
		+( 0.5d0*(VhX(i,j)+VhX(i,j-1))*NormalSx(i,j)+0.5d0*(VhY(i,j)+VhY(i,j-1))*NormalSy(i,j)) * CellAreaS(i,j)) &
		! diffusive term for energy - rewrite correctly
! 		+ ( (Uh(i+1,j)-Uh(i,j))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 & 
! 		+ y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0)*( (kappah(i+1,j)+kappah(i,j))/(Ch(i+1,j)+Ch(i,j)) )*CellAreaE(i,j) &
! 		- 1d0/(x(i,j)**2-2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j) + & 
! 		y(i-1,j)**2)**(0.5d0)*((kappah(i-1,j)+kappah(i,j))/(Ch(i-1,j)+Ch(i,j)))*(Uh(i,j)-Uh(i-1,j))*CellAreaW(i,j) & 
! 		+ 1d0*((kappah(i,j+1)+kappah(i,j))/(Ch(i,j+1)+Ch(i,j)))*(Uh(i,j+1)-Uh(i,j))*CellAreaN(i,j)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j) & 
! 		+x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j)+y(i,j)**2)**(0.5d0) & 
! 		-1d0*((kappah(i,j)+kappah(i,j-1))/(Ch(i,j)+Ch(i,j-1)))*(Uh(i,j)-Uh(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 & 
! 		-2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) &
		)*dt/CellVol(i,j) & 
		+( & 
		    0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappah(i,j)/Ch(i,j)+kappah(i+1,j)/Ch(i+1,j))*(Uh(i+1,j)-Uh(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
		  - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappah(i,j)/Ch(i,j)+kappah(i+1,j)/Ch(i+1,j))*(0.25d0*Uh(i+1,j+1)+0.25d0*Uh(i,j+1)-0.25d0*Uh(i+1,j-1)-0.25d0*Uh(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
		  + 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappah(i-1,j)/Ch(i-1,j)+kappah(i,j)/Ch(i,j))*(Uh(i,j)-Uh(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
		  - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappah(i-1,j)/Ch(i-1,j)+kappah(i,j)/Ch(i,j))*(0.25d0*Uh(i,j+1)+0.25d0*Uh(i-1,j+1)-0.25d0*Uh(i-1,j-1)-0.25d0*Uh(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
		  + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappah(i,j+1)/Ch(i,j+1)+kappah(i,j)/Ch(i,j))*(Uh(i,j+1)-Uh(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
		  - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappah(i,j+1)/Ch(i,j+1)+kappah(i,j)/Ch(i,j))*(0.25d0*Uh(i+1,j+1)+0.25d0*Uh(i+1,j)-0.25d0*Uh(i-1,j)-0.25d0*Uh(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
		  + 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappah(i,j-1)/Ch(i,j-1)+kappah(i,j)/Ch(i,j))*(Uh(i,j)-Uh(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
		  - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappah(i,j-1)/Ch(i,j-1)+kappah(i,j)/Ch(i,j))*(0.25d0*Uh(i+1,j)+0.25d0*Uh(i+1,j-1)-0.25d0*Uh(i-1,j-1)-0.25d0*Uh(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
! ! Korfiatis 2007 equation
! 		    0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappah(i,j)*diffusionH(i,j)/Nh(i,j)+kappah(i+1,j)*diffusionH(i+1,j)/Nh(i+1,j))*(Nh(i+1,j)-Nh(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
! 		  - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappah(i,j)*diffusionH(i,j)/Nh(i,j)+kappah(i+1,j)*diffusionH(i+1,j)/Nh(i+1,j))*(0.25d0*Nh(i+1,j+1)+0.25d0*Nh(i,j+1)-0.25d0*Nh(i+1,j-1)-0.25d0*Nh(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 		  + 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappah(i-1,j)*diffusionH(i-1,j)/Nh(i-1,j)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(Nh(i,j)-Nh(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
! 		  - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappah(i-1,j)*diffusionH(i-1,j)/Nh(i-1,j)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(0.25d0*Nh(i,j+1)+0.25d0*Nh(i-1,j+1)-0.25d0*Nh(i-1,j-1)-0.25d0*Nh(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 		  + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappah(i,j+1)*diffusionH(i,j+1)/Nh(i,j+1)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(Nh(i,j+1)-Nh(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
! 		  - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappah(i,j+1)*diffusionH(i,j+1)/Nh(i,j+1)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(0.25d0*Nh(i+1,j+1)+0.25d0*Nh(i+1,j)-0.25d0*Nh(i-1,j)-0.25d0*Nh(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
! 		  + 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappah(i,j-1)*diffusionH(i,j-1)/Nh(i,j-1)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(Nh(i,j)-Nh(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
! 		  - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappah(i,j-1)*diffusionH(i,j-1)/Nh(i,j-1)+kappah(i,j)*diffusionH(i,j)/Nh(i,j))*(0.25d0*Nh(i+1,j)+0.25d0*Nh(i+1,j-1)-0.25d0*Nh(i-1,j-1)-0.25d0*Nh(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
		  ) & 
		  *dt/CellVol(i,j) &
		+Uh(i,j)
		
	    end if
  !     else 
  !       TeNew(i,j) = ((-CouplingE(i,j)+SourceE(i,j))*((1d0/8d0)*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1))-(1d0/8d0)*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1))) &
  ! 		      +.2500000000*kb**2d0*((0.5d0)*Te(i+1,j)+(0.5d0)*Te(i,j))**2d0*(mobilityE(i+1,j)+mobilityE(i,j))*(FermiTableE(ColFermi1,FermiIndexE(i+1,j))+FermiTableE(ColFermi1,FermiIndexE(i,j)))*(FermiTableE(ColFermiHalf,FermiIndexE(i+1,j)) &
  ! 		      +FermiTableE(ColFermiHalf,FermiIndexE(i,j)))*(Ne(i+1,j)-Ne(i,j))*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)-2d0*x(i+1,j+1)*x(i+1,j-1)-2d0*x(i+1,j+1)*x(i,j-1)+x(i,j+1)**2 &
  ! 		      -2d0*x(i,j+1)*x(i+1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)-2d0*y(i+1,j+1)*y(i+1,j-1) &
  ! 		      -2d0*y(i+1,j+1)*y(i,j-1)+y(i,j+1)**2-2d0*y(i,j+1)*y(i+1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1) &
  ! 		      +y(i,j-1)**2)**(0.5)/ec/(FermiTableE(ColFermi0,FermiIndexE(i+1,j))+FermiTableE(ColFermi0,FermiIndexE(i,j)))/(FermiTableE(ColFermiMenusHalf,FermiIndexE(i+1,j))+FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 &
  ! 		      +y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5)-.2500000000*kb**2d0*((0.5d0)*Te(i-1,j)+(0.5d0)*Te(i,j))**2d0*(mobilityE(i-1,j)+mobilityE(i,j))*(FermiTableE(ColFermi1,FermiIndexE(i-1,j)) &
  ! 		      +FermiTableE(ColFermi1,FermiIndexE(i,j)))*(FermiTableE(ColFermiHalf,FermiIndexE(i-1,j))+FermiTableE(ColFermiHalf,FermiIndexE(i,j)))*(Ne(i,j)-Ne(i-1,j))*(x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)-2d0*x(i,j+1)*x(i-1,j-1) &
  ! 		      -2d0*x(i,j+1)*x(i,j-1)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i-1,j-1)-2d0*x(i-1,j+1)*x(i,j-1)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i,j+1)**2 &
  ! 		      +2d0*y(i,j+1)*y(i-1,j+1)-2d0*y(i,j+1)*y(i-1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i-1,j-1)-2d0*y(i-1,j+1)*y(i,j-1)+y(i-1,j-1)**2 &
  ! 		      +2d0*y(i-1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5)/ec/(FermiTableE(ColFermi0,FermiIndexE(i-1,j))+FermiTableE(ColFermi0,FermiIndexE(i,j)))/(FermiTableE(ColFermiMenusHalf,FermiIndexE(i-1,j))+FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))/(x(i,j)**2 &
  ! 		      -2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j)+y(i-1,j)**2)**(0.5)+.2500000000*kb**2d0*((0.5d0)*Te(i,j+1)+(0.5d0)*Te(i,j))**2d0*(mobilityE(i,j+1) &
  ! 		      +mobilityE(i,j))*(FermiTableE(ColFermi1,FermiIndexE(i,j+1))+FermiTableE(ColFermi1,FermiIndexE(i,j)))*(FermiTableE(ColFermiHalf,FermiIndexE(i,j+1))+FermiTableE(ColFermiHalf,FermiIndexE(i,j)))*(Ne(i,j+1)-Ne(i,j))*(x(i+1,j+1)**2 &
  ! 		      +2d0*x(i+1,j+1)*x(i+1,j)-2d0*x(i+1,j+1)*x(i-1,j)-2d0*x(i+1,j+1)*x(i-1,j+1)+x(i+1,j)**2-2d0*x(i+1,j)*x(i-1,j)-2d0*x(i+1,j)*x(i-1,j+1)+x(i-1,j)**2 &
  ! 		      +2d0*x(i-1,j)*x(i-1,j+1)+x(i-1,j+1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)-2d0*y(i+1,j+1)*y(i-1,j)-2d0*y(i+1,j+1)*y(i-1,j+1)+y(i+1,j)**2 &
  ! 		      -2d0*y(i+1,j)*y(i-1,j)-2d0*y(i+1,j)*y(i-1,j+1)+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5)/ec/(FermiTableE(ColFermi0,FermiIndexE(i,j+1)) & 
  ! 		      +FermiTableE(ColFermi0,FermiIndexE(i,j)))/(FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j+1))+FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j)+x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j) &
  ! 		      +y(i,j)**2)**(0.5)-.2500000000*kb**2d0*((0.5d0)*Te(i,j-1)+(0.5d0)*Te(i,j))**2d0*(mobilityE(i,j-1)+mobilityE(i,j))*(FermiTableE(ColFermi1,FermiIndexE(i,j-1)) &
  ! 		      +FermiTableE(ColFermi1,FermiIndexE(i,j)))*(FermiTableE(ColFermiHalf,FermiIndexE(i,j-1))+FermiTableE(ColFermiHalf,FermiIndexE(i,j)))*(Ne(i,j)-Ne(i,j-1))*(x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)-2d0*x(i+1,j)*x(i-1,j-1) &
  ! 		      -2d0*x(i+1,j)*x(i-1,j)+x(i+1,j-1)**2-2d0*x(i+1,j-1)*x(i-1,j-1)-2d0*x(i+1,j-1)*x(i-1,j)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)+x(i-1,j)**2+y(i+1,j)**2 & 
  ! 		      +2d0*y(i+1,j)*y(i+1,j-1)-2d0*y(i+1,j)*y(i-1,j-1)-2d0*y(i+1,j)*y(i-1,j)+y(i+1,j-1)**2-2d0*y(i+1,j-1)*y(i-1,j-1)-2d0*y(i+1,j-1)*y(i-1,j)+y(i-1,j-1)**2 & 
  ! 		      +2d0*y(i-1,j-1)*y(i-1,j)+y(i-1,j)**2)**(0.5)/ec/(FermiTableE(ColFermi0,FermiIndexE(i,j-1))+FermiTableE(ColFermi0,FermiIndexE(i,j)))/(FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j-1))+FermiTableE(ColFermiMenusHalf,FermiIndexE(i,j)))/(x(i,j)**2 & 
  ! 		      -2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5)+(0.25d0)/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2+y(i+1,j)**2-2d0*y(i+1,j)*y(i,j) & 
  ! 		      +y(i,j)**2)**(0.5)*kappae(i+(0.5d0),j)*(Te(i+1,j)-Te(i,j))*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)-2d0*x(i+1,j+1)*x(i+1,j-1)-2d0*x(i+1,j+1)*x(i,j-1) & 
  ! 		      +x(i,j+1)**2-2d0*x(i,j+1)*x(i+1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i+1,j-1)**2+2d0*x(i+1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1) & 
  ! 		      -2d0*y(i+1,j+1)*y(i+1,j-1)-2d0*y(i+1,j+1)*y(i,j-1)+y(i,j+1)**2-2d0*y(i,j+1)*y(i+1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1) & 
  ! 		      +y(i,j-1)**2)**(0.5)-(0.25d0)/(x(i,j)**2-2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j)+y(i-1,j)**2)**(0.5)*kappae(i-(0.5d0),j)*(Te(i,j) & 
  ! 		      -Te(i-1,j))*(x(i,j+1)**2+2d0*x(i,j+1)*x(i-1,j+1)-2d0*x(i,j+1)*x(i-1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i-1,j+1)**2-2d0*x(i-1,j+1)*x(i-1,j-1) & 
  ! 		      -2d0*x(i-1,j+1)*x(i,j-1)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)-2d0*y(i,j+1)*y(i-1,j-1) & 
  ! 		      -2d0*y(i,j+1)*y(i,j-1)+y(i-1,j+1)**2-2d0*y(i-1,j+1)*y(i-1,j-1)-2d0*y(i-1,j+1)*y(i,j-1)+y(i-1,j-1)**2+2d0*y(i-1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5) & 
  ! 		      +(0.25d0)*kappae(i,j+(0.5d0))*(Te(i,j+1)-Te(i,j))*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)-2d0*x(i+1,j+1)*x(i-1,j)-2d0*x(i+1,j+1)*x(i-1,j+1)+x(i+1,j)**2 &
  ! 		      -2d0*x(i+1,j)*x(i-1,j)-2d0*x(i+1,j)*x(i-1,j+1)+x(i-1,j)**2+2d0*x(i-1,j)*x(i-1,j+1)+x(i-1,j+1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j) & 
  ! 		      -2d0*y(i+1,j+1)*y(i-1,j)-2d0*y(i+1,j+1)*y(i-1,j+1)+y(i+1,j)**2-2d0*y(i+1,j)*y(i-1,j)-2d0*y(i+1,j)*y(i-1,j+1)+y(i-1,j)**2+2d0*y(i-1,j)*y(i-1,j+1) &
  ! 		      +y(i-1,j+1)**2)**(0.5)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j)+x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j)+y(i,j)**2)**(0.5) & 
  ! 		      -(0.25d0)*kappae(i,j-(0.5d0))*(Te(i,j)-Te(i,j-1))*(x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)-2d0*x(i+1,j)*x(i-1,j-1)-2d0*x(i+1,j)*x(i-1,j)+x(i+1,j-1)**2 & 
  ! 		      -2d0*x(i+1,j-1)*x(i-1,j-1)-2d0*x(i+1,j-1)*x(i-1,j)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i-1,j)+x(i-1,j)**2+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1) & 
  ! 		      -2d0*y(i+1,j)*y(i-1,j-1)-2d0*y(i+1,j)*y(i-1,j)+y(i+1,j-1)**2-2d0*y(i+1,j-1)*y(i-1,j-1)-2d0*y(i+1,j-1)*y(i-1,j)+y(i-1,j-1)**2 & 
  ! 		      +2d0*y(i-1,j-1)*y(i-1,j)+y(i-1,j)**2)**(0.5)/(x(i,j)**2-2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1) & 
  ! 		      +y(i,j-1)**2)**(0.5))/Ce(i,j)*dt/((1d0/8d0)*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1))-(1d0/8d0)*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1))) & 
  ! 		      +Te(i,j)
  !     end if
      

  !  scheme ready for enhanced conductivity and drift of vector (OmegaX, OmegaY). 
  ! 	if(ConductivityFix < 2) then
	if(HolesOff.eq.0) then
	    if(ConvectionEnergy.eq.0) then 
	    
! 	    ThNew(i,j) = (&
! 		  + 0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappah(i,j)+kappah(i+1,j))*(Th(i+1,j)-Th(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
! ! 		  - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappah(i,j)+kappah(i+1,j))*(0.25d0*Th(i+1,j+1)+0.25d0*Th(i,j+1)-0.25d0*Th(i+1,j-1)-0.25d0*Th(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 		  - 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappah(i-1,j)+kappah(i,j))*(Th(i,j)-Th(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
! ! 		  - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappah(i-1,j)+kappah(i,j))*(0.25d0*Th(i,j+1)+0.25d0*Th(i-1,j+1)-0.25d0*Th(i-1,j-1)-0.25d0*Th(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 		  + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappah(i,j+1)+kappah(i,j))*(Th(i,j+1)-Th(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
! ! 		  - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappah(i,j+1)+kappah(i,j))*(0.25d0*Th(i+1,j+1)+0.25d0*Th(i+1,j)-0.25d0*Th(i-1,j)-0.25d0*Th(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
! 		  - 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappah(i,j-1)+kappah(i,j))*(Th(i,j)-Th(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
! ! 		  - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappah(i,j-1)+kappah(i,j))*(0.25d0*Th(i+1,j)+0.25d0*Th(i+1,j-1)-0.25d0*Th(i-1,j-1)-0.25d0*Th(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! 		  + CrossCoeff*( &
! 		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappah(i,j)+kappah(i+1,j))*( ThDual(i,j) - ThDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappah(i-1,j)+kappah(i,j))*( ThDual(i-1,j-1) - ThDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappah(i,j+1)+kappah(i,j))*( ThDual(i-1,j) - ThDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
! 		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappah(i,j-1)+kappah(i,j))*( ThDual(i,j-1) - ThDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! 		  ) &
! 		  +(-CouplingH(i,j)+SourceH(i,j))*CellVol(i,j)) & 
! 		  /Ch(i,j)*dt/CellVol(i,j)+Th(i,j)
	    
!
		      
	else !convection scheme
	    UhNew(i,j) = ((SourceUh(i,j)-CouplingH(i,j))*CellVol(i,j)-( & 
		   (0.5d0*(VhX(i+1,j)+VhX(i,j))*NormalEx(i,j)+0.5d0*(VhY(i+1,j)+VhY(i,j))*NormalEy(i,j)) * CellAreaE(i,j) & 
		  +(0.5d0*(VhX(i,j)+VhX(i-1,j))*NormalWx(i,j)+0.5d0*(VhY(i,j)+VhY(i-1,j))*NormalWy(i,j)) * CellAreaW(i,j) & 
		  +(0.5d0*(VhX(i,j)+VhX(i,j+1))*NormalNx(i,j)+0.5d0*(VhY(i,j)+VhY(i,j+1))*NormalNy(i,j)) * CellAreaN(i,j) &
		  +(0.5d0*(VhX(i,j)+VhX(i,j-1))*NormalSx(i,j)+0.5d0*(VhY(i,j)+VhY(i,j-1))*NormalSy(i,j)) * CellAreaS(i,j)) &
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
		  )*dt/CellVol(i,j)+Uh(i,j)
	end if
	

    end if		      
! 	  TsNew(i,j) = (0.25d0/(x(i+1,j)-x(i,j))*(0.5d0*kappas(i,j)+0.5d0*kappas(i+1,j))*(Ts(i+1,j)- & 
! 		      Ts(i,j))*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i,j+1)-2d0*x(i+1,j+1)*x(i+1,j-1)- & 
! 		      2d0*x(i+1,j+1)*x(i,j-1)+x(i,j+1)**2-2d0*x(i,j+1)*x(i+1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i+1,j-1)**2 & 
! 		      +2d0*x(i+1,j-1)*x(i,j-1)+x(i,j-1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i,j+1)-2d0*y(i+1,j+1)*y(i+1,j-1) & 
! 		      -2d0*y(i+1,j+1)*y(i,j-1)+y(i,j+1)**2-2d0*y(i,j+1)*y(i+1,j-1)-2d0*y(i,j+1)*y(i,j-1) & 
! 		      +y(i+1,j-1)**2+2d0*y(i+1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5)- & 
! 		      0.25d0/(x(i,j)-x(i-1,j))*(0.5d0*kappas(i-1,j)+0.5d0*kappas(i,j))*(Ts(i,j)-Ts(i-1,j))*(x(i,j+1)**2 & 
! 		      +2d0*x(i,j+1)*x(i-1,j+1)-2d0*x(i,j+1)*x(i-1,j-1)-2d0*x(i,j+1)*x(i,j-1)+x(i-1,j+1)**2 & 
! 		      -2d0*x(i-1,j+1)*x(i-1,j-1)-2d0*x(i-1,j+1)*x(i,j-1)+x(i-1,j-1)**2+2d0*x(i-1,j-1)*x(i,j-1)+x(i,j-1)**2 & 
! 		      +y(i,j+1)**2+2d0*y(i,j+1)*y(i-1,j+1)-2d0*y(i,j+1)*y(i-1,j-1)-2d0*y(i,j+1)*y(i,j-1)+y(i-1,j+1)**2 & 
! 		      -2d0*y(i-1,j+1)*y(i-1,j-1)-2d0*y(i-1,j+1)*y(i,j-1)+y(i-1,j-1)**2 & 
! 		      +2d0*y(i-1,j-1)*y(i,j-1)+y(i,j-1)**2)**(0.5)+0.25d0/(y(i,j+1)-y(i,j))*(0.5d0*kappas(i,j+1) & 
! 		      +0.5d0*kappas(i,j))*(Ts(i,j+1)-Ts(i,j))*(x(i+1,j+1)**2+2d0*x(i+1,j+1)*x(i+1,j)-2d0*x(i+1,j+1)*x(i-1,j) & 
! 		      -2d0*x(i+1,j+1)*x(i-1,j+1)+x(i+1,j)**2-2d0*x(i+1,j)*x(i-1,j)-2d0*x(i+1,j)*x(i-1,j+1)+x(i-1,j)**2 & 
! 		      +2d0*x(i-1,j)*x(i-1,j+1)+x(i-1,j+1)**2+y(i+1,j+1)**2+2d0*y(i+1,j+1)*y(i+1,j)-2d0*y(i+1,j+1)*y(i-1,j) & 
! 		      -2d0*y(i+1,j+1)*y(i-1,j+1)+y(i+1,j)**2-2d0*y(i+1,j)*y(i-1,j)-2d0*y(i+1,j)*y(i-1,j+1)+y(i-1,j)**2 & 
! 		      +2d0*y(i-1,j)*y(i-1,j+1)+y(i-1,j+1)**2)**(0.5)-0.25d0/(y(i,j)-y(i,j-1))*(0.5d0*kappas(i,j-1) & 
! 		      +0.5d0*kappas(i,j))*(Ts(i,j)-Ts(i,j-1))*(x(i+1,j)**2+2d0*x(i+1,j)*x(i+1,j-1)-2d0*x(i+1,j)*x(i-1,j-1) & 
! 		      -2d0*x(i+1,j)*x(i-1,j)+x(i+1,j-1)**2-2d0*x(i+1,j-1)*x(i-1,j-1)-2d0*x(i+1,j-1)*x(i-1,j)+x(i-1,j-1)**2 & 
! 		      +2d0*x(i-1,j-1)*x(i-1,j)+x(i-1,j)**2+y(i+1,j)**2+2d0*y(i+1,j)*y(i+1,j-1)-2d0*y(i+1,j)*y(i-1,j-1) & 
! 		      -2d0*y(i+1,j)*y(i-1,j)+y(i+1,j-1)**2-2d0*y(i+1,j-1)*y(i-1,j-1)-2d0*y(i+1,j-1)*y(i-1,j)+y(i-1,j-1)**2 & 
! 		      +2d0*y(i-1,j-1)*y(i-1,j)+y(i-1,j)**2)**(0.5) & 
! 		      +(CouplingH(i,j)+CouplingE(i,j))*(1d0/8d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1) & 
! 		      -y(i+1,j-1))-1d0/8d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1))))/Cs(i,j)*dt/(1d0/8d0*(x(i+1,j+1) & 
! 		      -x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1))-1d0/8d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1))) & 
! 		      +Ts(i,j)

    ! version with cross-diffusion
!     TsNew(i,j) =  (&
! 		  + 0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*(Ts(i+1,j)-Ts(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
! ! 		  - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*(0.25d0*Ts(i+1,j+1)+0.25d0*Ts(i,j+1)-0.25d0*Ts(i+1,j-1)-0.25d0*Ts(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 		  - 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*(Ts(i,j)-Ts(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) & 
! ! 		  - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*(0.25d0*Ts(i,j+1)+0.25d0*Ts(i-1,j+1)-0.25d0*Ts(i-1,j-1)-0.25d0*Ts(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 		  + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*(Ts(i,j+1)-Ts(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) & 
! ! 		  - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*(0.25d0*Ts(i+1,j+1)+0.25d0*Ts(i+1,j)-0.25d0*Ts(i-1,j)-0.25d0*Ts(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) & 
! 		  - 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*(Ts(i,j)-Ts(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) & 
! ! 		  - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*(0.25d0*Ts(i+1,j)+0.25d0*Ts(i+1,j-1)-0.25d0*Ts(i-1,j-1)-0.25d0*Ts(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! 		  + CrossCoeff*( &
! 		  + 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappas(i,j)+kappas(i+1,j))*( TsDual(i,j) - TsDual(i,j-1) )/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) & 
! 		  + 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappas(i-1,j)+kappas(i,j))*( TsDual(i-1,j-1) - TsDual(i-1,j) )/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
! 		  + 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappas(i,j+1)+kappas(i,j))*( TsDual(i-1,j) - TsDual(i,j) )/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
! 		  + 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappas(i,j-1)+kappas(i,j))*( TsDual(i,j-1) - TsDual(i-1,j-1) )/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) & 
! 		  ) &
! 		  +(CouplingE(i,j)+CouplingH(i,j))*CellVol(i,j)) & 
! 		  /Cs(i,j)*dt/CellVol(i,j)+Ts(i,j)

    if(ConvectionEnergy.eq.1) then !define temperatures from energy
      TeNew(i,j) = Te(i,j) + ((UeNew(i,j) -  Ue(i,j))-1.5d0*kb*Te(i,j)*(NeNew(i,j) - Ne(i,j))*FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)) ) / Ce(i,j)
      ThNew(i,j) = Th(i,j) + ((UhNew(i,j) -  Uh(i,j))-1.5d0*kb*Th(i,j)*(NhNew(i,j) - Nh(i,j))*FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)) ) / Ch(i,j)
    end if
		      
    ! 	TsNew(i,j) = (.5*(kappas(i,j)+kappas(i+1,j))*(Ts(i+1,j)-Ts(i,j))*dy/dx &
    ! 		      -.5*(kappas(i-1,j)+kappas(i,j))*(Ts(i,j)-Ts(i-1,j))*dy/dx & 
    ! 		      + 0.5*(kappas(i,j+1)+kappas(i,j))*(Ts(i,j+1)-Ts(i,j))*dx/dy &
    ! 		      - .5*(kappas(i,j-1)+kappas(i,j))*(Ts(i,j)-Ts(i,j-1))*dx/dy &
    ! 		      + CouplingE(i,j)*(Te(i,j)-Ts(i,j))*dx*dy &
    ! 		      + CouplingH(i,j)*(Th(i,j)-Ts(i,j))*dx*dy &
    ! 		      )/Cs(i,j)/dx/dy*dt & 
    ! 		      + Ts(i,j)

	  CFLxT(i,j)=kappae(i,j)/Ce(i,j) * dt/(x(i,j)-x(i-1,j))**2
	  CFLyT(i,j)=kappae(i,j)/Ce(i,j) * dt/(y(i,j)-y(i,j-1))**2
	  CFLxTs(i,j)=kappas(i,j)/Cs(i,j) * dt/(x(i,j)-x(i-1,j))**2
	  CFLyTs(i,j)=kappas(i,j)/Cs(i,j) * dt/(y(i,j)-y(i,j-1))**2
	  
	end if
	
	CFLxN(i,j)=diffusionE(i,j)*dt/(x(i,j)-x(i-1,j))**2 !+dt/(x(i,j)-x(i-1,j))*mobilityE(i,j)*sqrt(Ex(i,j)**2+Ey(i,j)**2)
	CFLyN(i,j)=diffusionE(i,j)*dt/(y(i,j)-y(i,j-1))**2 !+dt/(y(i,j)-y(i,j-1))*mobilityE(i,j)*sqrt(Ex(i,j)**2+Ey(i,j)**2)
	
	TotalElectrons(i,j)=NeNew(i,j)*(1d0/8d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1))-1d0/8d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1)))
	TotalHoles(i,j)=NhNew(i,j)*(1d0/8d0*(x(i+1,j+1)-x(i-1,j-1))*(y(i-1,j+1)-y(i+1,j-1))-1d0/8d0*(x(i-1,j+1)-x(i+1,j-1))*(y(i+1,j+1)-y(i-1,j-1)))
	
	
	ThermalEnergy(i,j)=1.5d0*kb*Te(i,j)*Ne(i,j)+1.5d0*kb*Th(i,j)*Nh(i,j)+1.5d0*kb*Ts(i,j)*SiDensity
	LaserEnergy(i,j)=OnePhotonIonizationRate(lambda,epsilonInf)*intensity(i,j)/(1d0-reflectivity(i,j))+ & !energy loss by interband absorption
			  TwoPhotonIonizationRate(lambda)*intensity(i,j)**2/(1d0-reflectivity(i,j))**2+ & !energy loss by two photon absorption
			  (absorptionDrudeE(i,j)+absorptionDrudeH(i,j))*intensity(i,j)/(1d0-reflectivity(i,j)) !energy loss by carrrier heating
	

      end do
      !$OMP END PARALLEL DO
      
    end do
    !$OMP END DO

  ! $ OMP END PARALLEL !!end of parallel section
    
    !boundary conditions
    
    do i=1, M !North and South boundaries
      ! finite differences finite difference fashion
      if(DriftOn.eq.0) then
!  	NeNew(i,1)=NeNew(i,2)	! disable if implicit solvers are used
! 	NhNew(i,1)=NhNew(i,2)
! 	NeNew(i,N)=NeNew(i,N-1)
! 	NhNew(i,N)=NhNew(i,N-1)
      end if

	UeNew(i,1)=UeNew(i,2)
	UhNew(i,1)=UhNew(i,2)
! 	TeNew(i,1)=TeNew(i,2)
! 	ThNew(i,1)=ThNew(i,2)
! 	TsNew(i,1)=TsNew(i,2)

	UeNew(i,N)=UeNew(i,N-1)
	UhNew(i,N)=UhNew(i,N-1)
! 	TeNew(i,N)=TeNew(i,N-1)
! 	ThNew(i,N)=ThNew(i,N-1)
! 	TsNew(i,N)=TsNew(i,N-1)
	
! 	potential(i,1)=0d0 !(0d0,0d0)
!                potential(i,N)=0d0 !(0d0,0d0)
    end do
    
    do i=2,M-1
	
	 ! boundary condition v.n = 0 on boundaries. 
	! NORTH
	
	 GradNeX(i,N) = GradNeX(i,N-1)
	 GradNeY(i,N) = GradNeY(i,N-1)
	 
! !       	 NeNew(i,N) = ((GainsE(i,N)-LossesE(i,N))*CellVol(i,N) & 
! ! 	      -( &
! ! ! 		(0.5d0*(JeX(i+1,N)+JeX(i,N))*NormalEx(i,N)+0.5d0*(JeY(i+1,N)+JeY(i,N))*NormalEy(i,N))*CellAreaE(i,N) & 
! ! ! 	       +(0.5d0*(JeX(i,N)+JeX(i-1,N))*NormalWx(i,N)+0.5d0*(JeY(i,N)+JeY(i-1,N))*NormalWy(i,N))*CellAreaW(i,N) &
! ! ! 	       +(0.5d0*(JeX(i,N)+JeX(i,N-1))*NormalSx(i,N)+0.5d0*(JeY(i,N)+JeY(i,N-1))*NormalSy(i,N))*CellAreaS(i,N) &
! ! 		 (0.5d0*(JeX(i-1,N)+JeX(i,N))*TangentWx(i,N)+0.5d0*(JeY(i-1,N)+JeY(i,N))*TangentWy(i,N))*TangentWx(i,N)*NormalWx(i,N)*CellAreaW(i,N) &
! ! 		+(0.5d0*(JeX(i-1,N)+JeX(i,N))*TangentWx(i,N)+0.5d0*(JeY(i-1,N)+JeY(i,N))*TangentWy(i,N))*TangentWy(i,N)*NormalWy(i,N)*CellAreaW(i,N) &
! ! 		+(0.5d0*(JeX(i,N)+JeX(i+1,N))*TangentEx(i,N)+0.5d0*(JeY(i,N)+JeY(i+1,N))*TangentEy(i,N))*TangentEx(i,N)*NormalEx(i,N)*CellAreaE(i,N) &
! ! 		+(0.5d0*(JeX(i,N)+JeX(i+1,N))*TangentEx(i,N)+0.5d0*(JeY(i,N)+JeY(i+1,N))*TangentEy(i,N))*TangentEy(i,N)*NormalEy(i,N)*CellAreaE(i,N) &
! ! 		+(0.5d0*(JeX(i,N)+JeX(i,N-1))*TangentSx(i,N)+0.5d0*(JeY(i,N)+JeY(i,N-1))*TangentSy(i,N))*TangentSx(i,N)*NormalSx(i,N)*CellAreaS(i,N) &
! ! 		+(0.5d0*(JeX(i,N)+JeX(i,N-1))*TangentSx(i,N)+0.5d0*(JeY(i,N)+JeY(i,N-1))*TangentSy(i,N))*TangentSy(i,N)*NormalSy(i,N)*CellAreaS(i,N) &
! ! 		+(0.5d0*(JeX(i-1,N)+JeX(i,N))*NormalWx(i,N)+0.5d0*(JeY(i-1,N)+JeY(i,N))*NormalWy(i,N))*NormalWx(i,N)*NormalWx(i,N)*CellAreaW(i,N) &
! ! 		+(0.5d0*(JeX(i-1,N)+JeX(i,N))*NormalWx(i,N)+0.5d0*(JeY(i-1,N)+JeY(i,N))*NormalWy(i,N))*NormalWy(i,N)*NormalWy(i,N)*CellAreaW(i,N) &
! ! 		+(0.5d0*(JeX(i,N)+JeX(i+1,N))*NormalEx(i,N)+0.5d0*(JeY(i,N)+JeY(i+1,N))*NormalEy(i,N))*NormalEx(i,N)*NormalEx(i,N)*CellAreaE(i,N) &
! ! 		+(0.5d0*(JeX(i,N)+JeX(i+1,N))*NormalEx(i,N)+0.5d0*(JeY(i,N)+JeY(i+1,N))*NormalEy(i,N))*NormalEy(i,N)*NormalEy(i,N)*CellAreaE(i,N) &
! ! 		+(0.5d0*(JeX(i,N)+JeX(i,N-1))*NormalSx(i,N)+0.5d0*(JeY(i,N)+JeY(i,N-1))*NormalSy(i,N))*NormalSx(i,N)*NormalSx(i,N)*CellAreaS(i,N) &
! ! 		+(0.5d0*(JeX(i,N)+JeX(i,N-1))*NormalSx(i,N)+0.5d0*(JeY(i,N)+JeY(i,N-1))*NormalSy(i,N))*NormalSy(i,N)*NormalSy(i,N)*CellAreaS(i,N) &
! ! 	       ) &
! ! 	     )*dt/CellVol(i,N) & 
! ! 	     +(0.5d0*(NormalEx(i,N)**2+NormalEy(i,N)**2)*CellAreaE(i,N)*(diffusionE(i,N)+diffusionE(i+1,N))*(Ne(i+1,N)-Ne(i,N))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistE(i,N) & 
! ! 		     + 0.5d0*(NormalWx(i,N)**2+NormalWy(i,N)**2)*CellAreaW(i,N)*(diffusionE(i-1,N)+diffusionE(i,N))*(Ne(i,N)-Ne(i-1,N))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistW(i,N) & 
! ! 		     + 0.5d0*(NormalSx(i,N)**2+NormalSy(i,N)**2)*CellAreaS(i,N)*(diffusionE(i,N-1)+diffusionE(i,N))*(Ne(i,N)-Ne(i,N-1))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistS(i,N) & 
! ! 		      -(CurviEx(i,N)*TangentEx(i,N)+CurviEy(i,N)*TangentEy(i,N))*CellAreaE(i,N)*(0.75d0*diffusionE(i,N)+0.75d0*diffusionE(i+1,N)+0.25d0*diffusionE(i+1,N-1)+0.25d0*diffusionE(i,N-1))*(0.25d0*Ne(i+1,N)+0.25d0*Ne(i,N)-0.25d0*Ne(i+1,N-1)-0.25d0*Ne(i,N-1))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistDualE(i,N) &
! ! 		      -(CurviWx(i,N)*TangentWx(i,N)+CurviWy(i,N)*TangentWy(i,N))*CellAreaW(i,N)*(0.75d0*diffusionE(i-1,N)+0.75d0*diffusionE(i,N)+0.25d0*diffusionE(i-1,N-1)+0.25d0*diffusionE(i,N-1))*(0.25d0*Ne(i,N)+0.25d0*Ne(i-1,N)-0.25d0*Ne(i-1,N-1)-0.25d0*Ne(i,N-1))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistDualW(i,N) &
! ! 		      -(CurviNx(i,N)*TangentNx(i,N)+CurviNy(i,N)*TangentNy(i,N))*CellAreaN(i,N)*diffusionE(i,N)*(0.5d0*Ne(i+1,N)-0.5d0*Ne(i-1,N))/(CurviNx(i,N)*NormalNx(i,N)+CurviNy(i,N)*NormalNy(i,N))/DistDualN(i,N) & 
! ! 		      -0.5d0*(CurviSx(i,N)*TangentSx(i,N)+CurviSy(i,N)*TangentSy(i,N))*CellAreaS(i,N)*(diffusionE(i,N-1)+diffusionE(i,N))*(0.25d0*Ne(i+1,N)+0.25d0*Ne(i+1,N-1)-0.25d0*Ne(i-1,N-1)-0.25d0*Ne(i-1,N))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistDualS(i,N) & 
! ! 		      )/CellVol(i,N)*dt &
! ! 	     +Ne(i,N)
! 	     
! !       	 NhNew(i,N) = ((GainsH(i,N)-LossesH(i,N))*CellVol(i,N) & 
! ! 	      -( &
! ! 		 (0.5d0*(JhX(i-1,N)+JhX(i,N))*TangentWx(i,N)+0.5d0*(JhY(i-1,N)+JhY(i,N))*TangentWy(i,N))*TangentWx(i,N)*NormalWx(i,N)*CellAreaW(i,N) &
! ! 		+(0.5d0*(JhX(i-1,N)+JhX(i,N))*TangentWx(i,N)+0.5d0*(JhY(i-1,N)+JhY(i,N))*TangentWy(i,N))*TangentWy(i,N)*NormalWy(i,N)*CellAreaW(i,N) &
! ! 		+(0.5d0*(JhX(i,N)+JhX(i+1,N))*TangentEx(i,N)+0.5d0*(JhY(i,N)+JhY(i+1,N))*TangentEy(i,N))*TangentEx(i,N)*NormalEx(i,N)*CellAreaE(i,N) &
! ! 		+(0.5d0*(JhX(i,N)+JhX(i+1,N))*TangentEx(i,N)+0.5d0*(JhY(i,N)+JhY(i+1,N))*TangentEy(i,N))*TangentEy(i,N)*NormalEy(i,N)*CellAreaE(i,N) &
! ! 		+(0.5d0*(JhX(i,N)+JhX(i,N-1))*TangentSx(i,N)+0.5d0*(JhY(i,N)+JhY(i,N-1))*TangentSy(i,N))*TangentSx(i,N)*NormalSx(i,N)*CellAreaS(i,N) &
! ! 		+(0.5d0*(JhX(i,N)+JhX(i,N-1))*TangentSx(i,N)+0.5d0*(JhY(i,N)+JhY(i,N-1))*TangentSy(i,N))*TangentSy(i,N)*NormalSy(i,N)*CellAreaS(i,N) &
! ! 		+(0.5d0*(JhX(i-1,N)+JhX(i,N))*NormalWx(i,N)+0.5d0*(JhY(i-1,N)+JhY(i,N))*NormalWy(i,N))*NormalWx(i,N)*NormalWx(i,N)*CellAreaW(i,N) &
! ! 		+(0.5d0*(JhX(i-1,N)+JhX(i,N))*NormalWx(i,N)+0.5d0*(JhY(i-1,N)+JhY(i,N))*NormalWy(i,N))*NormalWy(i,N)*NormalWy(i,N)*CellAreaW(i,N) &
! ! 		+(0.5d0*(JhX(i,N)+JhX(i+1,N))*NormalEx(i,N)+0.5d0*(JhY(i,N)+JhY(i+1,N))*NormalEy(i,N))*NormalEx(i,N)*NormalEx(i,N)*CellAreaE(i,N) &
! ! 		+(0.5d0*(JhX(i,N)+JhX(i+1,N))*NormalEx(i,N)+0.5d0*(JhY(i,N)+JhY(i+1,N))*NormalEy(i,N))*NormalEy(i,N)*NormalEy(i,N)*CellAreaE(i,N) &
! ! 		+(0.5d0*(JhX(i,N)+JhX(i,N-1))*NormalSx(i,N)+0.5d0*(JhY(i,N)+JhY(i,N-1))*NormalSy(i,N))*NormalSx(i,N)*NormalSx(i,N)*CellAreaS(i,N) &
! ! 		+(0.5d0*(JhX(i,N)+JhX(i,N-1))*NormalSx(i,N)+0.5d0*(JhY(i,N)+JhY(i,N-1))*NormalSy(i,N))*NormalSy(i,N)*NormalSy(i,N)*CellAreaS(i,N) &
! ! ! 		(0.5d0*(JhX(i+1,N)+JhX(i,N))*NormalEx(i,N)+0.5d0*(JhY(i+1,N)+JhY(i,N))*NormalEy(i,N))*CellAreaE(i,N) & 
! ! ! 	       +(0.5d0*(JhX(i,N)+JhX(i-1,N))*NormalWx(i,N)+0.5d0*(JhY(i,N)+JhY(i-1,N))*NormalWy(i,N))*CellAreaW(i,N) &
! ! ! 	       +(0.5d0*(JhX(i,N)+JhX(i,N-1))*NormalSx(i,N)+0.5d0*(JhY(i,N)+JhY(i,N-1))*NormalSy(i,N))*CellAreaS(i,N) &
! ! 	       ) &
! ! 	     )*dt/CellVol(i,N) & 
! ! 	     +(0.5d0*(NormalEx(i,N)**2+NormalEy(i,N)**2)*CellAreaE(i,N)*(diffusionH(i,N)+diffusionH(i+1,N))*(Nh(i+1,N)-Nh(i,N))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistE(i,N) & 
! ! 		     + 0.5d0*(NormalWx(i,N)**2+NormalWy(i,N)**2)*CellAreaW(i,N)*(diffusionH(i-1,N)+diffusionH(i,N))*(Nh(i,N)-Nh(i-1,N))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistW(i,N) & 
! ! 		     + 0.5d0*(NormalSx(i,N)**2+NormalSy(i,N)**2)*CellAreaS(i,N)*(diffusionH(i,N-1)+diffusionH(i,N))*(Nh(i,N)-Nh(i,N-1))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistS(i,N) & 
! ! 		      -(CurviEx(i,N)*TangentEx(i,N)+CurviEy(i,N)*TangentEy(i,N))*CellAreaE(i,N)*(0.75d0*diffusionH(i,N)+0.75d0*diffusionH(i+1,N)+0.25d0*diffusionH(i+1,N-1)+0.25d0*diffusionH(i,N-1))*(0.25d0*Nh(i+1,N)+0.25d0*Nh(i,N)-0.25d0*Nh(i+1,N-1)-0.25d0*Nh(i,N-1))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistDualE(i,N) &
! ! 		      -(CurviWx(i,N)*TangentWx(i,N)+CurviWy(i,N)*TangentWy(i,N))*CellAreaW(i,N)*(0.75d0*diffusionH(i-1,N)+0.75d0*diffusionH(i,N)+0.25d0*diffusionH(i-1,N-1)+0.25d0*diffusionH(i,N-1))*(0.25d0*Nh(i,N)+0.25d0*Nh(i-1,N)-0.25d0*Nh(i-1,N-1)-0.25d0*Nh(i,N-1))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistDualW(i,N) &
! ! 		      -(CurviNx(i,N)*TangentNx(i,N)+CurviNy(i,N)*TangentNy(i,N))*CellAreaN(i,N)*diffusionH(i,N)*(0.5d0*Nh(i+1,N)-0.5d0*Nh(i-1,N))/(CurviNx(i,N)*NormalNx(i,N)+CurviNy(i,N)*NormalNy(i,N))/DistDualN(i,N) & 
! ! 		      -0.5d0*(CurviSx(i,N)*TangentSx(i,N)+CurviSy(i,N)*TangentSy(i,N))*CellAreaS(i,N)*(diffusionH(i,N-1)+diffusionH(i,N))*(0.25d0*Nh(i+1,N)+0.25d0*Nh(i+1,N-1)-0.25d0*Nh(i-1,N-1)-0.25d0*Nh(i-1,N))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistDualS(i,N) & 
! ! 		      )/CellVol(i,N)*dt &
! ! 	     +Nh(i,N)
! 	     
! 	 UeNew(i,N) = ((SourceUe(i,N)-CouplingE(i,N))*CellVol(i,N) & 
! 	      -( &
! ! 		(0.5d0*(JeX(i+1,N)+JeX(i,N))*NormalEx(i,N)+0.5d0*(JeY(i+1,N)+JeY(i,N))*NormalEy(i,N))*CellAreaE(i,N) & 
! ! 	       +(0.5d0*(JeX(i,N)+JeX(i-1,N))*NormalWx(i,N)+0.5d0*(JeY(i,N)+JeY(i-1,N))*NormalWy(i,N))*CellAreaW(i,N) &
! ! 	       +(0.5d0*(JeX(i,N)+JeX(i,N-1))*NormalSx(i,N)+0.5d0*(JeY(i,N)+JeY(i,N-1))*NormalSy(i,N))*CellAreaS(i,N) &
! 		 (0.5d0*(VeX(i-1,N)+VeX(i,N))*TangentWx(i,N)+0.5d0*(VeY(i-1,N)+VeY(i,N))*TangentWy(i,N))*TangentWx(i,N)*NormalWx(i,N)*CellAreaW(i,N) &
! 		+(0.5d0*(VeX(i-1,N)+VeX(i,N))*TangentWx(i,N)+0.5d0*(VeY(i-1,N)+VeY(i,N))*TangentWy(i,N))*TangentWy(i,N)*NormalWy(i,N)*CellAreaW(i,N) &
! 		+(0.5d0*(VeX(i,N)+VeX(i+1,N))*TangentEx(i,N)+0.5d0*(VeY(i,N)+VeY(i+1,N))*TangentEy(i,N))*TangentEx(i,N)*NormalEx(i,N)*CellAreaE(i,N) &
! 		+(0.5d0*(VeX(i,N)+VeX(i+1,N))*TangentEx(i,N)+0.5d0*(VeY(i,N)+VeY(i+1,N))*TangentEy(i,N))*TangentEy(i,N)*NormalEy(i,N)*CellAreaE(i,N) &
! 		+(0.5d0*(VeX(i,N)+VeX(i,N-1))*TangentSx(i,N)+0.5d0*(VeY(i,N)+VeY(i,N-1))*TangentSy(i,N))*TangentSx(i,N)*NormalSx(i,N)*CellAreaS(i,N) &
! 		+(0.5d0*(VeX(i,N)+VeX(i,N-1))*TangentSx(i,N)+0.5d0*(VeY(i,N)+VeY(i,N-1))*TangentSy(i,N))*TangentSy(i,N)*NormalSy(i,N)*CellAreaS(i,N) &
! 	       ) &
! 	     )*dt/CellVol(i,N)+Ue(i,N)     
! 	     
! 	UhNew(i,N) = ((SourceUh(i,N)-CouplingH(i,N))*CellVol(i,N) & 
! 	      -( &
! ! 		(0.5d0*(JeX(i+1,N)+JeX(i,N))*NormalEx(i,N)+0.5d0*(JeY(i+1,N)+JeY(i,N))*NormalEy(i,N))*CellAreaE(i,N) & 
! ! 	       +(0.5d0*(JeX(i,N)+JeX(i-1,N))*NormalWx(i,N)+0.5d0*(JeY(i,N)+JeY(i-1,N))*NormalWy(i,N))*CellAreaW(i,N) &
! ! 	       +(0.5d0*(JeX(i,N)+JeX(i,N-1))*NormalSx(i,N)+0.5d0*(JeY(i,N)+JeY(i,N-1))*NormalSy(i,N))*CellAreaS(i,N) &
! 		 (0.5d0*(VhX(i-1,N)+VhX(i,N))*TangentWx(i,N)+0.5d0*(VhY(i-1,N)+VhY(i,N))*TangentWy(i,N))*TangentWx(i,N)*NormalWx(i,N)*CellAreaW(i,N) &
! 		+(0.5d0*(VhX(i-1,N)+VhX(i,N))*TangentWx(i,N)+0.5d0*(VhY(i-1,N)+VhY(i,N))*TangentWy(i,N))*TangentWy(i,N)*NormalWy(i,N)*CellAreaW(i,N) &
! 		+(0.5d0*(VhX(i,N)+VhX(i+1,N))*TangentEx(i,N)+0.5d0*(VhY(i,N)+VhY(i+1,N))*TangentEy(i,N))*TangentEx(i,N)*NormalEx(i,N)*CellAreaE(i,N) &
! 		+(0.5d0*(VhX(i,N)+VhX(i+1,N))*TangentEx(i,N)+0.5d0*(VhY(i,N)+VhY(i+1,N))*TangentEy(i,N))*TangentEy(i,N)*NormalEy(i,N)*CellAreaE(i,N) &
! 		+(0.5d0*(VhX(i,N)+VhX(i,N-1))*TangentSx(i,N)+0.5d0*(VhY(i,N)+VhY(i,N-1))*TangentSy(i,N))*TangentSx(i,N)*NormalSx(i,N)*CellAreaS(i,N) &
! 		+(0.5d0*(VhX(i,N)+VhX(i,N-1))*TangentSx(i,N)+0.5d0*(VhY(i,N)+VhY(i,N-1))*TangentSy(i,N))*TangentSy(i,N)*NormalSy(i,N)*CellAreaS(i,N) &
! 	       ) &
! 	     )*dt/CellVol(i,N)+Uh(i,N)   
! 	     
! 	 ! diffusion on irregular mesh    
! 	 if(TeOff.ne.1) then
! 	 TeNew(i,N) = (0.5d0*(NormalEx(i,N)**2+NormalEy(i,N)**2)*CellAreaE(i,N)*(kappae(i,N)+kappae(i+1,N))*(Te(i+1,N)-Te(i,N))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistE(i,N) & 
! 		     + 0.5d0*(NormalWx(i,N)**2+NormalWy(i,N)**2)*CellAreaW(i,N)*(kappae(i-1,N)+kappae(i,N))*(Te(i,N)-Te(i-1,N))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistW(i,N) & 
! 		     + 0.5d0*(NormalSx(i,N)**2+NormalSy(i,N)**2)*CellAreaS(i,N)*(kappae(i,N-1)+kappae(i,N))*(Te(i,N)-Te(i,N-1))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistS(i,N) & 
! 		      -(CurviEx(i,N)*TangentEx(i,N)+CurviEy(i,N)*TangentEy(i,N))*CellAreaE(i,N)*(0.75d0*kappae(i,N)+0.75d0*kappae(i+1,N)+0.25d0*kappae(i+1,N-1)+0.25d0*kappae(i,N-1))*(0.25d0*Te(i+1,N)+0.25d0*Te(i,N)-0.25d0*Te(i+1,N-1)-0.25d0*Te(i,N-1))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistDualE(i,N) &
! 		      -(CurviWx(i,N)*TangentWx(i,N)+CurviWy(i,N)*TangentWy(i,N))*CellAreaW(i,N)*(0.75d0*kappae(i-1,N)+0.75d0*kappae(i,N)+0.25d0*kappae(i-1,N-1)+0.25d0*kappae(i,N-1))*(0.25d0*Te(i,N)+0.25d0*Te(i-1,N)-0.25d0*Te(i-1,N-1)-0.25d0*Te(i,N-1))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistDualW(i,N) &
! 		      -(CurviNx(i,N)*TangentNx(i,N)+CurviNy(i,N)*TangentNy(i,N))*CellAreaN(i,N)*kappae(i,N)*(0.5d0*Te(i+1,N)-0.5d0*Te(i-1,N))/(CurviNx(i,N)*NormalNx(i,N)+CurviNy(i,N)*NormalNy(i,N))/DistDualN(i,N) & 
! 		      -0.5d0*(CurviSx(i,N)*TangentSx(i,N)+CurviSy(i,N)*TangentSy(i,N))*CellAreaS(i,N)*(kappae(i,N-1)+kappae(i,N))*(0.25d0*Te(i+1,N)+0.25d0*Te(i+1,N-1)-0.25d0*Te(i-1,N-1)-0.25d0*Te(i-1,N))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistDualS(i,N) & 
! 		      +(-CouplingE(i,N)+SourceE(i,N))*CellVol(i,N) &
! 		      )/Ce(i,N)/CellVol(i,N)*dt+Te(i,N)
! 		      
! 		      
! 	ThNew(i,N) = (0.5d0*(NormalEx(i,N)**2+NormalEy(i,N)**2)*CellAreaE(i,N)*(kappah(i,N)+kappah(i+1,N))*(Th(i+1,N)-Th(i,N))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistE(i,N) & 
! 		     + 0.5d0*(NormalWx(i,N)**2+NormalWy(i,N)**2)*CellAreaW(i,N)*(kappah(i-1,N)+kappah(i,N))*(Th(i,N)-Th(i-1,N))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistW(i,N) & 
! 		     + 0.5d0*(NormalSx(i,N)**2+NormalSy(i,N)**2)*CellAreaS(i,N)*(kappah(i,N-1)+kappah(i,N))*(Th(i,N)-Th(i,N-1))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistS(i,N) & 
! 		      -(CurviEx(i,N)*TangentEx(i,N)+CurviEy(i,N)*TangentEy(i,N))*CellAreaE(i,N)*(0.75d0*kappah(i,N)+0.75d0*kappah(i+1,N)+0.25d0*kappah(i+1,N-1)+0.25d0*kappah(i,N-1))*(0.25d0*Th(i+1,N)+0.25d0*Th(i,N)-0.25d0*Th(i+1,N-1)-0.25d0*Th(i,N-1))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistDualE(i,N) &
! 		      -(CurviWx(i,N)*TangentWx(i,N)+CurviWy(i,N)*TangentWy(i,N))*CellAreaW(i,N)*(0.75d0*kappah(i-1,N)+0.75d0*kappah(i,N)+0.25d0*kappah(i-1,N-1)+0.25d0*kappah(i,N-1))*(0.25d0*Th(i,N)+0.25d0*Th(i-1,N)-0.25d0*Th(i-1,N-1)-0.25d0*Th(i,N-1))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistDualW(i,N) &
! 		      -(CurviNx(i,N)*TangentNx(i,N)+CurviNy(i,N)*TangentNy(i,N))*CellAreaN(i,N)*kappah(i,N)*(0.5d0*Th(i+1,N)-0.5d0*Th(i-1,N))/(CurviNx(i,N)*NormalNx(i,N)+CurviNy(i,N)*NormalNy(i,N))/DistDualN(i,N) & 
! 		      -0.5d0*(CurviSx(i,N)*TangentSx(i,N)+CurviSy(i,N)*TangentSy(i,N))*CellAreaS(i,N)*(kappah(i,N-1)+kappah(i,N))*(0.25d0*Th(i+1,N)+0.25d0*Th(i+1,N-1)-0.25d0*Th(i-1,N-1)-0.25d0*Th(i-1,N))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistDualS(i,N) & 
! 		      +(-CouplingH(i,N)+SourceH(i,N))*CellVol(i,N) &
! 		      )/Ch(i,N)/CellVol(i,N)*dt+Th(i,N)
! 	      
! ! 	 TeNew(i,N) =  ! bad one &
! ! 		       (0.5d0*(NormalEx(i,N)**2+NormalEy(i,N)**2)*CellAreaE(i,N)*(kappae( i,N )+kappae(i+1,N))*(Te(i+1,N)-Te(i,N)) / (CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistE(i,N) &
! ! 		      + 0.5d0*(NormalWx(i,N)**2+NormalWy(i,N)**2)*CellAreaW(i,N)*(kappae(i-1,N)+kappae( i,N))*(Te(i,N)-Te(i-1,N))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistW(i,N) & 
! ! 		      + 0.5d0*(NormalSx(i,N)**2+NormalSy(i,N)**2)*CellAreaS(i,N)*(kappae(i,N-1)+kappae( i,N))*(Te(i,N)-Te(i,N-1))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistS(i,N) &
! ! 		      - 0.5d0*(CurviEx(i,N)*TangentEx(i,N)+CurviEy(i,N)*TangentEy(i,N))*CellAreaE(i,N)*(kappae(i,N)+kappae(i+1,N))*(0.25d0*Te(i+1,N)+0.25d0*Te( i,N  )  - 0.25d0*Te(i+1,N-1)-0.25d0*Te(i,N-1))/(CurviEx(i,N)*NormalEx(i,N)+CurviEy(i,N)*NormalEy(i,N))/DistDualE(i,N) & 
! ! 		      - 0.5d0*(CurviWx(i,N)*TangentWx(i,N)+CurviWy(i,N)*TangentWy(i,N))*CellAreaW(i,N)*(kappae(i-1,N)+kappae(i,N))*(0.25d0*Te(i,N ) +0.25d0*Te( i-1,N ) - 0.25d0*Te(i-1,N-1)-0.25d0*Te(i,N-1))/(CurviWx(i,N)*NormalWx(i,N)+CurviWy(i,N)*NormalWy(i,N))/DistDualW(i,N) & 
! ! 		      - 0.5d0*(CurviSx(i,N)*TangentSx(i,N)+CurviSy(i,N)*TangentSy(i,N))*CellAreaS(i,N)*(kappae(i,N-1)+kappae(i,N))*(0.25d0*Te(i+1,N)+0.25d0*Te(i+1,N-1) - 0.25d0*Te(i-1,N-1)-0.25d0*Te(i-1,N))/(CurviSx(i,N)*NormalSx(i,N)+CurviSy(i,N)*NormalSy(i,N))/DistDualS(i,N) & 
! ! 		      - (CurviNx(i,N)*TangentNx(i,N)+CurviNy(i,N)*TangentNy(i,N))*CellAreaN(i,N)*kappae(i,N)*(0.5d0*Te(i+1,N)-0.5d0*Te(i-1,N))/(CurviNx(i,N)*NormalNx(i,N)+CurviNy(i,N)*NormalNy(i,N))/DistDualN(i,N) & 
! ! 		      +(-CouplingE(i,N)+SourceE(i,N))*CellVol(i,N) &
! ! 		      )/Ce(i,N)*dt/CellVol(i,N) & 
! ! 		      +Te(i,N)
! 	end if
! 
	! SOUTH
	 GradNeX(i,1) = GradNeX(i,2)
	 GradNeY(i,1) = GradNeY(i,2)
! 	 
! !       	 NeNew(i,1) = ((GainsE(i,1)-LossesE(i,1))*CellVol(i,1) & 
! ! 	      -( &
! ! ! 		(0.5d0*(JeX(i,1)+JeX(i,2))*NormalNx(i,1)+0.5d0*(JeY(i,2)+JeY(i,1))*NormalNy(i,1))*CellAreaN(i,1) & 
! ! ! 	       +(0.5d0*(JeX(i+1,1)+JeX(i,1))*NormalEx(i,1)+0.5d0*(JeY(i+1,1)+JeY(i,1))*NormalEy(i,1))*CellAreaE(i,1) & 
! ! ! 	       +(0.5d0*(JeX(i,1)+JeX(i-1,1))*NormalWx(i,1)+0.5d0*(JeY(i,1)+JeY(i-1,1))*NormalWy(i,1))*CellAreaW(i,1) &
! ! 		
! ! 		 (0.5d0*(JeX(i-1,1)+JeX(i,1))*TangentWx(i,1)+0.5d0*(JeY(i-1,1)+JeY(i,1))*TangentWy(i,1))*TangentWx(i,1)*NormalWx(i,1)*CellAreaW(i,1) &
! ! 		+(0.5d0*(JeX(i-1,1)+JeX(i,1))*TangentWx(i,1)+0.5d0*(JeY(i-1,1)+JeY(i,1))*TangentWy(i,1))*TangentWy(i,1)*NormalWy(i,1)*CellAreaW(i,1) &
! ! 		+(0.5d0*(JeX(i,1)+JeX(i+1,1))*TangentEx(i,1)+0.5d0*(JeY(i,1)+JeY(i+1,1))*TangentEy(i,1))*TangentEx(i,1)*NormalEx(i,1)*CellAreaE(i,1) &
! ! 		+(0.5d0*(JeX(i,1)+JeX(i+1,1))*TangentEx(i,1)+0.5d0*(JeY(i,1)+JeY(i+1,1))*TangentEy(i,1))*TangentEy(i,1)*NormalEy(i,1)*CellAreaE(i,1) &
! ! 		+(0.5d0*(JeX(i,1)+JeX(i,2  ))*TangentNx(i,1)+0.5d0*(JeY(i,1)+JeY(i,2  ))*TangentNy(i,1))*TangentNx(i,1)*NormalNx(i,1)*CellAreaN(i,1) &
! ! 		+(0.5d0*(JeX(i,1)+JeX(i,2  ))*TangentNx(i,1)+0.5d0*(JeY(i,1)+JeY(i,2  ))*TangentNy(i,1))*TangentNy(i,1)*NormalNy(i,1)*CellAreaN(i,1) &
! ! 		+(0.5d0*(JeX(i-1,1)+JeX(i,1))*NormalWx(i,1)+0.5d0*(JeY(i-1,1)+JeY(i,1))*NormalWy(i,1))*NormalWx(i,1)*NormalWx(i,1)*CellAreaW(i,1) &
! ! 		+(0.5d0*(JeX(i-1,1)+JeX(i,1))*NormalWx(i,1)+0.5d0*(JeY(i-1,1)+JeY(i,1))*NormalWy(i,1))*NormalWy(i,1)*NormalWy(i,1)*CellAreaW(i,1) &
! ! 		+(0.5d0*(JeX(i,1)+JeX(i+1,1))*NormalEx(i,1)+0.5d0*(JeY(i,1)+JeY(i+1,1))*NormalEy(i,1))*NormalEx(i,1)*NormalEx(i,1)*CellAreaE(i,1) &
! ! 		+(0.5d0*(JeX(i,1)+JeX(i+1,1))*NormalEx(i,1)+0.5d0*(JeY(i,1)+JeY(i+1,1))*NormalEy(i,1))*NormalEy(i,1)*NormalEy(i,1)*CellAreaE(i,1) &
! ! 		+(0.5d0*(JeX(i,1)+JeX(i,2  ))*NormalNx(i,1)+0.5d0*(JeY(i,1)+JeY(i,2  ))*NormalNy(i,1))*NormalNx(i,1)*NormalNx(i,1)*CellAreaN(i,1) &
! ! 		+(0.5d0*(JeX(i,1)+JeX(i,2  ))*NormalNx(i,1)+0.5d0*(JeY(i,1)+JeY(i,2  ))*NormalNy(i,1))*NormalNy(i,1)*NormalNy(i,1)*CellAreaN(i,1) &
! ! 	       ) &
! ! 	     )*dt/CellVol(i,1) & 
! ! 	     +(&
! ! 	       0.5d0*(NormalNx(i,1)**2+NormalNy(i,1)**2)*CellAreaN(i,1)*(diffusionE(i,1)+diffusionE(i  ,2))*(Ne(i,2  )-Ne(i,1))/(CurviNx(i,1)*NormalNx(i,1)+CurviNy(i,1)*NormalNy(i,1))/DistN(i,1) &
! ! 	      +0.5d0*(NormalEx(i,1)**2+NormalEy(i,1)**2)*CellAreaE(i,1)*(diffusionE(i,1)+diffusionE(i+1,1))*(Ne(i+1,1)-Ne(i,1))/(CurviEx(i,1)*NormalEx(i,1)+CurviEy(i,1)*NormalEy(i,1))/DistE(i,1) & 
! ! 	      +0.5d0*(NormalWx(i,1)**2+NormalWy(i,1)**2)*CellAreaW(i,1)*(diffusionE(i,1)+diffusionE(i-1,1))*(Ne(i,1)-Ne(i-1,1))/(CurviWx(i,1)*NormalWx(i,1)+CurviWy(i,1)*NormalWy(i,1))/DistW(i,1) &
! ! 	      -0.5d0*(CurviEx(i,1)*TangentEx(i,1)+CurviEy(i,1)*TangentEy(i,1))*CellAreaE(i,1)*(0.25d0*diffusionE(i+1,2)+0.75d0*diffusionE(i+1,1)+0.75d0*diffusionE(i,1)+0.25d0*diffusionE(i,2))*(0.25d0*Ne(i+1,2)-0.25d0*Ne(i+1,1)-0.25d0*Ne(i,1)+0.25d0*Ne(i,2))/(CurviEx(i,1)*NormalEx(i,1)+CurviEy(i,1)*NormalEy(i,1))/DistDualE(i,1) &
! ! 	      -0.5d0*(CurviWx(i,1)*TangentWx(i,1)+CurviWy(i,1)*TangentWy(i,1))*CellAreaW(i,1)*(0.75d0*diffusionE(i,1)+0.25d0*diffusionE(i,2)+0.75d0*diffusionE(i-1,1)+0.25d0*diffusionE(i-1,2))*(-0.25d0*Ne(i,1)+0.25d0*Ne(i,2)-0.25d0*Ne(i-1,1)+0.25d0*Ne(i-1,2))/(CurviWx(i,1)*NormalWx(i,1)+CurviWy(i,1)*NormalWy(i,1))/DistDualW(i,1) &
! ! 	      -0.5d0*(CurviNx(i,1)*TangentNx(i,1)+CurviNy(i,1)*TangentNy(i,1))*CellAreaN(i,1)*(0.25d0*diffusionE(i+1,2)+0.25d0*diffusionE(i+1,1)+0.5d0*diffusionE(i,1)+0.5d0*diffusionE(i,2)+0.25d0*diffusionE(i-1,1)+0.25d0*diffusionE(i-1,2))*(0.25d0*Ne(i+1,2)+0.25d0*Ne(i+1,1)-0.25d0*Ne(i-1,1)-0.25d0*Ne(i-1,2))/(CurviNx(i,1)*NormalNx(i,1)+CurviNy(i,1)*NormalNy(i,1))/DistDualN(i,1) &
! ! ! 	      -(CurviSx(i,1)*TangentSx(i,1)+CurviSy(i,1)*TangentSy(i,1))*CellAreaS(i,1)*(diffusionE(i,1)+0.5d0*diffusionE(i+1,1)+0.5d0*diffusionE(i-1,1))*(0.5d0*Ne(i+1,1)-0.5d0*Ne(i-1,1))/(CurviSx(i,1)*NormalSx(i,1)+CurviSy(i,1)*NormalSy(i,1))/DistDualS(i,1) &
! ! 	      )/CellVol(i,1)*dt &
! ! 	     +Ne(i,1)
! ! 	     
! ! 	 NhNew(i,1) = ((GainsH(i,1)-LossesH(i,1))*CellVol(i,1) & 
! ! 	      -( &
! ! 		 (0.5d0*(JhX(i-1,1)+JhX(i,1))*TangentWx(i,1)+0.5d0*(JhY(i-1,1)+JhY(i,1))*TangentWy(i,1))*TangentWx(i,1)*NormalWx(i,1)*CellAreaW(i,1) &
! ! 		+(0.5d0*(JhX(i-1,1)+JhX(i,1))*TangentWx(i,1)+0.5d0*(JhY(i-1,1)+JhY(i,1))*TangentWy(i,1))*TangentWy(i,1)*NormalWy(i,1)*CellAreaW(i,1) &
! ! 		+(0.5d0*(JhX(i,1)+JhX(i+1,1))*TangentEx(i,1)+0.5d0*(JhY(i,1)+JhY(i+1,1))*TangentEy(i,1))*TangentEx(i,1)*NormalEx(i,1)*CellAreaE(i,1) &
! ! 		+(0.5d0*(JhX(i,1)+JhX(i+1,1))*TangentEx(i,1)+0.5d0*(JhY(i,1)+JhY(i+1,1))*TangentEy(i,1))*TangentEy(i,1)*NormalEy(i,1)*CellAreaE(i,1) &
! ! 		+(0.5d0*(JhX(i,1)+JhX(i,2  ))*TangentNx(i,1)+0.5d0*(JhY(i,1)+JhY(i,2  ))*TangentNy(i,1))*TangentNx(i,1)*NormalNx(i,1)*CellAreaN(i,1) &
! ! 		+(0.5d0*(JhX(i,1)+JhX(i,2  ))*TangentNx(i,1)+0.5d0*(JhY(i,1)+JhY(i,2  ))*TangentNy(i,1))*TangentNy(i,1)*NormalNy(i,1)*CellAreaN(i,1) &
! ! 		+(0.5d0*(JhX(i-1,1)+JhX(i,1))*NormalWx(i,1)+0.5d0*(JhY(i-1,1)+JhY(i,1))*NormalWy(i,1))*NormalWx(i,1)*NormalWx(i,1)*CellAreaW(i,1) &
! ! 		+(0.5d0*(JhX(i-1,1)+JhX(i,1))*NormalWx(i,1)+0.5d0*(JhY(i-1,1)+JhY(i,1))*NormalWy(i,1))*NormalWy(i,1)*NormalWy(i,1)*CellAreaW(i,1) &
! ! 		+(0.5d0*(JhX(i,1)+JhX(i+1,1))*NormalEx(i,1)+0.5d0*(JhY(i,1)+JhY(i+1,1))*NormalEy(i,1))*NormalEx(i,1)*NormalEx(i,1)*CellAreaE(i,1) &
! ! 		+(0.5d0*(JhX(i,1)+JhX(i+1,1))*NormalEx(i,1)+0.5d0*(JhY(i,1)+JhY(i+1,1))*NormalEy(i,1))*NormalEy(i,1)*NormalEy(i,1)*CellAreaE(i,1) &
! ! 		+(0.5d0*(JhX(i,1)+JhX(i,2  ))*NormalNx(i,1)+0.5d0*(JhY(i,1)+JhY(i,2  ))*NormalNy(i,1))*NormalNx(i,1)*NormalNx(i,1)*CellAreaN(i,1) &
! ! 		+(0.5d0*(JhX(i,1)+JhX(i,2  ))*NormalNx(i,1)+0.5d0*(JhY(i,1)+JhY(i,2  ))*NormalNy(i,1))*NormalNy(i,1)*NormalNy(i,1)*CellAreaN(i,1) &
! ! ! 		(0.5d0*(JhX(i,1)+JhX(i,2))*NormalNx(i,1)+0.5d0*(JhY(i,2)+JhY(i,1))*NormalNy(i,1))*CellAreaN(i,1) & 
! ! ! 	       +(0.5d0*(JhX(i+1,1)+JhX(i,1))*NormalEx(i,1)+0.5d0*(JhY(i+1,1)+JhY(i,1))*NormalEy(i,1))*CellAreaE(i,1) & 
! ! ! 	       +(0.5d0*(JhX(i,1)+JhX(i-1,1))*NormalWx(i,1)+0.5d0*(JhY(i,1)+JhY(i-1,1))*NormalWy(i,1))*CellAreaW(i,1) &
! ! 	       ) &
! ! 	     )*dt/CellVol(i,1) & 
! ! 	     +(&
! ! 	       0.5d0*(NormalNx(i,1)**2+NormalNy(i,1)**2)*CellAreaN(i,1)*(diffusionH(i,1)+diffusionH(i  ,2))*(Nh(i,2  )-Nh(i,1))/(CurviNx(i,1)*NormalNx(i,1)+CurviNy(i,1)*NormalNy(i,1))/DistN(i,1) &
! ! 	      +0.5d0*(NormalEx(i,1)**2+NormalEy(i,1)**2)*CellAreaE(i,1)*(diffusionH(i,1)+diffusionH(i+1,1))*(Nh(i+1,1)-Nh(i,1))/(CurviEx(i,1)*NormalEx(i,1)+CurviEy(i,1)*NormalEy(i,1))/DistE(i,1) & 
! ! 	      +0.5d0*(NormalWx(i,1)**2+NormalWy(i,1)**2)*CellAreaW(i,1)*(diffusionH(i,1)+diffusionH(i-1,1))*(Nh(i,1)-Nh(i-1,1))/(CurviWx(i,1)*NormalWx(i,1)+CurviWy(i,1)*NormalWy(i,1))/DistW(i,1) &
! ! 	      -0.5d0*(CurviEx(i,1)*TangentEx(i,1)+CurviEy(i,1)*TangentEy(i,1))*CellAreaE(i,1)*(0.25d0*diffusionH(i+1,2)+0.75d0*diffusionH(i+1,1)+0.75d0*diffusionH(i,1)+0.25d0*diffusionH(i,2))*(0.25d0*Nh(i+1,2)-0.25d0*Nh(i+1,1)-0.25d0*Nh(i,1)+0.25d0*Nh(i,2))/(CurviEx(i,1)*NormalEx(i,1)+CurviEy(i,1)*NormalEy(i,1))/DistDualE(i,1) &
! ! 	      -0.5d0*(CurviWx(i,1)*TangentWx(i,1)+CurviWy(i,1)*TangentWy(i,1))*CellAreaW(i,1)*(0.75d0*diffusionH(i,1)+0.25d0*diffusionH(i,2)+0.75d0*diffusionH(i-1,1)+0.25d0*diffusionH(i-1,2))*(-0.25d0*Nh(i,1)+0.25d0*Nh(i,2)-0.25d0*Nh(i-1,1)+0.25d0*Nh(i-1,2))/(CurviWx(i,1)*NormalWx(i,1)+CurviWy(i,1)*NormalWy(i,1))/DistDualW(i,1) &
! ! 	      -0.5d0*(CurviNx(i,1)*TangentNx(i,1)+CurviNy(i,1)*TangentNy(i,1))*CellAreaN(i,1)*(0.25d0*diffusionH(i+1,2)+0.25d0*diffusionH(i+1,1)+0.5d0*diffusionH(i,1)+0.5d0*diffusionH(i,2)+0.25d0*diffusionH(i-1,1)+0.25d0*diffusionH(i-1,2))*(0.25d0*Nh(i+1,2)+0.25d0*Nh(i+1,1)-0.25d0*Nh(i-1,1)-0.25d0*Nh(i-1,2))/(CurviNx(i,1)*NormalNx(i,1)+CurviNy(i,1)*NormalNy(i,1))/DistDualN(i,1) &
! ! ! 	      -(CurviSx(i,1)*TangentSx(i,1)+CurviSy(i,1)*TangentSy(i,1))*CellAreaS(i,1)*(diffusionH(i,1)+0.5d0*diffusionH(i+1,1)+0.5d0*diffusionH(i-1,1))*(0.5d0*Nh(i+1,1)-0.5d0*Nh(i-1,1))/(CurviSx(i,1)*NormalSx(i,1)+CurviSy(i,1)*NormalSy(i,1))/DistDualS(i,1) &
! ! 	      )/CellVol(i,1)*dt &
! ! 	     +Nh(i,1)
! ! 	SourceE(i,1)=0d0 
! 	if(TeOff.ne.1) then
! 	TeNew(i,1) = (&
! 	       0.5d0*(NormalNx(i,1)**2+NormalNy(i,1)**2)*CellAreaN(i,1)*(kappae(i,1)+kappae(i  ,2))*(Te(i,2  )-Te(i,1))/(CurviNx(i,1)*NormalNx(i,1)+CurviNy(i,1)*NormalNy(i,1))/DistN(i,1) &
! 	      +0.5d0*(NormalEx(i,1)**2+NormalEy(i,1)**2)*CellAreaE(i,1)*(kappae(i,1)+kappae(i+1,1))*(Te(i+1,1)-Te(i,1))/(CurviEx(i,1)*NormalEx(i,1)+CurviEy(i,1)*NormalEy(i,1))/DistE(i,1) & 
! 	      +0.5d0*(NormalWx(i,1)**2+NormalWy(i,1)**2)*CellAreaW(i,1)*(kappae(i,1)+kappae(i-1,1))*(Te(i,1)-Te(i-1,1))/(CurviWx(i,1)*NormalWx(i,1)+CurviWy(i,1)*NormalWy(i,1))/DistW(i,1) &
! 	      -0.5d0*(CurviEx(i,1)*TangentEx(i,1)+CurviEy(i,1)*TangentEy(i,1))*CellAreaE(i,1)*(0.25d0*kappae(i+1,2)+0.75d0*kappae(i+1,1)+0.75d0*kappae(i,1)+0.25d0*kappae(i,2))*(0.25d0*Te(i+1,2)-0.25d0*Te(i+1,1)-0.25d0*Te(i,1)+0.25d0*Te(i,2))/(CurviEx(i,1)*NormalEx(i,1)+CurviEy(i,1)*NormalEy(i,1))/DistDualE(i,1) &
! 	      -0.5d0*(CurviWx(i,1)*TangentWx(i,1)+CurviWy(i,1)*TangentWy(i,1))*CellAreaW(i,1)*(0.75d0*kappae(i,1)+0.25d0*kappae(i,2)+0.75d0*kappae(i-1,1)+0.25d0*kappae(i-1,2))*(-0.25d0*Te(i,1)+0.25d0*Te(i,2)-0.25d0*Te(i-1,1)+0.25d0*Te(i-1,2))/(CurviWx(i,1)*NormalWx(i,1)+CurviWy(i,1)*NormalWy(i,1))/DistDualW(i,1) &
! 	      -0.5d0*(CurviNx(i,1)*TangentNx(i,1)+CurviNy(i,1)*TangentNy(i,1))*CellAreaN(i,1)*(0.25d0*kappae(i+1,2)+0.25d0*kappae(i+1,1)+0.5d0*kappae(i,1)+0.5d0*kappae(i,2)+0.25d0*kappae(i-1,1)+0.25d0*kappae(i-1,2))*(0.25d0*Te(i+1,2)+0.25d0*Te(i+1,1)-0.25d0*Te(i-1,1)-0.25d0*Te(i-1,2))/(CurviNx(i,1)*NormalNx(i,1)+CurviNy(i,1)*NormalNy(i,1))/DistDualN(i,1) &
! 	      -(CurviSx(i,1)*TangentSx(i,1)+CurviSy(i,1)*TangentSy(i,1))*CellAreaS(i,1)*(kappae(i,1)+0.5d0*kappae(i+1,1)+0.5d0*kappae(i-1,1))*(0.5d0*Te(i+1,1)-0.5d0*Te(i-1,1))/(CurviSx(i,1)*NormalSx(i,1)+CurviSy(i,1)*NormalSy(i,1))/DistDualS(i,1) &
! 	      +(-CouplingE(i,1)+SourceE(i,1))*CellVol(i,1) &
! 	      )/Ce(i,1)/CellVol(i,1)*dt+Te(i,1)
! 	      
! 	ThNew(i,1) = (&
! 	       0.5d0*(NormalNx(i,1)**2+NormalNy(i,1)**2)*CellAreaN(i,1)*(kappah(i,1)+kappah(i  ,2))*(Th(i,2  )-Th(i,1))/(CurviNx(i,1)*NormalNx(i,1)+CurviNy(i,1)*NormalNy(i,1))/DistN(i,1) &
! 	      +0.5d0*(NormalEx(i,1)**2+NormalEy(i,1)**2)*CellAreaE(i,1)*(kappah(i,1)+kappah(i+1,1))*(Th(i+1,1)-Th(i,1))/(CurviEx(i,1)*NormalEx(i,1)+CurviEy(i,1)*NormalEy(i,1))/DistE(i,1) & 
! 	      +0.5d0*(NormalWx(i,1)**2+NormalWy(i,1)**2)*CellAreaW(i,1)*(kappah(i,1)+kappah(i-1,1))*(Th(i,1)-Th(i-1,1))/(CurviWx(i,1)*NormalWx(i,1)+CurviWy(i,1)*NormalWy(i,1))/DistW(i,1) &
! 	      -0.5d0*(CurviEx(i,1)*TangentEx(i,1)+CurviEy(i,1)*TangentEy(i,1))*CellAreaE(i,1)*(0.25d0*kappah(i+1,2)+0.75d0*kappah(i+1,1)+0.75d0*kappah(i,1)+0.25d0*kappah(i,2))*(0.25d0*Th(i+1,2)-0.25d0*Th(i+1,1)-0.25d0*Th(i,1)+0.25d0*Th(i,2))/(CurviEx(i,1)*NormalEx(i,1)+CurviEy(i,1)*NormalEy(i,1))/DistDualE(i,1) &
! 	      -0.5d0*(CurviWx(i,1)*TangentWx(i,1)+CurviWy(i,1)*TangentWy(i,1))*CellAreaW(i,1)*(0.75d0*kappah(i,1)+0.25d0*kappah(i,2)+0.75d0*kappah(i-1,1)+0.25d0*kappah(i-1,2))*(-0.25d0*Th(i,1)+0.25d0*Th(i,2)-0.25d0*Th(i-1,1)+0.25d0*Th(i-1,2))/(CurviWx(i,1)*NormalWx(i,1)+CurviWy(i,1)*NormalWy(i,1))/DistDualW(i,1) &
! 	      -0.5d0*(CurviNx(i,1)*TangentNx(i,1)+CurviNy(i,1)*TangentNy(i,1))*CellAreaN(i,1)*(0.25d0*kappah(i+1,2)+0.25d0*kappah(i+1,1)+0.5d0*kappah(i,1)+0.5d0*kappah(i,2)+0.25d0*kappah(i-1,1)+0.25d0*kappah(i-1,2))*(0.25d0*Th(i+1,2)+0.25d0*Th(i+1,1)-0.25d0*Th(i-1,1)-0.25d0*Th(i-1,2))/(CurviNx(i,1)*NormalNx(i,1)+CurviNy(i,1)*NormalNy(i,1))/DistDualN(i,1) &
! 	      -(CurviSx(i,1)*TangentSx(i,1)+CurviSy(i,1)*TangentSy(i,1))*CellAreaS(i,1)*(kappah(i,1)+0.5d0*kappah(i+1,1)+0.5d0*kappah(i-1,1))*(0.5d0*Th(i+1,1)-0.5d0*Th(i-1,1))/(CurviSx(i,1)*NormalSx(i,1)+CurviSy(i,1)*NormalSy(i,1))/DistDualS(i,1) &
! 	      +(-CouplingH(i,1)+SourceH(i,1))*CellVol(i,1) &
! 	      )/Ch(i,1)/CellVol(i,1)*dt+Th(i,1)
! 	end if
     end do
    
    
      do j=2, N-1 !West and East boundaries
      ! finite differences bad fashion
	if(DriftOn.eq.0) then
! 		NeNew(1,j)=NeNew(2,j)
! 		NhNew(1,j)=NhNew(2,j)
	end if
	
	UeNew(1,j)=UeNew(2,j)
	UhNew(1,j)=UhNew(2,j)
! 	TeNew(1,j)=TeNew(2,j)
! 	ThNew(1,j)=ThNew(2,j)
! 	TsNew(1,j)=TsNew(2,j)
  !       potential(1,j)=0d0 !(0d0, 0d0)
  !       potential(M,j)=potential0 !(potential0, 0d0)
  
	! conditions on the cone base - most important	
	if(DriftOn.eq.0) then
! 		NeNew(M,j)=NeNew(M-1,j) !Ne0 replace by dynamic bnd condition with flux equal to the one of cell(M-1,j)
! 		NhNew(M,j)=NhNew(M-1,j) !Nh0 replace by dynamic bnd condition with flux equal to the one of cell(M-1,j)
	end if
	
	!outlet condition on density and energy
! 	write(*,*) CellAreaE(M,j), DistE(M-2,j)
	

! 	NeNew(M,j) = -0.5d0*(diffusionE(M-2,j)+diffusionE(M-1,j))*(Ne(M-1,j)-Ne(M-2,j))/DistE(M-2,j)/(-0.5d0*diffusionE(M-1,j)-0.5d0*diffusionE(M,j))/DistE(M-1,j)+Ne(M-1,j)
! 	NhNew(M,j) = -0.5d0*(diffusionH(M-2,j)+diffusionH(M-1,j))*(Nh(M-1,j)-Nh(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*diffusionH(M-1,j)-0.5d0*diffusionH(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Nh(M-1,j)
	
! 	UeNew(M,j)=UeNew(M-1,j)
! 	UhNew(M,j)=UhNew(M-1,j)
	
! 	TeNew(M,j)=Tout !TeNew(M-1,j) !Tout 
! 	ThNew(M,j)=Tout !ThNew(M-1,j) !Tout 
! 	TsNew(M,j)=Tout !TsNew(M-1,j) ! Tout !cooling by diffusion from outside, TsNew(M-1,j)
	
! 	TeNew(M,j) = -0.5d0*(kappae(M-2,j)+kappae(M-1,j))*(Te(M-1,j)-Te(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappae(M-1,j)-0.5d0*kappae(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Te(M-1,j)
! 	ThNew(M,j) = -0.5d0*(kappah(M-2,j)+kappah(M-1,j))*(Th(M-1,j)-Th(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappah(M-1,j)-0.5d0*kappah(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Th(M-1,j)
! 	TsNew(M,j) = -0.5d0*(kappas(M-2,j)+kappas(M-1,j))*(Ts(M-1,j)-Ts(M-2,j))*CellAreaE(M-1,j)/DistE(M-2,j)/(-0.5d0*kappas(M-1,j)-0.5d0*kappas(M,j))/CellAreaE(M,j)*DistE(M-1,j)+Ts(M-1,j)
	
      end do
      

!       if(DriftOn.eq.1) then
	  
	  do j=2,N-1
! 	
	  ! WEST
	  
	  GradNeX(1,j) = GradNeX(2,j)
	  GradNeY(1,j) = GradNeY(2,j)
! 	  
! 	  if(y(1,j)<0d0) then
! 	    TangentNx(1,j)=-TangentNx(1,j)
! 	    TangentNy(1,j)=-TangentNy(1,j)
! 	    TangentSx(1,j)=-TangentSx(1,j)
! 	    TangentSy(1,j)=-TangentSy(1,j)
! 	    TangentEx(1,j)=-TangentEx(1,j)
! 	    TangentEy(1,j)=-TangentEy(1,j)
! 	  end if
! 	  
! ! 	  NeNew(1,j) = ((GainsE(1,j)-LossesE(1,j))*CellVol(1,j) & 
! ! 	      -( &
! !     ! 	     (0.5d0*(JeX(2,j)+JeX(1,j))*NormalEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j))*NormalEy(1,j))*CellAreaE(1,j) &
! !     ! 	    +(0.5d0*(JeX(1,j)+JeX(1,j+1))*NormalNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1))*NormalNy(1,j))*CellAreaN(1,j) &
! !     ! 	    +(0.5d0*(JeX(1,j)+JeX(1,j-1))*NormalSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*NormalSy(1,j))*CellAreaS(1,j) &
! ! 		((0.5d0 * (JeX(1,j)+JeX(1,j+1)) * TangentNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1)) * TangentNy(1,j) ) * TangentNx(1,j) & 
! ! 		 +(0.5d0 * (JeX(1,j)+JeX(1,j+1)) *  NormalNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1)) * NormalNy(1,j) ) *  NormalNx(1,j) )* & 
! ! 		 NormalNx(1,j)*CellAreaN(1,j) &
! ! 		+((0.5d0*(JeX(1,j)+JeX(1,j+1))*TangentNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1))*TangentNy(1,j))*TangentNy(1,j) & 
! ! 		+(0.5d0*(JeX(1,j)+JeX(1,j+1))*NormalNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1))*NormalNy(1,j))*TangentNy(1,j)) & 
! ! 		 *NormalNy(1,j)*CellAreaN(1,j) &
! ! 		+((0.5d0*(JeX(1,j)+JeX(1,j-1))*TangentSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*TangentSy(1,j))*TangentSx(1,j) &
! ! 		+ (0.5d0*(JeX(1,j)+JeX(1,j-1))*NormalSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*NormalSy(1,j))*NormalSx(1,j)) & 
! ! 		 *NormalSx(1,j)*CellAreaS(1,j) &
! ! 		+((0.5d0*(JeX(1,j)+JeX(1,j-1))*TangentSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*TangentSy(1,j))*TangentSy(1,j) & 
! ! 		+ (0.5d0*(JeX(1,j)+JeX(1,j-1))*NormalSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*NormalSy(1,j))*NormalSy(1,j)) & 
! ! 		*NormalSy(1,j)*CellAreaS(1,j) & 
! ! 		+((0.5d0*(JeX(2,j)+JeX(1,j)  )*TangentEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j)  )*TangentEy(1,j))*TangentEx(1,j) &
! ! 		+ (0.5d0*(JeX(2,j)+JeX(1,j)  )*NormalEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j)  )*NormalEy(1,j))*NormalEx(1,j)) &
! ! 		*NormalEx(1,j)*CellAreaE(1,j) &
! ! 		+((0.5d0*(JeX(2,j)+JeX(1,j)  )*TangentEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j)  )*TangentEy(1,j))*TangentEy(1,j) & 
! ! 		+ (0.5d0*(JeX(2,j)+JeX(1,j)  )*NormalEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j)  )*NormalEy(1,j))*NormalEy(1,j)) & 
! ! 		*NormalEy(1,j)*CellAreaE(1,j) & 
! ! 		) &
! ! 	      )*dt/CellVol(1,j)+Ne(1,j)
! 	      
! 	  NeNew(1,j) = ((GainsE(1,j)-LossesE(1,j))*CellVol(1,j) & 
! 	      -( &
! 		 (0.5d0*(JeX(1,j)+JeX(1,j+1))*TangentNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1))*TangentNy(1,j))*TangentNx(1,j)*NormalNx(1,j)*CellAreaN(1,j) &
! 		+(0.5d0*(JeX(1,j)+JeX(1,j+1))*TangentNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1))*TangentNy(1,j))*TangentNy(1,j)*NormalNy(1,j)*CellAreaN(1,j) &
! 		+(0.5d0*(JeX(1,j)+JeX(1,j-1))*TangentSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*TangentSy(1,j))*TangentSx(1,j)*NormalSx(1,j)*CellAreaS(1,j) &
! 		+(0.5d0*(JeX(1,j)+JeX(1,j-1))*TangentSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*TangentSy(1,j))*TangentSy(1,j)*NormalSy(1,j)*CellAreaS(1,j) & 
! 		+(0.5d0*(JeX(2,j)+JeX(1,j)  )*TangentEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j)  )*TangentEy(1,j))*TangentEx(1,j)*NormalEx(1,j)*CellAreaE(1,j) &
! 		+(0.5d0*(JeX(2,j)+JeX(1,j)  )*TangentEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j)  )*TangentEy(1,j))*TangentEy(1,j)*NormalEy(1,j)*CellAreaE(1,j) &
! 		+(0.5d0*(JeX(1,j)+JeX(1,j+1))*NormalNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1))*NormalNy(1,j))*NormalNx(1,j)*NormalNx(1,j)*CellAreaN(1,j) &
! 		+(0.5d0*(JeX(1,j)+JeX(1,j+1))*NormalNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1))*NormalNy(1,j))*NormalNy(1,j)*NormalNy(1,j)*CellAreaN(1,j) &
! 		+(0.5d0*(JeX(1,j)+JeX(1,j-1))*NormalSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*NormalSy(1,j))*NormalSx(1,j)*NormalSx(1,j)*CellAreaS(1,j) &
! 		+(0.5d0*(JeX(1,j)+JeX(1,j-1))*NormalSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*NormalSy(1,j))*NormalSy(1,j)*NormalSy(1,j)*CellAreaS(1,j) & 
! 		+(0.5d0*(JeX(2,j)+JeX(1,j)  )*NormalEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j)  )*NormalEy(1,j))*NormalEx(1,j)*NormalEx(1,j)*CellAreaE(1,j) &
! 		+(0.5d0*(JeX(2,j)+JeX(1,j)  )*NormalEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j)  )*NormalEy(1,j))*NormalEy(1,j)*NormalEy(1,j)*CellAreaE(1,j) &
! ! 		(0.5d0*(JeX(2,j)+JeX(1,j))*NormalEx(1,j)+0.5d0*(JeY(2,j)+JeY(1,j))*NormalEy(1,j))*CellAreaE(1,j) &
! ! 		+(0.5d0*(JeX(1,j)+JeX(1,j+1))*NormalNx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j+1))*NormalNy(1,j))*CellAreaN(1,j) &
! ! 		+(0.5d0*(JeX(1,j)+JeX(1,j-1))*NormalSx(1,j)+0.5d0*(JeY(1,j)+JeY(1,j-1))*NormalSy(1,j))*CellAreaS(1,j) &
! 		) &
! 	      )*dt/CellVol(1,j) & 
! 	      +( 0.5d0*(NormalEx(1,j)**2+NormalEy(1,j)**2)*CellAreaE(1,j)*(diffusionE(2,j)+diffusionE(1,j))*(Ne(2,j)-Ne(1,j))/(CurviEx(1,j)*NormalEx(1,j)+CurviEy(1,j)*NormalEy(1,j))/DistE(1,j) & 
! 		-0.5d0*(CurviEx(1,j)*TangentEx(1,j)+CurviEy(1,j)*TangentEy(1,j))*CellAreaE(1,j)*(0.25d0*diffusionE(2,j+1)+0.5d0*diffusionE(2,j)+0.5d0*diffusionE(1,j)+0.25d0*diffusionE(1,j+1)+0.25d0*diffusionE(2,j-1)+0.25d0*diffusionE(1,j-1))*(0.25d0*Ne(2,j+1)+0.25d0*Ne(1,j+1)-0.25d0*Ne(2,j-1)-0.25d0*Ne(1,j-1))/(CurviEx(1,j)*NormalEx(1,j)+CurviEy(1,j)*NormalEy(1,j))/DistDualE(1,j) & 
! 		-0.5d0*(CurviWx(1,j)*TangentWx(1,j)+CurviWy(1,j)*TangentWy(1,j))*CellAreaW(1,j)*(diffusionE(1,j)+0.5d0*diffusionE(1,j+1)+0.5d0*diffusionE(1,j-1))*(0.5d0*Ne(1,j+1)-0.5d0*Ne(1,j-1))/(CurviWx(1,j)*NormalWx(1,j)+CurviWy(1,j)*NormalWy(1,j))/DistDualW(1,j) & 
! 		+0.5d0*(NormalNx(1,j)**2+NormalNy(1,j)**2)*CellAreaN(1,j)*(diffusionE(1,j)+diffusionE(1,j+1))*(Ne(1,j+1)-Ne(1,j))/(CurviNx(1,j)*NormalNx(1,j)+CurviNy(1,j)*NormalNy(1,j))/DistN(1,j) & 
! 		-0.5d0*(CurviNx(1,j)*TangentNx(1,j)+CurviNy(1,j)*TangentNy(1,j))*CellAreaN(1,j)*(0.25d0*diffusionE(2,j+1)+0.25d0*diffusionE(2,j)+0.75d0*diffusionE(1,j)+0.75d0*diffusionE(1,j+1))*(0.25d0*Ne(2,j+1)+0.25d0*Ne(2,j)-0.25d0*Ne(1,j)-0.25d0*Ne(1,j+1))/(CurviNx(1,j)*NormalNx(1,j)+CurviNy(1,j)*NormalNy(1,j))/DistDualN(1,j) & 
! 		+0.5d0*(NormalSx(1,j)**2+NormalSy(1,j)**2)*CellAreaS(1,j)*(diffusionE(1,j-1)+diffusionE(1,j))*(Ne(1,j)-Ne(1,j-1))/(CurviSx(1,j)*NormalSx(1,j)+CurviSy(1,j)*NormalSy(1,j))/DistS(1,j) & 
! 		-0.5d0*(CurviSx(1,j)*TangentSx(1,j)+CurviSy(1,j)*TangentSy(1,j))*CellAreaS(1,j)*(0.25d0*diffusionE(2,j)+0.25d0*diffusionE(2,j-1)+0.75d0*diffusionE(1,j)+0.75d0*diffusionE(1,j-1))*(0.25d0*Ne(2,j)+0.25d0*Ne(2,j-1)-0.25d0*Ne(1,j)-0.25d0*Ne(1,j-1))/(CurviSx(1,j)*NormalSx(1,j)+CurviSy(1,j)*NormalSy(1,j))/DistDualS(1,j) & 
! 		)*dt/CellVol(1,j) & 
! 		+Ne(1,j)    
! 	      
! 	  NhNew(1,j) = ((GainsH(1,j)-LossesH(1,j))*CellVol(1,j) & 
! 	      -( &
! 		 (0.5d0*(JhX(1,j)+JhX(1,j+1))*TangentNx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j+1))*TangentNy(1,j))*TangentNx(1,j)*NormalNx(1,j)*CellAreaN(1,j) &
! 		+(0.5d0*(JhX(1,j)+JhX(1,j+1))*TangentNx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j+1))*TangentNy(1,j))*TangentNy(1,j)*NormalNy(1,j)*CellAreaN(1,j) &
! 		+(0.5d0*(JhX(1,j)+JhX(1,j-1))*TangentSx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j-1))*TangentSy(1,j))*TangentSx(1,j)*NormalSx(1,j)*CellAreaS(1,j) &
! 		+(0.5d0*(JhX(1,j)+JhX(1,j-1))*TangentSx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j-1))*TangentSy(1,j))*TangentSy(1,j)*NormalSy(1,j)*CellAreaS(1,j) & 
! 		+(0.5d0*(JhX(2,j)+JhX(1,j)  )*TangentEx(1,j)+0.5d0*(JhY(2,j)+JhY(1,j)  )*TangentEy(1,j))*TangentEx(1,j)*NormalEx(1,j)*CellAreaE(1,j) &
! 		+(0.5d0*(JhX(2,j)+JhX(1,j)  )*TangentEx(1,j)+0.5d0*(JhY(2,j)+JhY(1,j)  )*TangentEy(1,j))*TangentEy(1,j)*NormalEy(1,j)*CellAreaE(1,j) &
! 		+(0.5d0*(JhX(1,j)+JhX(1,j+1))*NormalNx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j+1))*NormalNy(1,j))*NormalNx(1,j)*NormalNx(1,j)*CellAreaN(1,j) &
! 		+(0.5d0*(JhX(1,j)+JhX(1,j+1))*NormalNx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j+1))*NormalNy(1,j))*NormalNy(1,j)*NormalNy(1,j)*CellAreaN(1,j) &
! 		+(0.5d0*(JhX(1,j)+JhX(1,j-1))*NormalSx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j-1))*NormalSy(1,j))*NormalSx(1,j)*NormalSx(1,j)*CellAreaS(1,j) &
! 		+(0.5d0*(JhX(1,j)+JhX(1,j-1))*NormalSx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j-1))*NormalSy(1,j))*NormalSy(1,j)*NormalSy(1,j)*CellAreaS(1,j) & 
! 		+(0.5d0*(JhX(2,j)+JhX(1,j)  )*NormalEx(1,j)+0.5d0*(JhY(2,j)+JhY(1,j)  )*NormalEy(1,j))*NormalEx(1,j)*NormalEx(1,j)*CellAreaE(1,j) &
! 		+(0.5d0*(JhX(2,j)+JhX(1,j)  )*NormalEx(1,j)+0.5d0*(JhY(2,j)+JhY(1,j)  )*NormalEy(1,j))*NormalEy(1,j)*NormalEy(1,j)*CellAreaE(1,j) &
! ! 		(0.5d0*(JhX(2,j)+JhX(1,j))*NormalEx(1,j)+0.5d0*(JhY(2,j)+JhY(1,j))*NormalEy(1,j))*CellAreaE(1,j) &
! ! 		+(0.5d0*(JhX(1,j)+JhX(1,j+1))*NormalNx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j+1))*NormalNy(1,j))*CellAreaN(1,j) &
! ! 		+(0.5d0*(JhX(1,j)+JhX(1,j-1))*NormalSx(1,j)+0.5d0*(JhY(1,j)+JhY(1,j-1))*NormalSy(1,j))*CellAreaS(1,j) &
! 		) &
! 	      )*dt/CellVol(1,j) & 
! 	      +( 0.5d0*(NormalEx(1,j)**2+NormalEy(1,j)**2)*CellAreaE(1,j)*(diffusionH(2,j)+diffusionH(1,j))*(Nh(2,j)-Nh(1,j))/(CurviEx(1,j)*NormalEx(1,j)+CurviEy(1,j)*NormalEy(1,j))/DistE(1,j) & 
! 		-0.5d0*(CurviEx(1,j)*TangentEx(1,j)+CurviEy(1,j)*TangentEy(1,j))*CellAreaE(1,j)*(0.25d0*diffusionH(2,j+1)+0.5d0*diffusionH(2,j)+0.5d0*diffusionH(1,j)+0.25d0*diffusionH(1,j+1)+0.25d0*diffusionH(2,j-1)+0.25d0*diffusionH(1,j-1))*(0.25d0*Nh(2,j+1)+0.25d0*Nh(1,j+1)-0.25d0*Nh(2,j-1)-0.25d0*Nh(1,j-1))/(CurviEx(1,j)*NormalEx(1,j)+CurviEy(1,j)*NormalEy(1,j))/DistDualE(1,j) & 
! 		-0.5d0*(CurviWx(1,j)*TangentWx(1,j)+CurviWy(1,j)*TangentWy(1,j))*CellAreaW(1,j)*(diffusionH(1,j)+0.5d0*diffusionH(1,j+1)+0.5d0*diffusionH(1,j-1))*(0.5d0*Nh(1,j+1)-0.5d0*Nh(1,j-1))/(CurviWx(1,j)*NormalWx(1,j)+CurviWy(1,j)*NormalWy(1,j))/DistDualW(1,j) & 
! 		+0.5d0*(NormalNx(1,j)**2+NormalNy(1,j)**2)*CellAreaN(1,j)*(diffusionH(1,j)+diffusionH(1,j+1))*(Nh(1,j+1)-Nh(1,j))/(CurviNx(1,j)*NormalNx(1,j)+CurviNy(1,j)*NormalNy(1,j))/DistN(1,j) & 
! 		-0.5d0*(CurviNx(1,j)*TangentNx(1,j)+CurviNy(1,j)*TangentNy(1,j))*CellAreaN(1,j)*(0.25d0*diffusionH(2,j+1)+0.25d0*diffusionH(2,j)+0.75d0*diffusionH(1,j)+0.75d0*diffusionH(1,j+1))*(0.25d0*Nh(2,j+1)+0.25d0*Nh(2,j)-0.25d0*Nh(1,j)-0.25d0*Nh(1,j+1))/(CurviNx(1,j)*NormalNx(1,j)+CurviNy(1,j)*NormalNy(1,j))/DistDualN(1,j) & 
! 		+0.5d0*(NormalSx(1,j)**2+NormalSy(1,j)**2)*CellAreaS(1,j)*(diffusionH(1,j-1)+diffusionH(1,j))*(Nh(1,j)-Nh(1,j-1))/(CurviSx(1,j)*NormalSx(1,j)+CurviSy(1,j)*NormalSy(1,j))/DistS(1,j) & 
! 		-0.5d0*(CurviSx(1,j)*TangentSx(1,j)+CurviSy(1,j)*TangentSy(1,j))*CellAreaS(1,j)*(0.25d0*diffusionH(2,j)+0.25d0*diffusionH(2,j-1)+0.75d0*diffusionH(1,j)+0.75d0*diffusionH(1,j-1))*(0.25d0*Nh(2,j)+0.25d0*Nh(2,j-1)-0.25d0*Nh(1,j)-0.25d0*Nh(1,j-1))/(CurviSx(1,j)*NormalSx(1,j)+CurviSy(1,j)*NormalSy(1,j))/DistDualS(1,j) & 
! 		)*dt/CellVol(1,j) & 
! 	      
! 	      +Nh(1,j)
! 	  if(TeOff.ne.1) then
! 	  TeNew(1,j) = ( 0.5d0*(NormalEx(1,j)**2+NormalEy(1,j)**2)*CellAreaE(1,j)*(kappae(2,j)+kappae(1,j))*(Te(2,j)-Te(1,j))/(CurviEx(1,j)*NormalEx(1,j)+CurviEy(1,j)*NormalEy(1,j))/DistE(1,j) & 
! 			-0.5d0*(CurviEx(1,j)*TangentEx(1,j)+CurviEy(1,j)*TangentEy(1,j))*CellAreaE(1,j)*(0.25d0*kappae(2,j+1)+0.5d0*kappae(2,j)+0.5d0*kappae(1,j)+0.25d0*kappae(1,j+1)+0.25d0*kappae(2,j-1)+0.25d0*kappae(1,j-1))*(0.25d0*Te(2,j+1)+0.25d0*Te(1,j+1)-0.25d0*Te(2,j-1)-0.25d0*Te(1,j-1))/(CurviEx(1,j)*NormalEx(1,j)+CurviEy(1,j)*NormalEy(1,j))/DistDualE(1,j) & 
! 			-0.5d0*(CurviWx(1,j)*TangentWx(1,j)+CurviWy(1,j)*TangentWy(1,j))*CellAreaW(1,j)*(kappae(1,j)+0.5d0*kappae(1,j+1)+0.5d0*kappae(1,j-1))*(0.5d0*Te(1,j+1)-0.5d0*Te(1,j-1))/(CurviWx(1,j)*NormalWx(1,j)+CurviWy(1,j)*NormalWy(1,j))/DistDualW(1,j) & 
! 			+0.5d0*(NormalNx(1,j)**2+NormalNy(1,j)**2)*CellAreaN(1,j)*(kappae(1,j)+kappae(1,j+1))*(Te(1,j+1)-Te(1,j))/(CurviNx(1,j)*NormalNx(1,j)+CurviNy(1,j)*NormalNy(1,j))/DistN(1,j) & 
! 			-0.5d0*(CurviNx(1,j)*TangentNx(1,j)+CurviNy(1,j)*TangentNy(1,j))*CellAreaN(1,j)*(0.25d0*kappae(2,j+1)+0.25d0*kappae(2,j)+0.75d0*kappae(1,j)+0.75d0*kappae(1,j+1))*(0.25d0*Te(2,j+1)+0.25d0*Te(2,j)-0.25d0*Te(1,j)-0.25d0*Te(1,j+1))/(CurviNx(1,j)*NormalNx(1,j)+CurviNy(1,j)*NormalNy(1,j))/DistDualN(1,j) & 
! 			+0.5d0*(NormalSx(1,j)**2+NormalSy(1,j)**2)*CellAreaS(1,j)*(kappae(1,j-1)+kappae(1,j))*(Te(1,j)-Te(1,j-1))/(CurviSx(1,j)*NormalSx(1,j)+CurviSy(1,j)*NormalSy(1,j))/DistS(1,j) & 
! 			-0.5d0*(CurviSx(1,j)*TangentSx(1,j)+CurviSy(1,j)*TangentSy(1,j))*CellAreaS(1,j)*(0.25d0*kappae(2,j)+0.25d0*kappae(2,j-1)+0.75d0*kappae(1,j)+0.75d0*kappae(1,j-1))*(0.25d0*Te(2,j)+0.25d0*Te(2,j-1)-0.25d0*Te(1,j)-0.25d0*Te(1,j-1))/(CurviSx(1,j)*NormalSx(1,j)+CurviSy(1,j)*NormalSy(1,j))/DistDualS(1,j) & 
! 			+(-CouplingE(1,j)+SourceE(1,j))*CellVol(1,j) & 
! 			)/Ce(1,j)*dt/CellVol(1,j) & 
! 			+Te(1,j)
! 			
! 	  ThNew(1,j) = ( 0.5d0*(NormalEx(1,j)**2+NormalEy(1,j)**2)*CellAreaE(1,j)*(kappah(2,j)+kappah(1,j))*(Th(2,j)-Th(1,j))/(CurviEx(1,j)*NormalEx(1,j)+CurviEy(1,j)*NormalEy(1,j))/DistE(1,j) & 
! 			-0.5d0*(CurviEx(1,j)*TangentEx(1,j)+CurviEy(1,j)*TangentEy(1,j))*CellAreaE(1,j)*(0.25d0*kappah(2,j+1)+0.5d0*kappah(2,j)+0.5d0*kappah(1,j)+0.25d0*kappah(1,j+1)+0.25d0*kappah(2,j-1)+0.25d0*kappah(1,j-1))*(0.25d0*Th(2,j+1)+0.25d0*Th(1,j+1)-0.25d0*Th(2,j-1)-0.25d0*Th(1,j-1))/(CurviEx(1,j)*NormalEx(1,j)+CurviEy(1,j)*NormalEy(1,j))/DistDualE(1,j) & 
! 			-0.5d0*(CurviWx(1,j)*TangentWx(1,j)+CurviWy(1,j)*TangentWy(1,j))*CellAreaW(1,j)*(kappah(1,j)+0.5d0*kappah(1,j+1)+0.5d0*kappah(1,j-1))*(0.5d0*Th(1,j+1)-0.5d0*Th(1,j-1))/(CurviWx(1,j)*NormalWx(1,j)+CurviWy(1,j)*NormalWy(1,j))/DistDualW(1,j) & 
! 			+0.5d0*(NormalNx(1,j)**2+NormalNy(1,j)**2)*CellAreaN(1,j)*(kappah(1,j)+kappah(1,j+1))*(Th(1,j+1)-Th(1,j))/(CurviNx(1,j)*NormalNx(1,j)+CurviNy(1,j)*NormalNy(1,j))/DistN(1,j) & 
! 			-0.5d0*(CurviNx(1,j)*TangentNx(1,j)+CurviNy(1,j)*TangentNy(1,j))*CellAreaN(1,j)*(0.25d0*kappah(2,j+1)+0.25d0*kappah(2,j)+0.75d0*kappah(1,j)+0.75d0*kappah(1,j+1))*(0.25d0*Th(2,j+1)+0.25d0*Th(2,j)-0.25d0*Th(1,j)-0.25d0*Th(1,j+1))/(CurviNx(1,j)*NormalNx(1,j)+CurviNy(1,j)*NormalNy(1,j))/DistDualN(1,j) & 
! 			+0.5d0*(NormalSx(1,j)**2+NormalSy(1,j)**2)*CellAreaS(1,j)*(kappah(1,j-1)+kappah(1,j))*(Th(1,j)-Th(1,j-1))/(CurviSx(1,j)*NormalSx(1,j)+CurviSy(1,j)*NormalSy(1,j))/DistS(1,j) & 
! 			-0.5d0*(CurviSx(1,j)*TangentSx(1,j)+CurviSy(1,j)*TangentSy(1,j))*CellAreaS(1,j)*(0.25d0*kappah(2,j)+0.25d0*kappah(2,j-1)+0.75d0*kappah(1,j)+0.75d0*kappah(1,j-1))*(0.25d0*Th(2,j)+0.25d0*Th(2,j-1)-0.25d0*Th(1,j)-0.25d0*Th(1,j-1))/(CurviSx(1,j)*NormalSx(1,j)+CurviSy(1,j)*NormalSy(1,j))/DistDualS(1,j) & 
! 			+(-CouplingH(1,j)+SourceH(1,j))*CellVol(1,j) & 
! 			)/Ch(1,j)*dt/CellVol(1,j) & 
! 			+Th(1,j)
! 	  end if
	 
	 ! EAST
	  
	  GradNeX(M,j) = GradNeX(M-1,j)
	  GradNeY(M,j) = GradNeY(M-1,j)
! 
! ! 	  NeNew(M,j) = ((GainsE(M,j)-LossesE(M,j))*CellVol(M,j) & 
! ! ! 	      -( &
! ! ! 		 (0.5d0*(JeX(M-1,j)+JeX(M,j))*TangentWx(M,j)+0.5d0*(JeY(M-1,j)+JeY(M,j))*TangentWy(M,j))*TangentWx(M,j)*NormalWx(M,j)*CellAreaW(M,j) &
! ! ! 		+(0.5d0*(JeX(M-1,j)+JeX(M,j))*TangentWx(M,j)+0.5d0*(JeY(M-1,j)+JeY(M,j))*TangentWy(M,j))*TangentWy(M,j)*NormalWy(M,j)*CellAreaW(M,j) &
! ! ! 		+(0.5d0*(JeX(M,j)+JeX(M,j+1))*TangentNx(M,j)+0.5d0*(JeY(M,j)+JeY(M,j+1))*TangentNy(M,j))*TangentNx(M,j)*NormalNx(M,j)*CellAreaN(M,j) &
! ! ! 		+(0.5d0*(JeX(M,j)+JeX(M,j+1))*TangentNx(M,j)+0.5d0*(JeY(M,j)+JeY(M,j+1))*TangentNy(M,j))*TangentNy(M,j)*NormalNy(M,j)*CellAreaN(M,j) &
! ! ! 		+(0.5d0*(JeX(M,j)+JeX(M,j-1))*TangentSx(M,j)+0.5d0*(JeY(M,j)+JeY(M,j-1))*TangentSy(M,j))*TangentSx(M,j)*NormalSx(M,j)*CellAreaS(M,j) &
! ! ! 		+(0.5d0*(JeX(M,j)+JeX(M,j-1))*TangentSx(M,j)+0.5d0*(JeY(M,j)+JeY(M,j-1))*TangentSy(M,j))*TangentSy(M,j)*NormalSy(M,j)*CellAreaS(M,j) &
! ! ! 		+(0.5d0*(JeX(M-1,j)+JeX(M,j))*NormalWx(M,j)+0.5d0*(JeY(M-1,j)+JeY(M,j))*NormalWy(M,j))*NormalWx(M,j)*NormalWx(M,j)*CellAreaW(M,j) &
! ! ! 		+(0.5d0*(JeX(M-1,j)+JeX(M,j))*NormalWx(M,j)+0.5d0*(JeY(M-1,j)+JeY(M,j))*NormalWy(M,j))*NormalWy(M,j)*NormalWy(M,j)*CellAreaW(M,j) &
! ! ! 		+(0.5d0*(JeX(M,j)+JeX(M,j+1))*NormalNx(M,j)+0.5d0*(JeY(M,j)+JeY(M,j+1))*NormalNy(M,j))*NormalNx(M,j)*NormalNx(M,j)*CellAreaN(M,j) &
! ! ! 		+(0.5d0*(JeX(M,j)+JeX(M,j+1))*NormalNx(M,j)+0.5d0*(JeY(M,j)+JeY(M,j+1))*NormalNy(M,j))*NormalNy(M,j)*NormalNy(M,j)*CellAreaN(M,j) &
! ! ! 		+(0.5d0*(JeX(M,j)+JeX(M,j-1))*NormalSx(M,j)+0.5d0*(JeY(M,j)+JeY(M,j-1))*NormalSy(M,j))*NormalSx(M,j)*NormalSx(M,j)*CellAreaS(M,j) &
! ! ! 		+(0.5d0*(JeX(M,j)+JeX(M,j-1))*NormalSx(M,j)+0.5d0*(JeY(M,j)+JeY(M,j-1))*NormalSy(M,j))*NormalSy(M,j)*NormalSy(M,j)*CellAreaS(M,j) &
! ! ! 		) &
! ! 	      )*dt/CellVol(M,j) & 
! ! 	      +( & 
! ! 	      -0.5d0*(CurviEx(M,j)*TangentEx(M,j)+CurviEy(M,j)*TangentEy(M,j))*CellAreaE(M,j)*(0.5d0*diffusionE(M,j+1)+diffusionE(M,j)+0.5d0*diffusionE(M,j-1))*(0.5d0*Ne(M,j+1)-0.5d0*Ne(M,j-1))/(CurviEx(M,j)*NormalEx(M,j)+CurviEy(M,j)*NormalEy(M,j))/DistDualE(M,j) & 
! ! 	      +0.5d0*(NormalWx(M,j)**2+NormalWy(M,j)**2)*CellAreaW(M,j)*(diffusionE(M,j)+diffusionE(M-1,j))*(Ne(M,j)-Ne(M-1,j))/(CurviWx(M,j)*NormalWx(M,j)+CurviWy(M,j)*NormalWy(M,j))/DistW(M,j) & 
! ! 	      -0.5d0*(CurviWx(M,j)*TangentWx(M,j)+CurviWy(M,j)*TangentWy(M,j))*CellAreaW(M,j)*(0.5d0*diffusionE(M,j)+0.25d0*diffusionE(M,j+1)+0.5d0*diffusionE(M-1,j)+0.25d0*diffusionE(M-1,j+1)+0.25d0*diffusionE(M-1,j-1)+0.25d0*diffusionE(M,j-1))*(0.25d0*Ne(M,j+1)+0.25d0*Ne(M-1,j+1)-0.25d0*Ne(M-1,j-1)-0.25d0*Ne(M,j-1))/(CurviWx(M,j)*NormalWx(M,j)+CurviWy(M,j)*NormalWy(M,j))/DistDualW(M,j) & 
! ! 	      +0.5d0*(NormalNx(M,j)**2+NormalNy(M,j)**2)*CellAreaN(M,j)*(diffusionE(M,j+1)+diffusionE(M,j))*(Ne(M,j+1)-Ne(M,j))/(CurviNx(M,j)*NormalNx(M,j)+CurviNy(M,j)*NormalNy(M,j))/DistN(M,j) & 
! ! 	      -0.5d0*(CurviNx(M,j)*TangentNx(M,j)+CurviNy(M,j)*TangentNy(M,j))*CellAreaN(M,j)*(0.75d0*diffusionE(M,j+1)+0.75d0*diffusionE(M,j)+0.25d0*diffusionE(M-1,j)+0.25d0*diffusionE(M-1,j+1))*(0.25d0*Ne(M,j+1)+0.25d0*Ne(M,j)-0.25d0*Ne(M-1,j)-0.25d0*Ne(M-1,j+1))/(CurviNx(M,j)*NormalNx(M,j)+CurviNy(M,j)*NormalNy(M,j))/DistDualN(M,j) & 
! ! 	      +0.5d0*(NormalSx(M,j)**2+NormalSy(M,j)**2)*CellAreaS(M,j)*(diffusionE(M,j)+diffusionE(M,j-1))*(Ne(M,j)-Ne(M,j-1))/(CurviSx(M,j)*NormalSx(M,j)+CurviSy(M,j)*NormalSy(M,j))/DistS(M,j) & 
! ! 	      -0.5d0*(CurviSx(M,j)*TangentSx(M,j)+CurviSy(M,j)*TangentSy(M,j))*CellAreaS(M,j)*(0.25d0*diffusionE(M,j)+0.25d0*diffusionE(M,j-1)-0.25d0*diffusionE(M-1,j-1)-0.25d0*diffusionE(M-1,j))*(0.25d0*Ne(M,j)+0.25d0*Ne(M,j-1)-0.25d0*Ne(M-1,j-1)-0.25d0*Ne(M-1,j))/(CurviSx(M,j)*NormalSx(M,j)+CurviSy(M,j)*NormalSy(M,j))/DistDualS(M,j) & 
! ! 	      )/CellVol(M,j)*dt & 
! ! 	      +Ne(M,j)
! ! 	      
! ! 	  NhNew(M,j) = ((GainsH(M,j)-LossesH(M,j))*CellVol(M,j) & 
! ! 	  -( &
! ! 	      (0.5d0*(JhX(M-1,j)+JhX(M,j))*TangentWx(M,j)+0.5d0*(JhY(M-1,j)+JhY(M,j))*TangentWy(M,j))*TangentWx(M,j)*NormalWx(M,j)*CellAreaW(M,j) &
! ! 	    +(0.5d0*(JhX(M-1,j)+JhX(M,j))*TangentWx(M,j)+0.5d0*(JhY(M-1,j)+JhY(M,j))*TangentWy(M,j))*TangentWy(M,j)*NormalWy(M,j)*CellAreaW(M,j) &
! ! 	    +(0.5d0*(JhX(M,j)+JhX(M,j+1))*TangentNx(M,j)+0.5d0*(JhY(M,j)+JhY(M,j+1))*TangentNy(M,j))*TangentNx(M,j)*NormalNx(M,j)*CellAreaN(M,j) &
! ! 	    +(0.5d0*(JhX(M,j)+JhX(M,j+1))*TangentNx(M,j)+0.5d0*(JhY(M,j)+JhY(M,j+1))*TangentNy(M,j))*TangentNy(M,j)*NormalNy(M,j)*CellAreaN(M,j) &
! ! 	    +(0.5d0*(JhX(M,j)+JhX(M,j-1))*TangentSx(M,j)+0.5d0*(JhY(M,j)+JhY(M,j-1))*TangentSy(M,j))*TangentSx(M,j)*NormalSx(M,j)*CellAreaS(M,j) &
! ! 	    +(0.5d0*(JhX(M,j)+JhX(M,j-1))*TangentSx(M,j)+0.5d0*(JhY(M,j)+JhY(M,j-1))*TangentSy(M,j))*TangentSy(M,j)*NormalSy(M,j)*CellAreaS(M,j) &
! ! 	    +(0.5d0*(JhX(M-1,j)+JhX(M,j))*NormalWx(M,j)+0.5d0*(JhY(M-1,j)+JhY(M,j))*NormalWy(M,j))*NormalWx(M,j)*NormalWx(M,j)*CellAreaW(M,j) &
! ! 	    +(0.5d0*(JhX(M-1,j)+JhX(M,j))*NormalWx(M,j)+0.5d0*(JhY(M-1,j)+JhY(M,j))*NormalWy(M,j))*NormalWy(M,j)*NormalWy(M,j)*CellAreaW(M,j) &
! ! 	    +(0.5d0*(JhX(M,j)+JhX(M,j+1))*NormalNx(M,j)+0.5d0*(JhY(M,j)+JhY(M,j+1))*NormalNy(M,j))*NormalNx(M,j)*NormalNx(M,j)*CellAreaN(M,j) &
! ! 	    +(0.5d0*(JhX(M,j)+JhX(M,j+1))*NormalNx(M,j)+0.5d0*(JhY(M,j)+JhY(M,j+1))*NormalNy(M,j))*NormalNy(M,j)*NormalNy(M,j)*CellAreaN(M,j) &
! ! 	    +(0.5d0*(JhX(M,j)+JhX(M,j-1))*NormalSx(M,j)+0.5d0*(JhY(M,j)+JhY(M,j-1))*NormalSy(M,j))*NormalSx(M,j)*NormalSx(M,j)*CellAreaS(M,j) &
! ! 	    +(0.5d0*(JhX(M,j)+JhX(M,j-1))*NormalSx(M,j)+0.5d0*(JhY(M,j)+JhY(M,j-1))*NormalSy(M,j))*NormalSy(M,j)*NormalSy(M,j)*CellAreaS(M,j) &
! ! 	    ) &
! ! 	  )*dt/CellVol(M,j) & 
! ! 	  +( & 
! ! ! 	  -0.5d0*(CurviEx(M,j)*TangentEx(M,j)+CurviEy(M,j)*TangentEy(M,j))*CellAreaE(M,j)*(0.5d0*diffusionH(M,j+1)+diffusionH(M,j)+0.5d0*diffusionH(M,j-1))*(0.5d0*Nh(M,j+1)-0.5d0*Nh(M,j-1))/(CurviEx(M,j)*NormalEx(M,j)+CurviEy(M,j)*NormalEy(M,j))/DistDualE(M,j) & 
! ! 	  +0.5d0*(NormalWx(M,j)**2+NormalWy(M,j)**2)*CellAreaW(M,j)*(diffusionH(M,j)+diffusionH(M-1,j))*(Nh(M,j)-Nh(M-1,j))/(CurviWx(M,j)*NormalWx(M,j)+CurviWy(M,j)*NormalWy(M,j))/DistW(M,j) & 
! ! 	  -0.5d0*(CurviWx(M,j)*TangentWx(M,j)+CurviWy(M,j)*TangentWy(M,j))*CellAreaW(M,j)*(0.5d0*diffusionH(M,j)+0.25d0*diffusionH(M,j+1)+0.5d0*diffusionH(M-1,j)+0.25d0*diffusionH(M-1,j+1)+0.25d0*diffusionH(M-1,j-1)+0.25d0*diffusionH(M,j-1))*(0.25d0*Nh(M,j+1)+0.25d0*Nh(M-1,j+1)-0.25d0*Nh(M-1,j-1)-0.25d0*Nh(M,j-1))/(CurviWx(M,j)*NormalWx(M,j)+CurviWy(M,j)*NormalWy(M,j))/DistDualW(M,j) & 
! ! 	  +0.5d0*(NormalNx(M,j)**2+NormalNy(M,j)**2)*CellAreaN(M,j)*(diffusionH(M,j+1)+diffusionH(M,j))*(Nh(M,j+1)-Nh(M,j))/(CurviNx(M,j)*NormalNx(M,j)+CurviNy(M,j)*NormalNy(M,j))/DistN(M,j) & 
! ! 	  -0.5d0*(CurviNx(M,j)*TangentNx(M,j)+CurviNy(M,j)*TangentNy(M,j))*CellAreaN(M,j)*(0.75d0*diffusionH(M,j+1)+0.75d0*diffusionH(M,j)+0.25d0*diffusionH(M-1,j)+0.25d0*diffusionH(M-1,j+1))*(0.25d0*Nh(M,j+1)+0.25d0*Nh(M,j)-0.25d0*Nh(M-1,j)-0.25d0*Nh(M-1,j+1))/(CurviNx(M,j)*NormalNx(M,j)+CurviNy(M,j)*NormalNy(M,j))/DistDualN(M,j) & 
! ! 	  +0.5d0*(NormalSx(M,j)**2+NormalSy(M,j)**2)*CellAreaS(M,j)*(diffusionH(M,j)+diffusionH(M,j-1))*(Nh(M,j)-Nh(M,j-1))/(CurviSx(M,j)*NormalSx(M,j)+CurviSy(M,j)*NormalSy(M,j))/DistS(M,j) & 
! ! 	  -0.5d0*(CurviSx(M,j)*TangentSx(M,j)+CurviSy(M,j)*TangentSy(M,j))*CellAreaS(M,j)*(0.25d0*diffusionH(M,j)+0.25d0*diffusionH(M,j-1)-0.25d0*diffusionH(M-1,j-1)-0.25d0*diffusionH(M-1,j))*(0.25d0*Nh(M,j)+0.25d0*Nh(M,j-1)-0.25d0*Nh(M-1,j-1)-0.25d0*Nh(M-1,j))/(CurviSx(M,j)*NormalSx(M,j)+CurviSy(M,j)*NormalSy(M,j))/DistDualS(M,j) & 
! ! 	  )/CellVol(M,j)*dt & 
! ! 	  +Nh(M,j)
! 	  
! 	if(TeOff.ne.1) then
! 	TeNew(M,j) = ( & 
! 	-0.5d0*(CurviEx(M,j)*TangentEx(M,j)+CurviEy(M,j)*TangentEy(M,j))*CellAreaE(M,j)*(0.5d0*kappae(M,j+1)+kappae(M,j)+0.5d0*kappae(M,j-1))*(0.5d0*Te(M,j+1)-0.5d0*Te(M,j-1))/(CurviEx(M,j)*NormalEx(M,j)+CurviEy(M,j)*NormalEy(M,j))/DistDualE(M,j) & 
! 	+0.5d0*(NormalWx(M,j)**2+NormalWy(M,j)**2)*CellAreaW(M,j)*(kappae(M,j)+kappae(M-1,j))*(Te(M,j)-Te(M-1,j))/(CurviWx(M,j)*NormalWx(M,j)+CurviWy(M,j)*NormalWy(M,j))/DistW(M,j) & 
! 	-0.5d0*(CurviWx(M,j)*TangentWx(M,j)+CurviWy(M,j)*TangentWy(M,j))*CellAreaW(M,j)*(0.5d0*kappae(M,j)+0.25d0*kappae(M,j+1)+0.5d0*kappae(M-1,j)+0.25d0*kappae(M-1,j+1)+0.25d0*kappae(M-1,j-1)+0.25d0*kappae(M,j-1))*(0.25d0*Te(M,j+1)+0.25d0*Te(M-1,j+1)-0.25d0*Te(M-1,j-1)-0.25d0*Te(M,j-1))/(CurviWx(M,j)*NormalWx(M,j)+CurviWy(M,j)*NormalWy(M,j))/DistDualW(M,j) & 
! 	+0.5d0*(NormalNx(M,j)**2+NormalNy(M,j)**2)*CellAreaN(M,j)*(kappae(M,j+1)+kappae(M,j))*(Te(M,j+1)-Te(M,j))/(CurviNx(M,j)*NormalNx(M,j)+CurviNy(M,j)*NormalNy(M,j))/DistN(M,j) & 
! 	-0.5d0*(CurviNx(M,j)*TangentNx(M,j)+CurviNy(M,j)*TangentNy(M,j))*CellAreaN(M,j)*(0.75d0*kappae(M,j+1)+0.75d0*kappae(M,j)+0.25d0*kappae(M-1,j)+0.25d0*kappae(M-1,j+1))*(0.25d0*Te(M,j+1)+0.25d0*Te(M,j)-0.25d0*Te(M-1,j)-0.25d0*Te(M-1,j+1))/(CurviNx(M,j)*NormalNx(M,j)+CurviNy(M,j)*NormalNy(M,j))/DistDualN(M,j) & 
! 	+0.5d0*(NormalSx(M,j)**2+NormalSy(M,j)**2)*CellAreaS(M,j)*(kappae(M,j)+kappae(M,j-1))*(Te(M,j)-Te(M,j-1))/(CurviSx(M,j)*NormalSx(M,j)+CurviSy(M,j)*NormalSy(M,j))/DistS(M,j) & 
! 	-0.5d0*(CurviSx(M,j)*TangentSx(M,j)+CurviSy(M,j)*TangentSy(M,j))*CellAreaS(M,j)*(0.25d0*kappae(M,j)+0.25d0*kappae(M,j-1)-0.25d0*kappae(M-1,j-1)-0.25d0*kappae(M-1,j))*(0.25d0*Te(M,j)+0.25d0*Te(M,j-1)-0.25d0*Te(M-1,j-1)-0.25d0*Te(M-1,j))/(CurviSx(M,j)*NormalSx(M,j)+CurviSy(M,j)*NormalSy(M,j))/DistDualS(M,j) & 
! 	+(-CouplingE(M,j)+SourceE(M,j))*CellVol(M,j) & 
! 	)/Ce(M,j)/CellVol(M,j)*dt & 
! 	+Te(M,j)
! 	
! 	ThNew(M,j) = ( & 
! 	-0.5d0*(CurviEx(M,j)*TangentEx(M,j)+CurviEy(M,j)*TangentEy(M,j))*CellAreaE(M,j)*(0.5d0*kappah(M,j+1)+kappah(M,j)+0.5d0*kappah(M,j-1))*(0.5d0*Th(M,j+1)-0.5d0*Th(M,j-1))/(CurviEx(M,j)*NormalEx(M,j)+CurviEy(M,j)*NormalEy(M,j))/DistDualE(M,j) & 
! 	+0.5d0*(NormalWx(M,j)**2+NormalWy(M,j)**2)*CellAreaW(M,j)*(kappah(M,j)+kappah(M-1,j))*(Th(M,j)-Th(M-1,j))/(CurviWx(M,j)*NormalWx(M,j)+CurviWy(M,j)*NormalWy(M,j))/DistW(M,j) & 
! 	-0.5d0*(CurviWx(M,j)*TangentWx(M,j)+CurviWy(M,j)*TangentWy(M,j))*CellAreaW(M,j)*(0.5d0*kappah(M,j)+0.25d0*kappah(M,j+1)+0.5d0*kappah(M-1,j)+0.25d0*kappah(M-1,j+1)+0.25d0*kappah(M-1,j-1)+0.25d0*kappah(M,j-1))*(0.25d0*Th(M,j+1)+0.25d0*Th(M-1,j+1)-0.25d0*Th(M-1,j-1)-0.25d0*Th(M,j-1))/(CurviWx(M,j)*NormalWx(M,j)+CurviWy(M,j)*NormalWy(M,j))/DistDualW(M,j) & 
! 	+0.5d0*(NormalNx(M,j)**2+NormalNy(M,j)**2)*CellAreaN(M,j)*(kappah(M,j+1)+kappah(M,j))*(Th(M,j+1)-Th(M,j))/(CurviNx(M,j)*NormalNx(M,j)+CurviNy(M,j)*NormalNy(M,j))/DistN(M,j) & 
! 	-0.5d0*(CurviNx(M,j)*TangentNx(M,j)+CurviNy(M,j)*TangentNy(M,j))*CellAreaN(M,j)*(0.75d0*kappah(M,j+1)+0.75d0*kappah(M,j)+0.25d0*kappah(M-1,j)+0.25d0*kappah(M-1,j+1))*(0.25d0*Th(M,j+1)+0.25d0*Th(M,j)-0.25d0*Th(M-1,j)-0.25d0*Th(M-1,j+1))/(CurviNx(M,j)*NormalNx(M,j)+CurviNy(M,j)*NormalNy(M,j))/DistDualN(M,j) & 
! 	+0.5d0*(NormalSx(M,j)**2+NormalSy(M,j)**2)*CellAreaS(M,j)*(kappah(M,j)+kappah(M,j-1))*(Th(M,j)-Th(M,j-1))/(CurviSx(M,j)*NormalSx(M,j)+CurviSy(M,j)*NormalSy(M,j))/DistS(M,j) & 
! 	-0.5d0*(CurviSx(M,j)*TangentSx(M,j)+CurviSy(M,j)*TangentSy(M,j))*CellAreaS(M,j)*(0.25d0*kappah(M,j)+0.25d0*kappah(M,j-1)-0.25d0*kappah(M-1,j-1)-0.25d0*kappah(M-1,j))*(0.25d0*Th(M,j)+0.25d0*Th(M,j-1)-0.25d0*Th(M-1,j-1)-0.25d0*Th(M-1,j))/(CurviSx(M,j)*NormalSx(M,j)+CurviSy(M,j)*NormalSy(M,j))/DistDualS(M,j) & 
! 	+(-CouplingH(M,j)+SourceH(M,j))*CellVol(M,j) & 
! 	)/Ch(M,j)/CellVol(M,j)*dt & 
! 	+Th(M,j)
! 	
! 		! boundary condition for Diffusion equation in finite volumes, at the APEX
!     ! 	  TeNew(1,j) = 2d0*((1d0/8d0)*(kappae(1,j+1)+kappae(1,j))*(Te(1,j+1)-Te(1,j))*(x(2,j+1)**2+2d0*x(2,j+1)*x(2,j)-2d0*x(2,j+1)*x(1,j)-2d0*x(2,j+1)*x(1,j+1)+x(2,j)**2-2d0*x(2,j)*x(1,j) &
!     ! 		-2d0*x(2,j)*x(1,j+1)+x(1,j)**2+2d0*x(1,j+1)*x(1,j)+x(1,j+1)**2+y(2,j+1)**2+2d0*y(2,j+1)*y(2,j)-2d0*y(2,j+1)*y(1,j)-2d0*y(2,j+1)*y(1,j+1)+y(2,j)**2-2d0*y(2,j)*y(1,j) &
!     ! 		-2d0*y(2,j)*y(1,j+1)+y(1,j)**2+2d0*y(1,j+1)*y(1,j)+y(1,j+1)**2)**(0.5d0)/(x(1,j+1)**2-2d0*x(1,j+1)*x(1,j)+x(1,j)**2+y(1,j+1)**2-2d0*y(1,j+1)*y(1,j)+y(1,j)**2)**(0.5d0) &
!     ! 		-(1d0/8d0)*(kappae(1,j-1)+kappae(1,j))*(Te(1,j)-Te(1,j-1))*(x(2,j)**2+2d0*x(2,j)*x(2,j-1)-2d0*x(2,j)*x(1,j)-2d0*x(2,j)*x(1,j-1)+x(2,j-1)**2-2d0*x(2,j-1)*x(1,j) &
!     ! 		-2d0*x(2,j-1)*x(1,j-1)+x(1,j)**2+2d0*x(1,j)*x(1,j-1)+x(1,j-1)**2+y(2,j)**2+2d0*y(2,j)*y(2,j-1)-2d0*y(2,j)*y(1,j)-2d0*y(2,j)*y(1,j-1)+y(2,j-1)**2-2d0*y(2,j-1)*y(1,j) &
!     ! 		-2d0*y(2,j-1)*y(1,j-1)+y(1,j)**2+2d0*y(1,j)*y(1,j-1)+y(1,j-1)**2)**(0.5d0)/(x(1,j)**2-2d0*x(1,j)*x(1,j-1)+x(1,j-1)**2+y(1,j)**2-2d0*y(1,j)*y(1,j-1)+y(1,j-1)**2)**(0.5d0) &
!     ! 		+0.5d0*(-CouplingE(1,j)+SourceE(1,j))*abs((0.25d0*x(2,j+1)+0.25d0*x(2,j)-0.25d0*x(1,j)+0.25d0*x(1,j+1)-0.5d0*x(1,j-1))*(0.25d0*y(1,j)+0.5d0*y(1,j+1)-0.25d0*y(2,j) & 
!     ! 		-0.25d0*y(2,j-1)-0.25d0*y(1,j-1))-(0.25d0*x(1,j)+0.5d0*x(1,j+1)-0.25d0*x(2,j)-0.25d0*x(2,j-1)-0.25d0*x(1,j-1))*(0.25d0*y(2,j+1)+0.25d0*y(2,j)-0.25d0*y(1,j)+0.25d0*y(1,j+1) &
!     ! 		-0.5d0*y(1,j-1))))/Ce(1,j)/abs((0.25d0*x(2,j+1)+0.25d0*x(2,j)-0.25d0*x(1,j)+0.25d0*x(1,j+1)-0.5d0*x(1,j-1))*(0.25d0*y(1,j)+0.5d0*y(1,j+1)-0.25d0*y(2,j)-0.25d0*y(2,j-1) &
!     ! 		-0.25d0*y(1,j-1))-(0.25d0*x(1,j)+0.5d0*x(1,j+1)-0.25d0*x(2,j)-0.25d0*x(2,j-1)-0.25d0*x(1,j-1))*(0.25d0*y(2,j+1)+0.25d0*y(2,j)-0.25d0*y(1,j)+0.25d0*y(1,j+1)-0.5d0*y(1,j-1)))*dt &
!     ! 		+Te(1,j)
! 
!     ! 	  ThNew(1,j) = 2d0*((1d0/8d0)*(kappah(1,j+1)+kappah(1,j))*(Th(1,j+1)-Th(1,j))*(x(2,j+1)**2+2d0*x(2,j+1)*x(2,j)-2d0*x(2,j+1)*x(1,j)-2d0*x(2,j+1)*x(1,j+1)+x(2,j)**2-2d0*x(2,j)*x(1,j) & 
!     ! 		-2d0*x(2,j)*x(1,j+1)+x(1,j)**2+2d0*x(1,j+1)*x(1,j)+x(1,j+1)**2+y(2,j+1)**2+2d0*y(2,j+1)*y(2,j)-2d0*y(2,j+1)*y(1,j)-2d0*y(2,j+1)*y(1,j+1)+y(2,j)**2-2d0*y(2,j)*y(1,j) & 
!     ! 		-2d0*y(2,j)*y(1,j+1)+y(1,j)**2+2d0*y(1,j+1)*y(1,j)+y(1,j+1)**2)**(0.5d0)/(x(1,j+1)**2-2d0*x(1,j+1)*x(1,j)+x(1,j)**2+y(1,j+1)**2-2d0*y(1,j+1)*y(1,j)+y(1,j)**2)**(0.5d0) & 
!     ! 		-(1d0/8d0)*(kappah(1,j-1)+kappah(1,j))*(Th(1,j)-Th(1,j-1))*(x(2,j)**2+2d0*x(2,j)*x(2,j-1)-2d0*x(2,j)*x(1,j)-2d0*x(2,j)*x(1,j-1)+x(2,j-1)**2-2d0*x(2,j-1)*x(1,j) & 
!     ! 		-2d0*x(2,j-1)*x(1,j-1)+x(1,j)**2+2d0*x(1,j)*x(1,j-1)+x(1,j-1)**2+y(2,j)**2+2d0*y(2,j)*y(2,j-1)-2d0*y(2,j)*y(1,j)-2d0*y(2,j)*y(1,j-1)+y(2,j-1)**2-2d0*y(2,j-1)*y(1,j) & 
!     ! 		-2d0*y(2,j-1)*y(1,j-1)+y(1,j)**2+2d0*y(1,j)*y(1,j-1)+y(1,j-1)**2)**(0.5d0)/(x(1,j)**2-2d0*x(1,j)*x(1,j-1)+x(1,j-1)**2+y(1,j)**2-2d0*y(1,j)*y(1,j-1)+y(1,j-1)**2)**(0.5d0) & 
!     ! 		+0.5d0*(-CouplingH(1,j)+SourceH(1,j))*abs((0.25d0*x(2,j+1)+0.25d0*x(2,j)-0.25d0*x(1,j)+0.25d0*x(1,j+1)-0.5d0*x(1,j-1))*(0.25d0*y(1,j)+0.5d0*y(1,j+1)-0.25d0*y(2,j) & 
!     ! 		-0.25d0*y(2,j-1)-0.25d0*y(1,j-1))-(0.25d0*x(1,j)+0.5d0*x(1,j+1)-0.25d0*x(2,j)-0.25d0*x(2,j-1)-0.25d0*x(1,j-1))*(0.25d0*y(2,j+1)+0.25d0*y(2,j)-0.25d0*y(1,j)+0.25d0*y(1,j+1) & 
!     ! 		-0.5d0*y(1,j-1))))/Ch(1,j)/abs((0.25d0*x(2,j+1)+0.25d0*x(2,j)-0.25d0*x(1,j)+0.25d0*x(1,j+1)-0.5d0*x(1,j-1))*(0.25d0*y(1,j)+0.5d0*y(1,j+1)-0.25d0*y(2,j)-0.25d0*y(2,j-1) & 
!     ! 		-0.25d0*y(1,j-1))-(0.25d0*x(1,j)+0.5d0*x(1,j+1)-0.25d0*x(2,j)-0.25d0*x(2,j-1)-0.25d0*x(1,j-1))*(0.25d0*y(2,j+1)+0.25d0*y(2,j)-0.25d0*y(1,j)+0.25d0*y(1,j+1)-0.5d0*y(1,j-1)))*dt &
!     ! 		+Th(1,j)
! 
! 	  end if
! 		  
	end do
! 	
!     !     ! CORNERS
!     !     
! 	!North-East
	! 	
! 	NeNew(M,N) = ((GainsE(M,N)-LossesE(M,N))*CellVol(M,N) &
! 		  -( &
! 		    (0.5d0*(JeX(M,N)+JeX(M,N-1))*TangentSx(M,N)+0.5d0*(JeY(M,N)+JeY(M,N-1))*TangentSy(M,N))*TangentSx(M,N)*NormalSx(M,N)*CellAreaS(M,N) & 
! 		   +(0.5d0*(JeX(M,N)+JeX(M,N-1))*TangentSx(M,N)+0.5d0*(JeY(M,N)+JeY(M,N-1))*TangentSy(M,N))*TangentSy(M,N)*NormalSy(M,N)*CellAreaS(M,N) &
! 		   +(0.5d0*(JeX(M,N)+JeX(M-1,N))*TangentWx(M,N)+0.5d0*(JeY(M,N)+JeY(M-1,N))*TangentWy(M,N))*TangentWx(M,N)*NormalWx(M,N)*CellAreaW(M,N) &
! 		   +(0.5d0*(JeX(M,N)+JeX(M-1,N))*TangentWy(M,N)+0.5d0*(JeY(M,N)+JeY(M-1,N))*TangentWy(M,N))*TangentWy(M,N)*NormalWy(M,N)*CellAreaW(M,N) &
! 		   +(0.5d0*(JeX(M,N)+JeX(M,N-1))*NormalSx(M,N)+0.5d0*(JeY(M,N)+JeY(M,N-1))*NormalSy(M,N))*NormalSx(M,N)*NormalSx(M,N)*CellAreaS(M,N) & 
! 		   +(0.5d0*(JeX(M,N)+JeX(M,N-1))*NormalSx(M,N)+0.5d0*(JeY(M,N)+JeY(M,N-1))*NormalSy(M,N))*NormalSy(M,N)*NormalSy(M,N)*CellAreaS(M,N) &
! 		   +(0.5d0*(JeX(M,N)+JeX(M-1,N))*NormalWx(M,N)+0.5d0*(JeY(M,N)+JeY(M-1,N))*NormalWy(M,N))*NormalWx(M,N)*NormalWx(M,N)*CellAreaW(M,N) &
! 		   +(0.5d0*(JeX(M,N)+JeX(M-1,N))*NormalWy(M,N)+0.5d0*(JeY(M,N)+JeY(M-1,N))*NormalWy(M,N))*NormalWy(M,N)*NormalWy(M,N)*CellAreaW(M,N) &
! ! 		   (0.5d0*(JeX(M,N)+JeX(M-1,N))*NormalWx(M,N)+0.5d0*(JeY(M,N)+JeY(M-1,N))*NormalWy(M,N))*CellAreaW(M,N) &
! ! 		  +(0.5d0*(JeX(M,N)+JeX(M,N-1))*NormalSx(M,N)+0.5d0*(JeY(M,N)+JeY(M,N-1))*NormalSy(M,N))*CellAreaS(M,N) & 
! 		  ) &
! 		  )*dt/CellVol(M,N) & 
! 		+(-(CurviEx(M,N)*TangentEx(M,N)+CurviEy(M,N)*TangentEy(M,N))*CellAreaE(M,N)*diffusionE(M,N)*(0.5d0*Ne(M,N)-0.5d0*Ne(M,N-1))/(CurviEx(M,N)*NormalEx(M,N)+CurviEy(M,N)*NormalEy(M,N))/DistDualE(M,N) & 
! 		+0.5d0*(NormalWx(M,N)**2+NormalWy(M,N)**2)*CellAreaW(M,N)*(diffusionE(M-1,N)+diffusionE(M,N))*(Ne(M,N)-Ne(M-1,N))/(CurviWx(M,N)*NormalWx(M,N)+CurviWy(M,N)*NormalWy(M,N))/DistW(M,N) & 
! 		-0.5d0*(CurviWx(M,N)*TangentWx(M,N)+CurviWy(M,N)*TangentWy(M,N))*CellAreaW(M,N)*(diffusionE(M-1,N)+diffusionE(M,N))*(0.25d0*Ne(M-1,N)+0.25d0*Ne(M,N)-0.25d0*Ne(M-1,N-1)-0.25d0*Ne(M,N-1))/(CurviWx(M,N)*NormalWx(M,N)+CurviWy(M,N)*NormalWy(M,N))/DistDualW(M,N) & 
! 		-(CurviNx(M,N)*TangentNx(M,N)+CurviNy(M,N)*TangentNy(M,N))*CellAreaN(M,N)*diffusionE(M,N)*(0.5d0*Ne(M,N)-0.5d0*Ne(M,N-1))/(CurviNx(M,N)*NormalNx(M,N)+CurviNy(M,N)*NormalNy(M,N))/DistDualN(M,N) & 
! 		+0.5d0*(NormalSx(M,N)**2+NormalSy(M,N)**2)*CellAreaS(M,N)*(diffusionE(M,N-1)+diffusionE(M,N))*(Ne(M,N)-Ne(M,N-1))/(CurviSx(M,N)*NormalSx(M,N)+CurviSy(M,N)*NormalSy(M,N))/DistS(M,N) & 
! 		-0.5d0*(CurviSx(M,N)*TangentSx(M,N)+CurviSy(M,N)*TangentSy(M,N))*CellAreaS(M,N)*(diffusionE(M,N-1)+diffusionE(M,N))*(0.25d0*Ne(M,N-1)+0.25d0*Ne(M,N)-0.25d0*Ne(M-1,N-1)-0.25d0*Ne(M-1,N))/(CurviSx(M,N)*NormalSx(M,N)+CurviSy(M,N)*NormalSy(M,N))/DistDualS(M,N) & 
! 		)/CellVol(M,N)*dt &
! 		  +Ne(M,N)
! 			
! 			
! 	NhNew(M,N) = ((GainsH(M,N)-LossesH(M,N))*CellVol(M,N) &
! 		  -( &
! 		   (0.5d0*(JhX(M,N)+JhX(M,N-1))*TangentSx(M,N)+0.5d0*(JhY(M,N)+JhY(M,N-1))*TangentSy(M,N))*TangentSx(M,N)*NormalSx(M,N)*CellAreaS(M,N) & 
! 		   +(0.5d0*(JhX(M,N)+JhX(M,N-1))*TangentSx(M,N)+0.5d0*(JhY(M,N)+JhY(M,N-1))*TangentSy(M,N))*TangentSy(M,N)*NormalSy(M,N)*CellAreaS(M,N) &
! 		   +(0.5d0*(JhX(M,N)+JhX(M-1,N))*TangentWx(M,N)+0.5d0*(JhY(M,N)+JhY(M-1,N))*TangentWy(M,N))*TangentWx(M,N)*NormalWx(M,N)*CellAreaW(M,N) &
! 		   +(0.5d0*(JhX(M,N)+JhX(M-1,N))*TangentWy(M,N)+0.5d0*(JhY(M,N)+JhY(M-1,N))*TangentWy(M,N))*TangentWy(M,N)*NormalWy(M,N)*CellAreaW(M,N) &
! 		   +(0.5d0*(JhX(M,N)+JhX(M,N-1))*NormalSx(M,N)+0.5d0*(JhY(M,N)+JhY(M,N-1))*NormalSy(M,N))*NormalSx(M,N)*NormalSx(M,N)*CellAreaS(M,N) & 
! 		   +(0.5d0*(JhX(M,N)+JhX(M,N-1))*NormalSx(M,N)+0.5d0*(JhY(M,N)+JhY(M,N-1))*NormalSy(M,N))*NormalSy(M,N)*NormalSy(M,N)*CellAreaS(M,N) &
! 		   +(0.5d0*(JhX(M,N)+JhX(M-1,N))*NormalWx(M,N)+0.5d0*(JhY(M,N)+JhY(M-1,N))*NormalWy(M,N))*NormalWx(M,N)*NormalWx(M,N)*CellAreaW(M,N) &
! 		   +(0.5d0*(JhX(M,N)+JhX(M-1,N))*NormalWy(M,N)+0.5d0*(JhY(M,N)+JhY(M-1,N))*NormalWy(M,N))*NormalWy(M,N)*NormalWy(M,N)*CellAreaW(M,N) &
! ! 		  (0.5d0*(JhX(M,N)+JhX(M-1,N))*NormalWx(M,N)+0.5d0*(JhY(M,N)+JhY(M-1,N))*NormalWy(M,N))*CellAreaW(M,N) &
! ! 		  +(0.5d0*(JhX(M,N)+JhX(M,N-1))*NormalSx(M,N)+0.5d0*(JhY(M,N)+JhY(M,N-1))*NormalSy(M,N))*CellAreaS(M,N) &
! 		  ) &
! 		  )*dt/CellVol(M,N) & 
! 		 +(-(CurviEx(M,N)*TangentEx(M,N)+CurviEy(M,N)*TangentEy(M,N))*CellAreaE(M,N)*diffusionH(M,N)*(0.5d0*Nh(M,N)-0.5d0*Nh(M,N-1))/(CurviEx(M,N)*NormalEx(M,N)+CurviEy(M,N)*NormalEy(M,N))/DistDualE(M,N) & 
! 		+0.5d0*(NormalWx(M,N)**2+NormalWy(M,N)**2)*CellAreaW(M,N)*(diffusionH(M-1,N)+diffusionH(M,N))*(Nh(M,N)-Nh(M-1,N))/(CurviWx(M,N)*NormalWx(M,N)+CurviWy(M,N)*NormalWy(M,N))/DistW(M,N) & 
! 		-0.5d0*(CurviWx(M,N)*TangentWx(M,N)+CurviWy(M,N)*TangentWy(M,N))*CellAreaW(M,N)*(diffusionH(M-1,N)+diffusionH(M,N))*(0.25d0*Nh(M-1,N)+0.25d0*Nh(M,N)-0.25d0*Nh(M-1,N-1)-0.25d0*Nh(M,N-1))/(CurviWx(M,N)*NormalWx(M,N)+CurviWy(M,N)*NormalWy(M,N))/DistDualW(M,N) & 
! 		-(CurviNx(M,N)*TangentNx(M,N)+CurviNy(M,N)*TangentNy(M,N))*CellAreaN(M,N)*diffusionH(M,N)*(0.5d0*Nh(M,N)-0.5d0*Nh(M,N-1))/(CurviNx(M,N)*NormalNx(M,N)+CurviNy(M,N)*NormalNy(M,N))/DistDualN(M,N) & 
! 		+0.5d0*(NormalSx(M,N)**2+NormalSy(M,N)**2)*CellAreaS(M,N)*(diffusionH(M,N-1)+diffusionH(M,N))*(Nh(M,N)-Nh(M,N-1))/(CurviSx(M,N)*NormalSx(M,N)+CurviSy(M,N)*NormalSy(M,N))/DistS(M,N) & 
! 		-0.5d0*(CurviSx(M,N)*TangentSx(M,N)+CurviSy(M,N)*TangentSy(M,N))*CellAreaS(M,N)*(diffusionH(M,N-1)+diffusionH(M,N))*(0.25d0*Nh(M,N-1)+0.25d0*Nh(M,N)-0.25d0*Nh(M-1,N-1)-0.25d0*Nh(M-1,N))/(CurviSx(M,N)*NormalSx(M,N)+CurviSy(M,N)*NormalSy(M,N))/DistDualS(M,N) & 
! 		)/CellVol(M,N)*dt &
! 		+Nh(M,N)
! 	  if(TeOff.ne.1) then
!     	  TeNew(M,N) = (-(CurviEx(M,N)*TangentEx(M,N)+CurviEy(M,N)*TangentEy(M,N))*CellAreaE(M,N)*kappae(M,N)*(0.5d0*Te(M,N)-0.5d0*Te(M,N-1))/(CurviEx(M,N)*NormalEx(M,N)+CurviEy(M,N)*NormalEy(M,N))/DistDualE(M,N) & 
! 		+0.5d0*(NormalWx(M,N)**2+NormalWy(M,N)**2)*CellAreaW(M,N)*(kappae(M-1,N)+kappae(M,N))*(Te(M,N)-Te(M-1,N))/(CurviWx(M,N)*NormalWx(M,N)+CurviWy(M,N)*NormalWy(M,N))/DistW(M,N) & 
! 		-0.5d0*(CurviWx(M,N)*TangentWx(M,N)+CurviWy(M,N)*TangentWy(M,N))*CellAreaW(M,N)*(kappae(M-1,N)+kappae(M,N))*(0.25d0*Te(M-1,N)+0.25d0*Te(M,N)-0.25d0*Te(M-1,N-1)-0.25d0*Te(M,N-1))/(CurviWx(M,N)*NormalWx(M,N)+CurviWy(M,N)*NormalWy(M,N))/DistDualW(M,N) & 
! 		-(CurviNx(M,N)*TangentNx(M,N)+CurviNy(M,N)*TangentNy(M,N))*CellAreaN(M,N)*kappae(M,N)*(0.5d0*Te(M,N)-0.5d0*Te(M,N-1))/(CurviNx(M,N)*NormalNx(M,N)+CurviNy(M,N)*NormalNy(M,N))/DistDualN(M,N) & 
! 		+0.5d0*(NormalSx(M,N)**2+NormalSy(M,N)**2)*CellAreaS(M,N)*(kappae(M,N-1)+kappae(M,N))*(Te(M,N)-Te(M,N-1))/(CurviSx(M,N)*NormalSx(M,N)+CurviSy(M,N)*NormalSy(M,N))/DistS(M,N) & 
! 		-0.5d0*(CurviSx(M,N)*TangentSx(M,N)+CurviSy(M,N)*TangentSy(M,N))*CellAreaS(M,N)*(kappae(M,N-1)+kappae(M,N))*(0.25d0*Te(M,N-1)+0.25d0*Te(M,N)-0.25d0*Te(M-1,N-1)-0.25d0*Te(M-1,N))/(CurviSx(M,N)*NormalSx(M,N)+CurviSy(M,N)*NormalSy(M,N))/DistDualS(M,N) & 
! 		+(-CouplingE(M,N)+SourceE(M,N))*CellVol(M,N))/Ce(M,N)/CellVol(M,N)*dt & 
! 		+Te(M,N)
! 		
! 	  ThNew(M,N) = (-(CurviEx(M,N)*TangentEx(M,N)+CurviEy(M,N)*TangentEy(M,N))*CellAreaE(M,N)*kappah(M,N)*(0.5d0*Th(M,N)-0.5d0*Th(M,N-1))/(CurviEx(M,N)*NormalEx(M,N)+CurviEy(M,N)*NormalEy(M,N))/DistDualE(M,N) & 
! 		+0.5d0*(NormalWx(M,N)**2+NormalWy(M,N)**2)*CellAreaW(M,N)*(kappah(M-1,N)+kappah(M,N))*(Th(M,N)-Th(M-1,N))/(CurviWx(M,N)*NormalWx(M,N)+CurviWy(M,N)*NormalWy(M,N))/DistW(M,N) & 
! 		-0.5d0*(CurviWx(M,N)*TangentWx(M,N)+CurviWy(M,N)*TangentWy(M,N))*CellAreaW(M,N)*(kappah(M-1,N)+kappah(M,N))*(0.25d0*Th(M-1,N)+0.25d0*Th(M,N)-0.25d0*Th(M-1,N-1)-0.25d0*Th(M,N-1))/(CurviWx(M,N)*NormalWx(M,N)+CurviWy(M,N)*NormalWy(M,N))/DistDualW(M,N) & 
! 		-(CurviNx(M,N)*TangentNx(M,N)+CurviNy(M,N)*TangentNy(M,N))*CellAreaN(M,N)*kappah(M,N)*(0.5d0*Th(M,N)-0.5d0*Th(M,N-1))/(CurviNx(M,N)*NormalNx(M,N)+CurviNy(M,N)*NormalNy(M,N))/DistDualN(M,N) & 
! 		+0.5d0*(NormalSx(M,N)**2+NormalSy(M,N)**2)*CellAreaS(M,N)*(kappah(M,N-1)+kappah(M,N))*(Th(M,N)-Th(M,N-1))/(CurviSx(M,N)*NormalSx(M,N)+CurviSy(M,N)*NormalSy(M,N))/DistS(M,N) & 
! 		-0.5d0*(CurviSx(M,N)*TangentSx(M,N)+CurviSy(M,N)*TangentSy(M,N))*CellAreaS(M,N)*(kappah(M,N-1)+kappah(M,N))*(0.25d0*Th(M,N-1)+0.25d0*Th(M,N)-0.25d0*Th(M-1,N-1)-0.25d0*Th(M-1,N))/(CurviSx(M,N)*NormalSx(M,N)+CurviSy(M,N)*NormalSy(M,N))/DistDualS(M,N) & 
! 		+(-CouplingE(M,N)+SourceE(M,N))*CellVol(M,N))/Ce(M,N)/CellVol(M,N)*dt & 
! 		+Th(M,N)
! 	end if
! 	!South-East
	
! 	NeNew(M,1) = ((GainsE(M,1)-LossesE(M,1))*CellVol(M,1) &
! 		  -( &
! 		   (0.5d0*(JeX(M,1)+JeX(M-1,1))*TangentWx(M,1)+0.5d0*(JeY(M-1,1)+JeY(M,1))*TangentWy(M,1))*TangentWx(M,1)*NormalWx(M,1)*CellAreaW(M,1) &
! 		  +(0.5d0*(JeX(M,1)+JeX(M-1,1))*TangentWx(M,1)+0.5d0*(JeY(M-1,1)+JeY(M,1))*TangentWy(M,1))*TangentWy(M,1)*NormalWy(M,1)*CellAreaW(M,1) &
! 		  +(0.5d0*( JeX(M,1)+JeX(M,2) )*TangentNx(M,1)+0.5d0*( JeY(M,1)+JeY(M,2) )*TangentNy(M,1))*TangentNx(M,1)*NormalNx(M,1)*CellAreaN(M,1) &
! 		  +(0.5d0*( JeX(M,1)+JeX(M,2) )*TangentNx(M,1)+0.5d0*( JeY(M,1)+JeY(M,2) )*TangentNy(M,1))*TangentNy(M,1)*NormalNy(M,1)*CellAreaN(M,1) &
! 		  +(0.5d0*(JeX(M,1)+JeX(M-1,1))*NormalWx(M,1)+0.5d0*(JeY(M-1,1)+JeY(M,1))*NormalWy(M,1))*NormalWx(M,1)*NormalWx(M,1)*CellAreaW(M,1) &
! 		  +(0.5d0*(JeX(M,1)+JeX(M-1,1))*NormalWx(M,1)+0.5d0*(JeY(M-1,1)+JeY(M,1))*NormalWy(M,1))*NormalWy(M,1)*NormalWy(M,1)*CellAreaW(M,1) &
! 		  +(0.5d0*( JeX(M,1)+JeX(M,2) )*NormalNx(M,1)+0.5d0*( JeY(M,1)+JeY(M,2) )*NormalNy(M,1))*NormalNx(M,1)*NormalNx(M,1)*CellAreaN(M,1) &
! 		  +(0.5d0*( JeX(M,1)+JeX(M,2) )*NormalNx(M,1)+0.5d0*( JeY(M,1)+JeY(M,2) )*NormalNy(M,1))*NormalNy(M,1)*NormalNy(M,1)*CellAreaN(M,1) &
! ! 		  (0.5d0*(JeX(M,1)+JeX(M-1,1))*NormalWx(M,1)+0.5d0*(JeY(M,1)+JeY(M-1,1))*NormalWy(M,1))*CellAreaW(M,1) &
! ! 		  +(0.5d0*(JeX(M,1)+JeX(M,2))*NormalNx(M,1)+0.5d0*(JeY(M,1)+JeY(M,2))*NormalNy(M,1))*CellAreaN(M,1) &
! 		  ) &
! 		  )*dt/CellVol(M,1)+Ne(M,1)
! 			
! 			
! 	NhNew(M,1) = ((GainsH(M,1)-LossesH(M,1))*CellVol(M,1) &
! 		    -( &
! 		     (0.5d0*(JhX(M,1)+JhX(M-1,1))*TangentWx(M,1)+0.5d0*(JhY(M-1,1)+JhY(M,1))*TangentWy(M,1))*TangentWx(M,1)*NormalWx(M,1)*CellAreaW(M,1) &
! 		    +(0.5d0*(JhX(M,1)+JhX(M-1,1))*TangentWx(M,1)+0.5d0*(JhY(M-1,1)+JhY(M,1))*TangentWy(M,1))*TangentWy(M,1)*NormalWy(M,1)*CellAreaW(M,1) &
! 		    +(0.5d0*( JhX(M,1)+JhX(M,2) )*TangentNx(M,1)+0.5d0*( JhY(M,1)+JhY(M,2) )*TangentNy(M,1))*TangentNx(M,1)*NormalNx(M,1)*CellAreaN(M,1) &
! 		    +(0.5d0*( JhX(M,1)+JhX(M,2) )*TangentNx(M,1)+0.5d0*( JhY(M,1)+JhY(M,2) )*TangentNy(M,1))*TangentNy(M,1)*NormalNy(M,1)*CellAreaN(M,1) &
! 		    +(0.5d0*(JhX(M,1)+JhX(M-1,1))*NormalWx(M,1)+0.5d0*(JhY(M-1,1)+JhY(M,1))*NormalWy(M,1))*NormalWx(M,1)*NormalWx(M,1)*CellAreaW(M,1) &
! 		    +(0.5d0*(JhX(M,1)+JhX(M-1,1))*NormalWx(M,1)+0.5d0*(JhY(M-1,1)+JhY(M,1))*NormalWy(M,1))*NormalWy(M,1)*NormalWy(M,1)*CellAreaW(M,1) &
! 		    +(0.5d0*( JhX(M,1)+JhX(M,2) )*NormalNx(M,1)+0.5d0*( JhY(M,1)+JhY(M,2) )*NormalNy(M,1))*NormalNx(M,1)*NormalNx(M,1)*CellAreaN(M,1) &
! 		    +(0.5d0*( JhX(M,1)+JhX(M,2) )*NormalNx(M,1)+0.5d0*( JhY(M,1)+JhY(M,2) )*NormalNy(M,1))*NormalNy(M,1)*NormalNy(M,1)*CellAreaN(M,1) &
! ! 		    (0.5d0*(JhX(M,1)+JhX(M-1,1))*NormalWx(M,1)+0.5d0*(JhY(M,1)+JhY(M-1,1))*NormalWy(M,1))*CellAreaW(M,1) &
! ! 		    +(0.5d0*(JhX(M,1)+JhX(M,2))*NormalNx(M,1)+0.5d0*(JhY(M,1)+JhY(M,2))*NormalNy(M,1))*CellAreaN(M,1) & 
! 		    ) &
! 		    )*dt/CellVol(M,1)+Nh(M,1)
! 	if(TeOff.ne.1) then
! ! 	TeNew(M,1) = (-(CurviEx(M,1)*TangentEx(M,1)+CurviEy(M,1)*TangentEy(M,1))*CellAreaE(M,1)*kappae(M,1)*(-0.5d0*Te(M,1)+0.5d0*Te(M,2))/(CurviEx(M,1)*NormalEx(M,1)+CurviEy(M,1)*NormalEy(M,1))/DistDualE(M,1)+0.5d0*(NormalWx(M,1)**2+NormalWy(M,1)**2)*CellAreaW(M,1)*(kappae(M,1)+kappae(M-1,1))*(Te(M,1)-Te(M-1,1))/(CurviWx(M,1)*NormalWx(M,1)+CurviWy(M,1)*NormalWy(M,1))/DistW(M,1) & 
! ! 	      -0.5d0*(CurviWx(M,1)*TangentWx(M,1)+CurviWy(M,1)*TangentWy(M,1))*CellAreaW(M,1)*(kappae(M,1)+kappae(M-1,1))*(0.25d0*Te(M-1,2)-0.25d0*Te(M-1,1)+0.25d0*Te(M,2)-0.25d0*Te(M,1))/(CurviWx(M,1)*NormalWx(M,1)+CurviWy(M,1)*NormalWy(M,1))/DistDualW(M,1) & 
! ! 	      +0.5d0*(NormalNx(M,1)**2+NormalNy(M,1)**2)*CellAreaN(M,1)*(kappae(M,1)+kappae(M,2))*(Te(M,2)-Te(M,1))/(CurviNx(M,1)*NormalNx(M,1)+CurviNy(M,1)*NormalNy(M,1))/DistN(M,1) & 
! ! 	      -0.5d0*(CurviNx(M,1)*TangentNx(M,1)+CurviNy(M,1)*TangentNy(M,1))*CellAreaN(M,1)*(kappae(M,1)+kappae(M,2))*(0.25d0*Te(M,1)+0.25d0*Te(M,2)-0.25d0*Te(M-1,2)-0.25d0*Te(M-1,1))/(CurviNx(M,1)*NormalNx(M,1)+CurviNy(M,1)*NormalNy(M,1))/DistDualN(M,1) & 
! ! 	      -(CurviSx(M,1)*TangentSx(M,1)+CurviSy(M,1)*TangentSy(M,1))*CellAreaS(M,1)*kappae(M,1)*(0.5d0*Te(M,1)-0.5d0*Te(M-1,1))/(CurviSx(M,1)*NormalSx(M,1)+CurviSy(M,1)*NormalSy(M,1))/DistDualS(M,1) & 
! ! 	      +(-CouplingE(M,1)+SourceE(M,1))*CellVol(M,1) & 
! ! 	      )/Ce(M,1)*dt/CellVol(M,1)+Te(M,1)
! 	end if
! 	!South-West
	
! 	
! ! 	NeNew(1,1) = ((GainsE(1,1)-LossesE(1,1))*CellVol(1,1) &
! ! 		  -( &
! !     ! 	       (0.5d0*(JeX(2,1)+JeX(1,1))*NormalEx(1,1)+0.5d0*(JeY(2,1)+JeY(1,1))*NormalEy(1,1))*CellAreaE(1,1) &
! !     ! 	      +(0.5d0*(JeX(1,1)+JeX(1,2))*NormalNx(1,1)+0.5d0*(JeY(1,1)+JeY(1,2))*NormalNy(1,1))*CellAreaN(1,1) &
! ! 		   (0.5d0*(JeX(1,1)+JeX(1,2))*TangentNx(1,1)+(0.5d0*(JeY(1,1)+JeY(1,2)))*TangentNy(1,1))*TangentNx(1,1)*NormalNx(1,1)*CellAreaN(1,1) & 
! ! 		  +(0.5d0*(JeX(1,1)+JeX(1,2))*TangentNx(1,1)+(0.5d0*(JeY(1,1)+JeY(1,2)))*TangentNy(1,1))*TangentNy(1,1)*NormalNy(1,1)*CellAreaN(1,1) & 
! ! 		  +(0.5d0*(JeX(1,1)+JeX(2,1))*TangentEx(1,1)+(0.5d0*(JeY(1,1)+JeY(2,1)))*TangentEy(1,1))*TangentEx(1,1)*NormalEx(1,1)*CellAreaE(1,1) & 
! ! 		  +(0.5d0*(JeX(1,1)+JeX(2,1))*TangentEx(1,1)+(0.5d0*(JeY(1,1)+JeY(2,1)))*TangentEy(1,1))*TangentEy(1,1)*NormalEy(1,1)*CellAreaE(1,1) &
! ! 		  +(0.5d0*(JeX(1,1)+JeX(1,2))*NormalNx(1,1)+(0.5d0*(JeY(1,1)+JeY(1,2)))*NormalNy(1,1))*NormalNx(1,1)*NormalNx(1,1)*CellAreaN(1,1) & 
! ! 		  +(0.5d0*(JeX(1,1)+JeX(1,2))*NormalNx(1,1)+(0.5d0*(JeY(1,1)+JeY(1,2)))*NormalNy(1,1))*NormalNy(1,1)*NormalNy(1,1)*CellAreaN(1,1) & 
! ! 		  +(0.5d0*(JeX(1,1)+JeX(2,1))*NormalEx(1,1)+(0.5d0*(JeY(1,1)+JeY(2,1)))*NormalEy(1,1))*NormalEx(1,1)*NormalEx(1,1)*CellAreaE(1,1) & 
! ! 		  +(0.5d0*(JeX(1,1)+JeX(2,1))*NormalEx(1,1)+(0.5d0*(JeY(1,1)+JeY(2,1)))*NormalEy(1,1))*NormalEy(1,1)*NormalEy(1,1)*CellAreaE(1,1) &
! ! 		  ) &
! ! 		  )*dt/CellVol(1,1) & 
! ! 		  +(0.5d0*(NormalEx(1,1)**2+NormalEy(1,1)**2)*CellAreaE(1,1)*(diffusionE(2,1)+diffusionE(1,1))*(Ne(2,1)-Ne(1,1))/(CurviEx(1,1)*NormalEx(1,1)+CurviEy(1,1)*NormalEy(1,1))/DistE(1,1) & 
! ! 		  -0.5d0*(CurviEx(1,1)*TangentEx(1,1)+CurviEy(1,1)*TangentEy(1,1))*CellAreaE(1,1)*(diffusionE(2,1)+diffusionE(1,1))*(-0.25d0*Ne(1,1)+0.25d0*Ne(1,2)-0.25d0*Ne(2,1)+0.25d0*Ne(2,2))/(CurviEx(1,1)*NormalEx(1,1)+CurviEy(1,1)*NormalEy(1,1))/DistDualE(1,1) & 
! ! 		  -(CurviWx(1,1)*TangentWx(1,1)+CurviWy(1,1)*TangentWy(1,1))*CellAreaW(1,1)*diffusionE(1,1)*(-0.5d0*Ne(1,1)+0.5d0*Ne(1,2))/(CurviWx(1,1)*NormalWx(1,1)+CurviWy(1,1)*NormalWy(1,1))/DistDualW(1,1) & 
! ! 		  +0.5d0*(NormalNx(1,1)**2+NormalNy(1,1)**2)*CellAreaN(1,1)*(diffusionE(1,1)+diffusionE(1,2))*(Ne(1,2)-Ne(1,1))/(CurviNx(1,1)*NormalNx(1,1)+CurviNy(1,1)*NormalNy(1,1))/DistN(1,1) & 
! ! 		  -0.5d0*(CurviNx(1,1)*TangentNx(1,1)+CurviNy(1,1)*TangentNy(1,1))*CellAreaN(1,1)*(diffusionE(1,1)+diffusionE(1,2))*(-0.25d0*Ne(1,1)-0.25d0*Ne(1,2)+0.25d0*Ne(2,1)+0.25d0*Ne(2,2))/(CurviNx(1,1)*NormalNx(1,1)+CurviNy(1,1)*NormalNy(1,1))/DistDualN(1,1) & 
! ! 		  -(CurviSx(1,1)*TangentSx(1,1)+CurviSy(1,1)*TangentSy(1,1))*CellAreaS(1,1)*diffusionE(1,1)*(-0.5d0*Ne(1,1)+0.5d0*Ne(2,1))/(CurviSx(1,1)*NormalSx(1,1)+CurviSy(1,1)*NormalSy(1,1))/DistDualS(1,1) & 
! ! 		  +(-LossesE(1,1)+GainsE(1,1))*CellVol(1,1) & 
! ! 		    )*dt/CellVol(1,1) &
! ! 		  +Ne(1,1)
! 			
! 			
! ! 	NhNew(1,1) = ((GainsH(1,1)-LossesH(1,1))*CellVol(1,1) &
! ! 		  -( &
! ! 		   (0.5d0*(JhX(1,1)+JhX(1,2))*TangentNx(1,1)+(0.5d0*(JhY(1,1)+JhY(1,2)))*TangentNy(1,1))*TangentNx(1,1)*NormalNx(1,1)*CellAreaN(1,1) & 
! ! 		  +(0.5d0*(JhX(1,1)+JhX(1,2))*TangentNx(1,1)+(0.5d0*(JhY(1,1)+JhY(1,2)))*TangentNy(1,1))*TangentNy(1,1)*NormalNy(1,1)*CellAreaN(1,1) & 
! ! 		  +(0.5d0*(JhX(1,1)+JhX(2,1))*TangentEx(1,1)+(0.5d0*(JhY(1,1)+JhY(2,1)))*TangentEy(1,1))*TangentEx(1,1)*NormalEx(1,1)*CellAreaE(1,1) & 
! ! 		  +(0.5d0*(JhX(1,1)+JhX(2,1))*TangentEx(1,1)+(0.5d0*(JhY(1,1)+JhY(2,1)))*TangentEy(1,1))*TangentEy(1,1)*NormalEy(1,1)*CellAreaE(1,1) &
! ! 		  +(0.5d0*(JhX(1,1)+JhX(1,2))*NormalNx(1,1)+(0.5d0*(JhY(1,1)+JhY(1,2)))*NormalNy(1,1))*NormalNx(1,1)*NormalNx(1,1)*CellAreaN(1,1) & 
! ! 		  +(0.5d0*(JhX(1,1)+JhX(1,2))*NormalNx(1,1)+(0.5d0*(JhY(1,1)+JhY(1,2)))*NormalNy(1,1))*NormalNy(1,1)*NormalNy(1,1)*CellAreaN(1,1) & 
! ! 		  +(0.5d0*(JhX(1,1)+JhX(2,1))*NormalEx(1,1)+(0.5d0*(JhY(1,1)+JhY(2,1)))*NormalEy(1,1))*NormalEx(1,1)*NormalEx(1,1)*CellAreaE(1,1) & 
! ! 		  +(0.5d0*(JhX(1,1)+JhX(2,1))*NormalEx(1,1)+(0.5d0*(JhY(1,1)+JhY(2,1)))*NormalEy(1,1))*NormalEy(1,1)*NormalEy(1,1)*CellAreaE(1,1) &
! ! ! 		  (0.5d0*(JhX(2,1)+JhX(1,1))*NormalEx(1,1)+0.5d0*(JhY(2,1)+JhY(1,1))*NormalEy(1,1))*CellAreaE(1,1) &
! ! ! 		  +(0.5d0*(JhX(1,1)+JhX(1,2))*NormalNx(1,1)+0.5d0*(JhY(1,1)+JhY(1,2))*NormalNy(1,1))*CellAreaN(1,1) & 
! ! 		  ) &
! ! 		  )*dt/CellVol(1,1) & 
! ! 		  +(0.5d0*(NormalEx(1,1)**2+NormalEy(1,1)**2)*CellAreaE(1,1)*(diffusionH(2,1)+diffusionH(1,1))*(Nh(2,1)-Nh(1,1))/(CurviEx(1,1)*NormalEx(1,1)+CurviEy(1,1)*NormalEy(1,1))/DistE(1,1) & 
! ! 		  -0.5d0*(CurviEx(1,1)*TangentEx(1,1)+CurviEy(1,1)*TangentEy(1,1))*CellAreaE(1,1)*(diffusionH(2,1)+diffusionH(1,1))*(-0.25d0*Nh(1,1)+0.25d0*Nh(1,2)-0.25d0*Nh(2,1)+0.25d0*Nh(2,2))/(CurviEx(1,1)*NormalEx(1,1)+CurviEy(1,1)*NormalEy(1,1))/DistDualE(1,1) & 
! ! 		  -(CurviWx(1,1)*TangentWx(1,1)+CurviWy(1,1)*TangentWy(1,1))*CellAreaW(1,1)*diffusionH(1,1)*(-0.5d0*Nh(1,1)+0.5d0*Nh(1,2))/(CurviWx(1,1)*NormalWx(1,1)+CurviWy(1,1)*NormalWy(1,1))/DistDualW(1,1) & 
! ! 		  +0.5d0*(NormalNx(1,1)**2+NormalNy(1,1)**2)*CellAreaN(1,1)*(diffusionH(1,1)+diffusionH(1,2))*(Nh(1,2)-Nh(1,1))/(CurviNx(1,1)*NormalNx(1,1)+CurviNy(1,1)*NormalNy(1,1))/DistN(1,1) & 
! ! 		  -0.5d0*(CurviNx(1,1)*TangentNx(1,1)+CurviNy(1,1)*TangentNy(1,1))*CellAreaN(1,1)*(diffusionH(1,1)+diffusionH(1,2))*(-0.25d0*Nh(1,1)-0.25d0*Nh(1,2)+0.25d0*Nh(2,1)+0.25d0*Nh(2,2))/(CurviNx(1,1)*NormalNx(1,1)+CurviNy(1,1)*NormalNy(1,1))/DistDualN(1,1) & 
! ! 		  -(CurviSx(1,1)*TangentSx(1,1)+CurviSy(1,1)*TangentSy(1,1))*CellAreaS(1,1)*diffusionH(1,1)*(-0.5d0*Nh(1,1)+0.5d0*Nh(2,1))/(CurviSx(1,1)*NormalSx(1,1)+CurviSy(1,1)*NormalSy(1,1))/DistDualS(1,1) & 
! ! 		  +(-LossesH(1,1)+GainsH(1,1))*CellVol(1,1) & 
! ! 		    )*dt/CellVol(1,1) &
! ! 		  +Nh(1,1)
!       if(TeOff.ne.1) then
! !     TeNew(1,1) = (0.5d0*(NormalEx(1,1)**2+NormalEy(1,1)**2)*CellAreaE(1,1)*(kappae(2,1)+kappae(1,1))*(Te(2,1)-Te(1,1))/(CurviEx(1,1)*NormalEx(1,1)+CurviEy(1,1)*NormalEy(1,1))/DistE(1,1) & 
! ! 		-0.5d0*(CurviEx(1,1)*TangentEx(1,1)+CurviEy(1,1)*TangentEy(1,1))*CellAreaE(1,1)*(kappae(2,1)+kappae(1,1))*(-0.25d0*Te(1,1)+0.25d0*Te(1,2)-0.25d0*Te(2,1)+0.25d0*Te(2,2))/(CurviEx(1,1)*NormalEx(1,1)+CurviEy(1,1)*NormalEy(1,1))/DistDualE(1,1) & 
! ! 		-(CurviWx(1,1)*TangentWx(1,1)+CurviWy(1,1)*TangentWy(1,1))*CellAreaW(1,1)*kappae(1,1)*(-0.5d0*Te(1,1)+0.5d0*Te(1,2))/(CurviWx(1,1)*NormalWx(1,1)+CurviWy(1,1)*NormalWy(1,1))/DistDualW(1,1) & 
! ! 		+0.5d0*(NormalNx(1,1)**2+NormalNy(1,1)**2)*CellAreaN(1,1)*(kappae(1,1)+kappae(1,2))*(Te(1,2)-Te(1,1))/(CurviNx(1,1)*NormalNx(1,1)+CurviNy(1,1)*NormalNy(1,1))/DistN(1,1) & 
! ! 		-0.5d0*(CurviNx(1,1)*TangentNx(1,1)+CurviNy(1,1)*TangentNy(1,1))*CellAreaN(1,1)*(kappae(1,1)+kappae(1,2))*(-0.25d0*Te(1,1)-0.25d0*Te(1,2)+0.25d0*Te(2,1)+0.25d0*Te(2,2))/(CurviNx(1,1)*NormalNx(1,1)+CurviNy(1,1)*NormalNy(1,1))/DistDualN(1,1) & 
! ! 		-(CurviSx(1,1)*TangentSx(1,1)+CurviSy(1,1)*TangentSy(1,1))*CellAreaS(1,1)*kappae(1,1)*(-0.5d0*Te(1,1)+0.5d0*Te(2,1))/(CurviSx(1,1)*NormalSx(1,1)+CurviSy(1,1)*NormalSy(1,1))/DistDualS(1,1) & 
! ! 		+(-CouplingE(1,1)+SourceE(1,1))*CellVol(1,1) & 
! ! 		  )/Ce(1,1)*dt/CellVol(1,1)+Te(1,1)
!       end if
! 	!North-West
	 	
! 	!north face, tangent is applied on W, and S of those cells. 
! 	NeNew(1,N) = ((GainsE(1,N)-LossesE(1,N))*CellVol(1,N) &
! 		  -( &
!     ! 	       (0.5d0*(JeX(2,N)+JeX(1,N)  )*NormalEx(1,N)+0.5d0*(JeY(2,N)+JeY(1,  N))*NormalEy(1,N))*CellAreaE(1,N) & !with accumulation problem
!     ! 	      +(0.5d0*(JeX(1,N)+JeX(1,N-1))*NormalSx(1,N)+0.5d0*(JeY(1,N)+JeY(1,N-1))*NormalSy(1,N))*CellAreaS(1,N) & !with accumulation problem
! 		  ( (0.5d0*(JeX(2,N)+JeX(1,N))*TangentEx(1,N) + 0.5d0*(JeY(2,N)+JeY(1,N))*TangentEy(1,N)) * TangentEx(1,N) * NormalEx(1,N) & ! fixed, east point
! 		  + (0.5d0*(JeX(2,N)+JeX(1,N))*TangentEx(1,N) + 0.5d0*(JeY(2,N)+JeY(1,N))*TangentEy(1,N)) * TangentEy(1,N) * NormalEy(1,N) &
! 		  )*CellAreaE(1,N) &
! 		  +( (0.5d0*(JeX(1,N)+JeX(1,N-1))*TangentSx(1,N)+0.5d0*(JeY(1,N)+JeY(1,N-1))*TangentSy(1,N))*TangentSx(1,N) * NormalSx(1,N) &
! 		    +(0.5d0*(JeX(1,N)+JeX(1,N-1))*TangentSx(1,N)+0.5d0*(JeY(1,N)+JeY(1,N-1))*TangentSy(1,N))*TangentSy(1,N) * NormalSy(1,N) &
! 		  )*CellAreaS(1,N) &
! 		  +( (0.5d0*(JeX(2,N)+JeX(1,N))*NormalEx(1,N) + 0.5d0*(JeY(2,N)+JeY(1,N))*NormalEy(1,N)) * NormalEx(1,N) * NormalEx(1,N) & ! fixed, east point
! 		  +  (0.5d0*(JeX(2,N)+JeX(1,N))*NormalEx(1,N) + 0.5d0*(JeY(2,N)+JeY(1,N))*NormalEy(1,N)) * NormalEy(1,N) * NormalEy(1,N) &
! 		  )*CellAreaE(1,N) &
! 		  +( (0.5d0*(JeX(1,N)+JeX(1,N-1))*NormalSx(1,N)+0.5d0*(JeY(1,N)+JeY(1,N-1))*NormalSy(1,N))*NormalSx(1,N) * NormalSx(1,N) &
! 		    +(0.5d0*(JeX(1,N)+JeX(1,N-1))*NormalSx(1,N)+0.5d0*(JeY(1,N)+JeY(1,N-1))*NormalSy(1,N))*NormalSy(1,N) * NormalSy(1,N) &
! 		  )*CellAreaS(1,N) &
! 		  ))*dt/CellVol(1,N)+Ne(1,N)
! 			
! 			
! 	NhNew(1,N) = ((GainsH(1,N)-LossesH(1,N))*CellVol(1,N) &
! 		  -( &
!     ! 	       (0.5d0*(JhX(2,N)+JhX(1,N)  )*NormalEx(1,N)+0.5d0*(JhY(2,N)+JhY(1,N  ))*NormalEy(1,N))*CellAreaE(1,N) &
!     ! 	      +(0.5d0*(JhX(1,N)+JhX(1,N-1))*NormalSx(1,N)+0.5d0*(JhY(1,N)+JhY(1,N-1))*NormalSy(1,N))*CellAreaS(1,N) &
! 		  ((0.5d0*(JhX(2,N)+JhX(1,N))*TangentEx(1,N)+0.5d0*(JhY(2,N)+JhY(1,N))*TangentEy(1,N))*TangentEx(1,N)*NormalEx(1,N) &
! 		  +(0.5d0*(JhX(2,N)+JhX(1,N))*TangentEx(1,N)+0.5d0*(JhY(2,N)+JhY(1,N))*TangentEy(1,N))*TangentEy(1,N)*NormalEy(1,N) &
! 		  )*CellAreaE(1,N) &
! 		  +((0.5d0*(JhX(1,N)+JhX(1,N-1))*TangentSx(1,N)+0.5d0*(JhY(1,N)+JhY(1,N-1))*TangentSy(1,N))*TangentSx(1,N)*NormalSx(1,N)  & 
! 		  + (0.5d0*(JhX(1,N)+JhX(1,N-1))*TangentSx(1,N)+0.5d0*(JhY(1,N)+JhY(1,N-1))*TangentSy(1,N))*TangentSy(1,N)*NormalSy(1,N) &
! 		  )*CellAreaS(1,N) &
! 		  +((0.5d0*(JhX(2,N)+JhX(1,N))*NormalEx(1,N)+0.5d0*(JhY(2,N)+JhY(1,N))*NormalEy(1,N))*NormalEx(1,N)*NormalEx(1,N) &
! 		  +(0.5d0*(JhX(2,N)+JhX(1,N))*NormalEx(1,N)+0.5d0*(JhY(2,N)+JhY(1,N))*NormalEy(1,N))*NormalEy(1,N)*NormalEy(1,N) &
! 		  )*CellAreaE(1,N) &
! 		  +((0.5d0*(JhX(1,N)+JhX(1,N-1))*NormalSx(1,N)+0.5d0*(JhY(1,N)+JhY(1,N-1))*NormalSy(1,N))*NormalSx(1,N)*NormalSx(1,N)  & 
! 		  + (0.5d0*(JhX(1,N)+JhX(1,N-1))*NormalSx(1,N)+0.5d0*(JhY(1,N)+JhY(1,N-1))*NormalSy(1,N))*NormalSy(1,N)*NormalSy(1,N) &
! 		  )*CellAreaS(1,N) &
! 		  ))*dt/CellVol(1,N)+Nh(1,N)
! 	if(TeOff.ne.1) then
! ! 	TeNew(1,N) = &
! ! 		    (&
! ! 		      0.5d0*(NormalEx(1,N)**2+NormalEy(1,N)**2)*CellAreaE(1,N)*(kappae(1,N)+kappae(2,N))*(Te(2,N)-Te(1,N))/(NormalEx(1,N)*CurviEx(1,N)+NormalEy(1,N)*CurviEy(1,N))/DistE(1,N) & 
! ! 		     -0.5d0*(CurviEx(1,N)*TangentEx(1,N)+CurviEy(1,N)*TangentEy(1,N))*CellAreaE(1,N)*(kappae(1,N)+kappae(2,N))*(0.25d0*Te(1,N)+0.25d0*Te(2,N)-0.25d0*Te(2,N-1)-0.25d0*Te(1,N-1))/(NormalEx(1,N)*CurviEx(1,N)+NormalEy(1,N)*CurviEy(1,N))/DistDualE(1,N) & 
! ! 		     -(CurviWx(1,N)*TangentWx(1,N)+CurviWy(1,N)*TangentWy(1,N))*CellAreaW(1,N)*kappae(1,N)*(0.5d0*Te(1,N)-0.5d0*Te(1,N-1))/(NormalWx(1,N)*CurviWx(1,N)+NormalWy(1,N)*CurviWy(1,N))/DistDualW(1,N) & 
! ! 		     -(CurviNx(1,N)*TangentNx(1,N)+CurviNy(1,N)*TangentNy(1,N))*CellAreaN(1,N)*kappae(1,N)*(-0.5d0*Te(1,N)+0.5d0*Te(2,N))/(NormalNx(1,N)*CurviNx(1,N)+NormalNy(1,N)*CurviNy(1,N))/DistDualN(1,N) & 
! ! 		     +0.5d0*(NormalSx(1,N)**2+NormalSy(1,N)**2)*CellAreaS(1,N)*(kappae(1,N)+kappae(1,N-1))*(Te(1,N)-Te(1,N-1))/(NormalSx(1,N)*CurviSx(1,N)+NormalSy(1,N)*CurviSy(1,N))/DistS(1,N) & 
! ! 		     -0.5d0*(CurviSx(1,N)*TangentSx(1,N)+CurviSy(1,N)*TangentSy(1,N))*CellAreaS(1,N)*(kappae(1,N)+kappae(1,N-1))*(0.25d0*Te(2,N-1)+0.25d0*Te(2,N)-0.25d0*Te(1,N-1)-0.25d0*Te(1,N))/(NormalSx(1,N)*CurviSx(1,N)+NormalSy(1,N)*CurviSy(1,N))/DistDualS(1,N)& 
! ! 		     +(-CouplingE(1,N)+SourceE(1,N)&)*CellVol(1,N) &
! ! 		     )/Ce(1,N)/CellVol(1,N)*dt+Te(1,N)
! 	end if
! !   end if !end drift yes
	      
	      
	      
    ! writing the final results
  
    do i=1,M

      do j=1, N
      
	if((MaxHeating(i,j) < Ts(i,j)) .AND. (t > 100d0*tau)) then
	   MaxHeating(i,j)=Ts(i,j)
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
	
	if(maxTe < TeNew(i,j)) then 
	  maxTe=TeNew(i,j)
	end if
	
	if(minTe > TeNew(i,j)) then 
	  minTe=TeNew(i,j)
	end if
	
	if(maxTh < ThNew(i,j)) then 
	  maxTh=ThNew(i,j)
	end if	
	
	if(minTh > ThNew(i,j)) then 
	  minTh=ThNew(i,j)
	end if

	if(maxTs < TsNew(i,j)) then 
	  maxTs=TsNew(i,j)
	end if
	
	if(minTs > TsNew(i,j)) then 
	  minTs=TsNew(i,j)
	end if
	
	if(maxNe < NeNew(i,j)) then 
	  maxNe=NeNew(i,j)
	end if
	
	if(minNe > NeNew(i,j)) then 
	  minNe=NeNew(i,j)
	end if
	
	if(maxNh < NhNew(i,j)) then 
	  maxNh=NhNew(i,j)
	end if
	
	if(minNh > NhNew(i,j)) then 
	  minNh=NhNew(i,j)
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
	
! 	if(ConductivityFix.eq.-1) then
! 	   TeNew(i,j)=Tout; ThNew(i,j)=Tout; 
! 	end if
	
	TotalNumOfE=TotalNumOfE + TotalElectrons(i,j)
	TotalNumOfH=TotalNumOfH + TotalHoles(i,j)
	TotalThermalEnergy=TotalThermalEnergy+ThermalEnergy(i,j)
	TotalLaserEnergy=TotalLaserEnergy+LaserEnergy(i,j)
	
	if((Te(i,j)).eq."NaN") then 
	  write(95,*) "Divergence of Te at t=", t, "x,y=", x, y
	  Diverged=.true.
	end if
	if((Th(i,j)).eq."NaN") then 
	  write(95,*) "Divergence of Th at t=", t, "x,y=", x, y
	  Diverged=.true.
	end if
	if((Ts(i,j)).eq."NaN") then 
	  write(95,*) "Divergence of Ts at t=", t, "x,y=", x, y
	  Diverged=.true.
	end if
	if((Ne(i,j)).eq."NaN") then 
	  write(95,*) "Divergence of Ne at t=", t, "x,y=", x, y
	  Diverged=.true.
! 	  stop !doesnt works using OpenMP
	end if
	if((Nh(i,j)).eq."NaN") then 
	  write(95,*) "Divergence of Nh at t=", t, "x,y=", x, y
	  Diverged=.true.
! 	  stop !doesnt works using OpenMP
	end if

	if(maxCFLxT.gt.1 .OR. maxCFLyT.gt.1) then 
	  write(95,*) "Bad convergence for Te,Th. t=", t, "(CFLx,CFLy)=", maxCFLxT, maxCFLyT
! 	  Diverged=.true.
! 	  stop !doesnt works using OpenMP
	end if
	if(maxCFLxN.gt.1 .OR. maxCFLyN.gt.1) then 
	  write(95,*) "Bad convergence for Ne,Nh. t=", t, "(CFLx,CFLy)=", maxCFLxN, maxCFLyN
! 	  Diverged=.true.
! 	  stop !doesnt works using OpenMP
	end if
	
	if(FermiIndexE(i,j) > real(FermiMaxLines) .OR. FermiIndexE(i,j) < 1d0) then
	  write(*,*) "t,i,j,FermiIndexE(i,j)=", t,i,j,FermiIndexE(i,j)
	end if
	if(FermiIndexH(i,j) > real(FermiMaxLines) .OR. FermiIndexH(i,j) < 1d0) then
	  write(*,*) "t,i,j,FermiIndexH(i,j)=", t,i,j,FermiIndexH(i,j)
	end if
	
	if(FermiRatioE(i,j) < 0d0 .OR. FermiRatioH(i,j) < 0d0) then
	  write(*,*) "Problem in DOS or Ne. DOS(i,j)=", i,j,DOSe(i,j), DOSh(i,j), "Ne,h(i,j)=", Ne(i,j), Nh(i,j)
	end if
	
	call flush(95)
	
	if(Diverged==.true.) then 
	  stop
	end if

              ! lets change dt when fast reponse is finished in order to catch the long one. 
        if(AdaptativeTimeStep.eq.1) then
	  if((t>1d1/2d0*tau) .AND. (dt.eq.dt0) .AND. (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then 
	    dt=10d0*dt0
	  else if ((t > 1d2*tau/2d0) .AND. (dt.eq.10d0*dt0) .AND. (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
	    dt=1d2*dt0
	  else if ((t > 1d3*tau/2d0) .AND. (dt.eq.1d2*dt0) .AND. (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs < maxCFL)) then
	    dt=1d3*dt0
! 	  else if (maxCFLxN+maxCFLyN + maxCFLxT + maxCFLyT + maxCFLxTs+maxCFLyTs > maxCFL) then
! 	    dt=dt/1d1
	  end if
	end if
	
	diffNe(i,j)=(NeNew(i,j)-Ne(i,j))/dt
	diffNh(i,j)=(NhNew(i,j)-Nh(i,j))/dt
	
		
	if(mod(nbiter,iterOut).eq.0) then 
	  write(97,887, advance="yes") t, x(i,j), y(i,j), intensity(i,j), Te(i,j), & !5
			Th(i,j), Ts(i,j), Ne(i,j), Nh(i,j), reflectivity(i,j), & !10
			absorptionDrudeE(i,j), absorptionDrudeH(i,j), diffNe(i,j), diffNh(i,j), TotalElectrons(i,j), & !15
			TotalHoles(i,j), real(FermiIndexE(i,j)), REAL(FermiIndexH(i,j)), FermiRatioE(i,j), FermiRatioH(i,j), & !20
			SourceE(i,j), SourceH(i,j), GainsE(i,j), GainsH(i,j), LossesE(i,j), & !25
			LossesH(i,j), real(DielectricDrudeE(i,j)), aimag(DielectricDrudeE(i,j)), Egap(i,j), real(Dielectric(i,j)), & !30
			aimag(Dielectric(i,j)), MaxHeatingTime(i,j), MaxHeating(i,j), real(potentialNeedle(i,j)), Ex(i,j), & !35
			Ey(i,j), diffusionE(i,j), diffusionH(i,j), intensity2(i,j), GradNeX(i,j), & !40
			GradNeY(i,j), real(EintField(i,j)), aimag(EintField(i,j)), EintFieldR(i,j), EintFieldI(i,j) !45
			
  887	FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)	
!	  write(97,886, advance='yes')
	  call flush(97); 
	end if
	
      end do !on Y
      
       if(mod(nbiter,iterOut).eq.0) then 
	write(97,886, advance="yes")
886	FORMAT (3x)
       end if
    end do !on X

    ! writing result for the electrostatic calculations
    if(mod(nbiter,iterOut).eq.0) then 
      do i=1,Mp
	do j=1,Np
	  ! ecriture des donnees dans un fichier different
	  write(101,889, advance="yes") t, xP(i,j), yP(i,j), real(potential(i,j)), real(ExPoisson(i,j)), & !
					real(EyPoisson(i,j)), DielectricStatic(i,j), NeP(i,j), NhP(i,j)
  889	FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, &
		  1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5) 
	  call flush(101); 
	end do
      end do
      
!       write(91,'(100E14.5)') potential !Xvector !Amatrix
!     write(91,*)
      
    end if
       
    ! output to files
    call cpu_time(calc_time_3)
    if(mod(nbiter,iterOut).eq.0) then 
      
      cpuefficiency=real(nbiter)/(calc_time_3-calc_time_begin)*real(nthreads)
            
      write(98,888, advance="YES") t, maxTe, maxTh, maxTs, maxNe, &
		    maxNh, maxIntensity, TotalLaserEnergy, TotalThermalEnergy, cpuefficiency, &
		    maxSourceE, maxGainsE, maxSourceH, maxGainsH, maxGap, &
		    maxDiffNe, maxDiffNh, TotalNumOfE, TotalNumOfH, real(maxFermiIndexE), &
		    real(maxFermiIndexH), NeTotal, NhTotal 
		    
888	FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1F12.8, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
		
	write(94,884, advance="YES") t, Te(1,N/2), Th(1,N/2), Ts(1,N/2), Ne(1,N/2), &
	      Nh(1,N/2), intensity(1,N/2), TotalLaserEnergy, TotalThermalEnergy, cpuefficiency, &
	      SourceE(1,N/2), GainsE(1,N/2), SourceH(1,N/2), GainsH(1,N/2), Egap(1,N/2), &
	      diffNe(1,N/2), diffNh(1,N/2), real(FermiIndexE(1,N/2)), real(FermiIndexH(1,N/2))
	      
884	FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)

	write(93,883, advance="YES") t, Te(M/2,N), Th(M/2,N), Ts(M/2,N), Ne(M/2,N), &
	      Nh(M/2,N), intensity(M/2,N), TotalLaserEnergy, TotalThermalEnergy, cpuefficiency, &
	      SourceE(M/2,N), GainsE(M/2,N), SourceH(M/2,N), GainsH(M/2,N), Egap(M/2,N), &
	      diffNe(M/2,N), diffNh(M/2,N), real(FermiIndexE(M/2,N)), real(FermiIndexH(M/2,N))
	      
883	FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)

	write(92,882, advance="YES") t, Te(M/2,1), Th(M/2,1), Ts(M/2,1), Ne(M/2,1), &
	      Nh(M/2,1), intensity(M/2,1), TotalLaserEnergy, TotalThermalEnergy, cpuefficiency, &
	      SourceE(M/2,1), GainsE(M/2,1), SourceH(M/2,1), GainsH(M/2,1), Egap(M/2,1), &
	      diffNe(M/2,1), diffNh(M/2,1), real(FermiIndexE(M/2,1)), real(FermiIndexH(M/2,1))
	      
882	FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
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
				  xDualNE(i,j), yDualNE(i,j), xDualNW(i,j), yDualNW(i,j), xDual(i,j), & !50
				  yDual(i,j)
				  
881	FORMAT (I3, 3x, I3, 3x, 1E16.8, 3x, 1E16.8, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5)
	end do
      end do
    end if
    
    ! write the functions on Dual Mesh
    if(mod(nbiter,iterOut).eq.0) then
      do i=1,M-1
	do j=1,N-1
	  
	  write(103, 890, advance="YES") t, xDual(i,j), yDual(i,j), TeDual(i,j), ThDual(i,j), & !5
					  TsDual(i,j), NeDual(i,j), NhDual(i,j), intensityDual(i,j) !9
					  
	  890	FORMAT (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		  3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
		  
	end do
      end do	
    end if
    call flush(98); call flush(94); call flush(93); call flush(101)
  end do
  
  

  contains 
  
    function DensityOfStateE(Te)
    implicit none
    real(8) DensityOfStateE, Te
      DensityOfStateE=2d0*(meDOS*kb*Te/(2d0*pi*hbar**2))**(1.5)
      return
    end function DensityOfStateE
    
    function DensityOfStateH(Th)
    implicit none
    real(8) DensityOfStateH, Th
      DensityOfStateH=2d0*(mhDOS*kb*Th/(2d0*pi*hbar**2))**(1.5)
      return
    end function DensityOfStateH
  
    subroutine TabCreateFL !(FermiTableE, FermiTableH)
      implicit none
      integer :: unit1, unit2
      unit1=15; unit2=16
      open (unit1,file='FermiDatasE.dat')
      open (unit2,file='FermiDatasH.dat')
      read (unit1,*) FermiTableE(:,:) !, FermiTableE(2,:) !, FermiTableE(:,3), FermiTableE(:,4), &
! 	    FermiTableE(:,5), FermiTableE(:,6), FermiTableE(:,7), FermiTableE(:,8), &
! 	    FermiTableE(:,9)
      read (unit2,*) FermiTableH(:,:) !1), FermiTableH(:,2), FermiTableH(:,3), FermiTableH(:,4), &
!  	    FermiTableH(:,5), FermiTableH(:,6), FermiTableH(:,7), FermiTableH(:,8), &
! 	    FermiTableH(:,9) 
! 222	format (1F10.2, 3x, 1F10.2, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x)
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
      
      omegape=sqrt(ne*ec**2/me/epsilon0)
      DielectricFunction=epsilonInf-(omegape/omegaLaser)**2*1d0/(1d0+Imaginary*nuColl/omegaLaser)
      return
    end function DielectricFunction
    
    function DielectricFunctionDrude(density, Collision, mass)
      complex(8) :: DielectricFunctionDrude
      real(8) density, Collision, omegape, mass
      
      omegape=sqrt(density*ec**2/(mass*epsilon0))
      DielectricFunctionDrude=1d0-1d0*(omegape/omegaLaser)**2*1d0/(1d0+Imaginary*Collision/omegaLaser)
      return
    end function DielectricFunctionDrude
    
    function ephCollisionFrequency
      real(8) :: ephCollisionFrequency
      ephCollisionFrequency=4d12 ! ~240 fs coupling time
!       CollisionFrequency=1d14 ! 
      ! CollisionFrequency=1d13 !
      !CollisionFrequency=5d13 ! 
      return
    end function ephCollisionFrequency
  
    function CollisionFrequency
      real(8) :: CollisionFrequency    
      CollisionFrequency=1d15
      return
    end function CollisionFrequency
  
    function ImpactIonizationRate(Te, Ne, Ts)
      real(8) :: ImpactIonizationRate
      real(8) Te, Ne, Ts
      ImpactIonizationRate=3.6d10*exp(-1.5d0*EgapValue(Ne, Ts)/kb/Te)
      return
    end function ImpactIonizationRate
    
    function OnePhotonIonizationRate(lambda,epsilonLinear)
      real(8) :: OnePhotonIonizationRate
      real(8) lambda
      complex(8) epsilonLinear
      OnePhotonIonizationRate=4d0*pi/lambda*aimag(sqrt(epsilonLinear))
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
! 	TwoPhotonIonizationRate=0d0
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
!      EgapValue=ec*(1.16d0-(7.02d-4*Ts**2)/(Ts+1108d0)-1.5d-10*Ne**(0.33333d0)) !Driel 1987
       EgapValue=ec*(1.1692d0-4.9d-4*Ts**2/(Ts+655d0)-1.5d-10*Ne**(1d0/3d0)) !Korfiatis 2007
!       EgapValue=ec*(1.1692d0-4.9d-4*Ts**2/(Ts+655d0))
! 	EgapValue=ec*(1.1692d0) !-4.9d-4*Ts**2/(Ts+655d0))
      if(EgapValue < 0d0) then 
	EgapValue=0d0
      end if
      return
    end function EgapValue
     
    function FermiIndex(NeNc)
      implicit none
      ! Input: Value of density/DOS
      ! returns the index to take in the Fermi files
      real(8) NeNc, NeNc0, dNeNc
      integer(8) FermiIndex
      NeNc0=1d-38
      dNeNc=1.1d0 !NeNc=NeNc0*dNeNc**n
      
      FermiIndex=int(log10(NeNc/NeNc0)/log10(dNeNc)+1d0)
      if(FermiIndex < 1 .OR. FermiIndex > FermiMaxLines) then
	write(*,*) "FermiIndex problem: NeNc=", NeNc, "FermiIndex=", FermiIndex
      end if
      return
    end function FermiIndex

    function MieScattering(r, phi, radius, dielectric)
      implicit none
    
      complex(8) :: MieScattering, dielectric
      complex(8) total
		 
      real(8) :: r, phi, radius, k=2d0*pi/lambda
      real(8) ireal
      integer(8) i, j

      
      total=Zero
      ! Just test functions to validate
!       total=BesselJ(1d0, Unit*k*r) !test: success
!       total=BesselJ(-1d0, Unit*k*r) !test: success
! 	total=Hankel1(1d0, Unit*k*r) !test: succes, undefined for z=0
! 	total=Hankel1(-1d0, Unit*k*r) !test: success, undefined for z=0
! 	total=BesselJprime(1d0, Unit*k*r) !test: success
! 	total=BesselJprime(-1d0, Unit*k*r) !test: success
! 	total=Hankel1prime(1d0, Unit*k*r) !test: success
! 	total=Hankel1prime(-1d0, Unit*k*r) !test: failed, strong divergence while r->0

      do i=1, 2*maxBesselOrder+1
	ireal=real(i)-maxBesselOrder-1 !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
! 	write(*,*) i, ireal
	total=total+Imaginary**ireal * exp(Imaginary*sqrt(dielectric)*phi) * BesselJ(ireal, sqrt(dielectric)*k*r) * MieCoeff1(ireal, radius, dielectric)
      end do

      MieScattering=total
      return
    end function MieScattering

    function MieCoeff1(order, radius, dielectric)
      implicit none
      complex(8) MieCoeff1, dielectric
      real(8) :: k=2d0*pi/lambda, radius
      real(8) :: order
! 	MieCoeff1=Unit !debug
      MieCoeff1=(BesselJ(order, Unit*k*radius) - MieCoeff2(order, radius, dielectric) * Hankel1(order, Unit*k*radius)) / (BesselJ(order, k*radius*sqrt(dielectric)))
      return
    end function MieCoeff1

    function MieCoeff2(order, radius, dielectric)
      implicit none
      complex(8) :: MieCoeff2, dielectric
      real(8) :: k=2d0*pi/lambda, radius
      real(8) :: order
      complex(8) value1
	MieCoeff2=( (sqrt(dielectric) * BesselJprime(order, k*radius*sqrt(dielectric)) * BesselJ(order, Unit*k*radius) ) - (BesselJ(order,sqrt(dielectric)*k*radius)*BesselJprime(order, Unit*k*radius)) ) &
		    / (( sqrt(dielectric)*BesselJprime(order,k*radius*sqrt(dielectric))*Hankel1(order, Unit*k*radius) ) - ( BesselJ(order, sqrt(dielectric)*k*radius)*Hankel1prime(order, Unit*k*radius) ))

! 	value1=radius*sqrt(dielectric)
! 	MieCoeff2=BesselJprime(order, value1) !debug

! 	MieCoeff2=BesselJprime(order, Unit*k*radius)
!  	write(*,*) order, MieCoeff2
      return
    end function MieCoeff2

    function BesselJ(order, z)
      implicit none
      complex(8) :: z, BesselJ
      real(8) zR, zC
      real(8) :: order
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
      BesselJ=Unit*cyr(besselArray)+Imaginary*cyi(besselArray)

!       write(*,*) "Bessel", order
      
      if(order .lt. 0d0) then
	BesselJ=(-1d0)**(abs(order)) * BesselJ
      end if

      return
    end function BesselJ

    function BesselJprime(order, z)
    implicit none
      external zbesj
      complex(8) z, BesselJprime
      real(8) order
  
      BesselJprime=0.5d0*(BesselJ(order-1d0,z)-BesselJ(order+1d0,z)) !Abramovitz, Eq. (9.1.27)

!! other form of the relation
!       if(z .eq. Zero) then
! 	BesselJprime=Zero
!       else
! 	BesselJprime=order*BesselJ(order,z)/z-BesselJ(order+1d0,z) !other form (still given in Abramovitz)
!       end if
      
      !debug
!       BesselJprime=Unit
!       end if
    end function BesselJprime
    
    function Hankel1(order, z)
      implicit none
      complex(8) :: z, Hankel1
      real(8) zR, zC
      real(8) :: order
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

      if(z .eq. Zero) then
	Hankel1=Zero
      else
	if(ierr.ne.0) then
	  write(*,*) "Hankel1 is not well configured."
	  write(*,*) z, order, ierr !, ZABS(zR, zC)
	end if
      end if
      Hankel1=Unit*cyr(besselArray)+Imaginary*cyi(besselArray)

!       write(*,*) "Hankel", order, Hankel1

      if(order .lt. 0d0) then
	Hankel1=exp(Imaginary*abs(order)*pi) * Hankel1
      end if
    return
    end function Hankel1

    function Hankel1prime(order, z)
      implicit none
      external zbesh
      complex(8) z, Hankel1prime
      real(8) order
!       Hankel1prime=0.5d0*(Hankel1(order-1d0,z)-Hankel1(order+1d0,z))
! debug
	Hankel1prime=Unit
    end function Hankel1prime
    
    subroutine swap(a, b)
      real(8) temp, a, b
      temp=a
      a=b
      b=temp
    end subroutine swap
    
    function OneDindex(i,j)
      implicit none
      integer(8) i,j, OneDindex
      OneDindex=(i-1)*Np+j
      return
    end function OneDindex

    
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function FillMatrixCondition(i,j,kl,ku, Nsolve)
    ! returns boolean if position i,j has to be filled in the band diagonal matrix for Lapack solver DGBTRF/DGBTRS.
      implicit none
      integer(8):: i,j, 		&
		 ku,kl, 	&
		 Nsolve, 	& !size of the matrix A
		 FillMatrixCondition
		 
      FillMatrixCondition=((max(1,j-ku) <= i) .AND. (min(Nsolve,j+kl) >= i))
      return
    end function FillMatrixCondition
    
    function FillMatrixRow(i,j,kl,ku, Nsolve)
    ! returns row index to fill the band diagonal matrix for Lapack solver DGBTRF/DGBTRS. 
      implicit none
      integer(8):: i,j,		&!Input matrix indexes
		 ku,kl,		&!Number of upper and lower diagonals
		 Nsolve, 	&!number of rows of the full band diagonal matrix
		 FillMatrixRow	 !Output for ABsolve solve matrix
      if(FillMatrixCondition(i,j,kl,ku, Nsolve)) then
	FillMatrixRow=kl+ku+1+i-j
      else
	FillMatrixRow=-1 !induce a crash
      end if
      return
    end function FillMatrixRow
    
    function LaInv(A) result(Ainv)
    ! invert A matrix using direct inversion from Lapack
      real(8), dimension(:,:), intent(in) :: A
      real(8), dimension(size(A,1),size(A,2)) :: Ainv

      real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info, nb, ilaenv

      ! External procedures defined in LAPACK
      external ILAENV
      external DGETRF !LU decomposition preparation routine
      external DGETRI !LU inversion routine

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges (!! very slow for a sparse matrix! )
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
	stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
!       nb = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
	stop 'Matrix inversion failed!'
      end if
    end function LaInv
    
    ! Returns the inverse of a sparse matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    subroutine SolveSparse(A, X, B, KL, KU)
      real(8), dimension(:,:), intent(in) :: A
      real(8), dimension(size(A,1),1):: X, B
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info, ldb, KL, KU, i, j
      real(8), dimension(2*KL+KU+1, 2*KL+KU+1) :: AB
      integer dimAx, dimAy
      dimAx=size(A,1); dimAy=size(A,2)

      ! External procedures defined in LAPACK
      ! external DGBSV !LU inversion routine for sparse matrix

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      
      do i=1,dimAx
	do j=1,dimAy
	  if(max(1,j-KU) <= i .AND. min(M,j+KL) >= i) then
	    AB(KU+1+i-j,j)=A(i,j)
	  end if
	end do
      end do
!                   if(mod(nbiter,iterOut).eq.0) then
! 		  write(91,'(5E30.15)') AB !BsolveP !Amatrix
! 		  write(91,*)
! 	    end if
      !
      n = size(A,1)

      ! DGBSV computes the solution to a real system of linear equations
      ! A * X = B, where A is a band matrix of order N with KL subdiagonals
      ! and KU superdiagonals, and X and B are N-by-NRHS matrices.
       CALL DGBSV( n, KL, KU, 1, AB, 2*KL+KU+1, ipiv, B, n, INFO )

      if (info /= 0) then
	stop 'Matrix is numerically singular!'
      end if

      if (info /= 0) then
	stop 'Matrix inversion failed!'
      end if
    end subroutine SolveSparse

    subroutine InterpolateBiCubic(Z, xS, yS, xT, yT, dimXs, dimYs, dimXt, dimYt, Interpolated, DerivativeX, DerivativeY)
    ! use bicubic interpolation 
    
    integer(8) dimXs, dimYs, &				!dimension of the source mesh
		dimXt, dimYt				!dimensions of the target mesh
      
!       real(8) :: InterpolateBiCubic(1:dimXt,1:dimYt) 		! new value of Z
      real(8) Z(1:dimXs,1:dimYs), &				!source data
	      xS(1:dimXs,1:dimYs), yS(1:dimXs,1:dimYs), &
	      xT(1:dimXt, 1:dimYt), yT(1:dimXt, 1:dimYt), &
	      xSs(1:dimXs*dimYs), ySs(1:dimXs*dimYs), & !same in 1D matrixes
	      xTt(1:dimXt*dimYt), yTt(1:dimXt*dimYt), &
	      ShepardSource(1:dimXs*dimYs), &
	      ShepardX(1:dimXs*dimYs), ShepardY(1:dimXs*dimYs), &
	      ShepardWeight(1:dimXs*dimYs), ShepardRadius(1:dimXs*dimYs), &
	      ShepardA(1:9,dimXs*dimYs), &
	      ShepardValues(1:dimXt, 1:dimYt), &
	      ShepardNext(1:dimXs*dimYs), &
	      ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, CSTVAL, &
	      Interpolated(1:dimXt, 1:dimYt), &
	      DerivativeX(1:dimXt, 1:dimYt), &
	      DerivativeY(1:dimXt, 1:dimYt)
	      
      integer(8) iS, jS, iT, jT, &
		  ShepardCells(1:dimXt, 1:dimYt)
		  
      integer(4) ShepardError
	     
      ! 1/ transfer everything in vectors
      !$O M P PARALLEL DO
      do iS=1,dimXs
	!$O M P DO
	do jS=1,dimYs
	  ShepardX((iS-1)*dimXs+jS)=xS(iS,jS)
	  ShepardY((iS-1)*dimXs+jS)=yS(iS,jS)
	  ShepardSource((iS-1)*dimXs+jS)=Z(iS,jS)
	end do
	!$O M P END DO
      end do
      !$O M P END PARALLEL DO
      
      ! 2/ compute Cshep2
      CALL Cshep2(dimXs*dimYs, ShepardX, ShepardY, ShepardSource, 11, 15, dimXt, &
	     ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, &
	     ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA, ShepardError) 
      
      ! 3/ compute the interpolation
      if(ShepardError.eq.0) then
	!$O M P PARALLEL DO
	do i=1,dimXt
	  !$O M P DO
	  do j=1,dimYt
! 	    Interpolated(i,j)=CSTVAL(xT(i,j), yT(i,j), dimXs*dimYs, ShepardX, ShepardY, ShepardSource, dimXt, ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA)
	    call CS2GRD(xT(i,j), yT(i,j), dimXs*dimYs, ShepardX, ShepardY, ShepardSource, dimXt, ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA, Interpolated(i,j), DerivativeX(i,j), DerivativeY(i,j), ShepardError)
	  end do
	  !$O M P END DO
	end do
	!$O M P END PARALLEL DO
	
	! 4/ extract the result
! 	InterpolateBiCubic(:,:)=Interpolated(:,:)
      else 
	write(*,*) 'BicubicInterp. error code = ', ShepardError
      end if
	     
    end subroutine InterpolateBiCubic
    
    function Interpolate(Z, xS, yS, xT, yT, dimXs, dimYs, dimXt, dimYt)
    ! from Z expressed on a source mesh (xS,yS),
    ! returns Z expressed on a target mesh (xT, yT) by bilinear interapolation for convex irregular mesh
    ! [Selner & Wastermann, J. Comp. Phys 79, 1-11 (1988)]
    
      integer(8) dimXs, dimYs, &				!dimension of the source mesh
		dimXt, dimYt				!dimensions of the target mesh
      
      real(8) :: Interpolate(1:dimXt,1:dimYt) 		! new value of Z
      real(8) Z(1:dimXs,1:dimYs), &				!source data
	      xS(1:dimXs,1:dimYs), yS(1:dimXs,1:dimYs), &
	      xT(1:dimXt, 1:dimYt), yT(1:dimXt, 1:dimYt)
      real(8) alpha1, alpha2, AireElementVoisin, AireSum
      real(8) CoordMatrix(1:2,1:2),&			! small matrix for base change
	      CoordsT(1,1:2), NewCoordsT(1,1:2),&		! coordinate vectors of target point
	      CoordsS(1,1:2), NewCoordsS(1,1:2),&		! coordinate vectors of source points
	      p, q					! simple coefficients
      
      integer(8) iS, jS, iT, jT 				! source indexes, target indexes
      integer(8) Neighbours(1:4,1:2) !contain the index of the nearest points in source mesh
      
      !$O M P PARALLEL DO
      do iT=1, dimXt
	!$O M P DO
	do jT=1, dimYt
	  ! catch the neighbours of point (i,j) contained in source mesh
	  Neighbours=FindNeighbours(xT(iT,jT), yT(iT,jT), xS, yS, 4, dimXs, dimYs)
	  if(Neighbours(1,1).eq.0) then !point is outside of the found element
! 	    write(*,*) AireElementVoisin, "!=", AireSum
	    Interpolate(iT,jT)=0d0
	  else !interpolate for real
! 	    write(*,*) AireElementVoisin, "==", AireSum

	    ! the coordinate of target point is transformed into squared element
	    CoordsT(1,1)=xT(iT,jT); CoordsT(1,2)=yT(iT,jT)
! 	    CoordMatrix(1,1)=Neighbours(4,1); CoordMatrix(1,2)=Neighbours(4,2); 
! 	    CoordMatrix(2,1)=Neighbours(2,1); CoordMatrix(2,2)=Neighbours(2,2)
! il semble qu  il y a une erreur dans les deux lignes precedentes ! test
	    CoordMatrix(1,1)=Neighbours(2,1); CoordMatrix(1,2)=Neighbours(4,1)
	    CoordMatrix(2,1)=Neighbours(2,2); CoordMatrix(2,2)=Neighbours(4,2)
	    CoordMatrix=LaInv(CoordMatrix)
	    
	    NewCoordsT=transpose(matmul(CoordMatrix, transpose(CoordsT)))
	    CoordsS(1,1)=Neighbours(3,1); CoordsS(1,2)=Neighbours(3,2)
	    NewCoordsS=transpose(matmul(CoordMatrix, transpose(CoordsS)))
	    
	    if((NewCoordsS(1,1)-1d0)<=1d-15) then
	      alpha2=NewCoordsT(1,2)/(1d0+(NewCoordsT(1,1)*(NewCoordsS(1,2)-1d0)))
	    else 
	      p=0.5d0*(1d0+NewCoordsT(1,1)*(NewCoordsS(1,2)-1d0)-NewCoordsT(1,2)*(NewCoordsS(1,1)-1d0))
	      q=NewCoordsT(1,2)*(NewCoordsS(1,1)-1d0)
	      alpha2=(-p+(p**2+q)**0.5d0)/(NewCoordsS(1,1)-1d0)
	    end if
	    
	    alpha1=NewCoordsT(1,1)/(1d0+alpha2*(NewCoordsS(1,1)-1d0))
	    
	    ! use the 4 neighbours to interpolate
	    Interpolate(iT,jT)=(1d0-alpha1)*(1d0-alpha2)*Z(Neighbours(1,1),Neighbours(1,2)) &
				+ alpha1*(1d0-alpha2)*Z(Neighbours(2,1),Neighbours(2,2)) &
				+ alpha1*alpha2*Z(Neighbours(3,1),Neighbours(3,2)) &
				+ (1d0-alpha1)*alpha2*Z(Neighbours(4,1),Neighbours(4,2))
	  end if    
	end do 
	!$O M P END DO
      end do
      !$O M P END PARALLEL DO
      ! $ OMP BARRIER
    end function Interpolate
    
    function FindNeighbours(x, y, xM, yM, number, dimX, dimY)
    ! [T. Westermann, J. Comp. Phys. 101, 307-313 (1992)]
    ! returns an array(1:4,1:2) with neighboor points of (x,y)
      integer(8) 	i, j, dimX, dimY, & !source mesh dimensions
			number, & !number of nearest neighbours to return
			times
      integer(8) FindNeighbours(1:4, 1:2), minimum(1:1), InElement !neighbour index array
      real(8) xM(1:dimX,1:dimY), yM(1:dimX,1:dimY), distance(1:dimX,1:dimY) !source mesh dimensions
      real(8) x,y, & !target point to find neighbours in source mesh
	      Area1, Area2, Area3, Area4, Element1, Element2, Element3, Element4, &
	      Areas(1:4)
      times=0
      InElement=0
      distance(:,:)=1d10
      
      do i=2, dimX-1
	do j=2,dimY-1
	  ! compute distance of point (x,y) with each point of (xM, yM)
	  distance(i,j)=sqrt((x-xM(i,j))**2d0+(y-yM(i,j))**2d0)
	end do
      end do

      ! reprendre ici si jamais la geometrie a rendu caduque le theoreme "times" fois
      do while (InElement.eq.0 .AND. times<=4)
  !     write(*,*) "1st neighboor index", minloc(distance), minval(distance)
! 	write(*,'(10E12.4)') distance
	FindNeighbours(1,1:2)=minloc(distance) 
! 	write(*,*) "Direct neighbour=", FindNeighbours(1,1:2)
	!first neighboor define 4 possible cells (not on boundaries!!)
	i=FindNeighbours(1,1); j=FindNeighbours(1,2)
	
  !       if(i.eq.1 .OR. i.eq.dimX .OR. j.eq.1 .OR. j.eq.dimY) then
	! cas particuliers
  !       else 
	  ! area of the 4 possible elements calculated with triangles defined by the point to interpolate
	  Area1=AreaTri(xM(i-1,j-1),yM(i-1,j-1),xM(i,j-1),yM(i,j-1),x,y)+AreaTri(xM(i,j-1),yM(i,j-1),x,y,xM(i,j),yM(i,j)) &
	  +AreaTri(xM(i,j),yM(i,j),x,y,xM(i-1,j),yM(i-1,j)) + AreaTri(xM(i-1,j),yM(i-1,j),x,y,xM(i-1,j-1),yM(i-1,j-1))
	  
	  Area2=AreaTri(xM(i,j-1),yM(i,j-1),xM(i+1,j-1),yM(i+1,j-1),x,y)+AreaTri(xM(i+1,j-1),yM(i+1,j-1),x,y,xM(i+1,j),yM(i+1,j)) &
	  +AreaTri(xM(i+1,j),yM(i+1,j),x,y,xM(i,j),yM(i,j)) + AreaTri(xM(i,j),yM(i,j),x,y,xM(i,j-1),yM(i,j-1))
	  
	  Area3=AreaTri(xM(i,j),yM(i,j),xM(i+1,j),yM(i+1,j),x,y)+AreaTri(xM(i+1,j),yM(i+1,j),x,y,xM(i+1,j+1),yM(i+1,j+1)) &
	  +AreaTri(xM(i+1,j+1),yM(i+1,j+1),x,y,xM(i,j+1),yM(i,j+1)) + AreaTri(xM(i,j+1),yM(i,j+1),x,y,xM(i,j),yM(i,j))
		
	  Area4=AreaTri(xM(i-1,j),yM(i-1,j),xM(i,j),yM(i,j),x,y)+AreaTri(xM(i,j),yM(i,j),x,y,xM(i,j+1),yM(i,j+1)) &
	  +AreaTri(xM(i,j+1),yM(i,j+1),x,y,xM(i-1,j+1),yM(i-1,j+1)) + AreaTri(xM(i-1,j+1),yM(i-1,j+1),x,y,xM(i-1,j),yM(i-1,j))
	  
	  ! calcauler l'aire de chaque element. Celui où Area i == AreaElement(i) contient alors le noeud.
	  Element1=AreaElement(xM(i-1,j-1),yM(i-1,j-1),xM(i,j-1),yM(i,j-1),xM(i,j),yM(i,j),xM(i-1,j),yM(i-1,j))
	  Element2=AreaElement(xM(i,j-1),yM(i,j-1),xM(i+1,j-1),yM(i+1,j-1),xM(i+1,j),yM(i+1,j),xM(i,j),yM(i,j))
	  Element3=AreaElement(xM(i,j),yM(i,j),xM(i+1,j),yM(i+1,j),xM(i+1,j+1),yM(i+1,j+1),xM(i,j+1),yM(i,j+1))
	  Element4=AreaElement(xM(i-1,j),yM(i-1,j),xM(i,j),yM(i,j),xM(i,j+1),yM(i,j+1),xM(i-1,j+1),yM(i-1,j+1))
	  
	  if(abs(Area1-Element1)<1d-19) then 
	    InElement=1
	  else if(abs(Area2-Element2)<1d-19) then 
	    InElement=2
	  else if(abs(Area3-Element3)<1d-19) then 
	    InElement=3
	  else if(abs(Area4-Element4)<1d-19) then 
	    InElement=4
	  else !autrement, il y a deux solutions: 
		!SOIT le plus proche voisin est dans un element plus loin (cas des points inclus mais pas pris en compte! car maillage non regulier)
		!soit il est bel et bien externe au maillage
	    InElement=0 ! 
	    Times=Times+1
	    distance(i,j)=1d10
	  end if
	end do
	
	! debug
!  	if(i.eq.5 .AND. j.eq.20) then
! 	  write(*,*) "Sum of triangle-point areas", Area1, Area2, Area3, Area4
! 	  write(*,*) "Element areas", Element1, Element2, Element3, Element4
! 	  write(*,*) "Inside", InElement
!  	end if
      
! 	! chercher l'indice de Area tel que Area(i) = AreaElement()
      
!       minimum=minloc(Areas)
!       InElement=minimum(1) !type bugfix
      
      ! Now the element is found, we organize the summits      
      if(InElement.eq.1) then
	FindNeighbours(1,1)=i-1; 	FindNeighbours(1,2)=j-1; 
	FindNeighbours(2,1)=i; 		FindNeighbours(2,2)=j-1; 
	FindNeighbours(3,1)=i; 		FindNeighbours(3,2)=j; 
	FindNeighbours(4,1)=i-1; 	FindNeighbours(4,2)=j; 
      else if(InElement.eq.2) then
	FindNeighbours(1,1)=i; 		FindNeighbours(1,2)=j-1; 
	FindNeighbours(2,1)=i+1; 	FindNeighbours(2,2)=j-1; 
	FindNeighbours(3,1)=i+1; 	FindNeighbours(3,2)=j; 
	FindNeighbours(4,1)=i; 		FindNeighbours(4,2)=j;
      else if(InElement.eq.3) then
	FindNeighbours(1,1)=i; 		FindNeighbours(1,2)=j; 
	FindNeighbours(2,1)=i+1; 	FindNeighbours(2,2)=j; 
	FindNeighbours(3,1)=i+1; 	FindNeighbours(3,2)=j+1; 
	FindNeighbours(4,1)=i; 		FindNeighbours(4,2)=j+1; 
      else if(InElement.eq.4) then
	FindNeighbours(1,1)=i-1; 	FindNeighbours(1,2)=j;
	FindNeighbours(2,1)=i; 		FindNeighbours(2,2)=j; 
	FindNeighbours(3,1)=i; 		FindNeighbours(3,2)=j+1; 
	FindNeighbours(4,1)=i-1; 	FindNeighbours(4,2)=j+1; 
      else if(InElement.eq.0) then
	FindNeighbours(:,:)=0 !out of the mesh
      else 
	FindNeighbours(:,:)=-1
      end if
      
!       if(i.eq.2 .AND. j.eq.7) then
! 	write(*,*) InElement 
! 	write(*,*) 
! 	write(*,'(2I6.2)') FindNeighbours
! 	write(*,*) " " 
!       end if
      
    end function FindNeighbours
    
    function AreaElement(x1,y1,x2,y2,x3,y3,x4,y4)
    ! works with convex elements!
      real(8) x1, y1, x2, y2, x3,y3,x4,y4, AreaElement
!       AreaElement=0.5d0*abs((x3-x1)*(y4-y2)-(y3-y1)*(x4-x2))
      AreaElement=0.5d0*abs((x3-x1)*(y2-y4)-(y3-y1)*(x2-x4))
    end function AreaElement
    
    function AreaTri(xA, yA, xB, yB, xP, yP)
      real(8) xA, yA, xB, yB, xP, yP, AreaTri
      AreaTri=0.5d0*abs((xA-xP)*(yB-yP)-(xB-xP)*(yA-yP))
    end function AreaTri
    
    function Distance(x1, y1, x2, y2) 
      real(8) Distance, x1, y1, x2, y2
      Distance=sqrt((x2-x1)**2+(y2-y1)**2)
    end function Distance
    
    function Normal(x1,y1,x2,y2,axis)
    !compute the normal vector to the line defined by given distances in cartesian plane
      real(8) x1, y1, x2, y2
      real(8) NormalX, NormalY, Dx, Dy, Normal
      real(8) normVector
      integer(8) axis
      Dx=x1-x2; Dy=y1-y2
      normVector=sqrt(Dx**2+Dy**2)
      NormalX=-Dy/normVector
      NormalY=Dx/normVector
      if(axis.eq.1) then 
	Normal=NormalX
      else 
	Normal=NormalY
      end if
    end function Normal
    
    function Tangent(x1,y1,x2,y2,axis)
    !compute the tangent vector to the line defined by given distances in cartesian plane
      real(8) x1, y1, x2, y2
      real(8) TangentX, TangentY, Dx, Dy, Tangent
      real(8) normVector
      integer(8) axis
      Dx=x1-x2; Dy=y1-y2
      normVector=sqrt(Dx**2+Dy**2)
      TangentX=Dx/normVector
      TangentY=Dy/normVector
      if(axis.eq.1) then 
	Tangent=TangentX
      else 
	Tangent=TangentY
      end if
    end function Tangent

    function ContourYofX(x, radius, angle)
    ! returns the value of the cone radius as a function of position X
    ! assumption: cone is symmetrical by rotation around (Ox) axis
      implicit none

      real(8)::x, radius, angle, ContourYofX
      real(8) a, b
      a=radius/(tan(angle/2d0)**2)
      b=radius/(tan(angle/2d0))
      ContourYofX=b*tan(acos(a/(x+a)))
      
      return
    end function ContourYofX
    
end program Flaps

