!calculation of intrinsic hole and electron density in semiconductor
!for given temperature, band gap, effective masses and fermi level
!calculate fermi integrals of different orders and its relations

program fermi

implicit none
  integer, parameter::fid=5
  real(8), parameter::qe=1.6d-19, egap=1.156d0, phi=1.1, ni=6.33d15, kb=1.38d-23 
  real(8), parameter::me=9.1d-31, hbar=1.05457d-34, T=100.d0
  real(8) fi1, egapJ, pi, fi2, mp, mc, n0, gamm, Te, Th, ef, Te0, neMax, TeMax, mass
  real(8) fle, fi, fi0, flh, fl, eta, etah, etae, Nc, Nv, nh, ne,na, eta0, dNe
  real(8) fi3, fi4,dTe, Ce, deta, SigmaE, De, Ke, mu0, nuColl, DeChen
  real(8) Fermi0, Fermi1, Fermi2, FermiHalf, FermiThreeHalf, FermiMenusHalf
  integer i,j,k, PrevIndex, NbLignes
  real(8) NeNc, NeNcMax, NeNc0, dNeNc

  !Te=T
  !Th=T
  egapJ=egap*qe

  mc=0.36d0*me !http://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html
  mp=0.81d0*me

  pi = 4.d0*datan(1.d0)
  PrevIndex=0
  
!  ef=0.1
!  flh=-ef*qe!-egapJ/2! hole fermi level      =-egapJ/2 in thermal equilibrium
!  fle=-ef*qe!-egapJ/2 ! electron fermi level

!  etah=flh/(kb*Th) !reduced hole fermi-level
!  etae=fle/(kb*Te) !reduced electron fermi-level

!  Nv=2.d0*(mp*kb*Th/(2.d0*pi*hbar**2))**(3.0/2.0) !effective density of valence band states
!  Nc=2.d0*(mc*kb*Te/(2.d0*pi*hbar**2))**(3.0/2.0) !effective density of conduction band states
!  
!  nh=Nv*fermi_integral(0.5, etah) !intrinsic hole density
!  ne=Nc*fermi_integral(0.5, etae) !intrinsic electron density
!  
!  n0=dsqrt(nh*ne) !average density
!  print*, nh, ne, n0
!  
!  eta=etae
!  fi4=fermi_integral(1.5, eta)/fermi_integral(0.5, eta)
!  write(*,*) eta, fi4
  
open(fid,FILE='Ch.dat',access='sequential',status='unknown')
open(987, FILE='FermiDatas2.dat', access='sequential', status='unknown')

  eta0=-80.d0	!initial for reduced potential
  deta=0.1d0	!increase for the solution of reduced potential

!   dTe=10.d0	!temperature step
!   Te0=1.d0	!initial sweep for temperature
!   TeMax=1d4	!maximum electron temperature
! 
!   ne=1.d1 	!initial sweep of density
!   dNe=10.d0	!multiplication coefficient of density
!   neMax=1*5d28 !maxdensity

!   nuColl=1d14	! collision frequency
  mass=mp !mc : eletron, mp: holes
  
! Additionnal parameters
!   mu0=qe/(mass*nuColl)


! On va balayer directement les valeurs de Ne/Nc. Inutile de balayer les parametres   
! Ne/Nc = [1d-38:1d8] On va parcourir par une loi géométrique. 

NbLignes=100000
NeNc0=1d-38
dNeNc=1.1d0
NeNcMax=1d8


NeNc=NeNc0
do i=1, NbLignes
     ! fi0=ne/Nc
     eta=eta0
     !$OMP DO
     do j=1,Nblignes !Calcul de l'integrale
        eta=eta+deta
        fi1=fermi_integral(0.5, eta)
        if (fi1>=NeNc) then
            exit
        endif
     enddo
     !$OMP END DO
     Fermi0=fermi_integral(0., eta)
     Fermi1=fermi_integral(1.,eta)
     Fermi2=fermi_integral(2.,eta)
     FermiHalf=fermi_integral(0.5,eta)
     FermiThreeHalf=fermi_integral(1.5,eta)
     FermiMenusHalf=fermi_integral(-0.5,eta)
     
!      fi4=fermi_integral(1.5, eta)/fermi_integral(0.5, eta)-eta*(1d0-(fermi_integral(1.5,eta)/fermi_integral(0.5,eta))*(fermi_integral(-0.5,eta)/fermi_integral(0.5,eta)))
!      Ce=1.5d0*kb*ne*fi4
! 
!      SigmaE=qe*ne*mu0*(fermi_integral(0.,eta)/fermi_integral(0.5,eta))
!      Ke=kb**2*SigmaE*Te/qe**2*(6d0*fermi_integral(2.0,eta)/fermi_integral(0.,eta)-4d0*(fermi_integral(1.,eta)/fermi_integral(0.,eta))**2)
!      De=Ke/Ce
!      DeChen=-kb*ne*mu0*fermi_integral(0.,eta)/fermi_integral(0.5,eta)*(eta-2*fermi_integral(1.,eta)/fermi_integral(0.,eta))
! 
!      write(fid,999) & !'(e12.5, f12.5, f12.5, f12.5, e12.5, f12.5, e12.5, e12.5, e12.5)'
! 	     ne, Te, eta, fi4, Ce, &
! 	     Ce/(kb*ne), Ke, De, DeChen, fi0
! 999	format (1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
! 		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)

! ici on construit le fichier à importer dans le code
!    if(int(FermiIndex(ne/Nc)) <> PrevIndex) then !write a line only if one change of index
!    PrevIndex=int(FermiIndex(ne/Nc))
     write(987, 998) &
	FermiIndex(NeNc), NeNc, eta, Fermi0, Fermi1, &
	Fermi2, FermiHalf, FermiThreeHalf, FermiMenusHalf
	
998	format (I8.3, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, &
		3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5, 3x, 1E12.5)
!    end if
!   enddo

  NeNc=NeNc*dNeNc
  if(NeNc > NeNcMax) then
    exit
  end if
enddo

  close(fid)


  contains
  
    function FermiIndex(NeNc)
      implicit none
      ! Input: Value of density/DOS
      ! returns the index to take in the Fermi files
      real(8) NeNc, NeNc0, dNeNc
      integer(8) FermiIndex
      NeNc0=1d-38
      dNeNc=1.1d0 !NeNc=NeNc0*dNeNc**n
      
      FermiIndex=int(log(NeNc/NeNc0)/log(dNeNc)+1d0)
      return
    end function FermiIndex

  function fermi_integral(degree, chp)

  implicit none

	real(8) fermi_integral
        real(4) degree  
	real(8) chp, dE, Emax,f, t, dt, eta
        real(8) tmin,tmax,df,dfmin,dfmax, xmax, xmin,x
	integer i , imax
   
        imax=1000
        tmin=-5
        tmax=10
        dt=(tmax-tmin)/imax

  	fermi_integral=0.d0
	f=0.d0
        eta=chp
        !print*, chp, kb, Te, eta
        do i = 1,imax-1
          t=tmin+i*dt
          x=dexp(t-dexp(-t))
          df = x*(1d0+dexp(-t))*x**degree/(1d0+dexp(x-eta))
          f = f + df
        enddo
        
        xmin=dexp(tmin-dexp(-tmin))
        dfmin = xmin*(1d0+dexp(-tmin))*xmin**degree/(1d0+dexp(xmin-eta))

        xmax=dexp(tmax-dexp(-tmax))
        dfmax = xmax*(1d0+dexp(-tmax))*xmax**degree/(1d0+dexp(xmax-eta))
        
       !gamma function
        if (degree==0.d0.or.degree==1.d0) then
	    gamm=1.d0
        else if (degree==0.5) then
            gamm=dsqrt(pi)/2.d0
        else if (degree==-0.5) then
            gamm=dsqrt(pi)
        else if (degree==1.5) then
            gamm=3.d0*dsqrt(pi)/4.d0
        else if (degree==2.0) then
            gamm=2.d0
        end if      

        fermi_integral = dt*(0.5d0*(dfmin+dfmax)+f)/gamm
            
   end function fermi_integral

end program fermi
