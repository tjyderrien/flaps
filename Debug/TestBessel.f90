program TestBessel
  
  !Check that Bessel functions are well solved
  integer, parameter:: maxBesselOrder=1
  complex, parameter:: Imaginary=(0d0,1d0), Unit=(1d0,0d0), Zero=(0d0,0d0)		! complex unity

  external ZBESJ, ZBESH
  
  write(*,*) "BesselJ=", BesselJ(1d0,Unit)
  write(*,*) "HankelJ=", Hankel1(1d0,Unit)

  contains

    function BesselJ(order, z)

      complex :: z, BesselJ
      double precision zR, zC
      double precision order
      integer nz, ierr
      double precision cyr(1:maxBesselOrder), cyi(1:maxBesselOrder)

      cyr(:)=0.d0; cyi(:)=0.d0
      ierr=0; nz=0

      zR=real(z)
      zC=aimag(z)

!       write(*,*) (zR, zC)

      CALL zbesj(zR, zC, order, 1, maxBesselOrder, cyr, cyi, nz, ierr)
!       CALL zbesh(zR, zC, order, 1, 1, maxBesselOrder, cyr, cyi, nz, ierr)
!       CALL ZBESJ(zR, zC, 0d0, 2, maxBesselOrder, cyr, cyi, nz, ierr)



      if(ierr.ne.0) then
	write(*,*) "BesselJ is not well configured."
	write(*,*) z, cyr, cyi, ierr !, ZABS(zR, zC)
      end if
      BesselJ=Unit*cyr(maxBesselOrder)+Imaginary*cyi(maxBesselOrder)

      if(order < 0d0) then
	BesselJ=(-1d0)**(abs(order))*BesselJ
      end if

      return
    end function BesselJ


    function Hankel1(order, z)

      complex :: z, Hankel1
      double precision zR, zC
      double precision order
      integer nz, ierr
      double precision cyr(1:maxBesselOrder), cyi(1:maxBesselOrder)

      cyr(:)=0.d0; cyi(:)=0.d0
      ierr=0; nz=0

      zR=real(z)
      zC=aimag(z)

!       write(*,*) (zR, zC)

!       CALL zbesj(zR, zC, order, 1, maxBesselOrder, cyr, cyi, nz, ierr)
      CALL zbesh(zR, zC, order, 1, 1, maxBesselOrder, cyr, cyi, nz, ierr)
!       CALL ZBESJ(zR, zC, 0d0, 2, maxBesselOrder, cyr, cyi, nz, ierr)



      if(ierr.ne.0) then
	write(*,*) "Hankel1 is not well configured."
	write(*,*) z, cyr, cyi, ierr !, ZABS(zR, zC)
      end if
      Hankel1=Unit*cyr(maxBesselOrder)+Imaginary*cyi(maxBesselOrder)
! 
      if(order .lt. 0d0) then
	Hankel1=exp(Imaginary*abs(order)*pi) * Hankel1
      end if

      return
    end function Hankel1

end program