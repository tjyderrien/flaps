! FLAPS2D
! T.JY. Derrien, CNRS/BAM 2014
! Functions for Mie scattering calculations in 3D

!     function MieScattering(r, phi, radius, dielectric)
!       implicit none
! 
!       complex(8) :: MieScattering, dielectric
!       complex(8) total
! 
!       real(8) :: r, phi, radius, k
!       real(8) ireal
!       integer(8) i, j
! 
!       k=2d0*pi/lambda
!       
!       total=Zero
!       ! Just test functions to validate
! !       total=BesselJ(1d0, Unit*k*r) !test: success
! !       total=BesselJ(-1d0, Unit*k*r) !test: success
! ! 	total=Hankel1(1d0, Unit*k*r) !test: succes, undefined for z=0
! ! 	total=Hankel1(-1d0, Unit*k*r) !test: success, undefined for z=0
! ! 	total=BesselJprime(1d0, Unit*k*r) !test: success
! ! 	total=BesselJprime(-1d0, Unit*k*r) !test: success
! ! 	total=Hankel1prime(1d0, Unit*k*r) !test: success
! ! 	total=Hankel1prime(-1d0, Unit*k*r) !test: failed, strong divergence while r->0
! 
!       do i=1, 2*maxBesselOrder+1
! 	ireal=real(i)-maxBesselOrder-1 !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
! ! 	write(*,*) i, ireal
! 	total=total+Imaginary**ireal * exp(Imaginary*ireal*phi) * BesselJ(ireal, sqrt(dielectric)*k*r) * MieCoeff1(ireal, radius, dielectric)
!       end do
! 
!       MieScattering=total
!       return
!     end function MieScattering
! 
!     function MieScatteringTE1(r, phi, radius, dielectric)
!       implicit none
! 
!       complex(8) :: MieScatteringTE1, dielectric
!       complex(8) :: total
! 
!       real(8) :: r, phi, radius, k
!       real(8) :: ireal
!       integer(8) :: i, j
! 
!       k=2d0*pi/lambda
!       
!       !!Careful !! This function is very sensitive to noise.
! 
!       total=Zero
! 
!       do i=1, 2*maxBesselOrder+1
! 	ireal=real(i)-maxBesselOrder-1 !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
! ! 	write(*,*) i, ireal
! 	total=total + ( Imaginary**ireal * exp(Imaginary*ireal*phi) * BesselJ(ireal, sqrt(dielectric)*k*r) * ireal * MieCoeff3(ireal, radius, dielectric) ) !original !!
!       end do
!       MieScatteringTE1=-total/(dielectric*k*r)
! ! 	MieScatteringTE1=Zero
!       return
!     end function MieScatteringTE1
! 
!     function MieScatteringTE2(r, phi, radius, dielectric)
!       implicit none
! 
!       complex(8) :: MieScatteringTE2, dielectric
!       complex(8) :: total
! 
!       real(8) :: r, phi, radius, k
!       real(8) :: ireal
!       integer(8) :: i, j
! 
!       k=2d0*pi/lambda
!       
!       total=Zero
! 
!       do i=1, 2*maxBesselOrder+1
! 	ireal=real(i)-maxBesselOrder-1 !ireal is included in [-tmin;tmin], but fortran does not accept loops with negative index
! ! 	write(*,*) i, ireal
! 	total=total + ( Imaginary**ireal * exp(Imaginary*ireal*phi) * BesselJprime(ireal, sqrt(dielectric)*k*r) * MieCoeff3(ireal, radius, dielectric) )
!       end do
! 
!       MieScatteringTE2=-total*Imaginary/sqrt(dielectric)
! ! 	MieScatteringTE2=Zero !debug
!       return
!     end function MieScatteringTE2
! 
!     function MieCoeff1(order, radius, dielectric)
!       implicit none
!       complex(8) MieCoeff1, dielectric
!       real(8) :: radius, k
!       real(8) :: order
! 
!       k=2d0*pi/lambda
!       
! ! 	MieCoeff1=Unit !debug
!       MieCoeff1=(BesselJ(order, Unit*k*radius) - MieCoeff2(order, radius, dielectric) * Hankel1(order, Unit*k*radius)) / (BesselJ(order, k*radius*sqrt(dielectric)))
!       return
!     end function MieCoeff1
! 
!     function MieCoeff2(order, radius, dielectric)
!       implicit none
!       complex(8) :: MieCoeff2, dielectric
!       real(8) :: k, radius
!       real(8) :: order
!       complex(8) value1
! 	MieCoeff2= ( (sqrt(dielectric) * BesselJprime(order, k*radius*sqrt(dielectric)) * BesselJ(order, Unit*k*radius) ) - (BesselJ(order,sqrt(dielectric)*k*radius)*BesselJprime(order, Unit*k*radius)) ) &
! 		    / (( sqrt(dielectric)*BesselJprime(order,k*radius*sqrt(dielectric))*Hankel1(order, Unit*k*radius) ) - ( BesselJ(order, sqrt(dielectric)*k*radius)*Hankel1prime(order, Unit*k*radius) )) !original
! 
! 	k=2d0*pi/lambda
! 		    
! ! 	value1=radius*sqrt(dielectric)
! ! 	MieCoeff2=BesselJprime(order, value1) !debug
! 
! ! 	MieCoeff2=BesselJprime(order, Unit*k*radius)
! !  	write(*,*) order, MieCoeff2
!       return
!     end function MieCoeff2

!     function MieCoeff3(order, radius, dielectric)
!       implicit none
!       complex(8) MieCoeff3, dielectric
!       real(8) :: k, radius
!       real(8) :: order
! ! 	MieCoeff3=Unit !debug
!       k=2d0*pi/lambda
!       MieCoeff3=(BesselJ(order, Unit*k*radius) - MieCoeff4(order, radius, dielectric) * Hankel1(order, Unit*k*radius)) / (BesselJ(order, k*radius*sqrt(dielectric)))
!       return
!     end function MieCoeff3
! 
!     function MieCoeff4(order, radius, dielectric)
!       implicit none
!       complex(8) :: MieCoeff4, dielectric
!       real(8) :: k, radius
!       real(8) :: order
!       complex(8) value1
! 	k=2d0*pi/lambda
! 	MieCoeff4= ( ( BesselJprime(order, k*radius*sqrt(dielectric)) * BesselJ(order, Unit*k*radius) ) - sqrt(dielectric) * (BesselJ(order,sqrt(dielectric)*k*radius)*BesselJprime(order, Unit*k*radius)) ) &
! 		    / (( BesselJprime(order,k*radius*sqrt(dielectric)) * Hankel1(order, Unit*k*radius) ) - sqrt(dielectric) * ( BesselJ(order, sqrt(dielectric)*k*radius) * Hankel1prime(order, Unit*k*radius) )) !original
! ! 	MieCoeff4=Unit
!       return
!     end function MieCoeff4
! 
!     function BesselJ(order, z)
!       implicit none
!       complex(8) :: z, BesselJ
!       real(8) zR, zC
!       real(8) :: order
!       integer(8) nz, ierr
!       real(8) cyr(1:besselArray), cyi(1:besselArray)
! 
!       external ZBESJ
! !       external ZABS
! 
!       cyr(:)=0.d0; cyi(:)=0.d0
!       ierr=0; nz=0
! 
!       zR=real(z)
!       zC=aimag(z)
! 
! !       write(*,*) (zR, zC)
! 
! !       CALL zbesj(1.d0, 0.d0, 0.d0, 1, besselArray, cyr, cyi, nz, ierr)
! 
! !       write(*,*) "Bessel 1", order
!       CALL ZBESJ(zR, zC, abs(order), 1, besselArray, cyr, cyi, nz, ierr)
! 
!       if(ierr.ne.0) then
! 	write(*,*) "BesselJ is not well configured."
! 	write(*,*) z, cyr, cyi, ierr !, ZABS(zR, zC)
!       end if
!       BesselJ=Unit*cyr(besselArray)+Imaginary*cyi(besselArray)
! 
! !       write(*,*) "Bessel", order
! 
!       if(order .lt. 0d0) then
! 	BesselJ=(-1d0)**(abs(order)) * BesselJ
!       end if
! 
!       return
!     end function BesselJ
! 
!     function BesselJprime(order, z)
!     implicit none
!       external zbesj
!       complex(8) z, BesselJprime
!       real(8) order
! 
!       BesselJprime=0.5d0*(BesselJ(order-1d0,z)-BesselJ(order+1d0,z)) !Abramovitz, Eq. (9.1.27)
! 
! !! other form of the relation
! !       if(z .eq. Zero) then
! ! 	BesselJprime=Zero
! !       else
! ! 	BesselJprime=order*BesselJ(order,z)/z-BesselJ(order+1d0,z) !other form (still given in Abramovitz)
! !       end if
! 
!       !debug
! !       BesselJprime=Unit
! !       end if
!     end function BesselJprime
! 
!     function Hankel1(order, z)
!       implicit none
!       complex(8) :: z, Hankel1
!       real(8) zR, zC
!       real(8) :: order
!       integer(8) nz, ierr
!       real(8) cyr(1:besselArray), cyi(1:besselArray)
!       external ZBESH
! !       external ZABS
!       cyr(:)=0.d0; cyi(:)=0.d0
!       ierr=0; nz=0
!       zR=real(z)
!       zC=aimag(z)
! 
! !       write(*,*) zR, zC
! 
!       CALL ZBESH(zR, zC, abs(order), 1, 1, besselArray, cyr, cyi, nz, ierr)
! 
!       if(z .eq. Zero) then
! 	Hankel1=Zero
!       else
! 	if(ierr.ne.0) then
! 	  write(*,*) "Hankel1 is not well configured."
! 	  write(*,*) z, order, ierr !, ZABS(zR, zC)
! 	end if
!       end if
!       Hankel1=Unit*cyr(besselArray)+Imaginary*cyi(besselArray)
! 
! !       write(*,*) "Hankel", order, Hankel1
! 
!       if(order .lt. 0d0) then
! 	Hankel1=exp(Imaginary*abs(order)*pi) * Hankel1
!       end if
!     return
!     end function Hankel1
! 
!     function Hankel1prime(order, z)
!       implicit none
!       external zbesh
!       complex(8) z, Hankel1prime
!       real(8) order
!       Hankel1prime=0.5d0*(Hankel1(order-1d0,z)-Hankel1(order+1d0,z))
! ! debug
! ! 	Hankel1prime=Unit
!     end function Hankel1prime

    

!   write(*,*) 'Smoothing the obtained intensity profile to reduce mesh size...'

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