! 
!     subroutine InterpolateBiCubic(Z, xS, yS, xT, yT, dimXs, dimYs, dimXt, dimYt, Interpolated, DerivativeX, DerivativeY)
!     ! use bicubic interpolation
! 
!     integer(8) dimXs, dimYs, &				!dimension of the source mesh
! 		dimXt, dimYt				!dimensions of the target mesh
! 
! !       real(8) :: InterpolateBiCubic(1:dimXt,1:dimYt) 		! new value of Z
!       real(8) Z(1:dimXs,1:dimYs), &				!source data
! 	      xS(1:dimXs,1:dimYs), yS(1:dimXs,1:dimYs), &
! 	      xT(1:dimXt, 1:dimYt), yT(1:dimXt, 1:dimYt), &
! 	      xSs(1:dimXs*dimYs), ySs(1:dimXs*dimYs), & !same in 1D matrixes
! 	      xTt(1:dimXt*dimYt), yTt(1:dimXt*dimYt), &
! 	      ShepardSource(1:dimXs*dimYs), &
! 	      ShepardX(1:dimXs*dimYs), ShepardY(1:dimXs*dimYs), &
! 	      ShepardWeight(1:dimXs*dimYs), ShepardRadius(1:dimXs*dimYs), &
! 	      ShepardA(1:9,dimXs*dimYs), &
! 	      ShepardValues(1:dimXt, 1:dimYt), &
! 	      ShepardNext(1:dimXs*dimYs), &
! 	      ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, CSTVAL, &
! 	      Interpolated(1:dimXt, 1:dimYt), &
! 	      DerivativeX(1:dimXt, 1:dimYt), &
! 	      DerivativeY(1:dimXt, 1:dimYt)
! 
!       integer(8) iS, jS, iT, jT, &
! 		  ShepardCells(1:dimXt, 1:dimYt)
! 
!       integer(4) ShepardError
! 
!       ! 1/ transfer everything in vectors
!       !$O M P PARALLEL DO
!       do iS=1,dimXs
! 	!$O M P DO
! 	do jS=1,dimYs
! 	  ShepardX((iS-1)*dimXs+jS)=xS(iS,jS)
! 	  ShepardY((iS-1)*dimXs+jS)=yS(iS,jS)
! 	  ShepardSource((iS-1)*dimXs+jS)=Z(iS,jS)
! 	end do
! 	!$O M P END DO
!       end do
!       !$O M P END PARALLEL DO
! 
!       ! 2/ compute Cshep2
!       CALL Cshep2(dimXs*dimYs, ShepardX, ShepardY, ShepardSource, 11, 15, dimXt, &
! 	     ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, &
! 	     ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA, ShepardError)
! 
!       ! 3/ compute the interpolation
!       if(ShepardError.eq.0) then
! 	!$O M P PARALLEL DO
! 	do i=1,dimXt
! 	  !$O M P DO
! 	  do j=1,dimYt
! ! 	    Interpolated(i,j)=CSTVAL(xT(i,j), yT(i,j), dimXs*dimYs, ShepardX, ShepardY, ShepardSource, dimXt, ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA)
! 	    call CS2GRD(xT(i,j), yT(i,j), dimXs*dimYs, ShepardX, ShepardY, ShepardSource, dimXt, ShepardCells, ShepardNext, ShepardXmin, ShepardYmin, ShepardDx, ShepardDy, ShepardRadiusMax, ShepardRadius, ShepardA, Interpolated(i,j), DerivativeX(i,j), DerivativeY(i,j), ShepardError)
! 	  end do
! 	  !$O M P END DO
! 	end do
! 	!$O M P END PARALLEL DO
! 
! 	! 4/ extract the result
! ! 	InterpolateBiCubic(:,:)=Interpolated(:,:)
!       else
! 	write(*,*) 'BicubicInterp. error code = ', ShepardError
!       end if
! 
!     end subroutine InterpolateBiCubic
! 
! 
! ! 
! !     function AreaElement(x1,y1,x2,y2,x3,y3,x4,y4)
! !     ! works with convex elements!
! !       real(8) x1, y1, x2, y2, x3,y3,x4,y4, AreaElement
! ! !       AreaElement=M_HALF*abs((x3-x1)*(y4-y2)-(y3-y1)*(x4-x2))
! !       AreaElement=M_HALF*abs((x3-x1)*(y2-y4)-(y3-y1)*(x2-x4))
! !     end function AreaElement
! ! 
! !     function AreaTri(xA, yA, xB, yB, xP, yP)
! !       real(8) xA, yA, xB, yB, xP, yP, AreaTri
! !       AreaTri=M_HALF*abs((xA-xP)*(yB-yP)-(xB-xP)*(yA-yP))
! !     end function AreaTri
! ! 
! !     function Distance(x1, y1, x2, y2)
! !       real(8) Distance, x1, y1, x2, y2
! !       Distance=sqrt((x2-x1)**2+(y2-y1)**2)
! !     end function Distance
! ! 
! !     function Normal(x1,y1,x2,y2,axis)
! !     !compute the normal vector to the line defined by given distances in cartesian plane
! !       real(8) x1, y1, x2, y2
! !       real(8) NormalX, NormalY, Dx, Dy, Normal
! !       real(8) normVector
! !       integer(8) axis
! !       Dx=x1-x2; Dy=y1-y2
! !       normVector=sqrt(Dx**2+Dy**2)
! !       NormalX=-Dy/normVector
! !       NormalY=Dx/normVector
! !       if(axis.eq.1) then
! ! 	Normal=NormalX
! !       else
! ! 	Normal=NormalY
! !       end if
! !     end function Normal
! ! 
! !     function Tangent(x1,y1,x2,y2,axis)
! !     !compute the tangent vector to the line defined by given distances in cartesian plane
! !       real(8) x1, y1, x2, y2
! !       real(8) TangentX, TangentY, Dx, Dy, Tangent
! !       real(8) normVector
! !       integer(8) axis
! !       Dx=x1-x2; Dy=y1-y2
! !       normVector=sqrt(Dx**2+Dy**2)
! !       TangentX=Dx/normVector
! !       TangentY=Dy/normVector
! !       if(axis.eq.1) then
! ! 	Tangent=TangentX
! !       else
! ! 	Tangent=TangentY
! !       end if
! !     end function Tangent
! 
! 
! !   write(*,*) "[TEST] Area function"
! !   write(*,*) AreaElement(5.99d-7,-2.5d-7, 7.68d-7,-2.85d-7, 6.11d-7,-2.11d-7,2.71d-7, -1.66d-7)
! !   write(*,*) "[TEST] Triangle area function"
! !   write(*,*) "5.99d-7,-2.5d-7 ; 7.68d-7,-2.85d-7 ; 6.11d-7,-2.11d-7 ; 2.71d-7, -1.66d-7"
! !   write(*,*) "P 3.26d-7, 3.82d-8"
! !   write(*,*) AreaTri(5.99d-7,-2.5d-7,7.68d-7,-2.85d-7,5.74d-7,-2.07d-7)
! !   write(*,*) AreaTri(5.74d-7,-2.07d-7,7.68d-7,-2.85d-7,6.11d-7,-2.11d-7)
! !   write(*,*) AreaTri(5.74d-7,-2.07d-7,6.11d-7,-2.11d-7, 2.71d-7, -1.66d-7)
! !   write(*,*) AreaTri(5.74d-7,-2.07d-7,2.71d-7, -1.66d-7,5.99d-7,-2.5d-7)
! !   write(*,*) "Sum=", AreaTri(5.99d-7,-2.5d-7,7.68d-7,-2.85d-7,5.74d-7,-2.07d-7) + &
! !   AreaTri(5.74d-7,-2.07d-7,7.68d-7,-2.85d-7,6.11d-7,-2.11d-7) + &
! !   AreaTri(5.74d-7,-2.07d-7,6.11d-7,-2.11d-7, 2.71d-7, -1.66d-7) + &
! !   AreaTri(5.74d-7,-2.07d-7,2.71d-7, -1.66d-7,5.99d-7,-2.5d-7)
!   
!     
!     function ContourYofX(x, radius, angle)
!     ! returns the value of the cone radius as a function of position X
!     ! assumption: cone is symmetrical by rotation around (Ox) axis
!       implicit none
! 
!       real(8)::x, radius, angle, ContourYofX
!       real(8) a, b
!       a=radius/(tan(angle/2d0)**2)
!       b=radius/(tan(angle/2d0))
!       ContourYofX=b*tan(acos(a/(x+a)))
! 
!       return
!     end function ContourYofX
! 
! 
!     function ConeExp1(t)
!     ! returns the parametric function x(t) adjusted on experimental shape given by Angela Vella, GPM Rouen, Dec (2013).
!       implicit none
! 
!       real(8)::ConeExp1, t
!       real(8) p1, p2, p3, p4, p5, q1, q2, q3, q4, x0, tc
! 
!       p1=42.15d0; p2=-62.93d0; p3=66.13d0; p4=-15.28d0; p5=0.981d0
!       q1=-2.848d0; q2=12.53d0; q3=-9.082d0; q4=5.09d0
!       x0=0.0045712d0; tc=-0.14033d0
! 
!       if(t <= 1d0 .AND. t >= -1.05d0) then !reduced cone size for calculation acceleration
! !       if(t <= 2.471556d0 .AND. t >= -3.35d0) then
! 	ConeExp1=(p1*(t-tc)**4+p2*(t-tc)**3+p3*(t-tc)**2+p4*(t-tc)+p5)/((t-tc)**4 + q1*(t-tc)**3 + q2*(t-tc)**2 + q3*(t-tc) + q4) + x0
!       end if
! 
!     end function ConeExp1
! 
! !     function ConeExp1Projection(t,t0,x)
! !       implicit none
! !       ! analytical approximated (order 2) solution of the equation ConeExp1(y)=x for a given y0 (=t0) and a given x.
! !       ! This function allows to calculate the radius of the experimental cone at each point of the mesh
! !       real(8)::ConeExp1Projection
! !       real(8) t, t0, x
! !
! !       ConeExp1Projection =  4.000000000*1D-8*(-4.225456665*1D34+2.072821280*1D40*t0**6+1.789019053*1D35*t0-1.449577628*1D40*t0**3 &
! !       +2.141440746*1D40*t0**4-4.542903643*1D39*t0**9-2.112163083*1D40*t0**5+1.088559844*1D40*t0**8-1.646040210*1D40*t0**7 &
! !       +1.464301234*1D39*t0**2+2.141745000*1D38*t0**10+sqrt(4.341025390*1D77*t0**20+3.060282445*1D77*x-1.239546760*1D74*t0**2 &
! !       +4.691833578*1D70*t0-1.791486012*1D78*t0**3+2.588361669*1D79*t0**4+3.649006560*1D80*t0**6-1.173054202*1D80*t0**5 &
! !       -7.923040723*1D80*t0**7-2.003760305*1D81*t0**11+1.748451675*1D81*t0**12-3.244744576*1D78*t0**19-1.483645894*1D80*x*t0**9 &
! !       +1.007697569*1D80*x*t0**10+8.119253020*1D78*x*t0**2-1.371653295*1D81*t0**13+1.343038650*1D81*t0**8-1.809791033*1D81*t0**9 &
! !       +2.039034289*1D81*t0**10-5.001334715*1D80*t0**15+2.177401384*1D80*t0**16+9.127445662*1D80*t0**14-7.170993259*1D79*t0**17 &
! !       +1.876491723*1D79*t0**18-1.504739035*1D76*t0**21+5.468899487*1D79*t0**13*x-4.150194227*1D79*t0**14*x+2.185788784*1D79*t0**15*x &
! !       -8.624117400*1D78*t0**16*x+2.514874386*1D78*t0**17*x+4.535801676*1D79*x*t0**4-2.103090456*1D78*x*t0+8.975722800*1D76*t0**19*x- &
! !       5.817015539*1D77*t0**18*x-7.235696433*1D69-1.066059421*1D76*t0**20*x+3.569575000*1D74*t0**21*x+1.100628950*1D80*x*t0**6- &
! !       7.682804097*1D79*x*t0**5-2.215292736*1D79*x*t0**3-1.414557090*1D80*x*t0**7-2.414310149*1D79*x*t0**11-3.251570330*1D79*x*t0**12 &
! !       +1.598734131*1D80*x*t0**8))/(6.717133238*1D31*t0**2-1.937769587*1D32*t0-1.461808351*1D32*t0**3+3.899801935*1D32*t0**4 &
! !       +3.850883251*1D31*t0**6-2.125921915*1D32*t0**5+2.490160016*1D32*t0**7-1.313896238*1D32*t0**8+5.711320000*1D30*t0**9+7.344752664*1D31)
! !
! !     end function ConeExp1Projection
! 
!     function ConeExp2(t)
!     ! returns the parametric function x(t) adjusted on experimental shape given by Angela Vella, GPM Rouen, Dec (2013).
!       implicit none
! 
!       real(8):: ConeExp2, t
!       real(8) a0, a1, a2, a3, a4, a5, a6, a7, a8, w, x0
! 
!       a0=6.486d0; a1=-6.901d0; a2=1.874d0; a3=-0.7881d0; a4=-0.2992d0; a5=-0.05438d0; a6=-0.1377d0; a7=-0.08898d0; a8=-0.05379d0; w=2.162d0; x0=-0.040850d0
! 
!       if(t <= 1.4d0 .AND. t >= -1.4d0) then
! 	ConeExp2=a0+a1*cos(1d0*w*t)+a2*cos(2d0*w*t)+a3*cos(3d0*w*t)+a4*cos(4d0*w*t)+a5*cos(5d0*w*t)+a6*cos(6d0*w*t)+a7*cos(7d0*w*t)+a8*cos(8d0*w*t)+x0
!       end if
! 
!     end function ConeExp2
! 
! !     function ConeExp2Projection(t, t0, x)
! !       implicit none
! !
! !       ! analytical approximated (order 2) solution of the equation ConeExp2(y)=x for a given y0 (=t0) and a given x.
! !       ! This function allows to calculate the radius of the experimental cone at each point of the mesh
! !       real(8)::ConeExp2Projection
! !       real(8) t, t0, x
! !
! !       ConeExp2Projection = 0.9250693802D-3 * (0.25000000D8 * x + 0.172525000D9 * cos(0.2162000000D1 * t0) + 0.19702500D8 * cos(0.6486000000D1 * t0)  &
! !       + 0.2224500D7 * cos(0.1513400000D2 * t0) - 0.46850000D8 * cos(0.4324000000D1 * t0) + 0.7480000D7 * cos(0.8648000000D1 * t0) &
! !       + 0.1344750D7 * cos(0.1729600000D2 * t0) + 0.3342500D7 * cos(0.1297200000D2 * t0) + 0.1359500D7 * cos(0.1081000000D2 * t0) &
! !       + 0.14696195D8 * sin(0.1081000000D2 * t0) * t0 - 0.202579400D9 * sin(0.4324000000D1 * t0) * t0 + 0.23258796D8 * sin(0.1729600000D2 * t0) * t0  &
! !       + 0.33665583D8 * sin(0.1513400000D2 * t0) * t0 + 0.127790415D9 * sin(0.6486000000D1 * t0) * t0 + 0.372999050D9 * sin(0.2162000000D1 * t0) * t0 &
! !       + 0.43358910D8 * sin(0.1297200000D2 * t0) * t0 + 0.64687040D8 * sin(0.8648000000D1 * t0) * t0 - 0.161128750D9) / (-0.187400D6 * sin(0.4324000000D1 * t0) &
! !       + 0.21516D5 * sin(0.1729600000D2 * t0) + 0.31143D5 * sin(0.1513400000D2 * t0) + 0.118215D6 * sin(0.6486000000D1 * t0) &
! !       + 0.345050D6 * sin(0.2162000000D1 * t0) + 0.40110D5 * sin(0.1297200000D2 * t0) + 0.59840D5 * sin(0.8648000000D1 * t0) &
! !       + 0.13595D5 * sin(0.1081000000D2 * t0))
! !
! !     end function ConeExp2Projection
! 
! 
! 
! !     function ConeExp1Newton(t, value, step)
! !     ! finds the root of the function ConeExp1=value by Newton-Raphson iteration
! !       implicit none
! !       real(8) ConeExp1Newton, t, step, value
! !       ConeExp1Newton=t-ConeExp1Equation(t, value)/( (ConeExp1Equation(t+step, value)-ConeExp1Equation(t, value))/step )
! !     end function ConeExp1Newton
! !
! !     function ConeExp2Newton(t, value, step)
! !     ! finds the root of the function ConeExp2=value by Newton-Raphson iteration
! !       implicit none
! !       real(8) ConeExp2Newton, t, step, value
! !       ConeExp2Newton=t-ConeExp2Equation(t, value)/( (ConeExp2Equation(t+step, value)-ConeExp2Equation(t, value))/step )
! !     end function ConeExp2Newton
! 
!     function ConeExp1Equation(t,value)
!     ! sets the equation to solve for Newton algorithm
!       implicit none
!       real(8) ConeExp1Equation, t, value
!       integer(8) ConeExp1
!       ConeExp1Equation=ConeExp1(t)-value
!     end function ConeExp1Equation
! 
!     function ConeExp2Equation(t, value)
!     ! sets the equation to solve for Newton algorithm
!       implicit none
!       real(8) ConeExp2Equation, t, value
!       integer(8) ConeExp2
!       ConeExp2Equation=ConeExp2(t)-value
! !       ConeExp2Equation=t**3+t-value
! !       write(*,*) "Value", value
!     end function ConeExp2Equation
! 
!  function ConeExp1Radius(t,value,step)
!     ! execute Newton algorithm to find the radius of the equivalent cylinder
!       implicit none
!       real(8) ConeExp1Radius, temp1, temp2, t, step, value, zeroin
!       integer(8) i, ConeExp1Equation
! 
! !       do i=1,NewtonIterations
! ! 	temp=ConeExp2Newton(temp,value,step)
! ! 	write(*,*) 'Solving contour radius...'
! ! 	temp1=zeroin(0d0, 2.471556d0, ConeExp1Equation, 1d-15, value)
! ! 	temp2=zeroin(-3.35d0, 0d0, ConeExp1Equation, 1d-15, value)
! 	! reduced cone size to accelerate calculations
! 	temp1=zeroin(0d0, 1d0, ConeExp1Equation, 1d-15, value)
! 	temp2=zeroin(-1.05d0, 0d0, ConeExp1Equation, 1d-15, value)
! ! 	write(*,*) 'Solve:', temp
! !       end do
! 	if(abs(temp1) < abs(temp2)) then
! 	  ConeExp1Radius=abs(temp2)
! 	else
! 	  ConeExp1Radius=abs(temp1)
! 	end if
! ! 	ConeExp1Radius=M_HALF*(abs(temp1)+abs(temp2)) !not the average, since the case r>R is divergent
! !       ConeExp2Radius=M_HALF*(temp1+temp2)
!     end function ConeExp1Radius
! 
! 
!     function ConeExp2Radius(t,value,step)
!     ! execute Newton algorithm to find the radius of the equivalent cylinder
!       implicit none
!       real(8) ConeExp2Radius, temp1, temp2, t, step, value, zeroin
!       integer(8) i, ConeExp2Equation
! 
! !       do i=1,NewtonIterations
! ! 	temp=ConeExp2Newton(temp,value,step)
! ! 	write(*,*) 'Solving contour radius...'
! 	temp1=zeroin(0d0, 1.4d0, ConeExp2Equation, 1d-15, value)
! 	temp2=zeroin(-1.4d0, 0d0, ConeExp2Equation, 1d-15, value)
! 
! 	! choose the maximum value of radius to avoid the non-physical case r>R.
! 	if(abs(temp1) < abs(temp2)) then
! 	  ConeExp2Radius=abs(temp2)
! 	else
! 	  ConeExp2Radius=abs(temp1)
! 	end if
! ! 	write(*,*) 'Solve:', temp
! !       end do
! ! 	ConeExp2Radius=M_HALF*(abs(temp1)+abs(temp2))
! !       ConeExp2Radius=M_HALF*(temp1+temp2)
!     end function ConeExp2Radius
