!------------------------------------------------------------------------------
!> @file solver.f90
!
! DESCRIPTION:
!> @brief All the routine computing the physical quantities
!
!> @author
!> Nicolas Tancogne-Dejean
!
!> @date
!> 24 Jun 2016 - Initial Version from Thibault's code
!------------------------------------------------------------------------------


! This routine computes the electronic density for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeNe( newmesh, mesh, dual, dt, InvCellVol, GainsE, LossesE, diffusionE, &
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE2, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW2, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN2, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS2 )
  use Types_m
  implicit none

  type(MeshValues),                   intent(inout) :: newmesh
  type(MeshValues),                   intent(in)    :: mesh, dual
  real(8),                            intent(in)    :: dt
  real(8), dimension(mesh%M, mesh%n), intent(in) :: InvCellVol, GainsE, LossesE, diffusionE, &
                                                    NormalW2, NormalE2, NormalN2, NormalS2, &
                                                    ShapeFactorNormalE, ShapeFactorNormalW, &
                                                    ShapeFactorNormalS, ShapeFactorNormalN, &
                                                    ShapeFactorTangentE, ShapeFactorTangentW, &
                                                    ShapeFactorTangentS, ShapeFactorTangentN


  integer :: i, j

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1 !(optimized)
    do i=2, mesh%M-1
!       do i=2, M-1
!         do j=2, N-1

     ! diffusion is separated from drift
     newmesh%Ne(i,j) = mesh%Ne(i,j) + dt*( GainsE(i,j)-LossesE(i,j) )
     newmesh%Ne(i,j) = newmesh%Ne(i,j) + 0.5d0*dt*InvCellVol(i,j)*( &
               ShapeFactorNormalE(i,j) &
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
                  *(diffusionE(i,j-1)+diffusionE(i,j))*( dual%Ne(i,j-1) - dual%Ne(i-1,j-1) ) ) )


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

    end do
  end do

end subroutine computeNe




! This routine computes the holes density for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeNh( newmesh, mesh, dual, dt, InvCellVol, GainsH, LossesH, diffusionH, &
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE2, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW2, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN2, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS2 )
  use Types_m
  implicit none

  type(MeshValues),                   intent(inout) :: newmesh
  type(MeshValues),                   intent(in)    :: mesh, dual
  real(8),                            intent(in)    :: dt
  real(8), dimension(mesh%M, mesh%n), intent(in) :: InvCellVol, GainsH, LossesH, diffusionH, &
                                                    NormalW2, NormalE2, NormalN2, NormalS2, &
                                                    ShapeFactorNormalE, ShapeFactorNormalW, &
                                                    ShapeFactorNormalS, ShapeFactorNormalN, &
                                                    ShapeFactorTangentE, ShapeFactorTangentW, &
                                                    ShapeFactorTangentS, ShapeFactorTangentN

  integer :: i, j

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1 !(optimized)
    do i=2, mesh%M-1
!       do i=2, M-1
!         do j=2, N-1

      !TODO: This can be frther optimise
      newmesh%Nh(i,j) = mesh%Nh(i,j) + dt*( GainsH(i,j)-LossesH(i,j) )
      newmesh%Nh(i,j) = newmesh%Nh(i,j) + 0.5d0*dt*InvCellVol(i,j)*( &
              ! drift
              !TODO: Warning, the commented code will have a problem of a factror ofd 0.5
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
                NormalE2(i,j)*ShapeFactorNormalE(i,j)*(diffusionH(i,j)+diffusionH(i+1,j))*(mesh%Nh(i+1,j)-mesh%Nh(i,j)) &
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
              )
    end do
  end do

end subroutine computeNh



! This routine computes the electronic temperature for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeTe( newmesh, mesh, dual, dt, InvCellVol, kappae,  CouplingE, SourceE, Ce,&
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE2, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW2, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN2, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS2 )
  use Types_m
  implicit none

  type(MeshValues),                   intent(inout) :: newmesh
  type(MeshValues),                   intent(in)    :: mesh, dual
  real(8),                            intent(in)    :: dt
  real(8), dimension(mesh%M, mesh%n), intent(in) :: InvCellVol, kappae, CouplingE, SourceE, Ce, &
                                                    NormalW2, NormalE2, NormalN2, NormalS2, &
                                                    ShapeFactorNormalE, ShapeFactorNormalW, &
                                                    ShapeFactorNormalS, ShapeFactorNormalN, &
                                                    ShapeFactorTangentE, ShapeFactorTangentW, &
                                                    ShapeFactorTangentS, ShapeFactorTangentN

  integer :: i, j

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1 !(optimized)
    do i=2, mesh%M-1
!       do i=2, M-1
!         do j=2, N-1

      !TODO: This can be further optimise
      newmesh%Te(i,j) = mesh%Te(i,j) + dt * (-CouplingE(i,j)+SourceE(i,j))/ Ce(i,j) !TODO: Do we need Ce or can we compute only its inverse?
      newmesh%Te(i,j) = newmesh%Te(i,j) + 0.5d0 / Ce(i,j) * dt * InvCellVol(i,j)*( &
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
                    *( dual%Te(i,j-1) - dual%Te(i-1,j-1) ) ) !source
    end do
  end do

end subroutine computeTe

