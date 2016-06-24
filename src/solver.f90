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
