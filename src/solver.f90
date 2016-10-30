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
!> @file solver.f90
!
! DESCRIPTION:
!> @brief All the routines for computing the physical quantities
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
  use Maths_m
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
  do j=2, mesh%N-1
    do i=2, mesh%M-1

     ! diffusion is separated from drift
     newmesh%Ne(i,j) = mesh%Ne(i,j) + dt*( GainsE(i,j)-LossesE(i,j) )
     newmesh%Ne(i,j) = newmesh%Ne(i,j) + M_HALF*dt*InvCellVol(i,j)*( &
               ShapeFactorNormalE(i,j) &
              * ( &
                 ! direct diffusion operator over irregular mesh
                 ( diffusionE(i,j)+diffusionE(i+1,j) )*( mesh%Ne(i+1,j)-mesh%Ne(i,j) )* NormalE2(i,j) &

                 ! Cross-diffusion from
                 ! [S. Mathur, J. Murthy, A pressure-based method for unstructured
                 ! meshes, Numerical Heat Transfert, Part B 31 (1997) 195–215]

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
  !$OMP END DO

end subroutine computeNe




! This routine computes the holes density for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeNh( newmesh, mesh, dual, dt, InvCellVol, GainsH, LossesH, diffusionH, &
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE2, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW2, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN2, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS2 )
  use Maths_m
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
  do j=2, mesh%N-1
    do i=2, mesh%M-1

      !TODO: This can be further optimised

      newmesh%Nh(i,j) = mesh%Nh(i,j) + dt*( GainsH(i,j)-LossesH(i,j) )
      newmesh%Nh(i,j) = newmesh%Nh(i,j) + M_HALF*dt*InvCellVol(i,j)*( &
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
  !$OMP END DO

end subroutine computeNh



! This routine computes the electronic temperature for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeTe( newmesh, mesh, dual, dt, InvCellVol, kappae,  CouplingE, SourceE, invCe,&
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE2, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW2, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN2, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS2 )
  use Maths_m
  use Types_m
  implicit none

  type(MeshValues),                   intent(inout) :: newmesh
  type(MeshValues),                   intent(in)    :: mesh, dual
  real(8),                            intent(in)    :: dt
  real(8), dimension(mesh%M, mesh%N), intent(in) :: InvCellVol, kappae, CouplingE, SourceE, invCe, &
                                                    NormalW2, NormalE2, NormalN2, NormalS2, &
                                                    ShapeFactorNormalE, ShapeFactorNormalW, &
                                                    ShapeFactorNormalS, ShapeFactorNormalN, &
                                                    ShapeFactorTangentE, ShapeFactorTangentW, &
                                                    ShapeFactorTangentS, ShapeFactorTangentN

  integer :: i, j

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1
    do i=2, mesh%M-1

      !TODO: This can be further optimised
      newmesh%Te(i,j) = mesh%Te(i,j) + dt * (-CouplingE(i,j)+SourceE(i,j))*invCe(i,j)
      newmesh%Te(i,j) = newmesh%Te(i,j) + M_HALF *invCe(i,j) * dt * InvCellVol(i,j)*( &
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
  !$OMP END DO

end subroutine computeTe


! This routine computes the hole temperature for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeTh( newmesh, mesh, dual, dt, InvCellVol, kappah,  CouplingH, SourceH, invCh,&
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE2, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW2, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN2, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS2 )
  use Maths_m
  use Types_m
  implicit none

  type(MeshValues),                   intent(inout) :: newmesh
  type(MeshValues),                   intent(in)    :: mesh, dual
  real(8),                            intent(in)    :: dt
  real(8), dimension(mesh%M, mesh%n), intent(in) :: InvCellVol, kappah, CouplingH, SourceH, invCh, &
                                                    NormalW2, NormalE2, NormalN2, NormalS2, &
                                                    ShapeFactorNormalE, ShapeFactorNormalW, &
                                                    ShapeFactorNormalS, ShapeFactorNormalN, &
                                                    ShapeFactorTangentE, ShapeFactorTangentW, &
                                                    ShapeFactorTangentS, ShapeFactorTangentN

  integer :: i, j

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1
    do i=2, mesh%M-1

  !TODO: This can be further optimised
        newmesh%Th(i,j) = mesh%Th(i,j) + (-CouplingH(i,j)+SourceH(i,j))*invCh(i,j)*dt
        newmesh%Th(i,j) = newmesh%Th(i,j)+ M_HALF*( &
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
              )*invCh(i,j)*dt*InvCellVol(i,j)

    end do
  end do
  !$OMP END DO

end subroutine computeTh


! This routine computes the lattice temperature for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeTs( newmesh, mesh, dual, dt, InvCellVol, kappas,  CouplingH, CouplingE, &
                      h1, h2, h3, invCs, TsPrev, TsOld, CellVol, &
                      ShapeFactorNormalE, ShapeFactorTangentE, NormalE2, &
                      ShapeFactorNormalW, ShapeFactorTangentW, NormalW2, &
                      ShapeFactorNormalN, ShapeFactorTangentN, NormalN2, &
                      ShapeFactorNormalS, ShapeFactorTangentS, NormalS2 )
  use Maths_m
  use Types_m
  implicit none

  type(MeshValues),                   intent(inout) :: newmesh
  type(MeshValues),                   intent(in)    :: mesh, dual
  real(8),                            intent(in)    :: dt, h1, h2, h3
  real(8), dimension(mesh%M, mesh%N), intent(in) :: InvCellVol, CellVol, kappas, CouplingH, CouplingE, invCs, TsPrev, TsOld, &
                                                    NormalW2, NormalE2, NormalN2, NormalS2, &
                                                    ShapeFactorNormalE, ShapeFactorNormalW, &
                                                    ShapeFactorNormalS, ShapeFactorNormalN, &
                                                    ShapeFactorTangentE, ShapeFactorTangentW, &
                                                    ShapeFactorTangentS, ShapeFactorTangentN

  integer :: i, j
  real(8) :: w1, w2, w3, w4

  w1 = h2 * h3 / h1 / (-h3 + h1) / (-h2 + h1)
  w2 = h1 * h3 / (-h2 + h1) / h2 / (-h3 + h2)
  w3 = h1 * h2 / h3 / (h3 ** 2 - h1 * h3 - h2 * h3 + h1 * h2)
  w4 = h2 * h1 * h3 / (h1 * h2 + h1 * h3 + h2 * h3)

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1
    do i=2, mesh%M-1
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

        newmesh%Ts(i,j) = ((M_HALF*( (&
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
                    + M_TWO*(CouplingE(i,j)+CouplingH(i,j)) * CellVol(i,j) &
                ) * InvCellVol(i,j) &
!                     - ((h1 * h2 + h1 * h3 &
!                     + h2 * h3) / h2 / h1 / h3 * Cs(i,j) - h2 * h3 / h1 &
!                     /(-h3 + h1) / (-h2 + h1) * CsOld(i,j) + h1 * h3 / (-h2 &
!                     + h1) / h2 / (-h3 + h2) * CsPrev(i,j) - h1 * h2 / h3 &
!                     / (h3 ** 2 - h1 * h3 - h2 * h3 + h1 * h2) * CsPrev2(i,j))*0d0 & !20150426-Temporal variation of Cs is killed here.
!                     * Ts(i,j)
                    ) *invCs(i,j) &
                    + w1 * mesh%Ts(i,j) - w2 * TsOld(i,j) + w3 * TsPrev(i,j)) * w4
!                    ) / Cs(i,j) - h2 * h3 / h1 / (-h3 + h1) / dt2 * Ts(i,j) &
!                    - h1 * h3 / h2 / dt2 / dt3 * TsOld(i,j) &
!                    + h1 * h2 / (h3 * ( h3 - h1 - h2 ) + h1 * h2) * TsPrev(i,j)) &
!                    / (h1 * h2 + h1 * h3 + h2 * h3) * h2 * h1

                  !TODO: NTD: The mixing should not be done here! 
                  !      TJYD: Which mixing? Do you mean temporal integration?  

    end do
  end do
  !$OMP END DO

end subroutine computeTs

! This routine computes Ue for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeUe( mesh, dt, InvCellVol, kappae,  CouplingE, SourceUe, &
                      invCe, Ue, UeNew, VeX, VeY, CellVol, &
                      ShapeFactorNormalE, ShapeFactorTangentE, ShapeFactorNormalW, ShapeFactorTangentW, &
                      ShapeFactorNormalN, ShapeFactorTangentN, ShapeFactorNormalS, ShapeFactorTangentS, &
                      CellAreaE, CellAreaW, CellAreaN, CellAreaS, &
                      NormalN, NormalS, NormalE, NormalW  )
  use Maths_m
  use Types_m
  implicit none

  type(MeshValues),                   intent(in)    :: mesh
  real(8), dimension(mesh%M, mesh%N), intent(inout) :: UeNew
  real(8),                            intent(in)    :: dt
  real(8), dimension(mesh%M, mesh%N), intent(in)    :: InvCellVol, CellVol, kappae, CouplingE, SourceUe, &
                                                       invCe, Ue, VeX, VeY, &
                                                       ShapeFactorNormalE, ShapeFactorNormalW, &
                                                       ShapeFactorNormalS, ShapeFactorNormalN, &
                                                       ShapeFactorTangentE, ShapeFactorTangentW, &
                                                       ShapeFactorTangentS, ShapeFactorTangentN, &
                                                       CellAreaE, CellAreaW, CellAreaN, CellAreaS
  type(VectorField)                 , intent(in)    :: NormalN, NormalS, NormalE, NormalW

  integer :: i, j

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1
    do i=2, mesh%M-1

! form with bug corrected in derivatives and (OmegaX, OmegaY) drift transport included in finite volumes
!     if(TransportModel < 2) then
         UeNew(i,j) = Ue(i,j) + ((SourceUe(i,j)-CouplingE(i,j))*CellVol(i,j) &
                ! convective term for transport of the energy by the field
                -M_HALF*(((VeX(i+1,j)+VeX(i,j))*NormalE%x(i,j)                   &
                        +(VeY(i+1,j)+VeY(i,j))*NormalE%y(i,j)) * CellAreaE(i,j) &
                +       ((VeX(i,j)+VeX(i-1,j))*NormalW%x(i,j)                   &
                        +(VeY(i,j)+VeY(i-1,j))*NormalW%y(i,j)) * CellAreaW(i,j) &
                +       ((VeX(i,j)+VeX(i,j+1))*NormalN%x(i,j)                   &
                        +(VeY(i,j)+VeY(i,j+1))*NormalN%y(i,j)) * CellAreaN(i,j) &
                +       ((VeX(i,j)+VeX(i,j-1))*NormalS%x(i,j)                   &
                        +(VeY(i,j)+VeY(i,j-1))*NormalS%y(i,j)) * CellAreaS(i,j))&
                ! diffusive term for energy - rewrite correctly
!                 + ( (Ue(i+1,j)-Ue(i,j))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 &
!                 + y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0)*( (kappae(i+1,j)+kappae(i,j))/(Ce(i+1,j)+Ce(i,j)) )*CellAreaE(i,j) &
!                 - 1d0/(x(i,j)**2-2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j) + &
!                 y(i-1,j)**2)**(0.5d0)*((kappae(i-1,j)+kappae(i,j))/(Ce(i-1,j)+Ce(i,j)))*(Ue(i,j)-Ue(i-1,j))*CellAreaW(i,j) &
!                 + 1d0*((kappae(i,j+1)+kappae(i,j))/(Ce(i,j+1)+Ce(i,j)))*(Ue(i,j+1)-Ue(i,j))*CellAreaN(i,j)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j) &
!                 +x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j)+y(i,j)**2)**(0.5d0) &
!                 -1d0*((kappae(i,j)+kappae(i,j-1))/(Ce(i,j)+Ce(i,j-1)))*(Ue(i,j)-Ue(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 &
!                 -2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) &
                )*dt*InvCellVol(i,j)

          UeNew(i,j) = UeNew(i,j) + M_HALF*5d0/3d0*dt*InvCellVol(i,j)*( &
                  NormalE%N(i,j)*ShapeFactorNormalE(i,j)*(kappae(i,j)*invCe(i,j)+kappae(i+1,j)*invCe(i+1,j))*(Ue(i+1,j)-Ue(i,j))  &
                  !
                - ShapeFactorTangentE(i,j)*ShapeFactorNormalE(i,j)*(kappae(i,j) &
                          *invCe(i,j)+kappae(i+1,j)*invCe(i+1,j))*0.25d0*(Ue(i+1,j+1)+Ue(i,j+1)-Ue(i+1,j-1)-Ue(i,j-1)) &
                  !
                + NormalW%N(i,j)*ShapeFactorNormalW(i,j)*(kappae(i-1,j)*invCe(i-1,j)+kappae(i,j)*invCe(i,j))*(Ue(i,j)-Ue(i-1,j)) &
                  !
                - ShapeFactorTangentW(i,j)*ShapeFactorNormalW(i,j)*(kappae(i-1,j) &
                        *invCe(i-1,j)+kappae(i,j)*invCe(i,j))*0.25d0*(Ue(i,j+1)+Ue(i-1,j+1)-Ue(i-1,j-1)-Ue(i,j-1)) &
                  !
                + NormalN%N(i,j)*ShapeFactorNormalN(i,j)*(kappae(i,j+1)*invCe(i,j+1)+kappae(i,j)*invCe(i,j))*(Ue(i,j+1)-Ue(i,j)) &
                  !
                - ShapeFactorTangentN(i,j)*ShapeFactorNormalN(i,j)*(kappae(i,j+1) &
                        *invCe(i,j+1)+kappae(i,j)*invCe(i,j))*0.25d0*(Ue(i+1,j+1)+Ue(i+1,j)-Ue(i-1,j)-Ue(i-1,j+1)) &
                  !
                + NormalS%N(i,j)*ShapeFactorNormalS(i,j)*(kappae(i,j-1)*invCe(i,j-1)+kappae(i,j)*invCe(i,j))*(Ue(i,j)-Ue(i,j-1)) &
                  !
                - ShapeFactorTangentS(i,j)*ShapeFactorNormalS(i,j)*(kappae(i,j-1) &
                        *invCe(i,j-1)+kappae(i,j)*invCe(i,j))*0.25d0*(Ue(i+1,j)+Ue(i+1,j-1)-Ue(i-1,j-1)-Ue(i-1,j)) &
!                   Korfiatis 2007 equation
!                     0.5d0*(NormalEx(i,j)**2+NormalEy(i,j)**2)*CellAreaE(i,j)*(kappae(i,j)*diffusionE(i,j)/Ne(i,j)+kappae(i+1,j)*diffusionE(i+1,j)/Ne(i+1,j))*(Ne(i+1,j)-Ne(i,j)) / (CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistE(i,j) &
!                   - 0.5d0*(CurviEx(i,j)*TangentEx(i,j)+CurviEy(i,j)*TangentEy(i,j))*CellAreaE(i,j)*(kappae(i,j)*diffusionE(i,j)/Ne(i,j)+kappae(i+1,j)*diffusionE(i+1,j)/Ne(i+1,j))*(0.25d0*Ne(i+1,j+1)+0.25d0*Ne(i,j+1)-0.25d0*Ne(i+1,j-1)-0.25d0*Ne(i,j-1))/(CurviEx(i,j)*NormalEx(i,j)+CurviEy(i,j)*NormalEy(i,j))/DistDualE(i,j) &
!                   + 0.5d0*(NormalWx(i,j)**2+NormalWy(i,j)**2)*CellAreaW(i,j)*(kappae(i-1,j)*diffusionE(i-1,j)/Ne(i-1,j)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j)-Ne(i-1,j))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistW(i,j) &
!                   - 0.5d0*(CurviWx(i,j)*TangentWx(i,j)+CurviWy(i,j)*TangentWy(i,j))*CellAreaW(i,j)*(kappae(i-1,j)*diffusionE(i-1,j)/Ne(i-1,j)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i,j+1)+0.25d0*Ne(i-1,j+1)-0.25d0*Ne(i-1,j-1)-0.25d0*Ne(i,j-1))/(CurviWx(i,j)*NormalWx(i,j)+CurviWy(i,j)*NormalWy(i,j))/DistDualW(i,j) &
!                   + 0.5d0*(NormalNx(i,j)**2+NormalNy(i,j)**2)*CellAreaN(i,j)*(kappae(i,j+1)*diffusionE(i,j+1)/Ne(i,j+1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j+1)-Ne(i,j))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistN(i,j) &
!                   - 0.5d0*(CurviNx(i,j)*TangentNx(i,j)+CurviNy(i,j)*TangentNy(i,j))*CellAreaN(i,j)*(kappae(i,j+1)*diffusionE(i,j+1)/Ne(i,j+1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i+1,j+1)+0.25d0*Ne(i+1,j)-0.25d0*Ne(i-1,j)-0.25d0*Ne(i-1,j+1))/(CurviNx(i,j)*NormalNx(i,j)+CurviNy(i,j)*NormalNy(i,j))/DistDualN(i,j) &
!                   + 0.5d0*(NormalSx(i,j)**2+NormalSy(i,j)**2)*CellAreaS(i,j)*(kappae(i,j-1)*diffusionE(i,j-1)/Ne(i,j-1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(Ne(i,j)-Ne(i,j-1))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistS(i,j) &
!                   - 0.5d0*(CurviSx(i,j)*TangentSx(i,j)+CurviSy(i,j)*TangentSy(i,j))*CellAreaS(i,j)*(kappae(i,j-1)*diffusionE(i,j-1)/Ne(i,j-1)+kappae(i,j)*diffusionE(i,j)/Ne(i,j))*(0.25d0*Ne(i+1,j)+0.25d0*Ne(i+1,j-1)-0.25d0*Ne(i-1,j-1)-0.25d0*Ne(i-1,j))/(CurviSx(i,j)*NormalSx(i,j)+CurviSy(i,j)*NormalSy(i,j))/DistDualS(i,j) &
                  )

    end do
  end do
  !$OMP END DO

end subroutine computeUe


! This routine computes Ue for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeUh( mesh, dt, InvCellVol, kappah,  CouplingH, SourceUh, &
                      invCh, Uh, UhNew, VhX, VhY, CellVol, &
                      ShapeFactorNormalE, ShapeFactorTangentE, ShapeFactorNormalW, ShapeFactorTangentW, &
                      ShapeFactorNormalN, ShapeFactorTangentN, ShapeFactorNormalS, ShapeFactorTangentS, &
                      CellAreaE, CellAreaW, CellAreaN, CellAreaS, &
                      NormalN, NormalS, NormalE, NormalW )
  use Maths_m
  use Types_m
  implicit none


  type(MeshValues),                   intent(in)    :: mesh
  real(8), dimension(mesh%M, mesh%N), intent(inout) :: UhNew
  real(8),                            intent(in)    :: dt
  real(8), dimension(mesh%M, mesh%N), intent(in)    :: InvCellVol, CellVol, kappah, CouplingH, SourceUh, &
                                                       invCh, Uh, VhX, VhY,                                 &
                                                       ShapeFactorNormalE, ShapeFactorNormalW,           &
                                                       ShapeFactorNormalS, ShapeFactorNormalN,           &
                                                       ShapeFactorTangentE, ShapeFactorTangentW,         &
                                                       ShapeFactorTangentS, ShapeFactorTangentN,         &
                                                       CellAreaE, CellAreaW, CellAreaN, CellAreaS
  type(VectorField)                 , intent(in)    :: NormalN, NormalS, NormalE, NormalW


  integer :: i, j

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1
    do i=2, mesh%M-1

      UhNew(i,j) = Uh(i,j) + ((SourceUh(i,j)-CouplingH(i,j))*CellVol(i,j) &
                ! convective term for transport of the energy by the field
                -M_HALF*(((VhX(i+1,j)+VhX(i,j))*NormalE%x(i,j)+(VhY(i+1,j)+VhY(i,j))*NormalE%y(i,j)) * CellAreaE(i,j) &
                +( (VhX(i,j)+VhX(i-1,j))*NormalW%x(i,j)+(VhY(i,j)+VhY(i-1,j))*NormalW%y(i,j)) * CellAreaW(i,j) &
                +( (VhX(i,j)+VhX(i,j+1))*NormalN%x(i,j)+(VhY(i,j)+VhY(i,j+1))*NormalN%y(i,j)) * CellAreaN(i,j) &
                +( (VhX(i,j)+VhX(i,j-1))*NormalS%x(i,j)+(VhY(i,j)+VhY(i,j-1))*NormalS%y(i,j)) * CellAreaS(i,j)) &
                ! diffusive term for energy - rewrite correctly
!                 + ( (Uh(i+1,j)-Uh(i,j))/(x(i+1,j)**2-2d0*x(i+1,j)*x(i,j)+x(i,j)**2 &
!                 + y(i+1,j)**2-2d0*y(i+1,j)*y(i,j)+y(i,j)**2)**(0.5d0)*( (kappah(i+1,j)+kappah(i,j))/(Ch(i+1,j)+Ch(i,j)) )*CellAreaE(i,j) &
!                 - 1d0/(x(i,j)**2-2d0*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-2d0*y(i,j)*y(i-1,j) + &
!                 y(i-1,j)**2)**(0.5d0)*((kappah(i-1,j)+kappah(i,j))/(Ch(i-1,j)+Ch(i,j)))*(Uh(i,j)-Uh(i-1,j))*CellAreaW(i,j) &
!                 + 1d0*((kappah(i,j+1)+kappah(i,j))/(Ch(i,j+1)+Ch(i,j)))*(Uh(i,j+1)-Uh(i,j))*CellAreaN(i,j)/(x(i,j+1)**2-2d0*x(i,j+1)*x(i,j) &
!                 +x(i,j)**2+y(i,j+1)**2-2d0*y(i,j+1)*y(i,j)+y(i,j)**2)**(0.5d0) &
!                 -1d0*((kappah(i,j)+kappah(i,j-1))/(Ch(i,j)+Ch(i,j-1)))*(Uh(i,j)-Uh(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 &
!                 -2d0*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-2d0*y(i,j)*y(i,j-1)+y(i,j-1)**2)**(0.5d0)) &
                )*dt*InvCellVol(i,j)

     UhNew(i,j) = UhNew(i,j) +M_HALF*( &
             NormalE%N(i,j)*ShapeFactorNormalS(i,j)*(kappah(i,j)*invCh(i,j)+kappah(i+1,j)*invCh(i+1,j))*(Uh(i+1,j)-Uh(i,j))  &
                  !
              - ShapeFactorTangentE(i,j)*ShapeFactorNormalE(i,j)*(kappah(i,j)*invCh(i,j)     &
                      +kappah(i+1,j)*invCh(i+1,j))*0.25d0*(Uh(i+1,j+1)+Uh(i,j+1)-Uh(i+1,j-1)-Uh(i,j-1))        &
                !
              + NormalW%N(i,j)*ShapeFactorNormalW(i,j)*(kappah(i-1,j)*invCh(i-1,j)+kappah(i,j)*invCh(i,j))*(Uh(i,j)-Uh(i-1,j))  &
                !
              - ShapeFactorTangentW(i,j)*ShapeFactorNormalW(i,j)*(kappah(i-1,j)*invCh(i-1,j) &
                      +kappah(i,j)*invCh(i,j))*0.25d0*(Uh(i,j+1)+Uh(i-1,j+1)-Uh(i-1,j-1)-Uh(i,j-1)) &
                !
              + NormalN%N(i,j)*ShapeFactorNormalN(i,j)*(kappah(i,j+1)*invCh(i,j+1)+kappah(i,j)*invCh(i,j))*(Uh(i,j+1)-Uh(i,j))  &
                !
              - ShapeFactorTangentN(i,j)*ShapeFactorNormalN(i,j)*(kappah(i,j+1)*invCh(i,j+1) &
                      +kappah(i,j)*invCh(i,j))*0.25d0*(Uh(i+1,j+1)+Uh(i+1,j)-Uh(i-1,j)-Uh(i-1,j+1))            &
                !
              + NormalS%N(i,j)*ShapeFactorNormalS(i,j)*(kappah(i,j-1)*invCh(i,j-1)+kappah(i,j)*invCh(i,j))*(Uh(i,j)-Uh(i,j-1))  &
                !
              - ShapeFactorTangentS(i,j)*ShapeFactorNormalS(i,j)*(kappah(i,j-1)*invCh(i,j-1) &
                      +kappah(i,j)*invCh(i,j))*0.25d0*(Uh(i+1,j)+Uh(i+1,j-1)-Uh(i-1,j-1)-Uh(i-1,j)) &
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
                  *dt*InvCellVol(i,j)
    end do
  end do
  !$OMP END DO

end subroutine computeUh


! This routine computes Ue for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeUh_alt( mesh, dt, InvCellVol, kappah,  CouplingH, SourceUh, &
                      Ch, Uh, UhNew, VhX, VhY, CellVol, x, y, &
                      CellAreaE, NormalEx, NormalEy, CellAreaW, NormalWx, NormalWy, &
                      CellAreaN, NormalNx, NormalNy, CellAreaS, NormalSx, NormalSy  )
  use Maths_m
  use Types_m
  implicit none


  type(MeshValues),                   intent(in)    :: mesh
  real(8), dimension(mesh%M, mesh%N), intent(inout) :: UhNew
  real(8),                            intent(in)    :: dt
  real(8), dimension(mesh%M, mesh%N), intent(in) :: InvCellVol, CellVol, kappah, CouplingH, SourceUh, &
                                                    Ch, Uh, VhX, VhY, x, y, &
                                                    CellAreaE, NormalEx, NormalEy, CellAreaW, NormalWx, NormalWy, &
                                                    CellAreaN, NormalNx, NormalNy, CellAreaS, NormalSx, NormalSy

  integer :: i, j

  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1
    do i=2, mesh%M-1

     !TODO: This must be optmised !
        UhNew(i,j) = ((SourceUh(i,j)-CouplingH(i,j))*CellVol(i,j)-M_HALF*( &
               ((VhX(i+1,j)+VhX(i,j))*NormalEx(i,j)+(VhY(i+1,j)+VhY(i,j))*NormalEy(i,j)) * CellAreaE(i,j) &
              +((VhX(i,j)+VhX(i-1,j))*NormalWx(i,j)+(VhY(i,j)+VhY(i-1,j))*NormalWy(i,j)) * CellAreaW(i,j) &
              +((VhX(i,j)+VhX(i,j+1))*NormalNx(i,j)+(VhY(i,j)+VhY(i,j+1))*NormalNy(i,j)) * CellAreaN(i,j) &
              +((VhX(i,j)+VhX(i,j-1))*NormalSx(i,j)+(VhY(i,j)+VhY(i,j-1))*NormalSy(i,j)) * CellAreaS(i,j)) &
              + ( ((kappah(i+1,j)+kappah(i,j))/(Ch(i+1,j)+Ch(i,j))) * (Uh(i+1,j)-Uh(i,j)) &
              /(x(i+1,j)**2-M_TWO*x(i+1,j)*x(i,j)+x(i,j)**2 + y(i+1,j)**2-M_TWO*y(i+1,j)*y(i,j)+y(i,j)**2)**M_HALF &
              *CellAreaE(i,j) &
              - ((kappah(i-1,j)+kappah(i,j))/(Ch(i-1,j)+Ch(i,j))) * (Uh(i,j)-Uh(i-1,j)) &
              /(x(i,j)**2-M_TWO*x(i,j)*x(i-1,j)+x(i-1,j)**2+y(i,j)**2-M_TWO*y(i,j)*y(i-1,j) + y(i-1,j)**2)**M_HALF &
              *CellAreaW(i,j) &
              + ((kappah(i,j+1)+kappah(i,j))/(Ch(i,j+1)+Ch(i,j)))*(Uh(i,j+1)-Uh(i,j)) &
              /(x(i,j+1)**2-M_TWO*x(i,j+1)*x(i,j) + x(i,j)**2+y(i,j+1)**2-M_TWO*y(i,j+1)*y(i,j)+y(i,j)**2)**M_HALF &
              *CellAreaN(i,j) &
              -((kappah(i,j)+kappah(i,j-1))/(Ch(i,j)+Ch(i,j-1)))*(Uh(i,j)-Uh(i,j-1))*CellAreaS(i,j)/(x(i,j)**2 &
              -M_TWO*x(i,j)*x(i,j-1)+x(i,j-1)**2+y(i,j)**2-M_TWO*y(i,j)*y(i,j-1)+y(i,j-1)**2)**M_HALF) &
              )*dt*InvCellVol(i,j)+Uh(i,j)

    end do
  end do
  !$OMP END DO

end subroutine computeUh_alt



! This routine computes Ue for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine computeConvection( mesh, newmesh, UeNew, UhNew, Ue, Uh, invCe, invCh, &
                              FermiTableE, FermiTableH, FermiIndexE, FermiIndexH, &
                              ColFermiThreeHalf, ColFermiHalf, ColFermiMenusHalf, ColFermiEta )
  use Maths_m
  use Types_m
  implicit none


  type(MeshValues),                   intent(in)    :: mesh
  type(MeshValues),                   intent(inout) :: newmesh
  real(8), dimension(mesh%M, mesh%N), intent(in)    :: UeNew, UhNew, Ue, Uh, invCe, invCh
  real(8),           intent(in)    :: FermiTableE(mesh%M,mesh%N)
  real(8),           intent(in)    :: FermiTableH(mesh%M,mesh%N)
  integer(8),        intent(in)    :: FermiIndexE(mesh%M,mesh%N)
  integer(8),        intent(in)    :: FermiIndexH(mesh%M,mesh%N)
  integer(8),        intent(in)    :: ColFermiThreeHalf, ColFermiHalf, ColFermiMenusHalf, ColFermiEta

  integer :: i, j
  !$OMP DO COLLAPSE(2)
  do j=2, mesh%N-1
    do i=2, mesh%M-1
    newmesh%Te(i,j) = mesh%Te(i,j) + ((UeNew(i,j) -  Ue(i,j))-1.5d0*kb*mesh%Te(i,j)*(newmesh%Ne(i,j) - mesh%Ne(i,j)) &
            *FermiTableE(ColFermiThreeHalf,FermiIndexE(i,j))/FermiTableE(ColFermiHalf,FermiIndexE(i,j)) ) * invCe(i,j)
    newmesh%Th(i,j) = mesh%Th(i,j) + ((UhNew(i,j) -  Uh(i,j))-1.5d0*kb*mesh%Th(i,j)*(newmesh%Nh(i,j) - mesh%Nh(i,j)) &
            *FermiTableH(ColFermiThreeHalf,FermiIndexH(i,j))/FermiTableH(ColFermiHalf,FermiIndexH(i,j)) ) * invCh(i,j)
    end do
  end do
  !$OMP END DO

end subroutine computeConvection



! This routine computes Ue for the entire mesh
! We assume that we are in a OMP parallel environement
subroutine applyBoundaryConditions( newmesh, UeNew, UhNew, GradNeX, GradNeY, DriftOn )
  use Maths_m
  use Types_m
  implicit none


  type(MeshValues),                   intent(inout)       :: newmesh
  real(8), dimension(newmesh%M, newmesh%N), intent(inout) :: UeNew, UhNew, GradNeX, GradNeY
  integer(8),                         intent(in)          :: DriftOn

  integer :: i, j

    !BOUNDARY CONDITIONS
    !$OMP PARALLEL  DEFAULT(NONE) SHARED(newmesh, UeNew, UhNew, GradNeX, GradNeY, DriftOn)
    !$OMP DO
    do i=1, newmesh%M !North and South boundaries
      ! finite differences finite difference fashion
      if(DriftOn.eq.0) then
        newmesh%Ne(i,1)=newmesh%Ne(i,2)
        newmesh%Nh(i,1)=newmesh%Nh(i,2)
        newmesh%Ne(i,newmesh%N)=newmesh%Ne(i,newmesh%N-1)
        newmesh%Nh(i,newmesh%N)=newmesh%Nh(i,newmesh%N-1)
      end if

      UeNew(i,1)=UeNew(i,2)
      UhNew(i,1)=UhNew(i,2)
      newmesh%Te(i,1)=newmesh%Te(i,2)
      newmesh%Th(i,1)=newmesh%Th(i,2)
      newmesh%Ts(i,1)=newmesh%Ts(i,2)

      UeNew(i,newmesh%N)=UeNew(i,newmesh%N-1)
      UhNew(i,newmesh%N)=UhNew(i,newmesh%N-1)
      newmesh%Te(i,newmesh%N)=newmesh%Te(i,newmesh%N-1)
      newmesh%Th(i,newmesh%N)=newmesh%Th(i,newmesh%N-1)
      newmesh%Ts(i,newmesh%N)=newmesh%Ts(i,newmesh%N-1)
        ! includes also the corners... WHy are not they written?

    !         potential(i,1)=0d0 !(0d0,0d0)
    !                potential(i,N)=0d0 !(0d0,0d0)
    end do
    !$OMP END DO NOWAIT

    !TODO: Why i is restricted here, and to for  the West/East corners, for the same quantities
    !$OMP DO
    do i=2,newmesh%M-1
      ! boundary condition v.n = 0 on boundaries.
      ! NORTH
      GradNeX(i,newmesh%N) = GradNeX(i,newmesh%N-1)
      GradNeY(i,newmesh%N) = GradNeY(i,newmesh%N-1)
      !
      ! SOUTH
      GradNeX(i,1) = GradNeX(i,2)
      GradNeY(i,1) = GradNeY(i,2)
    end do
    !$OMP END DO NOWAIT

    !$OMP DO
    do j=2, newmesh%N-1 !West and East boundaries

      UeNew(1,j)=UeNew(2,j)
      UhNew(1,j)=UhNew(2,j)
      newmesh%Te(1,j)=newmesh%Te(2,j)
      newmesh%Th(1,j)=newmesh%Th(2,j)
      newmesh%Ts(1,j)=newmesh%Ts(2,j)

      newmesh%Te(newmesh%M,j)=newmesh%Te(newmesh%M-1,j) !Tout
      newmesh%Th(newmesh%M,j)=newmesh%Th(newmesh%M-1,j) !Tout
      newmesh%Ts(newmesh%M,j)=newmesh%Ts(newmesh%M-1,j) ! Tout !cooling by diffusion from outside, TsNew(M-1,j)

    !       potential(1,j)=0d0 !(0d0, 0d0)
    !       potential(M,j)=potential0 !(potential0, 0d0)

      if(DriftOn.eq.0) then
        newmesh%Ne(1,j)=newmesh%Ne(2,j)
        newmesh%Nh(1,j)=newmesh%Nh(2,j)
        ! conditions on the cone base - most important
        newmesh%Ne(newmesh%M,j)=newmesh%Ne(newmesh%M-1,j) !Ne0 replace by dynamic bnd condition with flux equal to the one of cell(M-1,j)
        newmesh%Nh(newmesh%M,j)=newmesh%Nh(newmesh%M-1,j) !Nh0 replace by dynamic bnd condition with flux equal to the one of cell(M-1,j)
      end if

      !TODO: Could we clean up these comments? TJYD (Oct 2, 2016): Were they already commented? 

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
      GradNeX(newmesh%M,j) = GradNeX(newmesh%M-1,j)
      GradNeY(newmesh%M,j) = GradNeY(newmesh%M-1,j)
    end do
    !$OMP END DO

    !$OMP END PARALLEL

end subroutine applyBoundaryConditions
