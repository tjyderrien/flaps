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
!> @file mesh.f90
!
! MODULE: Types
!
! DESCRIPTION:
!> This module defines various types used in the code
!------------------------------------------------------------------------------

module Types_m

  !> Parameters for defining the mesh
  type MeshValues
    real(8), allocatable, dimension(:) :: Te  !> electron temperature
    real(8), allocatable, dimension(:) :: Th  !> hole temperature
    real(8), allocatable, dimension(:) :: Ts  !> lattice temperature
    real(8), allocatable, dimension(:) :: Ne  !> electron density
    real(8), allocatable, dimension(:) :: Nh  !> hole density
    integer(8), allocatable, dimension(:,:) :: Map !> A mapping array

    integer :: M,N !> Mesh dimension
   end type MeshValues

  !> A vector field structure
  type VectorField
    real(8), allocatable, dimension(:,:) :: x, y !> x, y components of a vector
    real(8), allocatable, dimension(:,:) :: N    !> Norm of this vector
  end type VectorField

  !> This structure contains to the information for one output file
  type OutputData
    integer :: unit       !> The unit for input/output
  end type OutputData

  !> This structure contains input parameters, to be read from the input file
  type InputParameters
    integer :: UseMieScattering  ! 1: Enable Mie scattering analytic formula, 0: badly fitted FDTD input, -1: constant intensity
    integer :: NeOff, TeOff, HolesOff, TsOff
    integer :: M,N !> Mesh dimension
    real(8) :: TimeStep, TimeMax
    real(8) :: phiMie0
    integer :: PolarizationSource             ! Value of the Mie angle that will be distributed on various processors
    integer :: DrudeHeating          ! free-carrier absorption, 0: Drude heating OFF, 1: enabled (1-epsDrude)
    integer :: TransportModel        ! 2: Consider ambipolar diffusion in equations (but careful with boundary conditions)
                                     ! 1: consider Tritt particle transport (great expression), but Dumber field is needed !!! -> Poisson !
                                     ! 0: only fourier conductivity
                                     !-1: diffusion and conductivity OFF
    integer :: CouplingDebug         !0: e/h - lattice coupling enabled, 1: disabled
    integer :: AugerOff              !TODO: Comment this value
    integer :: OutputIter            !TODO: Comment this value
    integer :: ImpactOff             !TODO: Comment this value
    integer :: ConvectionEnergy      !0: work with Te, no convection. 1: work with Ue, convection
    integer :: CrossDiffusionOff     !TODO: Comment this value
    real(8) :: Text                  !external temperature (K)
    integer :: AdaptativeTimeStep
    integer :: RestartCalc
    integer :: DumpInterval
    real(8) :: Fluence
    real(8) :: PulseFWHM
    real(8) :: Wavelength
  end type

end module Types_m

