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
!> @file output.f90
!
! DESCRIPTION:
!> @brief All the routines for the output
!
!> @author
!> Nicolas Tancogne-Dejean
!
!> @date
!> 28 Jul 2016 - Initial Version
!------------------------------------------------------------------------------

module Output_m
  use Types_m

  implicit none

  !TODO: Add a description for each of these files
  type(OutputData) :: LaplaceConvergence ! Useful for debugging
  type(OutputData) :: LaplaceMatrix      ! Useful for debugging
  type(OutputData) :: TimeBottom         
  type(OutputData) :: TimeUp
  type(OutputData) :: TimeApex
  type(OutputData) :: Error
  type(OutputData) :: Parameters
  type(OutputData) :: Depth !TODO: This name is not self-documenting
  type(OutputData) :: TimeMax
  type(OutputData) :: MeshInfo
  type(OutputData) :: MeshVessel
  type(OutputData) :: DepthVessel
  type(OutputData) :: DualDepth          ! Useful for debugging
  type(OutputData) :: Field
  type(OutputData) :: EnergyConservation
  type(OutputData) :: Temperature
  type(OutputData) :: Density

  contains

  subroutine InitOutputs()
    implicit none


    !TODO: status unknow is not clean. Better to ensure new files
    !TODO: We should maybe check is file exists first

    !TODO: we have to hav here the mkdir output_dir

    ! opening files
    LaplaceConvergence%unit = 90
    open(LaplaceConvergence%unit,FILE='output/LaplaceConvergence.dat', access='sequential',status='unknown')

    LaplaceMatrix%unit = 91
    open(LaplaceMatrix%unit, FILE='output/LaplaceMatrix.dat',access='sequential',status='unknown')

    TimeBottom%unit = 92
    open(TimeBottom%unit, FILE='output/TimeBottom.dat', access='sequential', status='unknown')         ! format 882

    TimeUp%unit = 93
    open(TimeUp%unit, FILE='output/TimeUp.dat', access='sequential', status='unknown')         ! format 883

    TimeApex%unit = 94
    open(TimeApex%unit, FILE='output/TimeApex.dat', access='sequential', status='unknown')         ! format 884

    Error%unit = 95
    open(Error%unit,FILE='output/error.dat', access='sequential', status='unknown')

    Parameters%unit = 96
    open(Parameters%unit,FILE='output/parameters.dat', access='sequential', status='unknown')

    Depth%unit = 97
    open(Depth%unit,FILE='output/Depth.dat',access='sequential',status='unknown')                ! format 887

    TimeMax%unit = 98
    open(TimeMax%unit,FILE='output/TimeMax.dat',access='sequential',status='unknown')                 ! format 888

    MeshInfo%unit = 99
    open(MeshInfo%unit,FILE='output/mesh.dat',access='sequential',status='unknown')                ! format 885, 8852

    MeshVessel%unit = 100
    open(MeshVessel%unit,FILE='output/meshVessel.dat',access='sequential',status='unknown')

    DepthVessel%unit = 101
    open(DepthVessel%unit,FILE='output/DepthVessel.dat',access='sequential',status='unknown')         ! format 889

    DualDepth%unit = 103
    open(DualDepth%unit,FILE='output/DualDepth.dat', access='sequential', status='unknown') ! format 890

    Field%unit = 104
    open(Field%unit,FILE='output/Field.dat', access='sequential', status='unknown') ! format 891

    EnergyConservation%unit = 105
    open(EnergyConservation%unit, FILE='output/EnergyConservation.dat', access='sequential', status='unknown') !format 892
    
    Temperature%unit = 106
    open(Temperature%unit, FILE='output/Temperature.dat', access='sequential', status='unknown') !format 893
    
    Temperature%unit = 107
    open(Density%unit, FILE='output/Density.dat', access='sequential', status='unknown') !format 894

  end subroutine InitOutputs


end module Output_m



