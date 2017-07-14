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
!------------------------------------------------------------------------------

module Output_m
  use Types_m

  implicit none

  !TODO: Add a description for each of these files
  type(OutputData) :: LaplaceConvergence ! Useful for debugging
  type(OutputData) :: LaplaceMatrix      ! Useful for debugging
  type(OutputData) :: TimeBottom         ! Probing the variables at one given point. 
  type(OutputData) :: TimeUp             ! Probing the variables at one given point. 
  type(OutputData) :: TimeApex		 ! Probing the variables at one given point. 
  type(OutputData) :: ErrorFile
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

  subroutine InitOutputs(RestartCalc)
    integer :: RestartCalc


    !We create the directory for the outputs
    call system("mkdir output")

    ! opening files
    LaplaceConvergence%unit = 90
    LaplaceMatrix%unit = 91
    TimeBottom%unit = 92
    TimeUp%unit = 93
    TimeApex%unit = 94
    ErrorFile%unit = 95
    Parameters%unit = 96
    Depth%unit = 97
    TimeMax%unit = 98
    MeshInfo%unit = 99
    MeshVessel%unit = 100
    DepthVessel%unit = 101
    DualDepth%unit = 103
    Field%unit = 104
    EnergyConservation%unit = 105
    Temperature%unit = 106
    Temperature%unit = 107
  end subroutine InitOutputs

  subroutine CloseOutputs()

    close(LaplaceConvergence%unit)
    close(LaplaceMatrix%unit)
    close(TimeBottom%unit)
    close(TimeUp%unit)
    close(TimeApex%unit)
    close(ErrorFile%unit)
    close(Parameters%unit)
    close(Depth%unit)
    close(TimeMax%unit)
    close(MeshInfo%unit)
    close(MeshVessel%unit)
    close(DepthVessel%unit)
    close(DualDepth%unit)
    close(Field%unit)
    close(EnergyConservation%unit)
    close(Temperature%unit)
    close(Density%unit)
  end subroutine CloseOutputs

  subroutine output_open(unit, filename, append)
    integer,          intent(in) :: unit
    character(len=*), intent(in) :: filename
    logical,          intent(in) :: append

    logical :: exist

    inquire(file=trim(filename), exist=exist)
    if(exist) then
      if(append) then
        open(unit, file=trim(filename), status="old", position="append", action="write", access='sequential')
      else
        open(unit, file=trim(filename), status="old", position="rewind", action="write", access='sequential')
      end if
    else
      open(unit, file=trim(filename), status="new", action="write",access='sequential')
    end if

  end subroutine

end module Output_m



