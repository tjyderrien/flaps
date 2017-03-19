!! Copyright (C) 2017  N. Tancogne-Dejean
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
!> @file fermi.f90
!
! DESCRIPTION:
!> @brief This files contains routines related to the Fermi integrals
!------------------------------------------------------------------------------


module fermi_m
  implicit none

  private

  !Indices of the column in the table of Fermi integrals
  integer, public, parameter :: FERMI_NENC = 2;
  integer, public, parameter :: FERMI_ETA = 3;
  integer, public, parameter :: FERMI_0 = 4;
  integer, public, parameter :: FERMI_1 = 5;
  integer, public, parameter :: FERMI_2 = 6;
  integer, public, parameter :: FERMI_HALF = 7;
  integer, public, parameter :: FERMI_THREE_HALF = 8;
  integer, public, parameter :: FERMI_MINUS_HALF = 9;

end module fermi_m
