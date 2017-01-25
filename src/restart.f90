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
!> @file Restart.f90
!
! DESCRIPTION:
!> @brief Restart module to compute calculations on long times.
!------------------------------------------------------------------------------

module Restart_m
  use Types_m

  implicit none

  private

  public ::                &
    Restart_load,          &
    Restart_dump

  contains
  !------------------------------------------------------------------
  !> Load previously stored data
  !------------------------------------------------------------------
  subroutine Restart_load( mesh, UeNew, UhNew, TsOld, Ce, Ch, CsPrev, CsOld, Cs, t, iter, &
       IntensityEnergy, ElectronEnergy, ElectronKineticEnergy, ElectronPotentialEnergy, HoleEnergy, LatticeEnergy )
    type(MeshValues),                 intent(inout) :: mesh
    real(8), dimension(mesh%M, mesh%N), intent(out) :: UeNew, UhNew, TsOld, Ce, Ch, Cs, CsPrev, CsOld
    real(8),                            intent(out) :: t
    integer(8),                         intent(out) :: iter
    real(8),                            intent(out) :: IntensityEnergy, ElectronEnergy, ElectronKineticEnergy, &
                                                       ElectronPotentialEnergy, HoleEnergy, LatticeEnergy

    integer :: unit, ios, N, M, i, j
    logical :: exists

    print *, 'Start reading restart information'

    inquire( FILE='restart.flp', EXIST=exists )
    if( exists .eqv. .false. ) then
        print *,'The restart file does not exist.'
        call StopProgram()
    endif

    OPEN( UNIT=unit, FILE='restart.flp',  POSITION="rewind", &
             FORM="formatted", access='sequential', ACTION="read", IOSTAT=ios )
    if ( ios /= 0 ) then ! Probleme ea l'ouverture
        print *, 'Error opening restart.flp file'
        call StopProgram()
    else
      read(unit,*, end=999) N, M
      if(N /= mesh%N .or. M /= mesh%M) then
        print *, 'Restarting file restart.flp does not correspond to the current simulation.'
        call StopProgram()
      end if

      read(unit,*, end=999) t, iter

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) mesh%Te(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) mesh%Th(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) mesh%Ts(i,j)
        end do
      end do
      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) mesh%Ne(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) mesh%Nh(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) UeNew(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) UhNew(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) TsOld(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) Ce(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) Ch(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) CsPrev(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) CsOld(i,j)
        end do
      end do

       do j=1, mesh%N
        do i=1, mesh%M
          read(unit,*, end=999) Cs(i,j)
        end do
      end do

     read(unit,*,end=999) IntensityEnergy, ElectronEnergy, ElectronKineticEnergy, ElectronPotentialEnergy, HoleEnergy, LatticeEnergy

    close(unit)

    end if

    print *, 'End reading restart information'

    return

 999 print *, 'Error while reading restart.flp. The file probably corrupted.'
     call StopProgram()

  end subroutine Restart_load

    !------------------------------------------------------------------
  !> Load previously stored data
  !------------------------------------------------------------------
  subroutine Restart_dump(mesh, UeNew, UhNew, TsOld, Ce, Ch, CsPrev, CsOld, Cs, t, iter, &
      IntensityEnergy, ElectronEnergy, ElectronKineticEnergy, ElectronPotentialEnergy, HoleEnergy, LatticeEnergy )
    type(MeshValues),                    intent(in) :: mesh
    real(8), dimension(mesh%M, mesh%N), intent(in) :: UeNew, UhNew, TsOld, Ce, Ch, Cs, CsPrev, CsOld
    real(8),                            intent(in) :: t
    integer(8),                         intent(in) :: iter
    real(8),                            intent(in) :: IntensityEnergy, ElectronEnergy, ElectronKineticEnergy, &
                                                      ElectronPotentialEnergy, HoleEnergy, LatticeEnergy

    integer :: unit, ios, i, j

    print *, 'Start writing restart information'


    OPEN( UNIT=unit, FILE='restart.flp', POSITION="rewind", FORM="formatted", &
             access='sequential', ACTION="write", IOSTAT=ios )
    if ( ios /= 0 ) then ! Probleme ea l'ouverture
        print *, 'Error opening restart.flp file'
        call StopProgram()
    else
      write(unit,*) mesh%N, mesh%M
      write(unit,*) t, iter

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) mesh%Te(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) mesh%Th(i,j)
        end do
      end do

       do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) mesh%Ts(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) mesh%Ne(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) mesh%Nh(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) UeNew(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) UhNew(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) TsOld(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) Ce(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) Ch(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) CsPrev(i,j)
        end do
      end do

      do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) CsOld(i,j)
        end do
      end do

       do j=1, mesh%N
        do i=1, mesh%M
          write(unit,*) Cs(i,j)
        end do
      end do

      write(unit,*) IntensityEnergy, ElectronEnergy, ElectronKineticEnergy, &
                        ElectronPotentialEnergy, HoleEnergy, LatticeEnergy


      write(unit,*) 'End of file'

      close(unit)
    end if

    print *, 'End writing restart information'

  end subroutine Restart_dump

end module Restart_m
