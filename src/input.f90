!! Copyright (C) 2016 N. Tancogne-Dejean
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
!> @file input.f90
!
! DESCRIPTION:
!> @brief Provide the code with a parser for the input file.
!
!> @author
!> N. Tanconge-Dejean
!
!> @date
!> 26 Aug 2016 - Initial Version
!------------------------------------------------------------------------------

!> Initialize  Parameters
!> @param Params Structure containing input parameters
subroutine InitInputParameter( Params )
    use Types_m
    implicit none
    type(InputParameters) :: Params

    !Paramters of the simulation
    Params%UseMieScattering = -2
    Params%M = -1
    Params%N = -1
    Params%TimeStep = 0.d0
    Params%TimeMax  = 0.d0

    !Debug options
    Params%NeOff = 0
    Params%TeOff = 0
    Params%HolesOff = 0
    Params%TsOff = 0
    Params%AugerOff=0

    Params%DrudeHeating=1       ! free-carrier absorption, 0: Drude heating OFF, 1: enabled (1-epsDrude)
    Params%ConductivityFix = -1 ! 2: Consider ambipolar diffusion in equations (but careful with boundary conditions)
                                ! 1: consider Tritt particle transport (great expression), but Dumber field is needed !!! -> Poisson !
                                ! 0: only fourier conductivity
                                !-1: diffusion and conductivity OFF
    Params%CouplingDebug=0      !0: e/h - lattice coupling enabled, 1: disabled


    !Mie scattering
    Params%phiMie0 = 0
    Params%PolarizationSource = 0



end subroutine InitInputParameter

!> Destroy Input Parameters
!> @param Params Structure containing input parameters
subroutine ReleaseInputParameters( Params )
    use Types_m
    implicit none
    type(InputParameters) :: Params

end subroutine ReleaseInputParameters

!> Check validity of parameters and compute usefull values for the program
!> @param Params Structure containing input parameters
subroutine CheckValidityInputParameters( Params )
    use Types_m

    implicit none
    type(InputParameters) :: Params

    if( abs(Params%UseMieScattering) > 1) then
        print *, 'Bad value for UseMieScattering'
        call StopProgram()
    end if

    if( Params%M <= 3) then
        print *, 'Bad value for M'
        call StopProgram()
    end if

    if( Params%N <= 3) then
        print *, 'Bad value for N'
        call StopProgram()
    end if

    if( Params%NeOff < 0 .or. Params%NeOff > 1 ) then
      print *, 'Bad value for NeOff'
      call StopProgram()
    end if

    if( Params%TeOff < 0 .or. Params%TeOff > 1 ) then
      print *, 'Bad value for TeOff'
      call StopProgram()
    end if

    if( Params%HolesOff < 0 .or. Params%HolesOff > 1 ) then
      print *, 'Bad value for HolesOff'
      call StopProgram()
    end if

    if( Params%TsOff < 0 .or. Params%TsOff > 1 ) then
      print *, 'Bad value for TeOff'
      call StopProgram()
    end if

    if( Params%TimeStep <= 0.0d0 .or. Params%TimeMax < 0.0d0 ) then
      print *, 'Bad value for TimeStep or TimeMax'
      call StopProgram()
    end if

    if( Params%phiMie0 < 0.0d0 .or. Params%phiMie0 > 360.0d0 ) then
      print *, 'Bad value for phiMie0'
      call StopProgram()
    end if

    if( Params%PolarizationSource < 0.or. Params%PolarizationSource > 1 ) then
      print *, 'Bad value for PolarizationSource'
      call StopProgram()
    end if

    if( Params%DrudeHeating < 0 .or. Params%DrudeHeating > 1 ) then
      print *, 'Bad value for DrudeHeating'
      call StopProgram()
    end if

    if( Params%ConductivityFix < -1 .or. Params%ConductivityFix > 2 ) then
      print *, 'Bad value for ConductivityFix'
      call StopProgram()
    end if

    if( Params%CouplingDebug < 0 .or. Params%CouplingDebug > 1 ) then
      print *, 'Bad value for CouplingDebug'
      call StopProgram()
    end if

    if( Params%AugerOff < 0 .or. Params%AugerOff > 1 ) then
      print *, 'Bad value for AugerOff'
      call StopProgram()
    end if

end subroutine CheckValidityInputParameters

subroutine LoadInputParameters( filename, Params )
    use Types_m
    use String_m

    implicit none
    type(InputParameters) :: Params
    CHARACTER( LEN=* ), intent(in) :: filename

    CHARACTER( LEN=512 ) :: line, id
    LOGICAL :: IsOK, ReadLine
    integer ios, k, unit

    unit = 1

      OPEN( UNIT=unit, FILE=filename, FORM="formatted", &
        ACCESS="sequential", STATUS="old", ACTION="read", &
        POSITION="rewind", IOSTAT=ios )
    if ( ios /= 0 ) then ! Probleme ea l'ouverture
        print *, 'Error opening Input Parameters file'
        call StopProgram()
    else
         IsOK=ReadLine( unit, line )

         !TODO: change to select case, if working

         do while( IsOK .eqv. .true. )
              !Si la ligne est vide, on la passe
              if( len_trim(line) == 0 ) then
                goto 999
              endif

              !La ligne commence par un commentaire
              k = scan(line, '#!')
              if( k == 1 ) then
                goto 999
              endif

              k = index(line, ' ')
              if(k <= 1 ) then !En cas de problème, on ignore la ligne
                goto 999
              endif
              id = line(1:k-1) !Nom du paramètre

              !UseMieScattering
              if( id .equals. 'UseMieScattering' ) then
                call ParseInt( line, Params%UseMieScattering )
                goto 999
              endif

              !M
              if( id .equals. 'M' ) then
                call ParseInt( line, Params%M )
                goto 999
              endif

              !N
              if( id .equals. 'N' ) then
                call ParseInt( line, Params%N )
                goto 999
              endif

              !NeOff
              if( id .equals. 'NeOff' ) then
                call ParseInt( line, Params%NeOff )
                goto 999
              endif

              !TeOff
              if( id .equals. 'TeOff' ) then
                call ParseInt( line, Params%TeOff )
                goto 999
              endif

              !HolesOff
              if( id .equals. 'HolesOff' ) then
                call ParseInt( line, Params%HolesOff )
                goto 999
              endif

              !TsOff
              if( id .equals. 'TsOff' ) then
                call ParseInt( line, Params%TsOff )
                goto 999
              endif

              !TimeStep
              if( id .equals. 'TimeStep' ) then
                call ParseDouble( line, Params%TimeStep )
                goto 999
              endif

              !TimeDouble
              if( id .equals. 'TimeMax' ) then
                call ParseDouble( line, Params%TimeMax )
                goto 999
              endif

              !TsOff
              if( id .equals. 'PolarizationSource' ) then
                call ParseInt( line, Params%PolarizationSource )
                goto 999
              endif

              !phiMie0
              if( id .equals. 'phiMie0' ) then
                call ParseDouble( line, Params%phiMie0 )
                goto 999
              endif

              !DrudeHeating
              if( id .equals. 'DrudeHeating' ) then
                call ParseInt( line, Params%DrudeHeating )
                goto 999
              endif

              !ConductivityFix
              if( id .equals. 'ConductivityFix' ) then
                call ParseInt( line, Params%ConductivityFix )
                goto 999
              endif

              !CouplingDebug
              if( id .equals. 'CouplingDebug' ) then
                call ParseInt( line, Params%CouplingDebug )
                goto 999
              endif

              !AugerOff
              if( id .equals. 'AugerOff' ) then
                call ParseInt( line, Params%AugerOff )
                goto 999
              endif

! Some examples

!              !Type de la base d'orbitales
!              if( id .equals. 'sp3s*d5' ) then
!                POrb%nopa = 10
!                IMParams%ndim = 5
!                POrb%ndim = 5
!                goto 999
!              endif

!              !Es
!              if( id .equals. 'Es' ) then
!                call ParseTabDouble( line, POrb%Es, 2, 2, .false. )
!                goto 999
!              endif


              print *,'Error : unexpected line.'
              print *,'line : ',trim(line)
              call StopProgram()

999      IsOK=ReadLine( unit, line )
         end do
    endif
    CLOSE( UNIT=unit )

end subroutine LoadInputParameters

