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
!------------------------------------------------------------------------------

!> Initialize  Parameters
!> @param Params Structure containing input parameters
subroutine InitInputParameter( Params )
    use Types_m
    use Maths_m
    implicit none
    type(InputParameters) :: Params

    !Physical parameters of the simulation
    Params%TimeStep = M_ZERO
    Params%TimeMax  = M_ZERO
    Params%Text     = -M_ONE
    Params%AdaptativeTimeStep = 1
    Params%RestartCalc = 0
    Params%DumpInterval = 10000

    !Numerical parameters of the simulation
    Params%UseMieScattering = -2
    Params%M = -1
    Params%N = -1
    Params%OutputIter = 10000
    
    !Laser parameters for the simulation
    Params%fluence = 0e0
    Params%pulseFWHM = 0e0
    Params%wavelength = 0e0

    !Debug options
    Params%NeOff = 0
    Params%TeOff = 0
    Params%HolesOff = 0
    Params%TsOff = 0
    Params%AugerOff=0
    Params%CouplingDebug=0      !0: e/h - lattice coupling enabled, 1: disabled
    Params%ImpactOff = 0
    Params%CrossDiffusionOff = 0

    Params%DrudeHeating=1       ! free-carrier absorption, 0: Drude heating OFF, 1: enabled (1-epsDrude)
    Params%TransportModel = -1 ! 2: Consider ambipolar diffusion in equations (but careful with boundary conditions)
                                ! 1: consider Tritt particle transport (great expression), but Dumber field is needed !!! -> Poisson !
                                ! 0: only fourier conductivity
                                !-1: diffusion and conductivity OFF
    Params%ConvectionEnergy=0   !0: work with Te, no convection. 1: work with Ue, convection

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
    use Maths_m

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

    if( Params%TransportModel < -1 .or. Params%TransportModel > 2 ) then
      print *, 'Bad value for TransportModel'
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

    if( Params%ImpactOff < 0 .or. Params%ImpactOff > 1 ) then
      print *, 'Bad value for ImpactOff'
      call StopProgram()
    end if

    if( Params%ConvectionEnergy < 0 .or. Params%ConvectionEnergy > 1 ) then
      print *, 'Bad value for ConvectionEnergy'
      call StopProgram()
    end if

    if( Params%CrossDiffusionOff < 0 .or. Params%CrossDiffusionOff > 1 ) then
      print *, 'Bad value for CrossDiffusionOff'
      call StopProgram()
    end if

    if( Params%OutputIter <= 0 ) then
      print *, 'Bad value for OutputIter'
      call StopProgram()
    end if

    if( Params%Text <= M_ZERO ) then
      print *, 'Bad value for the external temperature Text'
      call StopProgram()
    end if

    if( Params%AdaptativeTimeStep < 0 .or. Params%AdaptativeTimeStep > 1 ) then
      print *, 'Bad value for AdaptativeTimeStep'
      call StopProgram()
    end if

    if( Params%RestartCalc < 0 .or. Params%RestartCalc > 1 ) then
      print *, 'Bad value for RestartCalc'
      call StopProgram()
    end if

    if( Params%DumpInterval < 0 ) then
      print *, 'Bad value for DumpInterval'
      call StopProgram()
    end if
    
    if( Params%Fluence < 0d0 ) then
      print *, 'Laser fluence cannot be negative'
      call StopProgram()
    end if
    
    if( (Params%Wavelength .le. 0e0) .or. ( (Params%Wavelength  .ne. 343e-9) .and. (Params%Wavelength  .ne. 515e-9) &
        .or. (Params%Wavelength  .ne. 800e-9 ) .and. ( Params%Wavelength  .ne. 1030e-9 ) ) ) then
      print *, 'Laser wavelength is for now limited to 343, 515, 800, and 1030 nm.'
      call StopProgram()
    end if
    
    if( Params%PulseFWHM .le. 0e0 ) then
      print *, 'Pulse duration cannot be null or negative. Optimal value is close to 50e-15 s = 50 fs.'
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

              !TransportModel
              if( id .equals. 'TransportModel' ) then
                call ParseInt( line, Params%TransportModel )
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

              !ImpactOff
              if( id .equals. 'ImpactOff' ) then
                call ParseInt( line, Params%ImpactOff )
                goto 999
              endif

              !OutputIter
              if( id .equals. 'OutputIter' ) then
                call ParseInt( line, Params%OutputIter )
                goto 999
              endif

              !ConvectionEnergy
              if( id .equals. 'ConvectionEnergy' ) then
                call ParseInt( line, Params%ConvectionEnergy )
                goto 999
              endif

              !CrossDiffusionOff
              if( id .equals. 'CrossDiffusionOff' ) then
                call ParseInt( line, Params%CrossDiffusionOff )
                goto 999
              endif

              !Text
              if( id .equals. 'Text' ) then
                call ParseDouble( line, Params%Text )
                goto 999
              endif

              !AdaptativeTimeStep
              if( id .equals. 'AdaptativeTimeStep' ) then
                call ParseInt( line, Params%AdaptativeTimeStep )
                goto 999
              endif

              !RestartCalc
              if( id .equals. 'RestartCalc' ) then
                call ParseInt( line, Params%RestartCalc )
                goto 999
              endif

              !DumpInterval
              if( id .equals. 'DumpInterval' ) then
                call ParseInt( line, Params%DumpInterval )
                goto 999
              endif
              
              !Wavelength
              if( id .equals. 'Wavelength' ) then
                call ParseDouble( line, Params%Wavelength )
                goto 999
              endif
              
              !pulseFWHM
              if( id .equals. 'PulseFWHM' ) then
                call ParseDouble( line, Params%PulseFWHM )
                goto 999
              endif
              
              !Fluence
              if( id .equals. 'Fluence' ) then
                call ParseDouble( line, Params%Fluence )
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

    call PrintInputParameters(Params)

end subroutine LoadInputParameters


subroutine PrintInputParameters(Params)
  use Types_m
  implicit none
  type(InputParameters), intent(in) :: Params

  integer :: unit, ios

  unit = 1

  OPEN( UNIT=unit, FILE='parser.log', STATUS="unknown",access='sequential', ACTION="write", IOSTAT=ios )
  if ( ios /= 0 ) then ! Probleme ea l'ouverture
        print *, 'Error opening parser.log'
        call StopProgram()
  else
    !Physical parameters of the simulation
    write(unit, '(a,e12.5)') 'TimeStep = ', Params%TimeStep
    write(unit, '(a,e12.5)') 'TimeMax = ', Params%TimeMax
    write(unit, '(a,e12.5)') 'Text = ', Params%Text
    write(unit, '(a,i2)') 'AdaptativeTimeStep = ', Params%AdaptativeTimeStep
    write(unit, '(a,i2)') 'RestartCalc = ', Params%RestartCalc
    write(unit, '(a,i2)') 'DumpInterval = ', Params%DumpInterval


    !Numerical parameters of the simulation
    write(unit, '(a,i2)') 'UseMieScattering = ', Params%UseMieScattering
    write(unit, '(a,i5)') 'M = ', Params%M
    write(unit, '(a,i5)') 'N = ', Params%N
    write(unit, '(a,i5)') 'OutputIter = ', Params%OutputIter

    !Laser parameters
    write(unit, '(a,e12.5)') 'pulse duration (FWHM) = ', Params%PulseFWHM
    write(unit, '(a,e12.5)') 'Laser wavelength (m) = ', Params%Wavelength
    write(unit, '(a,e12.5)') 'Laser Fluence (J/m2) = ', Params%Fluence
    
    !Debug options
    write(unit, '(a,i2)') 'NeOff = ', Params%NeOff
    write(unit, '(a,i2)') 'TeOff = ', Params%TeOff
    write(unit, '(a,i2)') 'HolesOff = ', Params%HolesOff
    write(unit, '(a,i2)') 'TsOff = ', Params%TsOff
    write(unit, '(a,i2)') 'AugerOff = ', Params%AugerOff
    write(unit, '(a,i2)') 'CouplingDebug = ', Params%CouplingDebug
    write(unit, '(a,i2)') 'ImpactOff = ', Params%ImpactOff
    write(unit, '(a,i2)') 'CrossDiffusionOff = ', Params%CrossDiffusionOff
    write(unit, '(a,i2)') 'DrudeHeating = ', Params%DrudeHeating
    write(unit, '(a,i2)') 'TransportModel = ', Params%TransportModel
    write(unit, '(a,i2)') 'ConvectionEnergy = ', Params%ConvectionEnergy

    !Mie scattering
    write(unit, '(a,e12.5)') 'phiMie0 = ', Params%phiMie0
    write(unit, '(a,i2)') 'PolarizationSource = ', Params%PolarizationSource
  end if

  CLOSE( UNIT=unit )

end subroutine PrintInputParameters

