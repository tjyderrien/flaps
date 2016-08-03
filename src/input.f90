!> Initialize  Parameters
!> @param Params Structure containing input parameters
subroutine InitInputParameter( Params )
    use Types_m
    implicit none
    type(InputParameters) :: Params

    Params%UseMieScattering = -2
    Params%M = -1
    Params%N = -1

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

