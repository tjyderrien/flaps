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
!> @file Utils.f90
!
! DESCRIPTION:
!> @brief Classe qui fournit les routines la gestion des fichiers ( ouverture et recherche de noms )
!------------------------------------------------------------------------------


!> Méthode qui génère un nom de fichier à partir d'un nom
!> initial en recherchant le premier fichier non deja existant
!> Si le fichier n'existe pas, on lui ajoute A, B, .., AA,
!> AB, ... jusqu'a touver un fichier non existant ou atteindre ZZZ
!
!> @param filename Nom du fichier a ouvrir
!> @param length Taille du nom du fichier
subroutine GetValidName( filename, length )
    implicit none
    INTEGER, INTENT(IN) :: length
    INTEGER :: pos
    INTEGER, DIMENSION(3) :: counter
    CHARACTER( LEN=(length + 3) ) :: filename
    CHARACTER(LEN=26), PARAMETER :: letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    LOGICAL :: exists

    !On teste si le fichier existe
    INQUIRE( FILE=trim(filename), EXIST=exists )
    if( exists .eqv. .false. ) then
        print *,'Find valid file name ', trim(filename)
        return
    endif
    pos = 1
    counter(:) = 1
    !Si le fichier n'existe pas, on lui ajoute A, B, .., AA,
    ! AB, ... jusqu'a touver un fichier non existant ou atteindre ZZZ
    do while( exists .eqv. .true. .and. pos < 4)
        print *, 'The file ', trim(filename), ' already exists.'
        filename( (length + 1):(length + 1) ) = letters(counter(1):counter(1))
        if( pos >= 2 ) then
            filename( (length + 2):(length + 2) ) = letters(counter(2):counter(2))
        endif
        if( pos == 3 ) then
            filename( (length + 3):(length + 3) ) = letters(counter(3):counter(3))
        endif

        counter(pos) = counter(pos)+1

        if( counter(3) == 27 ) then
            if( counter(2) == 26 .and. counter(3) == 26 ) then
                pos = pos+1
            endif
            counter(2) = counter(2)+1
            counter(3) = 1
        endif

        if( counter(2) == 27 ) then
            if( pos == 2 .and. counter(1) == 26 ) then
                pos = pos+1
            endif
            counter(1) = counter(1)+1
            counter(2) = 1
        endif

        if( counter(1) == 27 ) then
            if( pos == 1 ) then
                pos = pos +1
            endif
            counter(1) = 1
        endif

        INQUIRE( FILE=filename, EXIST=exists )
    end do

    if( pos == 4 ) then
        print *,'Erreur : pas de nom de fichier disponible'
        call StopProgram()
    endif
    print *,'Find valid file name ', trim(filename)
end subroutine GetValidName

!> Méthode qui ouvre le fichier filename pour l'unité unit en vérifiant que l'ouverture est correcte
!> Le fichier est ouvert en mode écriture
!
!> @param filename Nom du fichier a ouvrir
!> @param unit Unité associée à l'ouverture du fichier
subroutine OpenFileForWritting( filename, unit )
   implicit none
   CHARACTER( LEN=512 ) :: filename
   integer :: ios, unit


   call GetValidName( filename, LEN( trim(filename)) )
   OPEN( UNIT=unit, FILE=filename, FORM="formatted", ACCESS="sequential", &
   STATUS="new", ACTION="write", POSITION="rewind", IOSTAT=ios )
   if ( ios /= 0 ) then
       print *, 'Error opening file'
   call StopProgram()
   endif

end subroutine OpenFileForWritting

!> Méthode qui ouvre le fichier filename pour l'unité unit en vérifiant que l'ouverture est correcte
!> Le fichier est ouvert en mode lecture
!
!> @param filename Nom du fichier a ouvrir
!> @param unit Unité associée à l'ouverture du fichier
subroutine OpenFileForReading( filename, unit )
   implicit none
   CHARACTER( LEN=512 ) :: filename
   integer :: ios, unit

    OPEN( UNIT=unit, FILE=filename, FORM="formatted", ACCESS="sequential", &
    STATUS="old", ACTION="read", POSITION="rewind", IOSTAT=ios )
    if ( ios /= 0 ) then
        print *, 'Error opening file'
        call StopProgram()
    endif
end subroutine OpenFileForReading

!> Permet d'afficher la date et l'heure de début d'exécution du programme
subroutine PrintTime( )
    implicit none
    character (len=8) :: date
    character (len=10) :: time
    character (len=5) :: zone
    integer, dimension(8) :: values

    call date_and_time ( date, time, zone, values )
    print '(A20,1X,2(I2,A1),I4,A4,I2,A1,I2)',&
    'Program started the', values(3), '/', values(2), '/', values(1), &
    ' at ', values(5),':', values(6)
    print *, ' '
end subroutine PrintTime


!> Méthode qui écrit le nom de l'ordinateur de manière formaté
!> Le fichier est supposé déjà ouvert
!
!> @param unit Unité associée à l'ouverture du fichier
!> @param computername Nom de l'ordinateur à imprimer
subroutine PrintComputerName( unit, computername )
   implicit none
   CHARACTER( LEN=512 ) :: computername, format
   integer :: unit

   write ( format, '(A11,I2,A1)' ) '(A10, 1X, A', LEN_TRIM(computername), ')'

   write( unit, format) '#Computer ', TRIM(computername)
  write( unit, '(A1)') '#'

end subroutine PrintComputerName

!
!>Permet de fermer proprement le programme
!
SUBROUTINE StopProgram ()
    implicit none

    print *, 'Decision is taken to leave ... '

    stop
END SUBROUTINE StopProgram



