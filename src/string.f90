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
! MODULE: String
!
! DESCRIPTION:
!> Module qui fournit un opérateur .equals. pour les chaines de caractères
!------------------------------------------------------------------------------

module String_m

    interface operator (.equals.)
     MODULE PROCEDURE equals
    end interface

CONTAINS

    logical function equals(s1, s2)
        CHARACTER(LEN=*), INTENT(IN) :: s1, s2

        if( lle( s1, s2 ) .and. lge( s1, s2 ) ) then
             equals = .true.
        else
             equals = .false.
        endif
        return
    end function equals
end module String_m


LOGICAL FUNCTION ReadLine( unit, line )
    implicit none
    INTEGER, INTENT(IN) :: unit
    CHARACTER(len=512) :: line
    integer :: k
    read (unit, '(A)', end=999) line
    ReadLine = .true.
    line = trim(line)
    !On traite # et ! comme un commentaire
    k = scan(line, '#!')
    if( k > 1 ) then
       line = line(1:k-1)
       line(k:512) = ' '
    endif
    return
    !READ (unit, fmt=405, advance = ‘no’, size = l, eor=106) value
999 ReadLine = .false.
    return
END FUNCTION ReadLine



!> Parse un int
SUBROUTINE ParseInt( line, value )
    implicit none
    INTEGER :: i
    CHARACTER(len=512), INTENT(IN)  :: line
    INTEGER :: value
    value = 0

    i = index(line, ' ')
    read( line(i:), '(i8)' ) value
END SUBROUTINE ParseInt

!> Parse un float
SUBROUTINE ParseFloat( line, value )
    implicit none
    INTEGER :: i
    CHARACTER(len=512), INTENT(IN)  :: line
    REAL :: value
    value = 0

    i = index(line, ' ')
    read( line(i:), '(f8.4)' ) value
END SUBROUTINE ParseFloat

!> Parse double
SUBROUTINE ParseDouble( line, value )
    implicit none
    INTEGER :: i
    CHARACTER(len=512), INTENT(IN)  :: line
    DOUBLE PRECISION :: value
    value = 0

    i = index(line, ' ')
    read( line(i:), '(e12.5)' ) value
    !TODO: input like 10e0 are not well interpreted. 
END SUBROUTINE ParseDouble


!> Parse une liste de paramétre
!> Les paramètres sont : ligne, matrice d'arrivée (nxm), n, m, ligne de la matrice, nombre de valeurs a lire
SUBROUTINE ParseTabInt( line, Tab, n, nbValues, noID )
    implicit none
    INTEGER :: i
    INTEGER, INTENT(in) :: n, nbValues
    INTEGER, DIMENSION( n) :: Tab
    CHARACTER(len=512), INTENT(IN)  :: line
    LOGICAL, INTENT(IN) :: noID

    if( noID ) then
        i = 1
    else
        i = index(line, ' ')
    endif
    read( line(i:), * ) Tab(1:nbValues)
END SUBROUTINE ParseTabInt

!> Parse une liste de paramétre
!> Les paramètres sont : ligne, matrice d'arrivée (nxm), n, m, ligne de la matrice, nombre de valeurs a lire
SUBROUTINE ParseTabFloat( line, Tab, n, nbValues, noID )
    implicit none
    INTEGER :: i
    INTEGER, INTENT(in) :: n, nbValues
    REAL, DIMENSION( n) :: Tab
    CHARACTER(len=512), INTENT(IN)  :: line
    LOGICAL, INTENT(IN) :: noID

    if( noID ) then
        i = 1
    else
        i = index(line, ' ')
    endif
    read( line(i:), * ) Tab(1:nbValues)
END SUBROUTINE ParseTabFloat

!> Parse une liste de paramétre
SUBROUTINE ParseTabDouble( line, Tab, n, nbValues, noID )
    implicit none
    INTEGER :: i
    INTEGER, INTENT(in) :: n, nbValues
    DOUBLE PRECISION, DIMENSION( n) :: Tab
    CHARACTER(len=512), INTENT(IN)  :: line
    LOGICAL, INTENT(IN) :: noID

    if( noID ) then
        i = 1
    else
        i = index(line, ' ')
    endif
    read( line(i:), * ) Tab(1:nbValues)
END SUBROUTINE ParseTabDouble

!> Parse une liste de paramétre
!> Les paramètres sont : ligne, matrice d'arrivée (nxm), n, m, ligne de la matrice, nombre de valeurs a lire
SUBROUTINE ParseMatInt( line, Mat, n, m, u, nbValues, noID )
    implicit none
    INTEGER :: i
    INTEGER, INTENT(in) :: n, m, u, nbValues
    INTEGER, DIMENSION( n,m) :: Mat
    CHARACTER(len=512), INTENT(IN)  :: line
    LOGICAL, INTENT(IN) :: noID

    if( noID ) then
        i = 1
    else
        i = index(line, ' ')
    endif
    read( line(i:), * ) Mat(u, 1:nbValues)
END SUBROUTINE ParseMatInt

!> Parse une liste de paramétre
!> Les paramètres sont : ligne, matrice d'arrivée (nxm), n, m, ligne de la matrice, nombre de valeurs a lire
SUBROUTINE ParseMatFloat( line, Mat, n, m, u, nbValues, noID )
    implicit none
    INTEGER :: i
    INTEGER, INTENT(in) :: n, m, u, nbValues
    REAL, DIMENSION( n,m) :: Mat
    CHARACTER(len=512), INTENT(IN)  :: line
    LOGICAL, INTENT(IN) :: noID

    if( noID ) then
        i = 1
    else
        i = index(line, ' ')
    endif
    read( line(i:), * ) Mat(u, 1:nbValues)
END SUBROUTINE ParseMatFloat

!> Parse une liste de paramétre
!> Les paramètres sont : ligne, matrice d'arrivée (nxm), n, m, ligne de la matrice, nombre de valeurs a lire
SUBROUTINE ParseMatDouble( line, Mat, n, m, u, nbValues, noID )
    implicit none
    INTEGER :: i
    INTEGER, INTENT(in) :: n, m, u, nbValues
    DOUBLE PRECISION, DIMENSION( n,m) :: Mat
    CHARACTER(len=512), INTENT(IN)  :: line
    LOGICAL, INTENT(IN) :: noID

    if( noID ) then
        i = 1
    else
        i = index(line, ' ')
    endif
    read( line(i:), * ) Mat(u, 1:nbValues)

END SUBROUTINE ParseMatDouble

!> Parse une liste de paramétre
!> Les paramètres sont : ligne, matrice d'arrivée (nxm), n, m, ligne de la matrice, nombre de valeurs a lire
SUBROUTINE ParseMatDoubleReverse( line, Mat, n, m, u, nbValues, noID )
    implicit none
    INTEGER :: i
    INTEGER, INTENT(in) :: n, m, u, nbValues
    DOUBLE PRECISION, DIMENSION( n,m) :: Mat
    CHARACTER(len=512), INTENT(IN)  :: line
    LOGICAL, INTENT(IN) :: noID

    if( noID ) then
        i = 1
    else
        i = index(line, ' ')
    endif
    read( line(i:), * ) Mat(1:nbValues, u)

END SUBROUTINE ParseMatDoubleReverse
