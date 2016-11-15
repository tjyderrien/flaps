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
!> @file Timer.f90
!
! DESCRIPTION:
!> @brief Classe qui fournit les routines pour un timer\n
!> Voir le site de l'IDRIS http://www.idris.fr/su/Shared/fct_F95.html
!------------------------------------------------------------------------------

!> Methode qui initialise le timer
SUBROUTINE TimerInit  ( )
    INTEGER :: nb_periodes_initial,&
    nb_periodes_max,     & ! valeur maximale du compteur d'horloge
    nb_periodes_sec        ! nombre de periodes d'horloge par seconde
    COMMON /InitialPeriod/ nb_periodes_initial
    COMMON /PeriodMax/ nb_periodes_max
    COMMON /PeriodSec/ nb_periodes_sec

    ! Initialisation
    CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, COUNT_MAX=nb_periodes_max)
END SUBROUTINE TimerInit

!> Methode qui demarre le timer
SUBROUTINE TimerStart( )
    implicit none
    INTEGER :: nb_periodes_initial ! valeur initiale du compteur de periodes d'horloge
    COMMON /InitialPeriod/ nb_periodes_initial

    CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)
END SUBROUTINE TimerStart


!> Fonction qui retourne le temps en seconde ecoule depuis l'appel de la fonction Timer_Start
REAL FUNCTION ElapsedTime ( )
    implicit none
    INTEGER :: &
    nb_periodes_initial, & ! valeur initiale du compteur de periodes d'horloge
    nb_periodes_final,   & ! valeur finale   du compteur de periodes d'horloge
    nb_periodes_max,     & ! valeur maximale du compteur d'horloge
    nb_periodes_sec,     & ! nombre de periodes d'horloge par seconde
    nb_periodes            ! nombre de periodes d'horloge du code
    COMMON /InitialPeriod/ nb_periodes_initial
    COMMON /PeriodMax/ nb_periodes_max
    COMMON /PeriodSec/ nb_periodes_sec

    CALL SYSTEM_CLOCK(COUNT=nb_periodes_final)
    nb_periodes = nb_periodes_final - nb_periodes_initial
    IF (nb_periodes_final < nb_periodes_initial) &
        nb_periodes = nb_periodes + nb_periodes_max
    ElapsedTime   = REAL(nb_periodes) / nb_periodes_sec
END FUNCTION
