!------------------------------------------------------------------------------
!> @file Timer.f90
!
! DESCRIPTION:
!> @brief Classe qui fournit les routines pour un timer\n
!> Voir le site de l'IDRIS http://www.idris.fr/su/Shared/fct_F95.html
!
!> @author
!> Nicolas Tancogne-Dejean
!> @date
!> 01 Aug 2012 - Initial Version
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
