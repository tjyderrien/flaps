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
!> @brief A timer module \n
!> See IDRIS website http://www.idris.fr/su/Shared/fct_F95.html
!------------------------------------------------------------------------------

module Timer_m
  implicit none

  private

  public ::                &
    timer,                 &
    timer_init,            &
    timer_start,           &
    timer_elapsedtime

  type Timer
    integer :: nb_period_init
    integer :: nb_period_max
    integer :: nb_period_sec
  end type

  contains


  !> Init the timer
  subroutine timer_init( this )
    type(Timer), intent(inout) :: this

    ! Initialisation
    CALL SYSTEM_CLOCK(COUNT_RATE=this%nb_period_sec, COUNT_MAX=this%nb_period_max)
  end subroutine timer_init

  !> Start the timer
  subroutine timer_start( this )
    type(Timer), intent(inout) :: this

    CALL SYSTEM_CLOCK(COUNT=this%nb_period_init)
  end subroutine timer_start


  !>  Return the time elapsed from the last call of timer_start
  real function timer_elapsedtime ( this ) result(elapsed)
    type(Timer), intent(in) :: this

    integer :: nb_period_final,   & ! final value of the clock period counter
               nb_period            ! Number of clock period

    CALL SYSTEM_CLOCK(COUNT=nb_period_final)
    nb_period = nb_period_final - this%nb_period_init
    IF (nb_period_final < this%nb_period_init) &
        nb_period = nb_period + this%nb_period_max

    elapsed   = REAL(nb_period) / this%nb_period_sec
  end function timer_elapsedtime

end module Timer_m
