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
!> @file profiler.f90
!
! DESCRIPTION:
!> @brief A timer module \n
!> See IDRIS website http://www.idris.fr/su/Shared/fct_F95.html
!------------------------------------------------------------------------------

module Profiler_m
  use Timer_m
  implicit none

  private

  public ::                &
    Profiler,              &
    Profiler_init,         &
    Profiler_start,        &
    Profiler_stop,         &
    Profiler_global_init,  &
    Profiler_end_global


  type ptr
    type(Profiler), pointer :: p
  end type ptr

  type Profiler
    real           :: cumulative_time
    real           :: self_time
    integer        :: num_calls
    type(Timer)    :: timer
    logical        :: initialized = .false.
    character(256) :: name

    type(Profiler), pointer :: parent
    integer        :: nchild
    type(ptr)      :: children(512)
  end type


  type(Profiler), public  :: prof_init
  type(Profiler), public  :: prof_timeloop

  type(Profiler), target  :: prof_full
  type(Profiler), pointer :: prof_current

  contains

  !------------------------------------------------------------------
  !> Init global timer and the associated tree
  !------------------------------------------------------------------
  subroutine Profiler_global_init( )

    call Profiler_start(prof_full, "FULL")

    prof_full%nchild = 0
  end subroutine Profiler_global_init

  !------------------------------------------------------------------
  !> Init the timer
  !------------------------------------------------------------------
  subroutine Profiler_init( this, name )
    type(Profiler), target, intent(inout) :: this
    character(len=*),   intent(in) :: name

    this%cumulative_time = 0.0
    this%num_calls = 0
    call Timer_init( this%timer )
    this%initialized = .true.
    this%name = name

    !Tree related initialization
    this%parent => prof_current
    this%nchild = 0
    this%parent%nchild = this%parent%nchild +1
    this%parent%children(this%parent%nchild)%p => this
  end subroutine Profiler_init


  !------------------------------------------------------------------
  subroutine Profiler_start( this, name )
    type(Profiler) , target, intent(inout) :: this
    character(len=*),   intent(in) :: name

    if(.not.this%initialized) then
      call Profiler_init(this, name)
    end if

    if(name /= this%name) then
      print *, 'The profiler ', name, ' has already been used with a different name.'
      call StopProgram()
    endif

    call Timer_start( this%timer )
    prof_current => this

  end subroutine Profiler_start

  !------------------------------------------------------------------
  subroutine Profiler_stop( this )
    type(Profiler), intent(inout) :: this

    this%num_calls = this%num_calls + 1
    this%cumulative_time = this%cumulative_time + Timer_elapsedtime(this%timer)

    prof_current => this%parent
  end subroutine Profiler_stop

  !------------------------------------------------------------------
  subroutine Profiler_end_global( )

    integer :: iunit, ios

    call Profiler_stop(prof_full)

    iunit = 1
    OPEN( UNIT=iunit, FILE='profiling.log', STATUS="unknown",access='sequential', ACTION="write", IOSTAT=ios )
    if ( ios /= 0 ) then ! Probleme ea l'ouverture
      print *, 'Error opening profiling.log'
      call StopProgram()
    endif

    write(iunit, '(a,3x,a,3x,a,3x,a)') '# Name         ', '# of calls', 'Self time [s]', 'Cumulative time [s]', 'Time/call [s]'

   call addtorepport(iunit,prof_full)

   CLOSE( UNIT=iunit )
  end subroutine Profiler_end_global

   !------------------------------------------------------------------
   subroutine addtorepport(iunit, prof)
     integer,        intent(in) :: iunit
     type(Profiler), intent(inout) :: prof

     integer :: ichild

     !We compute the time which is spend in the profiled region, and that it is not spent in children
     prof%self_time = prof%cumulative_time
     do ichild = 1, prof%nchild
       prof%self_time = prof%self_time - prof%children(ichild)%p%cumulative_time
     end do

     write(iunit, '(a15,3x,i11,3x,f12.5,3x,f12.5,3x,f12.5)') trim(prof%name), prof%num_calls, prof%self_time, &
                       prof%cumulative_time, prof%cumulative_time/prof%num_calls

     !Recurive call
     do ichild = 1, prof%nchild
       call addtorepport(iunit, prof%children(ichild)%p)
     end do

   end subroutine addtorepport
end module Profiler_m
