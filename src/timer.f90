!=============================================================================!
!                            T I M E R                                        !
!=============================================================================!
!                                                                             !
! $Id: timer.f90,v 1.4 2011/08/02 12:56:47 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Contains routines to check time relative to a queue slot length. For        !
! example one may wish to run the code on a system with a maximum 120 hour    !
! run length. This module allows the code to run until some specified time    !
! before 120 hours is up, and then shutdown cleanly.
!-----------------------------------------------------------------------------!
module timer 

  use iso_c_binding
  use constants, only : dp,it
  implicit none

  private                ! Everything private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.
  public :: timer_init
  public :: timer_check_runtime
  public :: timer_elapsed_time

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  public :: timer_qtime
  public :: timer_closetime
  real(kind=dp),bind(c),save :: timer_qtime = 432000         ! Queue length
  real(kind=dp),bind(c),save :: timer_closetime = 3600       ! Time taken to shutdown

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  real(kind=dp),save     :: last_time,current_time,elapsed_time
  integer(kind=it),save  :: start_day
  logical,save           :: timer_initialised = .false.

contains

  subroutine timer_init() bind(c)
    !-------------------------------------------------------------------------!
    ! Initialises the timer                                                   !
    !-------------------------------------------------------------------------!
    ! D.Quigley March 2010                                                    !
    !-------------------------------------------------------------------------!
    implicit none

    character(12) :: dat,tim,zon
    integer(kind=it),dimension(8) :: info

    call date_and_time(dat,tim,zon,info)

    last_time = 3600_dp*real(info(5),kind=dp)+60.0_dp*real(info(6),kind=dp) &
               + real(info(7),kind=dp)+0.001_dp*real(info(8),kind=dp)

    start_day  = info(3)

    write(*,*)
    write(*,'("=====================")')
    write(*,'("| Timer initialised |")')
    write(*,'("=====================")')
    write(*,*)

    elapsed_time = 0.0_dp
    timer_initialised = .true.

    return

  end subroutine timer_init

  real(kind=dp) function timer_elapsed_time() bind(c)
    !-------------------------------------------------------------------------!
    ! Returns time since timer was initialised                                !
    !-------------------------------------------------------------------------!
    ! D.Quigley March 2010                                                    !
    !-------------------------------------------------------------------------!
    implicit none

    character(12) :: dat,tim,zon
    integer,dimension(8) :: info

    if (.not.timer_initialised) then
       stop 'Called timer_elasped time before initialising timer!'
    end if

    call date_and_time(dat,tim,zon,info)

    current_time = 3600_dp*real(info(5),kind=dp)+60.0_dp*real(info(6),kind=dp) &
                   + real(info(7),kind=dp)+0.001_dp*real(info(8),kind=dp)

    ! if the day has changed....
    if (start_day/=info(3)) then
        last_time = last_time - 86400_dp
       start_day  = info(3)
    end if

    ! Compute the elapsed time
    elapsed_time = elapsed_time + current_time - last_time
    last_time    = current_time
    timer_elapsed_time = elapsed_time

    return

  end function timer_elapsed_time

  subroutine timer_check_runtime(safe) bind(c)
    !-------------------------------------------------------------------------!
    ! Checks if we are within time_closetime of timer_qtime being reached.    !
    ! Returns safe=.false. if the code should now exit cleanly.               !
    !-------------------------------------------------------------------------!
    ! D.Quigley March 2010                                                    !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(out) :: safe
    real(kind=dp)       :: tmptime
    logical,save        :: lwarn = .false.

    safe = 0

    tmptime = timer_elapsed_time()

    if (timer_qtime - tmptime>timer_closetime) then
       safe = 1
    elseif (.not.lwarn) then
       write(*,'("! ====================================================================== ! ")')
       write(*,'("! WARNING less than timer_closetime remaining before timer_qtime reached ! ")')
       write(*,'("! ====================================================================== ! ")')
       lwarn = .true.
    end if

    return

  end subroutine timer_check_runtime

  integer(kind=it) function timer_check_continuation() bind (c)
    !-------------------------------------------------------------------------!
    ! Alternative to above for wrapping with Python                           !
    !-------------------------------------------------------------------------!
    ! D.Quigley Feb 2021                                                      !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it) :: safe

    call timer_check_runtime(safe)
    timer_check_continuation = safe

    return
    
  end function timer_check_continuation
    

end module timer
