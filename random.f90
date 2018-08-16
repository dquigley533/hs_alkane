! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                          R A N D O M                                        !
!=============================================================================!
!                                                                             !
! $Id: random.f90,v 1.2 2011/08/02 12:56:47 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Random number generation. Can be modified to act as a wrapper to your RNG   !
! choice. Here I'm just using the intrinsic Fortran90 RANDOM_NUMBER function. !!-----------------------------------------------------------------------------!
module Random

  use iso_c_binding
  use Constants , only : dp,it
  implicit none                                 ! Impose strong typing

  private                                       ! Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: random_set_random_seed               !... unless exposed here.
  public :: random_uniform_random
  public :: random_unit_vector

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  logical,save         :: random_initialised=.false.

contains

  subroutine Random_set_random_seed(seed) bind(c)
    !=========================================================================!
    ! Set internal random number generator seed as per GNU fortran manual.    !
    ! Use seed = 0 to initialise the RNG in its default state, otherwise      !
    ! seeded from the clock.                                                  !
    !-------------------------------------------------------------------------!
    ! D. Quigley, August 2018                                                 !
    !=========================================================================!
    implicit none
    integer(kind=it), intent(in) :: seed
    integer(kind=it) :: i, n, clock, ierr
    integer,dimension(:),allocatable :: seed_array

    ! If the RNG hasn't already been seeded...
    if (.not.random_initialised) then

      ! https://gcc.gnu.org/onlinedocs/gcc-4.5.1/gfortran/RANDOM_005fSEED.html
      if (seed == 0 ) then
        call random_seed()
      else

        call random_seed(size = N)
        allocate(seed_array(1:n),stat=ierr)
        if (ierr/=0) stop 'Error allocating random seed array in random.f90'

        call system_clock(count=clock)
        seed_array = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put=seed_array)

        deallocate(seed_array)

      end if

      random_initialised = .true.                         ! it has now

    end if

    return

  end subroutine random_set_random_seed

  function random_uniform_random() bind(c)
    !=========================================================================!
    ! Return a single random deviate ~ uniform [0,1].                         !
    !-------------------------------------------------------------------------!
    ! D. Quigley, August 2018                                                 !
    !=========================================================================!
    implicit none
    real(kind=dp)                 :: Random_uniform_random
    real(kind=dp) :: x

    call random_number(x)
    random_uniform_random = x

    return

  end function random_uniform_random

  function random_unit_vector()
    !-------------------------------------------------------------------------!
    ! Generates a random orientation vector on the unit sphere. Pinched from  !
    ! cbmc.f90 by M.P. Allen, orginal algorithm by Marsaglia, 1972.           !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp), dimension(3) :: random_unit_vector

    ! Local variables
    real(kind=dp), dimension(2) :: ran
    real(kind=dp)               :: ransq, ranh

    do
       ran(1) = random_uniform_random()
       ran(2) = random_uniform_random()
       ran(:) = 1.0_dp-2.0_dp*ran(:)
       ransq = sum ( ran(:)**2 )
       if ( ransq < 1.0_dp ) exit
    end do

    ranh = 2.0 * sqrt ( 1.0 - ransq )
    random_unit_vector= (/ ran(1)*ranh, ran(2)*ranh, (1.0_dp-2.0_dp*ransq) /)

    return

  end function random_unit_vector

end module random
