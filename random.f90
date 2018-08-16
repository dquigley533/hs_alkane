! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                          R A N D O M                                        !
!=============================================================================!
!                                                                             !
! $Id: random.f90,v 1.2 2011/08/02 12:56:47 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! This is the serial version of the random number generator written by Matt   !
! Probert and used in CASTEP. This module is basically a copy of the relevant !
! routines from the CASTEP Algor module.                                      !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: random.f90,v $
! Revision 1.2  2011/08/02 12:56:47  phseal
! Added C bindings to all procedures which should be callable externally
! when compiled as a library.
!
! Revision 1.1.1.1  2011/02/02 11:48:36  phseal
! Initial import from prototype code.
!
!
!=============================================================================!

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

  !For the uniform random number generator, we presume a 32-bit integer in order
  !to do bit-shuffling, so define and store here.
  integer, parameter        :: bit_32=selected_int_kind(9)
  integer(kind=bit_32),save :: ix,iy                     !intermediate randoms

  logical,save         :: random_initialised=.false.

contains

  subroutine Random_set_random_seed(seed) bind(c)
    !=========================================================================!
    ! Set module internal random number generator seed.                       !
    ! If seed=0, then use system clock to get a unique seed every call.       !
    ! If seed<>0 then use as deterministic seed.                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   seed intent(in)=> =0 random seed, else deterministic seed.            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   iseed, ix, iy and random_initialised are set here.                    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   comms for different iseed on different nodes iff seed=0               !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   iseed as starting 32-bit integer seed for all generators              !
    !   ix as seed 32-bit integer in Marsaglia generator                      !
    !   iy as seed 32-bit integer in Park-Miller generator                    !
    !-------------------------------------------------------------------------!
    ! Architecture dependencies:                                              !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   comms must be initialised first.                                      !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 01/07/2000                               !
    !=========================================================================!
    implicit none
    integer(kind=it), intent(in) :: seed                !may or may not be 32 bit...

    !32 bit integer version of seed
    integer(kind=bit_32) :: iseed              !random number seed

    !f90 intrinsic time call
    character(len=10) :: system_time           !length is crucial ...
    real (kind=dp)    :: rtime

    !We only want to set the seed if it has not already been set
    !and the flag is set to .false. by the compiler
    if (.not.random_initialised) then

       if (seed == 0 ) then
          !if (myRank.eq.0) then
             call date_and_time(time=system_time)  !character string hhmmss.xxx
             read (system_time,*) rtime            !convert to real
             rtime = rtime * 1000.0_dp             !0<rtime<235959999.0 which fits within huge(1)
          !end if                                   !and then copy
       else
          rtime = real(abs(seed),kind=dp)          !convert seed to real
       end if

       !Be careful to only gcopy default types 
       !call Comms_BcastReal(rtime,1)

       !and then convert to bit_32 size integer
       iseed = int(rtime,kind=bit_32)              !must fit within huge(1)
       !and make sure it is different on each node
       !iseed = iseed + int(MyRank,kind=bit_32) -1_bit_32
       iseed = iseed - 1_bit_32

       ix=ieor(777755555_bit_32,iseed)                   !Marsaglia generator
       iy=ior(ieor(888889999_bit_32,iseed),1_bit_32)     !Parks-Miller generator

       random_initialised=.true.                         !set flag
    end if

    return

  end subroutine random_set_random_seed

  function random_uniform_random() bind(c)
    !=========================================================================!
    ! Return a single random deviate ~ uniform [0,1].                         !
    ! Based on Park-Miller "minimal standard" generator with Schrage's method !
    !  to do 32-bit multiplication without requiring higher precision, plus   !
    !  additional Marsaglia shift to suppress any weaknesses & correlations.  !
    ! Using two independent methods greatly increases the period of the       !
    !  generator, s.t. resulting period ~2*10^18                              !
    ! NB Routine is only set to work with 32 bit integers!                    !
    !-------------------------------------------------------------------------!
    ! References:                                                             !
    !   S.K. Park and K.W. Miller, Commun. ACM, 31, p1192-1201 (1988)         !
    !   L. Schrage, ACM Trans. Math. Soft., 5, p132-138 (1979)                !
    !   G. Marsaglia, Linear Algebra and its Applications, 67, p147-156 (1985)!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Return value:                                                           !
    !   algor_uniform_random => required random deviate                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ix as next 32-bit integer in Marsaglia generator (updated)            !
    !   iy as next 32-bit integer in Park-Miller generator (updated)          !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 01/07/2000                               !
    !=========================================================================!
    implicit none
    real(kind=dp)                 :: Random_uniform_random

    !NB We use 3 logical XOR operations to do Marsaglia shift
    !=> hard-wire 3 shift values (not all triplets any good)
    !=> entire routine preset to only work with 32 bit integers.

    !local variables ...
    integer(kind=bit_32)            :: iy_tmp       !working value to force integer division
    integer(kind=bit_32), parameter :: iy_max=2147483647 !2^31-1
    real(kind=dp), parameter        :: inv_iy_max=1.0_dp/2147483647.0_dp

    !Catch uninitialised random number generator and set to random seed
    if (.not.random_initialised) call random_set_random_seed(0)

    !do Marsaglia shift sequence, period 2^32-1, to get new ix
    ix=ieor(ix,ishft(ix, 13_bit_32))
    ix=ieor(ix,ishft(ix,-17_bit_32))
    ix=ieor(ix,ishft(ix,  5_bit_32))

    !Park-Miller sequence, period iy_max-1, to get new iy
    iy_tmp=iy/127773_bit_32                         !NB integer division
    iy=16807_bit_32*(iy-iy_tmp*127773_bit_32)-iy_tmp*2836_bit_32  !next value of iy
    if (iy < 0_bit_32) iy=iy+iy_max                 !integer modulo iy_max

    !Combine ix and iy to get new random number, rescale onto range [0,1]
    !with masking to ensure non-zero
    random_uniform_random=inv_iy_max*ior(iand(iy_max,ieor(ix,iy)),1_bit_32)

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
