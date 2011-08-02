! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                                   M   C                                     !
!=============================================================================!
!                                                                             !
! $Id: mc_dummy.f90,v 1.2 2011/08/02 12:27:18 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Dummy MC module for use when compiling associated routines into a library.  !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: mc_dummy.f90,v $
! Revision 1.2  2011/08/02 12:27:18  phseal
! Updated all integers to use the integer type in constants.f90 where
! applicable. This allows the integer type it to be set to a C compatible
! type via the instrinsic iso_c_bindings module.
!
! Revision 1.1  2011/08/02 10:55:11  phseal
! Initial version for compilation as a library
!
!
!
!=============================================================================!
module mc

  use constants, only : dp
  implicit none

  private                ! Everything private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  !... unless exposed here.

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  ! Set by the IO module upon reading of input file.
  public :: mc_target_ratio
  public :: eq_adjust_mc
  public :: max_mc_cycles

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  real(kind=dp)    :: mc_target_ratio = 0.5_dp     ! Target acceptance ratio
  logical          :: eq_adjust_mc    = .true.     ! Do we adjust MC to reach it
  integer(kind=it) :: max_mc_cycles   = 500000000  ! How many cycles to perform


  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!
 
end module mc
