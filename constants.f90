! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                            C O N S T A N T S                                !
!=============================================================================!
!                                                                             !
! $Id: constants.f90,v 1.2 2011/08/02 12:27:18 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Various compile-time constants used in other modules.                       !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: constants.f90,v $
! Revision 1.2  2011/08/02 12:27:18  phseal
! Updated all integers to use the integer type in constants.f90 where
! applicable. This allows the integer type it to be set to a C compatible
! type via the instrinsic iso_c_bindings module.
!
! Revision 1.1.1.1  2011/02/02 11:48:36  phseal
! Initial import from prototype code.
!
!
!=============================================================================!
module constants

  use iso_c_binding
  implicit none

  !-------------------------------------------------------------------------!
  ! Type parameters                                                         !
  !-------------------------------------------------------------------------!
  ! Fortran style
  !integer,parameter   :: dp = 8  ! default precision - might try using single
  !integer,parameter   :: ep = 8  ! extra precision  - for accumulators
  !integer,parameter   :: it = 4  ! default integer kind

  ! For C interopability - use types in iso_c_binding intrinsic module
  integer,parameter    :: dp = c_double
  integer,parameter    :: ep = c_double
  integer,parameter    :: it = c_int

  !-------------------------------------------------------------------------!
  ! Fundamental constants                                                   !
  !-------------------------------------------------------------------------!

  ! pi
  real(kind=dp),parameter :: Pi = 3.141592653589793238462643383279502884197_dp
  real(kind=dp),parameter :: invPi = 1.0_dp/3.141592653589793238462643383279502884197_dp

  ! complex
  complex(kind=dp),parameter :: cmplx_0 = (0.0_dp,0.0_dp)
  complex(kind=dp),parameter :: cmplx_1 = (1.0_dp,0.0_dp)
  complex(kind=dp),parameter :: cmplx_i = (0.0_dp,1.0_dp)

  !-------------------------------------------------------------------------!
  ! Unit conversions                                                        !
  !-------------------------------------------------------------------------!

  ! 1/(4pi epsilon_0) in atomic units
  real(kind=dp),parameter :: inv4peps = 1.0_dp

  ! kB in atomic units of Hartree/Kelvin
  real(kind=dp),parameter :: kB = 1.0_dp/3.1577465e5_dp

  ! length conversions
  real(kind=dp),parameter :: bohr_to_ang = 0.5291772108_dp
  real(kind=dp),parameter :: ang_to_bohr = 1.0_dp/0.5291772108_dp 

  ! energy conversions
  real(kind=dp),parameter :: hart_to_dlpol = 2.625501E+05_dp
  real(kind=dp),parameter :: hart_to_eV    = 27.211396181_dp
  real(kind=dp),parameter :: hart_to_SI    = 4.3597482E-18_dp

  ! mass converstions
  real(kind=dp),parameter :: aum_to_SI = 1.66053886E-27_dp

  ! pressure conversions
  real(kind=dp),parameter :: aup_to_SI  = 2.942103918E13_dp
  real(kind=dp),parameter :: aup_to_atm = 2.90363081E8_dp
  real(kind=dp),parameter :: aup_to_MPa = 2.942103918E7_dp
  real(kind=dp),parameter :: aup_to_GPa = 2.942103918E4_dp
  
  ! density conversions
  real(kind=dp),parameter :: aud_to_kgm3 = 1.120587168E4_dp
  real(kind=dp),parameter :: kgm3_to_aud = 1.0_dp/1.120587168E4_dp

end module Constants
