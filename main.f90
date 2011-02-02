! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                       H  S    A  L  K  A  N  E                              !
!=============================================================================!
!                                                                             !
! $Id: main.f90,v 1.1 2011/02/02 11:48:36 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Simulation code to perform NPT simulations of hard-sphere chain alkanes.    !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: main.f90,v $
! Revision 1.1  2011/02/02 11:48:36  phseal
! Initial revision
!
!
!=============================================================================!


program hs_alkane

  !TODO - replace with minimal usage
  use constants, only : dp,ep,Pi
  use timer
  use alkane
  use random
  use box
  use vis
  use io
  use mc
  implicit none

  !========================================================!
  ! Local variables                                        !
  !========================================================!
  ! Timing
  integer :: t1,t2,rate
  logical :: safe 

  ! Loop counters / error flags
  integer :: ibead,ichain,ierr

  !--------------------------!
  ! Set timer module going   !
  !--------------------------!
  call timer_init()

  !------------------------------!
  ! Process input file           !
  !------------------------------!
  call io_read_input()

  !------------------------------!
  ! Initialise the alkane module !
  !------------------------------!
  call alkane_init()

  !------------------------------------!
  ! Set initial chains and initialise  !
  ! Monte-Carlo module.                !
  !------------------------------------!
  if (read_xmol) then
     call io_read_xmol()
     call mc_initialise_chains()
  else
     call mc_initialise_chains(grow_new_chains=.true.)
  end if

  !------------------------------------!
  ! Initialise I/O                     !
  !------------------------------------!
  call write_psf(nbeads,nchains)
  call write_dcd_header(Nchains,Nbeads)

  open(unit=25,file='density.dat',status='replace',iostat=ierr)
  if (ierr/=0) stop 'Error opening density.dat'

  !------------------------------------!
  ! Time measurement for cycle rate    !
  !------------------------------------!
  call system_clock(count=t1)

  !------------------------------------!
  ! Loop over MC cycles                !
  !------------------------------------!
  do mc_cycle_num = 1,max_mc_cycles

     ! Perform a MC cycle
     call mc_cycle()

     ! Write snapshot of system to DCD file every traj_output_int
     if (mod(mc_cycle_num,traj_output_int)==0) call write_dcd_snapshot(Nchains,nbeads,Rchain)

     ! Write density to density.dat every file_output_int
     if ((mod(mc_cycle_num,file_output_int)==0).and.pbc) then
        write(25,'(I10,F15.6)')mc_cycle_num,real(nchains,kind=dp)*(sigma**3)/box_compute_volume()
     end if

     ! Every 1000 cycles write the current cycle rate to stdout
     if (mod(mc_cycle_num,1000)==0) then

        call system_clock(count=t2,count_rate=rate) ! new count

        write(*,*)
        write(*,'("!=======================================!")')
        write(*,'("! MC cycles per second : ",F10.4, "     !")') &
             real(1000,kind=ep)/(real(t2-t1,kind=ep)/real(rate,kind=ep))
        write(*,'("!=======================================!")')
        write(*,*)

        call system_clock(count=t1,count_rate=rate)  ! reset clock

     end if

     ! Check if we should exit the program, i.e. are we within 
     ! time_closetime of reaching timer_qtime.
     call timer_check_runtime(safe)
     if (.not.safe) then
        write(*,*)
        write(*,'("!============================================!")')
        write(*,'("! Approaching end of queue time - stopping   !")')
        write(*,'("!============================================!")')
        write(*,*)
        exit
     else if (mod(mc_cycle_num,1000)==0) then
        write(*,*)
        write(*,'("!=======================================!")')
        write(*,'("! Time elapsed (s) : ",F18.6)')timer_elapsed_time()
        write(*,'("!=======================================!")')
        write(*,*)
     end if

  end do

  close(25) ! close density.dat

  !-------------------------------------!
  ! Write final snapshot in xmol format !
  !-------------------------------------!
  open(unit=25,file='final.xmol',status='replace',iostat=ierr)
  if (ierr/=0) stop 'Error opening final.xmol for writing'

  write(25,*)nbeads*nchains
  write(25,'("* ",9F15.6)')hmatrix
  do ichain = 1,nchains
     do ibead = 1,nbeads
        write(25,'("C ",3F15.6)')Rchain(:,ibead,ichain)
     end do
  end do

  close(25)

  !-------------------------------------!
  ! Release memory                      !
  !-------------------------------------!
  call alkane_destroy()
  call box_destroy_link_cells()


end program hs_alkane
