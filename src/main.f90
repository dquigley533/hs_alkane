!=============================================================================!
!                       H  S    A  L  K  A  N  E                              !
!=============================================================================!
!                                                                             !
! $Id: main.f90,v 1.8 2012/08/09 14:50:43 phrkao Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Simulation code to perform NPT simulations of hard-sphere chain alkanes.    !
!-----------------------------------------------------------------------------!
subroutine cleanexit()
  !-------------------------------------------------------------------------!
  ! Routine to handle clean shutdown of code in the event of SIGTERM being  !
  ! sent to the process by PBS or Condor etc.                               !
  !-------------------------------------------------------------------------!
  ! D.Quigley August 2011                                                   !
  !-------------------------------------------------------------------------!
  use box
  use alkane
  implicit none
  integer :: ibox

  write(0,'("Received signal SIGTERM - shutting down")')

  do ibox = 1,nboxes
     close(25+ibox-1) ! close density.dat
  end do

  !-------------------------------------!
  ! Release memory                      !
  !-------------------------------------!
  call alkane_destroy()
  call box_destroy_link_cells()
  call box_destroy()

end subroutine cleanexit

program hs_alkane

  !TODO - replace with minimal usage
  use iso_c_binding
  use constants, only : dp,ep,it,Pi
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
  integer(kind=it) :: t1,t2,rate
  integer(kind=it) :: safe

  ! Loop counters / error flags
  integer(kind=it) :: ibead,ichain,ierr,ibox

  ! Strings
  character(30) :: denfile
  character(3)  :: boxstring

  ! external signal handling subroutine
  external :: cleanexit

  io_standalone = .true.

  !--------------------------!
  ! Set timer module going   !
  !--------------------------!
  call timer_init()

  !------------------------------!
  ! Process input file           !
  !------------------------------!
  call io_read_input()

  !------------------------------!
  ! Initialise the box module    !
  !------------------------------!
  call box_initialise()

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
     call mc_initialise()

     ! Debug testing of alkane_get_external_overlaps
!!$     do ibox = 1,nboxes
!!$        do ichain = 1,nchains
!!$           call alkane_get_external_overlaps(ichain,ibox,noverlaps)
!!$           write(0,'("Chain ",I5," in box ",I5, " has ",I5," external overlaps.")')ichain,ibox,noverlaps
!!$           call alkane_get_internal_overlaps(ichain,ibox,noverlaps)
!!$           write(0,'("Chain ",I5," in box ",I5, " has ",I5," internal overlaps.")')ichain,ibox,noverlaps
!!$        end do
!!$     end do


  else
     call mc_initialise(grow_new_chains=.true.)
  end if

  !------------------------------------!
  ! Initialise I/O                     !
  !------------------------------------!
  call write_psf(nbeads,nchains)
  call write_dcd_header(Nchains,Nbeads)

  if (nboxes==1) then
     open(unit=25,file='density.dat',status='replace',iostat=ierr)
     if (ierr/=0) stop 'Error opening density.dat'
  else
     do ibox = 1,nboxes
        write(boxstring,'(".",I2.2)')ibox
        denfile = 'density.dat'//boxstring
        open (unit=25+ibox-1,file=trim(denfile),status='replace',iostat=ierr)
        if (ierr/=0) stop 'Error opening density output files'
     end do
  end if

  !------------------------------------!
  ! Setup a signal handler for SIGTERM !
  ! which is usually signal no. 15.    !
  !------------------------------------!
  call signal(15,cleanexit)

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
     if (mod(mc_cycle_num,traj_output_int)==0) call write_dcd_snapshot()

     ! Write density to density.dat every file_output_int
     do ibox = 1,nboxes
        if ((mod(mc_cycle_num,file_output_int)==0).and.pbc) then
           write(25+ibox-1,'(I10,F15.6)')mc_cycle_num,real(nchains,kind=dp)*(sigma**3)/box_compute_volume(ibox)
        end if
     end do

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
     if (.not.(safe==1)) then
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

  do ibox = 1,nboxes
     close(25+ibox-1) ! close density.dat
  end do

  !-------------------------------------!
  ! Write final snapshot in xmol format !
  !-------------------------------------!
  call io_write_xmol()


  !-------------------------------------!
  ! Release memory                      !
  !-------------------------------------!
  call alkane_destroy()
  call box_destroy_link_cells()
  call box_destroy()

end program hs_alkane
