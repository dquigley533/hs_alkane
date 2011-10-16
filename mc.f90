! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                                   M   C                                     !
!=============================================================================!
!                                                                             !
! $Id: mc.f90,v 1.11 2011/10/16 20:41:02 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Contains routines to perform a number of Monte-Carlo moves on hard-sphere   !
! chain models of alkanes. This module contains only the high level logic     !
! regarding selection and accept/reject mechanics. Actual moves and overlaps  !
! are performed within alkane.f90.                                            !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: mc.f90,v $
! Revision 1.11  2011/10/16 20:41:02  phseal
! Fixed lack of volume moves when using rigid chains
!
! Revision 1.10  2011/10/16 18:18:23  phseal
! Changed the minimum length to the side of a link cell to be an input
! parameter. Hence the second argument to box_construct_link_cells is
! no longer present, and link_cell_length is read from the input file.
!
! Revision 1.9  2011/09/01 16:55:13  phrkao
! Changed alkane_check_chain_geometry to be C compatible, had to change
! the argument "violate" from logical to integer and subsequently changed
! mc.f90 where this was used.
!
! Revision 1.8  2011/08/30 16:26:20  phseal
! Added ifail argument to alkane_grow_chain to indicate if at any step non
! of the ktrial segments had a viable weight and hence CBMC move should
! have zero probability of acceptance.
!
! Revision 1.7  2011/08/18 17:27:00  phrkao
! randomly had to add in i's back onto ichain and ibox
!
! Revision 1.6  2011/08/18 17:20:26  phrkao
! alkane.f90 has been updated to return quaternion and bond,angle information
! for use with lattice_switching code, bond_rotate and rotate_chain were changed.
! Dummy variables added to mc.f90 to account for this
!
! Revision 1.5  2011/08/02 12:56:47  phseal
! Added C bindings to all procedures which should be callable externally
! when compiled as a library.
!
! Revision 1.4  2011/08/02 12:31:54  phseal
! Removed use of optional arguments in routines which need to be called
! by C when compiled as a library.
!
! Revision 1.3  2011/08/02 12:27:18  phseal
! Updated all integers to use the integer type in constants.f90 where
! applicable. This allows the integer type it to be set to a C compatible
! type via the instrinsic iso_c_bindings module.
!
! Revision 1.2  2011/07/29 15:58:29  phseal
! Added multiple simulation box support.
!
! Revision 1.1.1.1  2011/02/02 11:48:36  phseal
! Initial import from prototype code.
!
!
!=============================================================================!
module mc

  use iso_c_binding
  use constants, only : dp,it
  implicit none

  private                ! Everything private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  !... unless exposed here.
  public :: mc_cycle                            ! Performs an MC sweep
  public :: mc_initialise                       ! Initialise this module

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  ! Set by the IO module upon reading of input file.
  public :: mc_target_ratio
  public :: eq_adjust_mc
  public :: max_mc_cycles
  public :: mc_cycle_num

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  real(kind=dp) :: mc_target_ratio = 0.5_dp     ! Target acceptance ratio
  logical       :: eq_adjust_mc    = .true.     ! Do we adjust MC to reach it
  integer(kind=it)      :: max_mc_cycles   = 500000000  ! How many cycles to perform
         
  integer(kind=it)      :: mc_cycle_num = 0             ! Current cycle number
  integer(kind=it)      :: last_move    = 0             ! Most recent move type

  integer(kind=it),parameter :: nmove_types = 5 ! Number of move types
  character(11),dimension(nmove_types) :: name  ! Move names

  ! Move attempt/accept counters
  integer(kind=it) :: mc_accepted_trans = 0, mc_attempted_trans = 0
  integer(kind=it) :: mc_accepted_rot   = 0, mc_attempted_rot   = 0
  integer(kind=it) :: mc_accepted_dih   = 0, mc_attempted_dih   = 0
  integer(kind=it) :: mc_accepted_box   = 0, mc_attempted_box   = 0
  integer(kind=it) :: mc_accepted_cbmc  = 0, mc_attempted_cbmc  = 0

  ! Backup chain used to restore configuration after failed MC moves
  real(kind=dp),allocatable,dimension(:,:) :: backup_chain

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!
 

contains

  subroutine mc_cycle
    !-------------------------------------------------------------------------!
    ! Performs nbeads x nchains trial moves selected at random, also updates  !
    ! MC move parameters to match the specified acceptance ratio (not to be   !
    ! used for production runs). Move types are volume/box moves, CBMC chain  !
    ! regrowth moves, chain translation moves, chain rotation moves, and      !
    ! internal chain torsion angle moves.                                     !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : dp,ep,Pi
    use random,    only : random_uniform_random
    use box      , only : pbc,pressure,box_construct_link_cells,nboxes
    use alkane   , only : Rchain,ktrial,mc_dh_max,mc_dr_max,mc_dt_max,mc_dv_max,&
                          nbeads,nchains,rigid,alkane_box_resize,alkane_grow_chain,&
                          alkane_bond_rotate,alkane_translate_chain,&
                          alkane_rotate_chain,alkane_check_chain_overlap,&
                          alkane_check_chain_geometry,alkane_update_linked_lists
    
    implicit none

    real(kind=dp) :: xi                  	   ! Random number on 0-1.
    real(kind=dp) :: acc_prob            	   ! Acceptance probability
    real(kind=dp) :: old_rb_factor       	   ! Old Rosenbluth factor
    real(kind=dp) :: new_rb_factor       	   ! New Rosenbluth factor
    real(kind=dp) :: new_boltz          	   ! New Boltmann factor
    real(kind=dp) :: acrat              	   ! Acceptance ratio
    real(kind=dp) :: dummy_dp			   ! Dummy variable
    real(kind=dp) :: volP                          ! Vol move prob

    real(kind=dp),dimension(4) :: dummy_quat  	   ! Dummy quarternion 

    logical :: overlap				   ! Check on overlaps
    integer(kind=it) :: violated           	   ! Check on geometry

    ! Loop counters
    integer(kind=it) :: ipass,ichain,ibead,k,ibox,dummy_int,ifail

    !-------------------------------------!
    ! nchains trial moves per cycle       !
    !-------------------------------------!
    do ipass = 1,nchains*nbeads*nboxes

       ! Generatate a random number to choose which box to operate on
       xi = random_uniform_random()
       ibox = int(xi*real(nboxes,kind=dp)) + 1
       ibox = min(nboxes,ibox)

       ! Generate a random number to choose a chain to operate on
       xi = random_uniform_random()
       ichain = int(xi*real(nchains,kind=dp))+1
       ichain = min(nchains,ichain)

       volP = 1.0_dp/real(nchains*nbeads,kind=dp)

       ! Generate a third random number to select a move type
       xi = random_uniform_random()
       if (nchains<2) xi = 0.5_dp*xi ! don't attempt upper half of moveset for single chains
       if (rigid)     then
          xi = 0.5_dp+0.5_dp*xi ! ..or lower half for rigid chains
          volP = volP + 0.5_dp
       end if

       if ( rigid.and.(nchains<2) ) stop 'No meaningfull moves in mc.f90'

       !----------------------------------------------------------!
       ! Box resize moves, these may be isotropic or anisotropic  !
       ! Attempt with prob 1/(nchains*nbeads). Don't attempt if   !
       ! if there are no periodic boundary conditions in use.     !
       ! ---------------------------------------------------------!
       if ( (xi<volP ).and.pbc ) then

          ! Resize the box, and return the new acceptance ratio
          ! Setting reset flag to zero
          call alkane_box_resize(pressure,ibox,acc_prob,0)
          mc_attempted_box = mc_attempted_box + 1

          ! Generate a random number
          xi = random_uniform_random()

          if ( xi < acc_prob ) then
             ! Accept the new config
             mc_accepted_box = mc_accepted_box + 1
             last_move = 1
          else
             ! Reset, call the routine with the reset flag.
             call alkane_box_resize(pressure,ibox,acc_prob,1)
          end if

          !----------------------------------------------------------!
          ! Configurational bias chain regrowth moves.               !
          ! Attempted with prob 0.4 - 1/(nchains*nbeads) unless we   !
          ! have a rigid/inflexible chains.                          !
          ! ---------------------------------------------------------!
       else if (xi < 0.4_dp) then

          mc_attempted_cbmc = mc_attempted_cbmc + 1

          ! 50 % probability of reversing the labelling of the
          ! chain, i.e. we can regrow from either end.
          xi = random_uniform_random()
          if ( xi > 0.5_dp ) then

             ! Reverse the chain and regrow from the other end
             k = 1
             backup_chain(:,:) = Rchain(:,:,ichain,ibox)   
             do ibead = nbeads,1,-1
                Rchain(:,k,ichain,ibox) = backup_chain(:,ibead)
                k = k + 1
             end do

             ! Update linked lists to reflect the renumbering
             do ibead = 1,nbeads
                call alkane_update_linked_lists(ibead,ichain,ibox, &
                     backup_chain(:,ibead),Rchain(:,ibead,ichain,ibox))
             end do

          end if

          ! Store the backup chain
          backup_chain(:,:) = Rchain(:,:,ichain,ibox)   

          ! Compute the Rosenbluth factor of the old chain
          ! i.e. call grow_chain with regrow = .false.
          call alkane_grow_chain(ichain,ibox,old_rb_factor,0,ifail)

          ! Compute the Rosenbluth factor of the new chain
          ! i.e. call grow_chain with regrow = .true.
          call alkane_grow_chain(ichain,ibox,new_rb_factor,1,ifail)

          ! Generate a random number and accept of reject the move
          xi = random_uniform_random()
          if ( xi < new_rb_factor/old_rb_factor ) then

             ! Accept the new config
             mc_accepted_cbmc = mc_accepted_cbmc + 1
             last_move = 2

             ! Update the linked lists
             do ibead = 1,nbeads
                call alkane_update_linked_lists(ibead,ichain,ibox, &
                     backup_chain(:,ibead),Rchain(:,ibead,ichain,ibox))
             end do

          else

             ! Reject and put the old chain back
             Rchain(:,:,ichain,ibox) = backup_chain(:,:)

          end if


          !----------------------------------------------------------!
          ! Internal torsion angle move, only attempted if chains    !
          ! are not constrained to be rigid.                         !
          ! ---------------------------------------------------------!
       elseif ( xi<0.5_dp ) then

          mc_attempted_dih = mc_attempted_dih + 1

          ! 50 % probability of reversing the labelling of the
          ! chain, i.e. we can rotate in either direction
          xi = random_uniform_random()
          if ( xi > 0.5_dp ) then

             ! reverse the chain and regrow from the other end
             k = 1

             backup_chain(:,:) = Rchain(:,:,ichain,ibox)   
             do ibead = nbeads,1,-1
                Rchain(:,k,ichain,ibox) = backup_chain(:,ibead)
                k = k + 1
             end do

             ! Update linked lists to reflect the renumbering
             do ibead = 1,nbeads
                call alkane_update_linked_lists(ibead,ichain,ibox, &
                     backup_chain(:,ibead),Rchain(:,ibead,ichain,ibox))
             end do

          end if

          ! Store the backup chain
          backup_chain(:,:) = Rchain(:,:,ichain,ibox)  

          ! Rotate a randomly selected torsion angle along the chain
          ! TODO - this will fail for model III with continuous
          ! torsion potentials.
          call alkane_bond_rotate(ichain,ibox,new_boltz,dummy_int,dummy_dp)

          ! Accept or reject move.
          xi = random_uniform_random()
          if ( xi < new_boltz ) then  

             ! Accept the new config
             mc_accepted_dih = mc_accepted_dih + 1
             last_move = 3

             ! Update the linked lists for the new config
             do ibead = 1,nbeads
                call alkane_update_linked_lists(ibead,ichain,ibox, &
                     backup_chain(:,ibead),Rchain(:,ibead,ichain,ibox))
             end do

          else

             ! reject and put the old position back
             Rchain(:,:,ichain,ibox) = backup_chain(:,:)

          end if

          !----------------------------------------------------------!
          ! Translate a single chain.                                !
          ! ---------------------------------------------------------!
       elseif ( xi < 0.75_dp ) then

          mc_attempted_trans = mc_attempted_trans + 1

          ! Store backup chain
          backup_chain(:,:) = Rchain(:,:,ichain,ibox)  

          ! Translate entire chain by a random vector. 
          ! new_boltz is the new Boltzmann factor
          call alkane_translate_chain(ichain,ibox,new_boltz)

          ! Accept or reject this trial move
          xi = random_uniform_random()
          if ( xi < new_boltz ) then

             ! Accept the new config
             mc_accepted_trans = mc_accepted_trans + 1
             last_move = 4

             ! Update linked lists for the new chain position
             do ibead = 1,nbeads
                call alkane_update_linked_lists(ibead,ichain,ibox, &
                     backup_chain(:,ibead),Rchain(:,ibead,ichain,ibox))
             end do

          else

             ! Reject and put the old position back
             Rchain(:,:,ichain,ibox) = backup_chain(:,:)

          end if

          !----------------------------------------------------------!
          ! Rotate a single chain about its centre of mass           !
          ! ---------------------------------------------------------!
       else 

          mc_attempted_rot = mc_attempted_rot + 1

          ! Store backup chain
          backup_chain(:,:) = Rchain(:,:,ichain,ibox)  

          ! Rotate by a random angle about a random axis
          ! new_boltz holds the new Boltzmann factor
          call alkane_rotate_chain(ichain,ibox,new_boltz,dummy_quat)

          ! Accept or reject.
          xi = random_uniform_random()
          if ( xi < new_boltz ) then

             ! Accept the new config
             mc_accepted_rot = mc_accepted_rot + 1
             last_move = 5

             ! Update linked-lists for the new chain 
             do ibead = 1,nbeads
                call alkane_update_linked_lists(ibead,ichain,ibox, &
                     backup_chain(:,ibead),Rchain(:,ibead,ichain,ibox))
             end do

          else

             ! Reject and put the old position back
             Rchain(:,:,ichain,ibox) = backup_chain(:,:)

          end if

       end if

       ! Check the chain geometry
       call alkane_check_chain_geometry(ichain,ibox,violated)
       if (violated == 1) then
          write(0,'("Last move was of type : ",A11)')name(last_move)
          stop
       end if

    end do

    ! Check for chain overlap periodically
    if (mod(mc_cycle_num,1000)==0) then
       overlap = .false.
       if (nchains>1) call alkane_check_chain_overlap(ibox,overlap)
       if (overlap) then
          write(0,'("Last move was of type : ",A11)')name(last_move)
          stop
       end if
    end if

    !--------------------------------------------------------------------!
    ! Check move acceptance rations and adjust move parameters if flag   !
    ! eq_adjust_mc is true. Should only be used to select moves params   !
    ! for subsequent simulations with fixed values.                      !
    !--------------------------------------------------------------------!
    if ( (mod(mc_cycle_num,1000)==0).and.eq_adjust_mc ) then

       write(*,'("!=======================================!")') 
       write(*,'("! Checking move acceptance ratios...    !")')
       write(*,'("!                                       !")')

       ! Enforce sensible minima here
       acrat  = real(mc_accepted_cbmc)/real(mc_attempted_cbmc)
       ktrial = max(5,int(real(ktrial,kind=dp)*mc_target_ratio/acrat))
       ktrial = min(10,ktrial)

       write(*,'("! Configurational bias moves   : ",I2" %   ")')nint(acrat*100.0_dp)

       acrat = real(mc_accepted_dih)/real(mc_attempted_dih)
       mc_dh_max = max(0.001,mc_dh_max*acrat/mc_target_ratio)
       mc_dh_max = min(mc_dh_max,Pi)

       write(*,'("! Dihedral angle moves         : ",I2," %   ")')nint(acrat*100.0_dp)

       if (nchains > 1) then

       acrat = real(mc_accepted_trans)/real(mc_attempted_trans)
       mc_dr_max = max(0.001,mc_dr_max*acrat/mc_target_ratio)

       write(*,'("! Molecule translation moves   : ",I2," %   ")')nint(acrat*100.0_dp)

       acrat = real(mc_accepted_rot)/real(mc_attempted_rot)
       mc_dt_max = max(0.001,mc_dt_max*acrat/mc_target_ratio)
       mc_dt_max = min(mc_dt_max,Pi)

       write(*,'("! Molecule rotation moves      : ",I2," %   ")')nint(acrat*100.0_dp)

       end if

       if (pbc) then

       acrat = real(mc_accepted_box)/real(mc_attempted_box)
       mc_dv_max = max(0.0001,mc_dv_max*acrat/mc_target_ratio)

       write(*,'("! Box moves                    : ",I2," %   ")')nint(acrat*100.0_dp)

       end if

       mc_accepted_trans = 0 ; mc_attempted_trans = 0
       mc_accepted_rot   = 0 ; mc_attempted_rot   = 0
       mc_accepted_dih   = 0 ; mc_attempted_dih   = 0
       mc_accepted_box   = 0 ; mc_attempted_box   = 0
       mc_accepted_cbmc  = 0 ; mc_attempted_cbmc  = 0

       write(*,'("!=======================================!")') 
       write(*,*)

       write(*,'("!=======================================!")')
       write(*,'("! Updated move parameters               !")')
       write(*,'("!                                       !")')
       write(*,'("! ktrial (CBMC)        : ",I10  , "     !")')ktrial
       write(*,'("! mc_dh_max            : ",F10.4, "     !")')mc_dh_max
       if (nchains>1) then
       write(*,'("! mc_dr_max            : ",F10.4, "     !")')mc_dr_max
       write(*,'("! mc_dt_max            : ",F10.4, "     !")')mc_dt_max
       end if
       if (pbc) then
       write(*,'("! mc_dv_max            : ",F10.4, "     !")')mc_dv_max
       end if
       write(*,'("!=======================================!")')
       write(*,*)

    end if


  end subroutine mc_cycle


  subroutine mc_initialise(grow_new_chains)
    !-------------------------------------------------------------------------!
    ! Populates the chain arrays in alkane.f90 with new chains grown by CBMC  !
    ! checks initial chain configurations and allocated memory for a backup   !
    ! chain used to reset after failed MC moves.                              !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use box,       only : box_construct_link_cells,nboxes
    use constants, only : ep,dp
    use alkane,    only : alkane_grow_chain,nchains,nbeads,&
                          alkane_check_chain_overlap,alkane_check_chain_geometry,&
                          alkane_construct_linked_lists
    implicit none
    logical,optional,intent(in) :: grow_new_chains

    logical :: overlap                  	 ! Check on overlaps 
    integer(kind=it) :: violated		 ! Check on geometry
    integer(kind=it) :: ichain,nccopy,ibox,ifail ! Loop counters etc
    real(kind=dp) :: rb_factor                   ! Rosenbluth factor for chain growth

    ! Loop over simulation boxes
    do ibox = 1,nboxes

       ! Grow new chains by CBMC
       if (present(grow_new_chains).and.grow_new_chains) then

          !------------------------------------!
          ! Grow initial chain configurations  !
          !------------------------------------!
          write(*,'("!=======================================!")')
          write(*,'("! Growing initial chain configurations  !")')
          nccopy = nchains
          do ichain = 1,nccopy
             nchains = ichain ! only check overlap for j < i
             do
                call alkane_grow_chain(ichain,ibox,rb_factor,1,ifail)
                if (rb_factor>tiny(1.0_ep)) exit
             end do
             write(*,'("! ...",I5)')ichain
          end do
          nchains = nccopy
          write(*,'("! ....done                              !")')
          write(*,'("!=======================================!")')
          write(*,*)

       end if

       ! Construct link cells and build linked-lists
       call box_construct_link_cells(ibox)
       call alkane_construct_linked_lists(ibox)

       ! Check initial configuration is sane
       overlap = .false.
       if (nchains>1) call alkane_check_chain_overlap(ibox,overlap)
       if (overlap) stop
       do ichain = 1,nchains
          call alkane_check_chain_geometry(ichain,ibox,violated)
       end do
       if (violated == 1) stop

    end do

    ! Allocate backup chain
    allocate(backup_chain(1:3,1:nbeads))

    ! Set move names
    name(1) = 'Box resize '
    name(2) = 'CBMC       '
    name(3) = 'Dihedral   '
    name(4) = 'Translation'
    name(5) = 'Rotation   '

    return

  end subroutine mc_initialise


end module mc
