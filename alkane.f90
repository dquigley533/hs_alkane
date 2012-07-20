! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                            A  L  K  A  N  E                                 !
!=============================================================================!
!                                                                             !
! $Id: alkane.f90,v 1.29 2012/07/20 08:58:32 phrkao Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Contains routines to store and manipulate (i.e. attempt trial MC moves) a   !
! system of one or more alkane chains within the various models of            !
! Malanoski, Monson and Cao. Bond lengths are held constant, as are angles.   !
! Uses the link-cell structure in box.f90, and returns Boltzmann factors for  !
! any trial moves attempted. These will be one or zero in most cases.         !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: alkane.f90,v $
! Revision 1.29  2012/07/20 08:58:32  phrkao
! removed first_bead=1 from alkane_grow_chain
!
! Revision 1.28  2012/07/05 13:04:42  phrkao
! merged
!
! Revision 1.27  2012/06/19 16:40:22  phrkao
! changed centre of mass to first bead
!
! Revision 1.26  2011/11/21 16:07:47  phseal
! Removed unused variables and refactored volume moves to prevent compiler warnings
!
! Revision 1.25  2011/11/21 12:53:41  phseal
! Stopped alkane_construct_linked_lists from trying to build before cells constructed
!
! Revision 1.24  2011/11/21 11:24:22  phseal
! Fixed out-of-bounds errors when one box had less link cells
!
! Revision 1.23  2011/11/04 16:12:20  phseal
! Added more informative check of link cell consistency
!
! Revision 1.22  2011/10/26 16:16:03  phrkao
! *** empty log message ***
!
! Revision 1.21  2011/10/26 16:14:57  phrkao
! made alkane_check_chain_overlap c compatible and corresponding mc.f90 calls adapted
!
! Revision 1.20  2011/10/19 16:21:42  phseal
! alkane_check_overlap will now check for overlaps both with and without the
! use of linked lists, and compare the result of the two methods.
!
! Revision 1.19  2011/10/18 09:26:46  phseal
! Reduced dihedral violation errors to warnings, as these were being triggered
! for model IV when the angle was violated by rounding errors, i.e. small differences
! in the result when computing the angle for construction and checking purposes.
!
! Revision 1.18  2011/10/16 18:18:23  phseal
! Changed the minimum length to the side of a link cell to be an input
! parameter. Hence the second argument to box_construct_link_cells is
! no longer present, and link_cell_length is read from the input file.
!
! Revision 1.17  2011/09/01 16:58:29  phrkao
! *** empty log message ***
!
! Revision 1.15  2011/08/30 16:26:19  phseal
! Added ifail argument to alkane_grow_chain to indicate if at any step non
! of the ktrial segments had a viable weight and hence CBMC move should
! have zero probability of acceptance.
!
! Revision 1.14  2011/08/30 16:00:15  phseal
! Trapped mostly harmless computation of NaN for segment selection probability
!
! Revision 1.13  2011/08/30 12:39:03  phseal
! Caught case where RNG returned 1.00000 in alkane_grow_chain
!
! Revision 1.12  2011/08/30 10:49:24  phseal
! Added routines to manipulate max_regrow from C
!
! Revision 1.11  2011/08/30 10:45:53  phseal
! Fixed various issues with number of beads regrown in cbmc
!
! Revision 1.10  2011/08/18 17:27:00  phrkao
! randomly had to add in i's back onto ichain and ibox
!
! Revision 1.9  2011/08/18 17:20:26  phrkao
! alkane.f90 has been updated to return quaternion and bond,angle information
! for use with lattice_switching code, bond_rotate and rotate_chain were changed.
! Dummy variables added to mc.f90 to account for this
!
! Revision 1.8  2011/08/02 12:56:47  phseal
! Added C bindings to all procedures which should be callable externally
! when compiled as a library.
!
! Revision 1.7  2011/08/02 12:30:53  phseal
! mc.f90
!
! Revision 1.6  2011/08/02 12:27:18  phseal
! Updated all integers to use the integer type in constants.f90 where
! applicable. This allows the integer type it to be set to a C compatible
! type via the instrinsic iso_c_bindings module.
!
! Revision 1.5  2011/08/02 10:04:18  phseal
! Added routines to manipulate inactive boxes from outside of the box and
! alkane module. Updates overlap counting routines to return only the
! total number of overlaps found rather than lists of overlapping atoms.
!
! Revision 1.4  2011/07/29 19:28:30  phseal
! Added experimental routines to count and store the number of intra-chain
! and inter-chain bead-bead overlaps. For use on inactive lattices with
! lattice-switching calculations.
!
! Revision 1.3  2011/07/29 15:58:29  phseal
! Added multiple simulation box support.
!
! Revision 1.2  2011/03/11 13:47:19  phseal
! Moved from single to double linked-lists
!
! Revision 1.1.1.1  2011/02/02 11:48:36  phseal
! Initial import from prototype code.
!
!
!=============================================================================!

module alkane

  use iso_c_binding
  use constants, only : dp,ep,it
  implicit none

  private                ! Everything private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.

  public :: alkane_init                       ! Initialise this module
  public :: alkane_destroy                    ! Shutdown this module
  public :: alkane_grow_chain                 ! Grow/regrow chain (CBMC)
  public :: alkane_translate_chain            ! Trial chain translation 
  public :: alkane_rotate_chain               ! Trial chain rotation
  public :: alkane_box_resize                 ! Trial box move
  public :: alkane_bond_rotate                ! Trial internal torsion move
  public :: alkane_check_chain_overlap        ! Check for overlapping chains
  public :: alkane_construct_neighbour_list   ! Build neighbour lists
  public :: alkane_construct_linked_lists     ! Build linked-lists
  public :: alkane_update_linked_lists        ! Update linked-lists
  public :: alkane_check_chain_geometry       ! Check chain constraints
  public :: alkane_box_scale                  ! Scale box anisotropically

  public :: alkane_get_external_overlaps      ! Count overlaps involving chain i
  public :: alkane_get_internal_overlaps      ! Count overlaps within chain i

  public :: alkane_get_dr_max,alkane_set_dr_max  ! Manipulate dr_max externally
  public :: alkane_get_dt_max,alkane_set_dt_max  ! Manipulate dt_max externally
  public :: alkane_get_dv_max,alkane_set_dv_max  ! Manipulate dv_max externally
  public :: alkane_get_dh_max,alkane_set_dh_max  ! Manipulate dh_max externally
  public :: alkane_get_ktrial,alkane_set_ktrial  ! Manipulate ktrial externally
  public :: alkane_get_max_regrow,alkane_set_max_regrow

  public :: alkane_get_nchains                   ! Query number of chain per box
  public :: alkane_get_nbeads                    ! Query number of beads per chain

  public :: alkane_get_chain                     ! Query coordinates of a single chain
  public :: alkane_set_chain                     ! Update coordinates of a single chain

  public :: alkane_change_box                    ! Updates the matrix of cell vectors
                                                 ! and scales chain C.O.M.s with it.

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  public :: chain_created                   ! Coordinates present for chain?
  public :: Rchain                          ! Chain coordinates
  public :: nchains                         ! Number of chains
  public :: nbeads                          ! No. beabs per chain
  public :: rigid                           ! Are chains rigid?
  public :: sigma                           ! Bead radius
  public :: L                               ! Bond length
  public :: model_type,torsion_type         ! Model specification
  public :: ktrial                          ! No of trial regrowths per bead
  public :: max_regrow                      ! No. beads to regrow in CBMC
  public :: mc_dr_max                       ! Maximum whole chain displacement
  public :: mc_dt_max                       ! Maximum whole chain rotation
  public :: mc_dv_max                       ! Maximum volume move
  public :: mc_dh_max                       ! Maximum torsion angle change


  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  ! Model parameters
  integer(kind=it),save  :: nchains = 35        ! Number of chains per box
  integer(kind=it),save  :: nbeads  = 4         ! length of alkane chains
  real(kind=dp),save     :: sigma   = 1.0_dp    ! hard sphere bead diameter
  real(kind=dp),save     :: L       = 0.4_dp    ! bond length

  !---------------------------------------------------------------!
  !                      Model specification                      !
  !---------------------------------------------------------------!
  !                                                               !
  ! model_type   = 1 : HS interactions i-j for j > i+3            !
  !              = 2 : HS interactions i-j for j > i+2            !
  !              = 3 : As 1 with explicit torsion potential       !
  !              = 4 : As 1 with 'hard' torsion potential         !
  !                                                               !
  ! torsion_type = 1 : Ryckaert-Bellemans torsion potenial        !
  !              = 2 : Padilla-Toxvaerd torsion potential         !
  !                                                               !
  !---------------------------------------------------------------!
  ! WARNING - this code adopts the convention that the zero of    !
  ! the dihedral angle is taken in the trans configurations,      !
  ! hence gauche conformations are found at +/- 120 degrees. This !
  ! matches the convention used in defining the three torsion     !
  ! potentials but differs from most other sources by 180 degrees !
  !---------------------------------------------------------------!
  integer(kind=it),save   :: model_type    = 4    ! default model type
  integer(kind=it),save   :: torsion_type  = 1    ! default torsion type
  integer(kind=it),save   :: nexclude             ! exclusion length
  
  ! Chain coordinates and flag for initial configuration
  real(kind=dp),allocatable,dimension(:,:,:,:),save :: Rchain
  logical,allocatable,dimension(:,:),save           :: chain_created

  ! Intermolecular neighbour list (if applicable)
  integer(kind=it),allocatable,dimension(:,:),save :: list,startinlist,endinlist

  ! Linked list arrays (if applicable) See F&S for details
  integer(kind=it),allocatable,dimension(:,:,:) ,save  :: head_of_cell
  integer(kind=it),allocatable,dimension(:,:,:,:),save :: linked_list
  logical,save       :: rigid     = .false.      ! Chains are inflexible

  ! CBMC parameters
  integer(kind=it),save :: ktrial = 5                     ! Number of trial segments
  integer(kind=it),save :: max_regrow = 3                 ! Maximum segments to regrow

  ! Other MC move parameters
  real(kind=dp),save :: mc_dr_max = 0.0421_dp    ! Maximum translation move
  real(kind=dp),save :: mc_dt_max = 0.3181_dp    ! Maximum rotation angle
  real(kind=dp),save :: mc_dv_max = 0.2958_dp    ! Maximum volume change
  real(kind=dp),save :: mc_dh_max = 0.0159_dp    ! Maximum dihedral change

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!


contains

  subroutine alkane_init() bind(c)
    !-------------------------------------------------------------------------!
    ! Initialises the alkane module                                           !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use box, only : nboxes
    implicit none

    integer(kind=it),dimension(3) :: ierr = 0

    ! Allocate chain position array
    allocate(Rchain(1:3,1:nbeads,1:nchains,1:nboxes),stat=ierr(1))
    if (any(ierr/=0)) stop 'Error allocating memory in alkane module'

    if ( model_type == 2 ) then
       nexclude = 3  ! 1-4 hard sphere interactions
    else
       nexclude = 4  ! no 1-4 hard sphere interactions
    end if

    allocate(chain_created(1:nchains,1:nboxes),stat=ierr(1))
    if (any(ierr/=0)) stop 'Error allocating memory in alkane module'

    ! Chains have not yet been created
    chain_created(1:nchains,1:nboxes) = .false.

    ! Neighbour list
    allocate(list(1:nbeads*nchains*100,1:nboxes)   ,stat=ierr(1))
    allocate(startinlist(1:nbeads*nchains,1:nboxes),stat=ierr(2))
    allocate(endinlist(1:nbeads*nchains,1:nboxes)  ,stat=ierr(3))
    if (any(ierr/=0)) stop 'Error allocating neighbour list arrays'
  
    return

  end subroutine alkane_init

  subroutine alkane_destroy() bind(c)
    !-------------------------------------------------------------------------!
    ! Performs a clean shutdown of the alkane module                          !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use box, only : use_link_cells
    implicit none
    
    integer(kind=it),dimension(2) :: ierr = 0

    deallocate(Rchain,chain_created,list,startinlist,endinlist,stat=ierr(1))
    if (use_link_cells) deallocate(head_of_cell,linked_list,stat=ierr(2))
    if (any(ierr/=0)) stop 'Error releasing memory in alkane module'

  end subroutine alkane_destroy

  subroutine alkane_translate_chain(ichain,ibox,new_boltz) bind(c)
    !-------------------------------------------------------------------------!
    ! Implements an MC trial move in which the specified chain is translated  !
    ! by a random vector. The new Boltzmann factor after the trial move is    !
    ! returned as new_boltz. The old Boltzmann factor will always be one.     !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use random, only : random_uniform_random
    implicit none
    integer(kind=it),intent(in)  :: ichain,ibox
    real(kind=dp),intent(out)    :: new_boltz
    real(kind=dp),dimension(3)   :: dr
    real(kind=dp)                :: tboltz
    integer(kind=it) :: ibead

    ! generate random move
    dr(1) = 2.0_dp*random_uniform_random() - 1.0_dp
    dr(2) = 2.0_dp*random_uniform_random() - 1.0_dp
    dr(3) = 2.0_dp*random_uniform_random() - 1.0_dp

    dr = dr*mc_dr_max
    
    ! translate the chain
    do ibead = 1,nbeads
       Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) + dr(:)
    end do

    tboltz    = alkane_chain_inter_boltz(ichain,ibox)
    new_boltz = tboltz 

    return


  end subroutine alkane_translate_chain

  subroutine alkane_rotate_chain(ichain,ibox,new_boltz,quat) bind(c)
    !-------------------------------------------------------------------------!
    ! Implements an MC trial move in which the specified chain is rotated     !
    ! about the first bead. The new Boltzmann factor after the trial move     !
    ! is returned as new_boltz. The old Boltzmann factor will always be one.  !
    ! Changed from rotation about Centre of Mass to about first bead          !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use random, only     : random_uniform_random
    use quaternion, only : quat_axis_angle_to_quat,quat_conjugate_q_with_v 
    implicit none
    integer(kind=it),intent(in)             :: ichain,ibox
    real(kind=dp),intent(out)               :: new_boltz
    real(kind=dp),dimension(3)              :: axis,first
    real(kind=dp),dimension(4),intent(out)  :: quat
    real(kind=dp)                           :: theta

    integer(kind=it) :: ibead

    ! generate random rotation axis
    axis(1) = 2.0_dp*random_uniform_random() - 1.0_dp
    axis(2) = 2.0_dp*random_uniform_random() - 1.0_dp
    axis(3) = 2.0_dp*random_uniform_random() - 1.0_dp
    axis(:) = axis(:)/sqrt(dot_product(axis,axis))

    ! generate random rotation angle
    theta = (2.0_dp*random_uniform_random() - 1.0_dp) * mc_dt_max

    call quat_axis_angle_to_quat(axis,theta,quat)

    ! chain center of mass
   ! rcom(:) = 0.0_dp
    !do ibead = 1,nbeads
     !  rcom(:) = rcom(:) + Rchain(:,ibead,ichain,ibox)
    !end do
    !rcom(:) = rcom(:)/real(nbeads,kind=dp)

    !first bead 
    first(:) = Rchain(:,1,ichain,ibox)

    !move so first bead is at the origin, rotate, move back
    do ibead = 1,nbeads
       Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) - first(:)
       Rchain(:,ibead,ichain,ibox) = quat_conjugate_q_with_v(quat,Rchain(:,ibead,ichain,ibox))
       Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) + first(:)
  end do

    ! rotate the chain
!    do ibead = 1,nbeads
       !Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) - rcom(:)
       !Rchain(:,ibead,ichain,ibox) = quat_conjugate_q_with_v(quat,Rchain(:,ibead,ichain,ibox))
       !Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) + rcom(:)
 !   end do

    new_boltz = alkane_chain_inter_boltz(ichain,ibox)

    return

  end subroutine alkane_rotate_chain

  subroutine alkane_box_resize(pressure,ibox,acc_prob,reset) bind(c)
    !-------------------------------------------------------------------------!
    ! Implements an MC trial move in which the size (and possibly shape) of   !
    ! the simulation box is altered by a random amount. The ratio  of new/old !
    ! Boltmann factors (i.e. the acceptance probability) is returned as       !
    ! new_boltz. The position of the FIRST BEAD is constant in fractional     !
    ! coordinates is constant during the move. Note that any torsional        !
    ! potential does not enter into the acceptance criteria.                  !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    use constants,only : invPi
    use box,      only : hmatrix,recip_matrix,box_construct_link_cells, &
         box_update_recipmatrix,box_compute_volume,isotropic
    use random,   only : random_uniform_random
    implicit none

    real(kind=dp),intent(in)     :: pressure
    integer(kind=it),intent(in)  :: ibox
    integer(kind=it),intent(in)  :: reset
    real(kind=dp),intent(out)    :: acc_prob
!    real(kind=dp),dimension(3)   :: oldcom,comchain,tmpcom
    real(kind=dp),dimension(3)   :: first,frac_first,first_chain,delta_first
    real(kind=dp),dimension(3,3),save :: old_hmatrix,new_hmatrix,delta_hmatrix
    real(kind=dp) :: old_volume,new_volume,delta_vol,x

    integer(kind=it) :: ichain,ibead,jdim,idim

    if ( reset==1 ) then

       do ichain = 1,nchains

          first(:) = 0.0_dp
          ibead = 1
          first(:) =  Rchain(:,ibead,ichain,ibox)

          
          ! Compute fractional first bead position using the current recip_matrix          
          frac_first(1) = recip_matrix(1,1,ibox)*first(1) + &
                         recip_matrix(2,1,ibox)*first(2) + &
                         recip_matrix(3,1,ibox)*first(3)
          frac_first(2) = recip_matrix(1,2,ibox)*first(1) + &
                         recip_matrix(2,2,ibox)*first(2) + &
                         recip_matrix(3,2,ibox)*first(3)  
          frac_first(3) = recip_matrix(1,3,ibox)*first(1) + &
                         recip_matrix(2,3,ibox)*first(2) + &
                         recip_matrix(3,3,ibox)*first(3) 

          frac_first = frac_first*0.5_dp*invPi 


          ! Scale to the previous cell
          first_chain(1) = old_hmatrix(1,1)*frac_first(1) + &
                           old_hmatrix(1,2)*frac_first(2) + &
                           old_hmatrix(1,3)*frac_first(3)

          first_chain(2) = old_hmatrix(2,1)*frac_first(1) + &
                           old_hmatrix(2,2)*frac_first(2) + &
                           old_hmatrix(2,3)*frac_first(3)

          first_chain(3) = old_hmatrix(3,1)*frac_first(1) + &
                           old_hmatrix(3,2)*frac_first(2) + &
                           old_hmatrix(3,3)*frac_first(3)

          delta_first(:) = first_chain(:) - first(:)

          !write(*,'(3F15.6)')delta_first

          do ibead = 1,nbeads
             Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) + delta_first(:)
          end do

       end do

       hmatrix(:,:,ibox) = old_hmatrix(:,:)
       call box_update_recipmatrix(ibox)

       ! Restoring the old configuration, acc_prob irrelevant
       acc_prob = 1.0_dp


    else


       old_hmatrix  = hmatrix(:,:,ibox)
       old_volume   = box_compute_volume(ibox)
       delta_vol    = 0.0_dp

       if (isotropic) then

          ! Uniform scaling of box dimensions
          new_volume  = old_volume+(2.0_dp*random_uniform_random()-1.0_dp)*mc_dv_max
          delta_vol   = new_volume - old_volume
          new_hmatrix = old_hmatrix*(new_volume/old_volume)**(1.0_dp/3.0_dp)
          hmatrix(:,:,ibox) = new_hmatrix

       else

          ! Tweak a random sub-diagonal element         
          x = random_uniform_random()
          idim = int(x*3.0_dp)+1
          x = random_uniform_random()
          jdim = int(x*real(idim,kind=dp))+1

          x = random_uniform_random()

          delta_hmatrix = 0.0_dp
          delta_hmatrix(jdim,idim) =  (2.0_dp*x-1.0_dp)*mc_dv_max

!!$          write(*,'(3F15.6)')delta_hmatrix(:,1)
!!$          write(*,'(3F15.6)')delta_hmatrix(:,2)
!!$          write(*,'(3F15.6)')delta_hmatrix(:,3)
!!$          write(*,*)

          new_hmatrix       = old_hmatrix + delta_hmatrix
          hmatrix(:,:,ibox) = new_hmatrix(:,:)
          new_volume        = box_compute_volume(ibox)
          delta_vol         = new_volume - old_volume
          

       end if

       do ichain = 1,nchains

          first(:) = 0.0_dp
          ibead = 1
          first(:) =  Rchain(:,ibead,ichain,ibox)
          
          ! Compute fractional first bead position using the current recip_matrix          
          frac_first(1) = recip_matrix(1,1,ibox)*first(1) + &
                          recip_matrix(2,1,ibox)*first(2) + &
                          recip_matrix(3,1,ibox)*first(3)
          frac_first(2) = recip_matrix(1,2,ibox)*first(1) + &
                          recip_matrix(2,2,ibox)*first(2) + &
                          recip_matrix(3,2,ibox)*first(3)  
          frac_first(3) = recip_matrix(1,3,ibox)*first(1) + &
                          recip_matrix(2,3,ibox)*first(2) + &
                          recip_matrix(3,3,ibox)*first(3) 

          frac_first = frac_first*0.5_dp*invPi 

          ! Scale to the new cell
          first_chain(1) = hmatrix(1,1,ibox)*frac_first(1) + &
                           hmatrix(1,2,ibox)*frac_first(2) + &
                           hmatrix(1,3,ibox)*frac_first(3)

          first_chain(2) = hmatrix(2,1,ibox)*frac_first(1) + &
                           hmatrix(2,2,ibox)*frac_first(2) + &
                           hmatrix(2,3,ibox)*frac_first(3)

          first_chain(3) = hmatrix(3,1,ibox)*frac_first(1) + &
                           hmatrix(3,2,ibox)*frac_first(2) + &
                           hmatrix(3,3,ibox)*frac_first(3)

          delta_first(:) = first_chain(:) - first(:)

          !write(*,'(3F15.6)')delta_first

          do ibead = 1,nbeads
             Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) + delta_first(:)
          end do

       end do

       !stop
       call box_update_recipmatrix(ibox)

       ! This is a trial move and we need to check for overlaps
       ! if the box has shrunk
       acc_prob = exp(-pressure*delta_vol + real(nchains,kind=dp)*log(new_volume/old_volume))

       do ichain = 1,nchains
          if ( alkane_chain_inter_boltz(ichain,ibox) < tiny(1.0_dp) ) then
             acc_prob = 0.0_dp
             return
          end if
       end do

    end if

    call box_construct_link_cells(ibox)
    call alkane_construct_linked_lists(ibox)

    return

  end subroutine alkane_box_resize

  subroutine alkane_box_scale(ibox,scaleA,scaleB,scaleC)  bind(c)
    !-------------------------------------------------------------------------!
    ! Implements an MC trial move in which the size (and possibly shape) of   !
    ! the simulation box is altered by a random amount. The ratio  of new/old !
    ! Boltmann factors (i.e. the acceptance probability) is returned as       !
    ! new_boltz. The first bead position of each chain in fractional          !
    ! coordinates is constant during the move. Note that any torsional        !
    ! potential does not enter into the acceptance criteria.                  !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    use constants,only : invPi
    use box,      only : hmatrix,recip_matrix,box_construct_link_cells, &
                         box_update_recipmatrix,box_compute_volume
    use random,   only : random_uniform_random
    implicit none

    integer(kind=it),intent(in)  :: ibox
    real(kind=dp),intent(in)     :: scaleA,scaleB,scaleC
!    real(kind=dp),dimension(3)   :: oldcom,comchain,tmpcom
    real(kind=dp),dimension(3)   :: first,frac_first,first_chain,delta_first
    real(kind=dp),dimension(3,3) :: old_hmatrix,new_hmatrix
    real(kind=dp) :: old_volume

    integer(kind=it) :: ichain,ibead


    ! Uniform scaling of box dimensions
    old_hmatrix = hmatrix(:,:,ibox)
    old_volume  = box_compute_volume(ibox)

    new_hmatrix(:,1) = old_hmatrix(:,1)*scaleA
    new_hmatrix(:,2) = old_hmatrix(:,2)*scaleB
    new_hmatrix(:,3) = old_hmatrix(:,3)*scaleC
    
    hmatrix(:,:,ibox)  = new_hmatrix(:,:)

    do ichain = 1,nchains

       ibead = 1
       first(:) = Rchain(:,ibead,ichain,ibox)

       ! Compute fractional com position using the current recip_matrix          
       frac_first(1) = recip_matrix(1,1,ibox)*first(1) + &
                   recip_matrix(2,1,ibox)*first(2) + &
                   recip_matrix(3,1,ibox)*first(3)
       frac_first(2) = recip_matrix(1,2,ibox)*first(1) + &
                   recip_matrix(2,2,ibox)*first(2) + &
                   recip_matrix(3,2,ibox)*first(3)  
       frac_first(3) = recip_matrix(1,3,ibox)*first(1) + &
                   recip_matrix(2,3,ibox)*first(2) + &
                   recip_matrix(3,3,ibox)*first(3) 

       frac_first = frac_first*0.5_dp*invPi 

       ! Scale to the new cell
       first_chain(1) = hmatrix(1,1,ibox)*frac_first(1) + &
                        hmatrix(1,2,ibox)*frac_first(2) + &
                        hmatrix(1,3,ibox)*frac_first(3)
                                
       first_chain(2) = hmatrix(2,1,ibox)*frac_first(1) + &
                        hmatrix(2,2,ibox)*frac_first(2) + &
                        hmatrix(2,3,ibox)*frac_first(3)
                                
       first_chain(3) = hmatrix(3,1,ibox)*frac_first(1) + &
                        hmatrix(3,2,ibox)*frac_first(2) + &
                        hmatrix(3,3,ibox)*frac_first(3)

       delta_first(:) = first_chain(:) - first(:)

       do ibead = 1,nbeads
          Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox ) + delta_first(:)
       end do

    end do

    call box_update_recipmatrix(ibox)
    call box_construct_link_cells(ibox)
    call alkane_construct_linked_lists(ibox)

    return

  end subroutine alkane_box_scale


  subroutine alkane_bond_rotate(ichain,ibox,new_boltz,ia,angle) bind(c)
    !-------------------------------------------------------------------------!
    ! Selects a random dihedral angle on the selected chain and alters it by  !
    ! a random angle. The Boltzmann factor after the move is returned as      !
    ! new_boltz. Small ( within a torsional P.E. well ) and large (between    !
    ! torsional wells) moves are attempted with equal probability.            !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    use constants,        only : Pi
    use random,           only : random_uniform_random
    use quaternion,       only : quat_axis_angle_to_quat,quat_conjugate_q_with_v
    implicit none
    integer(kind=it),intent(in) :: ichain,ibox
    real(kind=dp),intent(out)   :: new_boltz

    real(kind=dp),dimension(3)  :: axis,r12,r23,r34
    real(kind=dp),dimension(4)  :: quat

    real(kind=dp),intent(out) :: angle
    real(kind=dp)             :: xi

    integer(kind=it)             :: ibead
    integer(kind=it),intent(out) :: ia

    if (nbeads<4) stop 'Called alkane_bond_rotate with nbeads < 4'

    ! select a bead at random from 1 to nbeads-3
    ia = int(real(nbeads-3,kind=dp)*random_uniform_random()) + 1
    ia = min(ia,nbeads-3)

    ! we rotate about the vector r23
    axis(:) = Rchain(:,ia+2,ichain,ibox) - Rchain(:,ia+1,ichain,ibox)
    axis(:) = axis(:)/sqrt(dot_product(axis,axis))

    ! generate a random rotation
    angle = mc_dh_max*(2.0_dp*random_uniform_random() - 1.0_dp)

    ! Flip between square-well torsions if using model IV
    if ( model_type==4 ) then
       xi = random_uniform_random()
       if ( xi < 0.5_dp ) then ! flip between basins
          angle = sign(2.0_dp*Pi/3.0_dp,angle) + angle 
       end if
    end if

    ! get quaternion
    call quat_axis_angle_to_quat(axis,angle,quat)

    ! apply rotation
    do ibead = ia+3,nbeads
       Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) - Rchain(:,ia+2,ichain,ibox)
       Rchain(:,ibead,ichain,ibox) = quat_conjugate_q_with_v(quat,Rchain(:,ibead,ichain,ibox))
       Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox) + Rchain(:,ia+2,ichain,ibox)
    end do

    ! compute Boltzmann factor for acceptance
    r12(:) = Rchain(:,ia+1,ichain,ibox) - Rchain(:,ia,ichain,ibox)
    r23(:) = axis(:)*L
    r34(:) = Rchain(:,ia+3,ichain,ibox) - Rchain(:,ia+2,ichain,ibox)
    new_boltz = alkane_dihedral_boltz(r12,r23,r34)

    ! return now if already zero
    if (new_boltz<tiny(1.0_dp)) return

    ! check non-bonded interactions (including 1-4 for model II)
    do ibead = ia+3,nbeads
       new_boltz = new_boltz*alkane_nonbonded_boltz(ibead,ichain,ibox,Rchain(:,ibead,ichain,ibox))
       if (new_boltz<tiny(1.0_dp)) return
    end do

    return

  end subroutine alkane_bond_rotate

  function alkane_dihedral_boltz(b1,b2,b3) bind(c)
    !-------------------------------------------------------------------------!
    ! Return  the multiplicative contribution to the Boltzmann factor arising !
    ! from the torsional potential acting on the dihedral potentail defined   !
    ! by the three input vectors.                                             !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    use constants, only : Pi
    implicit none
    real(kind=dp),dimension(3),intent(in) :: b1,b2,b3
    real(kind=dp),dimension(3)            :: t1,t2
    real(kind=dp) :: arg1,arg2,angle
    real(kind=dp) :: alkane_dihedral_boltz

    t1 = cross_product(b2,b3)
    t2 = cross_product(b1,b2)

    arg1 = L*dot_product(b1,t1)
    arg2 =   dot_product(t2,t1)

    angle = atan2(arg1,arg2)

    ! convert to representation used in these models
    angle = -sign(Pi-abs(angle),angle)

    select case (model_type)

       case(1)

          ! no torsional potential
          alkane_dihedral_boltz = 1.0_dp
          return

       case(2)

          ! check for 1-4 overlap
          t1 = b1 + b2 + b3
          if ( dot_product(t1,t1) < sigma*sigma ) then
             alkane_dihedral_boltz = 0.0_dp
          else
             alkane_dihedral_boltz = 1.0_dp
          end if
          return

       case(3)
          alkane_dihedral_boltz = 0.0_dp
          stop 'Not implemented'
       case(4)

          ! Model 4 torsion potential. Allows 17.4 degrees
          ! either side of zero (trans), or 10 degrees either
          ! side of +/- 120 (gauche)
          if ( abs(angle) < Pi*17.4_dp/180.0_dp ) then
             alkane_dihedral_boltz = 1.0_dp
             !print*,"Accepting a new dihedral of ",angle*180.0_dp/Pi
             return
          elseif ( abs(abs(angle)-120.0_dp*Pi/180.0_dp) < Pi*10.0_dp/180.0_dp ) then
             alkane_dihedral_boltz = 1.0_dp
             !print*,"Accepting a new dihedral of ",angle*180.0_dp/Pi
             return
          else
             alkane_dihedral_boltz = 0.0_dp
             !print*,"rejecting a new dihedral of ",angle*180.0_dp/Pi
             return
          end if
          
       case default
          alkane_dihedral_boltz = 0.0_dp
             
    end select


  contains
          
    function cross_product(a,b)
      !-------------------------------------------------------------------------!
      ! Does exactly what is says on the tin.                                   !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none
      real(kind=dp),dimension(3),intent(in) :: a,b
      real(kind=dp),dimension(3) :: cross_product

      cross_product(1) = a(2)*b(3) - a(3)*b(2)
      cross_product(2) = -(a(1)*b(3)-a(3)*b(1))
      cross_product(3) = a(1)*b(2)-b(1)*a(2)

    end function cross_product

  end function alkane_dihedral_boltz

  subroutine alkane_grow_chain(ichain,ibox,rb_factor,new_conf,ifail) bind(c)
    !-------------------------------------------------------------------------!
    ! (Re)grows chain ichain. Call with ichain,rb_factor,.true. for each      !
    ! chain after calling alkane_init. Subsequently is used in CBMC moves.    !
    ! Calling with new_conf=.false. will choose a random bead and compute the !
    ! Rosenbluth factor for a subchain (random first_bead to nbeads) without  !
    ! modifying the current chain configuration. Calling with new_conf=.true. !
    ! will generate a new chain and Rosenbluth factor, replacing the current  !
    ! chain configuration from the same first_bead to nbeads.                 !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use box,        only : hmatrix,pbc
    use constants,  only : Pi
    use random,     only : random_unit_vector,random_uniform_random
    use quaternion, only : quat_get_minimum_arc_q,quat_conjugate_q_with_v, &
                           quat_axis_angle_to_quat
    implicit none

    integer(kind=it),intent(in)  :: ichain,ibox   ! chain number to grow/regrow
    real(kind=ep),intent(out)    :: rb_factor     ! Rosenbluth factor for chain 
    integer(kind=it),intent(in)  :: new_conf      ! old or new configuration?
    integer(kind=it),intent(out) :: ifail 

   ! Dihedral / angle calculation
    real(kind=dp),dimension(3) :: r12,r23,r34,tmpvect,axis
    real(kind=dp),dimension(4) :: quat
    real(kind=dp)              :: dih,theta

    ! Trial segments and weights
    real(kind=dp),allocatable,dimension(:,:) :: rtrial
    real(kind=dp),allocatable,dimension(:)   :: wtrial
    real(kind=dp),allocatable,dimension(:)   :: wset

    ! Loop counters / error flags
    integer(kind=it) :: ierr,n,j,ib,jl

    ! Rosenbluth variables
    real(kind=dp) :: wsum

    ! Bead from which to start growth
    integer(kind=it),save :: first_bead

    ! Allocate memory for trial segments
    allocate(rtrial(1:3,1:ktrial),stat=ierr)
    if ( ierr/=0 ) stop 'Error allocating Rtrial'

    allocate(wtrial(1:ktrial),stat=ierr)
    if ( ierr/=0 ) stop 'Error allocating wtrial'

    ! Set bead from which to (re)grow
    if ( .not.chain_created(ichain,ibox) ) then
       first_bead = 1  ! grow whole chain from scratch
       !write(0,'("alkane : Creating new chain from first bead")')
    elseif (new_conf==0) then
       first_bead = int(random_uniform_random()*real(max_regrow,kind=dp)) + 1 ! integer between 1 and max_regrow
       if (first_bead>max_regrow) first_bead = max_regrow
       first_bead = nbeads - max_regrow + first_bead                          ! integer between (nbeads-max_regrow + 1) and nbeads
       !write(0,'("alkane : Regrowing an old chain from bead", I5)')first_bead
    elseif (new_conf==1) then
       !write(0,'("alkane : Growing   an new chain from bead", I5)')first_bead
    end if

    allocate(wset(first_bead:nbeads),stat=ierr)
    if (ierr/=0) stop 'Error allocating wset'

    ! Set loop counters and accumulators
    if (new_conf==1) then
       ! All trial segments/weights are new
       jl   = 1
       wset = 0.0_dp
    else
       ! The first trial segment at each bead is the old config
       jl = 2
       do ib = first_bead,nbeads
          wset(ib) = alkane_nonbonded_boltz(ib,ichain,ibox,Rchain(:,ib,ichain,ibox))
       end do
       !write(0,'("alkane : Computing Rosenbluth factor for old configuration")')
    end if

    ! Loop over beads
    rb_factor = 1_ep       
    do ib = first_bead,nbeads

       !write(0,'("alkane : Loop at bead number ", I5)')ib

       !======================================================!
       ! First bead                                           !
       !======================================================!
       if ( ib==1 ) then

          !write(0,'("alkane : At first bead")')
          if (new_conf==1) then
             Rchain(1,1,ichain,ibox) = random_uniform_random()
             Rchain(2,1,ichain,ibox) = random_uniform_random()
             Rchain(3,1,ichain,ibox) = random_uniform_random()
             Rchain(:,1,ichain,ibox) = matmul(hmatrix(:,:,ibox),Rchain(:,1,ichain,ibox))
          end if

          rb_factor = real(alkane_nonbonded_boltz(1,ichain,ibox,Rchain(:,1,ichain,ibox)),kind=ep)

          !write(0,'(I5,3F15.6)')ib,Rchain(:,ib,ichain,ibox)

       !======================================================!
       ! Second bead                                          !
       !======================================================!
       elseif ( ib==2 ) then

          !write(0,'("alkane : At second bead")')

          wsum = wset(ib)
          do j = jl,ktrial
             ! All trial segments proposed with equal probability
             rtrial(:,j) = Rchain(:,1,ichain,ibox) + L*random_unit_vector()
             wtrial(j)   = alkane_nonbonded_boltz(2,ichain,ibox,rtrial(:,j))
             !print*,wtrial(j)
             wsum = wsum + wtrial(j) 
             !write(0,'(I5,3F15.6)')ib,rtrial(:,j)
          end do
          !stop

          rb_factor = rb_factor*real(wsum,kind=ep)
          if (new_conf==1) then
             call select_next_segment(n,ifail)         
             if (ifail/=0) return ! Rosenbluth factor is zero
             Rchain(:,2,ichain,ibox) = rtrial(:,n)
          end if



       !======================================================!
       ! Third bead                                           !
       !======================================================!
       elseif ( ib==3 ) then

          !write(0,'("alkane : At third bead")')

          r12 = Rchain(:,2,ichain,ibox) - Rchain(:,1,ichain,ibox)
          wsum = wset(ib)
          do j = jl,ktrial

             ! Position of third bead if in x-y plane
             rtrial(1,j) = L+L*cos(Pi-109.47_dp*Pi/180.0_dp)
             rtrial(2,j) =   L*sin(Pi-109.47_dp*Pi/180.0_dp)
             rtrial(3,j) = 0.0_dp

             ! Rotate random angle about r12
             tmpvect(:) = (/1.0_dp,0.0_dp,0.0_dp/)
             theta = 2.0_dp*random_uniform_random() - 1.0_dp
             call quat_axis_angle_to_quat(tmpvect,theta,quat)
             rtrial(:,j) = quat_conjugate_q_with_v(quat,rtrial(:,j))

             ! Rotate onto end of r12
             tmpvect(:) = (/L,0.0_dp,0.0_dp/)
             call quat_get_minimum_arc_q(tmpvect,r12,quat)
             rtrial(:,j) = quat_conjugate_q_with_v(quat,rtrial(:,j))      
             rtrial(:,j) = Rchain(:,1,ichain,ibox) + rtrial(:,j)

             wtrial(j)   = alkane_nonbonded_boltz(3,ichain,ibox,rtrial(:,j))
             wsum = wsum + wtrial(j) 

             !write(0,'(I5,3F15.6)')ib,rtrial(:,j)

          end do

          rb_factor = rb_factor*real(wsum,kind=ep)
          
          if (new_conf==1) then
             call select_next_segment(n,ifail)         
             if (ifail/=0) return ! Rosenbluth factor is zero
             Rchain(:,3,ichain,ibox) = rtrial(:,n)
          end if

          !write(0,'(I5,3F15.6)')new_conf,Rchain(:,ib,chain,box)

       else

          !======================================================!
          ! Fourth and subsequent beads                          !
          !======================================================!

         ! write(0,'("alkane : At fourth bead")')

          r12  = Rchain(:,ib-2,ichain,ibox) - Rchain(:,ib-3,ichain,ibox)
          r23  = Rchain(:,ib-1,ichain,ibox) - Rchain(:,ib-2,ichain,ibox)
          axis = r23/sqrt(dot_product(r23,r23))
          wsum = wset(ib)
          do j = jl,ktrial

             if ( model_type == 2 ) then

                ! Keep generating random dihedrals until we
                ! get one with no 1-4 overlap
                do

                   dih = alkane_random_dihedral()

                   ! Rotate r12 into r34 about r23
                   call quat_axis_angle_to_quat(axis,dih,quat)
                   r34 = quat_conjugate_q_with_v(quat,r12)
                   if ( dot_product(r12+r23+r34,r12+r23+r34) > sigma*sigma ) exit

                end do

             else

                ! Generate a random dihedral with probability
                ! given by the intramolecular pdf only.
                dih = alkane_random_dihedral()
                !write(0,'("Random dihedral angle =     ",F15.6)')dih*180.0_dp/Pi

                ! Rotate r12 into r34 about r23
                call quat_axis_angle_to_quat(axis,dih,quat)
                r34 = quat_conjugate_q_with_v(quat,r12)

             end if

             rtrial(:,j) = Rchain(:,ib-1,ichain,ibox) + r34
             !call alkane_check_dihedral(r12,r23,r34)

             wtrial(j)   = alkane_nonbonded_boltz(ib,ichain,ibox,rtrial(:,j))
             wsum = wsum + wtrial(j) 

             !write(0,'(I5,3F15.6)')ib,rtrial(:,j)

          end do

          rb_factor = rb_factor*real(wsum,kind=ep)

          if (new_conf==1) then
             call select_next_segment(n,ifail)         
             if (ifail/=0) return ! Rosenbluth factor is zero
             Rchain(:,ib,ichain,ibox) = rtrial(:,n)
          end if

          !write(0,'(I5,3F15.6)')new_conf,Rchain(:,ib,chain,box)

       end if

    end do

    ! This chain now exists if it didn't before
    chain_created(ichain,ibox) = .true.


    deallocate(rtrial,wtrial,wset,stat=ierr)
    if (ierr/=0) stop 'Error releasing memory in alkane_grow_chain'


    return

  contains


    subroutine select_next_segment(nout,ifail)
      !-------------------------------------------------------------------------!
      ! Randomly selects, according to the appropriate weight, from the parent  !
      ! set of ktrial trial segments.                                           !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none
      integer(kind=it),intent(out) :: nout,ifail
      real(kind=dp) :: sump,zeta,rws

      sump = 0.0_dp
      zeta = random_uniform_random()

      ! Trap divide by zero in case where none of the trial segments 
      ! have a non-zero weight
      if ( wsum < tiny(1.0_dp) ) then
         ifail = 1
         return
      end if

      rws  = 1.0_dp/wsum

      if ( zeta < 0.0_dp ) stop 'under'
      if ( zeta > 1.0_dp ) stop 'over'

      nout  = -1
      ifail =  0
      do j = 1,ktrial
         sump = sump + wtrial(j)
         if ( zeta < sump*rws ) then
            nout = j
            exit
         end if
      end do

      if (nout==-1) then
         ifail = 1
         write(0,'("Warning in alkane.f90 - failed to select a regrowth segment!")')
      end if

      return

    end subroutine select_next_segment

  end subroutine alkane_grow_chain

  function alkane_nonbonded_boltz(i,ichain,ibox,rbead) bind(c)
    !-------------------------------------------------------------------------!
    ! Computes the Boltzmann factor exp[-beta*U_ext(rbead)] which is either   !
    ! zero (hard-sphere overlap) or one (no hard-sphere overlap)              !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    use box, only :  box_minimum_image,use_link_cells,ncellx,ncelly,ncellz, &
                     lcellx,lcelly,lcellz,lcneigh,hmatrix,recip_matrix

    implicit none

    integer(kind=it),intent(in) :: i,ichain,ibox
    real(kind=dp),dimension(3),intent(in) :: rbead
    real(kind=dp),dimension(3)  :: rsep,sbead

    integer(kind=it) :: jchain,ix,iy,iz,icell,ni,jcell,jbead,tmpint
    real(kind=dp) :: sigma_sq
    real(kind=dp) :: rx,ry,rz
    real(kind=dp) :: sx,sy,sz

    real(kind=dp) :: alkane_nonbonded_boltz
    logical       :: overlap

    sigma_sq = sigma*sigma
    overlap  = .false.

!!$    ! Run up chain 
!!$    if ( i+nexclude < nbeads ) then 
!!$       do j = i+nexclude,nbeads
!!$          rsep(:) = rbead(:) - Rchain(:,j,ichain)
!!$          if ( dot_product(rsep,rsep) < sigma*sigma ) then
!!$             ! overlap - no point going any further
!!$             alkane_nonbonded_boltz = 0.0_dp
!!$             return
!!$          end if
!!$       end do
!!$    end if

    ! Run down previously placed beads in the current chain
    if ( i-nexclude > 0 ) then
       do jbead = i-nexclude,1,-1
          ! The compiler needs to inline this, use -ipo on Intel
          rsep(:) = box_minimum_image( Rchain(:,jbead,ichain,ibox), rbead(:),ibox )
          if ( dot_product(rsep,rsep) < sigma_sq ) then
             ! overlap - no point going any further
             alkane_nonbonded_boltz = 0.0_dp
             return
          end if
       end do
    end if

    ! Run over beads in other chains....
    if (nchains>1) then

       if (use_link_cells) then
          
          ! compute fractional coordinates sbead from rbead
          sbead(1) = recip_matrix(1,1,ibox)*rbead(1) + &
                     recip_matrix(2,1,ibox)*rbead(2) + &
                     recip_matrix(3,1,ibox)*rbead(3)
          sbead(2) = recip_matrix(1,2,ibox)*rbead(1) + &
                     recip_matrix(2,2,ibox)*rbead(2) + &
                     recip_matrix(3,2,ibox)*rbead(3)  
          sbead(3) = recip_matrix(1,3,ibox)*rbead(1) + &
                     recip_matrix(2,3,ibox)*rbead(2) + &
                     recip_matrix(3,3,ibox)*rbead(3) 

          sbead = sbead*0.5_dp*invPi 

          ! link cell containing rbead
          ix = floor(sbead(1)/lcellx(ibox))
          iy = floor(sbead(2)/lcelly(ibox))
          iz = floor(sbead(3)/lcellz(ibox))

          ix = modulo(ix,ncellx(ibox)) + 1
          iy = modulo(iy,ncelly(ibox)) + 1
          iz = modulo(iz,ncellz(ibox)) + 1

          icell = (iz-1)*ncellx(ibox)*ncelly(ibox) + (iy-1)*ncellx(ibox) + ix

          ! loop over link cells
          do ni = 1,27
             jcell   = lcneigh(ni,icell,ibox)
             jbead   = head_of_cell(1,jcell,ibox)
             jchain  = head_of_cell(2,jcell,ibox)
             do while ( jchain.ne.0 )

                !rsep(:) = box_minimum_image(Rchain(:,jbead,jchain),rbead(:))

                rx = rbead(1) - Rchain(1,jbead,jchain,ibox)
                ry = rbead(2) - Rchain(2,jbead,jchain,ibox)
                rz = rbead(3) - Rchain(3,jbead,jchain,ibox)

                sx = recip_matrix(1,1,ibox)*rx + &
                     recip_matrix(2,1,ibox)*ry + &
                     recip_matrix(3,1,ibox)*rz
                sy = recip_matrix(1,2,ibox)*rx + &
                     recip_matrix(2,2,ibox)*ry + &
                     recip_matrix(3,2,ibox)*rz  
                sz = recip_matrix(1,3,ibox)*rx + &
                     recip_matrix(2,3,ibox)*ry + &
                     recip_matrix(3,3,ibox)*rz 

                sx = sx*0.5_dp*invPi 
                sy = sy*0.5_dp*invPi
                sz = sz*0.5_dp*invPi 

                ! apply boundary conditions
                sx = sx - floor(sx+0.5_dp,kind=dp)
                sy = sy - floor(sy+0.5_dp,kind=dp)
                sz = sz - floor(sz+0.5_dp,kind=dp)

                ! scale back up
                rx = hmatrix(1,1,ibox)*sx + &
                     hmatrix(1,2,ibox)*sy + &
                     hmatrix(1,3,ibox)*sz
                                
                ry = hmatrix(2,1,ibox)*sx + &
                     hmatrix(2,2,ibox)*sy + &
                     hmatrix(2,3,ibox)*sz
                                
                rz = hmatrix(3,1,ibox)*sx + &
                     hmatrix(3,2,ibox)*sy + &
                     hmatrix(3,3,ibox)*sz

!!$                rx = rx - Lx*anint(rx*rLx)
!!$                ry = ry - Ly*anint(ry*rLy)
!!$                rz = rz - Lz*anint(rz*rLz)

                overlap = overlap.or.( ( rx*rx+ry*ry+rz*rz < sigma_sq ).and.(ichain/=jchain) )

                tmpint  = linked_list(1,jbead,jchain,ibox)
                jchain  = linked_list(2,jbead,jchain,ibox)
                jbead   = tmpint              

             end do
             if ( overlap ) then
                alkane_nonbonded_boltz = 0.0_dp
                return
             end if
          end do

       else ! search over all pairs

          if ( ichain > 1 ) then
             do jchain = 1,ichain-1
                do jbead = 1,nbeads
                   ! The compiler needs to inline this, use -ipo on Intel
                   rsep(:) = box_minimum_image( Rchain(:,jbead,jchain,ibox),rbead(:),ibox )
                   overlap = overlap.or.(dot_product(rsep,rsep) < sigma_sq)
                end do
                if ( overlap ) then
                   alkane_nonbonded_boltz = 0.0_dp
                   return
                end if
             end do
          end if

          if ( ichain < nchains ) then
             do jchain = ichain + 1,nchains
                do jbead = 1,nbeads
                   ! The compiler needs to inline this, use -ipo on Intel
                   rsep(:) = box_minimum_image( Rchain(:,jbead,jchain,ibox),rbead(:),ibox )
                   overlap = overlap.or.(dot_product(rsep,rsep) < sigma_sq)
                end do
                if ( overlap ) then
                   alkane_nonbonded_boltz = 0.0_dp
                   return
                end if
             end do
          end if

       end if

    end if

    ! if we've got this far then....
    alkane_nonbonded_boltz = 1.0_dp

    return

  end function alkane_nonbonded_boltz

  function alkane_chain_inter_boltz(ichain,ibox) bind(c)
    !-------------------------------------------------------------------------!
    ! Computes the Boltzmann factor arising from all non-bonded intermolecular!
    ! interactions with the currently selected chain. Used in whole-chain     !
    ! translation or rotation moves.                                          !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    use box,       only :  box_minimum_image,use_link_cells,ncellx,ncelly,ncellz, &
                           lcellx,lcelly,lcellz,lcneigh,hmatrix,recip_matrix
    implicit none

    integer(kind=it),intent(in) :: ichain,ibox
    real(kind=dp),dimension(3)  :: rbead
    real(kind=dp),dimension(3)  :: rsep,sbead
    real(kind=dp),dimension(9)  :: hcache

    integer(kind=it) :: j,jchain,ibead,icell,ix,iy,iz,jbead,jcell,ni,tmpint
    real(kind=dp) :: sigma_sq
    real(kind=dp) :: rx,ry,rz
    real(kind=dp) :: sx,sy,sz

    real(kind=dp) :: alkane_chain_inter_boltz
    logical       :: overlap

    sigma_sq = sigma*sigma
    overlap  = .false.

    iz = 1
    do ix = 1,3
       do iy = 1,3
          hcache(iz) = hmatrix(ix,iy,ibox) 
          iz = iz + 1
       end do
    end do

    ! Run over beads in other chains....
    if (nchains>1) then

       if ( use_link_cells ) then

  
          do ibead = 1,nbeads

             rbead(:) = Rchain(:,ibead,ichain,ibox)

             ! compute fractional coordinates sbead from rbead
             sbead(1) = recip_matrix(1,1,ibox)*rbead(1) + &
                        recip_matrix(2,1,ibox)*rbead(2) + &
                        recip_matrix(3,1,ibox)*rbead(3)
             sbead(2) = recip_matrix(1,2,ibox)*rbead(1) + &
                        recip_matrix(2,2,ibox)*rbead(2) + &
                        recip_matrix(3,2,ibox)*rbead(3)  
             sbead(3) = recip_matrix(1,3,ibox)*rbead(1) + &
                        recip_matrix(2,3,ibox)*rbead(2) + &
                        recip_matrix(3,3,ibox)*rbead(3) 

             sbead = sbead*0.5_dp*invPi 

             ! link cell containing rbead
             ix = floor(sbead(1)/lcellx(ibox))
             iy = floor(sbead(2)/lcelly(ibox))
             iz = floor(sbead(3)/lcellz(ibox))

             ix = modulo(ix,ncellx(ibox)) + 1
             iy = modulo(iy,ncelly(ibox)) + 1
             iz = modulo(iz,ncellz(ibox)) + 1

             icell = (iz-1)*ncellx(ibox)*ncelly(ibox) + (iy-1)*ncellx(ibox) + ix

             ! loop over link cells
             do ni = 1,27
                jcell   = lcneigh(ni,icell,ibox)
                jbead   = head_of_cell(1,jcell,ibox)
                jchain  = head_of_cell(2,jcell,ibox)
                do while ( jchain.ne.0 )

                   !rsep(:) = box_minimum_image(Rchain(:,jbead,jchain),rbead(:))
                   rx = rbead(1) - Rchain(1,jbead,jchain,ibox)
                   ry = rbead(2) - Rchain(2,jbead,jchain,ibox)
                   rz = rbead(3) - Rchain(3,jbead,jchain,ibox)

                   sx = recip_matrix(1,1,ibox)*rx + &
                        recip_matrix(2,1,ibox)*ry + &
                        recip_matrix(3,1,ibox)*rz
                   sy = recip_matrix(1,2,ibox)*rx + &
                        recip_matrix(2,2,ibox)*ry + &
                        recip_matrix(3,2,ibox)*rz  
                   sz = recip_matrix(1,3,ibox)*rx + &
                        recip_matrix(2,3,ibox)*ry + &
                        recip_matrix(3,3,ibox)*rz 

                   sx = sx*0.5_dp*invPi 
                   sy = sy*0.5_dp*invPi
                   sz = sz*0.5_dp*invPi 

                   ! apply boundary conditions
                   sx = sx - floor(sx+0.5_dp,kind=dp)
                   sy = sy - floor(sy+0.5_dp,kind=dp)
                   sz = sz - floor(sz+0.5_dp,kind=dp)

                   ! scale back up
                   rx = hcache(1)*sx + hcache(2)*sy + hcache(3)*sz
                   ry = hcache(4)*sx + hcache(5)*sy + hcache(6)*sz
                   rz = hcache(7)*sx + hcache(8)*sy + hcache(9)*sz
                   
!!$                   rx = hmatrix(1,1,ibox)*sx + &
!!$                        hmatrix(1,2,ibox)*sy + &
!!$                        hmatrix(1,3,ibox)*sz

                                   
!!$                   ry = hmatrix(2,1,ibox)*sx + &
!!$                        hmatrix(2,2,ibox)*sy + &
!!$                        hmatrix(2,3,ibox)*sz
                                   
!!$                   rz = hmatrix(3,1,ibox)*sx + &
!!$                        hmatrix(3,2,ibox)*sy + &
!!$                        hmatrix(3,3,ibox)*sz

                   !overlap = overlap.or.( (rx*rx+ry*ry+rz*rz < sigma_sq).and.(ichain/=jchain) )
                   overlap = ( (rx*rx+ry*ry+rz*rz < sigma_sq).and.(ichain/=jchain) )
                   if ( overlap ) then
                      alkane_chain_inter_boltz = 0.0_dp
                      return
                   end if

                   tmpint  = linked_list(1,jbead,jchain,ibox)
                   jchain  = linked_list(2,jbead,jchain,ibox)
                   jbead   = tmpint              

                end do
                !if ( overlap ) then
                !   alkane_chain_inter_boltz = 0.0_dp
                !   return
                !end if
             end do

          end do

       else

          if ( ichain > 1 ) then
             do ibead = 1,nbeads
                rbead(:) = Rchain(:,ibead,ichain,ibox)
                do jchain = 1,ichain-1
                   do j = 1,nbeads
                      ! The compiler needs to inline this, use -ipo on Intel
                      rsep(:) = box_minimum_image( Rchain(:,j,jchain,ibox),rbead(:), ibox )
                      overlap = overlap.or.(dot_product(rsep,rsep) < sigma_sq)
                   end do
                   if ( overlap ) then
                      alkane_chain_inter_boltz = 0.0_dp
                      return
                   end if
                end do
             end do
          end if

          if ( ichain < nchains ) then
             do ibead = 1,nbeads
                rbead(:) = Rchain(:,ibead,ichain,ibox)
                do jchain = ichain + 1,nchains
                   do j = 1,nbeads
                      ! The compiler needs to inline this, use -ipo on Intel
                      rsep(:) = box_minimum_image( Rchain(:,j,jchain,ibox),rbead(:),ibox )
                      overlap = overlap.or.(dot_product(rsep,rsep) < sigma_sq)
                   end do
                   if ( overlap ) then
                      alkane_chain_inter_boltz = 0.0_dp
                      return
                   end if
                end do
             end do
          end if

       end if

    end if

    ! if we've got this far then....
    alkane_chain_inter_boltz = 1.0_dp
    return

  end function alkane_chain_inter_boltz

  function alkane_random_dihedral() bind(c)
    !-------------------------------------------------------------------------!
    ! Generates a random dihedral potential to the Boltzmann distribution     !
    ! arising from a purely torsional
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : Pi
    use random,    only : random_uniform_random
    implicit none

    real(kind=dp) :: alkane_random_dihedral
    real(kind=dp) :: xi,phi

    select case (model_type)

       case(1:2)

          ! Generate uniformly on zero to 2*Pi
          ! will be rejected in calling routine
          ! if results in 1-4 overlap in model II
          xi = random_uniform_random()
          phi = xi*360.0_dp - 180.0_dp

       case(3)
          phi = 0.0_dp
          stop 'Not implemented'
       case(4)

          ! Model 4 torsion potential. Allows 17.4 degrees
          ! either side of zero (trans), or 10 degrees either
          ! side of +/- 120 (gauche)

          xi = random_uniform_random()
          if ( xi < 34.8_dp/74.8_dp ) then
             xi = random_uniform_random()
             phi = xi*34.8_dp-17.4_dp
          else
             xi = random_uniform_random()
             if ( xi < 0.5_dp ) then
                xi = random_uniform_random()
                phi =  120.0_dp + xi*20.0_dp - 10.0_dp
             else
                xi = random_uniform_random()
                phi = -120.0_dp + xi*20.0_dp - 10.0_dp
             end if
          end if

       case default
          phi = 0.0_dp

    end select

    alkane_random_dihedral = phi*Pi/180.0_dp

    return

  end function alkane_random_dihedral
    
  subroutine alkane_check_dihedral(b1,b2,b3,angle) bind(c)
    !-------------------------------------------------------------------------!
    ! Computes the dihedral angle formed as that between the two planes       !
    ! formed by b1xb2 and b2xb3 and shifts to be consistent with the          !
    ! convention employed in these models.                                    !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : Pi
    implicit none
    real(kind=dp),dimension(3),intent(in) :: b1,b2,b3
    real(kind=dp),intent(out) :: angle
    real(kind=dp),dimension(3)            :: t1,t2
    real(kind=dp) :: arg1,arg2

    t1 = cross_product(b2,b3)
    t2 = cross_product(b1,b2)

    arg1 = L*dot_product(b1,t1)
    arg2 =   dot_product(t2,t1)

    angle = atan2(arg1,arg2)*180_dp/Pi 

    ! convect to be consistent with the angle-origin define by the model
    angle = -sign(180.0_dp-abs(angle),angle)

    !write(0,'("Measured dihedral angle as : ",F15.6)')angle 

  contains

    function cross_product(a,b)
      !-------------------------------------------------------------------------!
      ! Does exactly what is says on the tin.                                   !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none
      real(kind=dp),dimension(3),intent(in) :: a,b
      real(kind=dp),dimension(3) :: cross_product

      cross_product(1) = a(2)*b(3) - a(3)*b(2)
      cross_product(2) = -(a(1)*b(3)-a(3)*b(1))
      cross_product(3) = a(1)*b(2)-b(1)*a(2)

    end function cross_product

  end subroutine alkane_check_dihedral

  subroutine alkane_check_chain_overlap(ibox,overlap) bind(c)
    !-------------------------------------------------------------------------!
    ! Sanity test for debugging. Checks if any two chains overlap. Does not   !
    ! check if a chain overlaps with itself.                                  !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use box, only : use_link_cells
    implicit none
    integer(kind=it),intent(in) :: ibox
    real(kind=dp) :: test,acc,acc_link
    integer(kind=it) :: ichain,num_chains_overlapping,num_chains_overlapping_link
    integer(kind=it),intent(out) :: overlap
    
    !write(*,*)"check chain ibox=",ibox

    if (nchains < 2) stop 'Called alkane_check_chain_overlap with one chain'
    
    ! We take one of two paths through this routine.   
    if (.not.use_link_cells) then
       
       ! Just check for overlaps 
       
       overlap = 0
       acc = 1.0_dp
       do ichain = 1,nchains
          test = alkane_chain_inter_boltz(ichain,ibox)
          acc  = acc*test
          if (test < tiny(1.0_dp) ) then
            ! write(0,'("Chain ",I5", overlaps with another chain")')ichain
          end if
       end do

       if (acc<tiny(1.0_dp)) then
          !write(0,'("Stopping")')
          !stop
          overlap = 1
       end if
        
       
    else
    
       ! Check for overlaps using link cells
       overlap = 0
       num_chains_overlapping_link = 0
       acc_link = 1.0_dp
       do ichain = 1,nchains
          test = alkane_chain_inter_boltz(ichain,ibox)
          acc_link  = acc_link*test
          if (test < tiny(1.0_dp) ) then
             !write(0,'("Chain ",I5", overlaps with another chain (computed using link cells)")')ichain
             num_chains_overlapping_link = num_chains_overlapping_link + 1
          end if
       end do
        
       ! Temporarily disable link cells
       use_link_cells = .false.

       overlap = 0
       num_chains_overlapping = 0
       acc = 1.0_dp
       do ichain = 1,nchains
          test = alkane_chain_inter_boltz(ichain,ibox)
          acc  = acc*test
          if (test < tiny(1.0_dp) ) then
             !write(0,'("Chain ",I5", overlaps with another chain (computed without link cells)")')ichain
             num_chains_overlapping = num_chains_overlapping + 1
          end if
       end do

       ! Re-enable link cells
       use_link_cells = .true.

       if (acc<tiny(1.0_dp)) then
          !write(0,'("Stopping")')
          !stop
          overlap = 1
       end if

       ! Check for consistency in the overall Boltzmann factor of the cell
       if ( abs(acc-acc_link) > epsilon(1.0_dp) ) then
          write(0,'("Error : Link cell calculation of overlaps disagrees with brute force method!")')
          write(0,'("Error : Boltzmann factors are ",F12.8," (link cell) ",F12.8," (brute force).")')acc_link,acc
       end if

       ! Check for consitency in the number of chains identified as overlapping with others
       if ( num_chains_overlapping /= num_chains_overlapping_link ) then
          write(0,'("Error : Link cell calculation of overlaps disagrees with brute force method!")')
          write(0,'("Error : Num. overlapping chains :",I5," (link cell) ",I5," (brute force).")') &
               num_chains_overlapping_link,num_chains_overlapping
       end if


    end if

    return

  end subroutine alkane_check_chain_overlap

  subroutine alkane_check_chain_geometry(ichain,ibox,violated) bind(c)
    !-------------------------------------------------------------------------!
    ! Sanity test for debugging. Checks bond lengths within a chain           !
    ! violated = 0 for no problems					      !
    ! violated = 1 if there are problems				      !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : Pi
    use box
    implicit none
    integer(kind=it),intent(in)  :: ichain
    integer(kind=it),intent(in)  :: ibox
    integer(kind=it),intent(out) :: violated 
    
    real(kind=dp),dimension(3) :: rsep,r12,r23,r34
    integer(kind=it) :: ibead,jbead
    real(kind=dp) :: boltzf,angle


    violated = 0 ! No problems found yet

    !==============================================!
    ! Check that all bond lengths are equal to L   !
    !==============================================!
    do ibead = 1,nbeads-1
       jbead = ibead+1
       rsep(:) = box_minimum_image( Rchain(:,jbead,ichain,ibox), Rchain(:,ibead,ichain,ibox), ibox )
       if ( dot_product(rsep,rsep) - L*L > 0.00001_dp ) then
          violated = 1
          write(0,'("Found a bond length of ",E15.6," between beads ",I5," and ",I5)') &
               sqrt(dot_product(rsep,rsep)),ibead,jbead
       end if
    end do

    if ( violated == 1) then
       write(0,'("Bond length violation for chain ",I5)')ichain
       write(0,'("Other chains may be affected")')
       return
    end if


    !=================================================!
    ! Check that bond angles are equal to 109.47 deg  !
    !=================================================!
    do ibead=1,nbeads-2

       r12(:) = box_minimum_image(Rchain(:,ibead,  ichain,ibox),Rchain(:,ibead+1,ichain,ibox), ibox)
       r23(:) = box_minimum_image(Rchain(:,ibead+2,ichain,ibox),Rchain(:,ibead+1,ichain,ibox), ibox)

       r12(:) = r12(:)/sqrt(dot_product(r12,r12))
       r23(:) = r23(:)/sqrt(dot_product(r23,r23))

       if ( abs(acos(dot_product(r12,r23)) - 109.47_dp*Pi/180.0_dp) > 0.0002_dp ) then
          violated = 1
          write(0,'("Found a bond angle of ",F15.6," involving beads ",3I5)') &
               acos(dot_product(r12,r23))*180.0_dp/Pi,ibead,ibead+1,ibead+2
       end if

    end do

    if ( violated == 1) then
       write(0,'("Bond angle violation for chain ",I5)')ichain
       write(0,'("Other chains may be affected")')
       return
    end if


    !==============================================!
    ! Check that torsions are sensible             !
    !==============================================!
    do ibead = 1,nbeads-3

       ! compute Boltzmann factor for acceptance
       r12(:) = box_minimum_image(Rchain(:,ibead,  ichain,ibox),Rchain(:,ibead+1,ichain,ibox),ibox)
       r23(:) = box_minimum_image(Rchain(:,ibead+1,ichain,ibox),Rchain(:,ibead+2,ichain,ibox),ibox)
       r34(:) = box_minimum_image(Rchain(:,ibead+2,ichain,ibox),Rchain(:,ibead+3,ichain,ibox),ibox)
       boltzf = alkane_dihedral_boltz(r12,r23,r34)

       if (boltzf<tiny(1.0_dp)) then
          write(*,*)"Warning: Found a bad dihedral angle on chain ",ichain,"box ",ibox
          write(*,*)"botltzman factor: ",boltzf
          call alkane_check_dihedral(r12,r23,r34,angle)
          write(*,*)"with angle ",angle
          ! DQ - this should only be a warning, as rounding errors often trigger this for dihedral
          ! angles slightly outside the allowed range when operations are performed in a different
          ! order to how they are generated.
          ! violated = 1 
       end if

    end do

    if ( violated == 1) then
       write(*,'("Dihedral angle violation for chain ",I5," in box ",I5)')ichain,ibox
       write(*,'("Other chains may be affected")')
       return
    end if

    !==============================================!
    ! Check that there are no overlaps             !
    !==============================================!
    do ibead = 1,nbeads-nexclude
       do jbead = ibead+nexclude,nbeads
          !write(0,'("Checking overlap between beads ",I5," and ",I5," on chain ",I5)')ibead,jbead,ichain
          r12(:) = box_minimum_image(Rchain(:,ibead, ichain, ibox),Rchain(:,jbead,ichain, ibox), ibox)
          !print*,sqrt(dot_product(r12,r12))
          if ( dot_product(r12,r12) < sigma*sigma ) then
             violated = 1
             write(0,'("Found an overlap between beads ",I5," and ",I5," on chain ",I5, " in box ",I5)')ibead,jbead,ichain,ibox
          end if
       end do
    end do

    if ( violated == 1) then
       write(0,'("Intra-chain bead overlap violation for chain ",I5, "in box ",I5)')ichain,ibox
       write(0,'("Other chains may be affected")')
       return
    else
       violated = 0   
    end if
   
   
    return

  end subroutine alkane_check_chain_geometry

  subroutine alkane_construct_neighbour_list(ibox) bind(c)
    !-------------------------------------------------------------------------!
    ! Constructs a Verlet neighbour list. All intra-chain interactions are    !
    ! excluded from this list and must therefore be included in a seperate    !
    ! search.                                                                 !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use box
    implicit none
    integer(kind=it),intent(in) :: ibox

    integer(kind=it) :: myint,logint,ibead,jbead,ichain,jchain,k,ierr,m
    integer(kind=it) :: t1,t2,rate
    logical :: lrange
    real(kind=dp) :: nl_range_sq 
    real(kind=dp),dimension(3) :: rbead,rsep

    integer(kind=it),allocatable,dimension(:) :: advance
    logical,save  :: firstpass = .true.
    integer(kind=it),save  :: lcnv

    nl_range_sq = (2.0_dp*sigma)**2

    if ( firstpass ) then

       ! Establish the manner in which the present compiler
       ! stores logical variables as integers and establish
       ! a conversion between this and "1".
       logint = transfer(.true.,myint)
       write(*,'("!=======================================!")')
       write(*,'("! This compiler stores .true. as : ",I5,"!")')logint
       write(*,'("!=======================================!")')

       lcnv = 1/logint

       firstpass = .false.

       allocate(advance(1:nchains*nbeads),stat=ierr)
       if (ierr/=0) stop 'Error allocating advance array'

    end if

    call system_clock(count=t1)

    m = 1
    do ichain = 1,nchains
       do ibead = 1,nbeads

          rbead(:) = Rchain(:,ibead,ichain,ibox)

          k = 1
          do jchain = 1,nchains
             do jbead = 1,nbeads
                ! Compiler need to in-line this next bit
                rsep(:) = box_minimum_image( Rchain(:,jbead,jchain,ibox),rbead(:),ibox )
                !rsep(:) = Rchain(:,jbead,jchain)-rbead(:)
                !rsep(1) = rsep(1) - Lx*anint(rsep(1)*rLx)
                !rsep(2) = rsep(2) - Ly*anint(rsep(2)*rLy)
                !rsep(3) = rsep(3) - Lz*anint(rsep(3)*rLz)
                lrange  = (dot_product(rsep,rsep)<nl_range_sq).and.(ichain/=jchain) 
                advance(k) = lcnv*transfer(lrange,myint)
                k = k + 1
             end do
          end do

          k = 1
          startinlist((ichain-1)*nbeads+ibead,ibox) = m
          do jchain = 1,nchains
             do jbead = 1,nbeads
                list(m,ibox) = k
                m = m + advance(k)
                k = k + 1
             end do
          end do
          endinlist((ichain-1)*nbeads+ibead,ibox) = m - 1

       end do
    end do

    deallocate(advance,stat=ierr)
    if (ierr/=0) stop 'Error releasing memory in alkane_construct_neighbour_list'

    call system_clock(count=t2,count_rate=rate)

    write(*,'("Neighbour list constructed in : ",F15.6," seconds")')real(t2-t1,kind=ep)/real(rate,kind=ep)
    !stop

    return

  end subroutine alkane_construct_neighbour_list

  subroutine alkane_construct_linked_lists(ibox) bind(c)
    !-------------------------------------------------------------------------!
    ! Assigns particles to link cells and constructs linked lists for use in  !
    ! searching over particle pairs. Closely follows procedure in appendix F  !
    ! of Frenkel and Smit.                                                    !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    use box,       only : ncellx,ncelly,ncellz,lcellx,lcelly,lcellz,use_link_cells, &
                          recip_matrix,nboxes,maxcells
    implicit none
    integer(kind=it),intent(in) :: ibox
    integer(kind=it) :: ichain,ibead,icell,ix,iy,iz,ierr
    integer(kind=it) :: jbead,jchain
    real(kind=dp)    :: rlcellx,rlcelly,rlcellz 
    real(kind=dp),dimension(3) :: rbead,sbead

    logical :: lrebuild_all_boxes = .false.
    integer(kind=it) :: ifirst,ilast
    integer(kind=it) :: jbox
    !integer(kind=it) :: t1,t2,rate,


    if (.not.use_link_cells) return

    ! Allocate head_of_cell array using this size if not allocated
    if ( .not.allocated(head_of_cell) ) then
        allocate(head_of_cell(1:2,1:maxcells,1:nboxes),stat=ierr)
       if (ierr/=0) stop 'Error allocating head_of_cell'
    end if

    ! Check that currently allocated head_of_cell arrray can hold all cells in each box
    if (size(head_of_cell,2)/=maxcells) then
       deallocate(head_of_cell,stat=ierr)
       if (ierr/=0) stop 'Error resizing head_of_cell'
       allocate(head_of_cell(1:2,1:maxcells,1:nboxes),stat=ierr)
       if (ierr/=0) stop 'Error allocating head_of_cell'
       lrebuild_all_boxes = .true.
    end if

    ! Allocate linked_list array itself
    if ( .not.allocated(linked_list) ) then
       allocate(linked_list(1:4,0:nbeads,0:nchains,1:nboxes),stat=ierr)
       if (ierr/=0) stop 'Error allocating linked_lists'
    end if

    if ( lrebuild_all_boxes ) then
       ifirst = 1 ; ilast = nboxes
    else
       ifirst = ibox;ilast = ibox
    end if

    !call system_clock(count=t1)

    do jbox = ifirst,ilast

       ! Don't try to build link cells if cell dimensions have not
       ! been set for jbox, i.e. at initialisation.
       if ( ncellx(jbox)*ncelly(jbox)*ncellz(jbox) == 0 ) cycle

       head_of_cell(:,:,jbox) = 0
       
       rlcellx = 1.0_dp/lcellx(jbox)
       rlcelly = 1.0_dp/lcelly(jbox)
       rlcellz = 1.0_dp/lcellz(jbox)
       
       do ichain = 1,nchains
          do ibead = 1,nbeads
             
             rbead(:) = Rchain(:,ibead,ichain,jbox)
             
             ! compute fractional coordinates sbead from rbead
             sbead(1) = recip_matrix(1,1,jbox)*rbead(1) + &
                        recip_matrix(2,1,jbox)*rbead(2) + &
                        recip_matrix(3,1,jbox)*rbead(3)
             sbead(2) = recip_matrix(1,2,jbox)*rbead(1) + &
                        recip_matrix(2,2,jbox)*rbead(2) + &
                        recip_matrix(3,2,jbox)*rbead(3)  
             sbead(3) = recip_matrix(1,3,jbox)*rbead(1) + &
                        recip_matrix(2,3,jbox)*rbead(2) + &
                        recip_matrix(3,3,jbox)*rbead(3) 

             sbead = sbead*0.5_dp*invPi 
             
             ! which link cell does this particle belong to
             ix = floor(sbead(1)*rlcellx)
             iy = floor(sbead(2)*rlcelly)
             iz = floor(sbead(3)*rlcellz)
             
             ix = modulo(ix,ncellx(jbox)) + 1
             iy = modulo(iy,ncelly(jbox)) + 1
             iz = modulo(iz,ncellz(jbox)) + 1
             
             icell = (iz-1)*ncellx(jbox)*ncelly(jbox) + (iy-1)*ncellx(jbox) + ix

             ! Bead and chain index for old head of cell
             ! (both zero if this is first atom to be added)
             jbead  = head_of_cell(1,icell,jbox)
             jchain = head_of_cell(2,icell,jbox)

             ! This bead points forward to the old head of cell,
             ! or to zero if it's the first bead to be added
             linked_list(1,ibead,ichain,jbox) = jbead
             linked_list(2,ibead,ichain,jbox) = jchain

             ! ..and becomes the new head of cell
             head_of_cell(1,icell,jbox) = ibead 
             head_of_cell(2,icell,jbox) = ichain 
          
             ! jbead, jchain points backward to ibead,ichain
             ! zero'th elements will be populated here for
             ! first bead added, and ignored.
             linked_list(3,jbead,jchain,jbox) = ibead
             linked_list(4,jbead,jchain,jbox) = ichain


          end do
       end do
       
    end do ! end loop over boxes

    !call system_clock(count=t2,count_rate=rate)
    !write(*,'("Linked lists rebuilt in : ",F15.6," seconds")')real(t2-t1,kind=ep)/real(rate,kind=ep)

    return

  end subroutine alkane_construct_linked_lists

  subroutine alkane_update_linked_lists(ibead,ichain,ibox,old_pos,new_pos) bind(c)
    !-------------------------------------------------------------------------!
    ! Attempts to correct the linked lists after a bead has moved from old    !
    ! pos to new pos (possibly crossing a link-cell boundary). As we will     !
    ! need to do this for several beads at a time it is possibly faster to    !
    ! reconstruct the entire list instead....                                 !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    use box,       only : ncellx,ncelly,ncellz,lcellx,lcelly,lcellz, &
                          use_link_cells,recip_matrix
    implicit none
    integer(kind=it),intent(in) :: ibead,ichain,ibox
    real(kind=dp),dimension(3),intent(in) :: old_pos,new_pos
    real(kind=dp),dimension(3) :: sbead
    integer(kind=it) :: ix,iy,iz,ncell,ocell,jchain,jbead
    integer(kind=it) :: kbead,kchain
    !integer(kind=it) :: t1,t2,rate


    if (.not.use_link_cells) return

    ! Compute sbead from old_pos
    sbead(1) = recip_matrix(1,1,ibox)*old_pos(1) + &
               recip_matrix(2,1,ibox)*old_pos(2) + &
               recip_matrix(3,1,ibox)*old_pos(3)
    sbead(2) = recip_matrix(1,2,ibox)*old_pos(1) + &
               recip_matrix(2,2,ibox)*old_pos(2) + &
               recip_matrix(3,2,ibox)*old_pos(3)  
    sbead(3) = recip_matrix(1,3,ibox)*old_pos(1) + &
               recip_matrix(2,3,ibox)*old_pos(2) + &
               recip_matrix(3,3,ibox)*old_pos(3) 

    sbead = sbead*0.5_dp*invPi 

    ! Compute link cell number for old and new positions
    ix = floor(sbead(1)/lcellx(ibox))
    iy = floor(sbead(2)/lcelly(ibox))
    iz = floor(sbead(3)/lcellz(ibox))

    ix = modulo(ix,ncellx(ibox)) + 1
    iy = modulo(iy,ncelly(ibox)) + 1
    iz = modulo(iz,ncellz(ibox)) + 1

    ocell = (iz-1)*ncellx(ibox)*ncelly(ibox) + (iy-1)*ncellx(ibox) + ix

    ! Compute sbead from new_pos
    sbead(1) = recip_matrix(1,1,ibox)*new_pos(1) + &
               recip_matrix(2,1,ibox)*new_pos(2) + &
               recip_matrix(3,1,ibox)*new_pos(3)
    sbead(2) = recip_matrix(1,2,ibox)*new_pos(1) + &
               recip_matrix(2,2,ibox)*new_pos(2) + &
               recip_matrix(3,2,ibox)*new_pos(3)  
    sbead(3) = recip_matrix(1,3,ibox)*new_pos(1) + &
               recip_matrix(2,3,ibox)*new_pos(2) + &
               recip_matrix(3,3,ibox)*new_pos(3) 

    sbead = sbead*0.5_dp*invPi 


    ix = floor(sbead(1)/lcellx(ibox))
    iy = floor(sbead(2)/lcelly(ibox))
    iz = floor(sbead(3)/lcellz(ibox))

    ix = modulo(ix,ncellx(ibox)) + 1
    iy = modulo(iy,ncelly(ibox)) + 1
    iz = modulo(iz,ncellz(ibox)) + 1

    ncell = (iz-1)*ncellx(ibox)*ncelly(ibox) + (iy-1)*ncellx(ibox) + ix

    ! ...nothing to see here
    if (ocell==ncell) return

    !call system_clock(count=t1)

    ! remove from old cell
!!$    do jchain = 1,nchains
!!$       do jbead = 1,nbeads
!!$          if ( (linked_list(1,jbead,jchain)==ibead).and.(linked_list(2,jbead,jchain)==ichain) ) then
!!$             ! This entry used to point to ibead,ichain. Make
!!$             ! it point to where ibead,ichains entry pointed to
!!$             linked_list(:,jbead,jchain)=linked_list(:,ibead,ichain)
!!$             exit
!!$          end if
!!$       end do
!!$    end do

!!$    write(0,'("Moving bead ",2I5," from link cell ",I5," to ",I5)')ibead,ichain,ocell,ncell

    !-----------------------------------!
    ! Remove from old cell linked lists !
    !-----------------------------------!
    jbead   = head_of_cell(1,ocell,ibox)
    jchain  = head_of_cell(2,ocell,ibox)
    if ( (jbead==ibead).and.(jchain==ichain) ) then

       ! ibead,ichain was the old head of cell

       !----------------------------------!
       !   Before             After       !
       !   ------             -----       !
       !  hoc->i->k->..      hoc->k->..   !
       !----------------------------------!
       
       ! k used to be pointed to by i
       kbead  = linked_list(1,ibead,ichain,ibox)
       kchain = linked_list(2,ibead,ichain,ibox)

       ! k is new head of cell
       head_of_cell(1,ocell,ibox) = kbead
       head_of_cell(2,ocell,ibox) = kchain

       ! k points backward to nothing
       linked_list(3,kbead,kchain,ibox) = 0
       linked_list(4,kbead,kchain,ibox) = 0

    else
       
       !-------------------------------!
       !    Before          After      !
       !    ------          -----      !
       ! ..j->i->k->...   ..j->k->..   !
       !-------------------------------!
       kbead  = linked_list(1,ibead,ichain,ibox)
       kchain = linked_list(2,ibead,ichain,ibox)

       jbead  = linked_list(3,ibead,ichain,ibox)
       jchain = linked_list(4,ibead,ichain,ibox)

       linked_list(1,jbead,jchain,ibox) = kbead
       linked_list(2,jbead,jchain,ibox) = kchain

       linked_list(3,kbead,kchain,ibox) = jbead
       linked_list(4,kbead,kchain,ibox) = jchain


!!$       do 
!!$
!!$          if ( (linked_list(1,jbead,jchain)==ibead).and.(linked_list(2,jbead,jchain)==ichain) ) then
!!$
!!$             ! This entry used to point to ibead,ichain. Make
!!$             ! it point to where ibead,ichains entry pointed to
!!$             linked_list(:,jbead,jchain)=linked_list(:,ibead,ichain)
!!$             exit
!!$
!!$          end if
!!$
!!$          tmpint  = linked_list(1,jbead,jchain)
!!$          jchain  = linked_list(2,jbead,jchain)
!!$          jbead   = tmpint    
!!$
!!$          if (jchain==0) then
!!$             write(0,'("Warning : Rebuild of link cell lists forced")')
!!$             call alkane_construct_linked_lists()
!!$             return
!!$          end if
!!$
!!$       end do

    end if

    !----------------------------------!
    ! Add into new cell's linked lists !
    !----------------------------------!
    !   Before             After       !
    !   ------             -----       !
    !  hoc->j->...     hoc->i->j->.... !
    !----------------------------------!

    ! Bead and chain index for old head of cell
    jbead  = head_of_cell(1,ncell,ibox)
    jchain = head_of_cell(2,ncell,ibox)

    ! add ibead,icheain to new cell as head of cell
    linked_list(1:2,ibead,ichain,ibox) = head_of_cell(:,ncell,ibox)
    linked_list(3:4,ibead,ichain,ibox) = 0
    head_of_cell(:,ncell,ibox) = (/ibead,ichain/)

    ! jbead, jchain points backward to ibead,ichain
    linked_list(3,jbead,jchain,ibox) = ibead
    linked_list(4,jbead,jchain,ibox) = ichain 


    !call system_clock(count=t2,count_rate=rate)
    !write(*,'("Linked lists updated in : ",F15.6," seconds")')real(t2-t1,kind=ep)/real(rate,kind=ep)

    return

  end subroutine alkane_update_linked_lists

  !subroutine alkane_get_internal_overlaps(ichain,ibox,mxoverlap,noverlap,poverlap)
  subroutine alkane_get_internal_overlaps(ichain,ibox,noverlap) bind(c)
    !-------------------------------------------------------------------------!
    ! Counts the number of internal intramolecular overlaps on a single chain !
    ! in a (presumably) inactive box/lattice. Useful for LSMC calculations.   !
    ! Note that a violation of a hard torsion potential (model IV) counts as  !
    ! a single overlap.                                                       !
    !-------------------------------------------------------------------------!
    ! D.Quigley July 2011                                                     !
    !-------------------------------------------------------------------------!
    use constants, only : Pi
    use box
    implicit none
    integer(kind=it),intent(in) :: ichain    ! Chain no. to check for internal overlaps
    integer(kind=it),intent(in) :: ibox      ! Box in which this chain resides
    !integer(kind=it),intent(in) :: mxoverlap ! Max. no. of overlaps

    integer(kind=it),intent(out) :: noverlap ! number of intra-chain overlaps
    !integer(kind=it),dimension(2,mxoverlap) :: poverlap ! list of overlapping pairs

    real(kind=dp),dimension(3) :: r12,r23,r34
    integer(kind=it) :: ibead,jbead
    real(kind=dp) :: boltzf

    noverlap = 0

    !==================================================================!
    ! Check torsions, and count a out-of-bounds torsion as an overlap. !
    !==================================================================!
    do ibead = 1,nbeads-3

       ! compute Boltzmann factor for acceptance
       r12(:) = box_minimum_image(Rchain(:,ibead,  ichain,ibox),Rchain(:,ibead+1,ichain,ibox),ibox)
       r23(:) = box_minimum_image(Rchain(:,ibead+1,ichain,ibox),Rchain(:,ibead+2,ichain,ibox),ibox)
       r34(:) = box_minimum_image(Rchain(:,ibead+2,ichain,ibox),Rchain(:,ibead+3,ichain,ibox),ibox)
       boltzf = alkane_dihedral_boltz(r12,r23,r34)

       if (boltzf<tiny(1.0_dp)) then
          noverlap = noverlap + 1  ! count this as an overlap
          !poverlap(1,noverlap) = ibead    ! Store the outer atoms of this torsion angle 
          !poverlap(2,noverlap) = ibead+3  ! ..as being those which overlap.
       end if

    end do

    !=======================================================!
    ! Check for other overlaps, discounting excluded pairs. ! 
    !=======================================================!
    do ibead = 1,nbeads-nexclude
       do jbead = ibead+nexclude,nbeads
          !write(0,'("Checking overlap between beads ",I5," and ",I5," on chain ",I5)')ibead,jbead,ichain
          r12(:) = box_minimum_image(Rchain(:,ibead, ichain, ibox),Rchain(:,jbead,ichain, ibox), ibox)
          !print*,sqrt(dot_product(r12,r12))
          if ( dot_product(r12,r12) < sigma*sigma ) then
             noverlap = noverlap + 1  ! count this as an overlap
             !poverlap(1,noverlap) = ibead  ! Store the two beads which overlap
             !poverlap(2,noverlap) = jbead  
          end if
       end do
    end do

    !if (noverlap> mxoverlap) stop 'Error mxoverlap exceeded in alkane_get_internal_overlaps'

    ! Debug - plz comment out when happy
   ! write(0,*)
    !write(0,'("Found ",I5," overlaps between beads on chain ",I5," in box ",I5)')noverlap,ichain,ibox
    !write(0,*)
!!$    do iovlp = 1,noverlap
!!$       write(0,'("Bead ",I5," overlaps with bead ",I5)')poverlap(1,iovlp),poverlap(2,iovlp)
!!$    end do

    return

  end subroutine alkane_get_internal_overlaps

  !subroutine alkane_get_external_overlaps(ichain,ibox,mxoverlap,noverlap,loverlap)
  subroutine alkane_get_external_overlaps(ichain,ibox,noverlap) bind(c)
    !-------------------------------------------------------------------------!
    ! Counts the number of overlaps between ichain and all other chains, and  !
    ! returns a list of the chains which overlap with ichain. Note that       !
    ! noverlap is the number of intermolecular bead-bead overlaps, not the    !
    ! number of chain-chain overlaps.                                         !
    !-------------------------------------------------------------------------!
    ! D.Quigley July 2011                                                     !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    use box,       only :  box_minimum_image,use_link_cells,ncellx,ncelly,ncellz, &
                           lcellx,lcelly,lcellz,lcneigh,hmatrix,recip_matrix
    implicit none
    integer(kind=it),intent(in) :: ichain    ! Chain no. to check for internal overlaps
    integer(kind=it),intent(in) :: ibox      ! Box in which this chain resides
    !integer(kind=it),intent(in) :: mxoverlap ! Max. no. of overlaps

    integer(kind=it),intent(out) :: noverlap ! number of inter-chain overlaps
    !integer(kind=it),dimension(2,mxoverlap) :: loverlap ! list of overlapping beads,chains


    real(kind=dp),dimension(3) :: rbead
    real(kind=dp),dimension(3) :: rsep,sbead

    integer(kind=it) :: j,jchain,ibead,icell,ix,iy,iz,jbead,jcell,ni,tmpint
    real(kind=dp) :: sigma_sq
    real(kind=dp) :: rx,ry,rz
    real(kind=dp) :: sx,sy,sz

    sigma_sq = sigma*sigma

    noverlap = 0 ! initialise


    ! Run over beads in other chains....
    if (nchains>1) then

       if ( use_link_cells ) then

          do ibead = 1,nbeads

             rbead(:) = Rchain(:,ibead,ichain,ibox)

             ! compute fractional coordinates sbead from rbead
             sbead(1) = recip_matrix(1,1,ibox)*rbead(1) + &
                        recip_matrix(2,1,ibox)*rbead(2) + &
                        recip_matrix(3,1,ibox)*rbead(3)
             sbead(2) = recip_matrix(1,2,ibox)*rbead(1) + &
                        recip_matrix(2,2,ibox)*rbead(2) + &
                        recip_matrix(3,2,ibox)*rbead(3)  
             sbead(3) = recip_matrix(1,3,ibox)*rbead(1) + &
                        recip_matrix(2,3,ibox)*rbead(2) + &
                        recip_matrix(3,3,ibox)*rbead(3) 

             sbead = sbead*0.5_dp*invPi 

             ! link cell containing rbead
             ix = floor(sbead(1)/lcellx(ibox))
             iy = floor(sbead(2)/lcelly(ibox))
             iz = floor(sbead(3)/lcellz(ibox))

             ix = modulo(ix,ncellx(ibox)) + 1
             iy = modulo(iy,ncelly(ibox)) + 1
             iz = modulo(iz,ncellz(ibox)) + 1

             icell = (iz-1)*ncellx(ibox)*ncelly(ibox) + (iy-1)*ncellx(ibox) + ix

             ! loop over link cells
             do ni = 1,27
                jcell   = lcneigh(ni,icell,ibox)
                jbead   = head_of_cell(1,jcell,ibox)
                jchain  = head_of_cell(2,jcell,ibox)
                do while ( jchain.ne.0 )

                   !rsep(:) = box_minimum_image(Rchain(:,jbead,jchain),rbead(:))
                   rx = rbead(1) - Rchain(1,jbead,jchain,ibox)
                   ry = rbead(2) - Rchain(2,jbead,jchain,ibox)
                   rz = rbead(3) - Rchain(3,jbead,jchain,ibox)

                   sx = recip_matrix(1,1,ibox)*rx + &
                        recip_matrix(2,1,ibox)*ry + &
                        recip_matrix(3,1,ibox)*rz
                   sy = recip_matrix(1,2,ibox)*rx + &
                        recip_matrix(2,2,ibox)*ry + &
                        recip_matrix(3,2,ibox)*rz  
                   sz = recip_matrix(1,3,ibox)*rx + &
                        recip_matrix(2,3,ibox)*ry + &
                        recip_matrix(3,3,ibox)*rz 

                   sx = sx*0.5_dp*invPi 
                   sy = sy*0.5_dp*invPi
                   sz = sz*0.5_dp*invPi 

                   ! apply boundary conditions
                   sx = sx - floor(sx+0.5_dp,kind=dp)
                   sy = sy - floor(sy+0.5_dp,kind=dp)
                   sz = sz - floor(sz+0.5_dp,kind=dp)

                   ! scale back up
                   rx = hmatrix(1,1,ibox)*sx + &
                        hmatrix(1,2,ibox)*sy + &
                        hmatrix(1,3,ibox)*sz
                                   
                   ry = hmatrix(2,1,ibox)*sx + &
                        hmatrix(2,2,ibox)*sy + &
                        hmatrix(2,3,ibox)*sz
                                   
                   rz = hmatrix(3,1,ibox)*sx + &
                        hmatrix(3,2,ibox)*sy + &
                        hmatrix(3,3,ibox)*sz

                   if ( (rx*rx+ry*ry+rz*rz < sigma_sq).and.(ichain/=jchain) ) then
                      ! ibead on ichain overlaps with jbead on jchain
                      noverlap = noverlap + 1
		      !write(*,*)'overlap between chain',ichain,'and',jchain,'separation',rx*rx+ry*ry+rz*rz
                      !loverlap(1,noverlap) = jbead   ! Store bead number
                      !loverlap(2,noverlap) = jchain  ! ..and chain number
                   end if

                   tmpint  = linked_list(1,jbead,jchain,ibox)
                   jchain  = linked_list(2,jbead,jchain,ibox)
                   jbead   = tmpint              

                end do

             end do

          end do

       else

          if ( ichain > 1 ) then
             do ibead = 1,nbeads
                rbead(:) = Rchain(:,ibead,ichain,ibox)
                do jchain = 1,ichain-1
                   do j = 1,nbeads
                      ! The compiler needs to inline this, use -ipo on Intel
                      rsep(:) = box_minimum_image( Rchain(:,j,jchain,ibox),rbead(:), ibox )
                      if (dot_product(rsep,rsep) < sigma_sq) then
                         ! ibead on ichain overlaps with jbead on jchain
                         noverlap = noverlap + 1
                         !loverlap(1,noverlap) = jbead   ! Store bead number
                         !loverlap(2,noverlap) = jchain  ! ..and chain number
                      end if
                   end do
                end do
             end do
          end if

          if ( ichain < nchains ) then
             do ibead = 1,nbeads
                rbead(:) = Rchain(:,ibead,ichain,ibox)
                do jchain = ichain + 1,nchains
                   do j = 1,nbeads
                      ! The compiler needs to inline this, use -ipo on Intel
                      rsep(:) = box_minimum_image( Rchain(:,j,jchain,ibox),rbead(:),ibox )
                      if (dot_product(rsep,rsep) < sigma_sq) then
                         ! ibead on ichain overlaps with jbead on jchain
                         noverlap = noverlap + 1
                         !loverlap(1,noverlap) = jbead   ! Store bead number
                         !loverlap(2,noverlap) = jchain  ! ..and chain number
                      end if
                   end do
                end do
             end do
          end if

       end if

    end if

!!$    if (noverlap> mxoverlap) stop 'Error mxoverlap exceeded in alkane_get_external_overlaps'

    ! Debug - plz comment out when happy
  !  write(0,*)
   ! write(0,'("Found ",I5," overlaps involving beads on chain ",I5," in box ",I5)')noverlap,ichain,ibox
    !write(0,*)
!!$    do iovlp = 1,noverlap
!!$       write(0,'("Bead ",I5," overlaps with bead ",I5," on chain ",I5)')ibead,loverlap(1,iovlp),loverlap(2,iovlp)
!!$    end do

    return

  end subroutine alkane_get_external_overlaps

  subroutine alkane_get_dr_max(dum_dr) bind(c)
    !-------------------------------------------------------------------------!
    ! Gets module level internal variable controlling the maximum molecule    !
    ! displacement during a translation move.                                 !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp),intent(out) :: dum_dr

    dum_dr = mc_dr_max

    return

  end subroutine alkane_get_dr_max

  subroutine alkane_set_dr_max(dum_dr) bind(c)
    !-------------------------------------------------------------------------!
    ! Sets module level internal variable controlling the maximum molecule    !
    ! displacement during a translation move.                                 !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp),intent(in) :: dum_dr

    mc_dr_max = dum_dr

    return

  end subroutine alkane_set_dr_max

  subroutine alkane_get_dt_max(dum_dt) bind(c)
    !-------------------------------------------------------------------------!
    ! Gets module level internal variable controlling the maximum molecule    !
    ! rotation during a chain rotation move.                                  !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp),intent(out) :: dum_dt

    dum_dt = mc_dt_max

    return

  end subroutine alkane_get_dt_max

  subroutine alkane_set_dt_max(dum_dt) bind(c)
    !-------------------------------------------------------------------------!
    ! Sets module level internal variable controlling the maximum molecule    !
    ! rotation during a chain rotation move.                                  !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp),intent(in) :: dum_dt

    mc_dt_max = dum_dt

    return

  end subroutine alkane_set_dt_max

  subroutine alkane_get_dv_max(dum_dv) bind(c)
    !-------------------------------------------------------------------------!
    ! Gets module level internal variable controlling the maximum volume      !
    ! change or cell vector displacement during a box resize/reshape move.    !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp),intent(out) :: dum_dv

    dum_dv = mc_dv_max

    return

  end subroutine alkane_get_dv_max

  subroutine alkane_set_dv_max(dum_dv) bind(c)
    !-------------------------------------------------------------------------!
    ! Sets module level internal variable controlling the maximum volume      !
    ! change or cell vector displacement during a box resize/reshape move.    !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp),intent(in) :: dum_dv

    mc_dv_max = dum_dv

    return

  end subroutine alkane_set_dv_max

  subroutine alkane_get_dh_max(dum_dh) bind(c)
    !-------------------------------------------------------------------------!
    ! Gets module level internal variable controlling the maximum change in   !
    ! dihedral angle during a torsion move.                                   !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp),intent(out) :: dum_dh

    dum_dh = mc_dh_max

    return

  end subroutine alkane_get_dh_max

  subroutine alkane_set_dh_max(dum_dh) bind(c)
    !-------------------------------------------------------------------------!
    ! Sets module level internal variable controlling the maximum change in   !
    ! dihedral angle during a torsion move.                                   !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp),intent(in) :: dum_dh

    mc_dh_max = dum_dh

    return

  end subroutine alkane_set_dh_max

   subroutine alkane_get_ktrial(dum_kt) bind(c)
    !-------------------------------------------------------------------------!
    ! Gets module level internal variable controlling the number of segment   !
    ! re-growths during a configurational bias MC move.                       !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(out) :: dum_kt

    dum_kt = ktrial

    return

  end subroutine alkane_get_ktrial

  subroutine alkane_set_ktrial(dum_kt) bind(c)
    !-------------------------------------------------------------------------!
    ! Sets module level internal variable controlling the number of segment   !
    ! re-growths during a configurational bias MC move.                       !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(in) :: dum_kt

    ktrial = dum_kt

    return

  end subroutine alkane_set_ktrial

   subroutine alkane_get_max_regrow(dum_max_regrow) bind(c)
    !-------------------------------------------------------------------------!
    ! Gets module level internal variable controlling the number of segment   !
    ! beads to regrow during a CBMC move                                      !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(out) :: dum_max_regrow

    dum_max_regrow = max_regrow

    return

  end subroutine alkane_get_max_regrow

  subroutine alkane_set_max_regrow(dum_max_regrow) bind(c)
    !-------------------------------------------------------------------------!
    ! Sets module level internal variable controlling the number of segment   !
    ! beads to regrow during a CBMC move                                      !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(in) :: dum_max_regrow

    max_regrow = dum_max_regrow

    return

  end subroutine alkane_set_max_regrow
 
  subroutine alkane_get_nchains(dumchains) bind(c)
    !-------------------------------------------------------------------------!
    ! Queries the number of chains per box in use by this module.             !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(out) :: dumchains

    dumchains = nchains

    return
    
  end subroutine alkane_get_nchains

  subroutine alkane_get_nbeads(dumbeads) bind(c)
    !-------------------------------------------------------------------------!
    ! Queries the number of beads per chain in use by this module.            !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(out) :: dumbeads

    dumbeads = nbeads

    return
    
  end subroutine alkane_get_nbeads

  subroutine alkane_set_chain(ichain,ibox,r) bind(c)
    !-------------------------------------------------------------------------!
    ! Overwrites the coordinates of a single chain in box ibox.               !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(in) :: ichain,ibox
    real(kind=dp),dimension(1:3,1:nbeads),intent(in) :: r

    Rchain(:,:,ichain,ibox) = r(:,:)

    return

  end subroutine alkane_set_chain

  subroutine alkane_get_chain(ichain,ibox,r) bind(c)
    !-------------------------------------------------------------------------!
    ! Returns the coordinates of a single chain in box ibox.                  !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(in) :: ichain,ibox
    real(kind=dp),dimension(1:3,1:nbeads),intent(out) :: r

    r(:,:) = Rchain(:,:,ichain,ibox) 

    return

  end subroutine alkane_get_chain

  subroutine alkane_change_box(ibox,delta_H) bind(c)
    !-------------------------------------------------------------------------!
    ! Implements a change in the matrix of cell vectors for box ibox, by      !
    ! the matrix delta_H. The box is changed, the chain first bead positions  !
    ! are scaled accordinly, and the matrix of reciprocal lattice vectors is  !
    ! updated for that box. For use in applying moves to an inactive box      !
    ! in lattice-switching calculations.                                      !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    use constants,only : invPi
    use box,      only : hmatrix,recip_matrix,box_construct_link_cells, &
                         box_update_recipmatrix,box_compute_volume
    use random,   only : random_uniform_random
    implicit none

    integer(kind=it),intent(in) :: ibox
    real(kind=dp),dimension(3,3),intent(in)     :: delta_H

!    real(kind=dp),dimension(3)   :: oldcom,comchain,tmpcom
    real(kind=dp),dimension(3)   :: first,frac_first,first_chain,delta_first
    real(kind=dp),dimension(3,3) :: old_hmatrix,new_hmatrix
    real(kind=dp) :: old_volume

    integer(kind=it) :: ichain,ibead


    ! Change hmatrix by delta_H
    old_hmatrix = hmatrix(:,:,ibox)
    old_volume  = box_compute_volume(ibox)

    new_hmatrix(:,:)  = old_hmatrix(:,:) + delta_H(:,:)
    hmatrix(:,:,ibox) = new_hmatrix(:,:)

    do ichain = 1,nchains

       ibead = 1
       first(:) = Rchain(:,ibead,ichain,ibox)

       ! Compute fractional first bead position using the current recip_matrix          
       frac_first(1) = recip_matrix(1,1,ibox)*first(1) + &
                       recip_matrix(2,1,ibox)*first(2) + &
                       recip_matrix(3,1,ibox)*first(3)
       frac_first(2) = recip_matrix(1,2,ibox)*first(1) + &
                       recip_matrix(2,2,ibox)*first(2) + &
                       recip_matrix(3,2,ibox)*first(3)  
       frac_first(3) = recip_matrix(1,3,ibox)*first(1) + &
                       recip_matrix(2,3,ibox)*first(2) + &
                       recip_matrix(3,3,ibox)*first(3) 

       frac_first = frac_first*0.5_dp*invPi 

       ! Scale to the new cell
       first_chain(1) = hmatrix(1,1,ibox)*frac_first(1) + &
                        hmatrix(1,2,ibox)*frac_first(2) + &
                        hmatrix(1,3,ibox)*frac_first(3)
                                
       first_chain(2) = hmatrix(2,1,ibox)*frac_first(1) + &
                        hmatrix(2,2,ibox)*frac_first(2) + &
                        hmatrix(2,3,ibox)*frac_first(3)
                                
       first_chain(3) = hmatrix(3,1,ibox)*frac_first(1) + &
                        hmatrix(3,2,ibox)*frac_first(2) + &
                        hmatrix(3,3,ibox)*frac_first(3)

       delta_first(:) = first_chain(:) - first(:)

       do ibead = 1,nbeads
          Rchain(:,ibead,ichain,ibox) = Rchain(:,ibead,ichain,ibox ) + delta_first(:)
       end do

    end do

    ! Book keeping
    call box_update_recipmatrix(ibox)
    call box_construct_link_cells(ibox)
    call alkane_construct_linked_lists(ibox)

    return

  end subroutine alkane_change_box



end module alkane
