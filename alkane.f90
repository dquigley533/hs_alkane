! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                            A  L  K  A  N  E                                 !
!=============================================================================!
!                                                                             !
! $Id: alkane.f90,v 1.2 2011/03/11 13:47:19 phseal Exp $
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
! Revision 1.2  2011/03/11 13:47:19  phseal
! Moved from single to double linked-lists
!
! Revision 1.1.1.1  2011/02/02 11:48:36  phseal
! Initial import from prototype code.
!
!
!=============================================================================!

module alkane

  use constants, only : dp,ep
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
  integer,save           :: nchains = 35        ! Number of chains
  integer,save           :: nbeads  = 4         ! length of alkane
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
  integer,save           :: model_type    = 4    ! default model type
  integer,save           :: torsion_type  = 1    ! default torsion type
  integer,save           :: nexclude             ! exclusion length
  

  ! CBMC parameters
  integer,save :: ktrial = 5                     ! Number of trial segments
  integer,save :: max_regrow = 3                 ! Maximum segments to regrow

  ! Other MC move parameters
  real(kind=dp),save :: mc_dr_max = 0.0421_dp    ! Maximum translation move
  real(kind=dp),save :: mc_dt_max = 0.3181_dp    ! Maximum rotation angle
  real(kind=dp),save :: mc_dv_max = 0.2958_dp    ! Maximum volume change
  real(kind=dp),save :: mc_dh_max = 0.0159_dp    ! Maximum dihedral change

  logical,save       :: rigid     = .false.      ! Chains are inflexible

  ! Chain coordinates and flag for initial configuration
  real(kind=dp),allocatable,dimension(:,:,:),save :: Rchain
  logical,allocatable,dimension(:),save           :: chain_created

  ! Intermolecular neighbour list (if applicable)
  integer,allocatable,dimension(:),save :: list,startinlist,endinlist

  ! Linked list arrays (if applicable) See F&S for details
  integer,allocatable,dimension(:,:),save   :: head_of_cell
  integer,allocatable,dimension(:,:,:),save :: linked_list

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!


contains

  subroutine alkane_init()
    !-------------------------------------------------------------------------!
    ! Initialises the alkane module                                           !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    implicit none

    integer,dimension(3) :: ierr = 0

    ! Allocate chain position array
    allocate(Rchain(1:3,1:nbeads,1:nchains),stat=ierr(1))
    if (any(ierr/=0)) stop 'Error allocating memory in alkane module'

    if ( model_type == 2 ) then
       nexclude = 3  ! 1-4 hard sphere interactions
    else
       nexclude = 4  ! no 1-4 hard sphere interactions
    end if

    allocate(chain_created(1:nchains),stat=ierr(1))
    if (any(ierr/=0)) stop 'Error allocating memory in alkane module'

    ! Chains have not yet been created
    chain_created = .false.

    ! Neighbour list
    allocate(list(1:nbeads*nchains*100),stat=ierr(1))
    allocate(startinlist(1:nbeads*nchains),stat=ierr(2))
    allocate(endinlist(1:nbeads*nchains),stat=ierr(3))
    if (any(ierr/=0)) stop 'Error allocating neighbour list arrays'
  
    return

  end subroutine alkane_init

  subroutine alkane_destroy()
    !-------------------------------------------------------------------------!
    ! Performs a clean shutdown of the alkane module                          !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use box, only : use_link_cells
    implicit none
    
    integer,dimension(2) :: ierr

    deallocate(Rchain,chain_created,list,startinlist,endinlist,stat=ierr(1))
    if (use_link_cells) deallocate(head_of_cell,linked_list,stat=ierr(2))
    if (any(ierr/=0)) stop 'Error releasing memory in alkane module'

  end subroutine alkane_destroy

  subroutine alkane_translate_chain(ichain,new_boltz)
    !-------------------------------------------------------------------------!
    ! Implements an MC trial move in which the specified chain is translated  !
    ! by a random vector. The new Boltzmann factor after the trial move is    !
    ! returned as new_boltz. The old Boltzmann factor will always be one.     !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use random, only : random_uniform_random
    implicit none
    integer,intent(in)         :: ichain
    real(kind=dp),intent(out)  :: new_boltz
    real(kind=dp),dimension(3) :: dr
    real(kind=dp)              :: tboltz
    integer :: ibead

    ! generate random move
    dr(1) = 2.0_dp*random_uniform_random() - 1.0_dp
    dr(2) = 2.0_dp*random_uniform_random() - 1.0_dp
    dr(3) = 2.0_dp*random_uniform_random() - 1.0_dp

    dr = dr*mc_dr_max
    
    ! translate the chain
    do ibead = 1,nbeads
       Rchain(:,ibead,ichain) = Rchain(:,ibead,ichain) + dr(:)
    end do

    tboltz    = alkane_chain_inter_boltz(ichain)
    new_boltz = tboltz 

    return


  end subroutine alkane_translate_chain

  subroutine alkane_rotate_chain(ichain,new_boltz)
    !-------------------------------------------------------------------------!
    ! Implements an MC trial move in which the specified chain is rotated     !
    ! about its centre of mass. The new Boltzmann factor after the trial move !
    ! is returned as new_boltz. The old Boltzmann factor will always be one.  !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use random, only     : random_uniform_random
    use quaternion, only : quat_axis_angle_to_quat,quat_conjugate_q_with_v 
    implicit none
    integer,intent(in)         :: ichain
    real(kind=dp),intent(out)  :: new_boltz
    real(kind=dp),dimension(3) :: axis,rcom
    real(kind=dp),dimension(4) :: quat
    real(kind=dp)              :: theta

    integer :: ibead

    ! generate random rotation axis
    axis(1) = 2.0_dp*random_uniform_random() - 1.0_dp
    axis(2) = 2.0_dp*random_uniform_random() - 1.0_dp
    axis(3) = 2.0_dp*random_uniform_random() - 1.0_dp
    axis(:) = axis(:)/sqrt(dot_product(axis,axis))

    ! generate random rotation angle
    theta = (2.0_dp*random_uniform_random() - 1.0_dp) * mc_dt_max

    call quat_axis_angle_to_quat(axis,theta,quat)

    ! chain center of mass
    rcom(:) = 0.0_dp
    do ibead = 1,nbeads
       rcom(:) = rcom(:) + Rchain(:,ibead,ichain)
    end do
    rcom(:) = rcom(:)/real(nbeads,kind=dp)

    ! rotate the chain
    do ibead = 1,nbeads
       Rchain(:,ibead,ichain) = Rchain(:,ibead,ichain) - rcom(:)
       Rchain(:,ibead,ichain) = quat_conjugate_q_with_v(quat,Rchain(:,ibead,ichain))
       Rchain(:,ibead,ichain) = Rchain(:,ibead,ichain) + rcom(:)
    end do

    new_boltz = alkane_chain_inter_boltz(ichain)

    return

  end subroutine alkane_rotate_chain

  subroutine alkane_box_resize(pressure,acc_prob,reset)
    !-------------------------------------------------------------------------!
    ! Implements an MC trial move in which the size (and possibly shape) of   !
    ! the simulation box is altered by a random amount. The ratio  of new/old !
    ! Boltmann factors (i.e. the acceptance probability) is returned as       !
    ! new_boltz. The centre-of-mass position of each chain in fractional      !
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
    logical,optional,intent(in)  :: reset
    real(kind=dp),intent(out)    :: acc_prob
    real(kind=dp),dimension(3)   :: oldcom,comchain,tmpcom
    real(kind=dp),dimension(3,3),save :: old_hmatrix,new_hmatrix,delta_hmatrix
    real(kind=dp) :: old_volume,new_volume,delta_vol,x

    integer :: ichain,ibead,jdim,idim

    if ( present(reset).and.reset ) then

       do ichain = 1,nchains

          comchain(:) = 0.0_dp
          do ibead = 1,nbeads
             comchain(:) = comchain(:) + Rchain(:,ibead,ichain)
          end do
          comchain(:) = comchain(:)/real(nbeads,kind=dp)
          oldcom(:)   = comchain(:)

          ! Compute fractional com position using the current recip_matrix          
          tmpcom(1) = recip_matrix(1,1)*oldcom(1) + &
                      recip_matrix(2,1)*oldcom(2) + &
                      recip_matrix(3,1)*oldcom(3)
          tmpcom(2) = recip_matrix(1,2)*oldcom(1) + &
                      recip_matrix(2,2)*oldcom(2) + &
                      recip_matrix(3,2)*oldcom(3)  
          tmpcom(3) = recip_matrix(1,3)*oldcom(1) + &
                      recip_matrix(2,3)*oldcom(2) + &
                      recip_matrix(3,3)*oldcom(3) 

          tmpcom = tmpcom*0.5_dp*invPi 


          ! Scale to the previous cell
          comchain(1) = old_hmatrix(1,1)*tmpcom(1) + &
                        old_hmatrix(1,2)*tmpcom(2) + &
                        old_hmatrix(1,3)*tmpcom(3)

          comchain(2) = old_hmatrix(2,1)*tmpcom(1) + &
                        old_hmatrix(2,2)*tmpcom(2) + &
                        old_hmatrix(2,3)*tmpcom(3)

          comchain(3) = old_hmatrix(3,1)*tmpcom(1) + &
                        old_hmatrix(3,2)*tmpcom(2) + &
                        old_hmatrix(3,3)*tmpcom(3)

          tmpcom(:) = comchain(:) - oldcom(:)

          !write(*,'(3F15.6)')tmpcom

          do ibead = 1,nbeads
             Rchain(:,ibead,ichain) = Rchain(:,ibead,ichain) + tmpcom(:)
          end do

       end do

       hmatrix = old_hmatrix(:,:)
       call box_update_recipmatrix()


    else


       old_hmatrix  = hmatrix
       old_volume   = box_compute_volume()

       if (isotropic) then

          ! Uniform scaling of box dimensions
          new_volume  = old_volume+(2.0_dp*random_uniform_random()-1.0_dp)*mc_dv_max
          delta_vol   = new_volume - old_volume
          new_hmatrix = old_hmatrix*(new_volume/old_volume)**(1.0_dp/3.0_dp)
          hmatrix     = new_hmatrix

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

          new_hmatrix  = old_hmatrix + delta_hmatrix
          hmatrix(:,:) = new_hmatrix(:,:)
          new_volume   = box_compute_volume()
          delta_vol    = new_volume - old_volume
          

       end if

       do ichain = 1,nchains

          comchain(:) = 0.0_dp
          do ibead = 1,nbeads
             comchain(:) = comchain(:) + Rchain(:,ibead,ichain)
          end do
          comchain(:) = comchain(:)/real(nbeads,kind=dp)
          oldcom(:)   = comchain(:)

          ! Compute fractional com position using the current recip_matrix          
          tmpcom(1) = recip_matrix(1,1)*oldcom(1) + &
                      recip_matrix(2,1)*oldcom(2) + &
                      recip_matrix(3,1)*oldcom(3)
          tmpcom(2) = recip_matrix(1,2)*oldcom(1) + &
                      recip_matrix(2,2)*oldcom(2) + &
                      recip_matrix(3,2)*oldcom(3)  
          tmpcom(3) = recip_matrix(1,3)*oldcom(1) + &
                      recip_matrix(2,3)*oldcom(2) + &
                      recip_matrix(3,3)*oldcom(3) 

          tmpcom = tmpcom*0.5_dp*invPi 

          ! Scale to the new cell
          comchain(1) = hmatrix(1,1)*tmpcom(1) + &
                        hmatrix(1,2)*tmpcom(2) + &
                        hmatrix(1,3)*tmpcom(3)

          comchain(2) = hmatrix(2,1)*tmpcom(1) + &
                        hmatrix(2,2)*tmpcom(2) + &
                        hmatrix(2,3)*tmpcom(3)

          comchain(3) = hmatrix(3,1)*tmpcom(1) + &
                        hmatrix(3,2)*tmpcom(2) + &
                        hmatrix(3,3)*tmpcom(3)

          tmpcom(:) = comchain(:) - oldcom(:)

          !write(*,'(3F15.6)')tmpcom

          do ibead = 1,nbeads
             Rchain(:,ibead,ichain) = Rchain(:,ibead,ichain) + tmpcom(:)
          end do

       end do

       !stop
       call box_update_recipmatrix()

    end if

    call box_construct_link_cells(1.001_dp)

    call alkane_construct_linked_lists()

    if ( present(reset).and.reset ) then
       ! Restoring the old configuration, acc_prob irrelevant
       acc_prob = 1.0_dp
       return
    else
       ! This is a trial move and we need to check for overlaps
       ! if the box has shrunk
       acc_prob = exp(-pressure*delta_vol + real(nchains,kind=dp)*log(new_volume/old_volume))

       !if ( any(boxscale<1.0_dp) ) then

       do ichain = 1,nchains
          if ( alkane_chain_inter_boltz(ichain) < tiny(1.0_dp) ) then
             acc_prob = 0.0_dp
             return
          end if
       end do

       !end if

    end if

    return

  end subroutine alkane_box_resize

  subroutine alkane_box_scale(scaleA,scaleB,scaleC)
    !-------------------------------------------------------------------------!
    ! Implements an MC trial move in which the size (and possibly shape) of   !
    ! the simulation box is altered by a random amount. The ratio  of new/old !
    ! Boltmann factors (i.e. the acceptance probability) is returned as       !
    ! new_boltz. The centre-of-mass position of each chain in fractional      !
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

    real(kind=dp),intent(in)     :: scaleA,scaleB,scaleC
    real(kind=dp),dimension(3)   :: oldcom,comchain,tmpcom
    real(kind=dp),dimension(3,3) :: old_hmatrix,new_hmatrix
    real(kind=dp) :: old_volume,new_volume,delta_vol

    integer :: ichain,ibead,jdim,idim


    ! Uniform scaling of box dimensions
    old_hmatrix = hmatrix
    old_volume  = box_compute_volume()

    new_hmatrix(:,1) = old_hmatrix(:,1)*scaleA
    new_hmatrix(:,2) = old_hmatrix(:,2)*scaleA
    new_hmatrix(:,3) = old_hmatrix(:,3)*scaleA
    
    hmatrix(:,:)  = new_hmatrix(:,:)

    do ichain = 1,nchains

       comchain(:) = 0.0_dp
       do ibead = 1,nbeads
          comchain(:) = comchain(:) + Rchain(:,ibead,ichain)
       end do
       comchain(:) = comchain(:)/real(nbeads,kind=dp)
       oldcom(:)   = comchain(:)

       ! Compute fractional com position using the current recip_matrix          
       tmpcom(1) = recip_matrix(1,1)*oldcom(1) + &
                   recip_matrix(2,1)*oldcom(2) + &
                   recip_matrix(3,1)*oldcom(3)
       tmpcom(2) = recip_matrix(1,2)*oldcom(1) + &
                   recip_matrix(2,2)*oldcom(2) + &
                   recip_matrix(3,2)*oldcom(3)  
       tmpcom(3) = recip_matrix(1,3)*oldcom(1) + &
                   recip_matrix(2,3)*oldcom(2) + &
                   recip_matrix(3,3)*oldcom(3) 

       tmpcom = tmpcom*0.5_dp*invPi 

       ! Scale to the new cell
       comchain(1) = hmatrix(1,1)*tmpcom(1) + &
                     hmatrix(1,2)*tmpcom(2) + &
                     hmatrix(1,3)*tmpcom(3)

       comchain(2) = hmatrix(2,1)*tmpcom(1) + &
                     hmatrix(2,2)*tmpcom(2) + &
                     hmatrix(2,3)*tmpcom(3)

       comchain(3) = hmatrix(3,1)*tmpcom(1) + &
                     hmatrix(3,2)*tmpcom(2) + &
                     hmatrix(3,3)*tmpcom(3)

       tmpcom(:) = comchain(:) - oldcom(:)

       do ibead = 1,nbeads
          Rchain(:,ibead,ichain) = Rchain(:,ibead,ichain) + tmpcom(:)
       end do

    end do

    call box_update_recipmatrix()

    call box_construct_link_cells(1.001_dp)

    call alkane_construct_linked_lists()

    return

  end subroutine alkane_box_scale


  subroutine alkane_bond_rotate(ichain,new_boltz)
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
    integer,intent(in)        :: ichain
    real(kind=dp),intent(out) :: new_boltz

    real(kind=dp),dimension(3) :: axis,r12,r23,r34
    real(kind=dp),dimension(4) :: quat

    real(kind=dp) :: angle,xi

    integer :: ibead,ia

    if (nbeads<4) stop 'Called alkane_bond_rotate with nbeads < 4'

    ! select a bead at random from 1 to nbeads-3
    ia = int(real(nbeads-3,kind=dp)*random_uniform_random()) + 1
    ia = min(ia,nbeads-3)

    ! we rotate about the vector r23
    axis(:) = Rchain(:,ia+2,ichain) - Rchain(:,ia+1,ichain)
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
       Rchain(:,ibead,ichain) = Rchain(:,ibead,ichain) - Rchain(:,ia+2,ichain)
       Rchain(:,ibead,ichain) = quat_conjugate_q_with_v(quat,Rchain(:,ibead,ichain))
       Rchain(:,ibead,ichain) = Rchain(:,ibead,ichain) + Rchain(:,ia+2,ichain)
    end do

    ! compute Boltzmann factor for acceptance
    r12(:) = Rchain(:,ia+1,ichain) - Rchain(:,ia,ichain)
    r23(:) = axis(:)*L
    r34(:) = Rchain(:,ia+3,ichain) - Rchain(:,ia+2,ichain)
    new_boltz = alkane_dihedral_boltz(r12,r23,r34)

    ! return now if already zero
    if (new_boltz<tiny(1.0_dp)) return

    ! check non-bonded interactions (including 1-4 for model II)
    do ibead = ia+3,nbeads
       new_boltz = new_boltz*alkane_nonbonded_boltz(ibead,ichain,Rchain(:,ibead,ichain))
       if (new_boltz<tiny(1.0_dp)) return
    end do

    return

  end subroutine alkane_bond_rotate

  function alkane_dihedral_boltz(b1,b2,b3)
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

  subroutine alkane_grow_chain(ichain,rb_factor,new_conf)
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

    integer,intent(in)        :: ichain        ! chain number to grow/regrow
    real(kind=ep),intent(out) :: rb_factor     ! Rosenbluth factor for chain 
    logical,intent(in)        :: new_conf      ! old or new configuration?

    ! Dihedral / angle calculation
    real(kind=dp),dimension(3) :: r12,r23,r34,tmpvect,axis
    real(kind=dp),dimension(4) :: quat
    real(kind=dp)              :: dih,theta

    ! Trial segments and weights
    real(kind=dp),allocatable,dimension(:,:) :: rtrial
    real(kind=dp),allocatable,dimension(:)   :: wtrial
    real(kind=dp),allocatable,dimension(:)   :: wset

    ! Loop counters / error flags
    integer :: ierr,n,j,ib,jl,ifail

    ! Rosenbluth variables
    real(kind=dp) :: wsum

    ! Bead from which to start growth
    integer,save :: first_bead

    ! Allocate memory for trial segments
    allocate(rtrial(1:3,1:ktrial),stat=ierr)
    if ( ierr/=0 ) stop 'Error allocating Rtrial'

    allocate(wtrial(1:ktrial),stat=ierr)
    if ( ierr/=0 ) stop 'Error allocating wtrial'

    ! Set bead from which to (re)grow
    if ( .not.chain_created(ichain) ) then
       first_bead = 1  ! grow whole chain from scratch
    elseif (.not.new_conf) then
       first_bead = int(random_uniform_random()*real(nbeads-max_regrow,kind=dp)) + 1 + max_regrow
    end if

    allocate(wset(first_bead:nbeads),stat=ierr)
    if (ierr/=0) stop 'Error allocating wset'

    ! Set loop counters and accumulators
    if (new_conf) then
       ! All trial segments/weights are new
       jl   = 1
       wset = 0.0_dp
    else
       ! The first trial segment at each bead is the old config
       jl = 2
       do ib = first_bead,nbeads
          wset(ib) = alkane_nonbonded_boltz(ib,ichain,Rchain(:,ib,ichain))
       end do
    end if

    ! Loop over beads
    rb_factor = 1_ep       
    do ib = first_bead,nbeads

       !======================================================!
       ! First bead                                           !
       !======================================================!
       if ( ib==1 ) then

          if (new_conf) then
             Rchain(1,1,ichain) = random_uniform_random()
             Rchain(2,1,ichain) = random_uniform_random()
             Rchain(3,1,ichain) = random_uniform_random()
             Rchain(:,1,ichain) = matmul(hmatrix,Rchain(:,1,ichain))
          end if

          rb_factor = real(alkane_nonbonded_boltz(1,ichain,Rchain(:,1,ichain)),kind=ep)

          !write(0,'(I5,3F15.6)')ib,Rchain(:,ib,ichain)

       !======================================================!
       ! Second bead                                          !
       !======================================================!
       elseif ( ib==2 ) then

          wsum = wset(ib)
          do j = jl,ktrial
             ! All trial segments proposed with equal probability
             rtrial(:,j) = Rchain(:,1,ichain) + L*random_unit_vector()
             wtrial(j)   = alkane_nonbonded_boltz(2,ichain,rtrial(:,j))
             !print*,wtrial(j)
             wsum = wsum + wtrial(j) 
          end do
          !stop

          rb_factor = rb_factor*real(wsum,kind=ep)
          call select_next_segment(n,ifail)         
          if (ifail/=0) return ! Rosenbluth factor is zero
          if (new_conf) Rchain(:,2,ichain) = rtrial(:,n)

          !write(0,'(I5,3F15.6)')ib,Rchain(:,ib,ichain)

       !======================================================!
       ! Third bead                                           !
       !======================================================!
       elseif ( ib==3 ) then

          r12 = Rchain(:,2,ichain) - Rchain(:,1,ichain)
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
             rtrial(:,j) = Rchain(:,1,ichain) + rtrial(:,j)

             wtrial(j)   = alkane_nonbonded_boltz(3,ichain,rtrial(:,j))
             wsum = wsum + wtrial(j) 

          end do

          rb_factor = rb_factor*real(wsum,kind=ep)
          call select_next_segment(n,ifail)         
          if (ifail/=0) return ! Rosenbluth factor is zero
          if (new_conf) Rchain(:,3,ichain) = rtrial(:,n)

       else

          !======================================================!
          ! Fourth and subsequent beads                          !
          !======================================================!

          r12  = Rchain(:,ib-2,ichain) - Rchain(:,ib-3,ichain)
          r23  = Rchain(:,ib-1,ichain) - Rchain(:,ib-2,ichain)
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

             rtrial(:,j) = Rchain(:,ib-1,ichain) + r34
             !call alkane_check_dihedral(r12,r23,r34)

             wtrial(j)   = alkane_nonbonded_boltz(ib,ichain,rtrial(:,j))
             wsum = wsum + wtrial(j) 

          end do

          rb_factor = rb_factor*real(wsum,kind=ep)
          call select_next_segment(n,ifail)         
          if (ifail/=0) return ! Rosenbluth factor is zero
          if (new_conf) Rchain(:,ib,ichain) = rtrial(:,n)

       end if

    end do

    ! This chain now exists if it didn't before
    chain_created(ichain) = .true.


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
      integer,intent(out) :: nout,ifail
      real(kind=dp) :: sump,zeta,rws

      sump = 0.0_dp
      zeta = random_uniform_random()
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

      if (nout==-1) ifail = 1


      return

    end subroutine select_next_segment

  end subroutine alkane_grow_chain

  function alkane_nonbonded_boltz(i,ichain,rbead)
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

    integer,intent(in) :: i,ichain
    real(kind=dp),dimension(3),intent(in) :: rbead
    real(kind=dp),dimension(3) :: rsep,sbead

    integer :: jchain,ix,iy,iz,icell,ni,jcell,jbead,tmpint
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
          rsep(:) = box_minimum_image( Rchain(:,jbead,ichain), rbead(:) )
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
          sbead(1) = recip_matrix(1,1)*rbead(1) + &
                     recip_matrix(2,1)*rbead(2) + &
                     recip_matrix(3,1)*rbead(3)
          sbead(2) = recip_matrix(1,2)*rbead(1) + &
                     recip_matrix(2,2)*rbead(2) + &
                     recip_matrix(3,2)*rbead(3)  
          sbead(3) = recip_matrix(1,3)*rbead(1) + &
                     recip_matrix(2,3)*rbead(2) + &
                     recip_matrix(3,3)*rbead(3) 

          sbead = sbead*0.5_dp*invPi 

          ! link cell containing rbead
          ix = floor(sbead(1)/lcellx)
          iy = floor(sbead(2)/lcelly)
          iz = floor(sbead(3)/lcellz)

          ix = modulo(ix,ncellx) + 1
          iy = modulo(iy,ncelly) + 1
          iz = modulo(iz,ncellz) + 1

          icell = (iz-1)*ncellx*ncelly + (iy-1)*ncellx + ix

          ! loop over link cells
          do ni = 1,27
             jcell   = lcneigh(ni,icell)
             jbead   = head_of_cell(1,jcell)
             jchain  = head_of_cell(2,jcell)
             do while ( jchain.ne.0 )

                !rsep(:) = box_minimum_image(Rchain(:,jbead,jchain),rbead(:))

                rx = rbead(1) - Rchain(1,jbead,jchain)
                ry = rbead(2) - Rchain(2,jbead,jchain)
                rz = rbead(3) - Rchain(3,jbead,jchain)

                sx = recip_matrix(1,1)*rx + &
                     recip_matrix(2,1)*ry + &
                     recip_matrix(3,1)*rz
                sy = recip_matrix(1,2)*rx + &
                     recip_matrix(2,2)*ry + &
                     recip_matrix(3,2)*rz  
                sz = recip_matrix(1,3)*rx + &
                     recip_matrix(2,3)*ry + &
                     recip_matrix(3,3)*rz 

                sx = sx*0.5_dp*invPi 
                sy = sy*0.5_dp*invPi
                sz = sz*0.5_dp*invPi 

                ! apply boundary conditions
                sx = sx - floor(sx+0.5_dp,kind=dp)
                sy = sy - floor(sy+0.5_dp,kind=dp)
                sz = sz - floor(sz+0.5_dp,kind=dp)

                ! scale back up
                rx = hmatrix(1,1)*sx + &
                     hmatrix(1,2)*sy + &
                     hmatrix(1,3)*sz

                ry = hmatrix(2,1)*sx + &
                     hmatrix(2,2)*sy + &
                     hmatrix(2,3)*sz

                rz = hmatrix(3,1)*sx + &
                     hmatrix(3,2)*sy + &
                     hmatrix(3,3)*sz

!!$                rx = rx - Lx*anint(rx*rLx)
!!$                ry = ry - Ly*anint(ry*rLy)
!!$                rz = rz - Lz*anint(rz*rLz)

                overlap = overlap.or.( ( rx*rx+ry*ry+rz*rz < sigma_sq ).and.(ichain/=jchain) )

                tmpint  = linked_list(1,jbead,jchain)
                jchain  = linked_list(2,jbead,jchain)
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
                   rsep(:) = box_minimum_image( Rchain(:,jbead,jchain),rbead(:) )
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
                   rsep(:) = box_minimum_image( Rchain(:,jbead,jchain),rbead(:) )
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

  function alkane_chain_inter_boltz(ichain)
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

    Integer,intent(in) :: ichain

    real(kind=dp),dimension(3) :: rbead
    real(kind=dp),dimension(3) :: rsep,sbead

    integer :: j,jchain,ibead,icell,ix,iy,iz,jbead,jcell,ni,tmpint
    real(kind=dp) :: sigma_sq
    real(kind=dp) :: rx,ry,rz
    real(kind=dp) :: sx,sy,sz

    real(kind=dp) :: alkane_chain_inter_boltz
    logical       :: overlap

    sigma_sq = sigma*sigma
    overlap  = .false.

    ! Run over beads in other chains....
    if (nchains>1) then

       if ( use_link_cells ) then

          do ibead = 1,nbeads

             rbead(:) = Rchain(:,ibead,ichain)

             ! compute fractional coordinates sbead from rbead
             sbead(1) = recip_matrix(1,1)*rbead(1) + &
                        recip_matrix(2,1)*rbead(2) + &
                        recip_matrix(3,1)*rbead(3)
             sbead(2) = recip_matrix(1,2)*rbead(1) + &
                        recip_matrix(2,2)*rbead(2) + &
                        recip_matrix(3,2)*rbead(3)  
             sbead(3) = recip_matrix(1,3)*rbead(1) + &
                        recip_matrix(2,3)*rbead(2) + &
                        recip_matrix(3,3)*rbead(3) 

             sbead = sbead*0.5_dp*invPi 

             ! link cell containing rbead
             ix = floor(sbead(1)/lcellx)
             iy = floor(sbead(2)/lcelly)
             iz = floor(sbead(3)/lcellz)

             ix = modulo(ix,ncellx) + 1
             iy = modulo(iy,ncelly) + 1
             iz = modulo(iz,ncellz) + 1

             icell = (iz-1)*ncellx*ncelly + (iy-1)*ncellx + ix

             ! loop over link cells
             do ni = 1,27
                jcell   = lcneigh(ni,icell)
                jbead   = head_of_cell(1,jcell)
                jchain  = head_of_cell(2,jcell)
                do while ( jchain.ne.0 )

                   !rsep(:) = box_minimum_image(Rchain(:,jbead,jchain),rbead(:))
                   rx = rbead(1) - Rchain(1,jbead,jchain)
                   ry = rbead(2) - Rchain(2,jbead,jchain)
                   rz = rbead(3) - Rchain(3,jbead,jchain)

                   sx = recip_matrix(1,1)*rx + &
                        recip_matrix(2,1)*ry + &
                        recip_matrix(3,1)*rz
                   sy = recip_matrix(1,2)*rx + &
                        recip_matrix(2,2)*ry + &
                        recip_matrix(3,2)*rz  
                   sz = recip_matrix(1,3)*rx + &
                        recip_matrix(2,3)*ry + &
                        recip_matrix(3,3)*rz 

                   sx = sx*0.5_dp*invPi 
                   sy = sy*0.5_dp*invPi
                   sz = sz*0.5_dp*invPi 

                   ! apply boundary conditions
                   sx = sx - floor(sx+0.5_dp,kind=dp)
                   sy = sy - floor(sy+0.5_dp,kind=dp)
                   sz = sz - floor(sz+0.5_dp,kind=dp)

                   ! scale back up
                   rx = hmatrix(1,1)*sx + &
                        hmatrix(1,2)*sy + &
                        hmatrix(1,3)*sz

                   ry = hmatrix(2,1)*sx + &
                        hmatrix(2,2)*sy + &
                        hmatrix(2,3)*sz

                   rz = hmatrix(3,1)*sx + &
                        hmatrix(3,2)*sy + &
                        hmatrix(3,3)*sz

                   overlap = overlap.or.( (rx*rx+ry*ry+rz*rz < sigma_sq).and.(ichain/=jchain) )

                   tmpint  = linked_list(1,jbead,jchain)
                   jchain  = linked_list(2,jbead,jchain)
                   jbead   = tmpint              

                end do
                if ( overlap ) then
                   alkane_chain_inter_boltz = 0.0_dp
                   return
                end if
             end do

          end do

       else

          if ( ichain > 1 ) then
             do ibead = 1,nbeads
                rbead(:) = Rchain(:,ibead,ichain)
                do jchain = 1,ichain-1
                   do j = 1,nbeads
                      ! The compiler needs to inline this, use -ipo on Intel
                      rsep(:) = box_minimum_image( Rchain(:,j,jchain),rbead(:) )
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
                rbead(:) = Rchain(:,ibead,ichain)
                do jchain = ichain + 1,nchains
                   do j = 1,nbeads
                      ! The compiler needs to inline this, use -ipo on Intel
                      rsep(:) = box_minimum_image( Rchain(:,j,jchain),rbead(:) )
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

  function alkane_random_dihedral()
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

    end select

    alkane_random_dihedral = phi*Pi/180.0_dp

    return

  end function alkane_random_dihedral
    
  subroutine alkane_check_dihedral(b1,b2,b3)
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
    real(kind=dp),dimension(3)            :: t1,t2
    real(kind=dp) :: arg1,arg2,angle

    t1 = cross_product(b2,b3)
    t2 = cross_product(b1,b2)

    arg1 = L*dot_product(b1,t1)
    arg2 =   dot_product(t2,t1)

    angle = atan2(arg1,arg2)*180_dp/Pi 

    ! convect to be consistent with the angle-origin define by the model
    angle = -sign(180.0_dp-abs(angle),angle)

    write(0,'("Measured dihedral angle as : ",F15.6)')angle 

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



  subroutine alkane_check_chain_overlap(overlap)
    !-------------------------------------------------------------------------!
    ! Sanity test for debugging. Checks if any two chains overlap. Does not   !
    ! check if a chain overlaps with itself.                                  !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    implicit none
    real(kind=dp) :: test,acc
    integer :: ichain
    logical,intent(out) :: overlap

    if (nchains < 2) stop 'Called alkane_check_chain_overlap with one chain'
    
    overlap = .false.

    acc = 1.0_dp
    do ichain = 1,nchains
       test = alkane_chain_inter_boltz(ichain)
       acc  = acc*test
       if (test < tiny(1.0_dp) ) then
          write(0,'("Chain ",I5", overlaps with another chain")')ichain
       end if
    end do

    if (acc<tiny(1.0_dp)) then
       !write(0,'("Stopping")')
       !stop
       overlap = .true.
    end if

    return

  end subroutine alkane_check_chain_overlap

  subroutine alkane_check_chain_geometry(ichain,violated)
    !-------------------------------------------------------------------------!
    ! Sanity test for debugging. Checks bond lengths within a chain           !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : Pi
    use box
    implicit none
    integer,intent(in)  :: ichain
    logical,intent(out) :: violated 
    real(kind=dp),dimension(3) :: rsep,r12,r23,r34
    integer :: ibead,jbead
    real(kind=dp) :: boltzf

    violated = .false. ! No problems found yet

    !==============================================!
    ! Check that all bond lengths are equal to L   !
    !==============================================!
    do ibead = 1,nbeads-1
       jbead = ibead+1
       rsep(:) = box_minimum_image( Rchain(:,jbead,ichain), Rchain(:,ibead,ichain) )
       if ( dot_product(rsep,rsep) - L*L > 0.00001_dp ) then
          violated = .true.
          write(0,'("Found a bond length of ",E15.6," between beads ",I5," and ",I5)') &
               sqrt(dot_product(rsep,rsep)),ibead,jbead
       end if
    end do

    if ( violated ) then
       write(0,'("Bond length violation for chain ",I5)')ichain
       write(0,'("Other chains may be affected")')
       return
    end if


    !=================================================!
    ! Check that bond angles are equal to 109.47 deg  !
    !=================================================!
    do ibead=1,nbeads-2

       r12(:) = box_minimum_image(Rchain(:,ibead,  ichain),Rchain(:,ibead+1,ichain))
       r23(:) = box_minimum_image(Rchain(:,ibead+2,ichain),Rchain(:,ibead+1,ichain))

       r12(:) = r12(:)/sqrt(dot_product(r12,r12))
       r23(:) = r23(:)/sqrt(dot_product(r23,r23))

       if ( abs(acos(dot_product(r12,r23)) - 109.47_dp*Pi/180.0_dp) > 0.0002_dp ) then
          violated = .true.
          write(0,'("Found a bond angle of ",F15.6," involving beads ",3I5)') &
               acos(dot_product(r12,r23))*180.0_dp/Pi,ibead,ibead+1,ibead+2
       end if

    end do

    if ( violated ) then
       write(0,'("Bond angle violation for chain ",I5)')ichain
       write(0,'("Other chains may be affected")')
       return
    end if


    !==============================================!
    ! Check that torsions are sensible             !
    !==============================================!
    do ibead = 1,nbeads-3

       ! compute Boltzmann factor for acceptance
       r12(:) = box_minimum_image(Rchain(:,ibead,  ichain),Rchain(:,ibead+1,ichain))
       r23(:) = box_minimum_image(Rchain(:,ibead+1,ichain),Rchain(:,ibead+2,ichain))
       r34(:) = box_minimum_image(Rchain(:,ibead+2,ichain),Rchain(:,ibead+3,ichain))
       boltzf = alkane_dihedral_boltz(r12,r23,r34)

       if (boltzf<tiny(1.0_dp)) then
          violated = .true.
          write(0,'("Found a bad dihedral angle involving beads ",I5,",",I5," and ",I5)') &
               ibead,ibead+1,ibead+2
          call alkane_check_dihedral(r12,r23,r34)
       end if

    end do

    if ( violated ) then
       write(0,'("Dihedral angle violation for chain ",I5)')ichain
       write(0,'("Other chains may be affected")')
       return
    end if

    !==============================================!
    ! Check that there are no overlaps             !
    !==============================================!
    do ibead = 1,nbeads-nexclude
       do jbead = ibead+nexclude,nbeads
          !write(0,'("Checking overlap between beads ",I5," and ",I5," on chain ",I5)')ibead,jbead,ichain
          r12(:) = box_minimum_image(Rchain(:,ibead, ichain),Rchain(:,jbead,ichain))
          !print*,sqrt(dot_product(r12,r12))
          if ( dot_product(r12,r12) < sigma*sigma ) then
             violated = .true.
             write(0,'("Found an overlap between beads ",I5," and ",I5," on chain ",I5)')ibead,jbead,ichain
          end if
       end do
    end do

    if ( violated ) then
       write(0,'("Intra-chain bead overlap violation for chain ",I5)')ichain
       write(0,'("Other chains may be affected")')
       return
    end if

    return

  end subroutine alkane_check_chain_geometry

  subroutine alkane_construct_neighbour_list()
    !-------------------------------------------------------------------------!
    ! Constructs a Verlet neighbour list. All intra-chain interactions are    !
    ! excluded from this list and must therefore be included in a seperate    !
    ! search.                                                                 !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use box
    implicit none

    integer :: myint,logint,ibead,jbead,ichain,jchain,k,ierr,m
    integer :: t1,t2,rate
    logical :: lrange
    real(kind=dp) :: nl_range_sq 
    real(kind=dp),dimension(3) :: rbead,rsep

    integer,allocatable,dimension(:) :: advance

    logical,save  :: firstpass = .true.
    integer,save  :: lcnv

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

          rbead(:) = Rchain(:,ibead,ichain)

          k = 1
          do jchain = 1,nchains
             do jbead = 1,nbeads
                ! Compiler need to in-line this next bit
                rsep(:) = box_minimum_image( Rchain(:,jbead,jchain),rbead(:) )
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
          startinlist((ichain-1)*nbeads+ibead) = m
          do jchain = 1,nchains
             do jbead = 1,nbeads
                list(m) = k
                m = m + advance(k)
                k = k + 1
             end do
          end do
          endinlist((ichain-1)*nbeads+ibead) = m - 1

       end do
    end do

    deallocate(advance,stat=ierr)
    if (ierr/=0) stop 'Error releasing memory in alkane_construct_neighbour_list'

    call system_clock(count=t2,count_rate=rate)

    write(*,'("Neighbour list constructed in : ",F15.6," seconds")')real(t2-t1,kind=ep)/real(rate,kind=ep)
    !stop

    return

  end subroutine alkane_construct_neighbour_list

  subroutine alkane_construct_linked_lists()
    !-------------------------------------------------------------------------!
    ! Assigns particles to link cells and constructs linked lists for use in  !
    ! searching over particle pairs. Closely follows procedue in appendix F   !
    ! of Frenkel and Smit.                                                    !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2010                                                  !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    use box,       only : ncellx,ncelly,ncellz,lcellx,lcelly,lcellz,use_link_cells, &
                          recip_matrix
    implicit none
    integer :: ichain,ibead,icell,ix,iy,iz,ierr
    integer :: t1,t2,rate,jbead,jchain
    real(kind=dp) :: rlcellx,rlcelly,rlcellz 
    real(kind=dp),dimension(3) :: rbead,sbead

    if (.not.use_link_cells) return

    if ( .not.allocated(head_of_cell) ) then
       allocate(head_of_cell(1:2,1:ncellx*ncelly*ncellz),stat=ierr)
       if (ierr/=0) stop 'Error allocating head_of_cell'
    end if
    if (size(head_of_cell,2)/=ncellx*ncelly*ncellz) then
       deallocate(head_of_cell,stat=ierr)
       if (ierr/=0) stop 'Error resizing head_of_cell'
       allocate(head_of_cell(1:2,1:ncellx*ncelly*ncellz),stat=ierr)
       if (ierr/=0) stop 'Error allocating head_of_cell'
    end if
    if ( .not.allocated(linked_list) ) then
       allocate(linked_list(1:4,0:nbeads,0:nchains),stat=ierr)
       if (ierr/=0) stop 'Error allocating linked_lists'
    end if

    !call system_clock(count=t1)

    head_of_cell(:,:) = 0

    rlcellx = 1.0_dp/lcellx
    rlcelly = 1.0_dp/lcelly
    rlcellz = 1.0_dp/lcellz

    do ichain = 1,nchains
       do ibead = 1,nbeads

          rbead(:) = Rchain(:,ibead,ichain)

          ! compute fractional coordinates sbead from rbead
          sbead(1) = recip_matrix(1,1)*rbead(1) + &
                     recip_matrix(2,1)*rbead(2) + &
                     recip_matrix(3,1)*rbead(3)
          sbead(2) = recip_matrix(1,2)*rbead(1) + &
                     recip_matrix(2,2)*rbead(2) + &
                     recip_matrix(3,2)*rbead(3)  
          sbead(3) = recip_matrix(1,3)*rbead(1) + &
                     recip_matrix(2,3)*rbead(2) + &
                     recip_matrix(3,3)*rbead(3) 

          sbead = sbead*0.5_dp*invPi 

          ! which link cell does this particle belong to
          ix = floor(sbead(1)*rlcellx)
          iy = floor(sbead(2)*rlcelly)
          iz = floor(sbead(3)*rlcellz)

          ix = modulo(ix,ncellx) + 1
          iy = modulo(iy,ncelly) + 1
          iz = modulo(iz,ncellz) + 1

          icell = (iz-1)*ncellx*ncelly + (iy-1)*ncellx + ix

          ! Bead and chain index for old head of cell
          ! (both zero if this is first atom to be added)
          jbead  = head_of_cell(1,icell)
          jchain = head_of_cell(2,icell)

          ! This bead points forward to the old head of cell,
          ! or to zero if it's the first bead to be added
          linked_list(1,ibead,ichain) = jbead
          linked_list(2,ibead,ichain) = jchain

          ! ..and becomes the new head of cell
          head_of_cell(1,icell) = ibead 
          head_of_cell(2,icell) = ichain 
          
          ! jbead, jchain points backward to ibead,ichain
          ! zero'th elements will be populated here for
          ! first bead added, and ignored.
          linked_list(3,jbead,jchain) = ibead
          linked_list(4,jbead,jchain) = ichain


       end do
    end do

    !call system_clock(count=t2,count_rate=rate)
    !write(*,'("Linked lists rebuilt in : ",F15.6," seconds")')real(t2-t1,kind=ep)/real(rate,kind=ep)

    return

  end subroutine alkane_construct_linked_lists

  subroutine alkane_update_linked_lists(ibead,ichain,old_pos,new_pos)
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
    integer,intent(in) :: ibead,ichain
    real(kind=dp),dimension(3),intent(in) :: old_pos,new_pos
    real(kind=dp),dimension(3) :: sbead
    integer :: ix,iy,iz,ncell,ocell,jchain,jbead
    integer :: t1,t2,rate,tmpint,kbead,kchain

    if (.not.use_link_cells) return

    ! Compute sbead from old_pos
    sbead(1) = recip_matrix(1,1)*old_pos(1) + &
               recip_matrix(2,1)*old_pos(2) + &
               recip_matrix(3,1)*old_pos(3)
    sbead(2) = recip_matrix(1,2)*old_pos(1) + &
               recip_matrix(2,2)*old_pos(2) + &
               recip_matrix(3,2)*old_pos(3)  
    sbead(3) = recip_matrix(1,3)*old_pos(1) + &
               recip_matrix(2,3)*old_pos(2) + &
               recip_matrix(3,3)*old_pos(3) 

    sbead = sbead*0.5_dp*invPi 

    ! Compute link cell number for old and new positions
    ix = floor(sbead(1)/lcellx)
    iy = floor(sbead(2)/lcelly)
    iz = floor(sbead(3)/lcellz)

    ix = modulo(ix,ncellx) + 1
    iy = modulo(iy,ncelly) + 1
    iz = modulo(iz,ncellz) + 1

    ocell = (iz-1)*ncellx*ncelly + (iy-1)*ncellx + ix

    ! Compute sbead from new_pos
    sbead(1) = recip_matrix(1,1)*new_pos(1) + &
               recip_matrix(2,1)*new_pos(2) + &
               recip_matrix(3,1)*new_pos(3)
    sbead(2) = recip_matrix(1,2)*new_pos(1) + &
               recip_matrix(2,2)*new_pos(2) + &
               recip_matrix(3,2)*new_pos(3)  
    sbead(3) = recip_matrix(1,3)*new_pos(1) + &
               recip_matrix(2,3)*new_pos(2) + &
               recip_matrix(3,3)*new_pos(3) 

    sbead = sbead*0.5_dp*invPi 


    ix = floor(sbead(1)/lcellx)
    iy = floor(sbead(2)/lcelly)
    iz = floor(sbead(3)/lcellz)

    ix = modulo(ix,ncellx) + 1
    iy = modulo(iy,ncelly) + 1
    iz = modulo(iz,ncellz) + 1

    ncell = (iz-1)*ncellx*ncelly + (iy-1)*ncellx + ix

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
    jbead   = head_of_cell(1,ocell)
    jchain  = head_of_cell(2,ocell)
    if ( (jbead==ibead).and.(jchain==ichain) ) then

       ! ibead,ichain was the old head of cell

       !----------------------------------!
       !   Before             After       !
       !   ------             -----       !
       !  hoc->i->k->..      hoc->k->..   !
       !----------------------------------!
       
       ! k used to be pointed to by i
       kbead  = linked_list(1,ibead,ichain)
       kchain = linked_list(2,ibead,ichain)

       ! k is new head of cell
       head_of_cell(1,ocell) = kbead
       head_of_cell(2,ocell) = kchain

       ! k points backward to nothing
       linked_list(3,kbead,kchain) = 0
       linked_list(4,kbead,kchain) = 0

    else
       
       !-------------------------------!
       !    Before          After      !
       !    ------          -----      !
       ! ..j->i->k->...   ..j->k->..   !
       !-------------------------------!
       kbead  = linked_list(1,ibead,ichain)
       kchain = linked_list(2,ibead,ichain)

       jbead  = linked_list(3,ibead,ichain)
       jchain = linked_list(4,ibead,ichain)

       linked_list(1,jbead,jchain) = kbead
       linked_list(2,jbead,jchain) = kchain

       linked_list(3,kbead,kchain) = jbead
       linked_list(4,kbead,kchain) = jchain


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
    jbead  = head_of_cell(1,ncell)
    jchain = head_of_cell(2,ncell)

    ! add ibead,icheain to new cell as head of cell
    linked_list(1:2,ibead,ichain) = head_of_cell(:,ncell)
    linked_list(3:4,ibead,ichain) = 0
    head_of_cell(:,ncell) = (/ibead,ichain/)

    ! jbead, jchain points backward to ibead,ichain
    linked_list(3,jbead,jchain) = ibead
    linked_list(4,jbead,jchain) = ichain 


    !call system_clock(count=t2,count_rate=rate)
    !write(*,'("Linked lists updated in : ",F15.6," seconds")')real(t2-t1,kind=ep)/real(rate,kind=ep)

    return

  end subroutine alkane_update_linked_lists



end module alkane
