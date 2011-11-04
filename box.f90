! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                               B  O  X                                       !
!=============================================================================!
!                                                                             !
! $Id: box.f90,v 1.8 2011/11/04 16:12:44 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Stores properties of the simulation 'box' (i.e. not the alkane chains) and  !
! routines to manipulate them. Also containts routines associated with the    !
! operation of link cells and periodic boundary conditions.                   !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: box.f90,v $
! Revision 1.8  2011/11/04 16:12:44  phseal
! Added box_get_num_boxes for C integration purposes.
!
! Revision 1.7  2011/10/19 16:20:17  phseal
! bypass_link_cells can now be read from the input file.
!
! Revision 1.6  2011/10/16 18:18:23  phseal
! Changed the minimum length to the side of a link cell to be an input
! parameter. Hence the second argument to box_construct_link_cells is
! no longer present, and link_cell_length is read from the input file.
!
! Revision 1.5  2011/08/02 12:56:47  phseal
! Added C bindings to all procedures which should be callable externally
! when compiled as a library.
!
! Revision 1.4  2011/08/02 12:27:18  phseal
! Updated all integers to use the integer type in constants.f90 where
! applicable. This allows the integer type it to be set to a C compatible
! type via the instrinsic iso_c_bindings module.
!
! Revision 1.3  2011/08/02 10:04:18  phseal
! Added routines to manipulate inactive boxes from outside of the box and
! alkane module. Updates overlap counting routines to return only the
! total number of overlaps found rather than lists of overlapping atoms.
!
! Revision 1.2  2011/07/29 15:58:29  phseal
! Added multiple simulation box support.
!
! Revision 1.1.1.1  2011/02/02 11:48:36  phseal
! Initial import from prototype code.
!
!
!=============================================================================!
module box

  use constants, only : dp,it
  implicit none

  private                ! Everything private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  !... unless exposed here.

  public :: box_initialise              ! Initialise this module
  public :: box_destroy                 ! Release memory used in this module
  public :: box_minimum_image           ! Minimum image routine
  public :: box_construct_link_cells    ! Build link cell structure
  public :: box_destroy_link_cells      ! Destroy link cell structure
  public :: box_update_recipmatrix      ! Compute reciprocal lattice
  public :: box_compute_volume          ! Compute volume

  public :: box_get_cell, box_set_cell  ! Manipulate hmatrix externally
  public :: box_get_num_boxes           ! Query number of boxes in use


  public :: box_cart_to_frac            ! Convert absolute coords to fractional
  public :: box_frac_to_cart            ! Convert fractional coords to absolute


  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  public :: nboxes                  ! Number of simulation boxes
  public :: hmatrix                 ! Matrix of cell vectors
  public :: recip_matrix            ! Matrix of reciprocal lattice vectors
  public :: pbc                     ! Periodic boundary conditions
  public :: use_link_cells          ! Use link cell algorithm
  public :: bypass_link_cells       ! Force link cell bypass
  public :: link_cell_length        ! Minimum length of link cell size
  public :: isotropic               ! Isotropic volume moves True/False
  public :: pressure                ! External pressure


  public :: ncellx,ncelly,ncellz    ! Number of link cells
  public :: lcellx,lcelly,lcellz    ! dimensions of link cells
  public :: lcneigh                 ! link cell topology

  public :: CellA, CellB, CellC     ! Temporary cell vectors to populate hmatrix

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  integer(kind=it),save       :: nboxes = 1        ! Number of boxes

  logical,save  :: pbc               = .true.      ! Use periodic bcs
  logical,save  :: use_link_cells    = .false.     ! Use link cells
  logical,save  :: bypass_link_cells = .false.     ! Force bypass of above

  real(kind=dp),save :: link_cell_length   = 1.5   ! Minimum link cell length in each dimension

  real(kind=dp),allocatable,dimension(:,:,:),save :: hmatrix       ! Matrix of cell vectors
  real(kind=dp),allocatable,dimension(:,:,:),save :: recip_matrix  ! Reciprocal lattice 

  logical,save       :: isotropic = .false. ! isotropic volume moves
  real(kind=dp),save :: pressure  = 6.0     ! external pressure

  ! Linked list accounting
  integer(kind=it),allocatable,dimension(:),save  :: ncellx,ncelly,ncellz   ! Number of link cells in each dir.
  real(kind=dp),allocatable,dimension(:),save     :: lcellx,lcelly,lcellz   ! Size of link cells in each dir.

  ! Array containing for each link cell, a list of the neighbouring link cells.
  integer(kind=it),allocatable,dimension(:,:,:),save :: lcneigh

  ! Temporary cell vectors to populate hmatrix
  real(kind=dp),dimension(3) :: CellA,CellB,CellC

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!
  logical :: box_initialised = .false.

contains
  
  subroutine box_initialise() bind(c)
    !------------------------------------------------------------------------------!
    ! Allocated memory for the basic properties of each simulation box.            !
    !------------------------------------------------------------------------------!
    ! D.Quigley July 2011                                                          !
    !------------------------------------------------------------------------------!   
    implicit none
    integer(kind=it),dimension(2) :: ierr
    integer(kind=it) :: ibox

    ! Allocate cell vectors and reciprocal lattice
    allocate(hmatrix(1:3,1:3,1:nboxes),stat=ierr(1))
    allocate(recip_matrix(1:3,1:3,1:nboxes),stat=ierr(2))
    if (any(ierr/=0)) stop 'Error allocating lattice vectors in box_initialise'

    ! Set all boxes the same if reading via this mechanism
    do ibox = 1,nboxes
       
       hmatrix(:,1,ibox) = CellA(:)
       hmatrix(:,2,ibox) = CellB(:)
       hmatrix(:,3,ibox) = CellC(:)

       if (pbc) then
          call box_update_recipmatrix(ibox)
       else
          recip_matrix = 0.0_dp
       end if
       
    end do

    ! Linked lists
    allocate(ncellx(1:nboxes),ncelly(1:nboxes),ncellz(1:nboxes),stat=ierr(1))
    allocate(lcellx(1:nboxes),lcelly(1:nboxes),lcellz(1:nboxes),stat=ierr(2))
    if (any(ierr/=0)) stop 'Error allocating link-cell sizes in box_initialise'

    box_initialised = .true.

    write(*,*)
    write(*,'("|=======================================|")')
    write(*,'("| Initialised ",I3," simulation box(es)      |")')nboxes
    write(*,'("|=======================================|")')
    write(*,*)

    return

  end subroutine box_initialise

  subroutine box_destroy() bind(c)
    !------------------------------------------------------------------------------!
    ! Allocated memory for the basic properties of each simulation box.            !
    !------------------------------------------------------------------------------!
    ! D.Quigley July 2011                                                          !
    !------------------------------------------------------------------------------!   
    implicit none

    deallocate(hmatrix,recip_matrix)
    deallocate(ncellx,ncelly,ncellz)
    deallocate(lcellx,lcelly,lcellz)

    return

  end subroutine box_destroy


  real(kind=dp) function box_compute_volume(ibox) bind(c)
    !------------------------------------------------------------------------------!
    ! Computes the determinant of a 3x3 matrix.                                    !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!     
    implicit none
    integer(kind=it),intent(in) :: ibox
    real(kind=dp) :: det

    Det =       hmatrix(1,1,ibox)*(hmatrix(2,2,ibox)*hmatrix(3,3,ibox) - &
                                   hmatrix(2,3,ibox)*hmatrix(3,2,ibox))
    Det = Det - hmatrix(1,2,ibox)*(hmatrix(2,1,ibox)*hmatrix(3,3,ibox) - &
                                   hmatrix(2,3,ibox)*hmatrix(3,1,ibox))
    Det = Det + hmatrix(1,3,ibox)*(hmatrix(2,1,ibox)*hmatrix(3,2,ibox) - &
                                   hmatrix(2,2,ibox)*hmatrix(3,1,ibox))

    box_compute_volume = abs(det)

    return

  end function box_compute_volume

  subroutine box_update_recipmatrix(ibox) bind(c)
    !------------------------------------------------------------------------------!
    ! Calculates the matrix of reciprocal lattive vectors from the hmatrix         !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------! 
    use constants, only : Pi
    implicit none
    integer(kind=it),intent(in) :: ibox
    real(kind=dp) :: vol

    ! invert hmatrix to get recip_matrix
    recip_matrix(1,1,ibox)=hmatrix(2,2,ibox)*hmatrix(3,3,ibox)-hmatrix(2,3,ibox)*hmatrix(3,2,ibox)
    recip_matrix(1,2,ibox)=hmatrix(2,3,ibox)*hmatrix(3,1,ibox)-hmatrix(2,1,ibox)*hmatrix(3,3,ibox)
    recip_matrix(1,3,ibox)=hmatrix(2,1,ibox)*hmatrix(3,2,ibox)-hmatrix(2,2,ibox)*hmatrix(3,1,ibox)

    recip_matrix(2,1,ibox)=hmatrix(1,3,ibox)*hmatrix(3,2,ibox)-hmatrix(1,2,ibox)*hmatrix(3,3,ibox)
    recip_matrix(2,2,ibox)=hmatrix(1,1,ibox)*hmatrix(3,3,ibox)-hmatrix(1,3,ibox)*hmatrix(3,1,ibox)
    recip_matrix(2,3,ibox)=hmatrix(1,2,ibox)*hmatrix(3,1,ibox)-hmatrix(1,1,ibox)*hmatrix(3,2,ibox)

    recip_matrix(3,1,ibox)=hmatrix(1,2,ibox)*hmatrix(2,3,ibox)-hmatrix(1,3,ibox)*hmatrix(2,2,ibox)
    recip_matrix(3,2,ibox)=hmatrix(1,3,ibox)*hmatrix(2,1,ibox)-hmatrix(1,1,ibox)*hmatrix(2,3,ibox)
    recip_matrix(3,3,ibox)=hmatrix(1,1,ibox)*hmatrix(2,2,ibox)-hmatrix(1,2,ibox)*hmatrix(2,1,ibox)
    
    ! Calculte cell volume
    vol =hmatrix(1,1,ibox)*recip_matrix(1,1,ibox) + &
         hmatrix(1,2,ibox)*recip_matrix(1,2,ibox) + &
         hmatrix(1,3,ibox)*recip_matrix(1,3,ibox)

    ! Scale reciprocal lattice by 2*pi/volume
    recip_matrix(:,:,ibox)=recip_matrix(:,:,ibox)*2.0_dp*Pi/vol

    return

  end subroutine box_update_recipmatrix


  function box_minimum_image(r1,r2,ibox)
    !-------------------------------------------------------------------------!
    ! Returns the shortest vector connecting images of particles located at   !
    ! r1 and r2. Should be in-lined where used for efficiency.                !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    implicit none
    real(kind=dp),dimension(3),intent(in) :: r1,r2
    integer(kind=it),intent(in) :: ibox
    real(kind=dp),dimension(3)  :: box_minimum_image
    real(kind=dp),dimension(3)  :: dr
    real(kind=dp) :: sx,sy,sz

    dr(1) = r2(1) - r1(1)
    dr(2) = r2(2) - r1(2)
    dr(3) = r2(3) - r1(3)

    if (.not.pbc) then
       box_minimum_image = dr(:)
       return   ! Use raw distance
    end if

    sx = recip_matrix(1,1,ibox)*dr(1) + &
         recip_matrix(2,1,ibox)*dr(2) + &
         recip_matrix(3,1,ibox)*dr(3)
    sy = recip_matrix(1,2,ibox)*dr(1) + &
         recip_matrix(2,2,ibox)*dr(2) + &
         recip_matrix(3,2,ibox)*dr(3)  
    sz = recip_matrix(1,3,ibox)*dr(1) + &
         recip_matrix(2,3,ibox)*dr(2) + &
         recip_matrix(3,3,ibox)*dr(3) 

    sx = sx*0.5_dp*invPi 
    sy = sy*0.5_dp*invPi
    sz = sz*0.5_dp*invPi 

    ! apply boundary conditions
    sx = sx - floor(sx+0.5_dp,kind=dp)
    sy = sy - floor(sy+0.5_dp,kind=dp)
    sz = sz - floor(sz+0.5_dp,kind=dp)
    
    ! scale back up
    dr(1) = hmatrix(1,1,ibox)*sx + &
            hmatrix(1,2,ibox)*sy + &
            hmatrix(1,3,ibox)*sz
                       
    dr(2) = hmatrix(2,1,ibox)*sx + &
            hmatrix(2,2,ibox)*sy + &
            hmatrix(2,3,ibox)*sz
                       
    dr(3) = hmatrix(3,1,ibox)*sx + &
            hmatrix(3,2,ibox)*sy + &
            hmatrix(3,3,ibox)*sz

!!$    dr(1) = dr(1) - Lx*anint(dr(1)*rLx)
!!$    dr(2) = dr(2) - Ly*anint(dr(2)*rLy)
!!$    dr(3) = dr(3) - Lz*anint(dr(3)*rLz)

    box_minimum_image = dr(:)
    
    return

  end function box_minimum_image

  subroutine box_destroy_link_cells() bind(c)
    !-------------------------------------------------------------------------!
    ! Released memory used by the link-cell algorithm                         !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it) :: ierr

    if (use_link_cells) then
       deallocate(lcneigh,stat=ierr)
    end if

  end subroutine box_destroy_link_cells

  subroutine box_construct_link_cells(ibox) bind(c)
    !-------------------------------------------------------------------------!
    ! Analyses the dimensions of the simulation cell and determines if use of !
    ! a link-cell algorithm is possible. The module level flag use_link_cells !
    ! is set appropriately. If link cells are possible then the number and    !
    ! dimension of cells is determined and data on the link-cell topology is  !
    ! constructed.                                                            !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(in) :: ibox
    real(kind=dp) :: Lx,Ly,Lz

    integer(kind=it) :: ix,iy,iz,jx,jy,jz,icell,jcell
    integer(kind=it) :: kx,ky,kz,jn
    integer(kind=it) :: ierr

    ! Sanity check
    if (.not.box_initialised) then
       write(0,'("Error in box_construct_link_cells : box module not initialised.")')
       stop
    endif

    use_link_cells = .true.

    ! Bomb out if not a periodic system
    if ((.not.pbc).or.bypass_link_cells) then
       use_link_cells = .false.
       return
    end if

    ! Lengths of the three cell vectors
    Lx = sqrt(dot_product(hmatrix(:,1,ibox),hmatrix(:,1,ibox)))
    Ly = sqrt(dot_product(hmatrix(:,2,ibox),hmatrix(:,2,ibox)))
    Lz = sqrt(dot_product(hmatrix(:,3,ibox),hmatrix(:,3,ibox)))

    ! Number of link cells in each direction
    ncellx(ibox) = int(Lx/link_cell_length)
    ncelly(ibox) = int(Ly/link_cell_length)
    ncellz(ibox) = int(Lz/link_cell_length)

    ! Fractional size of link-cells in each cell dimension
    lcellx(ibox) = 1.0_dp/real(ncellx(ibox),kind=dp)
    lcelly(ibox) = 1.0_dp/real(ncelly(ibox),kind=dp)
    lcellz(ibox) = 1.0_dp/real(ncellz(ibox),kind=dp)


    ! Bomb out if system is too small for link cells
    if ( (ncellx(ibox)<4).or.(ncelly(ibox)<4).or.(ncellz(ibox)<4) ) then
       use_link_cells = .false.
       write(0,'("System has become too small for link cells of minimum length ",F15.6)')link_cell_length
       return
    end if

    if (allocated(lcneigh)) deallocate(lcneigh)
    allocate(lcneigh(1:27,1:ncellx(ibox)*ncelly(ibox)*ncellz(ibox),1:nboxes),stat=ierr)
    if (ierr/=0) stop 'Error allocating link-cell neighbour array'

    ! Decide which cells are neighbours of each other.
    do iz = 1,ncellz(ibox)
       do iy = 1,ncelly(ibox)
          do ix = 1,ncellx(ibox)       

             icell = (iz-1)*ncellx(ibox)*ncelly(ibox) + (iy-1)*ncellx(ibox) + ix           

             jn = 1
             do jz = iz-1,iz+1
                do jy = iy-1,iy+1
                   do jx = ix-1,ix+1

                      kx = mod(jx,ncellx(ibox)) ; if (kx==0) kx=ncellx(ibox)
                      ky = mod(jy,ncelly(ibox)) ; if (ky==0) ky=ncelly(ibox)
                      kz = mod(jz,ncellz(ibox)) ; if (kz==0) kz=ncellz(ibox)                    

                      jcell = (kz-1)*ncellx(ibox)*ncelly(ibox) + (ky-1)*ncellx(ibox) + kx    

                      lcneigh(jn,icell,ibox) = jcell

                      jn = jn + 1

                   end do
                end do
             end do

          end do
       end do
    end do

    return

  end subroutine box_construct_link_cells

  subroutine box_get_cell(ibox,dumhmatrix) bind(c)
    !-------------------------------------------------------------------------!
    ! Queries the current matrix of cell vectors for box ibox, and returns    !
    ! via the dummy 3x3 matrix dumhmatrix.                                    !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(in) :: ibox
    real(kind=dp),dimension(3,3),intent(out) :: dumhmatrix

    if (ibox > nboxes ) stop 'Error in box_get_cell, ibox > nboxes'

    dumhmatrix = hmatrix(:,:,ibox)

    return

  end subroutine box_get_cell

  subroutine box_get_num_boxes(dumnb) bind (c)
    !-------------------------------------------------------------------------!
    ! Queries the number of simulation boxes in use as set when the input     !
    ! file is read. Provided for use from C when code compiled as library.    !
    !-------------------------------------------------------------------------!
    ! D.Quigley November 2011                                                 !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(out) :: dumnb

    dumnb = nboxes

    return

  end subroutine box_get_num_boxes

  subroutine box_set_cell(ibox,dumhmatrix) bind(c)
    !-------------------------------------------------------------------------!
    ! Sets the current matrix of cell vectors for box ibox. WARNING, this     !
    ! routine does not update the matrix of reciprocal lattice vectors, or    !
    ! rescale the coordinates of any atoms/molecules within the box.          !
    !                                                                         !
    ! Please consider as primarily for debugging purposes.                    !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    integer(kind=it),intent(in) :: ibox
    real(kind=dp),dimension(3,3),intent(in) :: dumhmatrix

    if (ibox > nboxes ) stop 'Error in box_get_cell, ibox > nboxes'

    hmatrix(:,:,ibox) = dumhmatrix

    return

  end subroutine box_set_cell

  subroutine box_cart_to_frac(ibox,in_vector,out_vector) bind(c)
    !-------------------------------------------------------------------------!
    ! Fox the specified box (ibox) convert in_vector in absolute cartesian    !
    ! coords into internal box-scaled coords.                                 !
    !                                                                         !
    ! Requirements : recip_matrix must be valid and current for ibox          !
    !                use box_update_recip_matrix(ibox) if not.                !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    implicit none
    integer(kind=it),intent(in) :: ibox
    real(kind=dp),dimension(3),intent(in)  :: in_vector
    real(kind=dp),dimension(3),intent(out) :: out_vector

    out_vector(1) = recip_matrix(1,1,ibox)*in_vector(1) + &
                    recip_matrix(2,1,ibox)*in_vector(2) + &
                    recip_matrix(3,1,ibox)*in_vector(3)
    out_vector(2) = recip_matrix(1,2,ibox)*in_vector(1) + &
                    recip_matrix(2,2,ibox)*in_vector(2) + &
                    recip_matrix(3,2,ibox)*in_vector(3)  
    out_vector(3) = recip_matrix(1,3,ibox)*in_vector(1) + &
                    recip_matrix(2,3,ibox)*in_vector(2) + &
                    recip_matrix(3,3,ibox)*in_vector(3) 

    out_vector(1) = out_vector(1)*0.5_dp*invPi 
    out_vector(2) = out_vector(2)*0.5_dp*invPi
    out_vector(3) = out_vector(3)*0.5_dp*invPi 

    return

  end subroutine box_cart_to_frac

  subroutine box_frac_to_cart(ibox,in_vector,out_vector) bind(c)
    !-------------------------------------------------------------------------!
    ! Fox the specified box (ibox) convert in_vector in box-scaled coords     !
    ! into out_vector in absolute cartesian coords.                           !
    !                                                                         !
    ! Requirements : hmatrix must be valid and current for ibox               !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2011                                                   !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    implicit none
    integer(kind=it),intent(in) :: ibox
    real(kind=dp),dimension(3),intent(in)  :: in_vector
    real(kind=dp),dimension(3),intent(out) :: out_vector

    out_vector(1) = hmatrix(1,1,ibox)*in_vector(1) + &
                    hmatrix(1,2,ibox)*in_vector(2) + &
                    hmatrix(1,3,ibox)*in_vector(3)
                       
    out_vector(2) = hmatrix(2,1,ibox)*in_vector(1) + &
                    hmatrix(2,2,ibox)*in_vector(2) + &
                    hmatrix(2,3,ibox)*in_vector(3)
                       
    out_vector(3) = hmatrix(3,1,ibox)*in_vector(1) + &
                    hmatrix(3,2,ibox)*in_vector(2) + &
                    hmatrix(3,3,ibox)*in_vector(3)

    return

  end subroutine box_frac_to_cart


end module box
