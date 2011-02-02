! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                               B  O  X                                       !
!=============================================================================!
!                                                                             !
! $Id: box.f90,v 1.1 2011/02/02 11:48:36 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Stores properties of the simulation 'box' (i.e. not the alkane chains) and  !
! routines to manipulate them. Also containts routines associated with the    !
! operation of link cells and periodic boundary conditions.                   !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: box.f90,v $
! Revision 1.1  2011/02/02 11:48:36  phseal
! Initial revision
!
!
!=============================================================================!
module box

  use constants, only : dp
  implicit none

  private                ! Everything private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  !... unless exposed here.


  public :: box_minimum_image           ! Minimum image routine
  public :: box_construct_link_cells    ! Build link cell structure
  public :: box_destroy_link_cells      ! Destroy link cell structure
  public :: box_update_recipmatrix      ! Compute reciprocal lattice
  public :: box_compute_volume          ! Compute volume


  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  public :: hmatrix                 ! Matrix of cell vectors
  public :: recip_matrix            ! Matrix of reciprocal lattice vectors
  public :: pbc                     ! Periodic boundary conditions
  public :: use_link_cells          ! Use link cell algorithm
  public :: isotropic               ! Isotropic volume moves True/False
  public :: pressure                ! External pressure


  public :: ncellx,ncelly,ncellz    ! Number of link cells
  public :: lcellx,lcelly,lcellz    ! dimensions of link cells
  public :: lcneigh                 ! link cell topology

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  logical,save       :: pbc               = .true.   ! Use periodic bcs
  logical,save       :: use_link_cells    = .false.  ! Use link cells
  logical,save       :: bypass_link_cells = .false.  ! Force bypass of above

  real(kind=dp),dimension(3,3),save :: hmatrix       ! Matrix of cell vectors
  real(kind=dp),dimension(3,3),save :: recip_matrix  ! Reciprocal lattice 
  real(kind=dp),save :: rcut

  logical,save       :: isotropic = .false. ! isotropic volume moves
  real(kind=dp),save :: pressure  = 6.0     ! external pressure

  ! Linked list accounting
  integer,save       :: ncellx,ncelly,ncellz   ! Number of link cells in each dir.
  real(kind=dp),save :: lcellx,lcelly,lcellz   ! Size of link cells in each dir.

  ! Array containing for each link cell, a list of the neighbouring link cells.
  integer,allocatable,dimension(:,:),save :: lcneigh

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

contains
  
  real(kind=dp) function box_compute_volume()
    !------------------------------------------------------------------------------!
    ! Computes the determinant of a 3x3 matrix.                                    !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!     
    implicit none
    real(kind=dp) :: det

    Det =       hmatrix(1,1)*(hmatrix(2,2)*hmatrix(3,3) - &
                              hmatrix(2,3)*hmatrix(3,2))
    Det = Det - hmatrix(1,2)*(hmatrix(2,1)*hmatrix(3,3) - &
                              hmatrix(2,3)*hmatrix(3,1))
    Det = Det + hmatrix(1,3)*(hmatrix(2,1)*hmatrix(3,2) - &
                              hmatrix(2,2)*hmatrix(3,1))

    box_compute_volume = abs(det)

    return

  end function box_compute_volume

  subroutine box_update_recipmatrix()
    !------------------------------------------------------------------------------!
    ! Calculates the matrix of reciprocal lattive vectors from the hmatrix         !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------! 
    use constants, only : Pi
    implicit none
    real(kind=dp) :: vol

    ! invert hmatrix to get recip_matrix
    recip_matrix(1,1)=hmatrix(2,2)*hmatrix(3,3)-hmatrix(2,3)*hmatrix(3,2)
    recip_matrix(1,2)=hmatrix(2,3)*hmatrix(3,1)-hmatrix(2,1)*hmatrix(3,3)
    recip_matrix(1,3)=hmatrix(2,1)*hmatrix(3,2)-hmatrix(2,2)*hmatrix(3,1)
    
    recip_matrix(2,1)=hmatrix(1,3)*hmatrix(3,2)-hmatrix(1,2)*hmatrix(3,3)
    recip_matrix(2,2)=hmatrix(1,1)*hmatrix(3,3)-hmatrix(1,3)*hmatrix(3,1)
    recip_matrix(2,3)=hmatrix(1,2)*hmatrix(3,1)-hmatrix(1,1)*hmatrix(3,2)
    
    recip_matrix(3,1)=hmatrix(1,2)*hmatrix(2,3)-hmatrix(1,3)*hmatrix(2,2)
    recip_matrix(3,2)=hmatrix(1,3)*hmatrix(2,1)-hmatrix(1,1)*hmatrix(2,3)
    recip_matrix(3,3)=hmatrix(1,1)*hmatrix(2,2)-hmatrix(1,2)*hmatrix(2,1)
    
    ! Calculte cell volume
    vol =hmatrix(1,1)*recip_matrix(1,1) + &
         hmatrix(1,2)*recip_matrix(1,2) + &
         hmatrix(1,3)*recip_matrix(1,3)

    ! Scale reciprocal lattice by 2*pi/volume
    recip_matrix(:,:)=recip_matrix(:,:)*2.0_dp*Pi/vol

    return

  end subroutine box_update_recipmatrix


  function box_minimum_image(r1,r2)
    !-------------------------------------------------------------------------!
    ! Returns the shortest vector connecting images of particles located at   !
    ! r1 and r2. Should be in-lined where used for efficiency.                !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    use constants, only : invPi
    implicit none
    real(kind=dp),dimension(3),intent(in) :: r1,r2
    real(kind=dp),dimension(3) :: box_minimum_image
    real(kind=dp),dimension(3) :: dr
    real(kind=dp) :: sx,sy,sz

    dr(1) = r2(1) - r1(1)
    dr(2) = r2(2) - r1(2)
    dr(3) = r2(3) - r1(3)

    if (.not.pbc) then
       box_minimum_image = dr(:)
       return   ! Use raw distance
    end if

    sx = recip_matrix(1,1)*dr(1) + &
         recip_matrix(2,1)*dr(2) + &
         recip_matrix(3,1)*dr(3)
    sy = recip_matrix(1,2)*dr(1) + &
         recip_matrix(2,2)*dr(2) + &
         recip_matrix(3,2)*dr(3)  
    sz = recip_matrix(1,3)*dr(1) + &
         recip_matrix(2,3)*dr(2) + &
         recip_matrix(3,3)*dr(3) 

    sx = sx*0.5_dp*invPi 
    sy = sy*0.5_dp*invPi
    sz = sz*0.5_dp*invPi 

    ! apply boundary conditions
    sx = sx - floor(sx+0.5_dp,kind=dp)
    sy = sy - floor(sy+0.5_dp,kind=dp)
    sz = sz - floor(sz+0.5_dp,kind=dp)
    
    ! scale back up
    dr(1) = hmatrix(1,1)*sx + &
            hmatrix(1,2)*sy + &
            hmatrix(1,3)*sz

    dr(2) = hmatrix(2,1)*sx + &
            hmatrix(2,2)*sy + &
            hmatrix(2,3)*sz

    dr(3) = hmatrix(3,1)*sx + &
            hmatrix(3,2)*sy + &
            hmatrix(3,3)*sz

!!$    dr(1) = dr(1) - Lx*anint(dr(1)*rLx)
!!$    dr(2) = dr(2) - Ly*anint(dr(2)*rLy)
!!$    dr(3) = dr(3) - Lz*anint(dr(3)*rLz)

    box_minimum_image = dr(:)
    
    return

  end function box_minimum_image

  subroutine box_destroy_link_cells()
    !-------------------------------------------------------------------------!
    ! Released memory used by the link-cell algorithm                         !
    !-------------------------------------------------------------------------!
    ! D.Quigley February 2010                                                 !
    !-------------------------------------------------------------------------!
    implicit none
    integer :: ierr

    if (use_link_cells) then
       deallocate(lcneigh,stat=ierr)
    end if

  end subroutine box_destroy_link_cells

  subroutine box_construct_link_cells(drcut)
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
    real(kind=dp),intent(in) :: drcut
    real(kind=dp) :: Lx,Ly,Lz

    integer :: ix,iy,iz,jx,jy,jz,icell,jcell
    integer :: kx,ky,kz,jn
    integer :: ierr

    ! Local copy
    rcut = drcut

    use_link_cells = .true.

    ! Bomb out if not a periodic system
    if ((.not.pbc).or.bypass_link_cells) then
       use_link_cells = .false.
       return
    end if

    ! Lengths of the three cell vectors
    Lx = sqrt(dot_product(hmatrix(:,1),hmatrix(:,1)))
    Ly = sqrt(dot_product(hmatrix(:,2),hmatrix(:,2)))
    Lz = sqrt(dot_product(hmatrix(:,3),hmatrix(:,3)))

    ! Number of link cells in each direction
    ncellx = int(Lx/rcut)
    ncelly = int(Ly/rcut)
    ncellz = int(Lz/rcut)

    ! Fractional size of link-cells in each cell dimension
    lcellx = 1.0_dp/real(ncellx,kind=dp)
    lcelly = 1.0_dp/real(ncelly,kind=dp)
    lcellz = 1.0_dp/real(ncellz,kind=dp)

    ! Bomb out if system is too small for link cells
    if ( (ncellx<4).or.(ncelly<4).or.(ncellz<4) ) then
       use_link_cells = .false.
       write(0,'("System has become too small for link cells")')
       return
    end if

    if (allocated(lcneigh)) deallocate(lcneigh)
    allocate(lcneigh(1:27,1:ncellx*ncelly*ncellz),stat=ierr)
    if (ierr/=0) stop 'Error allocating link-cell neighbour array'

    ! Decide which cells are neighbours of each other.
    do iz = 1,ncellz
       do iy = 1,ncelly
          do ix = 1,ncellx       

             icell = (iz-1)*ncellx*ncelly + (iy-1)*ncellx + ix           

             jn = 1
             do jz = iz-1,iz+1
                do jy = iy-1,iy+1
                   do jx = ix-1,ix+1

                      kx = mod(jx,ncellx) ; if (kx==0) kx=ncellx
                      ky = mod(jy,ncelly) ; if (ky==0) ky=ncelly
                      kz = mod(jz,ncellz) ; if (kz==0) kz=ncellz                    


                      jcell = (kz-1)*ncellx*ncelly + (ky-1)*ncellx + kx    

                      lcneigh(jn,icell) = jcell

                      jn = jn + 1

                   end do
                end do
             end do

          end do
       end do
    end do

    return

  end subroutine box_construct_link_cells

end module box
