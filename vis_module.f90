! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                               V  I  S                                       !
!=============================================================================!
!                                                                             !
! $Id: vis_module.f90,v 1.1 2011/02/02 11:48:36 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Routines to create psf and dcd files of alkane chains for visualisation.    !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: vis_module.f90,v $
! Revision 1.1  2011/02/02 11:48:36  phseal
! Initial revision
!
!
!=============================================================================!
module vis

  use constants, only : ep,dp
  implicit none

  private  ! Everything private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  !... unless exposed here.

  public :: write_psf            ! Create appropraite protein structure file
  public :: write_dcd_header     ! Write a header for a binary trajectory file
  public :: write_dcd_snapshot   ! Write a snapshot to a binary trajectory file

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------
  integer,save :: psf = 201
  integer,save :: dcd = 202

  contains

    subroutine write_psf(nbeads,nchains)
      !------------------------------------------------------!
      ! Writes a VMD compatible psf file for a linear chain  !
      ! of bonded beads. The symbol for carbon is used for   !
      ! each bead, this is an arbitary choice. The radius    !
      ! of the beads can be changed once loaded into vmd.    !
      !------------------------------------------------------!
      ! D.Quigley January 2010                               !
      !------------------------------------------------------! 
      implicit none
      integer,intent(in) :: nchains,nbeads

      integer,allocatable,dimension(:) :: arrbonds

      integer :: i,ierr,j,k,ichain

      ! open the psf file
      open(unit=psf,file='chain.psf',status='replace',iostat=ierr)
      if (ierr/=0) stop 'Error opening chain.psf for output'

      ! write the header
      write(psf,'(A3)')'PSF'
      write(psf,'("         1 !NTITLE")')
      write(psf,*)
      write(psf,'(I8,1x,"!NATOM")')nchains*nbeads

      ! format for x-plor style psf
10    format(I8,1x,A4,1x,I4,1x,A4,1x,A4,1x,A5,1x,F10.6,1x,5x,F8.4,10x,"0")

      ! Write the atoms information. The last entry is the mass
      ! which will be ignored by VMD anyway.
      k = 1
      do ichain = 1,nchains
         do j = 1,nbeads
            write(psf,10)k,"BULK",ichain,"UNK ","C","C",0.0,1.0
            k = k + 1
         end do
      end do

      ! Write the number of bonds
      write(psf,*)
      write(psf,'(I8,1x,"!NBOND: bonds")')nchains*(nbeads-1)


      ! Construct the array of bonding information
      allocate(arrbonds(1:2*nchains*(nbeads-1)),stat=ierr)
      if (ierr/=0) stop 'Error allocating bonds array in write_psf'
      k = 1
      do ichain = 1,nchains
         do i = 1,Nbeads-1
            arrbonds(k)   = (ichain-1)*nbeads+i
            arrbonds(k+1) = (ichain-1)*nbeads+i+1
            k = k + 2
         end do
      end do

      ! Write the bonding information in the correct format
30    format(1x,I7,1x,I7,1x,I7,1x,I7,1x,I7,1x,I7,1x,I7,1x,I7,1x)
      write(psf,30)arrbonds

      ! write irrelevant stuff - chimera will check for it and
      ! complain if it isn't there.
      write(psf,'(I8,1x,"!NTHETA: angles")')0
      write(psf,'(I8,1x,"!NPHI: torsions")')0
      write(psf,'(I8,1x,"!NIMPHI: torsions")')0    
      write(psf,'(I8,1x,"!NDON: donors")')0 
      write(psf,'(I8,1x,"!NACC: acceptors")')0 

      close(psf)
      deallocate(arrbonds)

      return

    end subroutine write_psf

  subroutine write_dcd_header(Nchains,Nbeads)
    !------------------------------------------------------!
    ! Writes the header of a VMD compatible dcd file for a !
    ! linear chain of bonded beads.                        !
    !------------------------------------------------------!
    ! D.Quigley January 2010                               !
    !------------------------------------------------------! 
    implicit none
    integer,intent(in) :: Nchains,Nbeads

    ! arrays for header
    integer,dimension(20) :: icntrl
    character(4) :: hdr='CORD'
    character*80,dimension(32) :: dcdtitle

    integer :: i,ierr

    open(unit=dcd,file='chain.dcd',status='replace',iostat=ierr,form='unformatted')
    if (ierr/=0) stop 'Error opening chain.dcd file - quitting'

    ! write the dcd header - most of this will be ignored
    icntrl(1)     = 1000              ! number of snapshots in history file
    icntrl(2)     = 0
    icntrl(3)     = 100               ! gap in steps between snapshots (doesn't matter)
    icntrl(4)     = 100*1000          ! total numbe of steps (VMD ignores this)
    icntrl(5:7)   = 0
    icntrl(8)     = 3*Nchains*Nbeads  ! Ndeg 
    icntrl(9)     = 0                 ! no fixed atoms
    icntrl(10)    = 0
    icntrl(11)    = 1                 ! 1/0 for unit cell presence
    icntrl(12:19) = 0
    icntrl(20)    = 24                ! Charmm version number (fixes dcd format)

    write(dcd)hdr,icntrl
    write(dcd)1,(dcdtitle(i),i=1,1)
    write(dcd)Nchains*nbeads

    close(dcd)

    return

  end subroutine write_dcd_header

  subroutine write_dcd_snapshot(Nchains,nbeads,r)
    !------------------------------------------------------!
    ! Writes a snapshot of the current positions to a the  !
    ! dcd file. Expects a 2D array r(1:3,1:nchain) in      !
    ! double precision which holds the coordinates.        !
    !------------------------------------------------------!
    ! D.Quigley January 2010                               !
    !------------------------------------------------------! 
    use constants, only : invPi
    use box,       only : hmatrix,pbc,recip_matrix
    implicit none

    integer,intent(in) :: Nchains,nbeads
    real(kind=dp),dimension(1:3,1:Nbeads,1:Nchains),intent(in) :: r
    real(kind=dp),allocatable,dimension(:,:,:) :: rcopy
    real(kind=dp),dimension(3) :: rcom,oldcom,tmpcom,comchain
    real(kind=dp),dimension(3) :: unita,unitb,unitc

    ! charmm style cell vector array
    real(kind=ep),dimension(6) :: xtlabc   

    integer :: ierr,j,i,ichain,ibead

    allocate(rcopy(1:3,1:nbeads,1:nchains),stat=ierr)

!!$    do ichain = 1,nchains
!!$       do j = 1,nbeads
!!$          rcopy(1,j,ichain) = r(1,j,ichain) - Lx*anint(r(1,j,ichain)*rLx)
!!$          rcopy(2,j,ichain) = r(2,j,ichain) - Ly*anint(r(2,j,ichain)*rLy)
!!$          rcopy(3,j,ichain) = r(3,j,ichain) - Lz*anint(r(3,j,ichain)*rLz)
!!$       end do
!!$    end do


!!$    rcopy = r

    do ichain = 1,nchains
       
       ! Find center of mass and the vector which translates
       ! it back inside the unit cell.
       comchain(:) = 0.0_dp
       do ibead = 1,nbeads
          comchain(:) = comchain(:) + r(:,ibead,ichain)
       end do
       comchain(:) = comchain(:)/real(nbeads,kind=dp)
       oldcom(:)   = comchain(:)

       tmpcom(1) = recip_matrix(1,1)*oldcom(1) + &
                   recip_matrix(2,1)*oldcom(2) + &
                   recip_matrix(3,1)*oldcom(3)
       tmpcom(2) = recip_matrix(1,2)*oldcom(1) + &
                   recip_matrix(2,2)*oldcom(2) + &
                   recip_matrix(3,2)*oldcom(3)  
       tmpcom(3) = recip_matrix(1,3)*oldcom(1) + &
                   recip_matrix(2,3)*oldcom(2) + &
                   recip_matrix(3,3)*oldcom(3) 

       tmpcom    = tmpcom*0.5_dp*invPi 
       
       oldcom(1) = - anint(tmpcom(1)) 
       oldcom(2) = - anint(tmpcom(2)) 
       oldcom(3) = - anint(tmpcom(3)) 

       tmpcom(1) = hmatrix(1,1)*oldcom(1) + &
                   hmatrix(1,2)*oldcom(2) + &
                   hmatrix(1,3)*oldcom(3)

       tmpcom(2) = hmatrix(2,1)*oldcom(1) + &
                   hmatrix(2,2)*oldcom(2) + &
                   hmatrix(2,3)*oldcom(3)

       tmpcom(3) = hmatrix(3,1)*oldcom(1) + &
                   hmatrix(3,2)*oldcom(2) + &
                   hmatrix(3,3)*oldcom(3)

       ! Apply to all beads in a chain
       do ibead = 1,nbeads
          rcopy(:,ibead,ichain) = r(:,ibead,ichain) + tmpcom(:)
       end do

    end do


    ! Special treatment for single chain with no pbc
    if ( (nchains==1).and.(.not.pbc) ) then
       rcom = 0.0_dp
       do i = 1,nbeads
          rcom(:) = rcom(:) + r(:,i,1)
       end do
       rcom(:) = rcom(:)/real(nbeads,kind=dp)
       do i = 1,nbeads
          rcopy(:,i,1) = r(:,i,1) - rcom(:)
       end do
    end if

    open(unit=dcd,file='chain.dcd',status='old',position='append',iostat=ierr,form='unformatted')
    if (ierr/=0) stop 'Error opening chain.dcd file - quitting'

    ! Specify an abritary periodic box. It doesn't matter that we aren't using
    ! periodic boundaries - vmd will ignore this information unless we ask it
    ! to do otherwise.
    xtlabc(1) = sqrt(dot_product(hmatrix(:,1),hmatrix(:,1)))
    xtlabc(3) = sqrt(dot_product(hmatrix(:,2),hmatrix(:,2)))
    xtlabc(6) = sqrt(dot_product(hmatrix(:,3),hmatrix(:,3)))

    unita(:)  = hmatrix(:,1)/xtlabc(1)
    unitb(:)  = hmatrix(:,2)/xtlabc(3)
    unitc(:)  = hmatrix(:,3)/xtlabc(6)

    xtlabc(2) = acos(dot_product(unita,unitb))*180.0*invPi
    xtlabc(4) = acos(dot_product(unita,unitc))*180.0*invPi
    xtlabc(5) = acos(dot_product(unitb,unitc))*180.0*invPi

    ! Write the information to file, note conversion to single
    ! precision to save on file size.
    write(dcd)xtlabc
    write(dcd)((-real(rcopy(1,i,j),kind=4),i=1,Nbeads),j=1,nchains)
    write(dcd)((-real(rcopy(2,i,j),kind=4),i=1,Nbeads),j=1,nchains)
    write(dcd)((real(rcopy(3,i,j),kind=4),i=1,Nbeads),j=1,nchains)

    !call flush(dcd)
    close(dcd)

    deallocate(rcopy,stat=ierr)
    

  end subroutine write_dcd_snapshot

end module vis
