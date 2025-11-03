!=============================================================================!
!                              I   O                                          !
!=============================================================================!
!                                                                             !
! $Id: io.f90,v 1.11 2014/05/06 14:39:29 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Holds routines to read the main input file, the xmol file containing        !
! initial coordinates, and variables relating to IO choices.                  !
!-----------------------------------------------------------------------------!
module io

  use iso_c_binding
  use constants, only : dp,ep,it
  implicit none

  private                ! Everything private


  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.

  public :: io_read_input                             ! Read input file
  public :: io_read_xmol                              ! Read initial coords
  public :: io_write_xmol                             ! Write current coords
  
  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  public :: read_xmol                                 ! Flag to read xmol
  public :: file_output_int                           ! Sample interval
  public :: traj_output_int                           ! Trajectory interval
  public :: io_standalone                             ! Not running as library?

  logical :: read_xmol = .false.   ! are we reading an xmol

  ! Sampling and snapshot intervals
  integer(kind=it) :: file_output_int = 25
  integer(kind=it) :: traj_output_int = 250

  logical,save :: io_standalone = .false.
  
  !---------------------------------------------------------------------------
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  ! Base name of xmol file from which structure is read
  character(len=60) :: basefilename = "chain.xmol"

  ! Base name of xmol file to write structure is written
  character(len=60) :: outfilename = "final.xmol"

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

contains

  subroutine io_read_input() bind(c)
    !-------------------------------------------------------------------------!
    ! Reads the input file specified as the first command line argument. Uses !
    ! Fortran nameslists, and populates user specified variables in box,      !
    ! alkane and timer.                                                       !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2011                                                  !
    !-------------------------------------------------------------------------!
    use alkane, only : nchains,nbeads,sigma,L,model_type,torsion_type,rigid,max_regrow
    use box,    only : pbc,isotropic,pressure, &
                       box_update_recipmatrix,nboxes,CellA,CellB,CellC, &
                       link_cell_length,bypass_link_cells,use_verlet_list
    use mc,     only : eq_mc_cycles,max_mc_cycles,eq_adjust_mc,mc_target_ratio,eq_adjust_int
    use timer,  only : timer_closetime,timer_qtime
    implicit none


    namelist/system/nboxes,nchains,nbeads,CellA,CellB,CellC,sigma,L,model_type, &
                    torsion_type,pbc,link_cell_length,bypass_link_cells, &
                    use_verlet_list,read_xmol,rigid,isotropic
    namelist/thermal/pressure
    namelist/bookkeeping/file_output_int,traj_output_int,timer_qtime, &
                    timer_closetime,eq_mc_cycles,max_mc_cycles,eq_adjust_mc, &
                    mc_target_ratio, eq_adjust_int



    ! command line data
    integer(kind=it)              :: iarg,idata
    integer(kind=it)              :: num_args
    character(30), dimension(0:10):: command_line
    character(30)                 :: file_name,seedname
    integer(kind=it)              :: last_dot
    !    integer(kind=it),  external ::  iargc
    !    external  getarg

    integer(kind=it) :: ierr  ! error flag
    logical :: lexist = .false.

    ! check that there is only one argument.
    num_args = iargc()

    if (num_args<1) then

       ! No arguments, which might mean we are operating in library mode,
       ! or we want to use the default input file hs_alkane.input
       inquire(file='hs_alkane.input',exist=lexist)
       if (lexist) then
          write(*,'("Using default input filename - hs_alkane.input")')
          file_name = 'hs_alkane.input'
       else
          write(*,*)
          write(*,*) '                 H S _ A L K A N E          '
          write(*,*)
          write(*,*) '            Usage: hs_alkane <input file>   '
          write(*,*)
          write(*,*) '        D. Quigley - University of Warwick  '
          write(*,*)
          write(*,*) '   Alternatively (or if operating in library mode) '
          write(*,*) 'ensure hs_alkane.input exists in the current directory'
          write(*,*)
          stop
       end if
    else

       do iarg = 1, num_args
          call getarg(iarg,command_line(iarg))
       end do

       ! find the last dot in the filename.
       last_dot = len(seedname)
       do idata = 1, len(seedname)
          if(command_line(1)(idata:idata)=='.') last_dot = idata
       end do

       ! set the seedname.
       seedname = command_line(1)(1:last_dot-1)

       ! set the file name
       file_name = trim(command_line(1))

    end if

    ! open it
    open (25,file=file_name,status='old',iostat=ierr)
    if(ierr/=0) stop 'Unable to open input file.'

    read(25,nml=system,iostat=ierr)
    if(ierr/=0) stop 'Error reading system namelist'

    ! enforce exclusivity of link cells and verlet list
    if ( (.not.bypass_link_cells).and.(use_verlet_list) ) then
       write(0,'("Error in system namelist - must bypass link cells if using ")')
       write(0,'("Verlet neighbour list instead.")')
       stop
    end if


    max_regrow = nbeads - 1

    read(25,nml=thermal,iostat=ierr)
    if (ierr/=0) stop 'Error reading thermal namelist'
    pressure = pressure/sigma**3

    read(25,nml=bookkeeping,iostat=ierr)
    if (ierr/=0) stop 'Error reading bookkeeping namelist'

    close(25)

  end subroutine io_read_input


  subroutine c_wrap_io_read_xmol(new_basefile) bind(c, name='io_read_xmol')
    !-------------------------------------------------------------------------!
    ! Interface to io_read_xmol for use when hs_alkane is compiled as a       !
    ! library. Overides the default filename to read from.                    !
    !-------------------------------------------------------------------------!
    ! D.Quigley August 2021                                                   !
    !-------------------------------------------------------------------------!
    implicit none
    
    character(kind=c_char), intent(in) :: new_basefile(*)
    character(len=60) :: oldbase
    integer :: term, i

    oldbase = basefilename

    basefilename = ""
    
    term = 1
    do i = 1,60

       if (new_basefile(i) /= c_null_char)  then
          basefilename(i:i) = new_basefile(i)
       else
          exit
       end if

    end do
    
    call io_read_xmol()
    basefilename = oldbase

    
    return
    
  end subroutine c_wrap_io_read_xmol
  
  subroutine io_read_xmol()
    !-------------------------------------------------------------------------!
    ! Reads nchains of nbeads each into the array Rchain in alkane module.    !
    ! Expects file chain.xmol to exist in present working directory. If this  !
    ! calculation uses more than one box, then files chain.xmol.xx will be    !
    ! read, where xx = 01,02,03.......nboxes.                                 !
    !-------------------------------------------------------------------------!
    ! D.Quigley January 2011                                                  !
    !-------------------------------------------------------------------------!
    use alkane, only : Rchain,nchains,nbeads,chain_created, &
                       alkane_construct_linked_lists,L,alkane_update_linked_lists
    use box   , only : box_update_recipmatrix,pbc,hmatrix,recip_matrix,nboxes, &
                       box_construct_link_cells,box_minimum_image
    implicit none

    integer(kind=it) :: ierr,ichain,ibead,dumint,ibox
    character(2)     :: dumchar
    character(5)     :: boxstring
    character(60)    :: filename

    real(kind=dp),dimension(3)    :: tmpvect,old_pos,new_pos
    
    ! Optional argument used only when called from C/Python
    do ibox = 1,nboxes

       filename = basefilename
       !write(0,'("Using : ",A60)')filename
       write(boxstring,'(".",I4.4)')ibox

       if ( nboxes > 1 ) filename = trim(filename)//boxstring

       open(unit=25,file=trim(filename),status='old',iostat=ierr)
       if (ierr/=0) then
          write(0,'("Error : Could not open local input file : ",A60)')filename
          if (io_standalone) then
             stop
          else
             chain_created(:,:) = .false.
             return
          end if
       end if

       read(25,*)dumint
       if (dumint/=nchains*nbeads) stop 'Wrong number of beads in chain.xmol'

       read(25,*)hmatrix(:,:,ibox)
       if (pbc) then
          call box_update_recipmatrix(ibox)
       else
          recip_matrix = 0.0_dp
       end if

       do ichain = 1,nchains
          do ibead = 1,nbeads
             read(25,*)dumchar,Rchain(:,ibead,ichain,ibox)
          end do
          chain_created(ichain,ibox) = .true.
       end do

       close(25)

       ! Update data structures which depend on cell vectors and positions
       call box_update_recipmatrix(ibox)

       ! Construct link cell and linked list data structures
       call box_construct_link_cells(ibox)
       call alkane_construct_linked_lists(ibox)

       !--------------------------!
       ! Enforce correct imaging  !
       !--------------------------!

       ! Position all beads relative to the first, i.e. if any part of a chain
       ! has been wrapped around the PBCs but it back in the correct position
       ! relative to the first bead in the chain. This will result in some beads
       ! being outside the box, but that's not a problem. All distances computed
       ! will still be correct.
       if (nbeads>1) then
          do ichain = 1,nchains
          
             do ibead = 2,nbeads
                old_pos = Rchain(:,ibead,ichain,ibox)
                tmpvect = box_minimum_image(Rchain(:,ibead-1,ichain,ibox),old_pos,ibox)
!                tmpvect = tmpvect/sqrt(dot_product(tmpvect,tmpvect))
                new_pos = Rchain(:,ibead-1,ichain,ibox) + tmpvect !*L
                Rchain(:,ibead,ichain,ibox) = new_pos
                call alkane_update_linked_lists(ibead,ichain,ibox,old_pos,new_pos)
             end do
             
          end do                 
       end if
       
       
    end do ! end loop over boxes

  end subroutine io_read_xmol

  subroutine c_wrap_io_write_xmol(new_outfile) bind(c, name='io_write_xmol')
    !-------------------------------------------------------------------------!
    ! Interface to io_write_xmol for use when hs_alkane is compiled as a      !
    ! library. Overides the default filename to write to.                     !
    !-------------------------------------------------------------------------!
    ! D.Quigley December 2021                                                 !
    !-------------------------------------------------------------------------!
    implicit none
    
    character(kind=c_char), intent(in) :: new_outfile(*)
    character(len=60) :: oldbase
    integer :: term, i

    oldbase = outfilename

    outfilename = ""
    
    term = 1
    do i = 1,60

       if (new_outfile(i) /= c_null_char)  then
          outfilename(i:i) = new_outfile(i)
       else
          exit
       end if

    end do
    
    call io_write_xmol()
    outfilename = oldbase

    
    return
    
  end subroutine c_wrap_io_write_xmol
  
  subroutine io_write_xmol()
    !-------------------------------------------------------------------------!
    ! Writes nchains of nbeads from the array Rchain in alkane module to      !
    ! file.  If this calculation uses more than one box, then files           !
    ! final.xmol.xx will be written, where xx = 01,02,03..nboxes.             !
    !-------------------------------------------------------------------------!
    ! D.Quigley December 2021                                                 !
    !-------------------------------------------------------------------------!
    use alkane, only : Rchain,nchains,nbeads,chain_created, &
                       alkane_construct_linked_lists
    use box   , only : box_update_recipmatrix,pbc,hmatrix,recip_matrix,nboxes, &
                       box_construct_link_cells
    implicit none

    integer(kind=it) :: ierr,ichain,ibead,ibox
    character(5)     :: boxstring
    character(60)    :: denfile
    
    if (nboxes==1) then
       open(unit=25,file=trim(outfilename),status='replace',iostat=ierr)
       if (ierr/=0) stop 'Error opening final.xmol'
    else
       do ibox = 1,nboxes
          write(boxstring,'(".",I4.4)')ibox
          denfile = trim(outfilename)//boxstring
          open (unit=25+ibox-1,file=trim(denfile),status='replace',iostat=ierr)
          if (ierr/=0) stop 'Error opening xmol files for output coordinates'
       end do
    end if

    do ibox = 1,nboxes
       write(25+ibox-1,*)nbeads*nchains
       write(25+ibox-1,'(9F20.12)')hmatrix(:,:,ibox)
       do ichain = 1,nchains
          do ibead = 1,nbeads
             write(25+ibox-1,'("C ",3F20.12)')Rchain(:,ibead,ichain,ibox)
          end do
       end do
    end do
    
    do ibox = 1,nboxes
       close(25+ibox-1) ! close final.xmol
    end do
    
    return
    
  end subroutine io_write_xmol

end module io
