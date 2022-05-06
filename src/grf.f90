module grf_module
  ! Intel
  use OMP_LIB


  ! Modules
  use constant_module
  use timing_module
  use input_module, only : input
  use sim_module  , only : sim
  

  ! Default
  implicit none
  public


  ! Types
  type grf_type
     ! Input
     integer(4)    :: Nproc,Nm1d,seed
     character(10) :: make
     character(80) :: dirin,dirout
     ! Variables
     integer(8)    :: Nmesh
     ! Pointers
     real(8), pointer, dimension(:,:,:) :: grf
  end type grf_type


  ! Objects
  type(grf_type) :: grf


contains


  subroutine grf_init
    ! Default
    implicit none


    ! Local variables


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from input
    grf%make   = input%grf_make
    grf%seed   = input%grf_seed
    grf%dirout = input%grf_dir
    grf%dirin  = input%sim_dirin
    grf%Nproc  = input%sim_Nproc
    grf%Nm1d   = input%mesh_Nm1d
    

    ! Number of cells/particles
    grf%Nmesh = int(grf%Nm1d,kind=8)**3


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : GRF init'
    return
  end subroutine grf_init


!------------------------------------------------------------------------------!
! IO
!------------------------------------------------------------------------------!


  subroutine grf_read
    ! Default
    implicit none


    ! Local variables
    integer(4)     :: un
    character(100) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Read Fourier transform of GRF
    un = 11
    fn = trim(grf%dirin)//'grf_'//trim(sim%Nstr)//'.dat'
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn,form='binary')
    read(un) grf%grf
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : GRF read'
    return
  end subroutine grf_read


  subroutine grf_write
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: un
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Write Fourier transform of GRF
    un = 11
    fn = trim(grf%dirout)//'grf_'//trim(sim%Nstr)//'.dat'
    write(*,*) 'Writing ',trim(fn)
    open(un,file=fn,form='binary')
    write(un) grf%grf
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : GRF write'
    return
  end subroutine grf_write


end module grf_module
