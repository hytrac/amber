module h21cm_module
  ! Intel
  use OMP_LIB


  ! Modules
  use constant_module
  use timing_module
  use cosmo_module, only : cosmo
  use input_module, only : input
  use reion_module, only : reion,xi_of_z
  use sim_module  , only : sim,unit


  ! Default
  implicit none
  public


  ! Types
  type ray_type
     real(8), dimension(3) :: z,a,r,xH
  end type ray_type

  
  type h21_type
     ! Input
     integer(4)    :: Nproc,Nm1d
     real(8)       :: zmin,zmax,zdel
     character(10) :: make,zspacing
     character(80) :: dirin,dirout
     ! Variables
     integer(4)    :: iz,Nz
     integer(8)    :: Nmesh
     ! Arrays
     real(8)       , allocatable, dimension(:) :: k,Phh
     type(ray_type), allocatable, dimension(:) :: ray
     ! Pointers
     real(8), pointer, dimension(:,:,:) :: Tb
  end type h21_type


  ! Objects
  type(h21_type) :: h21cm
  

contains


  subroutine h21cm_init
    ! Default
    implicit none


    ! Local variables
    integer(4) :: k


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from input
    h21cm%make     = input%h21cm_make
    h21cm%zmin     = input%h21cm_zmin
    h21cm%zmax     = input%h21cm_zmax
    h21cm%zdel     = input%h21cm_zdel
    h21cm%zspacing = input%h21cm_zspacing
    h21cm%dirout   = input%h21cm_dir
    h21cm%dirin    = input%sim_dirin
    h21cm%Nproc    = input%sim_Nproc
    h21cm%Nm1d     = input%mesh_Nm1d

    
    ! Number of cells/particles
    h21cm%Nmesh = int(h21cm%Nm1d,kind=8)**3


    ! Allocate
    allocate(h21cm%k(  h21cm%Nm1d))
    allocate(h21cm%Phh(h21cm%Nm1d))


    ! Init
    h21cm%Phh = 0

    do k=1,h21cm%Nm1d
       h21cm%k(k) = 2*pi/cosmo%Lbox*k
    enddo

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : 21cm init'
    return
  end subroutine h21cm_init


!------------------------------------------------------------------------------!
! IO
!------------------------------------------------------------------------------!  

  subroutine h21cm_write
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: i,un
    character(80) :: fn
    
    
    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : 21cm write'
    return
  end subroutine h21cm_write


!------------------------------------------------------------------------------!
! Functions
!------------------------------------------------------------------------------!


  function xH_of_z(z)
    ! Default
    implicit none

    ! Function arguments
    real(8) :: z
    real(8) :: xH_of_z

    ! Neutral fraction
    xH_of_z = 1 - xi_of_z(z)
    
    return
  end function xH_of_z
  

end module h21cm_module
