module mesh_module
  ! Intel
  use OMP_LIB


  ! Modules
  use constant_module
  use timing_module
  use cosmo_module, only : cosmo
  use input_module, only : input
  use sim_module  , only : sim


  ! Default
  implicit none
  public


  ! Types
  type mesh_type
     ! Input
     integer(4)    :: Nproc,Nm1d
     character(10) :: make
     character(80) :: dirin,dirout
     ! Variables
     integer(8)    :: Nmesh
     ! Arrays
     integer(8), allocatable, dimension(:,:)     :: proc
     real(4)   , allocatable, dimension(:,:,:)   :: rho1,rho2
     real(4)   , allocatable, dimension(:,:,:,:) :: vel1,mom2
     real(8)   , allocatable, dimension(:,:,:)   :: fft1,fft2,fft3
     ! Pointers
     real(4), pointer, dimension(:,:,:,:) :: mom1
     real(8), pointer, dimension(:,:,:)   :: d1,d2,p1,p2
  end type mesh_type
  

  ! Objects
  type(mesh_type), target :: mesh
  

contains


  subroutine mesh_init
    ! Default
    implicit none


    ! Local variables
    integer(4) :: iproc
    integer(4) :: k
    integer(8) :: Nmesh_proc
    

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from input
    mesh%make   = input%mesh_make
    mesh%Nm1d   = input%mesh_Nm1d
    mesh%dirout = input%mesh_dir
    mesh%dirin  = input%sim_dirout
    mesh%Nproc  = input%sim_Nproc


    ! Number of mesh cells
    mesh%Nmesh = int(mesh%Nm1d,kind=8)**3


    ! Number of mesh cells per processor
    Nmesh_proc = mesh%Nmesh/mesh%Nproc + min(1,mod(mesh%Nmesh,mesh%Nproc))


    ! Allocate arrays
    allocate(mesh%rho1(  mesh%Nm1d,mesh%Nm1d,mesh%Nm1d))
    allocate(mesh%rho2(  mesh%Nm1d,mesh%Nm1d,mesh%Nm1d))
    allocate(mesh%vel1(3,mesh%Nm1d,mesh%Nm1d,mesh%Nm1d))
    allocate(mesh%mom2(3,mesh%Nm1d,mesh%Nm1d,mesh%Nm1d))
    allocate(mesh%proc(2,mesh%Nproc))


    ! Allocate FFT arrays with padding
    allocate(mesh%fft1(mesh%Nm1d+2,mesh%Nm1d,mesh%Nm1d))
    allocate(mesh%fft2(mesh%Nm1d+2,mesh%Nm1d,mesh%Nm1d))
    allocate(mesh%fft3(mesh%Nm1d+2,mesh%Nm1d,mesh%Nm1d))


    ! First touch in parallel
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(iproc,k)
    !$omp do
    do iproc=1,mesh%Nproc
       mesh%proc(1,iproc) = 1 + (iproc-1)*Nmesh_proc
       mesh%proc(2,iproc) = min(iproc*Nmesh_proc,mesh%Nmesh)
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%rho1(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%rho2(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%vel1(:,:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%mom2(:,:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%fft1(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%fft2(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%fft3(:,:,k) = 0
    enddo
    !$omp end do
    !$omp end parallel


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : MESH init'
    return
  end subroutine mesh_init


!------------------------------------------------------------------------------!
! IO
!------------------------------------------------------------------------------!


  subroutine mesh_read
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: i,j,k
    integer(4)    :: un
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Read density fields
    un = 11
    fn = trim(mesh%dirout)//'rho1_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn,form='binary')
    read(un) mesh%rho1
    close(un)

    fn = trim(mesh%dirout)//'rho2_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn,form='binary')
    read(un) mesh%rho2
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : MESH read'
    return
  end subroutine mesh_read
  

  subroutine mesh_write
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: i,j,k
    integer(4)    :: un
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Write density fields
    un = 11
    fn = trim(mesh%dirout)//'rho1_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Writing ',trim(fn)
    open(un,file=fn,form='binary')
    write(un) mesh%rho1
    close(un)

    fn = trim(mesh%dirout)//'rho2_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Writing ',trim(fn)
    open(un,file=fn,form='binary')
    write(un) mesh%rho2
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : MESH write'
    return
  end subroutine mesh_write


end module mesh_module
