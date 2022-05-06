module lpt_module
  ! Intel
  use OMP_LIB


  ! Modules
  use constant_module
  use timing_module
  use cosmo_module, only : cosmo
  use input_module, only : input
  use sim_module  , only : sim,unit
  

  ! Default
  implicit none
  public


  ! Types
  type lpt_type
     ! Input
     integer(4)    :: Nproc,Nm1d
     character(10) :: make,order,assign
     character(80) :: dirin,dirout
     ! Variables
     integer(8)    :: Nmesh
     ! Arrays
     real(8), allocatable, dimension(:,:)     :: Pinit
     real(4), allocatable, dimension(:,:,:,:) :: gradphi
     ! Pointers
     real(8), pointer, dimension(:,:,:) :: delta1,delta2,gphi
  end type lpt_type


  ! Objects
  type(lpt_type) :: lpt


contains


  subroutine lpt_init
    ! Default
    implicit none


    ! Local variables
    integer(4) :: k


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from input
    lpt%make   = input%lpt_make
    lpt%order  = input%lpt_order
    lpt%assign = input%lpt_assign
    lpt%dirout = input%lpt_dir
    lpt%dirin  = input%sim_dirin
    lpt%Nproc  = input%sim_Nproc
    lpt%Nm1d   = input%mesh_Nm1d
    

    ! Number of cells/particles
    lpt%Nmesh = int(lpt%Nm1d,kind=8)**3


    ! Allocate arrays
    allocate(lpt%Pinit(3,lpt%Nm1d))
    allocate(lpt%gradphi(6,lpt%Nm1d,lpt%Nm1d,lpt%Nm1d))


    ! Init
    lpt%Pinit = 0


    ! First touch in parallel
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(k)
    !$omp do
    do k=1,lpt%Nm1d
       lpt%gradphi(:,:,:,k) = 0
    enddo
    !$omp end do
    !$omp end parallel


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : LPT init'
    return
  end subroutine lpt_init


!------------------------------------------------------------------------------!
! IO
!------------------------------------------------------------------------------!


  subroutine lpt_read
    ! Default
    implicit none


    ! Local variables
    integer(4)     :: un
    character(100) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Read Fourier transform of linear delta field
    un = 11
    fn = trim(lpt%dirin)//'delta_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn,form='binary')
    read(un) lpt%delta1
    close(un)

    
    ! Read gradphi array
    fn = trim(lpt%dirin)//'gradphi_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn,form='binary')
    read(un) lpt%gradphi
    close(un)

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : LPT read'
    return
  end subroutine lpt_read


  subroutine lpt_write
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: un
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Write Fourier transform of linear delta field
    un = 11
    fn = trim(lpt%dirout)//'delta_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Writing ',trim(fn)
    open(un,file=fn,form='binary')
    write(un) lpt%delta1
    close(un)

    
    ! Write gradphi array
    fn = trim(lpt%dirout)//'gradphi_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Writing ',trim(fn)
    open(un,file=fn,form='binary')
    write(un) lpt%gradphi
    close(un) 


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : LPT write'
    return
  end subroutine lpt_write


!------------------------------------------------------------------------------!
! Functions
!------------------------------------------------------------------------------!


  function x_lpt(i,j,k)
    ! Default
    implicit none


    ! Function arguments
    integer(4) :: i,j,k
    real(8)    :: x_lpt(3)


    ! Calculate position in grid units
    x_lpt = mod((/i,j,k/) - 0.5                   &
                - cosmo%D1*lpt%gradphi(1:3,i,j,k) &
                + cosmo%D2*lpt%gradphi(4:6,i,j,k) &
                + dble(lpt%Nm1d),dble(lpt%Nm1d))

    return
  end function x_lpt


  function v_lpt(i,j,k)
    ! Default
    implicit none


    ! Function arguments
    integer(4) :: i,j,k
    real(8)    :: v_lpt(3)


    ! Calculate velocity (dx/da) in grid units
    v_lpt = (-cosmo%v1*lpt%gradphi(1:3,i,j,k) &
             +cosmo%v2*lpt%gradphi(4:6,i,j,k))*unit%time

    return
  end function v_lpt


end module lpt_module
