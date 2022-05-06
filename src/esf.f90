module esf_module
  ! Intel
  use OMP_LIB


  ! Modules
  use constant_module
  use helper_module
  use timing_module
  use cosmo_module, only : cosmo
  use input_module, only : input
  use sim_module  , only : sim


  ! Default
  implicit none
  public


  ! Types
  type esf_type
     ! Input
     integer(4)    :: Nproc,Nm1d
     character(10) :: make,filter,assign
     character(80) :: dirin,dirout
     ! Variables
     integer(8) :: Nmesh
     real(8)    :: Mmin,Rmin,Smmin
     real(8)    :: Msmth,Rsmth,Ssmth
     real(8)    :: delcz,zmid
     ! Arrays
     real(4), allocatable, dimension(:,:,:) :: rho1,rho2
     ! Pointers
     real(8), pointer, dimension(:,:,:) :: delta,fcol,d1,d2
  end type esf_type


  ! Objects
  type(esf_type) :: esf


contains


  subroutine esf_init
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k
    

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from input
    esf%make   = input%esf_make
    esf%filter = input%esf_filter
    esf%assign = input%esf_assign 
    esf%dirout = input%esf_dir
    esf%dirin  = input%sim_dirin
    esf%Nproc  = input%sim_Nproc
    esf%Nm1d   = input%mesh_Nm1d
    

    ! Number of mesh cells
    esf%Nmesh = int(esf%Nm1d,kind=8)**3


    ! Minimum halo mass
    esf%Mmin = input%reion_Mmin


    ! Redshift midpoint
    esf%zmid = input%reion_zmid


    ! Allocate arrays
    allocate(esf%rho1(esf%Nm1d,esf%Nm1d,esf%Nm1d))
    allocate(esf%rho2(esf%Nm1d,esf%Nm1d,esf%Nm1d))


    ! First touch in parallel
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(k)
    !$omp do
    do k=1,esf%Nm1d
       esf%rho1(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,esf%Nm1d
       esf%rho2(:,:,k) = 0
    enddo
    !$omp end do
    !$omp end parallel
    

    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : ESF init'
    return
  end subroutine esf_init


!------------------------------------------------------------------------------!
! IO
!------------------------------------------------------------------------------!


  subroutine esf_read
    ! Default
    implicit none


    ! Local variables
    integer(4)     :: un
    character(100) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Read halo mass density field
    un = 11
    fn = trim(esf%dirin)//'rhoh1_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn,form='binary')
    read(un) esf%rho1
    close(un)

    fn = trim(esf%dirin)//'rhoh2_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn,form='binary')
    read(un) esf%rho2
    close(un)
    
    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : ESF read'
    return
  end subroutine esf_read


  subroutine esf_write
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: un
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Write halo mass density field
    un = 11
    fn = trim(esf%dirout)//'rhoh1_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Writing ',trim(fn)
    open(un,file=fn,form='binary')
    write(un) esf%rho1
    close(un)

    fn = trim(esf%dirout)//'rhoh2_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Writing ',trim(fn)
    open(un,file=fn,form='binary')
    write(un) esf%rho2
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : ESF write'
    return
  end subroutine esf_write
  

!------------------------------------------------------------------------------!
! Functions
!------------------------------------------------------------------------------!

  
  function M_of_T(T)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: T,z
    real(8) :: M_of_T


    ! Local variables
    real(8) :: mass,mu,rho


    ! M in Msun/h
    rho    = (18*pi**2)*cosmo%om*rhoc_cgs*cosmo%h**2*(1+cosmo%z)**3
    mu     = mH_cgs/(2*cosmo%XH + 3./4*cosmo%YHe)
    mass   = (2*k_cgs*T/(G_cgs*mu))**1.5/sqrt(4*pi/3*rho)
    M_of_T = mass/Msun_cgs*cosmo%h


    return
  end function M_of_T


  function M_of_R(R)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: R
    real(8) :: M_of_R


    ! Local variables
    real(8) :: M,rho0


    ! R in comoving Mpc/h
 

    ! Mean comoving matter density in units of (Msun/h)/(Mpc/h)^3
    rho0 = cosmo%om*rhoc_ast


    ! M in Msun/h
    M_of_R = 4*pi/3*rho0*R**3


    return
  end function M_of_R


  function R_of_M(M)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: M
    real(8) :: R_of_M


    ! Local variables
    real(8) :: R,rho0


    ! M in Msun/h
 

    ! Mean comoving matter density in units of (Msun/h)/(Mpc/h)^3
    rho0 = cosmo%om*rhoc_ast


    ! R in comoving Mpc/h
    R_of_M = (M/(4*pi/3*rho0))**(1./3)


    return
  end function R_of_M

  
  function S_of_R(R)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: R
    real(8) :: S_of_R


    ! Local variables
    integer(4) :: k
    real(8)    :: kr,dk,pow


    ! R in comoving Mpc/h


    ! Variance at z = 0
    S_of_R = 0

    do k=1,cosmo%Nk
       kr     = cosmo%Plin(1,k)
       dk     = cosmo%Plin(2,k)
       pow    = cosmo%Plin(3,k)
       S_of_R = S_of_R + pow*tophat_transform(kr*R)**2*(dk/kr)
    enddo


    return
  end function S_of_R

  
end module esf_module
