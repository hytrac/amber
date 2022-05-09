module cmb_module
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
     real(8), dimension(3) :: z,a,r,xe,tau
  end type ray_type

  
  type cmb_type
     ! Input
     integer(4)    :: Nproc,Nm1d,lmin,lmax
     real(8)       :: zmin,zmax,zdel
     character(10) :: make,zspacing
     character(80) :: dirin,dirout
     ! Variables
     integer(4)    :: iz,Nz
     integer(8)    :: Nmesh
     ! Arrays
     integer(4)    , allocatable, dimension(:) :: l
     real(8)       , allocatable, dimension(:) :: k,Pee,Pqq,Ctau,Cksz
     type(ray_type), allocatable, dimension(:) :: ray
     ! Pointers
     real(8), pointer, dimension(:,:,:) :: ne,qx,qy,qz
  end type cmb_type


  ! Objects
  type(cmb_type) :: cmb


contains


  subroutine cmb_init
    ! Default
    implicit none


    ! Local variables
    integer(4) :: k,l


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from input
    cmb%make     = input%cmb_make
    cmb%zmin     = input%cmb_zmin
    cmb%zmax     = input%cmb_zmax
    cmb%zdel     = input%cmb_zdel
    cmb%zspacing = input%cmb_zspacing
    cmb%lmin     = input%cmb_lmin
    cmb%lmax     = input%cmb_lmax
    cmb%dirout   = input%cmb_dir
    cmb%dirin    = input%sim_dirin
    cmb%Nproc    = input%sim_Nproc
    cmb%Nm1d     = input%mesh_Nm1d


    ! Warnings
    if (cmb%zmax < reion%zbeg) then
       write(*,*) 'Warning: zmax < zbeg'
    endif
    if (cmb%zmin > reion%zend) then
       write(*,*) 'Warning: zmin > zend'
    endif


    ! Number of cells/particles
    cmb%Nmesh = int(cmb%Nm1d,kind=8)**3


    ! Allocate
    allocate(cmb%k(  cmb%Nm1d))
    allocate(cmb%Pee(cmb%Nm1d))
    allocate(cmb%Pqq(cmb%Nm1d))
    allocate(cmb%l(   cmb%lmin:cmb%lmax))
    allocate(cmb%Ctau(cmb%lmin:cmb%lmax))
    allocate(cmb%Cksz(cmb%lmin:cmb%lmax))


    ! Init
    cmb%Pee  = 0
    cmb%Pqq  = 0
    cmb%Ctau = 0
    cmb%Cksz = 0

    do k=1,cmb%Nm1d
       cmb%k(k) = 2*pi/cosmo%Lbox*k
    enddo
    
    do l=cmb%lmin,cmb%lmax
       cmb%l(l) = l
    enddo
    

    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : CMB init'
    return
  end subroutine cmb_init


!------------------------------------------------------------------------------!
! IO
!------------------------------------------------------------------------------!


  subroutine cmb_write
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: l,un
    character(80) :: fn
    
    
    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! tau
    un = 11
    fn = trim(cmb%dirout)//'Cl_tau.txt'
    write(*,*) 'Writing ',trim(fn)
    
    open(un,file=fn)
    write(un,'(a6,a14)') 'l','C_l'
    
    do l=cmb%lmin,cmb%lmax
       write(un,'(i6,es14.6)') cmb%l(l),cmb%Ctau(l)
    enddo
    close(un)


    ! KSZ
    un = 11
    fn = trim(cmb%dirout)//'Cl_ksz.txt'
    write(*,*) 'Writing ',trim(fn)
    
    open(un,file=fn)
    write(un,'(a6,a14)') 'l','C_l'
    
    do l=cmb%lmin,cmb%lmax
       write(un,'(i6,es14.6)') cmb%l(l),cmb%Cksz(l)
    enddo
    close(un)

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : CMB write'
    return
  end subroutine cmb_write


!------------------------------------------------------------------------------!
! Functions
!------------------------------------------------------------------------------!


  function xe_of_z(z)
    ! Default
    implicit none

    
    ! Function arguments
    real(8) :: z
    real(8) :: xe_of_z


    ! Electron ionization fraction: xe = ne/ne_tot
    if (z > reion%zend) then
       ! EoR
       xe_of_z = xi_of_z(z) &
               * (1 + cosmo%YHe/(4*cosmo%XH))/(1 + cosmo%YHe/(2*cosmo%XH))
    else if (z > 3) then
       ! H ionized, He singly ionized
       xe_of_z = (1 + cosmo%YHe/(4*cosmo%XH))/(1 + cosmo%YHe/(2*cosmo%XH))
    else
       ! H ionized, He doubly ionized
       xe_of_z =  1
    endif

    
    return
  end function xe_of_z
  
  
  function tau_of_z(tau1,z1,z2)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: tau1,z1,z2
    real(8) :: tau_of_z


    ! Local variables
    integer(8) :: i,N
    real(8)    :: a1,a2,da


    ! Init
    ! Set |da| about 1E-6
    ! If z2 > z1, then da < 0, and -dtau > 0
    a1 = 1/(1 + z1)
    a2 = 1/(1 + z2)
    N  = 2*(ceiling(abs(a2-a1)/1E-6/2,kind=8))
    da = (a2 - a1)/N
    

    ! Integrate using Simpson's rule
    tau_of_z = tau1 - (dtau_da(a1) + dtau_da(a2))/3*da
    
    do i=1,N-1,2
       tau_of_z = tau_of_z - dtau_da(a1 + i*da)*(4./3)*da
    enddo
    do i=2,N-1,2
       tau_of_z = tau_of_z - dtau_da(a1 + i*da)*(2./3)*da
    enddo
    

    return
    

  contains


    function dtau_da(a)
      ! Function arguments
      real(8) :: a,dtau_da

      ! Local variables
      real(8) :: H,ne0,z
      
      ! Cosmo
      z = 1/a - 1
      H = cosmo%H0 &
        * sqrt(cosmo%or/a**4 + cosmo%om/a**3 + cosmo%ol/a**(3*(1+cosmo%w)))

      ! Comoving electron number density in g/cm^3
      ne0 = xe_of_z(z)*cosmo%ne0

      ! tau integrand
      dtau_da = sT_cgs*ne0*(c_cgs/H)/a**4

      return
    end function dtau_da
    
    
  end function tau_of_z

  
end module cmb_module
