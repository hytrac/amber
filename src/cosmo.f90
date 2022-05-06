module cosmo_module
  ! Intel
  use OMP_LIB


  ! Modules
  use constant_module
  use timing_module
  use input_module, only : input


  ! Default
  implicit none
  public


  ! Types
  type cosmo_type
     ! Input
     real(8)       :: Lbox
     real(8)       :: ob,om,ol,or
     real(8)       :: h,s8,ns,w
     real(8)       :: Tcmb0,XH,YHe
     character(80) :: file
     ! Variables
     real(8) :: a,z,t
     real(8) :: H0,Hz,fb,fc
     ! Densities
     real(8) :: rhocrit,rhocrit0
     real(8) :: rhob,rhob0,rhodm,rhodm0,rhom,rhom0
     real(8) :: nH,nH0,nHe,nHe0,ne,ne0
     ! Halo
     real(8) :: rho200a,rho200c,rho500c,rhovir     
     ! Growth factors
     real(8) :: D1,D2,delta0
     real(8) :: v1,v2
     ! Linear power spectrum
     integer(4) :: Nk
     real(8), allocatable, dimension(:,:) :: Plin
     ! IO
     character(80) :: dir
  end type cosmo_type


  ! Objects
  type(cosmo_type) :: cosmo


contains


  subroutine cosmo_init
    ! Default
    implicit none


    ! Local variables
    real(8), dimension(2) :: gf

    
    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from Input
    cosmo%Lbox   = input%cosmo_Lbox
    cosmo%om     = input%cosmo_om
    cosmo%ol     = input%cosmo_ol
    cosmo%ob     = input%cosmo_ob
    cosmo%or     = input%cosmo_or
    cosmo%h      = input%cosmo_h
    cosmo%s8     = input%cosmo_s8
    cosmo%ns     = input%cosmo_ns
    cosmo%w      = input%cosmo_w
    cosmo%Tcmb0  = input%cosmo_Tcmb0
    cosmo%XH     = input%cosmo_XH
    cosmo%YHe    = input%cosmo_YHe
    cosmo%file   = input%cosmo_file
    cosmo%dir    = input%cosmo_dir


    ! Redshift, scalefactor
    cosmo%a = 1
    cosmo%z = 0


    ! Growth factors
    gf = growthfactors_of_z(0D0)
    cosmo%delta0 = gf(1)

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : COSMO init'
    return
  end subroutine cosmo_init


!------------------------------------------------------------------------------!
! Functions
!------------------------------------------------------------------------------!

  
  function growthfactors_of_z(z)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: z
    real(8), dimension(2) :: growthfactors_of_z


    ! Local variables
    integer(8) :: i,N
    real(8)    :: a,a1,a2,aeq,da
    real(8), dimension(2) :: delta,dk,k1,k2,k3,k4


    ! Cosmo
    aeq      = cosmo%or/cosmo%om
    a1       = max(aeq,1E-6)
    a2       = 1/(1 + z)
    delta(1) = a1 + 2./3*aeq
    delta(2) = 1.


    ! Init
    ! Set da about 1E-4
    N  = ceiling(abs(a2 - a1)/1E-4,kind=8) 
    da = (a2 - a1)/N


    ! 4th-order Runge Kutta
    a = a1
    
    do i=1,N
       dk    = delta
       k1    = ode(a       ,dk)*da
       dk    = delta + k1/2
       k2    = ode(a + da/2,dk)*da
       dk    = delta + k2/2
       k3    = ode(a + da/2,dk)*da
       dk    = delta + k3
       k4    = ode(a + da  ,dk)*da
       delta = delta + (k1 + 2*k2 + 2*k3 + k4)/6
       a     = a + da
    enddo 


    ! delta and f = dlndelta/dlna
    growthfactors_of_z(1) = delta(1)
    growthfactors_of_z(2) = delta(2)*a/delta(1)


    return


  contains


    function ode(a,delta)
      ! Function arguments
      real(8) :: a
      real(8), dimension(2) :: delta,ode

      ! Local variables
      real(8) :: aw,H,q

      ! RHS of two 1st-order ODE from 2nd-order ODE
      aw     = a**(3*(1 + cosmo%w))
      H      = sqrt(cosmo%or/a**4 + cosmo%om/a**3 + cosmo%ol/aw)
      q      = (cosmo%or/a**4 + 0.5*cosmo%om/a**3 - cosmo%ol/aw)/H**2
      ode(1) = delta(2)
      ode(2) = 1.5*cosmo%om/a**5/H**2*delta(1) - (2-q)/a*delta(2)
      
      return
    end function ode

    
  end function growthfactors_of_z


  function dcom_of_z(d1,z1,z2)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: d1,z1,z2
    real(8) :: dcom_of_z


    ! Local variables
    integer(8) :: i,N
    real(8)    :: a1,a2,da


    ! Init
    ! Set |da| about 1E-6
    ! If z2 > z1, then da < 0 and -dr > 0
    a1 = 1/(1 + z1)
    a2 = 1/(1 + z2)
    N  = 2*(ceiling(abs(a2-a1)/1E-6/2,kind=8))
    da = (a2 - a1)/N

 
    ! Integrate using Simpson's rule
    dcom_of_z = d1 - (dcom_da(a1) + dcom_da(a2))/3*da
    
    do i=1,N-1,2
       dcom_of_z = dcom_of_z - dcom_da(a1 + i*da)*(4./3)*da
    enddo
    do i=2,N-1,2
       dcom_of_z = dcom_of_z - dcom_da(a1 + i*da)*(2./3)*da
    enddo

 
    return


  contains

    
    function dcom_da(a)
      ! Function arguments
      real(8) :: a,dcom_da
      
      ! Local variables
      real(8) :: H

      ! Cosmo
      H = 100*sqrt(cosmo%or/a**4 + cosmo%om/a**3 + cosmo%ol/a**(3*(1+cosmo%w)))

      ! dcom integrand
      dcom_da = (c_cgs/1E5)/H/a**2

      return
    end function dcom_da

    
  end function dcom_of_z


  function z_of_dcom(z1,d1,d2)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: z1,d1,d2
    real(8) :: z_of_dcom


    ! Local variables
    integer(8) :: i,N
    real(8)    :: r,dr,k1,k2,k3,k4


    ! Init
    ! Set |dr| about 1E-3 Mpc/h
    ! If d2 > d1, then dr > 0 and dz > 0
    N  = ceiling(abs(d2 - d1)/1E-3,kind=8)
    dr = (d2 - d1)/N

    
    ! 4th-order Runge Kutta
    z_of_dcom = z1
    r         = d1
    
    do i=1,N
       k1        = dz_dr(z_of_dcom       )*dr
       k2        = dz_dr(z_of_dcom + k1/2)*dr
       k3        = dz_dr(z_of_dcom + k2/2)*dr
       k4        = dz_dr(z_of_dcom + k3  )*dr
       z_of_dcom = z_of_dcom + (k1 + 2*k2 + 2*k3 + k4)/6
       r         = r + dr
    enddo 


    return


  contains


    function dz_dr(z)
      ! Function arguments
      real(8) :: z,dz_dr

      ! Local variables
      real(8) :: a,H

      ! Cosmo
      a = 1/(1 + z)
      H = 100*sqrt(cosmo%or/a**4 + cosmo%om/a**3 + cosmo%ol/a**(3*(1+cosmo%w)))

      ! ODE
      dz_dr = H/(c_cgs/1E5)

      return
    end function dz_dr
      
    
  end function z_of_dcom


end module cosmo_module
