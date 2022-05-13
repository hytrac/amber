module h21cmreion_module
  ! Intel
  use OMP_LIB


  ! Modules
  use h21cm_module
  use constant_module
  use mkl_module
  use timing_module
  use cosmo_module     , only : cosmo,dcom_of_z
  use cosmology_module , only : cosmo_calc
  use mesh_module      , only : mesh
  use meshmake_module  , only : mesh_density,mesh_velocity
  use reion_module     , only : reion
  use sim_module       , only : sim,unit
  use simulation_module, only : sim_calc 


  ! Default
  implicit none
  public


contains


  subroutine h21cm_make
    ! Default
    implicit none


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Pointers
    h21cm%rhom => mesh%fft1
    h21cm%rhoh => mesh%fft2
    h21cm%Tb   => mesh%fft3


    ! Global 21cm
    call h21cm_global

    
    ! Reionization
    call h21cm_reion

    
    ! IO
    if (h21cm%make == 'write') call h21cm_write

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : 21cm make'
    return
  end subroutine h21cm_make


!------------------------------------------------------------------------------!
! Global brightness temperature
!------------------------------------------------------------------------------!

  
  subroutine h21cm_global
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: iz,Nz,un
    real(8)       :: xH,Tb,T0,z,dz
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Assuming Ts >> Tcmb
    ! Only valid for xH < 0.75 (e.g. Santos et al 2008)


    ! Init
    dz = 0.01
    Nz = nint((reion%zbeg - reion%zend)/dz)


    ! 21cm
    ! e.g. Madau et al (1997)
    T0 = 28*(cosmo%ob*cosmo%h**2/0.022)*sqrt(0.15/(cosmo%om*cosmo%h**2))


    ! Write
    un = 11
    fn = trim(h21cm%dirout)//'21cm_global.txt'
    write(*,*) 'Writing ',trim(fn)

    open(un,file=fn)
    write(un,'(a5,2a14)') 'z','x_H','T_b'
    
    do iz=1,Nz
       ! Redshift
       z = (iz - 1)*dz + reion%zend
       
       ! HI fraction
       xH = xH_of_z(z)

       ! Brightness temperature
       Tb = T0*xH*sqrt((1+z)/10)

       ! IO
       write(un,'(f5.2,2es14.6)') z,xH,Tb
    enddo

    close(un)
    
    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : 21cm global'
    return
  end subroutine h21cm_global


!------------------------------------------------------------------------------!
! Reionization
!------------------------------------------------------------------------------!

  
  subroutine h21cm_reion
    ! Default
    implicit none


    ! Local variables
    integer(4) :: iz


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Redshift ray
    call h21cm_ray


    ! Loop over decreasing redshift
    do iz=h21cm%Nz,1,-1
       ! Write to screen
       write(*,*) 'iz : ',iz
       write(*,*) 'z  : ',real(h21cm%ray(iz)%z)
       write(*,*) 'a  : ',real(h21cm%ray(iz)%a)
       write(*,*) 'r  : ',real(h21cm%ray(iz)%r)
       write(*,*) 'xH : ',real(h21cm%ray(iz)%xH)

       
       ! Midpoint
       h21cm%iz = iz
       cosmo%z  = h21cm%ray(iz)%z(2)
       cosmo%a  = h21cm%ray(iz)%a(2)
       call cosmo_calc
       call sim_calc

       
       ! Make density and velocity fields at midpoint
       ! See meshmake.f90    
       call mesh_density
       call mesh_velocity


       ! Simple 21cm power spectrum
       ! Assume Ts >> Tcmb
       ! Ignore peculiar vel, lightcone effects
       call h21cm_simple
    enddo


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : 21cm reion'
    return


  contains


    subroutine h21cm_ray
      ! Default
      implicit none


      ! Local variables
      integer(4) :: iz
      real(8)    :: z,z1,z2
      real(8)    :: a,a1,a2
      real(8)    :: r,r1,r2
      real(8)    :: x,x1,x2


      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()

      
      ! z spacing
      if (h21cm%zspacing == 'lin') then
         ! zdel = dz
         h21cm%Nz   = ceiling((h21cm%zmax - h21cm%zmin)/h21cm%zdel)
      else
         ! zdel = dlg(1+z)
         h21cm%zdel = log10(1 + h21cm%zdel/(1 + h21cm%zmin))
         h21cm%Nz   = ceiling((log10((1+h21cm%zmax)/(1+h21cm%zmin)))/h21cm%zdel)
      endif


      ! Allocate
      allocate(h21cm%ray(h21cm%Nz))


      ! Init
      z1 = h21cm%zmin
      a1 = 1/(1 + z1)
      r1 = dcom_of_z(0D0,0D0,z1) 
      x1 = xH_of_z(z1)


      ! Loop over redshifts
      do iz=1,h21cm%Nz
         ! Redshift
         if (h21cm%zspacing == 'lin') then
            z  = z1 + h21cm%zdel/2
            z2 = z  + h21cm%zdel/2
         else
            z  = (1 + z1)*10**(h21cm%zdel/2) - 1
            z2 = (1 + z )*10**(h21cm%zdel/2) - 1
         endif

         ! Scalefactor
         a  = 1/(1 + z )
         a2 = 1/(1 + z2) 
         
         ! Comoving distance
         r  = dcom_of_z(r1,z1,z )
         r2 = dcom_of_z(r ,z ,z2)

         ! HI fraction
         x  = xH_of_z(z )
         x2 = xH_of_z(z2)

         ! Save
         h21cm%ray(iz)%z  = (/z1,z,z2/)
         h21cm%ray(iz)%a  = (/a1,a,a2/)
         h21cm%ray(iz)%r  = (/r1,r,r2/)
         h21cm%ray(iz)%xH = (/x1,x,x2/)

         ! Next
         z1 = z2
         a1 = a2
         r1 = r2
         x1 = x2
      enddo


      time2 = time()
      write(*,'(2a)') timing(time1,time2),' : 21cm ray'
      return
    end subroutine h21cm_ray


    subroutine h21cm_simple
      ! Default
      implicit none


      ! Local variables
      integer(4)    :: i,j,k,l,ip,kk,un
      integer(4)    :: imax,jmax,kmax
      real(8)       :: Tb,Ak,kx,ky,kz,kr
      real(8)       :: pmm,phh,ptt,phm,w
      character(80) :: fn
      real(8), allocatable, dimension(:,:) :: pow
      
      
      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()


      ! Assume Ts >> Tcmb
      ! Ignore peculiar vel, lightcone effects


      ! Use interlaced and deconvolved fields
      ! mesh%rho2


      ! Allocate
      allocate(pow(7,h21cm%Nm1d))
      pow = 0


      ! 21cm
      ! e.g. Madau et al (1997)
      Tb = 28*(cosmo%ob*cosmo%h**2/0.022)*sqrt(0.15/(cosmo%om*cosmo%h**2)) &
         * sqrt((1+cosmo%z)/10)

      
      ! FFT
      Ak   = 2*pi/cosmo%Lbox
      imax = h21cm%Nm1d   + 2
      jmax = h21cm%Nm1d/2 + 1
      kmax = jmax


      ! HI density
      !$omp parallel do     &
      !$omp default(shared) &
      !$omp private(i,j,k)
      do k=1,h21cm%Nm1d
         do j=1,h21cm%Nm1d
            do i=1,h21cm%Nm1d
               ! Matter
               h21cm%rhom(i,j,k) = mesh%rho2(i,j,k)
               
               ! Neutral?
               if (cosmo%z > reion%zre(i,j,k)) then
                  h21cm%rhoh      =    mesh%rho2(i,j,k)
                  h21cm%Tb(i,j,k) = Tb*mesh%rho2(i,j,k)
               else
                  h21cm%rhoh      = 0
                  h21cm%Tb(i,j,k) = 0
               endif
            enddo
         enddo
      enddo

      call fft_3d(h21cm%rhom,'f')
      call fft_3d(h21cm%rhoh,'f')
      call fft_3d(h21cm%Tb  ,'f')

      !$omp parallel do                            & 
      !$omp default(shared)                        &
      !$omp private(i,j,k,ip,kk)                   &
      !$omp private(kx,ky,kz,kr,pmm,phh,ptt,phm,w) &
      !$omp reduction(+:pow)
      do k=1,h21cm%Nm1d
         if (k <= kmax) then
            kz = Ak*(k-1)
         else
            kz = Ak*(k-1-h21cm%Nm1d)
         endif

         do j=1,h21cm%Nm1d
            if (j <= jmax) then
               ky = Ak*(j-1)
            else
               ky = Ak*(j-1-h21cm%Nm1d)
            endif

            do i=1,imax,2
               ip = i + 1
               kx = Ak*((i-1)/2)
               kr = sqrt(kx**2 + ky**2 + kz**2)

               ! Bin Fourier modes
               if (i > 1 .or. j > 1 .or. k > 1) then
                  ! FFT weight
                  if (i == 1 .or. i == h21cm%Nm1d + 1) then
                     w = 1
                  else
                     w = 2
                  endif

                  ! Power
                  pmm = sum(h21cm%rhom(i:ip,j,k)/mesh%Nmesh &
                           *h21cm%rhom(i:ip,j,k)/mesh%Nmesh)
                  phh = sum(h21cm%rhoh(i:ip,j,k)/mesh%Nmesh &
                           *h21cm%rhoh(i:ip,j,k)/mesh%Nmesh)
                  ptt = sum(h21cm%Tb(  i:ip,j,k)/mesh%Nmesh &
                           *h21cm%Tb(  i:ip,j,k)/mesh%Nmesh)
                  phm = sum(h21cm%rhoh(i:ip,j,k)/mesh%Nmesh &
                           *h21cm%rhom(i:ip,j,k)/mesh%Nmesh)

                  ! Add to bin
                  kk        = nint(kr/Ak)
                  pow(1,kk) = pow(1,kk) + w
                  pow(2,kk) = pow(2,kk) + w*pmm
                  pow(3,kk) = pow(3,kk) + w*phh
                  pow(4,kk) = pow(4,kk) + w*ptt
                  pow(5,kk) = pow(5,kk) + w*phm
               endif
            enddo
         enddo
      enddo
      !$omp end parallel do


      ! Power spectra
      do k=1,h21cm%Nm1d
         ! Divide by weights
         if (pow(1,k) > 0) then
            pow(2:5,k) = pow(2:5,k)/pow(1,k)
         else
            pow(2:5,k) = 0
         endif

         ! Bias and cc
         if (pow(1,k) > 0) then
            pow(6,k) = sqrt(pow(3,k)/pow(2,k))
            pow(7,k) = pow(5,k)/sqrt(pow(2,k)*pow(3,k))
         endif
 
         ! Save
         h21cm%Pmm(k) = pow(2,k)*cosmo%Lbox**3
         h21cm%Phh(k) = pow(3,k)*cosmo%Lbox**3
         h21cm%PTT(k) = pow(4,k)*cosmo%Lbox**3
      enddo


      ! Write
      un = 11
      fn = trim(h21cm%dirout)//'power_'//trim(sim%zstr)//'.txt'
      write(*,*) 'Writing ',trim(fn)
      
      open(un,file=fn)
      write(un,'(6a14)') 'k','P_mm','P_hh','P_TT','b_hm','r_hm'
      
      do k=1,h21cm%Nm1d/2
         write(un,'(6es14.6)') h21cm%k(k),h21cm%Pmm(k),h21cm%Phh(k), &
              h21cm%PTT(k),pow(6:7,k)
      enddo
      
      close(un)
      
      
      time2 = time()
      write(*,'(2a)') timing(time1,time2),' : 21cm simple'
      return
    end subroutine h21cm_simple

    
  end subroutine h21cm_reion


end module h21cmreion_module
