module cmbreion_module
  ! Intel
  use OMP_LIB

  
  ! Healpix
  use alm_tools
  use pix_tools


  ! Modules
  use cmb_module
  use constant_module
  use mkl_module
  use timing_module
  use cosmo_module     , only : cosmo,dcom_of_z
  use cosmology_module , only : cosmo_calc
  use lpt_module       , only : x_lpt,v_lpt
  use mesh_module      , only : mesh
  use meshmake_module  , only : mesh_density,mesh_velocity
  use reion_module     , only : reion
  use sim_module       , only : sim,unit
  use simulation_module, only : sim_calc 


  ! Default
  implicit none
  public


contains


  subroutine cmb_make
    ! Default
    implicit none


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Pointers
    cmb%rho => mesh%fft1
    cmb%ne  => mesh%fft2
    cmb%qx  => mesh%fft1
    cmb%qy  => mesh%fft2
    cmb%qz  => mesh%fft3


    ! Thomson optical depth
    call cmb_tau

    
    ! Reionization
    call cmb_reion

    
    ! IO
    if (cmb%make    == 'write') call cmb_write
    if (cmb%mapmake == 'write') call cmb_mapwrite

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : CMB make'
    return
  end subroutine cmb_make


!------------------------------------------------------------------------------!
! Thomson optical depth
!------------------------------------------------------------------------------!

  
  subroutine cmb_tau
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: iz,Nz,un
    real(8)       :: tau,xe,z,z1,z2,dz
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init
    z   = 0
    xe  = xe_of_z(z)
    tau = 0
    dz  = 0.01
    Nz  = nint(max(cmb%zmax,reion%zbeg)/dz)


    ! Write
    un = 11
    fn = trim(cmb%dirout)//'tau.txt'
    write(*,*) 'Writing ',trim(fn)
    
    open(un,file=fn)
    write(un,'(a5,2a14)') 'z','x_e','tau'
    write(un,'(f5.2,2es14.6)') z,xe,tau
    
    do iz=1,Nz
       ! Redshift
       z  = (iz - 0.5)*dz
       z1 = z - dz/2
       z2 = z + dz/2
       
       ! Electron ionization fraction: xe = ne/ne_tot
       xe = xe_of_z(z2)

       ! Thomson optical depth
       tau = tau_of_z(tau,z1,z2)

       ! IO
       write(un,'(f5.2,2es14.6)') z2,xe,tau
    enddo
    
    close(un)


    ! Write to screen
    write(*,*) 'tau : ',real(tau)
    

    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : CMB tau'
    return
  end subroutine cmb_tau


!------------------------------------------------------------------------------!
! Reionization
!------------------------------------------------------------------------------!

  
  subroutine cmb_reion
    ! Default
    implicit none


    ! Local variables
    integer(4) :: iz


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Redshift ray 
    call cmb_ray
    

    ! Loop over decreasing redshift
    do iz=cmb%Nz,1,-1
       ! Write to screen
       write(*,*) 'iz  : ',iz
       write(*,*) 'z   : ',real(cmb%ray(iz)%z)
       write(*,*) 'a   : ',real(cmb%ray(iz)%a)
       write(*,*) 'r   : ',real(cmb%ray(iz)%r)
       write(*,*) 'xe  : ',real(cmb%ray(iz)%xe)
       write(*,*) 'tau : ',real(cmb%ray(iz)%tau)

       
       ! Midpoint
       cmb%iz  = iz
       cosmo%z = cmb%ray(iz)%z(2)
       cosmo%a = cmb%ray(iz)%a(2)
       call cosmo_calc
       call sim_calc

       
       ! Density and velocity fields at midpoint
       ! See meshmake.f90    
       call mesh_density
       call mesh_velocity


       ! Patchy tau,ksz
       call cmb_powerspectrum
       call cmb_angularpowerspectrum


       ! Healpix maps
       if (cmb%mapmake == 'make' .or. cmb%mapmake == 'write') then
          call cmb_mapmake
       endif
    enddo


    ! Healpix spectra
    if (cmb%mapmake == 'make' .or. cmb%mapmake == 'write') then
       call cmb_mapspectrum
    endif

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : CMB reion'
    return


  contains


    subroutine cmb_ray
      ! Default
      implicit none


      ! Local variables
      integer(4) :: iz
      real(8)    :: z,z1,z2
      real(8)    :: a,a1,a2
      real(8)    :: r,r1,r2
      real(8)    :: x,x1,x2
      real(8)    :: tau,tau1,tau2


      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()

      
      ! z spacing
      if (cmb%zspacing == 'lin') then
         ! zdel = dz
         cmb%Nz   = ceiling((cmb%zmax - cmb%zmin)/cmb%zdel)
      else
         ! zdel = dlg(1+z)
         cmb%zdel = log10(1 + cmb%zdel/(1 + cmb%zmin))
         cmb%Nz   = ceiling((log10((1+cmb%zmax)/(1+cmb%zmin)))/cmb%zdel)
      endif


      ! Allocate
      allocate(cmb%ray(cmb%Nz))


      ! Init
      z1   = cmb%zmin
      a1   = 1/(1 + z1)
      r1   = dcom_of_z(0D0,0D0,z1) 
      x1   = xe_of_z(z1)
      tau1 = tau_of_z(0D0,0D0,z1)


      ! Loop over redshifts
      do iz=1,cmb%Nz
         ! Redshift
         if (cmb%zspacing == 'lin') then
            z  = z1 + cmb%zdel/2
            z2 = z  + cmb%zdel/2
         else
            z  = (1 + z1)*10**(cmb%zdel/2) - 1
            z2 = (1 + z )*10**(cmb%zdel/2) - 1
         endif

         ! Scalefactor
         a  = 1/(1 + z )
         a2 = 1/(1 + z2) 
         
         ! Comoving distance
         r  = dcom_of_z(r1,z1,z )
         r2 = dcom_of_z(r ,z ,z2)

         ! Electron ionization fraction
         x  = xe_of_z(z )
         x2 = xe_of_z(z2)

         ! Thomson optical depth
         tau  = tau_of_z(tau1,z1,z )
         tau2 = tau_of_z(tau ,z ,z2)

         ! Save
         cmb%ray(iz)%z   = (/z1,z,z2/)
         cmb%ray(iz)%a   = (/a1,a,a2/)
         cmb%ray(iz)%r   = (/r1,r,r2/)
         cmb%ray(iz)%xe  = (/x1,x,x2/)
         cmb%ray(iz)%tau = (/tau1,tau,tau2/)

         ! Next
         z1   = z2
         a1   = a2
         r1   = r2
         x1   = x2
         tau1 = tau2
      enddo


      time2 = time()
      write(*,'(2a)') timing(time1,time2),' : CMB ray'
      return
    end subroutine cmb_ray


    subroutine cmb_powerspectrum
      ! Default
      implicit none


      ! Local variables
      integer(4)    :: i,j,k,l,ip,kk,un
      integer(4)    :: imax,jmax,kmax
      real(8)       :: Ak,kx,kxh,ky,kyh,kz,kzh,kr
      real(8)       :: p,pmm,pee,pem,w,xe
      character(80) :: fn
      real(8), dimension(2,3) :: q
      real(8), allocatable, dimension(:,:) :: powe,powq
      
      
      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()


      ! Use interlaced and deconvolved fields
      ! mesh%rho2, mesh%mom2


      ! Allocate
      allocate(powe(6,cmb%Nm1d))
      allocate(powq(2,cmb%Nm1d))
      powe = 0
      powq = 0


      ! Electron ionization fraction, xe = ne/ne_tot
      ! H ionized, He singly ionized
      xe = (cosmo%XH + cosmo%YHe/4)/(cosmo%XH + cosmo%YHe/2)


      ! FFT
      Ak   = 2*pi/cosmo%Lbox
      imax = cmb%Nm1d   + 2
      jmax = cmb%Nm1d/2 + 1
      kmax = jmax
      

      ! Ionized electron density
      !$omp parallel do     &
      !$omp default(shared) &
      !$omp private(i,j,k)
      do k=1,cmb%Nm1d
         do j=1,cmb%Nm1d
            do i=1,cmb%Nm1d
               ! Matter
               cmb%rho(i,j,k) = mesh%rho2(i,j,k)
               
               ! Ionized?
               if (cosmo%z < reion%zre(i,j,k)) then
                  cmb%ne(i,j,k) = xe*mesh%rho2(i,j,k)
               else
                  cmb%ne(i,j,k) = 0
               endif
            enddo
         enddo
      enddo

      call fft_3d(cmb%rho,'f')
      call fft_3d(cmb%ne ,'f')

      !$omp parallel do                        & 
      !$omp default(shared)                    &
      !$omp private(i,j,k,ip,kk)               &
      !$omp private(kx,ky,kz,kr,pmm,pee,pem,w) &
      !$omp reduction(+:powe)
      do k=1,cmb%Nm1d
         if (k <= kmax) then
            kz = Ak*(k-1)
         else
            kz = Ak*(k-1-cmb%Nm1d)
         endif

         do j=1,cmb%Nm1d
            if (j <= jmax) then
               ky = Ak*(j-1)
            else
               ky = Ak*(j-1-cmb%Nm1d)
            endif

            do i=1,imax,2
               ip = i + 1
               kx = Ak*((i-1)/2)
               kr = sqrt(kx**2 + ky**2 + kz**2)

               ! Bin Fourier modes
               if (i > 1 .or. j > 1 .or. k > 1) then
                  ! FFT weight
                  if (i == 1 .or. i == cmb%Nm1d + 1) then
                     w = 1
                  else
                     w = 2
                  endif

                  ! Power
                  pmm = sum(cmb%rho(i:ip,j,k)/mesh%Nmesh &
                           *cmb%rho(i:ip,j,k)/mesh%Nmesh)
                  pee = sum(cmb%ne( i:ip,j,k)/mesh%Nmesh  &
                           *cmb%ne( i:ip,j,k)/mesh%Nmesh)
                  pem = sum(cmb%ne( i:ip,j,k)/mesh%Nmesh  &
                           *cmb%rho(i:ip,j,k)/mesh%Nmesh)

                  ! Add to bin
                  kk         = nint(kr/Ak)
                  powe(1,kk) = powe(1,kk) + w
                  powe(2,kk) = powe(2,kk) + w*pmm
                  powe(3,kk) = powe(3,kk) + w*pee
                  powe(4,kk) = powe(4,kk) + w*pem
               endif
            enddo
         enddo
      enddo
      !$omp end parallel do


      ! Ionized electron momentum
      !$omp parallel do     &
      !$omp default(shared) &
      !$omp private(i,j,k)
      do k=1,cmb%Nm1d
         do j=1,cmb%Nm1d
            do i=1,cmb%Nm1d
               ! Ionized?
               if (cosmo%z < reion%zre(i,j,k)) then
                  cmb%qx(i,j,k) = xe*mesh%mom2(1,i,j,k)
                  cmb%qy(i,j,k) = xe*mesh%mom2(2,i,j,k)
                  cmb%qz(i,j,k) = xe*mesh%mom2(3,i,j,k)
               else
                  cmb%qx(i,j,k) = 0
                  cmb%qy(i,j,k) = 0
                  cmb%qz(i,j,k) = 0
               endif
            enddo
         enddo
      enddo
      !$omp end parallel do

      call fft_3d(cmb%qx,'f')
      call fft_3d(cmb%qy,'f')
      call fft_3d(cmb%qz,'f')

      !$omp parallel do                            & 
      !$omp default(shared)                        &
      !$omp private(i,j,k,ip,kk)                   &
      !$omp private(kx,kxh,ky,kyh,kz,kzh,kr,p,q,w) &
      !$omp reduction(+:powq)
      do k=1,cmb%Nm1d
         if (k <= kmax) then
            kz = Ak*(k-1)
         else
            kz = Ak*(k-1-cmb%Nm1d)
         endif

         do j=1,cmb%Nm1d
            if (j <= jmax) then
               ky = Ak*(j-1)
            else
               ky = Ak*(j-1-cmb%Nm1d)
            endif

            do i=1,imax,2
               ip = i + 1
               kx = Ak*((i-1)/2)
               kr = sqrt(kx**2 + ky**2 + kz**2)

               ! Bin Fourier modes
               if (i > 1 .or. j > 1 .or. k > 1) then
                  ! FFT weight
                  if (i == 1 .or. i == cmb%Nm1d + 1) then
                     w = 1
                  else
                     w = 2
                  endif

                  ! Projection
                  ! e.g. Park et al (2013)
                  kxh    = kx/kr
                  kyh    = ky/kr
                  kzh    = kz/kr
                  q(:,1) = cmb%qx(i:ip,j,k) - kxh*(cmb%qx(i:ip,j,k)*kxh &
                                                  +cmb%qy(i:ip,j,k)*kyh &
                                                  +cmb%qz(i:ip,j,k)*kzh)
                  q(:,2) = cmb%qy(i:ip,j,k) - kyh*(cmb%qx(i:ip,j,k)*kxh &
                                                  +cmb%qy(i:ip,j,k)*kyh &
                                                  +cmb%qz(i:ip,j,k)*kzh)
                  q(:,3) = cmb%qz(i:ip,j,k) - kzh*(cmb%qx(i:ip,j,k)*kxh &
                                                  +cmb%qy(i:ip,j,k)*kyh &
                                                  +cmb%qz(i:ip,j,k)*kzh)

                  ! Power
                  p = sum((q/cmb%Nmesh)**2)

                  ! Bin
                  kk         = nint(kr/Ak)
                  powq(1,kk) = powq(1,kk) + w
                  powq(2,kk) = powq(2,kk) + w*p
               endif
            enddo
         enddo
      enddo
      !$omp end parallel do


      ! Power spectra
      ! Pmm(k) [(Mpc/h)^3]
      ! Pee(k) [(Mpc/h)^3]
      ! Pqq(k) [(Mpc/h)^3(km/s)^2]
      do k=1,cmb%Nm1d
         ! Divide by weights
         if (powe(1,k) > 0) then
            powe(2:4,k) = powe(2:4,k)/powe(1,k)
         endif
         if (powq(1,k) > 0) then
            powq(2  ,k) = powq(2  ,k)/powq(1,k)
         endif

         ! Bias and cc
         if (powe(2,k) > 0 .and. powe(3,k) > 0) then
            powe(5,k)   = sqrt(powe(3,k)/powe(2,k))
            powe(6,k)   = powe(4,k)/sqrt(powe(2,k)*powe(3,k))
         else
            powe(5:6,k) = 0
         endif

         ! Save
         cmb%Pmm(k) = powe(2,k)*cosmo%Lbox**3
         cmb%Pee(k) = powe(3,k)*cosmo%Lbox**3
         cmb%Pqq(k) = powq(2,k)*cosmo%Lbox**3
      enddo

      
      ! Write
      un = 11
      fn = trim(cmb%dirout)//'power_'//trim(sim%zstr)//'.txt'
      write(*,*) 'Writing ',trim(fn)
      
      open(un,file=fn)
      write(un,'(6a14)') 'k','P_mm','P_ee','P_qq','b_em','r_em'
      
      do k=1,cmb%Nm1d/2
         write(un,'(6es14.6)') cmb%k(k),cmb%Pmm(k),cmb%Pee(k),cmb%Pqq(k), &
                               powe(5:6,k)
      enddo
      close(un)
      

      time2 = time()
      write(*,'(2a)') timing(time1,time2),' : CMB power spectrum'
      return
    end subroutine cmb_powerspectrum


    subroutine cmb_angularpowerspectrum
      ! Default
      implicit none


      ! Local variables
      real(8) :: a,r,dr,tau
      real(8) :: Atau,Aksz
      real(8), allocatable, dimension(:) :: k,p


      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()


      ! Allocate
      allocate(k(cmb%lmin:cmb%lmax))
      allocate(p(cmb%lmin:cmb%lmax))


      ! Ray
      a   = cmb%ray(cmb%iz)%a(2)
      r   = cmb%ray(cmb%iz)%r(2)
      dr  = cmb%ray(cmb%iz)%r(3) - cmb%ray(cmb%iz)%r(1)
      tau = cmb%ray(cmb%iz)%tau(2)


      ! Patchy tau
      ! e.g. Dvorkin & Smith (2009)
      
      ! Limber approximation Pee(k = l/r)
      k = cmb%l/r
      call spline_cubic(cmb%k,cmb%Pee,k,p)

      ! (Mpc/h)^2: (Mpc/h)^3 from Pee, (Mpc/h)^-1 from dr/r^2
      Atau = (sT_cgs*cosmo%ne0)**2*(Mpc2cm/cosmo%h)**2

      ! Integrate/sum
      cmb%Ctau = cmb%Ctau + Atau*p/(r**2*a**4)*dr      

      
      ! Patchy KSZ
      ! e.g. Park et al (2013)

      ! Limber approximation Pqq(k = l/r)
      k = cmb%l/r
      call spline_cubic(cmb%k,cmb%Pqq,k,p)

      ! (Mpc/h)^2: (Mpc/h)^3 from Pqq, (Mpc/h)^-1 from dr/r^2
      ! (1E5)^2  : (km/s)^2  from Pqq
      Aksz = (sT_cgs*cosmo%ne0/c_cgs)**2*(Mpc2cm/cosmo%h)**2*(1E5)**2

      ! Integrate/sum
      cmb%Cksz = cmb%Cksz + Aksz*(p/2)*exp(-2*tau)/(r**2*a**4)*dr


      time2 = time()
      write(*,'(2a)') timing(time1,time2),' : CMB angular power spectrum'
      return
    end subroutine cmb_angularpowerspectrum

    
  end subroutine cmb_reion


  subroutine cmb_mapmake
    ! Default
    implicit none


    ! Local variables
    integer(4) :: iproc,k,kmin,kmax
    integer(8) :: ip
    real(8)    :: rbuf,rmax,kscale
    real(8), dimension(2) :: p,pavg,psig,pmax,pmin
  

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Redshift shell
    ! Buffer in Mpc/h for LPT particle displacement
    rbuf = 10
    rmax = cmb%ray(cmb%iz)%r(3) + rbuf
    kmax = 1 + int(rmax*unit%box_to_mesh)
    kmin = -kmax + 1


    ! Maps
    ! Loop over k latitudes
    ! Skip to avoid thread collisions
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(k)
    !$omp do schedule(dynamic,1)
    do k=kmin,kmax,2
       call map_make(k)
    enddo
    !$omp end do
    !$omp do schedule(dynamic,1)
    do k=kmin+1,kmax,2
       call map_make(k)
    enddo
    !$omp end do    
    !$omp end parallel


    ! Stats
    pavg    = 0
    psig    = 0
    pmax    = 0
    pmin(1) = huge(0.)
    pmin(2) = 0
    
    !$omp parallel do            &
    !$omp default(shared)        &
    !$omp private(iproc,ip,p)    &
    !$omp reduction(+:pavg,psig) &
    !$omp reduction(max:pmax)    &
    !$omp reduction(min:pmin)
    do iproc=1,cmb%Nproc
       do ip=cmb%proc(1,iproc),cmb%proc(2,iproc)
          p(1) = cmb%tau(ip)
          p(2) = cmb%ksz(ip)
          pavg = pavg + p
          psig = psig + p**2
          pmax = max(pmax,p)
          pmin = min(pmin,p)
       enddo
    enddo
    !$omp end parallel do


    ! Write to screen
    pavg   = pavg/cmb%Npix
    psig   = sqrt(psig/cmb%Npix - pavg**2)
    kscale = cosmo%Tcmb0*1E6
    write(*,*) 'tau : ',real((/pavg(1),psig(1),pmin(1),pmax(1)/))
    write(*,*) 'ksz : ',real((/pavg(2),psig(2),pmin(2),pmax(2)/)*kscale)

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : CMB map make'
    return


  contains


    subroutine map_make(k1)
      ! Default
      implicit none


      ! Subroutine arguments
      integer(4) :: k1


      ! Local variables
      integer(4) :: i,j,k
      integer(4) :: i1,i2,j1,j2,k2
      integer(4) :: imin,imax,jmax
      integer(8) :: ipix
      real(8)    :: a,z,dcom,dcom1,dcom2,dang
      real(8)    :: r,phi,theta,omega_str,area
      real(8)    :: dtau,dksz,mue,tau,tau1,tau2,vlos
      real(8)    :: rmin,rmax,rminsq,rmaxsq,xmin,xmax,ymax,ysq,zsq
      real(8), dimension(3) :: x,x1,x2,v


      ! Redshift shell
      a     = cmb%ray(cmb%iz)%a(2)
      z     = cmb%ray(cmb%iz)%z(2)
      dcom1 = cmb%ray(cmb%iz)%r(1)
      dcom  = cmb%ray(cmb%iz)%r(2)
      dcom2 = cmb%ray(cmb%iz)%r(3)
      tau1  = cmb%ray(cmb%iz)%tau(1)
      tau2  = cmb%ray(cmb%iz)%tau(3)      
      dang  = a*dcom


      ! Pixel
      ! area in proper cm^2
      ! sig_str in Msolar/steradian
      omega_str = 4*pi/cmb%Npix
      area      = omega_str*(dang/cosmo%h*Mpc2cm)**2


      ! tau, ksz
      mue  = mH_cgs/(cosmo%XH + cosmo%YHe/4)
      dtau =  sT_cgs*unit%mgas/mue/area
      dksz = -sT_cgs/c_cgs*unit%vel*unit%mgas/mue/area


      ! Block limits in part/mesh units
      ! Buffer in Mpc/h for LPT particle displacement
      rmin   = (dcom1 - rbuf)*unit%box_to_mesh
      rmax   = (dcom2 + rbuf)*unit%box_to_mesh
      rminsq = rmin**2
      rmaxsq = rmax**2


      ! Loop over block
      k2    = 1 + mod(mod(k1,cmb%Nm1d) - 1 + cmb%Nm1d,cmb%Nm1d)
      x1(3) = k1 - 0.5
      x2(3) = floor(x1(3)/cmb%Nm1d)*cmb%Nm1d
      zsq   = x1(3)**2
      ymax  = sqrt(max(rmaxsq - zsq,0.))
      jmax  = 1 + int(ymax)
      
      do j1=-jmax+1,jmax
         j2    = 1 + mod(mod(j1,cmb%Nm1d) - 1 + cmb%Nm1d,cmb%Nm1d)
         x1(2) = j1 - 0.5
         x2(2) = floor(x1(2)/cmb%Nm1d)*cmb%Nm1d
         ysq   = x1(2)**2
         xmax  = sqrt(max(rmaxsq - (ysq + zsq),0.))
         xmin  = sqrt(max(rminsq - (ysq + zsq),0.))
         imax  = 1 + int(xmax)
         imin  = 1 + int(xmin)
         
         do i1=-imax+1,imax
            if (i1 <= -imin+1 .or. i1 >= imin) then
               i2    = 1 + mod(mod(i1,cmb%Nm1d) - 1 + cmb%Nm1d,cmb%Nm1d)
               x1(1) = i1 - 0.5
               x2(1) = floor(x1(1)/cmb%Nm1d)*cmb%Nm1d
               
               ! Particle in part/mesh units
               x = x_lpt(i2,j2,k2)
               v = v_lpt(i2,j2,k2)
               i = 1 + int(x(1))
               j = 1 + int(x(2))
               k = 1 + int(x(3))
      
               ! Lightcone coordinates in Mpc/h
               x = (x + x2)*unit%mesh_to_box
               r = sqrt(sum(x**2))

               if (r > dcom1 .and. r < dcom2) then
                  ! Pixel
                  x = x/r
                  call vec2ang(x,theta,phi)
                  call ang2pix_ring(cmb%Nside,theta,phi,ipix)

                  ! Ionized
                  if (z < reion%zre(i,j,k)) then
                     vlos = sum(v*x)
                     tau  = tau1 + (tau2-tau1)*(r-dcom1)/(dcom2-dcom1)
                     cmb%tau(ipix) = cmb%tau(ipix) + dtau
                     cmb%ksz(ipix) = cmb%ksz(ipix) + dksz*vlos*exp(-tau)
                  endif
               endif
            endif
         enddo
      enddo

      return
    end subroutine map_make


  end subroutine cmb_mapmake


  subroutine cmb_mapspectrum
    ! Default
    implicit none


    ! Local variables
    integer(4)     :: l,un
    character(100) :: fn
    real(8), dimension(2) :: zb
    real(8), allocatable, dimension(:,:) :: rw,cl_tau,cl_ksz


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Allocate
    allocate(cl_tau(0:cmb%Nlmax,1:1))
    allocate(cl_ksz(0:cmb%Nlmax,1:1))
    allocate(rw(    2*cmb%Nside,1:1))


    ! C_l
    rw = 1D0
    zb = (/-1D0,1D0/)

    call map2alm(cmb%Nside,cmb%Nlmax,cmb%Nmmax,cmb%tau,cmb%alm,zb,rw)
    call alm2cl( cmb%Nlmax,cmb%Nmmax,cmb%alm,cl_tau)

    call map2alm(cmb%Nside,cmb%Nlmax,cmb%Nmmax,cmb%ksz,cmb%alm,zb,rw)
    call alm2cl( cmb%Nlmax,cmb%Nmmax,cmb%alm,cl_ksz)


    ! Write tau
    un = 11
    fn = trim(cmb%dirout)//'cl_tau_healpix.txt'
    write(*,*) 'Writing ',trim(fn)
    
    open(11,file=fn)
    write(un,'(a6,a14)') 'l','C_l'

    do l=1,cmb%Nlmax
       write(un,'(i6,es14.6)') l,cl_tau(l,1)
    enddo
    close(un)


    ! Write ksz
    un = 11
    fn = trim(cmb%dirout)//'cl_ksz_healpix.txt'
    write(*,*) 'Writing ',trim(fn)

    open(11,file=fn)
    write(un,'(a6,a14)') 'l','C_l'

    do l=1,cmb%Nlmax
       write(un,'(i6,es14.6)') l,cl_ksz(l,1)
    enddo
    close(un)
    

    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : CMB map spectrum'
    return
  end subroutine cmb_mapspectrum


end module cmbreion_module
