module cmbreion_module
  ! Intel
  use OMP_LIB


  ! Modules
  use cmb_module
  use constant_module
  use mkl_module
  use timing_module
  use cosmo_module       , only : cosmo,dcom_of_z
  use cosmology_module   , only : cosmo_calc
  use mesh_module        , only : mesh
  use meshmake_module    , only : mesh_density,mesh_velocity
  use reion_module       , only : reion
  use reionization_module, only : reion_ionfrac
  use sim_module         , only : sim,unit
  use simulation_module  , only : sim_calc 


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
    cmb%ne => mesh%fft1
    cmb%qx => mesh%fft1
    cmb%qy => mesh%fft2
    cmb%qz => mesh%fft3


    ! CMB
    if (cmb%make == 'make' .or. cmb%make =='write') then
       ! Thomson optical depth
       call cmb_tau

       ! Reionization
       call cmb_reion

       ! IO
       if (cmb%make == 'write') call cmb_write
    endif

    
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

       
       ! Make density and velocity fields at midpoint
       ! See meshmake.f90    
       call mesh_density
       call mesh_velocity


       ! Patchy tau,ksz
       call cmb_powerspectrum
       call cmb_angularpowerspectrum
    enddo

    
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
      real(8)       :: xe
      real(8)       :: Ak,kx,kxh,ky,kyh,kz,kzh,kr,p,w
      character(80) :: fn
      real(8), dimension(2,3) :: q
      real(8), allocatable, dimension(:,:) :: powe,powq
      
      
      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()


      ! Use interlaced and deconvolved fields
      ! mesh%rho2, mesh%mom2


      ! Allocate
      allocate(powe(2,cmb%Nm1d))
      allocate(powq(2,cmb%Nm1d))
      powe = 0
      powq = 0


      ! Electron ionization fraction, xe = ne/ne_tot
      ! H ionized, He singly ionized
      xe = (1 + cosmo%YHe/(4*cosmo%XH))/(1 + cosmo%YHe/(2*cosmo%XH))


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
               ! Ionized?
               if (cosmo%z < reion%zre(i,j,k)) then
                  cmb%ne(i,j,k) = xe*mesh%rho2(i,j,k)
               else
                  cmb%ne(i,j,k) = 0
               endif
            enddo
         enddo
      enddo

      call fft_3d(cmb%ne,'f')

      !$omp parallel do              & 
      !$omp default(shared)          &
      !$omp private(i,j,k,kk)        &
      !$omp private(kx,ky,kz,kr,p,w) &
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
                  p = sum(cmb%ne(i:i+1,j,k)/mesh%Nmesh &
                         *cmb%ne(i:i+1,j,k)/mesh%Nmesh)

                  ! Add to bin
                  kk         = nint(kr/Ak)
                  powe(1,kk) = powe(1,kk) + w
                  powe(2,kk) = powe(2,kk) + w*p
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
      ! Pee(k) [(Mpc/h)^3]
      ! Pqq(k) [(Mpc/h)^3(km/s)^2]
      do k=1,cmb%Nm1d
         if (powe(1,k) > 0) then
            cmb%Pee(k) = powe(2,k)/powe(1,k)*cosmo%Lbox**3
            cmb%Pqq(k) = powq(2,k)/powq(1,k)*cosmo%Lbox**3
         else
            cmb%Pee(k) = 0
            cmb%Pqq(k) = 0
         endif
      enddo

      
      ! Write
      un = 11
      fn = trim(cmb%dirout)//'power_'//trim(sim%zstr)//'.txt'
      write(*,*) 'Writing ',trim(fn)
      
      open(un,file=fn)
      write(un,'(3a14)') 'k','P_ee','P_qq'
      
      do k=1,cmb%Nm1d/2
         write(un,'(3es14.6)') cmb%k(k),cmb%Pee(k),cmb%Pqq(k)
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


end module cmbreion_module
