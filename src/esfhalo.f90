module esfhalo_module
  ! Intel
  use IFPORT
  use OMP_LIB
  

  ! Modules
  use esf_module
  use constant_module
  use helper_module
  use mkl_module
  use timing_module
  use cosmo_module     , only : cosmo
  use cosmology_module , only : cosmo_calc
  use lpt_module       , only : lpt,x_lpt
  use mesh_module      , only : mesh
  use reion_module     , only : reion
  use sim_module       , only : sim,domain,unit
  use simulation_module, only : sim_calc,domain_set,domain_end


  ! Default
  implicit none
  public


contains


  subroutine esf_make
    ! Default
    implicit none


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Pointers
    esf%d1    => mesh%fft1
    esf%d2    => mesh%fft2
    esf%delta => mesh%fft1
    esf%fcol  => mesh%fft3

    
    ! Reionization midpoint
    cosmo%z = esf%zmid
    cosmo%a = 1/(1 + cosmo%z)
    call cosmo_calc
    call sim_calc


    ! Minimum halo mass
    esf%Rmin  = R_of_M(esf%Mmin)
    esf%Smmin = S_of_R(esf%Rmin)


    ! Smoothing scale set to resolution limit
    esf%Rsmth = R_of_M(unit%Msunh)
    esf%Msmth = M_of_R(esf%Rsmth)
    esf%Ssmth = S_of_R(esf%Rsmth)


    ! Collapsed overdensity threshold
    esf%delcz = 1.686/cosmo%D1


    ! Write to scrren
    write(*,*) 'ESF'
    write(*,*) 'M       : ',real((/esf%Mmin ,esf%Msmth/))
    write(*,*) 'R       : ',real((/esf%Rmin ,esf%Rsmth/))
    write(*,*) 'S       : ',real((/esf%Smmin,esf%Ssmth/))
    write(*,*) 'delta_c : ',real(esf%delcz)

    
    ! ESF
    if (esf%make == 'read') then
       ! IO
       call esf_read
    else
       ! Lagrangian version
       call esf_delta
       call esf_collapsedfraction
       call esf_halodensity

       ! IO
       if (esf%make == 'write') call esf_write
    endif
    

    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : ESF make'
    return
  end subroutine esf_make


  subroutine esf_delta
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k,n,ip
    integer(4) :: imax,jmax,kmax
    real(8)    :: Ak,kr,kx,ky,kz,wsmth


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! esf%delta and lpt%delta1 share mesh%fft1
    ! lpt%delta1 already Fourier transformed

    
    ! Init
    Ak   = 2*pi/esf%Nm1d
    imax = esf%Nm1d   + 2
    jmax = esf%Nm1d/2 + 1
    kmax = jmax


    ! Smooth
    !$omp parallel do              &
    !$omp default(shared)          &
    !$omp private(i,j,k,ip)        &
    !$omp private(kx,ky,kz,kr,wsmth)
    do k=1,esf%Nm1d
       if (k <= kmax) then
          kz = Ak*(k-1)
       else
          kz = Ak*(k-1-esf%Nm1d)
       endif

       do j=1,esf%Nm1d
          if (j <= jmax) then
             ky = Ak*(j-1)
          else
             ky = Ak*(j-1-esf%Nm1d)
          endif

          do i=1,esf%Nm1d+2,2
             ip = i + 1
             kx = Ak*((i-1)/2)

             if (i > 1 .or. j > 1 .or. k > 1) then
                ! kr in comoving h/Mpc
                kr = sqrt(kx**2 + ky**2 + kz**2)/unit%mesh_to_box

                ! Filter
                select case (esf%filter)
                case ('sharpk')
                   wsmth = sharpk_filter(kr*esf%Rsmth)
                case ('tophat')
                   wsmth = tophat_transform(kr*esf%Rsmth)
                case default
                   wsmth = 1
                end select

                ! Complex multiply
                esf%delta(i:ip,j,k) = lpt%delta1(i:ip,j,k)*wsmth
             else
                esf%delta(i:ip,j,k) = 0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Inverse FFT delta field
    call fft_3d(esf%delta,'b')


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : ESF delta'
    return
  end subroutine esf_delta


  subroutine esf_collapsedfraction
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k
    real(8)    :: delta,x
    real(8)    :: fcol,favg,fsig,fmin,fmax


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init
    favg = 0
    fsig = 0
    fmax = 0
    fmin = huge(fmin)


    ! Collapsed fraction
    ! See Bond etal 1991, Lacey & Cole 1993
    !$omp parallel do            &
    !$omp default(shared)        &
    !$omp private(i,j,k)         &
    !$omp private(fcol,x)        &
    !$omp reduction(+:favg,fsig) &
    !$omp reduction(max:fmax)    &
    !$omp reduction(min:fmin)
    do k=1,esf%Nm1d
       do j=1,esf%Nm1d
          do i=1,esf%Nm1d
             ! Collapsed fraction
             x    = (esf%delcz - esf%delta(i,j,k)) &
                  / sqrt(2*(esf%Smmin - esf%Ssmth))
             fcol = erfc(max(x,0D0))
             esf%fcol(i,j,k) = fcol

             ! Stats
             favg = favg + fcol
             fsig = fsig + fcol**2
             fmax = max(fmax,fcol)
             fmin = min(fmin,fcol)
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Write to screen
    favg = favg/esf%Nmesh
    fsig = sqrt(fsig/esf%Nmesh - favg**2)
    write(*,*) 'fcol : ',real((/favg,fsig,fmin,fmax/))


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : ESF collapsed fraction'
    return
  end subroutine esf_collapsedfraction


  subroutine esf_halodensity
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k,n
    

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(k)
    !$omp do
    do k=1,esf%Nm1d
       esf%d1(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,esf%Nm1d
       esf%d2(:,:,k) = 0
    enddo
    !$omp end do
    !$omp end parallel


    ! Construct density fields from particles
    ! Use cubical domain decomposition, see domain.f90

    select case (esf%assign)
    case ('ngp')
       !$omp parallel do       &
       !$omp default(shared)   & 
       !$omp private(i,n)      &
       !$omp schedule(dynamic,1)
       do i=1,domain%Nd
          call domain_set(n)
          call ngp_assignment(domain%i(:,:,n))
          call domain_end(n)
       enddo
       !$omp end parallel do
    case ('cic')
       !$omp parallel do       &
       !$omp default(shared)   & 
       !$omp private(i,n)      &
       !$omp schedule(dynamic,1)
       do i=1,domain%Nd
          call domain_set(n)
          call cic_assignment(domain%i(:,:,n))
          call domain_end(n)
       enddo
       !$omp end parallel do
    case ('tsc')
       !$omp parallel do       &
       !$omp default(shared)   & 
       !$omp private(i,n)      &
       !$omp schedule(dynamic,1)
       do i=1,domain%Nd
          call domain_set(n)
          call tsc_assignment(domain%i(:,:,n))
          call domain_end(n)
       enddo
       !$omp end parallel do
    end select


    ! Interlace and average density fields
    call density_interlace

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : ESF halo density field'
    return

    
  contains

    
    subroutine tsc_assignment(indx)
      ! Default
      implicit none


      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx


      ! Local variable
      integer(4) :: a,b,c,l,m,n
      real(8)    :: lm1d,mass,mcol
      real(8)    :: dx,dy,dz
      integer(4), dimension(-1:1) :: i,j,k
      real(8),    dimension(-1:1) :: wx,wy,wz
      real(8),    dimension(3)    :: x


      ! Init
      lm1d = dble(esf%Nm1d)
      mass = 1


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x    = mod(x_lpt(l,m,n) - 0.5 + lm1d,lm1d)
         mcol = mass*esf%fcol(l,m,n)

         ! Indices and weights
         i( 0)  = 1 + int(x(1))
         i(-1)  = 1 + mod(i(0) - 2 + esf%Nm1d,esf%Nm1d)
         i( 1)  = 1 + mod(i(0)               ,esf%Nm1d)
         j( 0)  = 1 + int(x(2))
         j(-1)  = 1 + mod(j(0) - 2 + esf%Nm1d,esf%Nm1d)
         j( 1)  = 1 + mod(j(0)               ,esf%Nm1d)
         k( 0)  = 1 + int(x(3))
         k(-1)  = 1 + mod(k(0) - 2 + esf%Nm1d,esf%Nm1d)
         k( 1)  = 1 + mod(k(0)               ,esf%Nm1d)
         dx     = x(1) - (i(0) - 0.5)
         wx( 0) = 0.75 - dx**2
         wx(-1) = 0.5*(1.5 - abs(dx + 1))**2
         wx( 1) = 0.5*(1.5 - abs(dx - 1))**2
         dy     = x(2) - (j(0) - 0.5)
         wy( 0) = 0.75 - dy**2
         wy(-1) = 0.5*(1.5 - abs(dy + 1))**2
         wy( 1) = 0.5*(1.5 - abs(dy - 1))**2
         dz     = x(3) - (k(0) - 0.5)
         wz( 0) = 0.75 - dz**2
         wz(-1) = 0.5*(1.5 - abs(dz + 1))**2
         wz( 1) = 0.5*(1.5 - abs(dz - 1))**2

         ! Add collapsed mass
         do c=-1,1
         do b=-1,1
         do a=-1,1
            esf%d1(i(a),j(b),k(c)) = esf%d1(i(a),j(b),k(c)) &
                                   + mcol*wx(a)*wy(b)*wz(c)
         enddo
         enddo
         enddo


         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices and weights
         i( 0)  = 1 + int(x(1))
         i(-1)  = 1 + mod(i(0) - 2 + esf%Nm1d,esf%Nm1d)
         i( 1)  = 1 + mod(i(0)               ,esf%Nm1d)
         j( 0)  = 1 + int(x(2))
         j(-1)  = 1 + mod(j(0) - 2 + esf%Nm1d,esf%Nm1d)
         j( 1)  = 1 + mod(j(0)               ,esf%Nm1d)
         k( 0)  = 1 + int(x(3))
         k(-1)  = 1 + mod(k(0) - 2 + esf%Nm1d,esf%Nm1d)
         k( 1)  = 1 + mod(k(0)               ,esf%Nm1d)
         dx     = x(1) - (i(0) - 0.5)
         wx( 0) = 0.75 - dx**2
         wx(-1) = 0.5*(1.5 - abs(dx + 1))**2
         wx( 1) = 0.5*(1.5 - abs(dx - 1))**2
         dy     = x(2) - (j(0) - 0.5)
         wy( 0) = 0.75 - dy**2
         wy(-1) = 0.5*(1.5 - abs(dy + 1))**2
         wy( 1) = 0.5*(1.5 - abs(dy - 1))**2
         dz     = x(3) - (k(0) - 0.5)
         wz( 0) = 0.75 - dz**2
         wz(-1) = 0.5*(1.5 - abs(dz + 1))**2
         wz( 1) = 0.5*(1.5 - abs(dz - 1))**2

         ! Add collapsed mass
         do c=-1,1
         do b=-1,1
         do a=-1,1
            esf%d2(i(a),j(b),k(c)) = esf%d2(i(a),j(b),k(c)) &
                                   + mcol*wx(a)*wy(b)*wz(c)
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo


      return
    end subroutine tsc_assignment


    subroutine cic_assignment(indx)
      ! Default
      implicit none


      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx


      ! Local variables
      integer(4) :: i,j,k,l,m,n
      integer(4) :: i1,i2,j1,j2,k1,k2
      real(8)    :: lm1d,mass,mcol
      real(8)    :: dx1,dx2,dy1,dy2,dz1,dz2
      real(8), dimension(3) :: x
      

      ! Init
      lm1d = dble(esf%Nm1d)
      mass = 1


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x    = mod(x_lpt(l,m,n) - 0.5 + lm1d,lm1d)
         mcol = mass*esf%fcol(l,m,n)

         ! Indices and weights
         i1  = 1 + int(x(1))
         i2  = 1 + mod(i1,esf%Nm1d)
         dx1 = i1 - x(1)
         dx2 = 1  - dx1
         j1  = 1 + int(x(2))
         j2  = 1 + mod(j1,esf%Nm1d)
         dy1 = j1 - x(2)
         dy2 = 1  - dy1
         k1  = 1 + int(x(3))
         k2  = 1 + mod(k1,esf%Nm1d)
         dz1 = k1 - x(3)
         dz2 = 1  - dz1

         ! Add mass
         esf%d1(i1,j1,k1) = esf%d1(i1,j1,k1) + mcol*dx1*dy1*dz1
         esf%d1(i2,j1,k1) = esf%d1(i2,j1,k1) + mcol*dx2*dy1*dz1
         esf%d1(i1,j2,k1) = esf%d1(i1,j2,k1) + mcol*dx1*dy2*dz1
         esf%d1(i2,j2,k1) = esf%d1(i2,j2,k1) + mcol*dx2*dy2*dz1
         esf%d1(i1,j1,k2) = esf%d1(i1,j1,k2) + mcol*dx1*dy1*dz2
         esf%d1(i2,j1,k2) = esf%d1(i2,j1,k2) + mcol*dx2*dy1*dz2
         esf%d1(i1,j2,k2) = esf%d1(i1,j2,k2) + mcol*dx1*dy2*dz2
         esf%d1(i2,j2,k2) = esf%d1(i2,j2,k2) + mcol*dx2*dy2*dz2


         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices and weights
         i1  = 1 + int(x(1))
         i2  = 1 + mod(i1,esf%Nm1d)
         dx1 = i1 - x(1)
         dx2 = 1  - dx1
         j1  = 1 + int(x(2))
         j2  = 1 + mod(j1,esf%Nm1d)
         dy1 = j1 - x(2)
         dy2 = 1  - dy1
         k1  = 1 + int(x(3))
         k2  = 1 + mod(k1,esf%Nm1d)
         dz1 = k1 - x(3)
         dz2 = 1  - dz1

         ! Add mass
         esf%d2(i1,j1,k1) = esf%d2(i1,j1,k1) + mcol*dx1*dy1*dz1
         esf%d2(i2,j1,k1) = esf%d2(i2,j1,k1) + mcol*dx2*dy1*dz1
         esf%d2(i1,j2,k1) = esf%d2(i1,j2,k1) + mcol*dx1*dy2*dz1
         esf%d2(i2,j2,k1) = esf%d2(i2,j2,k1) + mcol*dx2*dy2*dz1
         esf%d2(i1,j1,k2) = esf%d2(i1,j1,k2) + mcol*dx1*dy1*dz2
         esf%d2(i2,j1,k2) = esf%d2(i2,j1,k2) + mcol*dx2*dy1*dz2
         esf%d2(i1,j2,k2) = esf%d2(i1,j2,k2) + mcol*dx1*dy2*dz2
         esf%d2(i2,j2,k2) = esf%d2(i2,j2,k2) + mcol*dx2*dy2*dz2
      enddo
      enddo
      enddo


      return
    end subroutine cic_assignment


    subroutine ngp_assignment(indx)
      ! Default
      implicit none
 

      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx


      ! Local variables
      integer(4) :: i,j,k,l,m,n
      real(8)    :: lm1d,mass,mcol
      real(8), dimension(3) :: x

      
      ! Init
      lm1d = dble(esf%Nm1d)
      mass = 1


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x    = mod(x_lpt(l,m,n),lm1d)
         mcol = mass*esf%fcol(l,m,n)

         ! Indices
         i = 1 + int(x(1))
         j = 1 + int(x(2))
         k = 1 + int(x(3))

         ! Add mass
         esf%d1(i,j,k) = esf%d1(i,j,k) + mcol

         
         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices
         i = 1 + int(x(1))
         j = 1 + int(x(2))
         k = 1 + int(x(3))

         ! Add mass
         esf%d2(i,j,k) = esf%d2(i,j,k) + mcol
      enddo
      enddo
      enddo


      return
    end subroutine ngp_assignment

    
    subroutine density_interlace
      ! Default
      implicit none


      ! Local variables
      integer(4) :: i,j,k,ip
      integer(4) :: imax,jmax,kmax
      real(8)    :: Ak,kx,ky,kz
      real(8)    :: dneg,dsum,dscale
      complex(8) :: c,t,w
      real(8), dimension(2) :: d,davg,dsig,dmax,dmin
      

      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()


      ! Init
      Ak   = 2*pi/esf%Nm1d
      imax = esf%Nm1d   + 2
      jmax = esf%Nm1d/2 + 1
      kmax = jmax


      ! Forward FFT density fields
      call fft_3d(esf%d1,'f')
      call fft_3d(esf%d2,'f')


      ! Shift in Fourier space
      !$omp parallel do           &
      !$omp default(shared)       &
      !$omp private(i,j,k,ip)     &
      !$omp private(kx,ky,kz,c,t,w)
      do k=1,esf%Nm1d
         if (k <= kmax) then
            kz = Ak*(k-1)
         else
            kz = Ak*(k-1-esf%Nm1d)
         endif

         do j=1,esf%Nm1d
            if (j <= jmax) then
               ky = Ak*(j-1)
            else
               ky = Ak*(j-1-esf%Nm1d)
            endif

            do i=1,esf%Nm1d+2,2
               ip = i + 1
               kx = Ak*((i-1)/2)

               ! Skip k=0 mode
               if (i > 1 .or. j > 1 .or. k > 1) then
                  ! Shift
                  t = cmplx(0.,(kx + ky + kz)/2)
                  c = exp(t)*cmplx(esf%d2(i,j,k),esf%d2(ip,j,k))
                  esf%d2(i ,j,k) = real(c)
                  esf%d2(ip,j,k) = imag(c)

                  ! Deconvolution
                  w = assignment_transform(kx,ky,kz,esf%assign)
                  esf%d1(i:ip,j,k) = esf%d1(i:ip,j,k)/w
                  esf%d2(i:ip,j,k) = esf%d2(i:ip,j,k)/w

                  ! Average
                  esf%d2(i:ip,j,k) = (esf%d1(i:ip,j,k) + esf%d2(i:ip,j,k))/2
               else
                  esf%d2(i:ip,j,k) = (esf%d1(i:ip,j,k) + esf%d2(i:ip,j,k))/2

               endif
            enddo
         enddo
      enddo
      !$omp end parallel do


      ! Forward FFT density fields
      call fft_3d(esf%d1,'b')
      call fft_3d(esf%d2,'b')


      ! Adjust for negative densities
      if (.true.) then
         dneg = 0
         dsum = 0

         !$omp parallel do          &
         !$omp default(shared)      &
         !$omp private(i,j,k)       &
         !$omp reduction(+:dneg,dsum)
         do k=1,esf%Nm1d
            do j=1,esf%Nm1d
               do i=1,esf%Nm1d
                  dsum = dsum + esf%d1(i,j,k)
                  
                  if (esf%d1(i,j,k) < 0) then
                     dneg = dneg + esf%d1(i,j,k)
                  endif
               enddo
            enddo
         enddo
         !$omp end parallel do

         dscale = dsum/(dsum + abs(dneg))

         !$omp parallel do     &
         !$omp default(shared) &
         !$omp private(i,j,k)
         do k=1,esf%Nm1d
            do j=1,esf%Nm1d
               do i=1,esf%Nm1d
                  esf%d1(i,j,k) = max(esf%d1(i,j,k),0D0)*dscale
               enddo
            enddo
         enddo
         !$omp end parallel do
      endif
      

      ! Stats
      davg = 0
      dsig = 0
      dmax = 0
      dmin = huge(dmin)
      
      !$omp parallel do            &
      !$omp default(shared)        &
      !$omp private(i,j,k,d)       &
      !$omp reduction(+:davg,dsig) &
      !$omp reduction(max:dmax)    &
      !$omp reduction(min:dmin)
      do k=1,esf%Nm1d
         do j=1,esf%Nm1d
            do i=1,esf%Nm1d
               ! Save
               d(1) = esf%d1(i,j,k)
               d(2) = esf%d2(i,j,k)
               esf%rho1(i,j,k) = d(1)
               esf%rho2(i,j,k) = d(2)

               ! Stats           
               davg = davg + d
               dsig = dsig + d*d
               dmax = max(dmax,d)
               dmin = min(dmin,d)
            enddo
         enddo
      enddo
      !$omp end parallel do


      ! Write to screen
      davg = davg/esf%Nmesh
      dsig = sqrt(dsig/esf%Nmesh - davg**2)
      write(*,*) 'rho 1 : ',real((/davg(1),dsig(1),dmin(1),dmax(1)/))
      write(*,*) 'rho 2 : ',real((/davg(2),dsig(2),dmin(2),dmax(2)/))

      
      time2 = time()
      write(*,'(2a)') timing(time1,time2),' : ESF density interlace'
      return
    end subroutine density_interlace
  
  
  end subroutine esf_halodensity


end module esfhalo_module
