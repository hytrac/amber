module meshmake_module
  ! Intel
  use OMP_LIB


  ! Modules
  use mesh_module
  use constant_module
  use helper_module
  use mkl_module
  use timing_module
  use cosmo_module     , only : cosmo
  use lpt_module       , only : lpt,x_lpt,v_lpt
  use reion_module     , only : reion
  use sim_module       , only : sim,domain,unit
  use simulation_module, only : domain_set,domain_end


  ! Default
  implicit none
  public


contains


  subroutine mesh_density
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k,n


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Pointers
    mesh%d1 => mesh%fft1
    mesh%d2 => mesh%fft2

    
    ! Init
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(k)
    !$omp do
    do k=1,mesh%Nm1d
       mesh%d1(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%d2(:,:,k) = 0
    enddo
    !$omp end do
    !$omp end parallel


    ! Construct density fields from particles
    ! Use cubical domain decomposition, see domain.f90

    select case (lpt%assign)
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
    write(*,'(2a)') timing(time1,time2),' : MESH density field'
    return


  contains


    subroutine tsc_assignment(indx)
      ! Default
      implicit none


      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx


      ! Local variable
      integer(4) :: a,b,c,l,m,n
      real(8)    :: lm1d,mass
      real(8)    :: dx,dy,dz
      integer(4), dimension(-1:1) :: i,j,k
      real(8),    dimension(-1:1) :: wx,wy,wz
      real(8),    dimension(3)    :: x
      

      ! Init
      lm1d = dble(mesh%Nm1d)
      mass = 1.0


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x = mod(x_lpt(l,m,n) - 0.5 + lm1d,lm1d)

         ! Indices and weights
         i( 0)  = 1 + int(x(1))
         i(-1)  = 1 + mod(i(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         i( 1)  = 1 + mod(i(0)                ,mesh%Nm1d)
         j( 0)  = 1 + int(x(2))
         j(-1)  = 1 + mod(j(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         j( 1)  = 1 + mod(j(0)                ,mesh%Nm1d)
         k( 0)  = 1 + int(x(3))
         k(-1)  = 1 + mod(k(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         k( 1)  = 1 + mod(k(0)                ,mesh%Nm1d)
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

         ! Add mass
         do c=-1,1
         do b=-1,1
         do a=-1,1
            mesh%d1(i(a),j(b),k(c)) = mesh%d1(i(a),j(b),k(c)) &
                                      + mass*wx(a)*wy(b)*wz(c)
         enddo
         enddo
         enddo


         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices and weights
         i( 0)  = 1 + int(x(1))
         i(-1)  = 1 + mod(i(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         i( 1)  = 1 + mod(i(0)                ,mesh%Nm1d)
         j( 0)  = 1 + int(x(2))
         j(-1)  = 1 + mod(j(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         j( 1)  = 1 + mod(j(0)                ,mesh%Nm1d)
         k( 0)  = 1 + int(x(3))
         k(-1)  = 1 + mod(k(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         k( 1)  = 1 + mod(k(0)                ,mesh%Nm1d)
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

         ! Add mass
         do c=-1,1
         do b=-1,1
         do a=-1,1
            mesh%d2(i(a),j(b),k(c)) = mesh%d2(i(a),j(b),k(c)) &
                                      + mass*wx(a)*wy(b)*wz(c)
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
      real(8)    :: lm1d,mass
      real(8)    :: dx1,dx2,dy1,dy2,dz1,dz2
      real(8), dimension(3) :: x
      

      ! Init
      lm1d = dble(mesh%Nm1d)
      mass = 1


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x = mod(x_lpt(l,m,n) - 0.5 + lm1d,lm1d)

         ! Indices and weights
         i1  = 1 + int(x(1))
         i2  = 1 + mod(i1,mesh%Nm1d)
         dx1 = i1 - x(1)
         dx2 = 1  - dx1
         j1  = 1 + int(x(2))
         j2  = 1 + mod(j1,mesh%Nm1d)
         dy1 = j1 - x(2)
         dy2 = 1  - dy1
         k1  = 1 + int(x(3))
         k2  = 1 + mod(k1,mesh%Nm1d)
         dz1 = k1 - x(3)
         dz2 = 1  - dz1

         ! Add mass
         mesh%d1(i1,j1,k1) = mesh%d1(i1,j1,k1) + mass*dx1*dy1*dz1
         mesh%d1(i2,j1,k1) = mesh%d1(i2,j1,k1) + mass*dx2*dy1*dz1
         mesh%d1(i1,j2,k1) = mesh%d1(i1,j2,k1) + mass*dx1*dy2*dz1
         mesh%d1(i2,j2,k1) = mesh%d1(i2,j2,k1) + mass*dx2*dy2*dz1
         mesh%d1(i1,j1,k2) = mesh%d1(i1,j1,k2) + mass*dx1*dy1*dz2
         mesh%d1(i2,j1,k2) = mesh%d1(i2,j1,k2) + mass*dx2*dy1*dz2
         mesh%d1(i1,j2,k2) = mesh%d1(i1,j2,k2) + mass*dx1*dy2*dz2
         mesh%d1(i2,j2,k2) = mesh%d1(i2,j2,k2) + mass*dx2*dy2*dz2


         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices and weights
         i1  = 1 + int(x(1))
         i2  = 1 + mod(i1,mesh%Nm1d)
         dx1 = i1 - x(1)
         dx2 = 1  - dx1
         j1  = 1 + int(x(2))
         j2  = 1 + mod(j1,mesh%Nm1d)
         dy1 = j1 - x(2)
         dy2 = 1  - dy1
         k1  = 1 + int(x(3))
         k2  = 1 + mod(k1,mesh%Nm1d)
         dz1 = k1 - x(3)
         dz2 = 1  - dz1

         ! Add mass
         mesh%d2(i1,j1,k1) = mesh%d2(i1,j1,k1) + mass*dx1*dy1*dz1
         mesh%d2(i2,j1,k1) = mesh%d2(i2,j1,k1) + mass*dx2*dy1*dz1
         mesh%d2(i1,j2,k1) = mesh%d2(i1,j2,k1) + mass*dx1*dy2*dz1
         mesh%d2(i2,j2,k1) = mesh%d2(i2,j2,k1) + mass*dx2*dy2*dz1
         mesh%d2(i1,j1,k2) = mesh%d2(i1,j1,k2) + mass*dx1*dy1*dz2
         mesh%d2(i2,j1,k2) = mesh%d2(i2,j1,k2) + mass*dx2*dy1*dz2
         mesh%d2(i1,j2,k2) = mesh%d2(i1,j2,k2) + mass*dx1*dy2*dz2
         mesh%d2(i2,j2,k2) = mesh%d2(i2,j2,k2) + mass*dx2*dy2*dz2
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
      real(8)    :: lm1d,mass
      real(8), dimension(3) :: x

      
      ! Init
      lm1d = dble(mesh%Nm1d)
      mass = 1


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x = mod(x_lpt(l,m,n),lm1d)

         ! Indices
         i = 1 + int(x(1))
         j = 1 + int(x(2))
         k = 1 + int(x(3))

         ! Add mass
         mesh%d1(i,j,k) = mesh%d1(i,j,k) + mass

         
         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices
         i = 1 + int(x(1))
         j = 1 + int(x(2))
         k = 1 + int(x(3))

         ! Add mass
         mesh%d2(i,j,k) = mesh%d2(i,j,k) + mass
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
      real(8)    :: dneg,dscale
      complex(8) :: c,t,w
      real(8), dimension(2) :: d,davg,dsig,dmax,dmin
      

      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()


      ! Init
      Ak   = 2*pi/mesh%Nm1d
      imax = mesh%Nm1d   + 2
      jmax = mesh%Nm1d/2 + 1
      kmax = jmax


      ! Forward FFT density fields
      call fft_3d(mesh%d1,'f')
      call fft_3d(mesh%d2,'f')


      ! Shift in Fourier space
      !$omp parallel do           &
      !$omp default(shared)       &
      !$omp private(i,j,k,ip)     &
      !$omp private(kx,ky,kz,c,t,w)
      do k=1,mesh%Nm1d
         if (k <= kmax) then
            kz = Ak*(k-1)
         else
            kz = Ak*(k-1-mesh%Nm1d)
         endif

         do j=1,mesh%Nm1d
            if (j <= jmax) then
               ky = Ak*(j-1)
            else
               ky = Ak*(j-1-mesh%Nm1d)
            endif

            do i=1,mesh%Nm1d+2,2
               ip = i + 1
               kx = Ak*((i-1)/2)

               if (i > 1 .or. j > 1 .or. k > 1) then
                  ! Shift
                  t = cmplx(0.,(kx + ky + kz)/2)
                  c = exp(t)*cmplx(mesh%d2(i,j,k),mesh%d2(ip,j,k))
                  mesh%d2(i ,j,k) = real(c)
                  mesh%d2(ip,j,k) = imag(c)

                  ! Deconvolution
                  w = assignment_transform(kx,ky,kz,lpt%assign)
                  mesh%d1(i:ip,j,k) = mesh%d1(i:ip,j,k)/w
                  mesh%d2(i:ip,j,k) = mesh%d2(i:ip,j,k)/w

                  ! Average
                  mesh%d2(i:ip,j,k) = (mesh%d1(i:ip,j,k) + mesh%d2(i:ip,j,k))/2
               else
                  mesh%d2(i:ip,j,k) = (mesh%d1(i:ip,j,k) + mesh%d2(i:ip,j,k))/2
               endif
            enddo
         enddo
      enddo
      !$omp end parallel do


      ! Inverse FFT density fields
      call fft_3d(mesh%d1,'b')
      call fft_3d(mesh%d2,'b')


      ! Adjust for negative densities
      if (.true.) then
         dneg = 0

         !$omp parallel do          &
         !$omp default(shared)      &
         !$omp private(i,j,k)       &
         !$omp reduction(+:dneg)
         do k=1,mesh%Nm1d
            do j=1,mesh%Nm1d
               do i=1,mesh%Nm1d
                  if (mesh%d1(i,j,k) < 0) then
                     dneg = dneg + mesh%d1(i,j,k)
                  endif
               enddo
            enddo
         enddo
         !$omp end parallel do

         dneg   = dneg/mesh%Nmesh
         dscale = 1/(1 + abs(dneg))

         !$omp parallel do     &
         !$omp default(shared) &
         !$omp private(i,j,k)
         do k=1,mesh%Nm1d
            do j=1,mesh%Nm1d
               do i=1,mesh%Nm1d
                  mesh%d1(i,j,k) = max(mesh%d1(i,j,k),0D0)*dscale
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
      do k=1,mesh%Nm1d
         do j=1,mesh%Nm1d
            do i=1,mesh%Nm1d
               ! Save
               d(1) = mesh%d1(i,j,k)
               d(2) = mesh%d2(i,j,k)
               mesh%rho1(i,j,k) = d(1)
               mesh%rho2(i,j,k) = d(2)

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
      davg = davg/mesh%Nmesh
      dsig = sqrt(dsig/mesh%Nmesh - davg**2)
      write(*,*) 'rho 1 : ',real((/davg(1),dsig(1),dmin(1),dmax(1)/))
      write(*,*) 'rho 2 : ',real((/davg(2),dsig(2),dmin(2),dmax(2)/))


      time2 = time()
      write(*,'(2a)') timing(time1,time2),' : MESH density interlace'
      return
    end subroutine density_interlace

    
  end subroutine mesh_density


!------------------------------------------------------------------------------!
! Velocity
!------------------------------------------------------------------------------!


  subroutine mesh_velocity
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k,n


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Pointers
    ! Do not overwrite mesh%d1=mesh%fft1
    mesh%mom1 => mesh%vel1
    mesh%d1   => mesh%fft1
    mesh%p1   => mesh%fft2
    mesh%p2   => mesh%fft3


    ! Init
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(k)
    !$omp do
    do k=1,mesh%Nm1d
       mesh%d1(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%mom1(:,:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       mesh%mom2(:,:,:,k) = 0
    enddo
    !$omp end do
    !$omp end parallel


    ! Construct velocity, momentum fields from particles
    ! Use cubical domain decomposition, see domain.f90

    select case (lpt%assign)
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


    ! Interlace and average velocity fields
    call velocity_interlace


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : MESH velocity field'
    return


  contains


    subroutine tsc_assignment(indx)
      ! Default
      implicit none


      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx


      ! Local variable
      integer(4) :: a,b,c,l,m,n
      real(8)    :: lm1d,mass
      real(8)    :: dx,dy,dz
      integer(4), dimension(-1:1) :: i,j,k
      real(8),    dimension(-1:1) :: wx,wy,wz
      real(8),    dimension(3)    :: x,v,mv
      

      ! Init
      lm1d = dble(mesh%Nm1d)
      mass = 1


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x  = mod(x_lpt(l,m,n) - 0.5 + lm1d,lm1d)
         v  = v_lpt(l,m,n)*unit%vel/1E5
         mv = mass*v

         ! Indices and weights
         i( 0)  = 1 + int(x(1))
         i(-1)  = 1 + mod(i(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         i( 1)  = 1 + mod(i(0)                ,mesh%Nm1d)
         j( 0)  = 1 + int(x(2))
         j(-1)  = 1 + mod(j(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         j( 1)  = 1 + mod(j(0)                ,mesh%Nm1d)
         k( 0)  = 1 + int(x(3))
         k(-1)  = 1 + mod(k(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         k( 1)  = 1 + mod(k(0)                ,mesh%Nm1d)
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

         ! Add mass and momentum
         do c=-1,1
         do b=-1,1
         do a=-1,1
            mesh%d1(    i(a),j(b),k(c)) = mesh%d1(    i(a),j(b),k(c)) &
                                        + mass*wx(a)*wy(b)*wz(c)
            mesh%mom1(:,i(a),j(b),k(c)) = mesh%mom1(:,i(a),j(b),k(c)) &
                                        + mv*wx(a)*wy(b)*wz(c)
         enddo
         enddo
         enddo


         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices and weights
         i( 0)  = 1 + int(x(1))
         i(-1)  = 1 + mod(i(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         i( 1)  = 1 + mod(i(0)                ,mesh%Nm1d)
         j( 0)  = 1 + int(x(2))
         j(-1)  = 1 + mod(j(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         j( 1)  = 1 + mod(j(0)                ,mesh%Nm1d)
         k( 0)  = 1 + int(x(3))
         k(-1)  = 1 + mod(k(0) - 2 + mesh%Nm1d,mesh%Nm1d)
         k( 1)  = 1 + mod(k(0)                ,mesh%Nm1d)
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

         ! Add momentum
         do c=-1,1
         do b=-1,1
         do a=-1,1
            mesh%mom2(:,i(a),j(b),k(c)) = mesh%mom2(:,i(a),j(b),k(c)) &
                                        + mv*wx(a)*wy(b)*wz(c)
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
      real(8)    :: lm1d,mass
      real(8)    :: dx1,dx2,dy1,dy2,dz1,dz2
      real(8), dimension(3) :: x,v,mv
      

      ! Init
      lm1d = dble(mesh%Nm1d)
      mass = 1


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x  = mod(x_lpt(l,m,n) - 0.5 + lm1d,lm1d)
         v  = v_lpt(l,m,n)*unit%vel/1E5
         mv = mass*v

         ! Indices and weights
         i1  = 1 + int(x(1))
         i2  = 1 + mod(i1,mesh%Nm1d)
         dx1 = i1 - x(1)
         dx2 = 1  - dx1
         j1  = 1 + int(x(2))
         j2  = 1 + mod(j1,mesh%Nm1d)
         dy1 = j1 - x(2)
         dy2 = 1  - dy1
         k1  = 1 + int(x(3))
         k2  = 1 + mod(k1,mesh%Nm1d)
         dz1 = k1 - x(3)
         dz2 = 1  - dz1

         ! Add mass
         mesh%d1(i1,j1,k1) = mesh%d1(i1,j1,k1) + mass*dx1*dy1*dz1
         mesh%d1(i2,j1,k1) = mesh%d1(i2,j1,k1) + mass*dx2*dy1*dz1
         mesh%d1(i1,j2,k1) = mesh%d1(i1,j2,k1) + mass*dx1*dy2*dz1
         mesh%d1(i2,j2,k1) = mesh%d1(i2,j2,k1) + mass*dx2*dy2*dz1
         mesh%d1(i1,j1,k2) = mesh%d1(i1,j1,k2) + mass*dx1*dy1*dz2
         mesh%d1(i2,j1,k2) = mesh%d1(i2,j1,k2) + mass*dx2*dy1*dz2
         mesh%d1(i1,j2,k2) = mesh%d1(i1,j2,k2) + mass*dx1*dy2*dz2
         mesh%d1(i2,j2,k2) = mesh%d1(i2,j2,k2) + mass*dx2*dy2*dz2

         ! Add momentum
         mesh%mom1(:,i1,j1,k1) = mesh%mom1(:,i1,j1,k1) + mv*dx1*dy1*dz1
         mesh%mom1(:,i2,j1,k1) = mesh%mom1(:,i2,j1,k1) + mv*dx2*dy1*dz1
         mesh%mom1(:,i1,j2,k1) = mesh%mom1(:,i1,j2,k1) + mv*dx1*dy2*dz1
         mesh%mom1(:,i2,j2,k1) = mesh%mom1(:,i2,j2,k1) + mv*dx2*dy2*dz1
         mesh%mom1(:,i1,j1,k2) = mesh%mom1(:,i1,j1,k2) + mv*dx1*dy1*dz2
         mesh%mom1(:,i2,j1,k2) = mesh%mom1(:,i2,j1,k2) + mv*dx2*dy1*dz2
         mesh%mom1(:,i1,j2,k2) = mesh%mom1(:,i1,j2,k2) + mv*dx1*dy2*dz2
         mesh%mom1(:,i2,j2,k2) = mesh%mom1(:,i2,j2,k2) + mv*dx2*dy2*dz2


         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices and weights
         i1  = 1 + int(x(1))
         i2  = 1 + mod(i1,mesh%Nm1d)
         dx1 = i1 - x(1)
         dx2 = 1  - dx1
         j1  = 1 + int(x(2))
         j2  = 1 + mod(j1,mesh%Nm1d)
         dy1 = j1 - x(2)
         dy2 = 1  - dy1
         k1  = 1 + int(x(3))
         k2  = 1 + mod(k1,mesh%Nm1d)
         dz1 = k1 - x(3)
         dz2 = 1  - dz1

         ! Add momentum
         mesh%mom2(:,i1,j1,k1) = mesh%mom2(:,i1,j1,k1) + mv*dx1*dy1*dz1
         mesh%mom2(:,i2,j1,k1) = mesh%mom2(:,i2,j1,k1) + mv*dx2*dy1*dz1
         mesh%mom2(:,i1,j2,k1) = mesh%mom2(:,i1,j2,k1) + mv*dx1*dy2*dz1
         mesh%mom2(:,i2,j2,k1) = mesh%mom2(:,i2,j2,k1) + mv*dx2*dy2*dz1
         mesh%mom2(:,i1,j1,k2) = mesh%mom2(:,i1,j1,k2) + mv*dx1*dy1*dz2
         mesh%mom2(:,i2,j1,k2) = mesh%mom2(:,i2,j1,k2) + mv*dx2*dy1*dz2
         mesh%mom2(:,i1,j2,k2) = mesh%mom2(:,i1,j2,k2) + mv*dx1*dy2*dz2
         mesh%mom2(:,i2,j2,k2) = mesh%mom2(:,i2,j2,k2) + mv*dx2*dy2*dz2
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
      real(8)    :: lm1d,mass
      real(8), dimension(3) :: x,v,mv

      
      ! Init
      lm1d = dble(mesh%Nm1d)
      mass = 1


      ! Loop over cells in domain
      do n=indx(1,3),indx(2,3)
      do m=indx(1,2),indx(2,2)
      do l=indx(1,1),indx(2,1)
         ! Particle
         x  = mod(x_lpt(l,m,n),lm1d)
         v  = v_lpt(l,m,n)*unit%vel/1E5
         mv = mass*v


         ! Indices
         i = 1 + int(x(1))
         j = 1 + int(x(2))
         k = 1 + int(x(3))

         ! Add mass
         mesh%d1(    i,j,k) = mesh%d1(    i,j,k) + mass
         mesh%mom1(:,i,j,k) = mesh%mom1(:,i,j,k) + mv

         
         ! Shift
         x = mod(x + 0.5,lm1d)

         ! Indices
         i = 1 + int(x(1))
         j = 1 + int(x(2))
         k = 1 + int(x(3))

         ! Add momentum
         mesh%mom2(:,i,j,k) = mesh%mom2(:,i,j,k) + mv
      enddo
      enddo
      enddo


      return
    end subroutine ngp_assignment


    subroutine velocity_interlace
      ! Default
      implicit none


      ! Local variables
      integer(4) :: i,j,k,n,ip
      integer(4) :: imax,jmax,kmax
      real(8)    :: Ak,kx,ky,kz
      complex(8) :: c,t,w
      real(8), dimension(3) :: v,vavg,vsig,vmin,vmax
      real(8), dimension(3) :: p,pavg,psig,pmin,pmax
    

      ! Timing variables
      integer(4) :: time1,time2
      time1 = time()


      ! Init
      Ak   = 2*pi/mesh%Nm1d
      imax = mesh%Nm1d   + 2
      jmax = mesh%Nm1d/2 + 1
      kmax = jmax


      ! Loop over directions
      do n=1,3
         ! Copy
         !$omp parallel do     &
         !$omp default(shared) &
         !$omp private(i,j,k)
         do k=1,mesh%Nm1d
            do j=1,mesh%Nm1d
               do i=1,mesh%Nm1d
                  mesh%p1(i,j,k) = mesh%mom1(n,i,j,k)
                  mesh%p2(i,j,k) = mesh%mom2(n,i,j,k)
               enddo
            enddo
         enddo
         !$omp end parallel do


         ! Forward FFT momentum fields
         call fft_3d(mesh%p1,'f')
         call fft_3d(mesh%p2,'f')


         ! Shift in Fourier space
         !$omp parallel do           &
         !$omp default(shared)       &
         !$omp private(i,j,k,ip)     &
         !$omp private(kx,ky,kz,c,t,w)
         do k=1,mesh%Nm1d
            if (k <= kmax) then
               kz = Ak*(k-1)
            else
               kz = Ak*(k-1-mesh%Nm1d)
            endif

            do j=1,mesh%Nm1d
               if (j <= jmax) then
                  ky = Ak*(j-1)
               else
                  ky = Ak*(j-1-mesh%Nm1d)
               endif

               do i=1,mesh%Nm1d+2,2
                  ip = i + 1
                  kx = Ak*((i-1)/2)

                  ! Skip k=0 mode
                  if (i > 1 .or. j > 1 .or. k > 1) then
                     ! Shift
                     t = cmplx(0.,(kx + ky + kz)/2)
                     c = exp(t)*cmplx(mesh%p2(i,j,k),mesh%p2(ip,j,k))
                     mesh%p2(i ,j,k) = real(c)
                     mesh%p2(ip,j,k) = imag(c)

                     ! Deconvolution
                     w = assignment_transform(kx,ky,kz,lpt%assign)

                     ! Average
                     mesh%p2(i:ip,j,k) = (mesh%p1(i:ip,j,k) &
                                       +  mesh%p2(i:ip,j,k))/(2*w)
                  else
                     mesh%p2(i:ip,j,k) = (mesh%p1(i:ip,j,k) &
                                       +  mesh%p2(i:ip,j,k))/2
                  endif
               enddo
            enddo
         enddo
         !$omp end parallel do


         ! Inverse FFT momentum field
         call fft_3d(mesh%p2,'b')


         ! Save mom2
         !$omp parallel do     &
         !$omp default(shared) &
         !$omp private(i,j,k)
         do k=1,mesh%Nm1d
            do j=1,mesh%Nm1d
               do i=1,mesh%Nm1d
                  mesh%mom2(n,i,j,k) = mesh%p2(i,j,k)
               enddo
            enddo
         enddo
         !$omp end parallel do
      enddo


      ! Velocity, momentum fields
      vavg = 0
      vsig = 0
      vmax = 0
      vmin = huge(vmin)
      pavg = 0
      psig = 0
      pmax = 0
      pmin = huge(pmin)

      !$omp parallel do                      &
      !$omp default(shared)                  &
      !$omp private(i,j,k,p,v)               &
      !$omp reduction(+:vavg,vsig,pavg,psig) &
      !$omp reduction(max:vmax,pmax)         &
      !$omp reduction(min:vmin,pmin)
      do k=1,mesh%Nm1d
         do j=1,mesh%Nm1d
            do i=1,mesh%Nm1d
               ! Divide by density
               if (mesh%d1(i,j,k) > 1E-6) then
                  v = mesh%mom1(:,i,j,k)/mesh%d1(i,j,k)
               else
                  v = 0
               endif
               mesh%vel1(:,i,j,k) = v
               p = mesh%mom2(:,i,j,k)

               ! Stats
               vavg = vavg + v
               vsig = vsig + v**2
               vmax = max(vmax,v)
               vmin = min(vmin,v)
               pavg = pavg + p
               psig = psig + p**2
               pmax = max(pmax,p)
               pmin = min(pmin,p)
            enddo
         enddo
      enddo


      ! Write to screen
      ! v [km/s]
      ! p = (1+delta)v [km/s]
      vavg = vavg/mesh%Nmesh
      vsig = sqrt(vsig/mesh%Nmesh - vavg**2)
      pavg = pavg/mesh%Nmesh
      psig = sqrt(psig/mesh%Nmesh - pavg**2)
      write(*,*) 'vx1 : ',real((/vavg(1),vsig(1),vmin(1),vmax(1)/))
      write(*,*) 'vy1 : ',real((/vavg(2),vsig(2),vmin(2),vmax(2)/))
      write(*,*) 'vz1 : ',real((/vavg(3),vsig(3),vmin(3),vmax(3)/))
      write(*,*) 'px2 : ',real((/pavg(1),psig(1),pmin(1),pmax(1)/))
      write(*,*) 'py2 : ',real((/pavg(2),psig(2),pmin(2),pmax(2)/))
      write(*,*) 'pz2 : ',real((/pavg(3),psig(3),pmin(3),pmax(3)/))


      time2 = time()
      write(*,'(2a)') timing(time1,time2),' : MESH velocity interlace'
      return
    end subroutine velocity_interlace

    
  end subroutine mesh_velocity


end module meshmake_module
