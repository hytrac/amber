module cosmology_module
  ! Intel
  use OMP_LIB


  ! Modules
  use cosmo_module
  use constant_module
  use helper_module
  use mkl_module
  use timing_module
  use esf_module  , only : esf
  use lpt_module  , only : lpt
  use mesh_module , only : mesh
  use reion_module, only : reion
  use sim_module  , only : sim
  

  ! Default
  implicit none
  public


contains


  subroutine cosmo_calc
    ! Default
    implicit none


    ! Local variables
    real(8) :: a,z,hsq,oma
    real(8) :: gf(2),f1,f2


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Cosmo
    a        = cosmo%a
    z        = cosmo%z
    hsq      = cosmo%or/a**4 + cosmo%om/a**3 + cosmo%ol
    oma      = cosmo%om/a**3/hsq**2
    cosmo%H0 = H0_cgs*cosmo%h
    cosmo%Hz = cosmo%H0*sqrt(hsq)
    cosmo%fb = cosmo%ob/cosmo%om
    cosmo%fc = 1 - cosmo%fb


    ! Comoving densities (cgs)
    cosmo%rhocrit0 = 3*cosmo%H0**2/(8*pi*G_cgs)
    cosmo%rhom0    = rhoc_cgs*cosmo%om*cosmo%h**2
    cosmo%rhob0    = cosmo%fb*cosmo%rhom0
    cosmo%rhodm0   = cosmo%fc*cosmo%rhom0
    cosmo%nH0      = cosmo%XH *cosmo%rhob0/mH_cgs
    cosmo%nHe0     = cosmo%YHe*cosmo%rhob0/mHe_cgs
    cosmo%ne0      = cosmo%nH0 + 2*cosmo%nHe0


    ! Proper densities (cgs)
    cosmo%rhocrit = 3*cosmo%Hz**2/(8*pi*G_cgs)
    cosmo%rhom    = cosmo%rhom0 /a**3
    cosmo%rhob    = cosmo%rhob0 /a**3
    cosmo%rhodm   = cosmo%rhodm0/a**3
    cosmo%nH      = cosmo%nH0   /a**3
    cosmo%nHe     = cosmo%nHe0  /a**3
    cosmo%ne      = cosmo%ne0   /a**3


    ! Halo variables
    cosmo%rho200a = 200*cosmo%rhom
    cosmo%rho200c = 200*cosmo%rhocrit
    cosmo%rho500c = 500*cosmo%rhocrit
    cosmo%rhovir  = (18*pi**2 + 82*(oma-1) - 39*(oma-1)**2)*cosmo%rhocrit


    ! Density and velocity growth factors
    gf       = growthfactors_of_z(z)
    cosmo%D1 = gf(1)/cosmo%delta0
    cosmo%D2 = -3./7*cosmo%D1**2*oma**(-1./143)
    f1       = gf(2)
    f2       = 2*oma**(6./11)
    cosmo%v1 = cosmo%D1*f1*cosmo%Hz
    cosmo%v2 = cosmo%D2*f2*cosmo%Hz


    time2 = time()
    !write(*,'(2a)') timing(time1,time2),' : Cosmo calc'
    return
  end subroutine cosmo_calc


!------------------------------------------------------------------------------!
! Linear power spectrum
!------------------------------------------------------------------------------!


  subroutine cosmo_linearpowerspectrum
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: k,k1,k2,un
    real(8)       :: kr,dlnk,pk,var
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Read dimensionless power spectrum
    ! First line assumed to be header    
    un = 11
    fn = trim(sim%dirin)//cosmo%file
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn)
 
    ! First loop through file to determine length
    k = 0
    do
       read(un,*,end=1)
       k = k + 1
    enddo
1   continue
    cosmo%Nk = k - 1


    ! Allocate power spectrum array
    allocate(cosmo%Plin(3,cosmo%Nk))


    ! Skip first line, read k and power
    rewind(un)
    read(un,*)

    do k=1,cosmo%Nk
       read(un,*) cosmo%Plin(1,k),cosmo%Plin(3,k)
    enddo
    close(un)


    ! Calculate dlnk, dk, sigma8
    var = 0

    do k=1,cosmo%Nk
       if (k == 1) then
          k1 = k
          k2 = k + 1
       else if (k == cosmo%Nk) then
          k1 = k - 1
          k2 = k
       else
          k1 = k - 1
          k2 = k + 1
       endif

       kr   = cosmo%Plin(1,k)
       pk   = cosmo%Plin(3,k)
       dlnk = log(cosmo%Plin(1,k2)/cosmo%Plin(1,k1))/(k2-k1)
       var  = var + pk*tophat_transform(kr*8)**2*dlnk
       cosmo%Plin(2,k) = kr*dlnk
    enddo


    ! Renormalize sigma8
    cosmo%Plin(3,:) = cosmo%Plin(3,:)*(cosmo%s8**2/var)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : COSMO linear power spectrum'
    return
  end subroutine cosmo_linearpowerspectrum


!------------------------------------------------------------------------------!
! Power spectrum
!------------------------------------------------------------------------------!


  subroutine cosmo_powerspectrum
    ! Default
    implicit none


    ! Local variables
    integer(4)     :: i,j,k,ip,kk,un
    integer(4)     :: imax,jmax,kmax
    real(8)        :: Ak,kr,kx,ky,kz
    real(8)        :: pmm,phh,pzz,phm,pzm,pzh
    real(8)        :: havg,zavg,w
    character(100) :: fn
    real(8), allocatable, dimension(:,:) :: pow


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Use interlaced and deconvolved fields
    ! mesh%rho2, esf%rho2


    ! Allocate
    allocate(pow(13,mesh%Nm1d))


    ! Init
    Ak   = 2*pi/cosmo%Lbox
    imax = mesh%Nm1d   + 2
    jmax = mesh%Nm1d/2 + 1
    kmax = jmax
    pow  = 0


    ! Average halo density and zre
    havg = 0
    zavg = 0
    
    !$omp parallel             &
    !$omp default(shared)      &
    !$omp private(i,j,k)       &
    !$omp reduction(+:havg,zavg)
    !$omp do
    do k=1,mesh%Nm1d
       do j=1,mesh%Nm1d
          do i=1,mesh%Nm1d
             havg = havg + esf%rho2(i,j,k)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       do j=1,mesh%Nm1d
          do i=1,mesh%Nm1d
             zavg = zavg + reion%zre(i,j,k)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    havg = havg/mesh%Nmesh
    zavg = zavg/mesh%Nmesh

    ! delta for mass, halo, zre
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(i,j,k)
    !$omp do
    do k=1,mesh%Nm1d
       do j=1,mesh%Nm1d
          do i=1,mesh%Nm1d
             mesh%fft1(i,j,k) = mesh%rho2(i,j,k) - 1
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       do j=1,mesh%Nm1d
          do i=1,mesh%Nm1d
             mesh%fft2(i,j,k) = esf%rho2(i,j,k)/havg - 1
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=1,mesh%Nm1d
       do j=1,mesh%Nm1d
          do i=1,mesh%Nm1d
             mesh%fft3(i,j,k) = (1+reion%zre(i,j,k))/(1+zavg) - 1
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel


    ! Forward FFT fields
    ! See fft.f90
    call fft_3d(mesh%fft1,'f')
    call fft_3d(mesh%fft2,'f')
    call fft_3d(mesh%fft3,'f')


    ! Auto and cross power
    !$omp parallel do                        & 
    !$omp default(shared)                    &
    !$omp private(i,j,k,ip,kk)               &
    !$omp private(kx,ky,kz,kr)               &
    !$omp private(pmm,phh,pzz,phm,pzm,pzh,w) &
    !$omp reduction(+:pow)
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

          do i=1,imax,2
             ip = ip + 1
             kx = Ak*((i-1)/2)
             kr = sqrt(kx**2 + ky**2 + kz**2)

             ! Bin Fourier modes
             if (i > 1 .or. j > 1 .or. k > 1) then
                ! FFT weight
                if (i == 1 .or. i == mesh%Nm1d + 1) then
                   w = 1
                else
                   w = 2
                endif

                ! Power
                pmm = sum(mesh%fft1(i:ip,j,k)/mesh%Nmesh &
                         *mesh%fft1(i:ip,j,k)/mesh%Nmesh)
                phh = sum(mesh%fft2(i:ip,j,k)/mesh%Nmesh &
                         *mesh%fft2(i:ip,j,k)/mesh%Nmesh)
                pzz = sum(mesh%fft3(i:ip,j,k)/mesh%Nmesh &
                         *mesh%fft3(i:ip,j,k)/mesh%Nmesh)                
                phm = sum(mesh%fft1(i:ip,j,k)/mesh%Nmesh &
                         *mesh%fft2(i:ip,j,k)/mesh%Nmesh)
                pzm = sum(mesh%fft1(i:ip,j,k)/mesh%Nmesh &
                         *mesh%fft3(i:ip,j,k)/mesh%Nmesh)
                pzh = sum(mesh%fft2(i:ip,j,k)/mesh%Nmesh &
                         *mesh%fft3(i:ip,j,k)/mesh%Nmesh)                

                ! Add to bin
                kk        = nint(kr/Ak)
                pow(1,kk) = pow(1,kk) + w
                pow(2,kk) = pow(2,kk) + w*pmm
                pow(3,kk) = pow(3,kk) + w*phh
                pow(4,kk) = pow(4,kk) + w*pzz
                pow(5,kk) = pow(5,kk) + w*phm
                pow(6,kk) = pow(6,kk) + w*pzm
                pow(7,kk) = pow(7,kk) + w*pzh
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Power and bias
    do k=1,mesh%Nm1d
       ! Divide by weights
       if (pow(1,k) > 0) then
          pow(2:7,k) = pow(2:7,k)/pow(1,k)
       else
          pow(2:7,k) = 0
       endif


       ! Bias and stochasticity
       if (pow(1,k) > 0) then
          ! matter-halo
          pow(8,k) = sqrt(pow(3,k)/pow(2,k))
          pow(9,k) = pow(5,k)/sqrt(pow(2,k)*pow(3,k))

          ! matter-zre
          pow(10,k) = sqrt(pow(4,k)/pow(2,k))
          pow(11,k) = pow(6,k)/sqrt(pow(2,k)*pow(4,k))

          ! halo-zre
          pow(12,k) = sqrt(pow(4,k)/pow(3,k))
          pow(13,k) = pow(7,k)/sqrt(pow(3,k)*pow(4,k))          
       endif


       ! Save
       pow(  1,k) = Ak*k
       pow(2:7,k) = 4*pi*dble(k)**3*pow(2:7,k)
    enddo


    ! Write to file
    un = 11
    fn = trim(cosmo%dir)//'power_'//trim(sim%zstr)//'.txt'
    write(*,*) 'Writing ',trim(fn)
    
    open(un,file=fn,recl=1000)
    write(un,'(12a14)') 'k','D^2_lin','D^2_lpt',   &
                        'D^2_mm','D^2_hh','D^2_zz', &
                        'b_hm','r_hm','b_zm','r_zm','b_zh','r_zh'
    
    do k=1,mesh%Nm1d/2
       write(un,'(12es14.6)') lpt%Pinit(1,k),lpt%Pinit(2:3,k)*cosmo%D1**2, &
                              pow(2:4,k),pow(8:13,k)
    enddo
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : COSMO power spectrum'
    return
  end subroutine cosmo_powerspectrum

  
end module cosmology_module
