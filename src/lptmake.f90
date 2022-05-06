module lptmake_module
  ! Intel
  use IFPORT
  use OMP_LIB
  

  ! Modules
  use lpt_module
  use constant_module
  use helper_module
  use mkl_module
  use timing_module
  use cosmo_module     , only : cosmo
  use cosmology_module , only : cosmo_calc
  use grf_module       , only : grf
  use mesh_module      , only : mesh
  use sim_module       , only : sim
  use simulation_module, only : sim_calc


  ! Default
  implicit none
  public


contains


  subroutine lpt_make
    ! Default
    implicit none


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Pointers
    ! lpt%delta1 and grf%grf share mesh%fft1
    lpt%delta1 => mesh%fft1
    lpt%delta2 => mesh%fft2
    lpt%gphi   => mesh%fft3


    ! Init comoving
    cosmo%a = 1
    cosmo%z = 0
    call cosmo_calc
    call sim_calc
    

    ! LPT
    if (lpt%make == 'read') then
       ! IO
       call lpt_read
    else
       ! Linear density field
       call lpt_delta

       ! 1LPT or 2LPT
       call lpt_gradphi1
       if (lpt%order == '2lpt') call lpt_gradphi2

       ! IO
       if (lpt%make == 'write') call lpt_write
    endif

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : LPT make'
    return
  end subroutine lpt_make


!------------------------------------------------------------------------------!
! Linear density field
!------------------------------------------------------------------------------!


  subroutine lpt_delta
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: i,j,k,kk,Nk,un
    integer(4)    :: imax,jmax,kmax
    real(8)       :: Ak,kr,kx,ky,kz,dk
    real(8)       :: p1,p2,w
    character(80) :: fn
    type(df_task) :: task
    real(8), dimension(1) :: x
    real(8), allocatable, dimension(:)   :: klin,plin
    real(8), allocatable, dimension(:,:) :: pow


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Allocate
    allocate(pow(3,lpt%Nm1d))


    ! Init
    Ak   = 2*pi/lpt%Nm1d
    imax = lpt%Nm1d   + 2
    jmax = lpt%Nm1d/2 + 1
    kmax = jmax
    pow  = 0


    ! Spline
    allocate(klin(cosmo%Nk))
    allocate(plin(cosmo%Nk))
    klin = cosmo%Plin(1,:)
    plin = cosmo%Plin(3,:)
    call spline_construct(task,klin,plin)


    ! Fourier modes delta(k)
    !$omp parallel do                    &
    !$omp default(shared)                &
    !$omp private(i,j,k,kk)              &
    !$omp private(kr,kx,ky,kz,p1,p2,w,x) &
    !$omp reduction(+:pow)
    do k=1,lpt%Nm1d
       if (k <= kmax) then
          kz = Ak*(k-1)
       else
          kz = Ak*(k-1-lpt%Nm1d)
       endif

       do j=1,lpt%Nm1d
          if (j <= jmax) then
             ky = Ak*(j-1)
          else
             ky = Ak*(j-1-lpt%Nm1d)
          endif

          do i=1,imax,2
             kx = Ak*((i-1)/2)
             kr = sqrt(kx**2 + ky**2 + kz**2)

             if (i > 1 .or. j > 1 .or. k > 1) then
                ! FFT weight
                if (i == 1 .or. i == lpt%Nm1d + 1) then
                   w = 1
                else
                   w = 2
                endif

                ! Interpolate
                x  = kr*lpt%Nm1d/cosmo%Lbox
                x  = spline_interp(task,x)
                p1 = (2*pi**2/kr**3)*x(1)

                ! Fourier mode
                lpt%delta1(i:i+1,j,k) = sqrt(p1)*grf%grf(i:i+1,j,k)
                p2 = sum((lpt%delta1(i:i+1,j,k)/lpt%Nmesh)**2)

                ! Bin power
                kk        = nint(kr/Ak)
                pow(1,kk) = pow(1,kk) + w
                pow(2,kk) = pow(2,kk) + w*p1/lpt%Nmesh
                pow(3,kk) = pow(3,kk) + w*p2
             else
                lpt%delta1(i:i+1,j,k) = 0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Linear power spectrum
    do k=1,lpt%Nm1d
       if (pow(1,k) > 0) then
          pow(2:3,k) = pow(2:3,k)/pow(1,k)
       else
          pow(2:3,k) = 0
       endif

       lpt%Pinit(1  ,k) = 2*pi*k/cosmo%Lbox
       lpt%Pinit(2:3,k) = 4*pi*dble(k)**3*pow(2:3,k)
    enddo


    ! Write to file
    un = 11
    fn = trim(lpt%dirout)//'power_linear.txt'
    write(*,*) 'Writing ',trim(fn)
    
    open(un,file=fn,recl=100)
    write(un,'(3a14)') 'k','D^2_lin','D^2_lpt'
    
    do k=1,lpt%Nm1d/2
       write(un,'(3es14.6)') lpt%Pinit(:,k)
    enddo
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : LPT delta field'
    return
  end subroutine lpt_delta


!------------------------------------------------------------------------------!
! 1LPT/Zel'dovich approximation
!------------------------------------------------------------------------------!


  subroutine lpt_gradphi1
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k,n
    integer(4) :: imax,jmax,kmax
    real(8)    :: Ak,kx,ky,kz,ksq,w


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init
    Ak   = 2*pi/lpt%Nm1d
    imax = lpt%Nm1d   + 2
    jmax = lpt%Nm1d/2 + 1
    kmax = jmax


    ! Calculate Grad phi_1 fields
    do n=1,3
       !$omp parallel do           &
       !$omp default(shared)       &
       !$omp private(i,j,k)        &
       !$omp private(kx,ky,kz,ksq,w)
       do k=1,lpt%Nm1d
          if (k <= kmax) then
             kz = Ak*(k-1)
          else
             kz = Ak*(k-1-lpt%Nm1d)
          endif

          do j=1,lpt%Nm1d
             if (j <= jmax) then
                ky = Ak*(j-1)
             else
                ky = Ak*(j-1-lpt%Nm1d)
             endif

             do i=1,imax,2
                kx  = Ak*((i-1)/2)
                ksq = kx**2 + ky**2 + kz**2

                ! Force kernels
                if (i > 1 .or. j > 1 .or. k > 1) then
                   select case (n)
                   case (1)
                      w = -kx/ksq
                   case (2)
                      w = -ky/ksq
                   case (3)
                      w = -kz/ksq
                   end select
                else
                   w = 0
                endif


                ! Complex multiply
                lpt%gphi(i  ,j,k) = -w*lpt%delta1(i+1,j,k)
                lpt%gphi(i+1,j,k) =  w*lpt%delta1(i  ,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do


       ! Inverse FFT Grad phi_1 fields
       ! See fft.f90
       call fft_3d(lpt%gphi,'b')


       ! Save
       !$omp parallel do     &
       !$omp default(shared) &
       !$omp private(i,j,k)
       do k=1,lpt%Nm1d
          do j=1,lpt%Nm1d
             do i=1,lpt%Nm1d
                lpt%gradphi(n,i,j,k) = lpt%gphi(i,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do
    enddo


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : LPT gradphi 1'
    return
  end subroutine lpt_gradphi1


!------------------------------------------------------------------------------!
! 2LPT
!------------------------------------------------------------------------------!


  subroutine lpt_gradphi2
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k,n
    integer(4) :: imax,jmax,kmax
    real(8)    :: Ak,kx,ky,kz,ksq,w


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init
    Ak   = 2*pi/lpt%Nm1d
    imax = lpt%Nm1d   + 2
    jmax = lpt%Nm1d/2 + 1
    kmax = jmax


    ! Calc Grad^2 phi_ii fields
    do n=1,3
       !$omp parallel do           &
       !$omp default(shared)       &
       !$omp private(i,j,k)        &
       !$omp private(kx,ky,kz,ksq,w)
       do k=1,lpt%Nm1d
          if (k <= kmax) then
             kz = Ak*(k-1)
          else
             kz = Ak*(k-1-lpt%Nm1d)
          endif

          do j=1,lpt%Nm1d
             if (j <= jmax) then
                ky = Ak*(j-1)
             else
                ky = Ak*(j-1-lpt%Nm1d)
             endif

             do i=1,imax,2
                kx  = Ak*((i-1)/2)
                ksq = kx**2 + ky**2 + kz**2

                ! Calc kernels
                if (i > 1 .or. j > 1 .or. k > 1) then
                   select case (n)
                   case (1)
                      w = -kx**2/ksq
                   case (2)
                      w = -ky**2/ksq
                   case (3)
                      w = -kz**2/ksq
                   end select
                else
                   w = 0
                endif


                ! Complex multiply
                lpt%gphi(i:i+1,j,k) = w*lpt%delta1(i:i+1,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do


       ! Inverse FFT Grad^2 phi_ii fields
       ! See fft.f90
       call fft_3d(lpt%gphi,'b')


       ! Save in gradphi array
       !$omp parallel do     &
       !$omp default(shared) &
       !$omp private(i,j,k)
       do k=1,lpt%Nm1d
          do j=1,lpt%Nm1d
             do i=1,lpt%Nm1d
                lpt%gradphi(n+3,i,j,k) = lpt%gphi(i,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do
    enddo


    ! Calc delta2
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k)
    do k=1,lpt%Nm1d
       do j=1,lpt%Nm1d
          do i=1,lpt%Nm1d
             lpt%delta2(i,j,k) = lpt%gradphi(4,i,j,k)*lpt%gradphi(5,i,j,k) &
                               + lpt%gradphi(5,i,j,k)*lpt%gradphi(6,i,j,k) &
                               + lpt%gradphi(6,i,j,k)*lpt%gradphi(4,i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Calculate Grad^2 phi_ij fields
    do n=1,3
       !$omp parallel do           &
       !$omp default(shared)       &
       !$omp private(i,j,k)        &
       !$omp private(kx,ky,kz,ksq,w)
       do k=1,lpt%Nm1d
          if (k <= kmax) then
             kz = Ak*(k-1)
          else
             kz = Ak*(k-1-lpt%Nm1d)
          endif

          do j=1,lpt%Nm1d
             if (j <= jmax) then
                ky = Ak*(j-1)
             else
                ky = Ak*(j-1-lpt%Nm1d)
             endif

             do i=1,imax,2
                kx  = Ak*((i-1)/2)
                ksq = kx**2 + ky**2 + kz**2

                ! Calc kernels
                if (i > 1 .or. j > 1 .or. k > 1) then
                   select case (n)
                   case (1)
                      w = -kx*ky/ksq
                   case (2)
                      w = -ky*kz/ksq
                   case (3)
                      w = -kz*kx/ksq
                   end select
                else
                   w = 0
                endif
                
                ! Complex multiply
                lpt%gphi(i:i+1,j,k) = w*lpt%delta1(i:i+1,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do


       ! Inverse FFT Grad^2 phi_ij fields
       call fft_3d(lpt%gphi,'b')


       ! Save in gradphi array
       !$omp parallel do     &
       !$omp default(shared) &
       !$omp private(i,j,k)
       do k=1,lpt%Nm1d
          do j=1,lpt%Nm1d
             do i=1,lpt%Nm1d
                lpt%gradphi(n+3,i,j,k) = lpt%gphi(i,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do
    enddo


    ! Add to delta2
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k)
    do k=1,lpt%Nm1d
       do j=1,lpt%Nm1d
          do i=1,lpt%Nm1d
             lpt%delta2(i,j,k) = lpt%delta2(i,j,k)       &
                               - lpt%gradphi(4,i,j,k)**2 &
                               - lpt%gradphi(5,i,j,k)**2 &
                               - lpt%gradphi(6,i,j,k)**2
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Forward FFT delta2 field
    ! See fft.f90
    call fft_3d(lpt%delta2,'f')


    ! Calculate Grad phi_2
    do n=1,3
       !$omp parallel do           &
       !$omp default(shared)       &
       !$omp private(i,j,k)        &
       !$omp private(kx,ky,kz,ksq,w)
       do k=1,lpt%Nm1d
          if (k <= kmax) then
             kz = Ak*(k-1)
          else
             kz = Ak*(k-1-lpt%Nm1d)
          endif

          do j=1,lpt%Nm1d
             if (j <= jmax) then
                ky = Ak*(j-1)
             else
                ky = Ak*(j-1-lpt%Nm1d)
             endif

             do i=1,imax,2
                kx  = Ak*((i-1)/2)
                ksq = kx**2 + ky**2 + kz**2

                ! Kernels
                if (i > 1 .or. j > 1 .or. k > 1) then
                   select case (n)
                   case (1)
                      w = -kx/ksq
                   case (2)
                      w = -ky/ksq
                   case (3)
                      w = -kz/ksq
                   end select
                else
                   w = 0
                endif

                ! Complex multiply
                lpt%gphi(i  ,j,k) = -w*lpt%delta2(i+1,j,k)
                lpt%gphi(i+1,j,k) =  w*lpt%delta2(i  ,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do


       ! Inverse FFT Grad phi_2 fields
       call fft_3d(lpt%gphi,'b')


       ! Save in gradphi array
       !$omp parallel do     &
       !$omp default(shared) &
       !$omp private(i,j,k)
       do k=1,lpt%Nm1d
          do j=1,lpt%Nm1d
             do i=1,lpt%Nm1d
                lpt%gradphi(n+3,i,j,k) = lpt%gphi(i,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do
    enddo


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : LPT gradphi 2'
    return
  end subroutine lpt_gradphi2


end module lptmake_module
