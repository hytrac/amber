module grfmake_module
  ! Intel
  use OMP_LIB
  

  ! Modules
  use grf_module
  use constant_module
  use helper_module
  use mkl_module
  use timing_module
  use cosmo_module, only : cosmo
  use mesh_module , only : mesh


  ! Default
  implicit none
  public


contains


  subroutine grf_make
    ! Default
    implicit none


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Pointers
    grf%grf => mesh%fft1

    
    ! Gaussian random field
    if (grf%make == 'read') then
       ! IO
       call grf_read
    else
       call grf_calc

       ! IO
       if (grf%make == 'write') call grf_write
    endif


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : GRF make'
    return
  end subroutine grf_make


!------------------------------------------------------------------------------!
! Gaussian random field
!------------------------------------------------------------------------------!


  subroutine grf_calc
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k,n
    integer(4) :: brng,method,errcode
    real(8)    :: a,b,sigma
    real(8)    :: d,davg,dsig,dmax,dmin
    type(VSL_STREAM_STATE) :: stream
    integer(4), allocatable, dimension(:) :: iseed
    real(8),    allocatable, dimension(:) :: xseed


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Allocate
    allocate(iseed(grf%Nm1d))
    allocate(xseed(grf%Nm1d))


    ! Init seeds
    iseed(1) = grf%seed
    a        = 0
    b        = 1
    brng     = VSL_BRNG_MT19937
    method   = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
    errcode  = vslnewstream(stream,brng,iseed(1))
    errcode  = vdrnguniform(method,stream,grf%Nm1d,xseed,a,b)
    iseed    = int(xseed*huge(0))


    ! Generate Gaussian random numbers
    brng   = VSL_BRNG_MT19937
    method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
    a      = 0
    sigma  = 1

    !$omp parallel do               &
    !$omp default(shared)           &
    !$omp private(j,k,stream,errcode)
    do k=1,grf%Nm1d
       ! Initialize stream
       errcode = vslnewstream(stream,brng,iseed(k))
       !call CheckVslError(errcode)

       ! Generate random numbers
       do j=1,grf%Nm1d
          errcode = vdrnggaussian(method,stream, &
                                  grf%Nm1d,grf%grf(1:grf%Nm1d,j,k),a,sigma)
          !call CheckVslError(errcode)
       enddo

       ! Delete stream
       errcode = vsldeletestream(stream)
       !call CheckVslError(errcode)
    enddo
    !$omp end parallel do


    ! Check gaussian random field
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
    do k=1,grf%Nm1d
       do j=1,grf%Nm1d
          do i=1,grf%Nm1d
             d    = grf%grf(i,j,k)
             davg = davg + d
             dsig = dsig + d**2
             dmax = max(dmax,d)
             dmin = min(dmin,d)
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Write to screen
    davg = davg/grf%Nmesh
    dsig = sqrt(dsig/grf%Nmesh - davg**2)
    write(*,*) 'GRF : ',real((/davg,dsig,dmin,dmax/))


    ! Forward FFT gaussian random field
    ! See fft.f90
    call fft_3d(grf%grf,'f')


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : GRF calc'
    return
  end subroutine grf_calc


end module grfmake_module
