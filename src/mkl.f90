! MKL
include 'mkl_df.f90'
include 'mkl_dfti.f90'
include 'mkl_vsl.f90'


module mkl_module
  ! Intel
  use IFPORT
  use MKL_DF
  use MKL_DFTI
  use MKL_VSL
  use OMP_LIB
  use, intrinsic :: iso_c_binding


  ! Modules
  use timing_module


  ! Default
  implicit none
  public


contains


!------------------------------------------------------------------------------!
! FFT
!------------------------------------------------------------------------------!
  

  subroutine fft_3d(a,c)
    ! Default
    implicit none


    ! Subroutine arguments
    character  :: c
    real(8), dimension(:,:,:) :: a


    ! Local variables
    integer(4)  :: n1,n2,n3,status
    integer(8)  :: Nfft
    real(8)     :: bscale
    type(C_PTR) :: fft_ptr
    integer, dimension(3) :: Lfft
    integer, dimension(4) :: strides_in,strides_out
    type(DFTI_DESCRIPTOR), pointer :: desc
    real(8), pointer, dimension(:) :: a_fft


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Size of array
    Lfft    = shape(a)
    n1      = Lfft(1) - 2
    n2      = Lfft(2)
    n3      = Lfft(3)
    Lfft(1) = n1


    ! Create pointer
    Nfft    = int(n1+2,kind=8)*int(n2,kind=8)*int(n3,kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr,a_fft,[Nfft])
    

    ! Choose direction
    select case (c)
    case ('f')
       ! Init
       strides_in  = (/0,1,n1+2,(n1+2)*n2/)
       strides_out = (/0,1,n1/2+1,(n1/2+1)*n2/)

       ! Init descriptor
       status = DftiCreateDescriptor(desc,DFTI_DOUBLE,DFTI_REAL,3,Lfft)
       status = DftiSetValue(desc,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)
       status = DftiSetValue(desc,DFTI_INPUT_STRIDES ,strides_in)
       status = DftiSetValue(desc,DFTI_OUTPUT_STRIDES,strides_out)
       status = DftiCommitDescriptor(desc)

       ! Compute transform
       status = DftiComputeForward(desc,a_fft)
       !write(*,*) DftiErrorMessage(status)

       ! Free descriptor
       !status = DftiFreeDescriptor(desc)
    case ('b')
       ! Init
       strides_in  = (/0,1,n1/2+1,(n1/2+1)*n2/)
       strides_out = (/0,1,n1+2,(n1+2)*n2/)
       bscale      = 1./(real(n1)*real(n2)*real(n3))

       ! Init descriptor
       status = DftiCreateDescriptor(desc,DFTI_DOUBLE,DFTI_REAL,3,Lfft)
       status = DftiSetValue(desc,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)
       status = DftiSetValue(desc,DFTI_INPUT_STRIDES ,strides_in)
       status = DftiSetValue(desc,DFTI_OUTPUT_STRIDES,strides_out)
       status = DftiSetValue(desc,DFTI_BACKWARD_SCALE,bscale)
       status = DftiCommitDescriptor(desc)

       ! Compute transform
       status = DftiComputeBackward(desc,a_fft)
       !write(*,*) DftiErrorMessage(status)

       ! Free descriptor
       !status = DftiFreeDescriptor(desc)
    end select


    time2 = time()
    !write(*,'(2a)') timing(time1,time2),' : MKL FFT'
    return
  end subroutine fft_3d


!------------------------------------------------------------------------------!
! Spline
!------------------------------------------------------------------------------!  

  subroutine spline_cubic(x,y,xs,ys)
    ! Default
    implicit none


    ! Subroutine arguments
    real(8), dimension(:) :: x,y,xs,ys


    ! Local variables
    integer(4)    :: status
    integer(4)    :: nx,xhint,ny,yhint
    integer(4)    :: sorder,stype,bctype,scoeffhint
    integer(4)    :: sformat,method,type,nsite,shint,rhint
    type(df_task) :: task
    real(8), allocatable, dimension(:) :: scoeff
    

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()

    
    ! Init task
    nx     = size(x)
    xhint  = DF_NO_HINT
    ny     = 1
    yhint  = DF_NO_HINT
    status = dfdnewtask1d(task,nx,x,xhint,ny,y,yhint)
    
    
    ! Edit task
    sorder     = DF_PP_CUBIC
    stype      = DF_PP_NATURAL
    bctype     = DF_BC_FREE_END
    scoeffhint = DF_NO_HINT
    allocate(scoeff(sorder*ny*(nx-1)))
    status     = dfdeditppspline1d(task, &
         sorder,stype,bctype,scoeff=scoeff,scoeffhint=scoeffhint)
 

    ! Construct task
    sformat = DF_PP_SPLINE
    method  = DF_METHOD_STD
    status  = dfdconstruct1d(task,sformat,method)


    ! Perform task
    type   = DF_INTERP
    method = DF_METHOD_PP
    nsite  = size(xs)
    shint  = DF_NO_HINT
    rhint  = DF_NO_HINT
    status = dfdinterpolate1d(task, &
         type,method,nsite,xs,shint,r=ys,rhint=rhint)


    ! Destroy task
    status = dfdeletetask(task)


    time2 = time()
    !write(*,'(2a)') timing(time1,time2),' : MKL spline cubic'
    return
  end subroutine spline_cubic


  subroutine spline_construct(task,x,y)
    ! Default
    implicit none


    ! Subroutine arguments
    type(df_task)         :: task
    real(8), dimension(:) :: x,y


    ! Local variables
    integer(4) :: status
    integer(4) :: nx,xhint,ny,yhint
    integer(4) :: sorder,stype,bctype,scoeffhint
    integer(4) :: sformat,method
    real(8), allocatable, dimension(:) :: scoeff
    

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()

    
    ! Init task
    nx     = size(x)
    xhint  = DF_NO_HINT
    ny     = 1
    yhint  = DF_NO_HINT
    status = dfdnewtask1d(task,nx,x,xhint,ny,y,yhint)
    
    
    ! Edit task
    sorder     = DF_PP_CUBIC
    stype      = DF_PP_NATURAL
    bctype     = DF_BC_FREE_END
    scoeffhint = DF_NO_HINT
    allocate(scoeff(sorder*ny*(nx-1)))
    status     = dfdeditppspline1d(task, &
         sorder,stype,bctype,scoeff=scoeff,scoeffhint=scoeffhint)
 

    ! Construct task
    sformat = DF_PP_SPLINE
    method  = DF_METHOD_STD
    status  = dfdconstruct1d(task,sformat,method)

    
    time2 = time()
    !write(*,'(2a)') timing(time1,time2),' : MKL spline construct'
    return
  end subroutine spline_construct


  function spline_interp(task,xs)
    ! Default
    implicit none


    ! Function arguments
    type(df_task) :: task
    real(8), dimension(:) :: xs
    real(8), allocatable, dimension(:) :: spline_interp


    ! Local variables
    integer(4) :: status,sformat,method,type,nsite,shint,rhint


    ! Allocate
    nsite = size(xs)
    allocate(spline_interp(nsite))

    
    ! Perform task
    type   = DF_INTERP
    method = DF_METHOD_PP
    shint  = DF_NO_HINT
    rhint  = DF_NO_HINT
    status = dfdinterpolate1d(task, &
         type,method,nsite,xs,shint,r=spline_interp,rhint=rhint)


    return
  end function spline_interp


end module mkl_module
