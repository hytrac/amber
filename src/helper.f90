module helper_module
  ! Intel
  use OMP_LIB
  use IFPORT


  ! Modules
  use constant_module


  ! Default
  implicit none
  public


contains


!------------------------------------------------------------------------------!
! Subroutines
!------------------------------------------------------------------------------!

  
  recursive subroutine quicksort(C,Ndiv,Nr,Ni,r4,i8)
    ! Default
    implicit none


    ! Subroutine arguments
    character  :: C
    integer(4) :: Ndiv
    integer(8) :: Nr,Ni
    real(4),    dimension(Nr) :: r4
    integer(8), dimension(Ni) :: i8      


    ! Local variables
    integer(4) :: Nd,Ndmax,Nimin
    integer(8) :: ip,N1,N2


    ! Maximum number of divisions
    Ndmax = 8*omp_get_max_threads()


    ! Mininum number of samples for recursion
    Nimin = 1024**2


    ! (Partition and continue recursion) or (qsort and end recursion)
    if (Ndiv < Ndmax .and. Ni > Nimin) then
       ! Partition
       ! ip is pivot index
       call partition(C,Nr,Ni,r4,i8,ip)

       ! Recursion
       Nd = 2*Ndiv
       N1 = ip - 1
       N2 = Ni - ip

       !$omp task          &
       !$omp default(shared)
       call quicksort(C,Nd,Nr,N1,r4,i8(:ip-1))
       !$omp end task

       !$omp task          &
       !$omp default(shared)
       call quicksort(C,Nd,Nr,N2,r4,i8(ip+1:))
       !$omp end task
       
       !$omp taskwait
    else
       ! Call IFPORT qsort         
       select case (C)
       case ('a')
          call qsort(i8,Ni,int(8,kind=8),ascending)
       case ('d')
          call qsort(i8,Ni,int(8,kind=8),descending)
       end select
    endif


    return


  contains


    subroutine partition(C,Nr,Ni,r4,i8,ip)
      ! Default
      implicit none


      ! Function arguments
      character  :: C
      integer(8) :: Nr,Ni,ip
      real(4),    dimension(Nr) :: r4
      integer(8), dimension(Ni) :: i8


      ! Local variables
      integer(8) :: i,j,n
      real(4)    :: pivot

      
      ! Init
      pivot = r4(i8(Ni))
      i     = 0

      ! Loop over i8
      select case (C)
      case ('a')
         do j=1,Ni-1
            if (r4(i8(j)) < pivot) then
               i     = i + 1
               n     = i8(i)
               i8(i) = i8(j)
               i8(j) = n
            endif
         enddo
      case ('d')
         do j=1,Ni-1
            if (r4(i8(j)) > pivot) then
               i     = i + 1
               n     = i8(i)
               i8(i) = i8(j)
               i8(j) = n
            endif
         enddo         
      end select


      ! Pivot
      ip     = i + 1
      n      = i8(ip)
      i8(ip) = i8(Ni)
      i8(Ni) = n


      return
    end subroutine partition


    function ascending(i1,i2)
      ! Function arguments
      integer(8) :: i1,i2
      integer(2) :: ascending

      if (r4(i1) <= r4(i2)) then
         ascending = -1
      else
         ascending =  1
      endif

      return
    end function ascending


    function descending(i1,i2)
      ! Function arguments
      integer(8) :: i1,i2
      integer(2) :: descending

      if (r4(i1) <= r4(i2)) then
         descending =  1
      else
         descending = -1
      endif

      return
    end function descending


  end subroutine quicksort


!------------------------------------------------------------------------------!
! Functions
!------------------------------------------------------------------------------!


  function assignment_transform(kx,ky,kz,scheme)
    ! Function arguments
    character(*) :: scheme
    real(8)      :: kx,ky,kz
    real(8)      :: assignment_transform

    ! Local variables
    integer(4) :: p

    select case (scheme)
    case ('ngp')
       p = 1
    case ('cic')
       p = 2
    case ('tsc')
       p = 3
    case default
       p = 0
    end select

    assignment_transform = (sinc(kx/2)*sinc(ky/2)*sinc(kz/2))**p

    return
  end function assignment_transform


  function sharpk_filter(x)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: x,sharpk_filter


    ! Local variables
    real(8) :: xk


    xk = (9*pi/2)**(1./3)

    if (x < xk) then
       sharpk_filter = 1
    else
       sharpk_filter = 0
    endif


    return
  end function sharpk_filter


  function sinc(x)
    ! Function arguments
    real(8) :: x,sinc


    if (abs(x) > 1E-6) then
       sinc = sin(x)/x
    else
       sinc = 1 - x**2/6
    endif


    return
  end function sinc


  function tophat_transform(x)
    ! Function arguments
    real(8) :: x,tophat_transform


    if (abs(x) > 1E-6) then
       tophat_transform = 3*(sin(x)-cos(x)*x)/x**3
    else
       tophat_transform = 1 - x**2/10
    endif


    return
  end function tophat_transform


end module helper_module
