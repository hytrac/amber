module reion_module
  ! Intel
  use OMP_LIB


  ! Modules
  use constant_module
  use timing_module
  use cosmo_module, only : cosmo
  use input_module, only : input
  use sim_module  , only : sim


  ! Default
  implicit none
  public


  ! Types
  type reion_type
     ! Input
     integer(4)    :: Nproc,Nm1d
     real(8)       :: zmid,zdel,zasy
     real(8)       :: xmid,xear,xlat
     real(8)       :: Mmin,mfp
     character(10) :: make
     character(80) :: dirin,dirout
     ! Variables
     integer(8)    :: Nmesh
     real(8)       :: zbeg,zear,zlat,zend
     real(8)       :: xbeg,xend
     real(8)       :: awb,bwb,cwb
     real(8)       :: xi,xim,xiv
     ! Arrays
     real(4), allocatable, dimension(:,:,:) :: zre
     ! Pointers
     integer(8), pointer, dimension(:)     :: isort
     real(4)   , pointer, dimension(:,:,:) :: rad
     real(8)   , pointer, dimension(:,:,:) :: d,m,r,s,w
  end type reion_type


  ! Objects
  type(reion_type), target :: reion


contains


  subroutine reion_init
    ! Default
    implicit none


    ! Local variables
    integer(4) :: k
    

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from input
    reion%make   = input%reion_make
    reion%zmid   = input%reion_zmid
    reion%zdel   = input%reion_zdel
    reion%zasy   = input%reion_zasy
    reion%xmid   = input%reion_xmid
    reion%xear   = input%reion_xear
    reion%xlat   = input%reion_xlat
    reion%Mmin   = input%reion_Mmin
    reion%mfp    = input%reion_mfp
    reion%dirout = input%reion_dir
    reion%dirin  = input%sim_dirin
    reion%Nproc  = input%sim_Nproc
    reion%Nm1d   = input%mesh_Nm1d


    ! Number of cells/particles
    reion%Nmesh = int(reion%Nm1d,kind=8)**3


    ! Allocate arrays
    allocate(reion%zre(reion%Nm1d,reion%Nm1d,reion%Nm1d))


    ! First touch in parallel
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(k)
    !$omp do
    do k=1,reion%Nm1d
       reion%zre(:,:,k) = 0
    enddo
    !$omp end do
    !$omp end parallel


    ! Redshifts
    reion%zear = reion%zmid + reion%zdel*reion%zasy/(1 + reion%zasy)
    reion%zlat = reion%zear - reion%zdel


    ! Weibull function
    call weibull_init
              

    ! Beginning and end
    reion%xbeg = 1D-4
    reion%zbeg = z_of_xi(reion%xbeg)
    reion%xend = 1D0
    reion%zend = z_of_xi(reion%xend)


    ! Write to screen
    write(*,*) 'Beginning : ',real((/reion%zbeg,reion%xbeg/))
    write(*,*) 'Early     : ',real((/reion%zear,reion%xear/))
    write(*,*) 'Midpoint  : ',real((/reion%zmid,reion%xmid/))
    write(*,*) 'Late      : ',real((/reion%zlat,reion%xlat/))
    write(*,*) 'End       : ',real((/reion%zend,reion%xend/))


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION init'
    return
  end subroutine reion_init


  subroutine weibull_init
    ! Default
    implicit none


    ! Local arguments
    integer(4) :: i,Nmax
    real(8)    :: c,cold,c1,c2
    real(8)    :: f,fold,dfdc,r,delta


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Weibull function
    ! xi(z) = exp[-((z-a)/b)^c]


    ! c is related to asymmetry?
    ! Newton's method with finite difference
    i     = 0
    delta = 1.0
    c     = 1.0
    cold  = 0.99*c
    fold  = fofc(cold)

    do while (abs(delta) > 1E-4 .and. i < 100)
       ! Increment
       i = i + 1

       ! Finite-difference
       f    = fofc(c)
       dfdc = (f - fold)/(c - cold)

       ! Save
       cold = c
       fold = f

       ! Guess
       c     = c - f/dfdc
       delta = (c - cold)/c
    enddo
    reion%cwb = c


    ! a is equal to z100
    r         = (log(reion%xear)/log(reion%xlat))**(1./reion%cwb)
    reion%awb = (reion%zear - reion%zlat*r)/(1 - r) 


    ! b is related to duration
    reion%bwb = (reion%zmid - reion%awb)/(-dlog(reion%xmid))**(1./reion%cwb)


    ! Write to screen
    write(*,*) 'Weibull : ',real((/reion%awb,reion%bwb,reion%cwb/))


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION Weibull init'
    return


  contains


    function fofc(c)
      ! Default
      implicit none

      ! Function arguments
      real(8) :: c,fofc

      ! Local arguments
      real(8) :: a,r

      r    = (log(reion%xlat)/log(reion%xmid))**(1./c)
      a    = (reion%zlat - reion%zmid*r)/(1 - r)
      fofc = (log(reion%xear)/log(reion%xmid))**(1./c) &
           - (reion%zear - a)/(reion%zmid - a)

      return
    end function fofc


  end subroutine weibull_init


!------------------------------------------------------------------------------!
! IO
!------------------------------------------------------------------------------!


  subroutine reion_read
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: un
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Read zre field
    un = 11
    fn = trim(reion%dirin)//'zre_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Reading ',trim(fn)
    open(un,file=fn,form='binary')
    read(un) reion%zre
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION read'
    return
  end subroutine reion_read

  
  subroutine reion_write
    ! Default
    implicit none


    ! Local variables
    integer(4)    :: un
    character(80) :: fn


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Write zre field
    un = 11
    fn = trim(reion%dirout)//'zre_'//trim(sim%fstr)//'.dat'
    write(*,*) 'Writing ',trim(fn)
    open(un,file=fn,form='binary')
    write(un) reion%zre
    close(un)


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION write'
    return
  end subroutine reion_write
  
  
!------------------------------------------------------------------------------!
! Functions
!------------------------------------------------------------------------------!


  function xi_of_z(z)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: z
    real(8) :: xi_of_z


    ! Weibull function
    if (z > reion%awb) then
       xi_of_z = exp(-((z - reion%awb)/reion%bwb)**reion%cwb)
    else
       xi_of_z = 1d0
    endif


    return
  end function xi_of_z


  function z_of_xi(xi)
    ! Default
    implicit none


    ! Function arguments
    real(8) :: xi
    real(8) :: z_of_xi


    ! Inverted Weibull function
    if (xi < 1d0) then
       z_of_xi = reion%awb + reion%bwb*(-dlog(xi))**(1./reion%cwb)
    else
       z_of_xi = reion%awb
    endif


    return
  end function z_of_xi


end module reion_module
