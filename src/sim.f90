module sim_module
  ! Intel
  use OMP_LIB


  ! Modules
  use timing_module
  use input_module, only : input


  ! Default
  implicit none
  public


  ! Types
  type sim_type
     ! Input
     integer(4)    :: Nproc,Nm1d
     real(8)       :: Lbox
     character(80) :: dirin,dirout
     ! Variables
     integer(8)    :: Nmesh
     ! IO
     character(10) :: Lstr,Nstr,zstr
     character(80) :: fstr
  end type sim_type


  type domain_type
     ! Variables
     integer(4) :: Nd,Nx,Ny,Nz
     ! Arrays
     integer(4), allocatable, dimension(:,:)   :: D
     integer(4), allocatable, dimension(:,:,:) :: i
  end type domain_type

  
  type unit_type
     ! Conversions
     real(8) :: len,rho,time
     real(8) :: mass,vel,acc,pres,temp
     real(8) :: mdm,mgas,Msun,Msunh
     real(8) :: mesh_to_box,box_to_mesh
  end type unit_type


  ! Objects
  type(sim_type)    :: sim
  type(unit_type)   :: unit
  type(domain_type) :: domain
  

contains


  subroutine sim_init
    ! Default
    implicit none


    ! Local variables


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init from input
    sim%Nproc  = input%sim_Nproc
    sim%dirin  = input%sim_dirin
    sim%dirout = input%sim_dirout
    sim%Lbox   = input%cosmo_Lbox
    sim%Nm1d   = input%mesh_Nm1d


    ! Domain decomposition
    call domain_init


    ! Unit conversions
    call unit_init
    
    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : SIM init'
    return
  end subroutine sim_init


  subroutine domain_init
    ! Default
    implicit none


    ! Local variables
    integer(4) :: iproc
    integer(4) :: a,b,c,i,j,k,l,m,n
    integer(4) :: i1,j1,k1,n1,n2,n3


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Calc domain subdivision
    n = nint(sim%Nproc**(1./3))

    do i=n,sim%Nm1d
       if (mod(sim%Nm1d,2*i) == 0) then
          n1          = 2*i
          n2          = n1
          n3          = n1
          domain%Nx   = sim%Nm1d/n1
          domain%Ny   = domain%Nx
          domain%Nz   = domain%Nx
          domain%Nd = n1*n2*n3
          exit
       endif
    enddo


    ! Allocate domain arrays
    allocate(domain%D(27 ,domain%Nd))
    allocate(domain%i(2,3,domain%Nd))


    ! Calculate domain information
    do n=1,domain%Nd
       i = 1 + mod((n-1),n1)
       j = 1 + mod((n-1)/n1,n2)
       k = 1 +    ((n-1)/n1/n2)

       ! Calculate indices
       domain%i(1,1,n) = 1 + (i-1)*domain%Nx
       domain%i(2,1,n) = i*domain%Nx
       domain%i(1,2,n) = 1 + (j-1)*domain%Ny
       domain%i(2,2,n) = j*domain%Ny
       domain%i(1,3,n) = 1 + (k-1)*domain%Nz
       domain%i(2,3,n) = k*domain%Nz


       ! Write domain info
       if (.false.) then
          write(*,*) n,i,domain%i(:,1,n)
          write(*,*) n,j,domain%i(:,2,n)
          write(*,*) n,k,domain%i(:,3,n)
       endif


       ! Init domain
       m             = 1
       domain%D(1,n) = 0

       ! Calculate neighbors
       do c=-1,1
          k1 = 1 + mod(k + c - 1 + n3,n3)
          do b=-1,1
             j1 = 1 + mod(j + b - 1 + n2,n2)
             do a=-1,1
                i1 = 1 + mod(i + a - 1 + n1,n1)
                l  = i1 + (j1-1)*n1 + (k1-1)*n1*n2
              
                if (l /= n) then
                   m             = m + 1
                   domain%D(m,n) = l
                endif
             enddo
          enddo
       enddo
    enddo


    ! Write to screen
    write(*,*) 'Domain : ',domain%Nd,n1,domain%Nx    


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : SIM domain init'
    return
  end subroutine domain_init


  subroutine unit_init
    ! Default
    implicit none


    ! Local variables


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Conversions
    unit%mesh_to_box = sim%Lbox/sim%Nm1d
    unit%box_to_mesh = 1/unit%mesh_to_box


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : SIM unit init'
    return
  end subroutine unit_init

  
end module sim_module
