module simulation_module
  ! Intel
  use OMP_LIB


  ! Modules
  use constant_module
  use sim_module
  use timing_module
  use cosmo_module, only : cosmo


  ! Default
  implicit none
  public


contains


  subroutine sim_calc
    ! Default
    implicit none


    ! Local variables
    integer(4) :: L
    

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Box length
    L = nint(cosmo%Lbox)
    select case (L)
    case (1:9)
       write(sim%Lstr,'(a2,i1)') 'L=',L
    case (10:99)
       write(sim%Lstr,'(a2,i2)') 'L=',L
    case (100:999)
       write(sim%Lstr,'(a2,i3)') 'L=',L
    case (1000:9999)
       write(sim%Lstr,'(a2,i4)') 'L=',L
    end select


    ! Nm1d
    select case (sim%Nm1d)
    case (1:9)
       write(sim%Nstr,'(a2,i1)') 'N=',sim%Nm1d
    case (10:99)
       write(sim%Nstr,'(a2,i2)') 'N=',sim%Nm1d
    case (100:999)
       write(sim%Nstr,'(a2,i3)') 'N=',sim%Nm1d
    case (1000:9999)
       write(sim%Nstr,'(a2,i4)') 'N=',sim%Nm1d
    end select


    ! Redshift
    if (cosmo%z < 10) then
       write(sim%zstr,'(a3,f4.2)') 'z=0',cosmo%z
    else if (cosmo%z < 100) then
       write(sim%zstr,'(a2,f5.2)') 'z=' ,cosmo%z
    else if (cosmo%z < 1000) then
       write(sim%zstr,'(a2,f5.1)') 'z=' ,cosmo%z
    endif


    ! File extension
    sim%fstr = trim(sim%Lstr)//'_'//trim(sim%Nstr)//'_'//trim(sim%zstr)


    ! Unit conversions
    call unit_calc


    time2 = time()
    !write(*,'(2a)') timing(time1,time2),' : SIM calc'
    return
  end subroutine sim_calc

  
!------------------------------------------------------------------------------!
! Domain
!------------------------------------------------------------------------------!

  
  subroutine domain_set(i)
    ! Default
    implicit none


    ! Subroutine arguments
    integer(4) :: i


    ! Local variables
    integer(4) :: l,m,n
    logical    :: loop,eval


    !$omp critical (access_domain)
    ! Repeat until eligible domain is found
    loop = .true.

    do while (loop)
       do n=1,domain%Nd
          if (domain%D(1,n) == 0) then
             eval = .true.
             do m=2,27
                l = domain%D(m,n)
                if (domain%D(1,l) == 1) then
                   eval = .false.
                   exit
                endif
             enddo

             if (eval) then
                i             = n
                domain%D(1,i) = 1
                loop          = .false.
                exit
             endif
          endif
       enddo
    enddo
    !$omp end critical (access_domain)


    return
  end subroutine domain_set


  subroutine domain_end(i)
    ! Default
    implicit none


    ! Subroutine arguments
    integer(4) :: i


    ! Domain evaluated
    domain%D(1,i) = 2


    ! Reset domains if all are evaluated
    if (sum(domain%D(1,:)) == 2*domain%Nd) domain%D(1,:) = 0


    return
  end subroutine domain_end


!------------------------------------------------------------------------------!
! Unit Conversion
!------------------------------------------------------------------------------!


  subroutine unit_calc
    ! Default
    implicit none


    ! Local variables
    real(8) :: a,z


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Cosmo
    a = cosmo%a
    z = cosmo%z


    ! Mesh/particle to physical values (cgs)
    unit%len  = (a*cosmo%Lbox/cosmo%h/sim%Nm1d)*Mpc2cm
    unit%rho  = rhoc_cgs*cosmo%om*cosmo%h**2/a**3
    unit%time = 2/(3*H0_cgs*cosmo%h*sqrt(cosmo%om))*a**2
    unit%mass = unit%rho*unit%len**3
    unit%vel  = unit%len/unit%time
    unit%acc  = unit%len/unit%time**2
    unit%temp = unit%vel**2
    unit%pres = unit%rho*unit%temp


    ! Mass
    unit%mdm   = cosmo%fc*unit%mass
    unit%mgas  = cosmo%fb*unit%mass
    unit%Msun  = unit%mass/Msun_cgs
    unit%Msunh = unit%Msun*cosmo%h


    time2 = time()
    !write(*,'(2a)') timing(time1,time2),' : UNIT calc'
    return
  end subroutine unit_calc

  
end module simulation_module
