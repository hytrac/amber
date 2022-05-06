module timing_module
  ! Intel
  use IFPORT


  ! Default
  implicit none


contains


  function timing(t1,t2)
    ! Default
    implicit none

    
    ! Function arguments
    integer(4)    :: t1,t2
    character(37) :: timing


    ! Local variables
    integer(4)   :: s,m,h,d,dt
    character(3) :: sc,mc,hc,dc


    ! Calc time
    dt = t2-t1
    s  = mod(dt      ,60)
    m  = mod(dt/60   ,60)
    h  = mod(dt/3600 ,24)
    d  =     dt/86400


    ! Write strings
    if (s < 10) then
       write(sc,'(a,i1,a)') '0',s,'s'
    else
       write(sc,'(  i2,a)')     s,'s'
    endif
    if (m < 10) then
       write(mc,'(a,i1,a)') '0',m,'m'
    else
       write(mc,'(  i2,a)')     m,'m'
    endif
    if (h < 10) then
       write(hc,'(a,i1,a)') '0',h,'h'
    else
       write(hc,'(  i2,a)')     h,'h'
    endif
    if (d < 10) then
       write(dc,'(a,i1,a)') '0',d,'d'
    else
       write(dc,'(  i2,a)')     d,'d'
    endif

    timing = ctime(t2)//' '//dc//hc//mc//sc


    return
  end function timing


end module timing_module
