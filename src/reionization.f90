module reionization_module
  ! Intel
  use IFPORT
  use OMP_LIB
  use, intrinsic :: iso_c_binding


  ! Modules
  use reion_module
  use constant_module
  use helper_module
  use mkl_module
  use timing_module
  use cosmo_module     , only : cosmo
  use cosmology_module , only : cosmo_calc,cosmo_powerspectrum
  use esf_module       , only : esf
  use mesh_module      , only : mesh
  use meshmake_module  , only : mesh_density,mesh_velocity
  use sim_module       , only : sim,unit
  use simulation_module, only : sim_calc


  ! Default
  implicit none
  public


contains


  subroutine reion_make
    ! Default
    implicit none


    ! Local variables
    type(C_ptr) :: ptr

    
    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Pointers
    reion%rad => reion%zre
    reion%d   => mesh%fft1
    reion%m   => mesh%fft2
    reion%r   => mesh%fft1
    reion%s   => mesh%fft1
    reion%w   => mesh%fft2


    ! I8 array reion%isort with R8 array mesh%fft3
    ptr = c_loc(mesh%fft3)
    call c_f_pointer(ptr,reion%isort,[reion%Nmesh])
    
    
    ! Reionization midpoint
    cosmo%z = reion%zmid
    cosmo%a = 1/(1 + cosmo%z)
    call cosmo_calc
    call sim_calc


    ! Make
    if (reion%make == 'read') then
       ! IO
       call reion_read
    else
       ! Radiation field
       call reion_radiation
       call reion_sort

       ! Density field
       ! See mesh.f90
       call mesh_density
       
       ! Abundance matching for zre
       call reion_redshift

       ! Ionization fraction
       call reion_ionfrac

       ! IO
       if (reion%make == 'write') call reion_write
    endif


    ! Power spectrum
    ! See cosmology.f90
    call cosmo_powerspectrum


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION make'
    return
  end subroutine reion_make


!------------------------------------------------------------------------------!
! Radiation Field
!------------------------------------------------------------------------------


  subroutine reion_radiation
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k
    real(8)    :: w,wsum
    real(8)    :: r,ravg,rsig,rmax,rmin
    

    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init
    wsum = 0
    ravg = 0
    rsig = 0
    rmax = 0
    rmin = huge(rmin)


    ! Source field and convolution window function
    !$omp parallel         &
    !$omp default(shared)  &
    !$omp private(i,j,k,w) &
    !$omp reduction(+:wsum)
    !$omp do
    do k=1,reion%Nm1d
       do j=1,reion%Nm1d
          reion%s(1:reion%Nm1d,j,k) = esf%rho1(:,j,k)
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=1,reion%Nm1d
       do j=1,reion%Nm1d
          do i=1,reion%Nm1d
             ! Kernel
             w    = radiation_kernel(i,j,k)
             wsum = wsum + w
             reion%w(i,j,k) = w
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel


    ! Forward FFT fields
    ! See fft.f90
    call fft_3d(reion%s,'f')
    call fft_3d(reion%w,'f')


    ! Convolution in Fourier space
    ! reion%r, reion%s share mesh%fft1
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k,w)
    do k=1,reion%Nm1d
       do j=1,reion%Nm1d
          do i=1,reion%Nm1d+2,2
             w = reion%w(i,j,k)/wsum
             reion%r(i:i+1,j,k) = reion%s(i:i+1,j,k)*w
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Inverse FFT radiation field
    ! See fft.f90
    call fft_3d(reion%r,'b')


    ! Stats
    !$omp parallel do            &
    !$omp default(shared)        &
    !$omp private(i,j,k,r)       &
    !$omp reduction(+:ravg,rsig) &
    !$omp reduction(min:rmin)    &
    !$omp reduction(max:rmax)    
    do k=1,reion%Nm1d
       do j=1,reion%Nm1d
          do i=1,reion%Nm1d
             ! Save
             r = reion%r(i,j,k)
             reion%rad(i,j,k) = r

             ! Stats
             ravg = ravg + r
             rsig = rsig + r**2
             rmax = max(rmax,r)
             rmin = min(rmin,r)
          enddo
       enddo
    enddo
    !$omp end parallel do   


    ! Write to screen
    ravg = ravg/reion%Nmesh
    rsig = sqrt(rsig/reion%Nmesh - ravg**2)
    write(*,*) 'Radiation : ',real((/ravg,rsig,rmin,rmax/))
    

    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION radiation field'
    return


  contains


    function radiation_kernel(i,j,k)
      ! Default
      implicit none


      ! Function arguments
      integer(4) :: i,j,k
      real(8)    :: radiation_kernel


      ! Local parameters
      integer(4) :: Ndiv = 4


      ! Local arguments
      integer(4) :: l,m,n
      real(8)    :: mfp,r,rsq,w
      real(8)    :: x,y,z,x1,y1,z1


      ! Left boundaries of cell
      if (i <= (reion%Nm1d/2+1)) then
         x1 = i - 1.5
      else
         x1 = i - 1.5 - reion%Nm1d
      endif
      if (j <= (reion%Nm1d/2+1)) then
         y1 = j - 1.5
      else
         y1 = j - 1.5 - reion%Nm1d
      endif
      if (k <= (reion%Nm1d/2+1)) then
         z1 = k - 1.5
      else
         z1 = k - 1.5 - reion%Nm1d
      endif


      ! Mean free path in mesh units
      mfp = reion%mfp*unit%box_to_mesh


      ! Subdivide cell
      w = 0 
      do n=1,Ndiv
         z = z1 + (n - 0.5)/Ndiv
         do m=1,Ndiv
            y = y1 + (m - 0.5)/Ndiv
            do l=1,Ndiv
               x   = x1 + (l - 0.5)/Ndiv
               rsq = x**2 + y**2 + z**2
               r   = sqrt(rsq)
               w   = w + exp(-r/mfp)/rsq
            enddo
         enddo
      enddo

      ! Save
      radiation_kernel = w/Ndiv**3

      
      return
    end function radiation_kernel


  end subroutine reion_radiation


  subroutine reion_sort
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k
    integer(8) :: i8,Nm1d


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Init
    Nm1d = reion%Nm1d
    

    ! Indexing
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k,i8)
    do k=1,reion%Nm1d
       do j=1,reion%Nm1d
          do i=1,reion%Nm1d
             i8 = i + (j-1)*Nm1d + (k-1)*Nm1d**2
             reion%isort(i8) = i8
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Parallel index quicksort
    ! See helper.f90

    ! Recursive quicksort
    !$omp parallel      &
    !$omp default(shared)
    !$omp single
    call quicksort('d',1,reion%Nmesh,reion%Nmesh,reion%rad,reion%isort)
    !$omp end single
    !$omp end parallel


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION sort'
    return
  end subroutine reion_sort


!------------------------------------------------------------------------------!
! Reionization-redshift field
!------------------------------------------------------------------------------!


  subroutine reion_redshift
    ! Default
    implicit none


    ! Local variables
    integer(4) :: iproc
    real(8), dimension(5) :: xi


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()

    
    ! Init
    xi = 0
    

    ! Required for mass-weighted ionization fraction
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(iproc)
    do iproc=1,reion%Nproc
       call mass_fraction(iproc)
    enddo
    !$omp end parallel do


    ! Abundance match reionization history
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(iproc)  &
    !$omp reduction(max:xi)
    do iproc=1,reion%Nproc
       call abund_match(iproc,xi)
    enddo
    !$omp end parallel do


    ! Write to screen
    write(*,*) 'Beginning : ',real((/reion%zbeg,reion%xbeg,xi(1)/))
    write(*,*) 'Early     : ',real((/reion%zear,reion%xear,xi(2)/))
    write(*,*) 'Midpoint  : ',real((/reion%zmid,reion%xmid,xi(3)/))
    write(*,*) 'Late      : ',real((/reion%zlat,reion%xlat,xi(4)/))
    write(*,*) 'End       : ',real((/reion%zend,reion%xend,xi(5)/))


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION redshift field'
    return


  contains


    subroutine mass_fraction(iproc)
      ! Default
      implicit none


      ! Subroutine arguments
      integer(4) :: iproc


      ! Local variables
      integer(4) :: i,j,k,ip
      integer(8) :: i8,n8
      real(8)    :: d,m


      ! Cumulative mass fraction
      ! Loop over Nmesh cells, skip buffer cells
      m = 0
      d = 0
   
      do n8=mesh%proc(1,iproc),mesh%proc(2,iproc)
         ! Indices
         i8 = reion%isort(n8)
         i  = 1 + mod((i8-1)           ,reion%Nm1d)
         j  = 1 + mod((i8-1)/reion%Nm1d,reion%Nm1d)
         k  = 1 +     (i8-1)/reion%Nm1d**2

         ! Sum and save mass and delta
         m = m +  mesh%rho1(i,j,k)     /reion%Nmesh
         d = d + (mesh%rho1(i,j,k) - 1)/reion%Nmesh
         reion%m(i,j,k) = m
         reion%d(i,j,k) = d
      enddo

      
      !$omp barrier

      ! Offsets for cumulative mass fraction
      m = 0
      d = 0
      do ip=1,iproc-1
         n8 = mesh%proc(2,ip)
         i8 = reion%isort(n8)
         i  = 1 + mod((i8-1)           ,reion%Nm1d)
         j  = 1 + mod((i8-1)/reion%Nm1d,reion%Nm1d)
         k  = 1 +     (i8-1)/reion%Nm1d**2
         m  = m + reion%m(i,j,k)
         d  = d + reion%d(i,j,k)
      enddo

      !$omp barrier

      
      if (iproc > 1) then
         do n8=mesh%proc(1,iproc),mesh%proc(2,iproc)
            ! Indices
            i8 = reion%isort(n8)
            i  = 1 + mod((i8-1)           ,reion%Nm1d)
            j  = 1 + mod((i8-1)/reion%Nm1d,reion%Nm1d)
            k  = 1 +     (i8-1)/reion%Nm1d**2

            ! Add offset
            reion%m(i,j,k) = reion%m(i,j,k) + m
            reion%d(i,j,k) = reion%d(i,j,k) + d
         enddo
      endif
      

      return
    end subroutine mass_fraction


    subroutine abund_match(iproc,xi)
      ! Default
      implicit none


      ! Subroutine arguments
      integer(4)            :: iproc
      real(8), dimension(5) :: xi


      ! Local variables
      integer(4) :: i,j,k,n
      integer(8) :: i8,n8
      real(8)    :: g,x,z

      
      ! Loop over sorted field in descending order
      do n8=mesh%proc(1,iproc),mesh%proc(2,iproc)
         ! Cell indices
         i8 = reion%isort(n8)
         i  = 1 + mod((i8-1)           ,reion%Nm1d)
         j  = 1 + mod((i8-1)/reion%Nm1d,reion%Nm1d)
         k  = 1 +     (i8-1)/reion%Nm1d**2

         ! Iterative guess
         g = 1
         do n=1,3
            x = (reion%m(i,j,k) - reion%d(i,j,k)) + g*reion%d(i,j,k)
            z = z_of_xi(x)
            g = (1 + cosmo%z)/(1 + z)
         enddo

         ! Save redshift
         reion%zre(i,j,k) = z

         ! Save xi
         if (z >= reion%zbeg) xi(1) = max(xi(1),x)
         if (z >= reion%zear) xi(2) = max(xi(2),x)
         if (z >= reion%zmid) xi(3) = max(xi(3),x)
         if (z >= reion%zlat) xi(4) = max(xi(4),x)
         if (z >= reion%zend) xi(5) = max(xi(5),x)
      enddo


      return
    end subroutine abund_match


  end subroutine reion_redshift


!------------------------------------------------------------------------------!
! Ionization fraction
!------------------------------------------------------------------------------!


  subroutine reion_ionfrac
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i,j,k
    real(8)    :: xim,xiv


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! Mass-weighted and volume-weighted ionization fractions
    xim = 0
    xiv = 0

    !$omp parallel do        &
    !$omp default(shared)    &
    !$omp private(i,j,k)     &
    !$omp reduction(+:xim,xiv)
    do k=1,reion%Nm1d
       do j=1,reion%Nm1d
          do i=1,reion%Nm1d
             ! Sum only ionized
             if (cosmo%z <= reion%zre(i,j,k)) then
                xim = xim + mesh%rho1(i,j,k)
                xiv = xiv + 1
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    ! Save and write to screen
    reion%xi  = xi_of_z(cosmo%z)
    reion%xim = xim/reion%Nmesh
    reion%xiv = xiv/reion%Nmesh
    write(*,*) 'xi : ',real((/reion%xi,reion%xim,reion%xiv/))

    
    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : REION ionization fraction'
    return
  end subroutine reion_ionfrac

  
end module reionization_module
