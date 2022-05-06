module input_module
  ! Intel
  use OMP_LIB


  ! Modules
  use timing_module


  ! Default
  implicit none
  public


  ! Types
  type input_type
     ! Simulation
     integer(4)    :: sim_Nproc
     character(80) :: sim_dirin,sim_dirout
     ! Cosmo
     real(8)       :: cosmo_Lbox
     real(8)       :: cosmo_ob,cosmo_om,cosmo_ol,cosmo_or
     real(8)       :: cosmo_h,cosmo_s8,cosmo_ns,cosmo_w
     real(8)       :: cosmo_Tcmb0,cosmo_XH,cosmo_YHe
     character(80) :: cosmo_file,cosmo_dir
     ! Reionization
     character(10) :: reion_make
     real(8)       :: reion_zmid,reion_zdel,reion_zasy
     real(8)       :: reion_xmid,reion_xear,reion_xlat
     real(8)       :: reion_Mmin,reion_mfp
     character(80) :: reion_dir
     ! GRF
     character(10) :: grf_make
     integer(4)    :: grf_seed
     character(80) :: grf_dir
     ! LPT
     character(10) :: lpt_make,lpt_order,lpt_assign
     character(80) :: lpt_dir
     ! ESF
     character(10) :: esf_make,esf_filter,esf_assign
     character(80) :: esf_dir
     ! Mesh
     character(10) :: mesh_make
     integer(4)    :: mesh_Nm1d
     character(80) :: mesh_dir
     ! CMB
     character(10) :: cmb_make
     integer(4)    :: cmb_lmin,cmb_lmax
     real(8)       :: cmb_zmin,cmb_zmax,cmb_zdel
     character(10) :: cmb_zspacing
     character(80) :: cmb_dir
     ! H21cm
     character(10) :: h21cm_make
     real(8)       :: h21cm_zmin,h21cm_zmax,h21cm_zdel
     character(10) :: h21cm_zspacing
     character(80) :: h21cm_dir     
  end type input_type


  ! Objects
  type(input_type) :: input


contains


  subroutine input_read
    ! Default
    implicit none


    ! Timing variables
    integer(4) :: time1,time2
    time1 = time()


    ! ./amber.x < input.txt
    write(*,*) 'Reading input parameters'


    ! Simulation
    read( *,*)
    read( *,*)
    read( *,*) input%sim_Nproc
    read( *,*) input%sim_dirin
    read( *,*) input%sim_dirout
    write(*,*) 'SIMULATION'
    write(*,*) 'Nproc          = ',input%sim_Nproc
    write(*,*) 'dir in         = ',trim(input%sim_dirin)
    write(*,*) 'dir out        = ',trim(input%sim_dirout)


    ! Cosmo
    read( *,*)
    read( *,*)
    read( *,*) input%cosmo_Lbox
    read( *,*) input%cosmo_om
    read( *,*) input%cosmo_ol
    read( *,*) input%cosmo_ob
    read( *,*) input%cosmo_or
    read( *,*) input%cosmo_h
    read( *,*) input%cosmo_s8
    read( *,*) input%cosmo_ns
    read( *,*) input%cosmo_w
    read( *,*) input%cosmo_Tcmb0
    read( *,*) input%cosmo_XH
    read( *,*) input%cosmo_YHe
    read( *,*) input%cosmo_file
    read( *,*) input%cosmo_dir
    write(*,*) 'COSMO'
    write(*,*) 'Box length     = ',real(input%cosmo_Lbox)
    write(*,*) 'Omega_m        = ',real(input%cosmo_om)
    write(*,*) 'Omega_l        = ',real(input%cosmo_ol)
    write(*,*) 'Omega_b        = ',real(input%cosmo_ob)
    write(*,*) 'Omega_r        = ',real(input%cosmo_or)
    write(*,*) 'h0             = ',real(input%cosmo_h)
    write(*,*) 'sigma_8        = ',real(input%cosmo_s8)
    write(*,*) 'n_s            = ',real(input%cosmo_ns)
    write(*,*) 'w_de           = ',real(input%cosmo_w)
    write(*,*) 'T_cmb          = ',real(input%cosmo_Tcmb0)
    write(*,*) 'X_hydrogen     = ',real(input%cosmo_XH)
    write(*,*) 'Y_helium       = ',real(input%cosmo_YHe)
    write(*,*) 'Plin file      = ',trim(input%cosmo_file)
    write(*,*) 'Dir            = ',trim(input%cosmo_dir)

    
    ! Reionization
    read( *,*)
    read( *,*)
    read( *,*) input%reion_make
    read( *,*) input%reion_zmid
    read( *,*) input%reion_zdel
    read( *,*) input%reion_zasy
    read( *,*) input%reion_xear
    read( *,*) input%reion_xmid
    read( *,*) input%reion_xlat
    read( *,*) input%reion_Mmin
    read( *,*) input%reion_mfp
    read( *,*) input%reion_dir
    write(*,*) 'REIONIZATION'
    write(*,*) 'Make           = ',trim(input%reion_make)
    write(*,*) 'z mid          = ',real(input%reion_zmid)
    write(*,*) 'z delta        = ',real(input%reion_zdel)
    write(*,*) 'z asymmetry    = ',real(input%reion_zasy)
    write(*,*) 'xi early       = ',real(input%reion_xear)
    write(*,*) 'xi mid         = ',real(input%reion_xmid)
    write(*,*) 'xi late        = ',real(input%reion_xlat)
    write(*,*) 'Halo min mass  = ',real(input%reion_Mmin)
    write(*,*) 'mean free path = ',real(input%reion_mfp)
    write(*,*) 'Dir            = ',trim(input%reion_dir)


    ! GRF
    read( *,*)
    read( *,*)
    read( *,*) input%grf_make
    read( *,*) input%grf_seed
    read( *,*) input%grf_dir
    write(*,*) 'GRF'
    write(*,*) 'Make           = ',trim(input%grf_make)
    write(*,*) 'Seed           = ',input%grf_seed
    write(*,*) 'Dir            = ',trim(input%grf_dir)


    ! LPT
    read( *,*)
    read( *,*)
    read( *,*) input%lpt_make
    read( *,*) input%lpt_order
    read( *,*) input%lpt_assign
    read( *,*) input%lpt_dir
    write(*,*) 'LPT'
    write(*,*) 'Make           = ',trim(input%lpt_make)
    write(*,*) 'Order          = ',trim(input%lpt_order)
    write(*,*) 'Assignment     = ',trim(input%lpt_assign)
    write(*,*) 'Dir            = ',trim(input%lpt_dir)


    ! ESF
    read( *,*)
    read( *,*)
    read( *,*) input%esf_make
    read( *,*) input%esf_filter
    read( *,*) input%esf_assign
    read( *,*) input%esf_dir
    write(*,*) 'ESF'
    write(*,*) 'Make           = ',trim(input%esf_make)
    write(*,*) 'Filter         = ',trim(input%esf_filter)
    write(*,*) 'Assignment     = ',trim(input%esf_assign)
    write(*,*) 'Dir            = ',trim(input%esf_dir)


    ! Mesh
    read( *,*)
    read( *,*)
    read( *,*) input%mesh_make
    read( *,*) input%mesh_Nm1d
    read( *,*) input%mesh_dir
    write(*,*) 'MESH'
    write(*,*) 'Make           = ',trim(input%mesh_make)
    write(*,*) 'Nm1d           = ',input%mesh_Nm1d
    write(*,*) 'Dir            = ',trim(input%mesh_dir)


    ! CMB
    read( *,*)
    read( *,*)
    read( *,*) input%cmb_make
    read( *,*) input%cmb_zmin
    read( *,*) input%cmb_zmax
    read( *,*) input%cmb_zdel
    read( *,*) input%cmb_zspacing
    read( *,*) input%cmb_lmin
    read( *,*) input%cmb_lmax
    read( *,*) input%cmb_dir
    write(*,*) 'CMB'
    write(*,*) 'Make           = ',trim(input%cmb_make)
    write(*,*) 'z min          = ',real(input%cmb_zmin)
    write(*,*) 'z max          = ',real(input%cmb_zmax)
    write(*,*) 'z del          = ',real(input%cmb_zdel)
    write(*,*) 'z spacing      = ',trim(input%cmb_zspacing)
    write(*,*) 'l min          = ',input%cmb_lmin
    write(*,*) 'l max          = ',input%cmb_lmax
    write(*,*) 'Dir            = ',trim(input%cmb_dir)


    ! H21cm
    read( *,*)
    read( *,*)
    read( *,*) input%h21cm_make
    read( *,*) input%h21cm_zmin
    read( *,*) input%h21cm_zmax
    read( *,*) input%h21cm_zdel
    read( *,*) input%h21cm_zspacing
    read( *,*) input%h21cm_dir
    write(*,*) 'H21cm'
    write(*,*) 'Make           = ',trim(input%h21cm_make)
    write(*,*) 'z min          = ',real(input%h21cm_zmin)
    write(*,*) 'z max          = ',real(input%h21cm_zmax)
    write(*,*) 'z del          = ',real(input%h21cm_zdel)
    write(*,*) 'z spacing      = ',trim(input%h21cm_zspacing)
    write(*,*) 'Dir            = ',trim(input%h21cm_dir)

    
    ! Redefine directories
    input%sim_dirin  = trim(input%sim_dirin )//'/'
    input%sim_dirout = trim(input%sim_dirout)//'/'
    input%cmb_dir    = trim(input%sim_dirout)//trim(input%cmb_dir  )//'/'
    input%cosmo_dir  = trim(input%sim_dirout)//trim(input%cosmo_dir)//'/'
    input%esf_dir    = trim(input%sim_dirout)//trim(input%esf_dir  )//'/'
    input%grf_dir    = trim(input%sim_dirout)//trim(input%grf_dir  )//'/'
    input%h21cm_dir  = trim(input%sim_dirout)//trim(input%h21cm_dir)//'/'
    input%lpt_dir    = trim(input%sim_dirout)//trim(input%lpt_dir  )//'/'
    input%mesh_dir   = trim(input%sim_dirout)//trim(input%mesh_dir )//'/'
    input%reion_dir  = trim(input%sim_dirout)//trim(input%reion_dir)//'/'


    time2 = time()
    write(*,'(2a)') timing(time1,time2),' : Input read'
    return
  end subroutine input_read


end module input_module
