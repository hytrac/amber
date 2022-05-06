! AMBER  :  Abundance Matching Box for the Epoch of Reionization
! Compile:  make clean; make
! Execute:  ./amber.x < input/input.txt > output/log


program main
  ! Intel
  use OMP_LIB


  ! Modules
  use cmb_module
  use cmbreion_module
  use constant_module
  use cosmo_module
  use cosmology_module
  use esf_module
  use esfhalo_module
  use grf_module
  use grfmake_module
  use h21cm_module
  use h21cmreion_module
  use helper_module
  use input_module
  use lpt_module
  use lptmake_module
  use mesh_module
  use meshmake_module
  use mkl_module
  use reion_module
  use reionization_module
  use sim_module
  use simulation_module
  use timing_module


  ! Default
  implicit none


  ! Timing variables
  integer(4) :: time1,time2
  time1 = time()


  ! Read input params
  ! See input.f90
  call input_read


  ! Init Openmp
  call OMP_SET_NUM_THREADS(input%sim_Nproc)


  ! Init
  ! See corresponding .f90 files
  call sim_init
  call cosmo_init
  call mesh_init


  ! Init
  cosmo%z = 0
  cosmo%a = 1/(1 + cosmo%z)
  call cosmo_calc
  call sim_calc

  
  ! Cosmo
  ! See cosmo.f90, cosmology.f90
  call cosmo_linearpowerspectrum

  
  ! GRF
  ! See grf.f90, grfmake.f90
  call grf_init
  call grf_make


  ! LPT
  ! See lpt.f90, lptmake.f90
  call lpt_init
  call lpt_make


  ! ESF
  ! See esf.f90, esfhalo.f90
  call esf_init
  call esf_make

  
  ! Reionization
  ! See reion.f90, reionization.f90
  call reion_init
  call reion_make


  ! CMB
  ! See cmb.f90, cmbreion.f90
  call cmb_init
  call cmb_make


  ! Hydrogen 21cm
  ! See h21cm.f90, h21cmreion.f90
  call h21cm_init
  call h21cm_make  


  time2 = time()
  write(*,'(2a)') timing(time1,time2),' : AMBER completed'


end program main



