! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_cice

use        types_mod, only : r8
!use    utilities_mod, only : initialize_utilities, finalize_utilities, &
!                             find_namelist_in_file, check_namelist_read, &
!                             file_exist, error_handler, E_ERR, E_MSG, to_upper
!use netcdf

implicit none

character(len=256) :: dart_to_cice_input_file = 'dart_restart.nc'
character(len=256) :: original_cice_input_file = 'restart_state.nc'

real(r8), allocatable :: aicen_org(:),vicen_org(:),vsnon_org(:)
real(r8), allocatable :: aicen(:),vicen(:),vsnon(:)
real(r8), allocatable :: Tsfcn(:)
real(r8), allocatable :: qice(:,:),sice(:,:),qsno(:,:)
real(r8), allocatable :: aice



end program dart_to_cice



