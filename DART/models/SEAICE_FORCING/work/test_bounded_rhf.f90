program test_bounded_rhf

! Test specific cases for obs_increment_bounded_norm_rhf
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

use     obs_utilities_mod, only : add_obs_to_seq, create_3d_obs, &
                                  getdimlen, getvar_int, set_missing_name,&
                                  getvar_real

use adaptive_inflate_mod, only : update_inflation, inflate_ens

use types_mod,         only : r8,metadatalength
use assim_tools_mod,   only : obs_increment_bounded_norm_rhf, &
                              obs_increment_rank_histogram, &
                              obs_increment_eakf, &
                              get_truncated_normal_like, &
                              update_from_obs_inc, &
                              assim_tools_init
                              
use random_seq_mod,    only : random_seq_type, init_random_seq, random_gaussian,random_beta, &
                              random_trunc_gaussian

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, nc_check, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode
use         utilities_mod, only : register_module, error_handler,               &
                                  E_ERR, E_MSG, nmlfileunit, get_unit,  &
                                  do_output, to_upper, logfileunit,   &
                                  find_namelist_in_file, check_namelist_read,   &
                                  file_exist, find_textfile_dims, file_to_text, &
                                  do_nml_file, do_nml_term
use default_model_mod,     only : init_time, init_conditions, adv_1step, &
                                  nc_write_model_vars

use beta_dist_functions,   only : betacdf,betaincinv

use netcdf
implicit none

real(r8), parameter :: truth = 0.9899_r8
integer,  parameter :: ens_size = 99
integer,  parameter :: num_steps = 5000
integer,  parameter :: bounds_case = 4         ! Case 4 is doubly bounded
real(r8), parameter :: bound(2) = (/0.0_r8, 1.0_r8/)
logical             :: is_bounded(4, 2)

character(len=25) :: in_file  = 'seaice_snow_forcings.nc'
character(len=19) :: out_file_control  = 'control_output.nc'
character(len=50) :: out_file_assim  = 'assim_output.nc'

integer :: ncid,iunit,io
integer :: filter_type = 1
integer :: obs_err_method = 1
integer :: obs_err_dist = 1
real(r8) :: const_err_val = 0.01_r8
logical :: inflate_prior = .true.
integer :: inflation_type = 1
real(r8) :: inflate_val = 2.0_r8
real(r8), allocatable :: control(:,:),ICs(:),input(:,:),preassim(:,:),output(:,:),post(:),forcings(:,:),prior(:)
real(r8), allocatable :: control_sic(:,:),control_sit(:,:),control_snow(:,:)
real(r8), allocatable :: prior_sic(:),prior_sit(:),prior_snow(:)
real(r8), allocatable :: post_sic(:),post_sit(:),post_snow(:)
real(r8), allocatable :: ICs_sic(:),ICs_sit(:),ICs_snow(:)
real(r8), allocatable :: forcings_sic(:,:),forcings_sit(:,:),forcings_snow(:,:)
real(r8), allocatable :: input_sic(:,:),input_sit(:,:),input_snow(:,:)
real(r8), allocatable :: preassim_sic(:,:),preassim_sit(:,:),preassim_snow(:,:)
real(r8), allocatable :: output_sic(:,:),output_sit(:,:),output_snow(:,:)

real(r8), allocatable :: assim(:,:),obs(:),obs_var(:),obs_inc(:)
real(r8), allocatable :: obs_sic(:),obs_sit(:),obs_snow(:)
real(r8), allocatable :: obsvar_sic(:),obsvar_sit(:),obsvar_snow(:)
real(r8), allocatable :: hold_FO_sic(:,:),hold_FO_sit(:,:),hold_FO_snow(:,:)

real(r8), allocatable :: truth_mem_sic(:), truth_mem_sit(:), truth_mem_snow(:),state_inc(:)

integer,parameter       :: n_sic = 35
real(r8), allocatable :: sicvals(:),sicerr(:),coefsicb(:),coefsicc(:),coefsicd(:)


real(r8)              :: a,sd_inflate,sic_err,hold_obs
real(r8)              :: varying_ss_inflate(99), varying_ss_inflate_sd(99),likelihood(ens_size),ens(ens_size),FO_SIT(ens_size),FO_SIC(ens_size),FO_SNOW(ens_size)
real(r8)              :: obs_std, prior_mean, prior_var,obs_place,initial_inflate,initial_inflate_sd
real(r8)              :: prior2(99),prior_var2,like_sum,ratio,reg_coef
type(random_seq_type) :: r_seq
integer               :: i, j,nmembers,ntimes,varid,nlocs,t,dt,fid,rc,nTimesDimID,nMembersDimID,truth_mem_ind
integer               :: save_index
character(len=512) :: msgstring
character(len=metadatalength) :: model_state_vars(3 * 1 ) = ' '
character(len=metadatalength) :: obs_kinds(3 * 1 ) = ' '
character(len=NF90_MAX_NAME), allocatable :: state_var_table(:),obs_kind_table(:)
character(len=NF90_MAX_NAME) :: varname

namelist /seaice_nml/  &
   filter_type, &
   obs_err_method, &
   obs_err_dist, &
   out_file_control, &
   out_file_assim, &
   const_err_val, &
   inflate_prior, &
   inflation_type, &
   inflate_val, &
   model_state_vars, &
   obs_kinds


!!!! INPUT MODS 
!filter_type = 1
!obs_err_method = 1
!obs_err_dist = 1

is_bounded(1:4, 1) = (/.false., .true., .false., .true./)
is_bounded(1:4, 2) = (/.false., .false., .true., .true./)

!initial_inflate = 1.0
!initial_inflate_sd = 0.6
!inflate_prior = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Initialize for error handler
call initialize_mpi_utilities('test_bounded_rhf')

call find_namelist_in_file('input.nml', 'seaice_nml', iunit)
read(iunit, nml = seaice_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')


if (model_state_vars(1) == ' ') then
  print*,'NO STATE MODEL VARS!'
  STOP
else
  do i =1, 3
    if (model_state_vars(i) == ' ') then
      save_index = i-1
      exit
    endif
  enddo
  allocate(state_var_table(save_index))
  do i = 1, 3
    if (model_state_vars(i) == ' ') exit
    varname = trim(model_state_vars(i))
    call to_upper(varname)
    state_var_table(i) = varname
  enddo
endif
if (obs_kinds(1) == ' ') then
  print*,'NO OBS TO BE ASSIMILATED'
  STOP
else
  do i =1, 3
    if (obs_kinds(i) == ' ') then
      save_index = i-1
      exit
    endif
  enddo
  allocate(obs_kind_table(save_index))
  do i = 1, 3
    if (obs_kinds(i) == ' ') exit
    varname = trim(obs_kinds(i))
    call to_upper(varname)
    obs_kind_table(i) = varname
  enddo
endif
print*,'Model State Variables'
print*,'----------------------'
do i=1,size(state_var_table)
  print*,state_var_table(i)
enddo
print*,'Assimilated Observation Types'
print*,'----------------------'
do i=1,size(obs_kind_table)
  print*,obs_kind_table(i)
enddo


print*,'-------------------------------------------------'
print*,'------------EXPERIMENT SETUP INFO----------------'
if (filter_type == 1) then
  print*,'FILTER TYPE: EaKF'
else if (filter_type == 2) then
  print*,'FILTER TYPE: RHF'
else if (filter_type == 3) then
  print*,'FILTER TYPE: Bounded RHF'
Else if (filter_type == 4) then
  print*,'FILTER TYPE: Particle Filter'
endif

if (obs_err_method == 1) then
  print*,'OBS ERR = 0.15*Truth'
else if (obs_err_method == 2) then
  print*,'OBS ERR = 0.15 - (0.99_r8*0.15)'
else if (obs_err_method == 3) then
  print*,'Beta maximum variance bell shape - sqrt((truth_mem_data(i)*(1.0_r8-truth_mem_data(i)))/ratio)'
else if (obs_err_method == 4) then
  print*,'Use spline function to fit sic errors - decrease with increasing sic'
else
  print*,'Use constant obs error', const_err_val
endif 

if (obs_err_dist == 1) then
  print*,'OBS ERR DIST is Gaussian'
else if (obs_err_dist == 2) then
  print*,'OBS ERR DIST is Beta'
endif

if (inflate_prior) then
  print*,'Prior inflation is turned on with an inflate val: ',inflate_val
  if (inflation_type == 1) then
    print*,'Prior inflation assuming Gaussian dist'
  else if (inflation_type == 2) then
    print*,'Prior inflation assuming Beta dist'
  endif
else
  print*,'Prior inflation is turned off!'
endif

print*,'------------------------------------'




! Initialize a repeating random sequence, seed fixed at 1
call init_random_seq(r_seq, 1)

call nc_check( nf90_open(in_file, nf90_nowrite, ncid), &
               'reading_in_data', 'opening file '//trim(in_file))

call getdimlen(ncid, "Members", nmembers)
call getdimlen(ncid, "Times", ntimes)


!allocate(control_sic(ntimes+1,nmembers))
!allocate(control_sit(ntimes+1,nmembers))
!allocate(control_snow(ntimes+1,nmembers))

!allocate(assim(ntimes+1,nmembers-1))
!allocate(input(ntimes+1,nmembers-1))
!allocate(preassim(ntimes+1,nmembers-1))
!allocate(output(ntimes+1,nmembers-1))
!allocate(obs(ntimes+1))
!allocate(obs_var(ntimes+1))


allocate(forcings_sic(ntimes,nmembers))
allocate(ICs_sic(nmembers))
allocate(forcings_sit(ntimes,nmembers))
allocate(ICs_sit(nmembers))
allocate(forcings_snow(ntimes,nmembers))
allocate(ICs_snow(nmembers))

call nc_check( nf90_inq_varid(ncid, 'aice_ics', varid), &
               'getvar_real', 'inquire var '// trim('aice_ics'))

call nc_check( nf90_get_var(ncid, varid, ICs_sic ), &
               'getvar_real', 'getting var '// trim('aice_ics'))

call nc_check( nf90_inq_varid(ncid, 'aice_forcings', varid), &
               'getvar_real', 'inquire var '// trim('aice_forcings'))

call nc_check( nf90_get_var(ncid, varid, forcings_sic ), &
               'getvar_real', 'getting var '// trim('aice_forcings'))
!!!!!!!!!!!!!!!!!!!!!!!!
call nc_check( nf90_inq_varid(ncid, 'hi_ics', varid), &
               'getvar_real', 'inquire var '// trim('hi_ics'))

call nc_check( nf90_get_var(ncid, varid, ICs_sit ), &
               'getvar_real', 'getting var '// trim('aice_ics'))!

call nc_check( nf90_inq_varid(ncid, 'hi_forcings', varid), &
               'getvar_real', 'inquire var '// trim('hi_forcings'))

call nc_check( nf90_get_var(ncid, varid, forcings_sit ), &
               'getvar_real', 'getting var '// trim('hi_forcings'))
!!!!!!!!!!!!!!!!!!!!!!!!!
call nc_check( nf90_inq_varid(ncid, 'hs_ics', varid), &
               'getvar_real', 'inquire var '// trim('hs_ics'))

call nc_check( nf90_get_var(ncid, varid, ICs_snow ), &
               'getvar_real', 'getting var '// trim('hs_ics'))

call nc_check( nf90_inq_varid(ncid, 'hs_forcings', varid), &
               'getvar_real', 'inquire var '// trim('hs_forcings'))

call nc_check( nf90_get_var(ncid, varid, forcings_snow ), &
               'getvar_real', 'getting var '// trim('hs_forcings'))

call nc_check( nf90_close(ncid), &
               'idealize_seaice', 'closing file '//trim(in_file))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Create the free forecasts for the different sea ice variables

allocate(prior_sic(nmembers))
allocate(post_sic(nmembers))
allocate(prior_sit(nmembers))
allocate(post_sit(nmembers))
allocate(prior_snow(nmembers))
allocate(post_snow(nmembers))

allocate(control_sic(ntimes+1,nmembers))
allocate(control_sit(ntimes+1,nmembers))
allocate(control_snow(ntimes+1,nmembers))

print*,'----------------------------------------------------'
print*,'Creating control ensemble from ICs and Forcings....'
do t=1, ntimes+1
  if (t==1) then
    print*,'Setting ICs'
    prior_sic(:) = ICs_sic
    prior_sit(:) = ICs_sit
    prior_snow(:) = ICs_snow
  else
    prior_sic = update_mem(post_sic,forcings_sic(t-1,:),1,r_seq,'SIC')
    prior_sit = update_mem(post_sit,forcings_sit(t-1,:),1,r_seq,'SIT')
    prior_snow = update_mem(post_snow,forcings_snow(t-1,:),1,r_seq,'SNO')
  endif
  !########################
  ! ASSIMILATION STEP!!!!
  !########################
  post_sic = prior_sic
  post_sit = prior_sit
  post_snow = prior_snow
  !########################
  ! FORECAST STEP !!!!
  !########################
  control_sic(t,:) = post_sic
  control_sit(t,:) = post_sit
  control_snow(t,:) = post_snow

end do
!########################################
rc = nf90_create('output_control.nc', NF90_CLOBBER, fid)
call nc_check(rc, 'create_output_file', 'creating "'//trim(out_file_control)//'"')

msgstring = 'Output from control ensemble run'
call set_global_char_att(fid, 'Data output:', msgstring)

call setup_sec_dim(fid, 'Times', nTimes+1, nTimesDimID)
call setup_sec_dim(fid, 'Members', 100, nMembersDimID)

call setup_sec_data_real(fid, 'aice', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'aice', 'description','output from test case')

call setup_sec_data_real(fid, 'hi', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'hi', 'description','output from test case')

call setup_sec_data_real(fid, 'hs', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'hs', 'description','output from test case')
rc = nf90_enddef(ncid)

do t=1,ntimes+1
  call write_sec_data_real(fid, t, 'aice', control_sic(t,:))
  call write_sec_data_real(fid, t, 'hi', control_sit(t,:))
  call write_sec_data_real(fid, t, 'hs', control_snow(t,:))
end do
rc = nf90_close(fid)
call nc_check(rc, 'close_output_file', 'closing '//trim(out_file_control))
!###########################################

deallocate(prior_sic, post_sic,prior_sit, post_sit,prior_snow, post_snow)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> END CREATING FREE FORECASTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Creating Observations for OSSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,'----------------------------------------------------'
print*,'Creating synthentic observations from truth members...'

truth_mem_ind = 100
truth_mem_sic = control_sic(:,truth_mem_ind)
truth_mem_sit = control_sit(:,truth_mem_ind)
truth_mem_snow = control_snow(:,truth_mem_ind)


if (obs_err_method == 4) then
  allocate(sicvals(n_sic))
  allocate(sicerr(n_sic))
  allocate(coefsicb(n_sic))
  allocate(coefsicc(n_sic))
  allocate(coefsicd(n_sic))
  call spline_sic_err(sicvals,sicerr,coefsicb,coefsicc,coefsicd)
endif


allocate(obsvar_sic(ntimes+1))
allocate(obs_sic(ntimes+1))
do i=1,ntimes+1
     obs_place = 99.0_r8
     call set_sic_err(truth_mem_sic(i),obs_err_method,obsvar_sic(i))
     do while(obs_place > 1.0_r8 .or. obs_place < 0.0_r8)
       if (obs_err_dist == 1) then
         obs_place = random_gaussian(r_seq, truth_mem_sic(i), sqrt(obsvar_sic(i)))
        else if (obs_err_dist == 2) then 
          obs_place = random_beta(r_seq, truth_mem_sic(i), sqrt(obsvar_sic(i)))
        endif
     end do
     if (filter_type/=3 .and. obs_err_method/=5) then
       if (i==1) then
         print*,'RESET OBS ERROR FROM OBS VALUE!!!!'
       endif
       call set_sic_err(obs_place,obs_err_method,obsvar_sic(i))
     endif
     obs_sic(i) = obs_place
end do



allocate(obsvar_sit(ntimes+1))
allocate(obs_sit(ntimes+1))
do i=1,ntimes+1
    !if (truth_mem_sic(i)<0.01_r8) then
    !  obs_sit(i) = -999.
    !  cycle
    ! endif
     obs_place = -99.0_r8
     obsvar_sit(i) = 0.1**2
     if (truth_mem_sic(i) == 0.0_r8) then
       hold_obs = 0.0_r8
     else
       hold_obs = truth_mem_sit(i)/truth_mem_sic(i)
     endif
     do while(obs_place < 0.0_r8)
       obs_place = random_gaussian(r_seq, hold_obs, sqrt(obsvar_sit(i)))
     end do
     obs_sit(i) = obs_place
end do

allocate(obsvar_snow(ntimes+1))
allocate(obs_snow(ntimes+1))
do i=1,ntimes+1
     obs_place = -99.0_r8
     obsvar_snow(i) = 0.1**2
     do while(obs_place < 0.0_r8)
       obs_place = random_gaussian(r_seq, truth_mem_snow(i), sqrt(obsvar_snow(i)))
     end do
     obs_snow(i) = obs_place
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> END MAKING OBSERVATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>RUNNING ASSIMILATION STEP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,'----------------------------------------------------'
print*,'Running assimilation cycling from ICs and Forcings....'
call assim_tools_init()


allocate(prior_sic(nmembers-1))
allocate(post_sic(nmembers-1))
allocate(prior_sit(nmembers-1))
allocate(post_sit(nmembers-1))
allocate(prior_snow(nmembers-1))
allocate(post_snow(nmembers-1))
allocate(obs_inc(nmembers-1))
allocate(state_inc(nmembers-1))


allocate(input_sic(ntimes+1,nmembers-1))
allocate(input_sit(ntimes+1,nmembers-1))
allocate(input_snow(ntimes+1,nmembers-1))
allocate(output_sic(ntimes+1,nmembers-1))
allocate(output_sit(ntimes+1,nmembers-1))
allocate(output_snow(ntimes+1,nmembers-1))

allocate(hold_FO_sic(ntimes+1,nmembers-1))
allocate(hold_FO_sit(ntimes+1,nmembers-1))
allocate(hold_FO_snow(ntimes+1,nmembers-1))

do t=1, ntimes+1
  !print*,t
  if (t==1) then
    print*,'Setting ICs'
    prior_sic(:) = ICs_sic(:ens_size)
    prior_sit(:) = ICs_sit(:ens_size)
    prior_snow(:) = ICs_snow(:ens_size)
  else
    prior_sic = update_mem(post_sic,forcings_sic(t-1,:ens_size),1,r_seq,'SIC')
    prior_sit = update_mem(post_sit,forcings_sit(t-1,:ens_size),1,r_seq,'SIT')
    prior_snow = update_mem(post_snow,forcings_snow(t-1,:ens_size),1,r_seq,'SNO')
    !!!!!!!
    prior_sit = compare_sic_sit(post_sic,prior_sic,prior_sit)
    !!!!!!!
  endif
  input_sic(t,:) = prior_sic
  input_sit(t,:) = prior_sit
  input_snow(t,:) = prior_snow
  !########################
  ! ASSIMILATION STEP!!!!
  !########################
  !!!!
  post_sic = prior_sic
  post_sit = prior_sit
  post_snow = prior_snow
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> ADD PRIOR INFLATION CODE HERE!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !> COMPUTE FOs if necessary
  do i=1,size(prior_sic)
    if (prior_sic(i) == 0.0_r8) then
      FO_SIC(i) = 0.0_r8
      FO_SIT(i) = 0.0_r8
      FO_SNOW(i) = 0.0_r8
    else
      FO_SIC(i) = prior_sic(i)
      FO_SIT(i) = prior_sit(i)/prior_sic(i)
      FO_SNOW(i) = prior_snow(i)/prior_sic(i)
    end if
  end do
!  FO_SIC = prior_sic
!  FO_SIT = prior_sit
!  FO_SNOW = prior_snow
  do i=1,size(obs_kind_table)
    if (obs_kind_table(i)=='AICE') then
      prior = FO_SIC
      prior_mean = sum(prior) / 99.0
      prior_var  = sum((prior - prior_mean)**2) / (99.0 - 1)
      obs_inc(:) = 0.0_r8
      call call_filter(prior,99,prior_mean,prior_var,obs_sic(t),obsvar_sic(t),obs_inc,filter_type,'SIC',obs_err_method)
    else if (obs_kind_table(i)=='HI') then
      prior = FO_SIT
      prior_mean = sum(prior) / 99.0
      prior_var  = sum((prior - prior_mean)**2) / (99.0 - 1)
      obs_inc(:) = 0.0_r8
      call call_filter(prior,99,prior_mean,prior_var,obs_sit(t),obsvar_sit(t),obs_inc,filter_type,'SIT',obs_err_method)
    else if (obs_kind_table(i)=='HS') then
      prior = FO_SNOW
      prior_mean = sum(prior) / 99.0
      prior_var  = sum((prior - prior_mean)**2) / (99.0 - 1)
      obs_inc(:) = 0.0_r8
      call call_filter(prior,99,prior_mean,prior_var,obs_snow(t),obsvar_snow(t),obs_inc,filter_type,'SNO',obs_err_method)
    endif
    if (obs_kind_table(i)=='AICE') then
      if (any(state_var_table(:)=='AICE')) then
        !print*,'Update SIC with SIC obs....'
        !call compute_regress(prior_sic,prior_sic,99,reg_coef)
        state_inc(:) = 0.0_r8
        call update_from_obs_inc(prior_sic,prior_mean,prior_var,obs_inc,prior_sic,99,state_inc,reg_coef,1.0_r8)
        post_sic = post_sic + state_inc
      endif
      if (any(state_var_table(:)=='HI')) then
        !print*,'Update SIT with SIC obs....'
        !call compute_regress(prior_sit,prior_sic,99,reg_coef)
        state_inc(:) = 0.0_r8
        call update_from_obs_inc(prior_sic,prior_mean,prior_var,obs_inc,prior_sit,99,state_inc,reg_coef,1.0_r8)
        post_sit = prior_sit + state_inc
      endif
      if (any(state_var_table(:)=='HS')) then
        call compute_regress(prior_snow,prior_sic,99,reg_coef)
        post_snow = prior_snow + (obs_inc*reg_coef)
        post_snow = max(0.0_r8,post_snow)
      endif
    else if (obs_kind_table(i)=='HI') then
      if (any(state_var_table(:)=='HI')) then
        !print*,'Update SIC with SIT obs...'
        !call compute_regress(prior_sic,prior_sit,99,reg_coef)
        state_inc(:) = 0.0_r8
        call update_from_obs_inc(prior_sit,prior_mean,prior_var,obs_inc,prior_sit,99,state_inc,reg_coef,1.0_r8)
        post_sit = prior_sit + state_inc
      endif
      if (any(state_var_table(:)=='AICE')) then
        !call compute_regress(prior_sit,prior_sit,99,reg_coef)
        state_inc(:) = 0.0_r8
        call update_from_obs_inc(prior_sit,prior_mean,prior_var,obs_inc,prior_sic,99,state_inc,reg_coef,1.0_r8)
        post_sic = prior_sic + state_inc
        !post_sit = max(0.0_r8,post_sit)
      endif
      if (any(state_var_table(:)=='HS')) then
        call compute_regress(prior_snow,prior_sit,99,reg_coef)
        post_snow = prior_snow + (obs_inc*reg_coef)
        !post_snow = max(0.0_r8,post_snow)
      endif
    else if (obs_kind_table(i)=='HS') then
      if (any(state_var_table(:)=='AICE')) then
        call compute_regress(prior_sic,prior_snow,99,reg_coef)
        post_sic = prior_sic + (obs_inc*reg_coef)
        !post_sic = max(0.0_r8,post_sic)
        !post_sic = min(1.0_r8,post_sic)
      endif
      if (any(state_var_table(:)=='HI')) then
        call compute_regress(prior_sit,prior_snow,99,reg_coef)
        post_sit = prior_sit + (obs_inc*reg_coef)
        !post_sit = max(0.0_r8,post_sit)
      endif
      if (any(state_var_table(:)=='HS')) then
        call compute_regress(prior_snow,prior_snow,99,reg_coef)
        post_snow = prior_snow + (obs_inc*reg_coef)
        !post_snow = max(0.0_r8,post_snow)
      endif
    endif
    prior_sic = post_sic
    prior_sit = post_sit
    prior_snow = post_snow
  end do
  !########################
  ! FORECAST STEP !!!!
  !########################
  hold_FO_sic(t,:) = FO_SIC
  hold_FO_sit(t,:) = FO_SIT
  hold_FO_snow(t,:) = FO_SNOW
  !!!!!!!
  !!!!!!!
  post_sit = compare_sic_sit(input_sic(t,:),post_sic,post_sit)
  !!!!!!!
  post_sic = max(0.0_r8,post_sic)
  post_sic = min(1.0_r8,post_sic)
  post_sit = max(0.0_r8,post_sit)
  post_snow = max(0.0_r8,post_snow)
  !!!!!!!
  output_sic(t,:) = post_sic
  output_sit(t,:) = post_sit
  output_snow(t,:) = post_snow
end do
!########################################
rc = nf90_create(out_file_assim, NF90_CLOBBER, fid)
call nc_check(rc, 'create_output_file', 'creating "'//trim(out_file_assim)//'"')

msgstring = 'Output from assim ensemble run'
call set_global_char_att(fid, 'Data output:', msgstring)

call setup_sec_dim(fid, 'Times', nTimes+1, nTimesDimID)
call setup_sec_dim(fid, 'Members', 100-1, nMembersDimID)

call setup_sec_data_real(fid, 'input_aice', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'input_aice', 'description','output from test case')

call setup_sec_data_real(fid, 'input_hi', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'input_hi', 'description','output from test case')

call setup_sec_data_real(fid, 'input_hs', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'input_hs', 'description','output from test case')

call setup_sec_data_real(fid, 'output_aice', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'output_aice', 'description','output from test case')

call setup_sec_data_real(fid, 'output_hi', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'output_hi', 'description','output from test case')

call setup_sec_data_real(fid, 'output_hs', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'output_hs', 'description','output from test case')

call setup_sec_data_real1d(fid, 'obs_aice', nTimesDimID)
call set_var_char_att(fid, 'obs_aice', 'description','obs from assimilation')

call setup_sec_data_real1d(fid, 'obs_hi', nTimesDimID)
call set_var_char_att(fid, 'obs_hi', 'description','obs from assimilation')

call setup_sec_data_real1d(fid, 'obs_hs', nTimesDimID)
call set_var_char_att(fid, 'obs_hs', 'description','obs from assimilation')

call setup_sec_data_real1d(fid, 'obsvar_aice', nTimesDimID)
call set_var_char_att(fid, 'obsvar_aice', 'description','obs error from test case')

call setup_sec_data_real1d(fid, 'obsvar_hi', nTimesDimID)
call set_var_char_att(fid, 'obsvar_hi', 'description','obs error from test case')

call setup_sec_data_real1d(fid, 'obsvar_hs', nTimesDimID)
call set_var_char_att(fid, 'obsvar_hs', 'description','obs error from test case')

call setup_sec_data_real(fid, 'FO_SIC', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'FO_SIC', 'description','output from test case')

call setup_sec_data_real(fid, 'FO_SIT', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'FO_SIT', 'description','output from test case')

call setup_sec_data_real(fid, 'FO_SNOW', nMembersDimID, nTimesDimID)
call set_var_char_att(fid, 'FO_SNOW', 'description','output from test case')

rc = nf90_enddef(ncid)

do t=1,ntimes+1
  call write_sec_data_real(fid, t, 'input_aice', input_sic(t,:))
  call write_sec_data_real(fid, t, 'input_hi', input_sit(t,:))
  call write_sec_data_real(fid, t, 'input_hs', input_snow(t,:))
  call write_sec_data_real(fid, t, 'output_aice', output_sic(t,:))
  call write_sec_data_real(fid, t, 'output_hi', output_sit(t,:))
  call write_sec_data_real(fid, t, 'output_hs', output_snow(t,:))
  call write_sec_data_real1d(fid, t, 'obs_aice', obs_sic(t))
  call write_sec_data_real1d(fid, t, 'obsvar_aice', obsvar_sic(t))
  call write_sec_data_real1d(fid, t, 'obs_hi', obs_sit(t))
  call write_sec_data_real1d(fid, t, 'obsvar_hi', obsvar_sit(t))
  call write_sec_data_real1d(fid, t, 'obs_hs', obs_snow(t))
  call write_sec_data_real1d(fid, t, 'obsvar_hs', obsvar_snow(t))
  call write_sec_data_real(fid, t, 'FO_SIC', hold_FO_sic(t,:))
  call write_sec_data_real(fid, t, 'FO_SIT', hold_FO_sit(t,:))
  call write_sec_data_real(fid, t, 'FO_SNOW', hold_FO_snow(t,:))
end do
rc = nf90_close(fid)
call nc_check(rc, 'close_output_file', 'closing '//trim(out_file_assim))
!###########################################



call finalize_mpi_utilities()

!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function update_mem(posterior,forcings,dt,r,type)
  real(r8), intent(in)  :: posterior(:),forcings(:)
  integer,  intent(in) :: dt
  character(len=3), intent(in) :: type
  real(r8), dimension(size(posterior))  :: update_mem
  integer :: i,m
  type(random_seq_type), intent(inout) :: r
  real(r8) :: noise,diff

  do i=1, size(posterior)
    if (trim(type) == 'SIC') then
      update_mem(i) = posterior(i) + forcings(i)*dt
      if (update_mem(i) > 1.0_r8) then
        update_mem(i) = 0.999999998_r8
      else if (update_mem(i) < 0.0_r8) then
        update_mem(i) = 0.0_r8
      endif
    else
      update_mem(i) = posterior(i) + forcings(i)*dt
      if (update_mem(i) < 0.0_r8) then
        update_mem(i) = 0.0_r8
      endif
    endif
  end do
end function update_mem
!!!!!!!!!!!!!!!!!!!!!!!!
function compare_sic_sit(old_sic,new_sic,new_sit)
  real(r8), intent(in) :: old_sic(99),new_sic(99),new_sit(99)
  real(r8), dimension(99)  :: compare_sic_sit
  integer :: m
  
  compare_sic_sit(:) = 0.0_r8
  do m=1, size(new_sic)
    if (old_sic(m) /= 0.0_r8 .and. new_sic(m) == 0.0_r8) then
      compare_sic_sit(m) = 0.0_r8
    else if (old_sic(m) == 0.0_r8 .and. new_sic(m)/=0.0_r8) then
      if (new_sit(m) == 0.0_r8) then
        compare_sic_sit(m) = new_sic(m)*0.32_r8
      else
        !print*,new_sit(m)
        compare_sic_sit(m) = new_sit(m)
      endif
    else
      if (new_sic(m) == 0.0_r8 .and. new_sit(m) /= 0.0) then
        compare_sic_sit(m) = 0.0_r8
      else if (new_sic(m) /= 0.0_r8 .and. new_sit(m) == 0.0) then
        compare_sic_sit(m) = new_sic(m)*0.32_r8
      else
        compare_sic_sit(m) = new_sit(m)
      endif
    endif
  end do
end function compare_sic_sit
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_sic_err(sic_val,obs_err_method,obsvar)
  real(r8), intent(in) :: sic_val
  integer, intent(in) :: obs_err_method
  real(r8), intent(out) :: obsvar

  if (obs_err_method == 1) then
    if (sic_val<0.01_r8) then
      obsvar = (0.01_r8*0.15)**2
    else
      obsvar = (sic_val*0.15_r8)**2
    endif
  else if (obs_err_method == 2) then
    if (sic_val>0.99_r8) then
      obsvar = (0.15 - (0.99_r8*0.15))**2
    else
      obsvar = (0.15 - (sic_val*0.15))**2
    endif
  else if (obs_err_method == 3) then
    ratio = 0.24997550249974496/0.0225
    if (sic_val < 0.01_r8) then
      obsvar = 0.01_r8*0.15_r8
    else if (sic_val > 0.99_r8) then
      obsvar = (0.15_r8 - (0.99_r8*0.15_r8))
    else
      obsvar = (sic_val*(1.0_r8-sic_val))/ratio
    endif
  else if (obs_err_method == 4) then
    call compute_spline(sic_val,sicvals,sicerr,coefsicb,coefsicc,coefsicd,n_sic,sic_err)
    obsvar = sic_err*sic_err
  else if (obs_err_method == 5) then
    obsvar_sic(i) = const_err_val**2
  endif
  return
end subroutine set_sic_err
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine call_filter(prior,ens_size,prior_mean,prior_var,obs,obsvar,obs_inc,filter_type,obs_kind,obs_err_method)
  real(r8), intent(in) :: prior_mean,prior_var,obs,obsvar
  real(r8), intent(in) :: prior(ens_size)
  character(len=3), intent(in) :: obs_kind
  real(r8), intent(inout) :: obs_inc(ens_size)
  integer, intent(in) :: ens_size,filter_type,obs_err_method
  real(r8) :: a,obsvar_hold
  integer :: i
  real(r8) :: likelihood(ens_size)
  real(r8) :: like_sum,bound(2)
  logical :: is_bounded(2)
  a = 0.0_r8

  if (filter_type == 1) then
    call obs_increment_eakf(prior,ens_size,prior_mean,prior_var,obs,obsvar, &
      obs_inc,a)
    return
  else if (filter_type == 2) then
    call obs_increment_rank_histogram(prior, ens_size, prior_var, obs, obsvar, &
      obs_inc)
    return
  else if (filter_type == 3) then
    ! NEW REQUIRMENTS FOR BOUNDED RHF
    ! Bounded normal RHF
    if (obs_kind == 'SIC') then
      is_bounded = (/.true.,.true./)
      bound = (/0.0_r8, 1.0_r8/)
    else
      is_bounded = (/.true.,.false./)
      bound = (/0.0_r8, 1.0_r8/)
    endif
    ! Test bounded normal likelihood; Could use an arbitrary likelihood
    do i = 1, ens_size
       if (obs_kind == 'SIC' .and. obs_err_method /= 5) then
         call set_sic_err(prior(i),obs_err_method,obsvar_hold)
       else
         obsvar_hold = obsvar
       endif
       likelihood(i) = get_truncated_normal_like(prior(i), obs, obsvar_hold, is_bounded, bound)
    end do
    
    ! Normalize the likelihood here
    like_sum = sum(likelihood)
    ! If likelihood underflow, assume flat likelihood, so no increments
    if(like_sum <= 0.0_r8) then
      obs_inc = 0.0_r8
      return
    else
      likelihood = likelihood / like_sum
    endif
    !call obs_increment_bounded_norm_rhf(prior, 99, prior_var, obs(t), obs_var(t), &
    !  obs_inc, is_bounded(bounds_case, :), bound)
    call obs_increment_bounded_norm_rhf(prior, likelihood, ens_size, prior_var, &
       obs_inc, is_bounded, bound)
    return
  endif
end subroutine call_filter
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_regress(state,obsspace,ens_size,reg_coef)
  real(r8), intent(in) :: state(ens_size),obsspace(ens_size)
  integer, intent(in) :: ens_size
  real(r8), intent(out) :: reg_coef
  real(r8) :: state_mean,obsspace_mean,obs_prior_var,obs_state_cov
  
  state_mean = sum(state)/ens_size
  obsspace_mean = sum(obsspace)/ens_size
  obs_prior_var = sum((obsspace - obsspace_mean)**2) / (ens_size - 1)

  obs_state_cov = sum((state-state_mean)*(obsspace-obsspace_mean))/(ens_size-1)
  if (obs_prior_var > 0.0_r8) then
    reg_coef = obs_state_cov/obs_prior_var
  else
    reg_coef = 0.0_r8
  endif
  return
end subroutine compute_regress
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_global_char_att(fid, attname, attvalue)

integer,          intent(in) :: fid
character(len=*), intent(in) :: attname
character(len=*), intent(in) :: attvalue

integer :: rc

rc = nf90_put_att(fid, NF90_GLOBAL, attname, attvalue)
call nc_check(rc, 'set_global_char_att', 'adding global attribute "'//trim(attname)//'"')

end subroutine set_global_char_att
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
subroutine setup_sec_unlimdim(ncid, c1, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(out) :: id1

integer :: rc

rc = nf90_def_dim(ncid, name=c1, len=NF90_UNLIMITED, dimid=id1)
call nc_check(rc, 'setup_sec_unlimdim', 'adding dimension "'//trim(c1)//'"')

end subroutine setup_sec_unlimdim
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
subroutine write_sec_data_real(ncid, col, c1, a1)

integer,          intent(in) :: ncid
integer,          intent(in) :: col
character(len=*), intent(in) :: c1
real(r8),         intent(in) :: a1(:)

integer :: rc, id1

rc = nf90_inq_varid(ncid, c1, id1)
call nc_check(rc, 'write_sec_data_real', 'querying variable "'//trim(c1)//'"')

rc = nf90_put_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'write_sec_data_real', 'writing variable "'//trim(c1)//'"')

end subroutine write_sec_data_real
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_var_char_att(fid, varname, attname, attvalue)

integer,          intent(in) :: fid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
character(len=*), intent(in) :: attvalue

integer :: rc, id1

rc = nf90_inq_varid(fid, varname, id1)
call nc_check(rc, 'set_var_char_att', 'inquiring variable id "'//trim(varname)//'"')

rc = nf90_put_att(fid, id1, attname, attvalue)
call nc_check(rc, 'set_var_char_att', 'adding attribute "'//trim(attname)//'" to variable "'//&
                                      trim(varname)//'"')

end subroutine set_var_char_att
!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setup_sec_dim(ncid, c1, n1, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: n1
integer,          intent(out) :: id1

integer :: rc

rc = nf90_def_dim(ncid, name=c1, len=n1, dimid=id1)
call nc_check(rc, 'setup_sec_dim', 'adding dimension "'//trim(c1)//'"')

end subroutine setup_sec_dim
!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setup_sec_data_real(ncid, c1, d1, d2)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: c1
integer,                    intent(in)  :: d1, d2

integer :: rc, id1

rc = nf90_def_var(ncid, name=c1, xtype=nf90_double, dimids=(/ d1, d2 /), varid=id1)
call nc_check(rc, 'setup_sec_data_real', 'defining variable "'//trim(c1)//'"')

end subroutine setup_sec_data_real
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_sec_data_real1d(ncid, col, c1, a1)

integer,          intent(in) :: ncid
integer,          intent(in) :: col
character(len=*), intent(in) :: c1
real(r8),         intent(in) :: a1

integer :: rc, id1

rc = nf90_inq_varid(ncid, c1, id1)
call nc_check(rc, 'write_sec_data_int1d', 'querying variable "'//trim(c1)//'"')

rc = nf90_put_var(ncid, id1, a1, start=(/ col /))
call nc_check(rc, 'write_sec_data_int1d', 'writing variable "'//trim(c1)//'"')

end subroutine write_sec_data_real1d
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setup_sec_data_real1d(ncid, c1, d1)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: c1
integer,                    intent(in)  :: d1

integer :: rc, id1

rc = nf90_def_var(ncid, name=c1, xtype=nf90_double, dimids=(/ d1 /), varid=id1)
call nc_check(rc, 'setup_sec_data_int1d', 'defining variable "'//trim(c1)//'"')

end subroutine setup_sec_data_real1d
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflate_prior_beta(inflate_val,ens_size,prior_mean,prior_var,prior)
  real(r8), intent(inout) :: inflate_val,prior_mean,prior_var
  integer, intent(in) :: ens_size
  real(r8), intent(inout) :: prior(ens_size)
  real(r8) :: theo_max_var,inflate_prior_var,a,b,inflate_a,inflate_b
  real(r8) :: cdf_val,inflate_prior_mean
  logical :: change_mean = .false.
  
  if (change_mean) then
    !inflate_var = prior_var*inflate_val
    !org_theo_max = org_mean*(1.0_r8-org_mean)
    !upper_bound = 0.5_r8*(sqrt(1.0_r8 - 4.0_r8*org_var) + 1.0_r8)
    !diff_mean = upper_bound - org_mean
    !diff_var = org_theo_max - org_var
    !hypot = np.sqrt(diff_mean**2 + diff_var**2)
    !inflate_upper_bound = 0.5_r8*(np.sqrt(1.0_r8 - 4.0_r8*inflate_var) + 1.0_r8)
    !new_mean = inflate_upper_bound - diff_mean
    !if new_mean < 0.5:
  !new_mean = 0.5

    !a = prior_mean*(((prior_mean*(1.0-prior_mean))/prior_var) - 1.0)
    !b = (1.0-prior_mean)*(((prior_mean*(1.0-prior_mean))/prior_var) - 1.0)
    print*,prior_mean,prior_var
    !inflate_a = inflate_prior_mean*(((inflate_prior_mean*(1.0-inflate_prior_mean))/inflate_prior_var) - 1.0)
    !inflate_b = (1.0-inflate_prior_mean)*(((inflate_prior_mean*(1.0-inflate_prior_mean))/inflate_prior_var) - 1.0)
    !print*,inflate_a,inflate_b
    !print*,prior_mean,inflate_prior_var
    !do i=1, ens_size
    !  !print*,prior(i)
    !  call betacdf(prior(i),a,b,cdf_val)
    !  call betaincinv(cdf_val,inflate_a,inflate_b,prior(i))
    !end do  
    !endif

  else 
    if (prior_var*inflate_val > prior_mean*(1.0_r8-prior_mean)) then
      theo_max_var = prior_mean*(1.0_r8 - prior_mean)
      if (theo_max_var > prior_var) then 
        !print*,'Theo var is used'
        !inflate_prior_var = theo_max_var - 1.0d-5
        return
      else
        !print*,'Can not inflate prior ensemble....continue'
        return
      endif  
    else
      !print*,'Using full inflate_val'
      inflate_prior_var = prior_var*inflate_val
    endif

    a = prior_mean*(((prior_mean*(1.0-prior_mean))/prior_var) - 1.0)
    b = (1.0-prior_mean)*(((prior_mean*(1.0-prior_mean))/prior_var) - 1.0)
    !print*,prior_mean,prior_var
    inflate_a = prior_mean*(((prior_mean*(1.0-prior_mean))/inflate_prior_var) - 1.0)
    inflate_b = (1.0-prior_mean)*(((prior_mean*(1.0-prior_mean))/inflate_prior_var) - 1.0)
    !print*,prior_mean,inflate_prior_var
    do i=1, ens_size
      !print*,prior(i)
      call betacdf(prior(i),a,b,cdf_val)
      call betaincinv(cdf_val,inflate_a,inflate_b,prior(i))
    end do  
  endif
end subroutine inflate_prior_beta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SPLINE FUNCTION
!> subroutine to initialize sic error splines
subroutine spline_sic_err(x,y,b,c,d)
integer, parameter :: n=35
real(r8), intent(out), dimension (n) :: b(n),c(n),d(n),x(n),y(n)
real(r8) :: h
integer :: i,j,gap,rcio

open(unit=1, file='sic_err_list.txt',form='formatted',action='read',status='old')
do i=1,n
  read(1,*,iostat=rcio)x(i),y(i)
  if (rcio .ne. 0) exit
!  print*,x(i),y(i)
end do

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if

!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do

!
! step 2: end conditions
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if

!
! step 3: forward elimination
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do

!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline_sic_err
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> subroutine to compute spline sit error val
subroutine compute_spline(u, x, y, b, c, d, n,ispline)
integer, intent(in) :: n
real(r8), intent(in), dimension (n) :: b(n),c(n),d(n),x(n),y(n)
real(r8), intent(in) :: u
real(r8), intent(out) :: ispline

integer :: i,j,k
real(r8) :: dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do

!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end subroutine compute_spline
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program test_bounded_rhf
