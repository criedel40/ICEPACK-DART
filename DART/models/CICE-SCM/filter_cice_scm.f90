program filter_cice_scm

! Test specific cases for obs_increment_bounded_norm_rhf
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

use     obs_utilities_mod, only : add_obs_to_seq, create_3d_obs, &
                                  getdimlen, getvar_int, set_missing_name,&
                                  getvar_real
use time_manager_mod,     only : generate_seed, time_type, set_time

use adaptive_inflate_mod, only : update_inflation, inflate_ens

use types_mod,         only : r8,metadatalength, MAX_FILES, MAX_NUM_DOMS
use assim_tools_mod,   only : obs_increment_bounded_norm_rhf, &
                              obs_increment_rank_histogram, &
                              obs_increment_eakf, &
                              get_truncated_normal_like, &
                              update_from_obs_inc, &
                              assim_tools_init
                              
use random_seq_mod,    only : random_seq_type, init_random_seq, random_gaussian,random_beta, &
                              random_trunc_gaussian,random_truncnorm_a,random_truncnorm_ab

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, nc_check, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode
use         utilities_mod, only : register_module, error_handler,               &
                                  E_ERR, E_MSG, get_unit,  &
                                  do_output, to_upper,    &
                                  find_namelist_in_file, check_namelist_read,   &
                                  file_exist, find_textfile_dims, file_to_text, &
                                  do_nml_file, do_nml_term, set_multiple_filename_lists

use default_model_mod,     only : init_time, init_conditions, adv_1step, &
                                  nc_write_model_vars

use beta_dist_functions,   only : betacdf,betaincinv

use netcdf
implicit none

integer :: ncid,iunit,io,ncid2
integer,  parameter :: ens_size = 79
integer :: filter_type = 1
integer :: obs_err_method = 1
integer :: obs_err_dist = 1
real(r8) :: const_err_val = 0.01_r8
logical :: inflate_prior = .false.
integer :: inflation_type = 1
real(r8) :: inflate_val = 2.0_r8
integer :: cyclenum = 1

logical :: check_state_var(3) = (/.true.,.true.,.true./)

real(r8), allocatable :: prior(:)
real(r8), allocatable :: input_state_vector(:,:,:),output_state_vector(:,:,:),work(:,:),work1D(:,:),preassim_state_vector(:,:,:)
real(r8), allocatable :: orginal_state_vector(:,:,:,:)
real(r8), allocatable :: aicen_org(:,:,:),vicen_org(:,:,:),vsnon_org(:,:,:)
real(r8), allocatable :: obs_inc(:),state_inc(:)
real(r8), allocatable :: truth_values(:), obs_values(:), obs_errs(:)
real(r8), allocatable :: state_vector_prior_FO(:,:),state_vector_post_FO(:,:)
real(r8), allocatable :: input_qice(:,:,:,:),input_sice(:,:,:,:),input_qsnow(:,:,:,:)
real(r8), allocatable :: input_tsfcn(:,:,:)!,output_tsfcn(:,:,:)
character(len=31), allocatable :: write_out_vars(:)
character(len=NF90_MAX_NAME), allocatable :: state_var_table(:),obs_kind_table(:)
character(len=256), allocatable :: file_array_input(:), file_array_output(:)

integer,parameter       :: n_sic = 35
real(r8), allocatable :: sicvals(:),sicerr(:),coefsicb(:),coefsicc(:),coefsicd(:)
character (len=3) :: nchar

type(time_type)       :: ens_time
real(r8)              :: sic_err,hold_mean,hold_stddev
real(r8)              :: prior_mean, prior_var,obs_place
real(r8)              :: ratio,reg_coef,time,time_forc,hold_err
type(random_seq_type) :: r_seq
integer               :: i,varid,varid2,f,dimid,nc,ob,dimid_obsnum,dimid_strlen
integer               :: aicen_index,vicen_index,vsnon_index
integer               :: save_index,iflag,seed,timesec
integer               :: icepack_data_index
integer               :: snow_enthapy_levs = 8
integer               :: ice_enthapy_sal_levs = 3
integer               :: NCAT,status,NX,istep0,istep1,dimid_ni,dimid_ncat,dimid_enssize
integer :: start2(2),count2(2)
character(len=metadatalength) :: model_state_vars(3 * 1 ) = ' '
character(len=metadatalength) :: obs_kinds(3 * 1 ) = ' '
character(len=NF90_MAX_NAME) :: varname

character(len=256) ::  input_state_files(MAX_FILES) = ''
character(len=256) :: output_state_files(MAX_FILES) = ''

character(len=256) ::  input_state_file_list(MAX_NUM_DOMS) = ''
character(len=256) :: output_state_file_list(MAX_NUM_DOMS) = ''
character(len=256) :: truth_mem_data = 'truth_mem_data.nc'
character(len=256) :: filename

namelist /seaice_nml/  &
   filter_type, &
   obs_err_method, &
   obs_err_dist, &
   const_err_val, &
   inflate_prior, &
   inflation_type, &
   inflate_val, &
   model_state_vars, &
   obs_kinds, &
   input_state_file_list, &
   output_state_file_list, &
   truth_mem_data, &
   icepack_data_index, &
   snow_enthapy_levs, &
   ice_enthapy_sal_levs, &
   cyclenum
   

! Initialize for error handler
call initialize_mpi_utilities('filter_cice_scm')

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
    !call to_upper(varname)
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
  allocate(write_out_vars(save_index))
  do i = 1, 3
    if (obs_kinds(i) == ' ') exit
    write_out_vars(i) = obs_kinds(i)
    varname = trim(obs_kinds(i))
    !call to_upper(varname)
    obs_kind_table(i) = varname
  enddo
endif
print*,'Model State Variables'
print*,'----------------------'
do i=1,size(state_var_table)
  write(6,*),trim(state_var_table(i))
  if (trim(state_var_table(i))=='aicen') then
    check_state_var(1) = .false.
  else if (trim(state_var_table(i))=='vicen') then
    check_state_var(2) = .false.
  else if (trim(state_var_table(i))=='vsnon') then
    check_state_var(3) = .false.
  endif
enddo
print*,'----------------------'
print*,'Assimilated Observation Types'
print*,'----------------------'
do i=1,size(obs_kind_table)
  write(6,*),trim(obs_kind_table(i))
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


call set_multiple_filename_lists(input_state_files(:), &
                                 input_state_file_list(:), &
                                 1, &
                                 ens_size,     &
                                 'filter_cice_scm','input_state_files','input_state_file_list')

call set_multiple_filename_lists(output_state_files(:), &
                                 output_state_file_list(:), &
                                 1, &
                                 ens_size,     &
                                 'filter_cice_scm','output_state_files','output_state_file_list')

allocate(file_array_input(ens_size), file_array_output(ens_size))

file_array_input  = RESHAPE(input_state_files,  (/ens_size/))
file_array_output = RESHAPE(output_state_files, (/ens_size/))

print*,'------------------------------------'
do i=1,size(file_array_input)
  filename = trim(file_array_input(i))
  if (filename == trim(truth_mem_data)) then
    write(6,*) "Truth member data read into state vector...STOP!"
    STOP
  endif
  write(6,*) 'Using restart dump = ', trim(filename)
  status = nf90_open(trim(filename), nf90_nowrite, ncid)
  if (status /= nf90_noerr) then
    write(6,*) "Netcdf file was not open correctly..STOP!"
    STOP
  endif
  if (i == 1) then
    status = nf90_inq_dimid(ncid,"ncat",dimid)
    status = nf90_inquire_dimension(ncid,dimid,len=NCAT)
    status = nf90_inq_dimid(ncid,"ni",dimid)
    status = nf90_inquire_dimension(ncid,dimid,len=NX)
    status = nf90_get_att(ncid, nf90_global, 'istep1', istep0)
    istep1 = istep0
    status = nf90_get_att(ncid, nf90_global, 'time', time)
    status = nf90_get_att(ncid, nf90_global, 'time_forc', time_forc)
    allocate(input_state_vector(size(state_var_table),NCAT,ens_size))
    allocate(output_state_vector(size(state_var_table),NCAT,ens_size))
    allocate(preassim_state_vector(size(state_var_table),NCAT,ens_size))
    allocate(orginal_state_vector(size(state_var_table),NX,NCAT,ens_size))
    allocate(input_qice(NX,NCAT,ice_enthapy_sal_levs,ens_size))
    allocate(input_sice(NX,NCAT,ice_enthapy_sal_levs,ens_size))
    allocate(input_qsnow(NX,NCAT,snow_enthapy_levs,ens_size))
    allocate(input_tsfcn(NX,NCAT,ens_size))
    allocate(work(NX,NCAT))
    allocate(aicen_org(NX,NCAT,ens_size))
    allocate(vicen_org(NX,NCAT,ens_size))
    allocate(vsnon_org(NX,NCAT,ens_size))
    !work(:,:) = 0.0_r8
    start2(1) = 1
    count2(1) = NX
    start2(2) = 1
    count2(2) = NCAT
    allocate(state_vector_prior_FO(size(obs_kind_table),ens_size))
  allocate(state_vector_post_FO(size(obs_kind_table),ens_size))
  endif
  do f=1,size(state_var_table)
    status = nf90_inq_varid(ncid, state_var_table(f), VarId)
    status = nf90_get_var(ncid, VarId, work, &
               start=start2, &
               count=count2 )
    orginal_state_vector(f,:,:,i) = work
    input_state_vector(f,:,i) = work(icepack_data_index,:)
  enddo
  do f=1, ice_enthapy_sal_levs
    write(nchar,'(i3.3)') f
    status = nf90_inq_varid(ncid, 'qice'//trim(nchar), VarId)
    status = nf90_get_var(ncid, VarId, input_qice(:,:,f,i), &
               start=start2, &
               count=count2 )
    status = nf90_inq_varid(ncid, 'sice'//trim(nchar), VarId)
    status = nf90_get_var(ncid, VarId, input_sice(:,:,f,i), &
               start=start2, &
               count=count2 )
  enddo
  do f=1, snow_enthapy_levs
    write(nchar,'(i3.3)') f
    status = nf90_inq_varid(ncid, 'qsno'//trim(nchar), VarId)
    status = nf90_get_var(ncid, VarId, input_qsnow(:,:,f,i), &
               start=start2, &
               count=count2 )
  enddo
  status = nf90_inq_varid(ncid, 'Tsfcn', VarId)
  status = nf90_get_var(ncid, VarId, input_tsfcn(:,:,i), &
               start=start2, &
               count=count2 )
  !!!!!!!! KEEP COPY OF ORGINALS IF VARIABLE
  !!!!!!!! IS NOT A STATE VECTOR
  status = nf90_inq_varid(ncid, 'aicen', VarId)
  status = nf90_get_var(ncid, VarId, aicen_org(:,:,i), &
               start=start2, &
               count=count2 )
  status = nf90_inq_varid(ncid, 'vicen', VarId)
  status = nf90_get_var(ncid, VarId, vicen_org(:,:,i), &
               start=start2, &
               count=count2 )
  status = nf90_inq_varid(ncid, 'vsnon', VarId)
  status = nf90_get_var(ncid, VarId, vsnon_org(:,:,i), &
               start=start2, &
               count=count2 )
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do f=1,size(obs_kind_table)
    if (obs_kind_table(f)=='si_concentration') then
      call FO_SIC(ncid,icepack_data_index,state_vector_prior_FO(f,i),hold_err)
    else if (obs_kind_table(f)=='si_thickness') then
      call FO_SIT(ncid,icepack_data_index,state_vector_prior_FO(f,i),hold_err)
    else
      write(6,*) 'Obs type is not supported...STOP!'
      STOP
    endif
  enddo
  status = nf90_close(ncid)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Creating Observations for OSSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,'------------------------------------'
! Initialize a repeating random sequence, seed fixed at 1
print*,'Creating synthentic observations from truth members...'
filename = trim(truth_mem_data)
write(6,*) 'Using restart dump for truth = ', trim(filename)
status = nf90_open(trim(filename), nf90_nowrite, ncid)
if (status /= nf90_noerr) then
  write(6,*) "Netcdf file was not open correctly..STOP!"
  STOP
endif
status = nf90_get_att(ncid, nf90_global, 'time', timesec)
ens_time = set_time(0,timesec)
seed = generate_seed(ens_time)
call init_random_seq(r_seq, seed)

allocate(truth_values(size(obs_kind_table)))
allocate(obs_values(size(obs_kind_table)))
allocate(obs_errs(size(obs_kind_table)))

do f=1,size(obs_kind_table)
  if (trim(obs_kind_table(f))=='si_concentration') then
    call FO_SIC(ncid,icepack_data_index,truth_values(f),obs_errs(f))
    obs_place = -99.0_r8
    do while(obs_place > 1.0_r8 .or. obs_place < 0.0_r8)
      obs_place = random_gaussian(r_seq, truth_values(f), sqrt(obs_errs(f)))
    enddo
    obs_values(f) = obs_place
    !obs_values(f) = random_truncnorm_ab(r_seq,truth_values(f),sqrt(obs_errs(f)),0.0_r8,1.0_r8)
  else if (trim(obs_kind_table(f))=='si_thickness') then
    call FO_SIT(ncid,icepack_data_index,truth_values(f),obs_errs(f))
    obs_place = -99.0_r8
    do while(obs_place < 0.0_r8)
      obs_place = random_gaussian(r_seq, truth_values(f), sqrt(obs_errs(f)))
    end do
    obs_values(f) = obs_place
    !obs_values(f) = random_truncnorm_a(r_seq,truth_values(f),sqrt(obs_errs(f)),0.0_r8)
  else
    write(6,*) 'Obs type is not supported...STOP!'
    STOP
  endif
enddo
status = nf90_close(ncid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> END MAKING OBSERVATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>RUNNING ASSIMILATION STEP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,'----------------------------------------------------'
print*,'Running assimilation cycling from ICs and Forcings....'
print*,'----------------------------------------------------'
call assim_tools_init()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ADD PRIOR INFLATION CODE HERE!
if (inflate_prior) then
  ! FOR NOW NO PRIOR INFLATE TURNED ON
  preassim_state_vector(:,:,:) = input_state_vector(:,:,:)
  output_state_vector(:,:,:) = preassim_state_vector(:,:,:)
else
  preassim_state_vector(:,:,:) = input_state_vector(:,:,:)
  output_state_vector(:,:,:) = preassim_state_vector(:,:,:)
endif


allocate(obs_inc(ens_size))
allocate(state_inc(ens_size))

do ob=1, size(obs_kind_table)
  print*,'----------------------------------------------------'
  write(6,*) 'Assimilating Ob => ',trim(obs_kind_table(ob))
  !########################
  ! ASSIMILATION STEP!!!!
  !#######################
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (obs_kind_table(ob)=='si_concentration') then
    prior = state_vector_prior_FO(ob,:)
    prior_mean = sum(prior) / ens_size
    prior_var  = sum((prior - prior_mean)**2) / (ens_size - 1)
    obs_inc(:) = 0.0_r8
    call call_filter(prior,ens_size,prior_mean,prior_var,obs_values(ob),obs_errs(ob),obs_inc,filter_type,'SIC',obs_err_method)
    state_vector_post_FO(ob,:) = prior + obs_inc
  else if (obs_kind_table(ob)=='si_thickness') then
    prior = state_vector_prior_FO(ob,:)
    prior_mean = sum(prior) / ens_size
    prior_var  = sum((prior - prior_mean)**2) / (ens_size - 1)
    obs_inc(:) = 0.0_r8
    call call_filter(prior,ens_size,prior_mean,prior_var,obs_values(ob),obs_errs(ob),obs_inc,filter_type,'SIT',obs_err_method)
    state_vector_post_FO(ob,:) = prior + obs_inc
  else if (obs_kind_table(ob)=='snow_depth') then
    prior = state_vector_prior_FO(ob,:)
    prior_mean = sum(prior) / ens_size
    prior_var  = sum((prior - prior_mean)**2) / (ens_size - 1)
    obs_inc(:) = 0.0_r8
    call call_filter(prior,ens_size,prior_mean,prior_var,obs_values(ob),obs_errs(ob),obs_inc,filter_type,'SNO',obs_err_method)
    state_vector_post_FO(ob,:) = prior + obs_inc
  endif
  !print*,"HELLO HELLO HELLO"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Update all state vectors with incs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do f=1,size(state_var_table)
    write(6,*) 'Updating state field => ',trim(state_var_table(f))
    do nc=1, NCAT
      state_inc(:) = 0.0_r8
      call update_from_obs_inc(prior,prior_mean,prior_var,obs_inc,output_state_vector(f,nc,:),ens_size,state_inc,reg_coef,1.0_r8)
      output_state_vector(f,nc,:) = output_state_vector(f,nc,:) + state_inc
    enddo
  end do
end do

! call clenaup subroutine
aicen_index = -1
if (any(state_var_table=='aicen')) then
 call findindex1(state_var_table,size(state_var_table),trim('aicen'),aicen_index)
endif
vicen_index = -1
if (any(state_var_table=='vicen')) then
  call findindex1(state_var_table,size(state_var_table),trim('vicen'),vicen_index)
endif
  vsnon_index = -1
if (any(state_var_table=='vsnon')) then
  call  findindex1(state_var_table,size(state_var_table),trim('vsnon'),vsnon_index)
endif
if (vicen_index /= -1 .and. vsnon_index /= -1) then !.and. present(vsnon)) then
    do i=1, ens_size   
      call cice_cleanup(output_state_vector(aicen_index,:,i),output_state_vector(vicen_index,:,i),output_state_vector(vsnon_index,:,i), &
               input_qice(icepack_data_index,:,:,i),input_sice(icepack_data_index,:,:,i),input_qsnow(icepack_data_index,:,:,i),&
               input_tsfcn(icepack_data_index,:,i),NCAT,ice_enthapy_sal_levs,snow_enthapy_levs,aicen_org=aicen_org(icepack_data_index,:,i),vsnon_org=vsnon_org(icepack_data_index,:,i))
    enddo
else if (vicen_index /= -1 .and. vsnon_index == -1) then
  do i=1, ens_size
    call cice_cleanup(aicen=output_state_vector(aicen_index,:,i),vicen=output_state_vector(vicen_index,:,i),vsnon=vsnon_org(icepack_data_index,:,i), &
               input_qice=input_qice(icepack_data_index,:,:,i),input_sice=input_sice(icepack_data_index,:,:,i),input_qsno=input_qsnow(icepack_data_index,:,:,i),&
               input_tsfcn=input_tsfcn(icepack_data_index,:,i),NCAT=NCAT,QICE=ice_enthapy_sal_levs,QSNO=snow_enthapy_levs,aicen_org=aicen_org(icepack_data_index,:,i),vsnon_org=vsnon_org(icepack_data_index,:,i))
  enddo
else if (vicen_index == -1 .and. vsnon_index /= -1) then
  do i=1, ens_size
    call cice_cleanup(aicen=output_state_vector(aicen_index,:,i),vicen=vicen_org(icepack_data_index,:,i),vsnon=output_state_vector(vsnon_index,:,i), &
               input_qice=input_qice(icepack_data_index,:,:,i),input_sice=input_sice(icepack_data_index,:,:,i),input_qsno=input_qsnow(icepack_data_index,:,:,i),&
               input_tsfcn=input_tsfcn(icepack_data_index,:,i),NCAT=NCAT,QICE=ice_enthapy_sal_levs,QSNO=snow_enthapy_levs,aicen_org=aicen_org(icepack_data_index,:,i),vsnon_org=vsnon_org(icepack_data_index,:,i))
  enddo
endif
print*,'------------------------------------'
print*,'WRITE OUT ADJUST STATE TO NETCDF FILES'
do i=1,size(file_array_output)
  filename = trim(file_array_output(i))
  write(6,*) 'Restart dump = ', trim(filename)
  if (filename == trim(truth_mem_data)) then
    write(6,*) "Truth member data read into state vector...STOP!"
    STOP
  endif
  iflag = nf90_clobber
  status = nf90_create(trim(filename), iflag, ncid)
  status = nf90_put_att(ncid,nf90_global,'istep1',istep1)
  status = nf90_put_att(ncid,nf90_global,'time',time)
  status = nf90_put_att(ncid,nf90_global,'time_forc',time_forc)
  status = nf90_def_dim(ncid,'ni',nx,dimid_ni)
  status = nf90_def_dim(ncid,'ncat',ncat,dimid_ncat)
  do f=1,size(state_var_table)
    status = nf90_def_var(ncid,trim(state_var_table(f)),nf90_double,(/dimid_ni,dimid_ncat/),varid)
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (check_state_var(1)) then
    status = nf90_def_var(ncid,'aicen',nf90_double,(/dimid_ni,dimid_ncat/),varid)
  endif
  if (check_state_var(2)) then
    status = nf90_def_var(ncid,'vicen',nf90_double,(/dimid_ni,dimid_ncat/),varid)
  endif
  if (check_state_var(3)) then
    status = nf90_def_var(ncid,'vsnon',nf90_double,(/dimid_ni,dimid_ncat/),varid)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do f=1, ice_enthapy_sal_levs
    write(nchar,'(i3.3)') f
    status = nf90_def_var(ncid,'qice'//trim(nchar),nf90_double,(/dimid_ni,dimid_ncat/),varid)
    status = nf90_def_var(ncid,'sice'//trim(nchar),nf90_double,(/dimid_ni,dimid_ncat/),varid)
  enddo
  do f=1, snow_enthapy_levs
    write(nchar,'(i3.3)') f
    status = nf90_def_var(ncid,'qsno'//trim(nchar),nf90_double,(/dimid_ni,dimid_ncat/),varid)
  enddo
  status = nf90_def_var(ncid,'Tsfcn',nf90_double,(/dimid_ni,dimid_ncat/),varid)
  status = nf90_enddef(ncid)
  do f=1,size(state_var_table)
    orginal_state_vector(f,icepack_data_index,:,i) = output_state_vector(f,:,i)
    status = nf90_inq_varid(ncid,trim(state_var_table(f)),varid)
    status = nf90_put_var(ncid, varid, orginal_state_vector(f,:,:,i), &
               start=start2, &
               count=count2)
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (check_state_var(1)) then
    status = nf90_inq_varid(ncid,'aicen',varid)
    status = nf90_put_var(ncid, varid, aicen_org(:,:,i), &
               start=start2, &
               count=count2)
  endif
  if (check_state_var(2)) then
    status = nf90_inq_varid(ncid,'vicen',varid)
    status = nf90_put_var(ncid, varid, vicen_org(:,:,i), &
               start=start2, &
               count=count2)
  endif
  if (check_state_var(3)) then
    status = nf90_inq_varid(ncid,'vsnon',varid)
    status = nf90_put_var(ncid, varid, vsnon_org(:,:,i), &
               start=start2, &
               count=count2)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do f=1, ice_enthapy_sal_levs
    write(nchar,'(i3.3)') f 
    status = nf90_inq_varid(ncid,'qice'//trim(nchar),varid)
    status = nf90_put_var(ncid, varid, input_qice(:,:,f,i), &
               start=start2, &
               count=count2)
    status = nf90_inq_varid(ncid,'sice'//trim(nchar),varid)
    status = nf90_put_var(ncid, varid, input_sice(:,:,f,i), &
               start=start2, &
               count=count2)
  enddo
  do f=1, snow_enthapy_levs
    write(nchar,'(i3.3)') f
    status = nf90_inq_varid(ncid,'qsno'//trim(nchar),varid)
    status = nf90_put_var(ncid, varid, input_qsnow(:,:,f,i), &
               start=start2, &
               count=count2)
  enddo
  status = nf90_inq_varid(ncid,'Tsfcn',varid)
  status = nf90_put_var(ncid, varid, input_tsfcn(:,:,i), &
               start=start2, &
               count=count2)
  status = nf90_close(ncid)
enddo

allocate(work1D(NCAT,ens_size+2))
iflag = nf90_clobber
status = nf90_create(trim('input_mean.nc'), iflag, ncid)
status = nf90_put_att(ncid,nf90_global,'istep1',istep1)
status = nf90_put_att(ncid,nf90_global,'time',time)
status = nf90_put_att(ncid,nf90_global,'time_forc',time_forc)
status = nf90_def_dim(ncid,'ncat',ncat,dimid_ncat)
status = nf90_def_dim(ncid,'ens_size',ens_size+2,dimid_enssize)
do f=1,size(state_var_table)
  status = nf90_def_var(ncid,trim(state_var_table(f)),nf90_double,(/dimid_ncat,dimid_enssize/),varid)
enddo
status = nf90_enddef(ncid)
!!!!!!!!!!
status = nf90_create(trim('output_mean.nc'), iflag, ncid2)
status = nf90_put_att(ncid2,nf90_global,'istep1',istep1)
status = nf90_put_att(ncid2,nf90_global,'time',time)
status = nf90_put_att(ncid2,nf90_global,'time_forc',time_forc)
status = nf90_def_dim(ncid2,'ncat',ncat,dimid_ncat)
status = nf90_def_dim(ncid2,'ens_size',ens_size+2,dimid_enssize)
do f=1,size(state_var_table)
  status = nf90_def_var(ncid2,trim(state_var_table(f)),nf90_double,(/dimid_ncat,dimid_enssize/),varid)
enddo
status = nf90_enddef(ncid2)
!!!!!!!!!!
do i=1,size(state_var_table)
  status = nf90_inq_varid(ncid,trim(state_var_table(i)),varid)
  status = nf90_inq_varid(ncid2,trim(state_var_table(i)),varid2)
  work1D(:,:) = 0.0_r8
  do f=1, NCAT
    work1D(f,:ens_size) = input_state_vector(i,f,:)
    hold_mean = sum(input_state_vector(i,f,:)) / ens_size
    hold_stddev = sum((input_state_vector(i,f,:) - hold_mean)**2) / (ens_size - 1)
    work1D(f,ens_size+1) = hold_mean
    work1D(f,ens_size+2) = hold_stddev
  enddo
  status = nf90_put_var(ncid, varid, work1D, &
               start=(/1,1/), &
               count=(/NCAT,ens_size+2/))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  work1D(:,:) = 0.0_r8
  do f=1, NCAT
    work1D(f,:ens_size) = output_state_vector(i,f,:)
    hold_mean = sum(output_state_vector(i,f,:)) / ens_size
    hold_stddev = sum((output_state_vector(i,f,:) - hold_mean)**2) / (ens_size - 1)
    work1D(f,ens_size+1) = hold_mean
    work1D(f,ens_size+2) = hold_stddev
  enddo
  status = nf90_put_var(ncid2, varid2, work1D, &
               start=(/1,1/), &
               count=(/NCAT,ens_size+2/))
 
enddo
status = nf90_close(ncid)
status = nf90_close(ncid2)

status = nf90_create(trim('obs_seq.final.nc'), iflag, ncid)
status = nf90_put_att(ncid,nf90_global,'istep1',istep1)
status = nf90_put_att(ncid,nf90_global,'time',time)
status = nf90_put_att(ncid,nf90_global,'time_forc',time_forc)
status = nf90_def_dim(ncid,'num_obs',size(obs_kind_table),dimid_obsnum)
status = nf90_def_dim(ncid,'ens_size',ens_size,dimid_enssize)
status = nf90_def_dim(ncid,'stringlen',31,dimid_strlen)
!!!!!!
status = nf90_def_var(ncid,'obs_types',nf90_char,(/dimid_strlen,dimid_obsnum/),varid)
status = nf90_def_var(ncid,'prior_ens',nf90_double,(/dimid_enssize,dimid_obsnum/),varid)
status = nf90_def_var(ncid,'post_ens',nf90_double,(/dimid_enssize,dimid_obsnum/),varid)
status = nf90_def_var(ncid,'truth',nf90_double,(/dimid_obsnum/),varid)
status = nf90_def_var(ncid,'observation',nf90_double,(/dimid_obsnum/),varid)
status = nf90_def_var(ncid,'obs_var',nf90_double,(/dimid_obsnum/),varid)
!!!!!!
status = nf90_enddef(ncid)
!!!!!
status = nf90_inq_varid(ncid,'obs_types',varid)
status = nf90_put_var(ncid, varid, write_out_vars, &
               start=(/1,1/), &
               count=(/31,size(write_out_vars)/))
status = nf90_inq_varid(ncid,'prior_ens',varid)
status = nf90_put_var(ncid, varid, transpose(state_vector_prior_FO), &
               start=(/1,1/), &
               count=(/ens_size,size(write_out_vars)/))
status = nf90_inq_varid(ncid,'post_ens',varid)
status = nf90_put_var(ncid, varid, transpose(state_vector_post_FO), &
               start=(/1,1/), &
               count=(/ens_size,size(write_out_vars)/))
status = nf90_inq_varid(ncid,'truth',varid)
status = nf90_put_var(ncid, varid, truth_values, &
               start=(/1/), &
               count=(/size(write_out_vars)/))
status = nf90_inq_varid(ncid,'observation',varid)
status = nf90_put_var(ncid, varid, obs_values, &
               start=(/1/), &
               count=(/size(write_out_vars)/))
status = nf90_inq_varid(ncid,'obs_var',varid)
status = nf90_put_var(ncid, varid, obs_errs, &
               start=(/1/), &
               count=(/size(write_out_vars)/))
status = nf90_close(ncid)

deallocate(prior,input_state_vector,output_state_vector,work,work1D,preassim_state_vector)
deallocate(orginal_state_vector,obs_inc,state_inc,truth_values,obs_values,obs_errs)
deallocate(state_vector_prior_FO,state_vector_post_FO,input_qice,input_sice,input_qsnow)
deallocate(input_tsfcn,write_out_vars,state_var_table,obs_kind_table)
deallocate(file_array_input,file_array_output)

call finalize_mpi_utilities()

open(10,file='filter_done.txt',status='replace')
close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cice_cleanup(aicen,vicen,vsnon,input_qice,input_sice,input_qsno,input_tsfcn,NCAT,QICE,QSNO,aicen_org,vsnon_org)
  integer, intent(in) :: NCAT,QICE,QSNO
  real(r8), intent(inout) :: aicen(NCAT),input_qice(NCAT,QICE),input_sice(NCAT,QICE),input_qsno(NCAT,QSNO),input_tsfcn(NCAT)
  real(r8), intent(inout), optional :: vicen(NCAT),vsnon(NCAT)
  real(r8), intent(in) :: aicen_org(NCAT),vsnon_org(NCAT) 
 
  real(r8), parameter :: Tsmelt = 0._r8
  real(r8), parameter :: c1  = 1.0_r8
  real(r8), parameter :: &     
     phi_init = 0.75_r8, &   
     dSin0_frazil = 3.0_r8   

  real(r8), parameter :: sss = 34.7_r8 
  real(r8) :: aice,aice_temp
  real(r8) :: vice,vice_temp
  real(r8) :: vsno,vsno_temp
  real(r8) :: squeeze,cc1,cc2,cc3,x1
  real(r8) :: Si0new,qsno_hold,qi0new,Ti
  real(r8), allocatable ::  hin_max(:)
  real(r8), allocatable ::  hcat_midpoint(:)
  integer :: n
  
  !!!!!! SAVE ORGINAL COPY

  input_qice = min(0.0_r8,input_qice)
  input_sice = max(0.0_r8,input_sice)
  input_qsno = min(0.0_r8,input_qsno)
  aicen = min(1.0_r8,aicen)
  input_tsfcn = min(Tsmelt,input_tsfcn)
  !!!!!! 
  aice = sum(aicen) 
  vice = sum(vicen)
  vsno = sum(vsnon)
  !!!
  aicen = max(0.0_r8,aicen)
  vicen   = max(0.0_r8,vicen)
  vsnon = max(0.0_r8,vsnon)
  !!!
  aice_temp = sum(aicen)
  vice_temp = sum(vicen)
  vsno_temp = sum(vsnon)
  !!!
  if (aice<0.0_r8) then
    aicen(:) = 0.0_r8 
    vicen(:) = 0.0_r8
    vsnon(:) = 0.0_r8
  endif
  !!!!!!
  do n=1, NCAT
    if (aice_temp > 0._r8 .and. aice>0._r8) then
      aicen(n) = aicen(n) - (aice_temp-aice)*aicen(n)/aice_temp
    endif
    if (vice_temp > 0._r8 .and. vice>0._r8) then
      vicen(n) = vicen(n) - (vice_temp-vice)*vicen(n)/vice_temp
    endif
    if (vsno_temp > 0._r8 .and. vsno > 0._r8) then
      vsnon(n) = vsnon(n) - (vsno_temp-vsno)*vsnon(n)/vsno_temp
    endif
  enddo
  !!!!!!!!!
  if (aice>1.0_r8) then
    squeeze = 1.0_r8/aice
    aicen(:) = aicen(:)*squeeze
  endif
  !!!!!!!!
  where(aicen==-999) aicen = 0.0_r8  
  !!!!!!!!
  cc1 = 3._r8/real(Ncat,kind=r8)
  cc2 = 15.0_r8*cc1
  cc3 = 3._r8
  allocate( hin_max(0:Ncat) )
  allocate( hcat_midpoint(Ncat) )
  hin_max(0) = 0._r8

  do n = 1, NCAT
      x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
      hin_max(n) = hin_max(n-1) &
             + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
      hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!
  do n=1,NCAT
    if (aicen(n) > 0.0_r8 .and. aicen_org(n) > 0.0_r8) then
      if (vicen(n) == 0.0_r8) then
        vicen(n) = aicen(n)*hcat_midpoint(n)
      endif
    endif
    if (aicen(n) == 0.0_r8 .and. aicen_org(n) > 0.0_r8) then
      vicen(n) = 0.0_r8
      input_qice(n,:) = 0.0_r8
      input_sice(n,:) = 0.0_r8
      input_qsno(n,:) = 0.0_r8
      vsnon(n) = 0.0_r8
      input_tsfcn(n) = -1.8_r8
    else if (aicen(n)>0.0_r8 .and. aicen_org(n) == 0.0_r8) then
      if (vicen(n) == 0.0_r8) vicen(n) =  aicen(n) * hcat_midpoint(n)
      Si0new = sss - dSin0_frazil
      input_sice(n,:) = Si0new
      Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
      qi0new = enthalpy_mush(Ti, Si0new)
      input_qice(n,:) = qi0new
      if (vsnon(n) == 0.0_r8 .and. vsnon_org(n) > 0.0_r8) then
        input_qsno(n,:) = 0.0_r8
      else if (vsnon(n) > 0.0_r8 .and. vsnon_org(n) == 0.0_r8) then
        qsno_hold = snow_enthaply(Ti)
        input_qsno(n,:) = qsno_hold
      endif
      input_tsfcn(n) = Ti
    endif
    if (aicen(n) == 0.0_r8) then
      vicen(n) = 0.0_r8
      vsnon(n) = 0.0_r8
    endif
  enddo
end subroutine cice_cleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
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
    obsvar = const_err_val**2
  endif
  return
end subroutine set_sic_err
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine call_filter(prior,ens_size,prior_mean,prior_var,obs,obsvar,obs_inc,filter_type,obs_kind,obs_err_method)
  integer, intent(in) :: ens_size,filter_type,obs_err_method
  real(r8), intent(in) :: prior_mean,prior_var,obs,obsvar
  real(r8), intent(in) :: prior(ens_size)
  character(len=3), intent(in) :: obs_kind
  real(r8), intent(inout) :: obs_inc(ens_size)
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
  integer, intent(in) :: ens_size
  real(r8), intent(in) :: state(ens_size),obsspace(ens_size)
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
  real(r8) :: cdf_val
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
! Reads a single restart field
! author Chris Newman, LANL

      subroutine read_restart_field(nu,work,ndim,nx)
      integer, intent(in) :: &
         nu            , & ! unit number (not used for netcdf)
         ndim,&             ! number of dimensions
         nx
      real(r8), dimension(nx,ndim), intent(inout) :: &
         work              ! input array (real, 8-byte)

      ! local variables
      integer :: &
         n, i               ! loop indices

      real(r8), dimension(nx) :: &
         work2              ! input array (real, 8-byte)

      real (r8) :: &
        minw, maxw, sumw    ! diagnostics

      character(len=*), parameter :: subname='(read_restart_field)'
      
      do n = 1, ndim
         print*,n
         read(nu) (work2(i), i=1,nx)
         print*,n
         work(:,n) = work2(:)
      enddo

      minw = minval(work)
      maxw = maxval(work)
      sumw = sum(work)
      write(6,*) subname, minw, maxw, sumw

      end subroutine read_restart_field


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FO_SIC(ncid,scm_index,sic_ob,sic_err)
  integer, intent(in) :: &
         ncid,scm_index
  real(r8), intent(out) :: sic_ob,sic_err
  real(r8), allocatable :: work(:,:)
  integer :: start(2),count(2)
  
  status = nf90_inq_dimid(ncid,"ncat",dimid)
  status = nf90_inquire_dimension(ncid,dimid,len=NCAT)
  status = nf90_inq_dimid(ncid,"ni",dimid)
  status = nf90_inquire_dimension(ncid,dimid,len=NX)
  allocate(work(NX,NCAT))
  start(1) = 1
  count(1) = NX
  start(2) = 1
  count(2) = NCAT
  status = nf90_inq_varid(ncid, 'aicen', VarId)
  status = nf90_get_var(ncid, VarId, work, &
               start=start, &
               count=count )
  sic_ob = sum(work(scm_index,:))
  sic_err = 0.075_r8**2
  deallocate(work)
end subroutine FO_SIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FO_SIT(ncid,scm_index,sit_ob,sit_err)
  integer, intent(in) :: &
         ncid,scm_index
  real(r8), intent(out) :: sit_ob,sit_err
  real(r8) :: obs_sic,obs_siv
  real(r8), allocatable :: work_sic(:,:),work_sit(:,:)
  integer :: start(2),count(2)

  status = nf90_inq_dimid(ncid,"ncat",dimid)
  status = nf90_inquire_dimension(ncid,dimid,len=NCAT)
  status = nf90_inq_dimid(ncid,"ni",dimid)
  status = nf90_inquire_dimension(ncid,dimid,len=NX)
  allocate(work_sic(NX,NCAT))
  allocate(work_sit(NX,NCAT))
  start(1) = 1
  count(1) = NX
  start(2) = 1
  count(2) = NCAT
  status = nf90_inq_varid(ncid, 'aicen', VarId)
  status = nf90_get_var(ncid, VarId, work_sic, &
               start=start2, &
               count=count2 )
  obs_sic = sum(work_sic(scm_index,:))
  status = nf90_inq_varid(ncid, 'vicen', VarId)
  status = nf90_get_var(ncid, VarId, work_sit, &
               start=start, &
               count=count )
  obs_siv = sum(work_sit(scm_index,:))
  sit_ob = obs_siv/obs_sic
  sit_err = 0.1_r8**2
  deallocate(work_sic,work_sit)
end subroutine FO_SIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function enthalpy_mush(zTin, zSin) result(zqin)

    ! enthalpy of mush from mush temperature and bulk salinity

    real(r8), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real(r8) :: &
         zqin    ! ice layer enthalpy (J m-3)

    real(r8) :: &
         phi     ! ice liquid fraction

! from shr_const_mod.F90
    real(r8),parameter :: SHR_CONST_CPSW  = 3.996e3_R8   ! specific heat of sea water ~ J/kg/K
    real(R8),parameter :: SHR_CONST_CPICE = 2.11727e3_R8 ! specific heat of fresh ice ~ J/kg/K
    real(R8),parameter :: SHR_CONST_RHOSW = 1.026e3_R8   ! density of sea water ~ kg/m^3
    real(R8),parameter :: SHR_CONST_RHOICE= 0.917e3_R8   ! density of ice        ~ kg/m^3
    real(R8),parameter :: SHR_CONST_LATICE= 3.337e5_R8   ! latent heat of fusion ~ J/kg


! from cice/src/drivers/cesm/ice_constants.F90
    real(r8) :: cp_ocn, cp_ice, rhoi, rhow, Lfresh

    cp_ice    = SHR_CONST_CPICE  ! specific heat of fresh ice (J/kg/K)
    cp_ocn    = SHR_CONST_CPSW   ! specific heat of ocn    (J/kg/K)
    rhoi      = SHR_CONST_RHOICE ! density of ice (kg/m^3)
    rhow      = SHR_CONST_RHOSW  ! density of seawater (kg/m^3)
    Lfresh    = SHR_CONST_LATICE ! latent heat of melting of fresh ice (J/kg)

    phi = liquid_fraction(zTin, zSin)

    zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
           rhoi * cp_ice * zTin - (1._r8 - phi) * rhoi * Lfresh

  end function enthalpy_mush
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function liquid_fraction(zTin, zSin) result(phi)

    ! liquid fraction of mush from mush temperature and bulk salinity

    real(r8), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real(r8) :: &
         phi , & ! liquid fraction
         Sbr     ! brine salinity (ppt)

    real (r8), parameter :: puny = 1.0e-11_r8 ! cice/src/drivers/cesm/ice_constants.F90

    Sbr = max(liquidus_brine_salinity_mush(zTin),puny)
    phi = zSin / max(Sbr, zSin)

  end function liquid_fraction
  !=======================================================================
  function snow_enthaply(Ti) result(qsno)
    real(r8), intent(in) :: Ti

    real(r8),parameter :: rhos = 330.0_r8, &
                        Lfresh = 2.835e6_r8 - 2.501e6_r8, &
                        cp_ice = 2106._r8
    real(r8) :: qsno

    qsno = -rhos*(Lfresh - cp_ice*min(0.0_r8,Ti))
  end function snow_enthaply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function liquidus_brine_salinity_mush(zTin) result(Sbr)

    ! liquidus relation: equilibrium brine salinity as function of temperature
    ! based on empirical data from Assur (1958)

    real(r8), intent(in) :: &
         zTin         ! ice layer temperature (C)

    real(r8) :: &
         Sbr          ! ice brine salinity (ppt)

    real(r8) :: &
         t_high   , & ! mask for high temperature liquidus region
         lsubzero     ! mask for sub-zero temperatures

    !constant numbers from ice_constants.F90
    real(r8), parameter :: &
         c1      = 1.0_r8 , &
         c1000   = 1000_r8

    ! liquidus relation - higher temperature region
    real(r8), parameter :: &
         az1_liq = -18.48_r8 ,&
         bz1_liq =   0.0_r8

    ! liquidus relation - lower temperature region
    real(r8), parameter :: &
         az2_liq = -10.3085_r8,  &
         bz2_liq =  62.4_r8

    ! liquidus break
    real(r8), parameter :: &
         Tb_liq = -7.6362968855167352_r8

    ! basic liquidus relation constants
    real(r8), parameter :: &
         az1p_liq = az1_liq / c1000, &
         bz1p_liq = bz1_liq / c1000, &
         az2p_liq = az2_liq / c1000, &
         bz2p_liq = bz2_liq / c1000

    ! temperature to brine salinity
    real(r8), parameter :: &
       J1_liq = bz1_liq / az1_liq         , &
       K1_liq = c1 / c1000                , &
       L1_liq = (c1 + bz1p_liq) / az1_liq , &
       J2_liq = bz2_liq  / az2_liq        , &
       K2_liq = c1 / c1000                , &
       L2_liq = (c1 + bz2p_liq) / az2_liq

    t_high   = merge(1._r8, 0._r8, (zTin > Tb_liq))
    lsubzero = merge(1._r8, 0._r8, (zTin <= 1._r8))

    Sbr = ((zTin + J1_liq) / (K1_liq * zTin + L1_liq)) * t_high + &
          ((zTin + J2_liq) / (K2_liq * zTin + L2_liq)) * (1._r8 - t_high)

    Sbr = Sbr * lsubzero

  end function liquidus_brine_salinity_mush

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================

  function liquidus_temperature_mush(Sbr) result(zTin)

    ! liquidus relation: equilibrium temperature as function of brine salinity
    ! based on empirical data from Assur (1958)

    real(r8), intent(in) :: &
         Sbr    ! ice brine salinity (ppt)

    real(r8) :: &
         zTin   ! ice layer temperature (C)

    real(r8) :: &
         t_high ! mask for high temperature liquidus region

    ! liquidus break
    real(r8), parameter :: &
       Sb_liq =  123.66702800276086_r8    ! salinity of liquidus break

    ! constant numbers from ice_constants.F90
    real(r8), parameter :: &
         c1      = 1.0_r8 , &
         c1000   = 1000_r8

    ! liquidus relation - higher temperature region
    real(r8), parameter :: &
         az1_liq = -18.48_r8 ,&
         bz1_liq =   0.0_r8

    ! liquidus relation - lower temperature region
    real(r8), parameter :: &
         az2_liq = -10.3085_r8,  &
         bz2_liq =  62.4_r8

    ! basic liquidus relation constants
    real(r8), parameter :: &
         az1p_liq = az1_liq / c1000, &
         bz1p_liq = bz1_liq / c1000, &
         az2p_liq = az2_liq / c1000, &
         bz2p_liq = bz2_liq / c1000

  ! brine salinity to temperature
    real(r8), parameter :: &
       M1_liq = az1_liq            , &
       N1_liq = -az1p_liq          , &
       O1_liq = -bz1_liq / az1_liq , &
       M2_liq = az2_liq            , &
       N2_liq = -az2p_liq          , &
       O2_liq = -bz2_liq / az2_liq

    t_high = merge(1._r8, 0._r8, (Sbr <= Sb_liq))

    zTin = ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + &
          ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (1._r8 - t_high)

  end function liquidus_temperature_mush
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findindex1(array,dimlen,varname,findex)
  integer, intent(in) :: dimlen
  character(len=NF90_MAX_NAME), intent(in) :: array(dimlen)
  character(len=5), intent(in) :: varname
  integer, intent(out) :: findex
  integer :: i  

  do i=1, dimlen
    if (trim(array(i)) == trim(varname)) then
      findex = i
      EXIT
    else
      CYCLE
    endif
  enddo
end subroutine findindex1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program filter_cice_scm
