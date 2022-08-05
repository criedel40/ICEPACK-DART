!> Run test cases for beta dist

program run_test_cases

use types_mod, only : r8
use truncnorm_dist_functions, only : truncnorm_cdf,truncnorm_cdf_inv

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, nc_check, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

use default_model_mod,     only : init_time, init_conditions, adv_1step, &
                                  nc_write_model_vars

use random_seq_mod,    only : random_seq_type, init_random_seq, random_gaussian,random_beta, &
                              random_trunc_gaussian

use sort_mod, only : sort

use netcdf
implicit none

real(r8) :: output,obs_place
integer :: i,j
!real(r8) :: cdfs(nparams,nx),inv_cdf(nparams,nx)
integer :: ncid,iunit,io,fid,nXDimID,nParamDimID,rc,t
character(len=512) :: msgstring
character(len=19) :: out_file_control  = 'output_beta_stats.nc'
logical, dimension(2), parameter :: bound_flags = (/.true.,.true./)
real(r8),dimension(2), parameter :: bounds = (/-1.0_r8,1.0_r8/)
integer, parameter :: num = 10000
real(r8),dimension(num) :: truth,sorted_truth
type(random_seq_type) :: r_seq

real(r8),dimension(num) :: cdfs,inv_cdfs

! Initialize for error handler
call initialize_mpi_utilities('run_test_case')
call init_random_seq(r_seq, 1)

do i=1,num
  obs_place = -99.0
  do while(obs_place < bounds(1) .or. obs_place>bounds(2))
    obs_place = random_gaussian(r_seq, 0.0_r8, 1.0_r8)
  enddo
  truth(i) = obs_place
enddo
sorted_truth = sort(truth)


do i=1, num
    cdfs(i) = truncnorm_cdf(sorted_truth(i),0.0_r8,1.0_r8,bound_flags,bounds)
    inv_cdfs(i) = truncnorm_cdf_inv(cdfs(i),0.0_r8,1.0_r8,bound_flags,bounds)
end do


rc = nf90_create('output_beta_stats.nc', NF90_CLOBBER, fid)
call nc_check(rc, 'create_output_file', 'creating "'//trim(out_file_control)//'"')

msgstring = 'Output from beta alg cdf and inverse cdf'
call set_global_char_att(fid, 'Data output:', msgstring)

call setup_sec_dim(fid, 'nx', num, nXDimID)

call setup_sec_data_real1d(fid, 'cdfs', nXDimID)
call set_var_char_att(fid, 'cdfs', 'description','output from test case')

call setup_sec_data_real1d(fid, 'inv_cdfs', nXDimID)
call set_var_char_att(fid, 'inv_cdfs', 'description','output from test case')

call setup_sec_data_real1d(fid, 'x', nXDimID)
call set_var_char_att(fid, 'x', 'description','output from test case')
rc = nf90_enddef(fid)

do i=1,num
  call write_sec_data_real1d(fid, i, 'x', sorted_truth(i))
  call write_sec_data_real1d(fid, i, 'cdfs', cdfs(i))
  call write_sec_data_real1d(fid, i, 'inv_cdfs', inv_cdfs(i))
enddo
rc = nf90_close(fid)
call nc_check(rc, 'close_output_file', 'closing '//trim(out_file_control))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

end program run_test_cases