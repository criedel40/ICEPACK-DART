!> Run test cases for beta dist

program run_test_cases

use types_mod, only : r8
use beta_dist_functions, only : log1p,betacdf,betaln,gammal,betaincinv

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, nc_check, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

use default_model_mod,     only : init_time, init_conditions, adv_1step, &
                                  nc_write_model_vars

use netcdf
implicit none

real(r8) :: output
integer :: i,j
integer, parameter :: nx = 1001
Integer, parameter :: nparams = 5 
real(r8) :: x(nx) = (/(i, i=0,1000, 1)/)*0.001
real(r8) :: a_array(nparams) = (/0.5_r8,1.0_r8,8.0_r8,3.0_r8,5.0_r8/)
real(r8) :: b_array(nparams) = (/5.0_r8,3.0_r8,8.0_r8,1.0_r8,0.5_r8/)
real(r8) :: cdfs(nparams,nx),inv_cdf(nparams,nx)
integer :: ncid,iunit,io,fid,nXDimID,nParamDimID,rc,t
character(len=512) :: msgstring
character(len=19) :: out_file_control  = 'output_beta_stats.nc'

! Initialize for error handler
call initialize_mpi_utilities('test_bounded_rhf')


do i=1, nparams
  do j=1, nx
    call betacdf(x(j),a_array(i), b_array(i), cdfs(i,j))
    call betaincinv(x(j),a_array(i), b_array(i), inv_cdf(i,j))
  end do
end do


rc = nf90_create('output_beta_stats.nc', NF90_CLOBBER, fid)
call nc_check(rc, 'create_output_file', 'creating "'//trim(out_file_control)//'"')

msgstring = 'Output from beta alg cdf and inverse cdf'
call set_global_char_att(fid, 'Data output:', msgstring)

call setup_sec_dim(fid, 'nparams', nparams, nParamDimID)
call setup_sec_dim(fid, 'nx', nx, nXDimID)

call setup_sec_data_real(fid, 'cdfs', nXDimID, nParamDimID)
call set_var_char_att(fid, 'cdfs', 'description','output from test case')

call setup_sec_data_real(fid, 'inv_cdfs', nXDimID, nParamDimID)
call set_var_char_att(fid, 'inv_cdfs', 'description','output from test case')

call setup_sec_data_real1d(fid, 'x', nXDimID)
call set_var_char_att(fid, 'x', 'description','output from test case')
rc = nf90_enddef(fid)

do t=1,nparams
  call write_sec_data_real(fid, t, 'cdfs', cdfs(t,:))
  call write_sec_data_real(fid, t, 'inv_cdfs', inv_cdf(t,:))
end do
do i=1,nx
  call write_sec_data_real1d(fid, i, 'x', x(i))
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