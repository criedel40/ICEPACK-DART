program test_bounded_rhf

! Test specific cases for obs_increment_bounded_norm_rhf
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use types_mod,         only : r8
use assim_tools_mod,   only : obs_increment_bounded_norm_rhf, &
                              obs_increment_rank_histogram
use random_seq_mod,    only : random_seq_type, init_random_seq, random_gaussian

implicit none

real(r8), parameter :: truth = 0.9899_r8
integer,  parameter :: ens_size = 80
integer,  parameter :: num_steps = 5000
integer,  parameter :: bounds_case = 4         ! Case 4 is doubly bounded
real(r8), parameter :: bound(2) = (/0.0_r8, 1.0_r8/)

real(r8)              :: ens(ens_size), obs_inc(ens_size) 
real(r8)              :: obs, obs_var, obs_std, prior_mean, prior_var
logical               :: is_bounded(4, 2)
type(random_seq_type) :: r_seq
integer               :: i, j

! Initialize for error handler
call initialize_mpi_utilities('test_bounded_rhf')

! Initialize a repeating random sequence, seed fixed at 1
call init_random_seq(r_seq, 1)

! The 4 possible bounding cases: none, left, right, both
is_bounded(1:4, 1) = (/.false., .true., .false., .true./)
is_bounded(1:4, 2) = (/.false., .false., .true., .true./)

! Set the observation error variance
if(truth > 0.1_r8) then 
   obs_std = 0.15 * truth
else
   ! NOTE: Go back to 0.15 for small truth to compare 0 and 1 truth cases
   obs_std = 0.15;
endif
obs_var = obs_std**2

! Initial ensemble is uninformative over [0 1]
do i = 1, ens_size
   ens(i) = (i * 1.0_r8) / (ens_size + 1.0_r8)
end do

! Loop through a sequence of observations
do i = 1, num_steps

   ! Generate truncated normal observation
   obs = 99.0_r8
   do while(obs >= 1.0_r8 .or. obs <= 0.0_r8)
      obs = random_gaussian(r_seq, truth, sqrt(obs_var))
   end do

   ! Compute prior variance and mean from sample
   prior_mean = sum(ens) / ens_size
   prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

   ! Naive output of ensemble mean and variance
   write(41, *) prior_mean, prior_var, obs
   
   ! Output of entire ensemble
   do j = 1, ens_size
      if (j == 1) then
        write(42, *) '##############'
      Endif
      write(42, *) j, ens(j)
   end do
   !print*,ens
   !print*,obs
   ! Use the bounded rhf
   !call obs_increment_bounded_norm_rhf(ens, ens_size, prior_var, obs, obs_var, &
   !   obs_inc, is_bounded(bounds_case, :), bound)

   ! Baseline check for original rhf
   call obs_increment_rank_histogram(ens, ens_size, prior_var, obs, obs_var, &
      obs_inc)

   ! Update the ensemble
   ens = ens + obs_inc
   !print*,ens

end do

! Metadata for entire ensemble
write(42, *) ens_size, truth

call finalize_mpi_utilities()

end program test_bounded_rhf
