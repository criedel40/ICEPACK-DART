!> Functions to compute did beta distribution quantities 
!> required additional functions are included so no outside code is needed
!> Currently working beta distributions functions
!> betacdf => compute CDFs for beta distribution using incomplete beta function
!> betapdf => compute PDFs for beta distrubution
!> betainv => compute inverse CDFs for beta distribution


module truncnorm_dist_functions

use types_mod, only: digits12, i8, r8, PI
use  utilities_mod,       only : error_handler, E_ERR, E_MSG

implicit none
private

public :: truncnorm_cdf, &
          truncnorm_cdf_inv

character(len=512)     :: msgstring
character(len=*), parameter :: source = 'truncnorm_dist_functions.f90'

contains



function norm_cdf(x)
  real(r8), intent(in) :: x
  real(r8) :: norm_cdf
  real(r8), parameter :: a1 = 0.398942280444D+00
  real(r8), parameter :: a2 = 0.399903438504D+00
  real(r8), parameter :: a3 = 5.75885480458D+00
  real(r8), parameter :: a4 = 29.8213557808D+00
  real(r8), parameter :: a5 = 2.62433121679D+00
  real(r8), parameter :: a6 = 48.6959930692D+00
  real(r8), parameter :: a7 = 5.92885724438D+00
  real(r8), parameter :: b0 = 0.398942280385D+00
  real(r8), parameter :: b1 = 3.8052D-08
  real(r8), parameter :: b2 = 1.00000615302D+00
  real(r8), parameter :: b3 = 3.98064794D-04
  real(r8), parameter :: b4 = 1.98615381364D+00
  real(r8), parameter :: b5 = 0.151679116635D+00
  real(r8), parameter :: b6 = 5.29330324926D+00
  real(r8), parameter :: b7 = 4.8385912808D+00
  real(r8), parameter :: b8 = 15.1508972451D+00
  real(r8), parameter :: b9 = 0.742380924027D+00
  real(r8), parameter :: b10 = 30.789933034D+00
  real(r8), parameter :: b11 = 3.99019417011D+00

  real(r8) :: y,q

  if (abs(x) < 1.28_r8) then
    y = 0.5_r8 * x * x
    q = 0.5_r8 - abs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7) ) ) )
  else if (abs(x) < 12.7_r8) then
    y = 0.5_r8 * x * x
    q = exp( -y ) * b0 / ( abs ( x ) - b1 + b2  / ( abs ( x ) + b3 + b4  / ( abs ( x ) - b5 + b6  / ( abs ( x ) + b7 - b8  / ( abs ( x ) + b9 + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
  else
    q = 0.0_r8
  endif
  
  if (x < 0.0_r8) then
    norm_cdf = q
  else
    norm_cdf = 1.0_r8 - q
  endif

end function norm_cdf

function poly_solver(coeffs, x)
    real(r8), dimension(:), intent (in) :: coeffs
    real(r8), intent (in) :: x
    real(r8) :: res,poly_solver
    integer :: i
 
    res = 0.0
    do i = size (coeffs), 1, -1
      res = res * x + coeffs (i)
    end do
    poly_solver = res
end function poly_solver

function norm_cdf_inv(p)
  real(r8), intent(in) :: p
  real(r8) :: norm_cdf_inv
  real(r8), parameter :: const1 = 0.180625D+00
  real(r8), parameter :: const2 = 1.6D+00
  real(r8), parameter :: split1 = 0.425D+00
  real(r8), parameter :: split2 = 5.0D+00
  real(r8), dimension(8), parameter :: a = (/3.3871328727963666080D+00, 1.3314166789178437745D+02, &
                                            1.9715909503065514427D+03,1.3731693765509461125D+04, &
                                            4.5921953931549871457D+04,6.7265770927008700853D+04, &   
                                            3.3430575583588128105D+04,2.5090809287301226727D+03/)
  real(r8), dimension(8), parameter :: b = (/1.0D+00,4.2313330701600911252D+01, &
                                            6.8718700749205790830D+02,5.3941960214247511077D+03, &
                                            2.1213794301586595867D+04,3.9307895800092710610D+04, &
                                            2.8729085735721942674D+04,5.2264952788528545610D+03 /)
  real(r8), dimension(8), parameter :: c = (/1.42343711074968357734D+00,4.63033784615654529590D+00, &
                                             5.76949722146069140550D+00,3.64784832476320460504D+00, &
                                             1.27045825245236838258D+00,2.41780725177450611770D-01, &
                                             2.27238449892691845833D-02,7.74545014278341407640D-04 /)
  real(r8), dimension(8), parameter :: d = (/1.0D+00,2.05319162663775882187D+00, &
                                             1.67638483018380384940D+00,6.89767334985100004550D-01, &
                                             1.48103976427480074590D-01,1.51986665636164571966D-02, &
                                             5.47593808499534494600D-04,1.05075007164441684324D-09 /)
  real(r8), dimension(8), parameter :: e = (/6.65790464350110377720D+00,5.46378491116411436990D+00, &
                                             1.78482653991729133580D+00,2.96560571828504891230D-01, &
                                             2.65321895265761230930D-02,1.24266094738807843860D-03, &
                                             2.71155556874348757815D-05,2.01033439929228813265D-07 /)
  real(r8), dimension(8), parameter :: f = (/1.0D+00,5.99832206555887937690D-01, &
                                             1.36929880922735805310D-01,1.48753612908506148525D-02, &
                                             7.86869131145613259100D-04,1.84631831751005468180D-05, &
                                             1.42151175831644588870D-07,2.04426310338993978564D-15 /)
  real(r8) :: q,r,x

  if (p < 0.0_r8 .or. p > 1.0_r8) then
    msgstring = 'Input quantile has to be between 0 and 1'
    call error_handler(E_ERR, 'norm_cdf_inv', msgstring, source)
  endif

  q = p - 0.5_r8
  if (abs(q) < split1) then
    r = const1 - q*q
    x = q*poly_solver(a,r)/poly_solver(b,r)
  else
    if (q < 0.0_r8) then
      r = p
    else
      r = 1.0_r8 - p
    endif
    if (r < 0.0_r8) then
      msgstring = 'Bad value (r) when computing inverse CDF'
      call error_handler(E_ERR, 'norm_cdf_inv', msgstring, source)
    else
      r = sqrt(-log(r))
      if (r<split2) then
        r = r - const2
        x = poly_solver(c,r)/poly_solver(d,r)
      else
        r = r - split2
        x = poly_solver(e,r)/poly_solver(f,r)
      endif
    endif
    if (q < 0.0_r8) then
      x = -x
    endif
  endif
  norm_cdf_inv = x
end function norm_cdf_inv

function truncnorm_cdf(x,mu,sigma,bound_flags,bounds)
  real(r8), intent(in) :: x,mu,sigma
  logical, dimension(2), intent(in) :: bound_flags
  real(r8), dimension(2), intent(in) :: bounds
  real(r8) :: truncnorm_cdf

  real(r8) :: cdf,alpha,beta,alpha_cdf,beta_cdf,xi,xi_cdf

  if (bound_flags(1) .and. .not. bound_flags(2)) then
    if (x<bounds(1)) then
      cdf = 0.0_r8
    else
      alpha = (bounds(1) - mu)/sigma
      xi = (x-mu)/sigma
      alpha_cdf = norm_cdf(alpha)
      xi_cdf = norm_cdf(xi)
      cdf = (xi_cdf - alpha_cdf)/(1.0_r8 - alpha_cdf)
    endif
  else if (.not. bound_flags(1) .and. bound_flags(2)) then
    if (x<bounds(2)) then
      beta = (bounds(2) - mu)/sigma
      xi = (x-mu)/sigma
      beta_cdf = norm_cdf(beta)
      xi_cdf = norm_cdf(xi)
      cdf = xi_cdf/beta_cdf
    else
      cdf = 1.0_r8
    endif
  else if (bound_flags(1) .and. bound_flags(2)) then
    if (x < bounds(1)) then
      cdf = 0.0_r8
    else if (x < bounds(2)) then
      alpha = (bounds(1) - mu)/sigma
      beta = (bounds(2) - mu)/sigma
      xi = (x-mu)/sigma
      alpha_cdf = norm_cdf(alpha)
      beta_cdf = norm_cdf(beta)
      xi_cdf = norm_cdf(xi)
      cdf = (xi_cdf - alpha_cdf)/(beta_cdf - alpha_cdf)
    else
      cdf = 1.0_r8
    endif
  else
    msgstring = 'No bounds were provide so use normal function for norm_cdf'
    call error_handler(E_ERR, 'truncnorm_cdf', msgstring, source)
  endif
  truncnorm_cdf = cdf
end function truncnorm_cdf

function truncnorm_cdf_inv(p,mu,sigma,bound_flags,bounds)
  real(r8), intent(in) :: p,mu,sigma
  logical, dimension(2), intent(in) :: bound_flags
  real(r8), dimension(2), intent(in) :: bounds
  real(r8) :: truncnorm_cdf_inv
  real(r8) :: cdf,alpha,beta,alpha_cdf,beta_cdf,xi,xi_cdf,x

  if (p < 0.0_r8 .or. p > 1.0_r8) then
    msgstring = 'Input quantile has to be between 0 and 1'
    call error_handler(E_ERR, 'norm_cdf_inv', msgstring, source)
  endif

  if (bound_flags(1) .and. .not. bound_flags(2)) then
    alpha = (bounds(1) - mu)/sigma
    alpha_cdf = norm_cdf(alpha)
    xi_cdf = (1.0_r8 - alpha_cdf) * p+alpha_cdf
    xi = norm_cdf_inv(xi_cdf)
    x = mu + sigma*xi
  else if (.not. bound_flags(1) .and. bound_flags(2)) then
    beta = (bounds(2) - mu)/sigma
    beta_cdf = norm_cdf(beta)
    xi_cdf = beta_cdf*p
    xi = norm_cdf_inv(xi_cdf)
    x = mu + sigma*xi
  else if (bound_flags(1) .and. bound_flags(2)) then
    alpha = (bounds(1) - mu)/sigma
    beta = (bounds(2) - mu)/sigma
    alpha_cdf = norm_cdf(alpha)
    beta_cdf = norm_cdf(beta)
    xi_cdf = (beta_cdf - alpha_cdf) * p+alpha_cdf
    xi = norm_cdf_inv(xi_cdf)
    x = mu + sigma*xi
  else
    msgstring = 'No bounds were provide so use normal function for norm_cdf_inv'
    call error_handler(E_ERR, 'truncnorm_cdf_inv', msgstring, source)
  endif
  truncnorm_cdf_inv = x
end function truncnorm_cdf_inv

end module truncnorm_dist_functions