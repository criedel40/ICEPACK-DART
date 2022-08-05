!> Functions to compute did beta distribution quantities 
!> required additional functions are included so no outside code is needed
!> Currently working beta distributions functions
!> betacdf => compute CDFs for beta distribution using incomplete beta function
!> betapdf => compute PDFs for beta distrubution
!> betainv => compute inverse CDFs for beta distribution


module beta_dist_functions

use types_mod, only: digits12, i8, r8, PI
use  utilities_mod,       only : error_handler, E_ERR, E_MSG

implicit none
private

public :: log1p, &
          gammal, &
          betaln, &
          incbeta, &
          betacdf, &
          betapdf, &
          betaincinv
character(len=512)     :: msgstring
character(len=*), parameter :: source = 'beta_dist_functions.f90'

contains

function log1p(x)
  real(r8), intent(in) :: x
  real(r8) :: log1p
  real(r8) :: careful,p1
  if (abs(x) < 1.0d-5) then
    p1 = x + 1.0_r8
    careful = x*log(p1)/(p1-1.0_r8)
    if (p1 == 1.0_r8) then
      log1p = x
    else
      log1p = careful
    endif
  else
    log1p = log(1.0_r8 + x)
  endif
end function log1p

function gammal(xx)
  real(r8), intent(in) :: xx
  real(r8) :: gammal
  real(r8), parameter :: cov(6) = (/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, -.5395239384953d-5/)
  real(r8), parameter :: stp = 2.5066282746310005d0
  real(r8) :: x,y,tmp,ser
  integer :: j

  x = xx
  y = x

  tmp = x+5.5d0
  tmp = (x+0.5d0)*log(tmp) - tmp
  ser = 1.000000000190015d0
  do j=1, 6
    y = y + 1.0_r8
    ser = ser + cov(j)/y  
  end do 
  gammal = tmp + log(stp*ser/x)
end function gammal

subroutine betaln(a,b,output)
  real(r8), intent(in) :: a,b
  real(r8), intent(out) :: output

  output = gammal(a) + gammal(b) - gammal(a+b)

end subroutine betaln

function incbeta(a,b,x)
  real(r8), intent(in) :: a,b,x
  real(r8) :: incbeta
  real(r8), parameter :: TINY = 1.0e-30
  real(r8), parameter :: STOP = 1.0e-8
  real(r8) :: lbeta_ab,front,f,c,d,numerator,cd
  integer :: m,i
  integer, parameter :: bot = 2

  if (x < 0 .or. x > 1) then
    msgstring = 'Input value for x is not between 0 - 1'
    call error_handler(E_ERR, 'incbeta', msgstring, source)
  endif

  call betaln(a,b,lbeta_ab)    
  front = exp(log(x)*a + log(1.0_r8-x)*b - lbeta_ab) / a
  f = 1.0_r8;c=1.0_r8;d=0.0_r8
  do i=0, 200
    m = floor(i/2.0_r8)
    if (i == 0) then
      numerator = 1.0_r8     
    else if (mod(i,2) == 0) then
      numerator = (m*(b-m)*x)/((a+2.0_r8*m-1.0)*(a+2.0_r8*m))
    else
      numerator = -((a+m)*(a+b+m)*x)/((a+2.0_r8*m)*(a+2.0_r8*m+1))
    end if   

    d = 1.0_r8 + (numerator * d)
    if (abs(d) < TINY) d = TINY
    d = 1.0_r8/d

    c = 1.0_r8 + (numerator/c)
    if (abs(c) < TINY) c = TINY

    cd = c*d
    
    f = cd*f

    if (abs(1.0_r8 - cd) < STOP) then
      incbeta = front * (f-1.0_r8)
      return
    end if
  end do
  msgstring = 'Alg. did not converge'
  call error_handler(E_ERR, 'incbeta', msgstring, source)
end function incbeta


subroutine betacdf(x,a,b,output)
  real(r8),intent(in) :: x,a,b
  real(r8),intent(out) :: output

  if (a < 0 .or. b < 0) then
    msgstring = 'Bad input beta parameters'
    call error_handler(E_ERR, 'betacdf', msgstring, source)
  endif

  if (x > (a+1.0_r8)/(a+b+2.0_r8)) then
    output = (1.0_r8 - incbeta(b,a,1.0_r8-x))
    return
  else
    output = incbeta(a,b,x)
    return
  endif
end subroutine betacdf

subroutine betapdf(x,a,b,output)
  real(r8),intent(in) :: x,a,b
  real(r8),intent(out) :: output

  real(r8) :: loga,logb,beta_ln

  if (a == 1 .and. x == 0) then
    output = b
    return
  else if (b == 1.0_r8 .and. x == 1.0_r8) then
    output = a
    return
  else if (a < 1.0_r8 .and. x == 0.0_r8) then
    output = -9999._r8
    return
  else if (b < 1.0_r8 .and. x == 1.0_r8) then
    output = -9999._r8
    return
  endif

  if (a<= 0.0 .or. b <= 0.0) then
    msgstring = 'Bad input beta parameters'
    call error_handler(E_ERR, 'incbeta', msgstring, source)
  endif

  loga = (a-1.0_r8)*log(x)

   if (x < 0.1_r8) then
     logb = (b-1.0_r8)*log1p(-x)
   else
     logb = (b-1.0_r8)*log(1.0_r8-x)
   endif
   call betaln(a,b,beta_ln)
   output = exp(loga+logb - beta_ln)
   return
end subroutine betapdf

subroutine betaincinv(p,a,b,output)
  real(r8), intent(in) :: p,a,b
  real(r8), intent(out) :: output

  real(r8) :: q,reltol,dF,F,ff,h,qNew,dFold
  integer :: maxiter,iter,j

  if (a < 0.0_r8 .or. b < 0.0_r8) then
    msgstring = 'Bad input beta parameters'
    call error_handler(E_ERR, 'betaincinv', msgstring, source)
  endif

  if (p < 0.0_r8 .or. p > 1.0_r8) then
    msgstring = 'Bad input quantile value'
    call error_handler(E_ERR, 'betaincinv', msgstring, source)
  endif

  if (p == 0.0_r8) then
    output = 0.0_r8
    return
  else if (p == 1.0_r8) then
    output = 1.0_r8
    return
  endif
  
  !Using Newton's Method to find a root of gamcdf(X,A,B) = P
  q = a/(a+b)
  maxiter = 500
  reltol = (EPSILON(q))**(3./4.)
  q = max(reltol, min(1.0_r8-reltol,q))
  
  call betacdf(q,a,b,f)
  
  dF = f - p

  do iter=1, maxiter
    call betapdf(q,a,b,ff)
    h = dF/ff

    qNew = max(q/10._r8,min(1.0_r8-(1.0_r8-q)/10._r8,q-h))

    if (abs(h) <= reltol*q) then
      output = qNew
      return
    endif
    
    dFold = dF
    call betacdf(qNew,a,b,f)
    do j=1, 25
      dF = f - p
      if (abs(dF) < abs(dFold)) then
        EXIT
      endif
      qNew = (q + qNew)/2.0_r8
      call betacdf(qNew,a,b,f)
    end do
    q = qNew
  end do
  output = qNew
  return
end subroutine betaincinv













end module beta_dist_functions