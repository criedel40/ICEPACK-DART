import numpy as np


def norm_cdf(x):
  a1 = 0.398942280444e+00
  a2 = 0.399903438504e+00
  a3 = 5.75885480458e+00
  a4 = 29.8213557808e+00
  a5 = 2.62433121679e+00
  a6 = 48.6959930692e+00
  a7 = 5.92885724438e+00
  b0 = 0.398942280385e+00
  b1 = 3.8052e-08
  b2 = 1.00000615302e+00
  b3 = 3.98064794e-04
  b4 = 1.98615381364e+00
  b5 = 0.151679116635e+00
  b6 = 5.29330324926e+00
  b7 = 4.8385912808e+00
  b8 = 15.1508972451e+00
  b9 = 0.742380924027e+00
  b10 = 30.789933034e+00
  b11 = 3.99019417011e+00
  
  
  if abs(x) < 1.28e0:
    y = 0.5e0 * x * x
    q = 0.5e0 - abs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7) ) ) )
  elif abs(x) < 12.7e0:
    y = 0.5e0 * x * x
    q = np.exp( -y ) * b0 / ( abs ( x ) - b1 \
         + b2  / ( abs ( x ) + b3 \
         + b4  / ( abs ( x ) - b5 \
         + b6  / ( abs ( x ) + b7 \
         - b8  / ( abs ( x ) + b9 \
         + b10 / ( abs ( x ) + b11 ) )
      ) ) ) )
  else:
    q = 0.0e0
    
  if x < 0.0e0:
    cdf = q
  else:
    cdf = 1.0e0 - q
    
  return cdf

def horner(coeffs, x):
    acc = 0
    for c in reversed(coeffs):
        acc = acc * x + c
    return acc
  


def norm_cdf_inv(p):
  const1 = 0.180625e+00
  const2 = 1.6e+00
  split1 = 0.425e+00
  split2 = 5.0e+00
  
  a = np.array([3.3871328727963666080e+00,
     1.3314166789178437745e+02,
     1.9715909503065514427e+03,
     1.3731693765509461125e+04,
     4.5921953931549871457e+04,
     6.7265770927008700853e+04,
     3.3430575583588128105e+04,
     2.5090809287301226727e+03])
  b = np.array([1.0e+00,
     4.2313330701600911252e+01,
     6.8718700749205790830e+02,
     5.3941960214247511077e+03,
     2.1213794301586595867e+04,
     3.9307895800092710610e+04,
     2.8729085735721942674e+04,
     5.2264952788528545610e+03])
  c = np.array([1.42343711074968357734e+00,
     4.63033784615654529590e+00,
     5.76949722146069140550e+00,
     3.64784832476320460504e+00,
     1.27045825245236838258e+00,
     2.41780725177450611770e-01,
     2.27238449892691845833e-02,
     7.74545014278341407640e-04])
  d = np.array([1.0e+00,
     2.05319162663775882187e+00,
     1.67638483018380384940e+00,
     6.89767334985100004550e-01,
     1.48103976427480074590e-01,
     1.51986665636164571966e-02,
     5.47593808499534494600e-04,
     1.05075007164441684324e-09])
  e = np.array([6.65790464350110377720e+00,
     5.46378491116411436990e+00,
     1.78482653991729133580e+00,
     2.96560571828504891230e-01,
     2.65321895265761230930e-02,
     1.24266094738807843860e-03,
     2.71155556874348757815e-05,
     2.01033439929228813265e-07])
  f = np.array([1.0e+00,
     5.99832206555887937690e-01,
     1.36929880922735805310e-01,
     1.48753612908506148525e-02,
     7.86869131145613259100e-04,
     1.84631831751005468180e-05,
     1.42151175831644588870e-07,
     2.04426310338993978564e-15])

  if p < 0.0e0:
    print('Can not have a negative qunatile...STOP!')
    sys.exit()
  if p > 1.0e0:
    print('Can not have a quantile greater than 1...STOP!')
    sys.exit()
    
  q = p - 0.5e0
  
  if abs(q) < split1:
    r = const1 - q * q
    x = q * horner(a,r)/horner(b,r)
  else:
    if q < 0.0e0:
      r = p
    else:
      r = 1.0e0 - p
    if r < 0.0e0:
      print('Bad value when computing inverse cdf...STOP')
      sys.exit()
    else:
      r = np.sqrt(-np.log(r))
      if r < split2:
        r = r - const2
        x = horner(c,r)/horner(d,r)
      else:
        r = r - split2
        x = horner(e,r)/horner(f,r)
    if q < 0.0e0:
      x = -x
  return x
     
def truncnorm_cdf(x,mu,sigma,a,b):
  if np.isinf(b):
    if x < a:
      cdf = 0.0e0
    else:
      alpha = (a-mu)/sigma
      xi = (x-mu)/sigma
      alpha_cdf = norm_cdf(alpha)
      xi_cdf = norm_cdf(xi)
      cdf = (xi_cdf - alpha_cdf)/(1.0e0 - alpha_cdf)
  elif np.isinf(a):
    if x < b:
      beta = (b-mu)/sigma
      xi = (x-mu)/sigma
      beta_cdf = norm_cdf(beta)
      xi_cdf = norm_cdf(xi)
      cdf = xi_cdf/beta_cdf
    else:
      cdf = 1.0e0
  else:
    if x < a:
      cdf = 0.0e0
    elif x < b:
      alpha = (a-mu)/sigma
      beta = (b-mu)/sigma
      xi = (x-mu)/sigma
      alpha_cdf = norm_cdf(alpha)
      beta_cdf = norm_cdf(beta)
      xi_cdf = norm_cdf(xi)
      cdf = (xi_cdf - alpha_cdf)/(beta_cdf - alpha_cdf)
    else:
      cdf = 1.0e0
      
  return cdf

def truncnorm_cdf_inv(p,mu,sigma,a,b):
  if p < 0.0 or p > 1.0:
    print('Invalid value for the cdf(p)')
    print('0<p<1.0')
    sys.exit()
  if np.isinf(b):
    if p == 0.0:
      x = a
      return
    alpha = (a-mu)/sigma
    alpha_cdf = norm_cdf(alpha)
    xi_cdf = (1.0e0 - alpha_cdf) * p+alpha_cdf
    xi = norm_cdf_inv(xi_cdf)
    x = mu + sigma*xi
  elif np.isinf(a):
    if p == 1.0:
      return b
    beta = (b - mu)/sigma
    beta_cdf = norm_cdf(beta)
    xi_cdf = beta_cdf*p
    xi = norm_cdf_inv(xi_cdf)
    x = mu + sigma*xi
  else:
    if p == 0.0:
      return a
    elif p == 1.0:
      return b
    alpha = (a-mu)/sigma
    beta = (b-mu)/sigma
    alpha_cdf = norm_cdf(alpha)
    beta_cdf = norm_cdf(beta)
    xi_cdf = ( beta_cdf - alpha_cdf ) * p + alpha_cdf
    xi = norm_cdf_inv(xi_cdf)
    x = mu + sigma*xi
  return x
