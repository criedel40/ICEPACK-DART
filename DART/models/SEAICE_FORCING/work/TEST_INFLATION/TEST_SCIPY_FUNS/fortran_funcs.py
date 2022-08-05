import numpy as np
import sys

def log1p(x):
  if abs(x) < 1.0e-5:
    p1 = x + 1.0
    careful = x*np.log(p1)/(p1-1.0)
    if p1 == 1:
      careful = x
    y = careful
  else:
    y = np.log(1.0 + x)
  return y
  
def gammal(xx):
  cof = np.array([76.18009172947146e0,-86.50532032941677e0, 24.01409824083091e0,-1.231739572450155e0,.1208650973866179e-2, -.5395239384953e-5])
  stp = 2.5066282746310005e0
  x = xx
  y=x
  
  tmp = x+5.5e0
  tmp = (x+0.5e0)*np.log(tmp)-tmp
 
  ser = 1.000000000190015e0
  for j in range(0,6):
    y = y + 1.0e0
    ser = ser + cof[j]/y
  gammln = tmp + np.log(stp*ser/x)
  return gammln


def betaln(a,b):
  return gammal(a) + gammal(b) - gammal(a+b)

def incbeta(a,b,x):
  if x < 0 or x >1:
    return np.nan
  TINY = 1.0e-30
  STOP = 1.0e-8
  #lbeta_ab = betaln(a,b)
  lbeta_ab = betaln(a,b)
  front = np.exp(np.log(x)*a+np.log(1.0-x)*b-lbeta_ab)/a
  f = 1.0;c=1.0;d=0.0
  for i in range(0,200):
    m = np.floor(i/2)
    if (i == 0):
      numerator = 1.0
    elif (i % 2 == 0):
      numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m))
    else:
      numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1))
   
    d = 1.0 + (numerator * d)
    if abs(d) < TINY: d = TINY
    d = 1.0/d

    c = 1.0 + (numerator/c)
    if abs(c) < TINY: c = TINY
    if i != 0:
      del cd
    cd = c*d
    
    f *= cd
    
    if abs(1.0-cd) < STOP:
      return front * (f-1.0)
      
  return np.nan

def betacdf(x,a,b):
  if a < 0 or b < 0:
    print('Incorrect beta parameters!!')
    return np.nan
  if x > (a+1.0)/(a+b+2.0):
    return (1.0 - incbeta(b,a,1.0-x))
  else:
    return incbeta(a,b,x)

def betapdf(x,a,b):
  if a == 1 and x ==0:
    y = b
  elif b == 1 and x == 1:
    y = a
  elif a<1 and x == 0:
    y = np.inf
  elif b<1 and x ==1:
    y = np.inf
    
  if a<=0:
    np.nan
  elif b<=0:
    np.nan
    
  loga = (a-1.0)*np.log(x)
  if x<0.1:
    logb = (b-1.0)*np.log1p(-x)
  else:
    logb = (b-1.0)*np.log(1-x)
  y = np.exp(loga+logb-betaln(a,b))
  return y

def betaincinv(p,a,b):
  if a<0 or b<0:
    print('Bad input beta parameters')
    return np.nan
  if p<0 or p>1:
    print('bad input quantile value')
    return np.nan
  if p == 0:
    return 0
  elif p == 1:
    return 1
  ## Using Newton's Method to find a root of gamcdf(X,A,B) = P
  q = a/(a+b)
  #print(q)
  maxiter = 500
  reltol = sys.float_info.epsilon**(3/4)
  #print(reltol)
  q = max(reltol,min(1.0-reltol,q))
  #print(q)
  F = betacdf(q,a,b)
  #print(F)
  dF = F - p
  #print(dF)
  for iter in range(1,maxiter+1):
    #print(iter)
    f = betapdf(q,a,b)
    #print(f)
    h = dF/f
    #print(h)
    qNew = max(q/10.0,min(1-(1-q)/10.,q-h))
    #print(qNew)
    if abs(h) <= reltol*q:
      q = qNew
      break
    dFold = dF
    F = betacdf(qNew,a,b)
    #print(F)
    for j in range(1,25):
      dF = F - p
      #print(abs(dF),abs(dFold))
      if abs(dF) < abs(dFold):
        break
      qNew = (q + qNew)/2
      F = betacdf(qNew,a,b)
    q = qNew
  badcdf = (abs(dF/F) > reltol**(2./3.))
  if iter>maxiter or badcdf:
    print('Alg. did not converage for a={0} b={1} c={2}'.format(a,b,p))
    return np.nan
  return q


