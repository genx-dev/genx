'''File: paratt_cython.py an implementation of the paratt algorithm in cython
should yield a speed increase > 2.
Programmed by: Matts Bjorck
Last changed: 2008 11 23
'''

from numpy import *
paratt_ext_built = False
debug = True
try:
    import genx.models.lib.paratt_ext as paratt_ext
    paratt_ext_built = True
except ImportError as e:
    import genx.models.lib.paratt_ext
    paratt_ext_built = False
    raise e

#if not paratt_ext_built and debug:
#    import subprocess
#    subprocess.run(['pwd'])
#    subprocess.run(['python3', '../build_ext.py', 'build_ext', '--inplace'])


def Refl(theta,lamda,n,d,sigma):
    theta = theta.astype(float64)
    lamda = float(lamda)
    n = n.astype(complex128)
    d = d.astype(float64)
    sigma = sigma.astype(float64)
    R = paratt_ext.refl(theta, lamda, n, d, sigma)
    return R

def ReflQ(Q,lamda,n,d,sigma):
    Q = Q.astype(float64)
    lamda = float(lamda)
    n = n.astype(complex128)
    d = d.astype(float64)
    sigma = sigma.astype(float64)
    R = paratt_ext.reflq(Q, lamda, n, d, sigma)
    return R

def Refl_nvary2(theta,lamda,n,d,sigma):
    theta=array(theta,dtype=float64)
    n = n.astype(complex128)
    d = d.astype(float64)
    sigma = sigma.astype(float64)
    lamda = lamda.astype(float64)
    R = paratt_ext.refl_nvary2(theta, lamda, n, d, sigma)
    return R

def Refl_nvary2_nosigma(theta, lamda, n, d):
    # Length of k-vector in vaccum
    #print n.shape, theta.shape, d.shape
    theta=array(theta,dtype=float64)
    n = n.astype(complex128)
    d = d.astype(float64)
    lamda = lamda.astype(float64)
    #print n.shape, theta.shape, d.shape
    R = paratt_ext.refl_nvary2_nosigma(theta, lamda, n, d)
    return R


if __name__=='__main__':
    import genx.models.lib.paratt as paratt
    import time
    import pylab as pl
    theta=arange(0,10,0.01)+1e-12
    q = 4*math.pi/1.54*sin(theta*pi/180.0)
    rep = 1000
    n = array([1-7.57e-6+1.73e-7j] + [1-2.24e-5+2.89e-6j,1-7.57e-6+1.73e-7j,1-2.24e-5+2.89e-6j,1-7.57e-6+1.73e-7j,1-2.24e-5+2.89e-6j,1-7.57e-6+1.73e-7j]*rep +[1])
    d = array([1] + [80,20,80,20,80,20]*rep + [1])*1.0
    sigma = array([0] + [0,0,0,0,0,0]*rep + [0])*1.0
#   print(n.shape)
    t1=time.clock()
    c1=paratt.Refl_nvary2(theta, 1.54*ones(theta.shape), n[:, newaxis]*ones(theta.shape), d,sigma)
#   c2=paratt.Refl_nvary2(theta, 1.54*ones(theta.shape), n[:, newaxis]*ones(theta.shape), d,sigma*0)
#   c3 = paratt.Refl(theta, 1.54, n, d, sigma)
    t2=time.clock()
    c4 = Refl_nvary2(theta, 1.54*ones(theta.shape), n[:, newaxis]*ones(theta.shape), d,sigma) 
#   c5 = Refl_nvary2_nosigma(theta, 1.54*ones(theta.shape), n[:, newaxis]*ones(theta.shape), d)
#   c6 = Refl(theta, 1.54, n, d, sigma)
    t3=time.clock()
    print(t2-t1,t3-t2)
    pl.plot(theta,log10(c1),'x',theta,log10(c4))
    pl.show()
