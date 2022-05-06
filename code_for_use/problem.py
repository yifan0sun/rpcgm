import numpy as np
from PIL import Image


#Problems encoded here

class LeastSquares():
    def __init__(self,m,n,x,eta = .1,seed = None):
        self.m = m
        self.n = n
        if seed is not None:
            np.random.seed(seed+10)
        self.A =  np.random.randn(m,n)/(n+0.)
        self.x_true = x
        self.b = np.dot(self.A,self.x_true) + eta * np.random.randn(m)
        self.L = self.get_L()
        if m > n:
            self.AA = np.dot(self.A.T,self.A)
            self.Ab = np.dot(self.A.T,self.b)
            self.bb = np.dot(self.b, self.b)
        print m,n
        
    def get_L(self):
        A = self.A + 0.
            
        #L smoothness w.r.t. 1 norm     
        Acolnorm = np.sum(A*A,0)
        return  np.max(Acolnorm)/self.m


    def get_grad(self,x):
        if self.m > self.n:
            return (np.dot(self.AA,x)-self.Ab)/self.m
        return np.dot(self.A.T,np.dot(self.A,x)-self.b)/self.m

    def get_obj(self,x):
        if self.m > self.n:
            z = self.bb + np.dot(x.T,np.dot(self.AA,x)) - 2.*np.dot(x,self.Ab)
            z = z / self.m/2.
            return z
        z = np.dot(self.A,x)-self.b
        return np.dot(z,z)/2./self.m
        
        
    
class Sensing(LeastSquares):
    def __init__(self,m,n,s,eta = .1,seed = None):
        np.random.seed(seed)
        self.x_true = np.zeros(n)
        idx = np.argsort(np.random.rand(n) )[:s]
        self.x_true[idx] = np.random.randn(s)
        LeastSquares.__init__(self,m,n,self.x_true,eta = .1,seed = None)
        
class Denoising(LeastSquares):
    def __init__(self,m,n,s,eta = .1,seed = None):
        if seed is not None: np.random.seed(seed)
        u = np.zeros(n)
        idx = np.argsort(np.random.rand(n-1) )[:s]
        u[idx+1] = np.random.randn(s)
        u[0] = np.random.randn()
        self.x_true = np.cumsum(u)
        LeastSquares.__init__(self,m,n,self.x_true,eta = .1,seed = None)
     
    
class GroupSparsity(LeastSquares):
    def __init__(self,m,n,s,groups,eta = .1,seed = None):
        if seed is not None: np.random.seed(seed)
        self.x_true = np.zeros(n)
        idx = np.argsort(np.random.rand(len(groups)) )[:s]
        for i in idx:
            self.x_true[groups[i]] = np.random.randn(len(groups[i]))
        self.groups = groups
        self.group_truth = np.zeros(len(groups))
        for i in idx:
            self.group_truth[i] = 1.
        
        LeastSquares.__init__(self,m,n,self.x_true,eta = .1,seed = None)

class SimpleDenoising():
    def __init__(self,im, noiselevel, seed):
        self.noiselevel = noiselevel
        self.seed = seed
        self.im = im
        im = np.asarray(im,np.float)
        im = np.mean(im,2)

        self.imshape = im.shape
        
        self.n0 = np.prod(self.imshape)
        self.n = 2*self.n0
        self.x0 = np.reshape(im,self.n0)
        self.normfact = np.mean(self.x0)#np.linalg.norm(self.x0,2)
        
        np.random.seed(seed)
        self.x0 = self.x0 + 10.*np.random.randn(self.n0)*np.less(np.random.rand(self.n0),noiselevel) 
        self.x0 = self.x0 / self.normfact
        
        self.x0 = np.reshape(self.x0,im.shape)
        self.x0 = np.hstack([ np.reshape(self.x0,self.n0), np.reshape(self.x0.T,self.n0)])
        
        self.L = 1.

    def get_L(self):
        return 1.


    def get_grad(self,x):
        return (x - self.x0)

    def get_obj(self,x):
        return np.sum(np.power(x - self.x0,2.))/2.
    
    def array2im(self,x):
        
        x = x*self.normfact
        
        
        x1 = x[:self.n0]
        x2 = x[self.n0:]
        x1 = np.reshape(x1,self.imshape)
        x2 = np.reshape(x2,(self.imshape[1],self.imshape[0]))
        
        im = Image.fromarray(np.uint8((x1+x2.T)/2.))
        return im
        
    