import numpy as np


#Phi and gamma functions encoded here. Use Convex, Nonconvex, and LocallyConvex classes to form the object that is fed into the method.

class MonomialAbs():
    def __init__(self, p):
        self.p = p
        self.rmin = 0.
        self.rmax = np.infty
    def get_fval(self,xi):
        if (type(xi) == np.float64) or (type(xi) == float):
            if xi == 0: return xi
            return np.power(np.abs(xi),self.p+0.) / (self.p+0.)
        else: 
            idx = np.not_equal(xi,0)
            xi = xi + 0.
            xi[idx] = np.power(np.abs(xi[idx]),self.p+0.) / (self.p+0.)
            return xi
        
    def get_grad(self,xi):
        if self.p == 2:
            return np.abs(xi) + 0.
        
        
        if (type(xi) == np.float64) or (type(xi) == float):
            if xi == 0: return 0.
            if self.p == 1.: return 1.
            return np.power(np.abs(xi),self.p-1.)
        else: 
            idx = np.not_equal(xi,0)
            xi = xi + 0.
            
            
            if self.p == 1. : xi[idx] = 1.
            else:
                xi[idx] = np.power(np.abs(xi[idx]),self.p-1.)
            return xi
    
    def get_conj_grad(self,nu):
        nu = np.abs(nu)
        if self.p == 2:
            return nu + 0.
    
        if (type(nu) == np.float64) or (type(nu) == float):
            if nu == 0: return 0.
            
            if self.p == 1.: return 1.
            return np.power(nu,1./(self.p-1.))
        else: 
            idx = np.not_equal(nu,0)
            nu = nu + 0.
            if self.p == 1. : nu[idx] = 1.
            else: nu[idx] = np.power(nu[idx],1./(self.p-1.))
            return nu
        
        


class MonomialAbsEps(MonomialAbs):
    def __init__(self, p, eps):
        self.p = p
        self.eps = eps
        self.rmax = self.get_grad(0.)
        self.rmin = 0.
        
    def get_fval(self,xi):
        return MonomialAbs.get_fval(self,xi+self.eps)
        
    def get_grad(self,xi):
        return MonomialAbs.get_grad(self,xi+self.eps)
        
        

class SCAD():
    def __init__(self, c, theta):
        self.c = c+0.
        self.theta = theta+0.
        assert theta > 2., 'theta argument must be > 2'
        assert c > 0., 'c argument must be > 0'
        self.rmax = self.get_grad(0.)
        self.rmin = 0.
    def get_fval(self,xi):
        xi = np.abs(xi)
    
    
        if (type(xi) == np.float64) or (type(xi) == float):
            if xi <= self.c:
                return self.c*xi
            elif xi <= self.theta*self.c:
                return (-(xi**2.) + 2.*self.theta*self.c*xi-self.c**2.)/(2.*(self.theta-1.))
            else:
                return (self.theta+1.)*(self.c**2.)/2.
            
        else: 
            idx1 = np.less_equal(xi,self.c)
            idx2 = np.logical_and(np.greater(xi,self.c), np.less_equal(xi,self.theta*self.c))
            idx3 = np.greater(xi,self.c*self.theta)
            
            g = xi + 0.
            g[idx1] = self.c*xi[idx1]
            g[idx2] = (-(xi[idx2]**2.) + 2.*self.theta*self.c*xi[idx2]-self.c**2.)/(2.*(self.theta-1.))
            g[idx3] = (self.theta+1.)*(self.c**2.)/2.
        
            return g
        
        
        
        
    def get_grad(self,xi):
        xi = np.abs(xi)
    
        if (type(xi) == np.float64) or (type(xi) == float):
            if xi <= self.c:
                return self.c
            elif xi <= self.theta*self.c:
                return (-xi + self.theta*self.c)/(self.theta-1.)
            else:
                return 0.
            
        else: 
            idx1 = np.less_equal(xi,self.c)
            idx2 = np.logical_and(np.greater(xi,self.c), np.less_equal(xi,self.theta*self.c))
            idx3 = np.greater(xi,self.c*self.theta)
            
            g = xi + 0.
            g[idx1] = self.c
            g[idx2] = (-xi[idx2] + self.theta*self.c)/(self.theta-1.)
            g[idx3] = 0.
        
            return g
        
        

class LSP():    
    def __init__(self, theta):
        self.theta = theta+0.
        assert theta > 0., 'theta argument must be > 2'
        self.rmax = self.get_grad(0.)
        self.rmin = 0.
    def get_fval(self,xi):
        xi = np.abs(xi)
    
        return np.log(1.+xi/self.theta)
            
        
    def get_grad(self,xi):    
        return 1./(self.theta+xi)
    



        
#########################################
class Nonconvex():  
    def __init__(self, fn = None):
        self.fn = fn
    def get_fval(self,xi):
        if self.fn is None: return np.abs(xi)
        return self.fn.get_fval(np.abs(xi))
        
    def get_grad(self,xi):
        if self.fn is None: return np.ones(xi.shape)
        g = self.fn.get_grad(np.abs(xi))
        
        return g
####################################3
class LocallyNonconvex(Nonconvex):   


    def __init__(self, fn,  boundary):
        Nonconvex.__init__(self, fn)
        self.boundary = boundary
        self.bval = self.fn.get_fval(boundary)
        self.bgrad = self.fn.get_grad(boundary)
        self.rmin = self.bgrad
        self.rmax = self.fn.get_grad(0)
        
    def get_fval(self,xi):
        xi = np.abs(xi)
    
        if type(xi) == np.float64:
            if xi <= self.boundary: return Nonconvex.get_fval(self,np.abs(xi))
            else: return (xi - self.boundary)*self.bgrad + self.bval
            
        idx = np.less_equal(np.abs(xi),self.boundary)
        nidx = np.logical_not(idx)
        f = xi * 0.
        f[idx] = Nonconvex.get_fval(self,np.abs(xi[idx]))
        f[nidx] = (xi[nidx] - self.boundary)*self.bgrad + self.bval
        
        return f
        
        
    def get_grad(self,xi):
        xi = np.abs(xi)
        if type(xi) == np.float64:
            if xi <= self.boundary: return Nonconvex.get_grad(self,np.abs(xi))
            else: return self.bgrad
            
        idx = np.less_equal(np.abs(xi),self.boundary)
        g = xi * 0.
        g[idx] = Nonconvex.get_grad(self,np.abs(xi[idx]))
        g[np.logical_not(idx)] = self.bgrad
        
        return g
        
        
class Convex():  
    def __init__(self, fn):
        self.fn = fn
    def get_fval(self,xi):
        return self.fn.get_fval(np.abs(xi))
        
    def get_conj_grad(self,xi):
        return self.fn.get_conj_grad(np.abs(xi))
        