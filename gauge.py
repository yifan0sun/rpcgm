import numpy as np
import scipy.sparse as ss

#define custom gauges here

class Gauge():
    def screen(self,zmax, z, resbnd):
        res = zmax - z
        return np.greater(res, resbnd)
    
        
class OneNorm(Gauge):
    def __init__(self,n):
        self.n = n
        self.Lfactor = 1.
        
    def LMO(self,z, weights = None):
        
        if weights is not None:
            z /= weights
        k = np.argmax(np.abs(z))
        sgn = np.sign(z[k])
        s = np.zeros(self.n)
        s[k] = sgn
        return s, np.abs(z[k])
    
    def get_fval(self,x, nonconvex = None, weights = None):
        x = np.abs(x)
        if weights is not None:
            x = x * weights
        if nonconvex is not None:
            x = nonconvex.get_val(x)
            
        return np.sum(x)
    
    def get_dval(self,z, weights = None):
        if weights is not None:
            z = z/ weights
        return np.max(np.abs(z))
    
    
    def screen(self,z, resbnd):
        zmax = np.max(np.abs(z))
        z = np.abs(z)
        return Gauge.screen(self,zmax,z,resbnd)

    def get_atom_weights(self,x):
        return np.abs(x)
    
    def miny(self,f,x,y):
        return y
    
class TVNorm(Gauge):
    def __init__(self,n):
        self.n = n
        #V = np.hstack([np.ones(n-1), -np.ones(n-1)])
        #I = np.hstack([range(n-1),range(n-1)])
        #J = np.hstack([range(n-1),range(1,n)])
        #self.P = ss.coo_matrix((V, (I,J)), shape = (n-1,n))
        #self.P = self.P.todense()
        
        self.Lfactor = n+0.
        
    def adjointDiff(self,z):
        
        u = np.zeros(self.n-1)
        u[0] = z[0]
        for k in xrange(1,self.n-1):
            u[k]  = z[k] +  u[k-1]
        u = np.array(u)
        return u
    
    def LMO(self,z, weights = None):
        z = z - np.mean(z)
        u = self.adjointDiff(z)
        
        
        if weights is not None:
            u /= weights
        
        
        k = np.argmax(np.abs(u))
        sgn = np.sign(u[k])
        s = np.zeros(self.n)
        s[:k] = sgn
        #s[(k+1):] = -sgn
        s = s - np.mean(s)
        
        
        #print np.linalg.norm(np.dot(self.P.T,u)-z), np.abs(np.dot(z,s)), np.max(np.abs(u)), u[k]
        return s, np.max(np.abs(u))#np.abs(np.dot(z,s))
    
    def get_fval(self,x, weights = None, nonconvex = None):
        x = np.abs(x[1:] - x[:-1])
        if weights is not None:
            x = x * w
        if nonconvex is not None:
            x = nonconvex.get_val(x)
            
        return np.sum(x)
    
    def get_dval(self,z, weights = None):
        u = self.adjointDiff(z)
        if weights is not None:
            u /= weights
        return np.max(np.abs(u))
    
    
    
    def screen(self,z, resbnd):
        u = self.adjointDiff(z)
        umax = np.max(np.abs(u))
        u = np.abs(u)
        return Gauge.screen(self,umax,u,resbnd)
    

    def get_atom_weights(self,x):   
        return np.abs(x[1:]-x[:-1])
    def miny(self,f,x,y):
        stepsize = 1./f.L
        for iter in xrange(1000):
            g = np.mean(f.get_grad(x+y))
            if np.abs(g) < .000001*np.linalg.norm(x)/(self.n+0.): 
                break
            y -= stepsize*g*np.ones(self.n)
           
        return y
    
class GroupNorm(Gauge):
    def __init__(self,n, groups):
        self.n = n
        self.groups = groups
        self.Lfactor = 1.
        
        self.form_proxy_matrices()
        self.update_gauge = False
        
    def form_proxy_matrices(self):
        I,J,V = [],[],[]
        offset = 0
        for idx in self.groups:
            I.extend(idx)
            J.extend(range(offset, offset + len(idx)))
            V.extend([1 for i in idx])
            offset += len(idx)


        P = ss.coo_matrix((V,(I,J)), shape = (self.n, sum([len(idx) for idx in self.groups])))
        P2 = np.dot(P,P.T)

        V = []
        offset = 0
        for idx in self.groups:
            V.extend([1./P2[i,i] for i in idx])
            offset += len(idx)

        Proj = ss.coo_matrix((V,(I,J)), shape = (self.n, sum([len(idx) for idx in self.groups])))

        self.P = P
        self.Proj = Proj.T
    
    
        
    def LMO(self,z, weights = None):
        zi = [np.linalg.norm(z[idx]) for idx in self.groups]
        if weights is not None:
            zi /= weights
        k = np.argmax(zi)
        
        s = np.zeros(self.n)
        s[self.groups[k]] = z[self.groups[k]]
        s[self.groups[k]]  = s[self.groups[k]]  / np.linalg.norm(s[self.groups[k]] )
        return s, zi[k]
    
    def get_dval(self,z, weights = None):
        zi = [np.linalg.norm(z[idx]) for idx in self.groups]
        if weights is not None:
            zi /= weights
        return np.max(zi)
    
    def solve_gauge_problem(self,x,xstart = None):
        atoms, val, res, affres = get_latent_groupnorm(x,self.groups,xstart=None, P = self.P, Proj = self.Proj)
        self.atoms = atoms
        self.gauge_val = val
        self.update_gauge = True
    
    def get_fval(self,x):
        if not self.update_gauge:
            if hasattr(self,'atoms'):
                xstart = self.atoms
            else:
                xstart = None
            self.solve_gauge_problem(x, xstart)
        return self.gauge_val
    
    
    def get_atom_weights(self,x):
        if np.linalg.norm(x) == 0.:
            return np.zeros(len(self.groups))
        
        
        if not self.update_gauge:
            if hasattr(self,'atoms'):
                xstart = self.atoms
            else:
                xstart = None
            self.solve_gauge_problem(x, xstart)
        return self.atoms
    
    def screen(self,z, resbnd):
        zi = [np.linalg.norm(z[idx]) for idx in self.groups]
        zgmax = np.max(zi)
        return Gauge.screen(self,zgmax,zi,resbnd)
    
    def miny(self,f,x,y):
        return y
    
    
    

def get_latent_groupnorm(b,groups,xstart, P, Proj):
    
    def shrinkage_scalar(u,t):
        if u < -t: return u + t
        if u > t: return u - t
        return 0.
    
    n = len(b)
    n_latent = sum(len(idx) for idx in groups)
    if xstart is None:
        z = np.zeros(n_latent)
    else:
        z = xstart
    t = 1.
    rho = 1.
    x = z + 0.
    for iter in xrange(100):
        
        #x = proxf(z,t)
        offset = 0
        for k in xrange(len(groups)):
            zn = np.linalg.norm(z[offset:offset + len(groups[k])])
            zn_new = shrinkage_scalar(zn,t)
            if zn_new == 0.: 
                zn = 1.
            x[offset:offset + len(groups[k])] = z[offset:offset + len(groups[k])] / zn * zn_new
            offset += len(groups[k])
        
       
        
        #y = proj(2*x-z)
        y = 2.*x-z
        y = y - ss.coo_matrix.dot(Proj ,ss.coo_matrix.dot(P,y)-b)
        
        
        z = z + rho*(y-x) 
    
        res = np.linalg.norm(y-x) 
        affres = np.linalg.norm(ss.coo_matrix.dot(P,x)-b)
        obj = sum(np.linalg.norm(x[idx]) for idx in groups)
        #print iter, obj,res ,affres
        if res < 1.e-8:
            break
    return x, obj, res,affres
