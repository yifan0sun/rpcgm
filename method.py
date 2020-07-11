import numpy as np
import matplotlib.pyplot as plt

#P-CGM and RP-CGM coded here

class Method():
    def __init__(self,prob,opts, tracklist, xstart = None):
        if xstart is None:
            self.xstart = np.zeros(prob.n)
        else:
            self.xstart = xstart
            
        
        self.prob = prob
        self.maxiter = opts['maxiter']
        self.track_period = opts['track_period']
        if 'debug_fn' in opts.keys():
            self.debug_fn = opts['debug_fn']
        else: self.debug_fn = None
        
        self.track = {'iter':[]}
        for t in tracklist:
            self.track[t] = []
            
    def run(self):
        self.vars['x'] = self.xstart + 0.
        self.vars['y'] = np.zeros(self.xstart.shape)

        for iter in xrange(self.maxiter):
            self.LMO_step()
            self.screen(iter)
            self.merge_step(iter)
            
            if iter % self.track_period == 0:
                self.update_tracking(iter)
            if self.debug_fn is not None:
                self.debug_fn(self, iter)
        
        
        for k,v in self.track.items():
            self.track[k] = np.array(v)
        

    def update_tracking(self,iter):
        self.track['iter'].append(iter)
        if 'obj' in self.track.keys():
            self.set_full_obj()
            self.track['obj'].append(self.obj)
        if 'nnz' in self.track.keys():
            self.track['nnz'].append(np.sum(np.greater(np.abs(self.vars['x']),1.e-8)))
            

        
class PCGM(Method):
    def __init__(self,prob,gauge,convexphi, lam,opts, tracklist, xstart = None):
        Method.__init__(self,prob,opts, tracklist, xstart)
            
        self.gauge = gauge
        self.convexphi = convexphi
        self.lam = lam
        self.L = prob.get_L() * gauge.Lfactor
        self.vars = {}
        
        self.screened = self.gauge.get_atom_weights(self.xstart*0.)
        self.screened[:] = False
        
    
    def update_tracking(self,iter):
        Method.update_tracking(self,iter)

        if 'gap' in self.track.keys():
            self.track['gap'].append(self.gap)
        if 'resbnd' in self.track.keys():
            self.track['resbnd'].append(self.resbnd)
            
        if 'nscreened' in self.track.keys():
            self.track['nscreened'].append(np.sum(self.screened))
            
        if 'pr_screen' in self.track.keys():
            
            guess = np.not_equal(self.screened,1.)
            truth = np.not_equal(self.prob.x_true,0.)
            overlap = np.logical_and(guess,truth)
            precision = sum(overlap)/(sum(guess)+0.)
            recall = sum(overlap)/(sum(truth)+0.)
            
            self.track['pr_screen'].append(np.vstack([precision,recall]))
            
            
        
    
            
        
    def h(self,x):
        return self.convexphi.get_fval(self.gauge.get_fval(x))
    
    
    def set_full_obj(self):
        self.obj = self.prob.get_obj(self.vars['x'] + self.vars['y']) + self.lam*self.h(self.vars['x'])
        

        
    def LMO_step(self):
        x = self.vars['x']
        y = self.vars['y']
        g = self.prob.get_grad(x+y)
        s,nu = self.gauge.LMO(-g)
        xi = self.convexphi.get_conj_grad(nu/self.lam)
        s *=  xi
        
        
        self.vars['s'] = s
        self.vars['g'] = g
        self.vars['nu'] = nu
        self.vars['xi'] = xi
        

    def merge_step(self,iter):
        t = 2./(1.+iter)
        self.vars['x']  = t * self.vars['s'] + (1.-t) * self.vars['x']
        x,y = self.vars['x'],self.vars['y']
        self.vars['y'] = self.gauge.miny(self.prob,x,y)
        
        if hasattr(self.gauge, 'update_gauge'):
            self.gauge.update_gauge = False
        
   
    def screen(self, iter):
        self.set_gap()
        self.set_resbnd()        
        #if self.gap <= 0.: 
        #    print 'warning, gap value %d < 0' % self.gap
        screened = self.gauge.screen(-self.vars['g'],self.resbnd)
        self.screened[screened] = True 

    
    
    def set_gap(self):
        g,s = self.vars['g'],self.vars['s']
        x = self.vars['x']
        self.gap = np.dot(-g,(s-x)) +self.lam*self.h(x) - self.lam*self.h(s)
        
        
    def set_resbnd(self):
        self.resbnd = 2.*np.sqrt(self.L*self.gap)
         
            
        
class RPCGM(PCGM, Method):
    def __init__(self,prob,gauge,convexphi,nonconvex, lam,opts, tracklist, xstart = None):
        PCGM.__init__(self,prob,gauge,convexphi, lam,opts, tracklist, xstart)
        self.nonconvex = nonconvex
        
        self.screened = self.gauge.get_atom_weights(self.xstart*0.)
        self.screened[:] = False
       

    def update_tracking(self,iter):
        PCGM.update_tracking(self,iter)
        
        if 'res' in self.track.keys():
            self.track['res'].append(self.res)
        if 'nscreened' in self.track.keys():
            self.track['nscreened'].append(np.sum(self.screened))
            
        if 'pr_screen' in self.track.keys():
            guess = np.not_equal(self.screened,1.)
            truth = np.not_equal(self.prob.x_true,0.)
            overlap = np.logical_and(guess,truth)
            precision = sum(overlap)/(sum(guess)+0.)
            recall = sum(overlap)/(sum(truth)+0.)
            self.track['pr_screen'].append(np.vstack([precision,recall]))
    

            
    def r(self,x):
        return np.sum(self.nonconvex.get_fval(self.gauge.get_atom_weights(x)))
    def h(self,x):
        return self.convexphi.get_fval(self.r(x))
    
    
    
    def rlin(self,x, xbar):
        atom_weights_xbar = self.gauge.get_atom_weights(xbar)
        nonconvex_slope_xbar = self.nonconvex.get_grad(atom_weights_xbar)
        atom_weights_x = self.gauge.get_atom_weights(x)
        return np.dot(nonconvex_slope_xbar, atom_weights_x)
    
    def hlin(self,x, xbar):
        return self.convexphi.get_fval(self.r(xbar) - self.rlin(xbar,xbar) + self.rlin(x,xbar))
        
        
    def LMO_step(self):
        x = self.vars['x']
        y = self.vars['y']
        g = self.prob.get_grad(x+y)
        
        a = self.gauge.get_atom_weights(x)
        w = self.nonconvex.get_grad(a)
        
        s,nu = self.gauge.LMO(-g,weights=w)
        xi = self.convexphi.get_conj_grad(nu/self.lam)
        s *=  xi
        
        
        self.vars['s'] = s
        self.vars['g'] = g
        self.vars['nu'] = nu
        self.vars['xi'] = xi
        

   
    def screen(self,iter):
        self.set_res()
        self.set_resbnd()        
        self.set_resbnd()        
        #if self.gap <= 0.: 
        #    print 'warning, gap value %d < 0' % self.gap
        #if self.res <= 0.: 
        #    print 'warning, res value %d < 0' % self.res
        
        screened = self.gauge.screen(-self.vars['g'],self.resbnd)
        self.screened[screened] = True 

    
        if iter > 5:
            screened = self.gauge.screen(-self.vars['g'],self.resbnd)
            self.screened[screened] = True 
    
    def set_gap(self):
        g,s = self.vars['g'],self.vars['s']
        x = self.vars['x']
        self.gap = np.dot(-g,(s-x)) +self.lam*self.hlin(x,x) - self.lam*self.hlin(s,x)
    
    def set_res(self):
        self.set_gap()
        x = self.vars['x']
        g = self.vars['g']
        
        a = self.gauge.get_atom_weights(x)
        w = self.nonconvex.get_grad(a)
        
        nu = self.gauge.get_dval(-g, weights = w)
        
        self.res = self.gap + (self.rlin(x,x) - self.r(x))*nu
        
       
                
        
    def set_resbnd(self):
        self.resbnd = 2.*np.sqrt(self.L*self.res)*self.nonconvex.rmax
    