function [x, gap_track,   xtrack,ztrack, itervec] = PCGM(n,maxiter, rho, smoothgrad,  atom_type,p,iter0)


    q = p/(p-1);
    phi = @(s,p)(s^p/p);
    phiconj = @(nu,p)(nu^q/q);
        
    x = zeros(n,1);
    
    
    itervec = [];
    gap_track = zeros(min(maxiter,10000),1);
    xtrack = zeros(min(maxiter,10000),n);
    ztrack = zeros(min(maxiter,10000),n);
    
    
    
    
    for iter = 1:maxiter
        if mod(iter,round(maxiter/10)) == 0
            iter
        end
        grad = smoothgrad(x);
        [s,nu] = lmo(-grad,atom_type);
        xi = (nu/rho)^(1/(p-1));

        theta = 2/(2+iter+iter0);
        x = x*(1-theta)+xi*s*theta;
        
        
        if iter == 10 || iter == 100 || iter == 1000 || iter == 10000 || iter == 30 || iter == 3 || iter == 300 || iter == 3000 || iter== 30000 || iter == 3e6
            itervec = [itervec,iter];
            i = length(itervec);
            
            
            
            gap = x'*grad+phi(sum(abs(diff(x))),p)*rho + rho*phiconj(nu/rho,p);%%%
            gap_track(i) = gap;
            
            
            xtrack(i,:) = x;
            ztrack(i,:) = -grad;
            
            
            
        end
        
    end
    L = length(itervec);
    gap_track = gap_track(1:L,:);
    xtrack = xtrack(1:L,:);
    ztrack = ztrack(1:L,:);
        
end