function [x,obj_track, train, test, gap_track] = PCGM(n,maxiter, rho, get_trainobj, smoothgrad,  atom_type,p,iter0)



    q = p/(p-1);
    phi = @(s,p)(s^p/p);
    phiconj = @(nu,p)(nu^q/q);
        
    x = zeros(n,1);
    
    obj_track = zeros(maxiter,1);
    gap_track = zeros(maxiter,1);
    
    train = get_metrics_obj(maxiter);
    test = get_metrics_obj(maxiter);
    
    for iter = 1:maxiter
        if mod(iter,100) == 0
            iter
        end
        
        grad = smoothgrad(x);
        [s,nu] = lmo(-grad,atom_type);
        xi = (nu/rho)^(1/(p-1));

        theta = 2/(2+iter+iter0);
        x = x*(1-theta)+xi*s*theta;
        
        [train, test] = update_metrics(train, test, x, iter);

        gap = x'*grad+phi(sum(abs(x)),p)*rho + rho*phiconj(nu/rho,p);
        gap_track(iter) = gap;
        obj_track(iter) = get_trainobj(x);
    end
        
end


