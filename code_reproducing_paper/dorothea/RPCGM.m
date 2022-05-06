function [x,obj_track, train, test, gap_track] = RPCGM(n,maxiter, rho, theta, cbar, get_trainobj, smoothgrad,  atom_type,p,iter0)
    gamma = @(x) (sum(log(1+abs(x)./theta).*(abs(x)<=cbar) + (1./(theta+cbar)*(abs(x)-cbar)+log(1+cbar/theta)).*(abs(x)>cbar)));
    gammagrad = @(x) (1./(theta+abs(x)).*(abs(x)<=cbar) + (1./(theta+cbar)).*(abs(x)>cbar));
    gammalin = @(x,xbar)(gammagrad(abs(xbar))'*abs(x));

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
        xbar = x;
        w = gammagrad(xbar);
        r0 = gamma(xbar) - gammalin(xbar,xbar);

        grad = smoothgrad(x);
        v = -grad./ w;
        v(w==0) = 0;
        
        [sw,nu] = lmo(v,atom_type);
        s = sw ./ w;
        
        
        xi = (nu/rho)^(1/(p-1))-r0;

        theta = 2/(2+iter+iter0);
        x = x*(1-theta)+xi*s*theta;
        [train, test] = update_metrics(train, test, x, iter);

        
        gap = x'*grad+phi(sum(gamma(xbar)-gammalin(xbar,xbar)+gammalin(x,xbar)),p)*rho + rho*(phiconj(nu/rho,p)-nu*r0/rho);
        gap_track(iter) = gap;
        obj_track(iter) = get_trainobj(x);
    end

end








