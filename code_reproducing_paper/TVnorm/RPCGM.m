function [x, gap_track,   xtrack,ztrack,itervec] = RPCGM(n, maxiter, rho, theta, cbar, smoothgrad,  atom_type,p,iter0)

    gamma = @(x) ((log(1+abs(x)./theta).*(abs(x)<=cbar) + (1./(theta+cbar)*(abs(x)-cbar)+log(1+cbar/theta)).*(abs(x)>cbar)));
    gammagrad = @(x) (1./(theta+abs(x)).*(abs(x)<=cbar) + (1./(theta+cbar)).*(abs(x)>cbar));
    gammalin = @(x,xbar)(gammagrad(abs(xbar)).*abs(x));

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
        xbar = x;
        w = gammagrad(xbar);
        r0 = gamma(xbar) - gammalin(xbar,xbar);

        grad = smoothgrad(x);
        v = -grad./ w;
        v(w==0) = 0;

        [sw,nu] = lmo(v,atom_type);
        s = sw ./ w;


        xi = (nu/rho)^(1/(p-1))-sum(r0);

        theta = 2/(2+iter+iter0);
        x = x*(1-theta)+xi*s*theta;


        if iter == 10 || iter == 100 || iter == 1000 || iter == 10000 || iter == 30 || iter == 3 || iter == 300 || iter == 3000 || iter== 30000 || iter == 3e6
            itervec = [itervec,iter];
            i = length(itervec);



            gap = x'*grad+phi(sum(gamma(abs(diff(xbar)))-gammalin(abs(diff(xbar)),diff(xbar))+gammalin(diff(x),diff(xbar))),p)*rho + rho*(phiconj(nu/rho,p))-nu*sum(r0);

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





    