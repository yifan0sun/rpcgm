function [x,supp_err, gap_track, supp, screen, xtrack,ztrack,itervec] = RPCGM(n,idx0,maxiter, rho, theta, cbar, smoothgrad,  atom_type,p,iter0,L,screen_eps)

    gamma = @(x) ((log(1+abs(x)./theta).*(abs(x)<=cbar) + (1./(theta+cbar)*(abs(x)-cbar)+log(1+cbar/theta)).*(abs(x)>cbar)));
    gammagrad = @(x) (1./(theta+abs(x)).*(abs(x)<=cbar) + (1./(theta+cbar)).*(abs(x)>cbar));
    gammalin = @(x,xbar)(gammagrad(abs(xbar)).*abs(x));

    q = p/(p-1);
    phi = @(s,p)(s^p/p);
    phiconj = @(nu,p)(nu^q/q);

    x = zeros(n,1);

    itervec = [];
    supp_err = zeros(min(maxiter,10000),1);
    gap_track = zeros(min(maxiter,10000),1);
    supp = zeros(min(maxiter,10000),n);
    screen = zeros(min(maxiter,10000),n);
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
        idx = x ~= 0;


        if (iter < 100) || ((iter < 1e4) && (mod(iter , 10) == 0))|| ((iter < 1e6) && (mod(iter , 100) == 0))
            itervec = [itervec,iter];
            i = length(itervec);

            supp_err(i) =sum(abs(idx-idx0))/n;


            gap = x'*grad+phi(sum(gamma(xbar)-gammalin(xbar,xbar)+gammalin(x,xbar)),p)*rho + rho*(phiconj(nu/rho,p))-nu*sum(r0);

            gap_track(i) = gap;

            xtrack(i,:) = x;
            ztrack(i,:) = -grad;
            supp(i,:) = x ~= 0;
            screen(i,:) =( (max(abs(grad))) - abs(grad) )< (screen_eps + (2*sqrt(L*gap + screen_eps)));

        end


    end

end


