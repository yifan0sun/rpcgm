function [x,rec_err, sing_rec_err, gap_track, obj_track] = RPCGM(U,V,maxiter, rho, theta, cbar, smoothgrad, smoothobj, atom_type,p,iter0)




gamma = @(x) ((log(1+abs(x)./theta).*(abs(x)<=cbar) + (1./(theta+cbar)*(abs(x)-cbar)+log(1+cbar/theta)).*(abs(x)>cbar)));
gammagrad = @(x) (1./(theta+abs(x)).*(abs(x)<=cbar) + (1./(theta+cbar)).*(abs(x)>cbar));
gammalin = @(x,xbar)(gammagrad(abs(xbar)).*abs(x));

q = p/(p-1);
phi = @(s)(s^p/p);
phiconj = @(nu)(nu^q/q);



m = size(U,1);
n = size(V,1);
x = zeros(m,n);


[U,~,~] = svd(U);
[V,~,~] = svd(V);

rec_err = zeros(maxiter,1);
gap_track = zeros(maxiter,1);
sing_rec_err = zeros(maxiter,m);
obj_track = zeros(maxiter,1);
for iter = 1:maxiter
    xbar = x;
    
    grad = smoothgrad(x);
    [uz,sz,vz] = svd(-grad);
    sbar = diag(uz'*xbar*vz);
    [ux,sx,vx] = svd(x);
    
    w = gammagrad(sbar);
    r0 = sum(gamma(sbar) - gammalin(sbar,sbar));
    v = diag(sz)./w;
    v(w==0) = 0;
    [sw,nu] = lmo(v,atom_type);
    s = sw ./ w;
    xi = (nu/rho)^(1/(p-1))-r0;
    theta = 2/(2+iter+iter0);
    
    x = x*(1-theta)+xi*(uz*diag(s)*vz')*theta;
    
    if sum(isnan(x(:))) == 0 && sum(isinf(x(:))) == 0
        rec_err(iter) = norm(x - U*V','fro')/m/n/2;
        
        
        ss = diag(uz'*x*vz);
      
        sing_rec_err(iter,:) = svd(x);


        gap = trace(x'*grad)+...
            phi(r0+sum(gammalin(ss,sbar)))*rho + rho*phiconj(nu/rho)-nu*r0;
        
% [gap,trace(x'*grad),diag(sz)'*(ss), phi(sum(r0+gammalin(ss,sbar)),p)*rho ,phi(sum(gamma(ss)),p)*rho ,  rho*phiconj(nu/rho,p),-nu*sum(r0)]
        gap_track(iter) = gap;
        obj_track(iter) = smoothobj(x);
        gap
    else
        
        
        rec_err(iter:end) = nan;
        sing_rec_err(iter:end,:) = nan;
        
        gap_track(iter:end) = nan;
        
        obj_track(iter:end) = nan;
        break
    end
    
    
end

end

















