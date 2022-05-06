function [x,rec_err, sing_rec_err, gap_track,obj_track] = PCGM(U,V,maxiter, rho, smoothgrad, smoothobj, atom_type,p,iter0)
m = size(U,1);
n = size(V,1);
x = zeros(m,n);


[U,~,~] = svd(U);
[V,~,~] = svd(V);

q = p/(p-1);
phi = @(s,p)(s^p/p);
phiconj = @(nu,p)(nu^q/q);

rec_err = zeros(maxiter,1);
gap_track = zeros(maxiter,1);
sing_rec_err = zeros(maxiter,m);
obj_track = zeros(maxiter,1);
for iter = 1:maxiter
    
    grad = smoothgrad(x);
    [s,nu] = lmo(-grad,atom_type);
    xi = (nu/rho)^(1/(p-1));
    
    theta = 2/(2+iter+iter0);
    x = x*(1-theta)+xi*s*theta;
    
    if sum(isnan(x(:))) == 0 && sum(isinf(x(:))) == 0
        [ux,sx,vx] = svd(x);
        
        
        rec_err(iter) = norm(x - U*V','fro')/m/n/2;
        sing_rec_err(iter,:) = diag(sx);
        
        gap = trace(x'*grad)+phi(sum(sum(sx)),p)*rho +rho* phiconj(nu/rho,p);
        gap_track(iter) = gap;
        
        obj_track(iter) = smoothobj(x);
    else
        
        
        rec_err(iter:end) = nan;
        sing_rec_err(iter:end,:) = nan;
        
        gap_track(iter:end) = nan;
        
        obj_track(iter:end) = nan;
        break
    end
    
end

end





