function [x,rec_err,sing_rec_err, gap_track,obj_track] = CGM(U,V,maxiter, alpha, smoothgrad, smoothobj, atom_type)

m = size(U,1);
n = size(V,1);
x = zeros(m,n);

[U,~,~] = svd(U);
[V,~,~] = svd(V);

rec_err = zeros(maxiter,1);
sing_rec_err = zeros(maxiter,m);
gap_track = zeros(maxiter,1);
obj_track = zeros(maxiter,1);


for iter = 1:maxiter
    grad = smoothgrad(x);
    
    
    
    [s,~] = lmo(-grad,atom_type);
    
    theta = 2/(2+iter);
    x = x*(1-theta)+alpha*s*theta;
    
    [ux,sx,vx] = svd(x);
    
    rec_err(iter) = norm(x - U*V','fro')/m/n/2;
    sing_rec_err(iter,:) = diag(sx);
    
    gap = trace((x-alpha*s)'*grad);
    gap_track(iter) = gap;
    
    obj_track(iter) = smoothobj(x);
end
end