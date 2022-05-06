function [x,supp_err, gap_track, supp, screen, xtrack,ztrack, itervec] = CGM(n,groups, active_groups,maxiter, alpha, smoothgrad,  atom_type,iter0,L)

    
    x = zeros(n,1);
    idx0 = x*0;
    for k = active_groups
        idx0(groups{k}) = 1;
    end
    itervec = [];
    supp_err = zeros(min(maxiter,10000),1);
    gap_track = zeros(min(maxiter,10000),1);
    supp = zeros(min(maxiter,10000),n);
    screen = zeros(min(maxiter,10000),length(groups));
    xtrack = zeros(min(maxiter,10000),n);
    ztrack = zeros(min(maxiter,10000),n);
    
    nz = zeros(length(groups),1);
    for iter = 1:maxiter
        if mod(iter,round(maxiter/10)) == 0
            iter
        end
        
        grad = smoothgrad(x);


        [s,~] = lmo(-grad,atom_type, groups);

        theta = 2/(2+iter+iter0);
        x = x*(1-theta)+alpha*s*theta;
        idx = x ~= 0;
        
        
        
        if (iter < 100) || ((iter < 1e4) && (mod(iter , 10) == 0))|| ((iter < 1e6) && (mod(iter , 100) == 0)) || iter == 30 || iter == 3 || iter == 300 || iter == 3000 || iter== 30000 || iter == 3e6
            itervec = [itervec,iter];
            i = length(itervec);
            supp_err(i) =sum(abs(idx-idx0))/n;
            gap = (x-alpha*s)'*grad;
            gap_track(i) = gap;

            xtrack(i,:) = x;
            supp(i,:) = x ~= 0;
            
            for k = 1:length(groups)
                nz(k) = norm(grad(groups{k}));
            end
            ztrack(i,:) = -grad;
            
            screen(i,:) = (max(nz) - nz) <  (2*sqrt(L*gap));
        end
        
    end
    
    
    L = length(itervec);
    supp_err = supp_err(1:L,:);
    gap_track = gap_track(1:L,:);
    screen = screen(1:L,:);
    supp = supp(1:L,:);
    xtrack = xtrack(1:L,:);
    ztrack = ztrack(1:L,:);
end