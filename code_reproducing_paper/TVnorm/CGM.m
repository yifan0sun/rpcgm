function [x, gap_track,  xtrack,ztrack, itervec] = CGM(n, maxiter, alpha, smoothgrad,  atom_type,iter0)

    
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


        [s,~] = lmo(-grad,atom_type);

        theta = 2/(2+iter+iter0);
        x = x*(1-theta)+alpha*s*theta;
        idx = x ~= 0;
        
        
        
        if (iter < 100) || ((iter < 1e4) && (mod(iter , 10) == 0))|| ((iter < 1e6) && (mod(iter , 100) == 0)) || iter == 30 || iter == 3 || iter == 300 || iter == 3000 || iter== 30000 || iter == 3e6
            itervec = [itervec,iter];
            i = length(itervec);
            
            
            gap = (x-alpha*s)'*grad;
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