function [x,supp_err, gap_track, supp, screen, xtrack,ztrack, itervec] = CGM(n,idx0,maxiter, alpha, smoothgrad,  atom_type,iter0,L)


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
        
        grad = smoothgrad(x);


        [s,~] = lmo(-grad,atom_type);

        theta = 2/(2+iter+iter0);
        x = x*(1-theta)+alpha*s*theta;
        idx = x ~= 0;
        
        
        
        if (iter < 100) || ((iter < 1e4) && (mod(iter , 10) == 0))|| ((iter < 1e6) && (mod(iter , 100) == 0))
            itervec = [itervec,iter];
            i = length(itervec);
            supp_err(i) =sum(abs(idx-idx0))/n;
            gap = (x-alpha*s)'*grad;
            gap_track(i) = gap;

            xtrack(i,:) = x;
            ztrack(i,:) = -grad;
            supp(i,:) = x ~= 0;
            screen(i,:) = (max(abs(grad))) - abs(grad) <  (2*sqrt(L*gap));
        end
        
    end
end