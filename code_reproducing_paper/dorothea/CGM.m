function [x,obj_track, train, test, gap_track] = ...
    CGM(n,maxiter, alpha, get_trainobj, smoothgrad,  atom_type)


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


        [s,~] = lmo(-grad,atom_type);

        theta = 2/(2+iter);
        x = x*(1-theta)+alpha*s*theta;
        
        
        [train, test] = update_metrics(train, test, x, iter);

    
        gap = (x-alpha*s)'*grad;
        gap_track(iter) = gap;
        
        obj_track(iter) = get_trainobj(x);
    end
end





