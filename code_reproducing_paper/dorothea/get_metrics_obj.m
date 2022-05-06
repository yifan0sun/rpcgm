function train = get_metrics_obj(maxiter)

    train.P = zeros(maxiter,1);
    train.R = zeros(maxiter,1);
    train.f1 = zeros(maxiter,1);
    train.acc = zeros(maxiter,1);
end