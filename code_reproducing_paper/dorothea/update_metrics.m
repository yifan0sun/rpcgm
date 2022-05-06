function [train, test] = update_metrics(train, test, theta, iter)


        [tr, te] = perf_eval(theta);

        train.P(iter) = tr.P;
        train.R(iter) = tr.R;
        train.f1(iter) = tr.f1;
        train.acc(iter) = tr.acc;
        
        
        test.P(iter) = te.P;
        test.R(iter) = te.R;
        test.f1(iter) = te.f1;
        test.acc(iter) = te.acc;
        
end