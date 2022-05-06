function [train, test] = perf_eval(theta)


load data/dorothea_train.mat X y
load data/dorothea_valid.mat Xv yv

yhat = sign(X*theta);
yvhat = sign(Xv*theta);

pos = sum((yhat > 0).*(y > 0));
train.P = pos / sum(y>0);
train.R = pos / sum(yhat > 0);
train.f1 = train.P*train.R/(train.P+train.R);
train.acc = mean((yhat > 0)==(y > 0));
if pos == 0
    train.P = 0;
    train.R = 0;
    train.f1 = 0;
end

pos = sum((yvhat > 0).*(yv > 0));
test.P = pos / sum(yv>0);
test.R = pos / sum(yvhat > 0);
test.f1 = test.P*test.R/(test.P+test.R);
if pos == 0
    test.P = 0;
    test.R = 0;
    test.f1 = 0;
end
test.acc = mean((yvhat > 0)==(yv > 0));
end