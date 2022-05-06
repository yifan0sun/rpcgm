clc
clear

load data/dorothea_train.mat X y
load data/dorothea_valid.mat Xv yv

%%
[m,n] = size(X);
sigmoid = @(x)(1./(1+exp(-x)));
[I,J] = find(X);
Z = sparse(I,J,y(I));
get_trainobj = @(theta)(-sum(log(sigmoid(Z*theta)))/m);
smoothgrad = @(theta)(Z'*(sigmoid(Z*theta)-1)/m);

%%
maxiter = 1000;
track_obj = zeros(maxiter,1);
track_trainacc = zeros(maxiter,1);
track_testacc = zeros(maxiter,1);
track_trainf1 = zeros(maxiter,1);
track_testf1 = zeros(maxiter,1);
theta = zeros(n,1);
for iter = 1:maxiter
    if mod(iter,100) == 0
        [iter,maxiter]
    end
    theta = theta - smoothgrad(theta)/10;
    track_obj(iter) = get_trainobj(theta);
    
    
    [train, test] = perf_eval(theta);
    track_trainacc(iter) = train.acc;
    track_testacc(iter) = test.acc;   
    track_trainf1(iter) = train.f1;
    track_testf1(iter) = test.f1;     
end

%%
figure(1)
clf
subplot(1,3,1)
plot(track_obj,'linewidth',2)
ylabel('train loss')
xlabel('iterations')

subplot(1,3,2)
plot(track_trainacc,'linewidth',2)
hold on
plot(track_testacc,'linewidth',2)
ylabel('accuracy')
xlabel('iterations')
legend('train','test')

subplot(1,3,3)
plot(track_trainf1,'linewidth',2)
hold on
plot(track_testf1,'linewidth',2)
ylabel('f1')
xlabel('iterations')
legend('train','test')








