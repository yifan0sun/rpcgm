clc
clear

%%
% splevel = 0.1;
% noise_variance = 0;
% for trial = 1:100
%     for m = [10,25,50,100,250,500,1000,2500,5000,10000]
%
%         n = 100;
%         A = randn(m,n);
%
%         x0 = zeros(n,1);
%         idx0 = rand(n,1) < splevel;
%         x0(idx0) = 1;
%         b = A*x0+noise_variance*randn(m,1);
%
%         savefilename = sprintf('small_experiments/problem_m%d_n%d_trial%d_splevel%f_noise%f.mat',m,n,trial,splevel,noise_variance);
%         save(savefilename,'A','b','x0','idx0')
%     end
% end

%%

maxiter = 1000000;
n = 100;
splevel = 0.1;
noise_variance = 0;
iter0 = 0;
m = 100;
trial = 1;

trial
m
loadfilename = sprintf('saved_data/problem_m%d_n%d_trial%d_splevel%f_noise%f.mat',m,n,trial,splevel,noise_variance);
data = load(loadfilename);
A = data.A;
b = data.b;
idx0 = data.idx0;
x0 = data.x0;

smoothobj = @(x)( sum((A*x-b).^2)/m/2 );
smoothgrad = @(x)(A'*(A*x-b)/m);

L = norm(A'*A,2);
%%
for alpha = [1,2,5,10]
    alpha
    filename = sprintf('saved_data/CGM_m%d_n%d_trial%d_splevel%f_noise%f_alpha%f_iter0%d.mat',m,n,trial,splevel,noise_variance,alpha,iter0);

        [x,supp_err, gap_track, supp, screen, xtrack,ztrack,itervec] = CGM(n,idx0,maxiter, alpha, smoothgrad, '1norm',iter0,L);
        save(filename,'x','gap_track','supp','supp_err','screen', 'xtrack','ztrack','itervec')


end
    
    
    
    
%%

for p = [1.25,2,5]
    if p == 1.25
        rhovec = [.1,.5,1,2,5,10];
        iter0 = 1000;
    elseif p == 2
        rhovec = [.001,.002,.005,.01,.02,.05,.1,.2,.5,1];
        iter0 = 1000;
    elseif p == 5
        rhovec = [.00001,.00002,.00005,.0001,.0002,.0005,.001,.002,.005,.01,.1];
        iter0 = 1000;
    end
    for rho = rhovec
        [p,rho]
        filename = sprintf('saved_data/PCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,iter0);

            [x,supp_err, gap_track, supp, screen, xtrack,ztrack, itervec] = PCGM(n,idx0,maxiter, rho, smoothgrad,  '1norm',p,iter0,L);
            save(filename,'x','gap_track','supp','supp_err','screen', 'xtrack','ztrack','itervec')
    end
end
%%
maxiter = 100
for iter0 = [1000]
    for p = [2]
        for theta = [1/8,1/4,1/2,3/4,7/8]
            for rho = 10%[.0001,.001,.01,.1,1,10,100]
                [theta,rho,p]
                cbar = 1-theta;
                filename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,cbar,theta,iter0);

                    [x,supp_err, gap_track, supp, screen, xtrack,ztrack,itervec] = RPCGM(n,idx0,maxiter, rho, theta, cbar, smoothgrad,  '1norm',p,iter0,L,0);

                figure(1)
                clf
                plot(gap_track)
            end
        end
    end
end











