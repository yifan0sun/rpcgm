clc
clear

%%
noise_variance = 0;
trial = 1;
n = 100;
m = 1000;
num_edges = 4;

% offset = randi(n,1,num_edges)
% x0 = zeros(n,1);
% 
% for k = 1:num_edges
% x0(offset(k)+1:end) = x0(offset(k)+1:end) + sign(randn());
% end
% figure(1)
% clf
% plot(x0)
% 
% 
% A = randn(m,n);
% 
% b = A*x0+noise_variance*randn(m,1);
% 
% savefilename = sprintf('saved_data/problem_m%d_n%d_trial%d_num_edges%d_noise%f.mat',m,n,trial,num_edges,noise_variance)
% save(savefilename,'A','b','x0','offset')

%%

iter0 = 0;

trial
m
loadfilename = sprintf('saved_data/problem_m%d_n%d_trial%d_num_edges%d_noise%f.mat',m,n,trial,num_edges,noise_variance);
data = load(loadfilename);
A = data.A;
b = data.b;
x0 = data.x0;



smoothobj = @(x)( sum((A*x-b).^2)/m/2 );
smoothgrad = @(x)(A'*(A*x-b)/m);

L = norm(A'*A,2);


%%

maxiter = 10001;
figure(1)
clf
for  alpha = [1]
    alpha
    filename = sprintf('saved_data/CGM_m%d_n%d_trial%d_num_edges%d_noise%f_alpha%f_iter0%d.mat',m,n,trial, num_edges,noise_variance,alpha,iter0);
    try
        
        data = load(filename);
        gap_track = data.gap_track;
        x = data.x;
        itervec = data.itervec;
    catch
        [x, gap_track,  xtrack,ztrack, itervec] = CGM(n, maxiter, alpha, smoothgrad,  'tv',iter0);
%         [x,supp_err, gap_track, supp, screen, xtrack,ztrack,itervec] = CGM(n,groups, active_groups,maxiter, alpha, smoothgrad, 'groupnorm',iter0,L);
        save(filename,'x','gap_track','xtrack','ztrack','itervec')
    end
    subplot(2,1,1)
    loglog(gap_track)
    hold on
    subplot(2,1,2)
    plot(x0,'k')
    hold on
    stem(x)
end



%%

maxiter = 10001;
figure(2)
clf
for p = [2]%,5]%,2,5]
    if p == 1.25
        rhovec = 5;%[.1,.5,1,2,5,10];
%         iter0 = 1000;
    elseif p == 2
        rhovec = .1;%[.001,.002,.005,.01,.02,.05,.1,.2,.5,1];
%         iter0 = 1000;
    elseif p == 5
        rhovec = .01;%[.00001,.00002,.00005,.0001,.0002,.0005,.001,.002,.005,.01,.1];
%         iter0 = 1000;
    end
    rhovec = [.1];
    for rho = rhovec
        [p,rho]
        filename = sprintf('saved_data/PCGM_m%d_n%d_trial%d_num_edges%d_noise%f_rho%f_p%f_iter0%d.mat',m,n,trial,num_edges,noise_variance,rho,p,iter0);
        
            

        try
            
            data = load(filename);
            gap_track = data.gap_track;
            x = data.x;
            itervec = data.itervec;
        catch


            [x, gap_track,   xtrack,ztrack, itervec] = PCGM(n,maxiter, rho, smoothgrad,  'tv',p,iter0);
            save(filename,'x','gap_track', 'xtrack','ztrack','itervec')
        end
    subplot(2,1,1)
    loglog(itervec,gap_track)
    hold on
    subplot(2,1,2)
    plot(x0,'k')
    hold on
    stem(x)
    end
end



%%

maxiter = 10001;
figure(2)
clf
% maxiter = 100
for iter0 = [100]
    for p = [2]
        for theta = 3/4%[1/8,1/4,1/2,3/4,7/8]
            for rho = .25
                [theta,rho,p]
                cbar = 1-theta;
                filename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_num_edges%d_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,num_edges, noise_variance,rho,p,cbar,theta,iter0);

                try
                    
                    data = load(filename);
                    gap_track = data.gap_track;
                    x = data.x;
                    supp_err = data.supp_err;
                    itervec = data.itervec;
                    supp = data.supp;
                    screen = data.screen;
                catch

                    [x, gap_track,   xtrack,ztrack,itervec] = RPCGM(n, maxiter, rho, theta, cbar, smoothgrad,  'tv',p,iter0);
        save(filename,'x','gap_track', 'xtrack','ztrack','itervec')
                end
    subplot(2,1,1)
    loglog(itervec,gap_track)
    hold on
    subplot(2,1,2)
    plot(x0,'k')
    hold on
    stem(x)
    x
            end
        end
    end
end











