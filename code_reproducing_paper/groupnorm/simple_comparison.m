clc
clear

%%
pulsewidth = 10;
noise_variance = 0;
trial = 1;
n = 100;
m = 100;
num_pulses = 3;
% 
% offset = randi(n-pulsewidth,1,num_pulses)
% x0 = zeros(n,1);
% for k = 1:num_pulses
% x0(offset(k)+1:offset(k)+pulsewidth) = 1;
% end
% x0 = x0.*randn(n,1);
% figure(1)
% clf
% plot(x0)
% 
% 
% A = randn(m,n);
% 
% b = A*x0+noise_variance*randn(m,1);
% 
% savefilename = sprintf('saved_data/problem_m%d_n%d_trial%d_numpulses%d_pulsewidth%d_noise%f.mat',m,n,trial,num_pulses, pulsewidth,noise_variance)
% save(savefilename,'A','b','x0','offset')


%%

iter0 = 0;

trial
m
loadfilename = sprintf('saved_data/problem_m%d_n%d_trial%d_numpulses%d_pulsewidth%d_noise%f.mat',m,n,trial,num_pulses, pulsewidth,noise_variance);
data = load(loadfilename);
A = data.A;
b = data.b;
active_groups = data.offset;
x0 = data.x0;


for k = 1:(n-pulsewidth)
    groups{k} = k+1:k+pulsewidth;
end

smoothobj = @(x)( sum((A*x-b).^2)/m/2 );
smoothgrad = @(x)(A'*(A*x-b)/m);

L = norm(A'*A,2);


%%

maxiter = 10001;

for alpha = [5]
    alpha
    filename = sprintf('saved_data/CGM_m%d_n%d_trial%d_numpulses%d_pulsewidth%d_noise%f_alpha%f_iter0%d.mat',m,n,trial,num_pulses, pulsewidth,noise_variance,alpha,iter0);
    try
        
        data = load(filename);
        gap_track = data.gap_track;
        x = data.x;
        supp_err = data.supp_err;
        itervec = data.itervec;
        supp = data.supp;
        screen = data.screen;
    catch
        [x,supp_err, gap_track, supp, screen, xtrack,ztrack,itervec] = CGM(n,groups, active_groups,maxiter, alpha, smoothgrad, 'groupnorm',iter0,L);
        save(filename,'x','gap_track','supp','supp_err','screen', 'xtrack','ztrack','itervec')
    end
    
    
    
end



%%

maxiter = 10001;
for p = [2,5]%,2,5]
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
    for rho = rhovec
        [p,rho]
        filename = sprintf('saved_data/PCGM_m%d_n%d_trial%d_numpulses%d_pulsewidth%d_noise%f_rho%f_p%f_iter0%d.mat',m,n,trial,num_pulses, pulsewidth,noise_variance,rho,p,iter0);
        
            

        try
            
            data = load(filename);
            gap_track = data.gap_track;
            x = data.x;
            supp_err = data.supp_err;
            itervec = data.itervec;
            supp = data.supp;
            screen = data.screen;
        catch
            [x,supp_err, gap_track, supp, screen, xtrack,ztrack, itervec] = PCGM(n,groups, active_groups,maxiter, rho, smoothgrad,  'groupnorm',p,iter0,L);
            save(filename,'x','gap_track','supp','supp_err','screen', 'xtrack','ztrack','itervec')
        end
    end
end



%%
% maxiter = 100
for iter0 = [100]
    for p = [2]
        for theta = 1/2%[1/8,1/4,1/2,3/4,7/8]
            for rho = .01
                [theta,rho,p]
                cbar = 1-theta;
                filename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_numpulses%d_pulsewidth%d_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,num_pulses, pulsewidth,noise_variance,rho,p,cbar,theta,iter0);

                try
                    
                    data = load(filename);
                    gap_track = data.gap_track;
                    x = data.x;
                    supp_err = data.supp_err;
                    itervec = data.itervec;
                    supp = data.supp;
                    screen = data.screen;
                catch
                    [x,supp_err, gap_track, supp, screen, xtrack,ztrack,itervec] = RPCGM(n,groups, active_groups,maxiter, rho, theta, cbar, smoothgrad,  'groupnorm',p,iter0,L,0);
        save(filename,'x','gap_track','supp','supp_err','screen', 'xtrack','ztrack','itervec')
                end
                figure(1)
                clf
                plot(gap_track)
            end
        end
    end
end











