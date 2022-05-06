clc
clear


%%

n = 100;
splevel = 0.1;
noise_variance = 0;
ntrials = 10;

% m = [10,25,50,100,250,500,1000,2500,5000,10000]
%
m = 100;

trial = 1;
%%
iter0 = 0;
figure(1)
clf

for alpha = [1,2,5]
    loadfilename = sprintf('saved_data/CGM_m%d_n%d_trial%d_splevel%f_noise%f_alpha%f_iter0%d.mat',m,n,trial,splevel,noise_variance,alpha,iter0);
    data = load(loadfilename);
    x = data.x;
    gap_track = data.gap_track;
    supp_err = data.supp_err;
    itervec = data.itervec;
    
    
    supp = data.supp;
    screen = data.screen;
    
    final_supp = supp(end,:);
    
    if sum(final_supp) <= 1
        fprintf('alpha %f too low\n', alpha)
        continue
    end
    
    screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp)),2);
    
    
    subplot(4,3,1)
    loglog(itervec,gap_track,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('gap')
    legend('$\alpha  = 1$','$\alpha = 2$','$\alpha = 5$','$\alpha = 10$','interpreter','latex','location','southwest')
    
    subplot(4,3,2)
    semilogx(itervec,supp_err,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('support error')
    title('CGM')
    
    subplot(4,3,3)
    semilogx(itervec, screen_err,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('screen error')
end


%%
iter0 = 1000;
p = 1.25;

rhovec = [.5,1,2];

for rho = rhovec
    
    loadfilename = sprintf('saved_data/PCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,iter0);
    data = load(loadfilename);
    x = data.x;
    gap_track = data.gap_track;
    supp_err = data.supp_err;
    itervec = data.itervec;
    supp = data.supp;
    screen = data.screen;
    final_supp = supp(end,:);
    if sum(final_supp) <= 1
        fprintf('rho %f too high\n', rho)
        continue
    end
    screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp)),2);
    
    
    
    
    subplot(4,3,4)
    loglog(itervec,gap_track,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('gap')
    legend('$\rho  = 0.1$','$\rho = 1$','$\rho = 10$','interpreter','latex','location','southwest')
    
    subplot(4,3,5)
    semilogx(itervec,supp_err,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('support error')
    title('PCGM, $p=1.25$','interpreter','latex')
    
    
    subplot(4,3,6)
    semilogx(itervec, screen_err,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('screen error')
    
    
end
%%
iter0 = 1000;
rhovec = [.01,.1,1.];

p = 2;
for rho = rhovec
    loadfilename = sprintf('saved_data/PCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,iter0);
    data = load(loadfilename);
    x = data.x;
    gap_track = data.gap_track;
    supp_err = data.supp_err;
    itervec = data.itervec;
    supp = data.supp;
    screen = data.screen;
    final_supp = supp(end,:);
    if sum(final_supp) <= 1
        fprintf('rho %f too high\n', rho)
        continue
    end
    screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp)),2);
    
    subplot(4,3,7)
    loglog(itervec,gap_track,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('gap')
    legend('$\rho  = 0.1$','$\rho = 0.2$','$\rho = 0.5$','interpreter','latex','location','southwest')
    
    subplot(4,3,8)
    semilogx(itervec,supp_err,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('support error')
    title('PCGM, $p=2$','interpreter','latex')
    
    subplot(4,3,9)
    semilogx(itervec, screen_err,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('screen error')
end

%%
iter0 = 1000;
p = 5;

rhovec = [.001,.01,.1];

for rho = rhovec
    loadfilename = sprintf('saved_data/PCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,iter0);
    data = load(loadfilename);
    x = data.x;
    gap_track = data.gap_track;
    supp_err = data.supp_err;
    itervec = data.itervec;
    supp = data.supp;
    screen = data.screen;
    final_supp = supp(end,:);
    if sum(final_supp) <= 1
        fprintf('rho %f too high\n', rho)
        continue
    end
    screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp)),2);
    
    
    subplot(4,3,10)
    loglog(itervec,gap_track,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('gap')
    legend('$\rho  = 0.1$','$\rho = 0.2$','$\rho = 0.5$','interpreter','latex','location','southwest')
    
    subplot(4,3,11)
    semilogx(itervec,supp_err,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('support error')
    title('PCGM, $p=5$','interpreter','latex')
    
    subplot(4,3,12)
    semilogx(itervec, screen_err,'linewidth',2)
    hold on
    xlabel('iterations')
    ylabel('screen error')
end

%%


figure(1)
clf
iter0 = 1000;
p = 1.25;

for theta = [1/8,1/4,1/2,3/4]
    
    
    for rho = [.5]
        cbar = 1-theta;
        
        
        
        
        loadfilename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,cbar,theta,iter0);
        data = load(loadfilename);
        x = data.x;
        gap_track = data.gap_track;
        supp_err = data.supp_err;
        itervec = data.itervec;
        supp = data.supp;
        screen = data.screen;
        final_supp = supp(end,:);
            if sum(final_supp) <= 1
                fprintf('rho %f too high\n', rho)
%                 continue
            end
        screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp)),2);
        
        
        subplot(4,3,1)
        loglog(itervec,gap_track,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('gap')
        legend('$\theta  = 1/8$','$\theta = 1/4$','$\theta = 1/2$','$\theta = 3/4$','interpreter','latex','location','southwest')
        
        subplot(4,3,2)
        semilogx(itervec,supp_err,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('support error')
        title('RPCGM, $p=1.25$','interpreter','latex')
        
        subplot(4,3,3)
        semilogx(itervec, screen_err,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('screen error')
        
        
        
    end
    
end

%%



iter0 = 1000;
p = 2;

for theta = [1/8,1/4,1/2,3/4]
    
    
    for rho = [.1]
        cbar = 1-theta;
        
        
        
        
        loadfilename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,cbar,theta,iter0);
        data = load(loadfilename);
        x = data.x;
        gap_track = data.gap_track;
        supp_err = data.supp_err;
        itervec = data.itervec;
        supp = data.supp;
        screen = data.screen;
        final_supp = supp(end,:);
            if sum(final_supp) <= 1
                fprintf('rho %f too high\n', rho)
        %         continue
            end
        screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp)),2);
        
        
        subplot(4,3,4)
        loglog(itervec,gap_track,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('gap')
        legend('$\theta  = 1/8$','$\theta = 1/4$','$\theta = 1/2$','$\theta = 3/4$','interpreter','latex','location','southwest')
        
        subplot(4,3,5)
        semilogx(itervec,supp_err,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('support error')
        title('RPCGM, $p=2$','interpreter','latex')
        
        subplot(4,3,6)
        semilogx(itervec, screen_err,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('screen error')
        
        
        
    end
    
end

%%

p = 5;

for theta = [1/8,1/4,1/2,3/4]
    
    
    for rho = [.01]
        cbar = 1-theta;
        

 loadfilename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,cbar,theta,iter0);
        data = load(loadfilename);
        x = data.x;
        gap_track = data.gap_track;
        supp_err = data.supp_err;
        itervec = data.itervec;
        supp = data.supp;
        screen = data.screen;
        final_supp = supp(end,:);
            if sum(final_supp) <= 1
                fprintf('rho %f too high\n', rho)
        %         continue
            end
        screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp)),2);
        
        
        subplot(4,3,7)
        loglog(itervec,gap_track,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('gap')
        legend('$\theta  = 1/8$','$\theta = 1/4$','$\theta = 1/2$','$\theta = 3/4$','interpreter','latex','location','southwest')

        subplot(4,3,8)
        semilogx(itervec,supp_err,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('support error')
        title('RPCGM, $p=5$','interpreter','latex')
        
        subplot(4,3,9)
        semilogx(itervec, screen_err,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('screen error')
        
        
        
    end
    
end

%%

p = 10;

for theta = [1/8,1/4,1/2,3/4]
    
    
    for rho = [.001]
        cbar = 1-theta;
        
        
        loadfilename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_splevel%f_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,splevel,noise_variance,rho,p,cbar,theta,iter0);
        data = load(loadfilename);
        x = data.x;
        gap_track = data.gap_track;
        supp_err = data.supp_err;
        itervec = data.itervec;
        supp = data.supp;
        screen = data.screen;
        final_supp = supp(end,:);
        if sum(final_supp) <= 1
            fprintf('rho %f too high\n', rho)
            %         continue
        end
        screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp)),2);
        
        
        subplot(4,3,10)
        loglog(itervec,gap_track,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('gap')
        legend('$\theta  = 1/8$','$\theta = 1/4$','$\theta = 1/2$','$\theta = 3/4$','interpreter','latex','location','southwest')

        subplot(4,3,11)
        semilogx(itervec,supp_err,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('support error')
        title('RPCGM, $p=10$','interpreter','latex')
        
        subplot(4,3,12)
        semilogx(itervec, screen_err,'linewidth',2)
        hold on
        xlabel('iterations')
        ylabel('screen error')
        
        
        
    end
    
end