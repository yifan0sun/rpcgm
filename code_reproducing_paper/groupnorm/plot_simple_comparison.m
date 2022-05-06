clc
clear

%%
pulsewidth = 10;
noise_variance = 0;
trial = 1;
n = 100;
m = 100;
num_pulses = 3;

% offset = randi(n-pulsewidth,1,num_pulses)
% x0 = zeros(n,1);
% for k = 1:num_pulses
% x0(offset(k)+1:offset(k)+pulsewidth) = 1;
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

figure(1)
clf
x0 = x0 / norm(x0,1);
for k = 1:6
    subplot(2,6,k)
    plot(x0,'k')
    hold on
end
figk = 0;
for alpha = 5%[1,2,5,10]
    alpha
    filename = sprintf('saved_data/CGM_m%d_n%d_trial%d_numpulses%d_pulsewidth%d_noise%f_alpha%f_iter0%d.mat',m,n,trial,num_pulses, pulsewidth,noise_variance,alpha,iter0);
    
    data = load(filename);
    gap_track = data.gap_track;
    x = data.x;
    supp_err = data.supp_err;
    itervec = data.itervec;
    supp = data.supp;
    screen = data.screen;
    ztrack = data.ztrack;
    xtrack = data.xtrack;
    
    
    
    
    
    final_supp = zeros(length(groups),1);
    for k = 1:length(groups)
        g = groups{k};
        if norm(supp(end,k)) > 0
            final_supp(k) = 1;
        end
    end
    
    if sum(final_supp) <= 1
        fprintf('alpha %f too low\n', alpha)
        continue
    end
    
    screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp')),2);
    
    
    figure(1)
    
    subplot(2,4,1)
    x = xtrack(itervec==3,:);
    x = x/ norm(x,1);
    stem(x)
    %     ylim([-2,2])
    
    subplot(2,4,2)
    x = xtrack(itervec==10,:);
    x = x/ norm(x,1);
    stem(x)
    %     ylim([-2,2])
    
    subplot(2,4,3)
    x = xtrack(itervec==30,:);
    x = x/ norm(x,1);
    stem(x)
    %     ylim([-2,2])
    
    subplot(2,4,4)
    x = xtrack(itervec==100,:);
    x = x/ norm(x,1);
    stem(x)
    %     ylim([-2,2])
    
    %
    
    subplot(2,4,5)
    z = abs(ztrack(itervec==3,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    subplot(2,4,6)
    z = abs(ztrack(itervec==10,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    
    
    subplot(2,4,7)
    z = abs(ztrack(itervec==30,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    
    subplot(2,4,8)
    z = abs(ztrack(itervec==100,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    
end


%%



figure(2)
clf
x0 = x0 / max(abs(x0));
for k = 1:6
    subplot(2,6,k)
    plot(x0,'k')
    hold on
end
figk = 0;




for p = [2]%,2,5]
    
    if p == 1.25
        rho = 5;
    elseif p == 2
        rho = .1;
    elseif p == 5
        rho = .01;
    end
    [p,rho]
    filename = sprintf('saved_data/PCGM_m%d_n%d_trial%d_numpulses%d_pulsewidth%d_noise%f_rho%f_p%f_iter0%d.mat',m,n,trial,num_pulses, pulsewidth,noise_variance,rho,p,iter0);
    
    
    
    
    data = load(filename);
    gap_track = data.gap_track;
    x = data.x;
    supp_err = data.supp_err;
    itervec = data.itervec;
    supp = data.supp;
    screen = data.screen;
    ztrack = data.ztrack;
    xtrack = data.xtrack;
    
    
    
    
    
    final_supp = zeros(length(groups),1);
    for k = 1:length(groups)
        g = groups{k};
        if norm(supp(end,k)) > 0
            final_supp(k) = 1;
        end
    end
    
    if sum(final_supp) <= 1
        fprintf('alpha %f too low\n', alpha)
        %         continue
    end
    
    screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp')),2);
    
    
    figure(2)
    subplot(2,6,1)
    x = xtrack(itervec==3,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    
    subplot(2,6,2)
    x = xtrack(itervec==10,:);
    x = x/ max(abs(x))
    stem(x)
    ylim([-2,2])
    
    subplot(2,6,3)
    x = xtrack(itervec==30,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    
    subplot(2,6,4)
    x = xtrack(itervec==100,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    
    subplot(2,6,5)
    x = xtrack(itervec==300,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    
    subplot(2,6,6)
    x = xtrack(itervec==1000,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    
    %
    
    subplot(2,6,7)
    z = abs(ztrack(itervec==3,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    subplot(2,6,8)
    z = abs(ztrack(itervec==10,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    
    
    subplot(2,6,9)
    z = abs(ztrack(itervec==30,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    
    subplot(2,6,10)
    z = abs(ztrack(itervec==100,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    
    subplot(2,6,11)
    z = abs(ztrack(itervec==300,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    
    subplot(2,6,12)
    z = abs(ztrack(itervec==1000,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    %
end




%%


iter0 = 100;

figure(3)
clf
x0 = x0 / max(abs(x0));
for k = 1:6
    subplot(2,6,k)
    plot(x0,'k')
    hold on
end
figk = 0;


theta = 1/2;
for p = [2]%,2,5]
    
    rho =.01;
    [p,rho]
    cbar = 1-theta;
    filename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_numpulses%d_pulsewidth%d_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,num_pulses, pulsewidth,noise_variance,rho,p,cbar,theta,iter0);
    
    
    
    [theta,rho,p]
    
    data = load(filename);
    gap_track = data.gap_track;
    x = data.x;
    supp_err = data.supp_err;
    itervec = data.itervec;
    supp = data.supp;
    screen = data.screen;
    ztrack = data.ztrack;
    xtrack = data.xtrack;
    
    
    
    
    
    final_supp = zeros(length(groups),1);
    for k = 1:length(groups)
        g = groups{k};
        if norm(supp(end,k)) > 0
            final_supp(k) = 1;
        end
    end
    
    if sum(final_supp) <= 1
        fprintf('alpha %f too low\n', alpha)
        %         continue
    end
    
    screen_err = mean((screen ~= (ones(size(screen,1),1)*final_supp')),2);
    
    
    figure(3)
    subplot(2,6,1)
    x = xtrack(itervec==3,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    title('t = 3')
    ylabel('x')
    axis tight
    
    subplot(2,6,2)
    x = xtrack(itervec==10,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    title('t = 10')
    ylabel('x')
    axis tight
    
    subplot(2,6,3)
    x = xtrack(itervec==30,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    title('t = 30')
    ylabel('x')
    axis tight
    
    subplot(2,6,4)
    x = xtrack(itervec==100,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    title('t = 100')
    ylabel('x')
    axis tight
    
    subplot(2,6,5)
    x = xtrack(itervec==300,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    title('t = 300')
    ylabel('x')
    axis tight
    
    subplot(2,6,6)
    x = xtrack(itervec==1000,:);
    x = x/ max(abs(x));
    stem(x)
    ylim([-2,2])
    title('t = 1000')
    ylabel('x')
    legend('$x_0$','$x^{(t)}$','interpreter','latex')
    axis tight
    %
    
    subplot(2,6,7)
    z = abs(ztrack(itervec==3,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 3')
    
    axis tight
    subplot(2,6,8)
    z = abs(ztrack(itervec==10,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 10')
    
    subplot(2,6,9)
    z = abs(ztrack(itervec==30,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 30')
    
    axis tight
    subplot(2,6,10)
    z = abs(ztrack(itervec==100,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 100')
    
    subplot(2,6,11)
    z = abs(ztrack(itervec==300,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 300')
    
    axis tight
    subplot(2,6,12)
    z = abs(ztrack(itervec==1000,:));
    nz = zeros(length(groups),1);
    for k = 1:length(groups)
        nz(k) = norm(z(groups{k}));
    end
    stem(nz)
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 1000')
    
    axis tight
    
end









