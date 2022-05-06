clc
clear


%%
noise_variance = 0;
trial = 1;
n = 100;
m = 1000;
num_edges = 4;

loadfilename = sprintf('saved_data/problem_m%d_n%d_trial%d_num_edges%d_noise%f.mat',m,n,trial,num_edges,noise_variance);
data = load(loadfilename);
A = data.A;
b = data.b;
x0 = data.x0;
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


theta = 3/4;
for p = [2]%,2,5]
    
    rho =.25;
    [p,rho]
    cbar = 1-theta;
    filename = sprintf('saved_data/RPCGM_m%d_n%d_trial%d_num_edges%d_noise%f_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,n,trial,num_edges, noise_variance,rho,p,cbar,theta,iter0);
    
    
    
    [theta,rho,p]
    
    data = load(filename);
    gap_track = data.gap_track;
    x = data.x;
    itervec = data.itervec;
    ztrack = data.ztrack;
    xtrack = data.xtrack;
    
    
    
    
    
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
    z = (ztrack(itervec==3,:));
    z = z - mean(z);
%     z = cumsum(z)
    stem(cumsum(z))
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 3')
    
    axis tight
    subplot(2,6,8)
    z = (ztrack(itervec==10,:));
    z = z - mean(z);
   
    stem(cumsum(z))
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 10')
    
    subplot(2,6,9)
    z = (ztrack(itervec==30,:));
    z = z - mean(z);
   
    stem(cumsum(z))
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 30')
    
    axis tight
    subplot(2,6,10)
    z = (ztrack(itervec==100,:));
    z = z - mean(z);
   
    stem(cumsum(z))
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 100')
    
    subplot(2,6,11)
    z = (ztrack(itervec==300,:));
    z = z - mean(z);
   
    stem(cumsum(z))
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 300')
    
    axis tight
    subplot(2,6,12)
    z = (ztrack(itervec==1000,:));
    z = z - mean(z);
    
    stem(cumsum(z))
    ylabel('$\|z_G\|_2$','interpreter','latex')
    title('t = 1000')
    
    axis tight
    
end









