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
maxiter = 1000;

%%
figure(1)
clf
for alpha = [100]
    alpha
    filename = sprintf('saved_data/CGM_alpha%f.mat',alpha);
    
    data = load(filename);
    gap_track = data.gap_track;
    obj_track = data.obj_track;
    train_track = data.train_track;
    test_track = data.test_track;
    
    
    
    subplot(3,4,1)
    loglog(gap_track)
    hold on
    title('gap')
    
    subplot(3,4,2)
    plot(obj_track)
    hold on
    title('obj')
    ylim([0,2])
    
    subplot(3,4,3)
    plot(train_track.P)
    hold on
    title('train P')
    
    subplot(3,4,4)
    plot(train_track.R)
    hold on
    title('train R')
    
    subplot(3,4,5)
    plot(train_track.acc)
    hold on
    title('train acc')
    
    subplot(3,4,6)
    plot(train_track.f1)
    hold on
    title('train f1')
    
    
    subplot(3,4,7)
    plot(test_track.P)
    hold on
    title('test P')
    
    subplot(3,4,8)
    plot(test_track.R)
    hold on
    title('test R')
    
    subplot(3,4,9)
    plot(test_track.acc)
    hold on
    title('test acc')
    
    subplot(3,4,10)
    plot(test_track.f1)
    hold on
    title('test f1')
    
    drawnow
end
%%

figure(2)
clf
for iter0 = [1000]%[0]%,10,100,1000]
    for p = [1.25,2,5,10]
        if p == 1.25
            rhovec = [.0001,]
        elseif p == 2
            rhovec = [1e-5]
        elseif p == 5
            rhovec = [1e-13,1e-10]
        elseif p == 10
            rhovec = [1e-26]
            
        end
        for rho  = rhovec
            
            if rho < 1e-5
                filename = sprintf('saved_data/PCGM_rho%e_p%f_iter0%d.mat',rho,p,iter0)
            else
                filename = sprintf('saved_data/PCGM_rho%f_p%f_iter0%d.mat',rho,p,iter0);
            end
            %                 try
            data = load(filename);
            gap_track = data.gap_track;
            obj_track = data.obj_track;
            train_track = data.train_track;
            test_track = data.test_track;
            %                 catch
            %                     continue
            %                 end
            
            [iter0,p,rho]
            
            
            
            subplot(3,4,1)
            loglog(gap_track)
            hold on
            title('gap')
            
            subplot(3,4,2)
            plot(obj_track)
            hold on
            title('obj')
            ylim([0,10])
            
            subplot(3,4,3)
            plot(train_track.P)
            hold on
            title('train P')
            
            subplot(3,4,4)
            plot(train_track.R)
            hold on
            title('train R')
            
            subplot(3,4,5)
            plot(train_track.acc)
            hold on
            title('train acc')
            
            subplot(3,4,6)
            plot(train_track.f1)
            hold on
            title('train f1')
            
            
            subplot(3,4,7)
            plot(test_track.P)
            hold on
            title('test P')
            
            subplot(3,4,8)
            plot(test_track.R)
            hold on
            title('test R')
            
            subplot(3,4,9)
            plot(test_track.acc)
            hold on
            title('test acc')
            
            subplot(3,4,10)
            plot(test_track.f1)
            hold on
            title('test f1')
            
            drawnow
            
            
        end
    end
end

%%

figure(3)
clf
figure(4)
clf
figure(5)
clf

figure(6)
clf

% for iter0 = [0,10,100,1000]
    for p =[1.25,2,5]%,10]
        rhovec = [1e-26,1e-14,1e-13,1e-12,1e-11,1e-5,1e-4,1e-3];%
        thvec = [1/8]%,1/4,1/2,3/4];
        if p == 1.25
            rhovec = [.01];
            thvec = [1/2,3/4];
        elseif p == 2
            rhovec = [1e-4];
            thvec = [1/8,1/4,1/2,3/4];
        elseif p == 5
            thvec = [1/8,1/4,1/2,3/4]
        elseif p == 10
%             rhovec = [1e-10];%
        end
        
        for th = thvec
            
            if p == 5
                if th == 1/8
                rhovec = [1e-10]
                iter0 = 1000;
                elseif th == 1/4
                    rhovec = [1e-10];
                    iter0 = 10;
                elseif th == 1/2
                    rhovec = [1e-10];
                    iter0 = 0;
                elseif th == 3/4
                    rhovec = [1e-10];
                    iter0 = 0;
                else
                    rhovec = []
                end
                
            end
%             if p == 1.25
%                 figure(3)
%             elseif p == 2
%                 figure(4)
%             elseif p == 5
%                 figure(5)
%             elseif p == 10
%                 figure(6)
%             end
            for rho = rhovec
                cbar = 1-th;
                
                %                 try
                if rho < 1e-5
                filename = sprintf('saved_data/RPCGM_rho%e_p%f_cbar%f_theta%f_iter0%d.mat',rho,p,cbar,th,iter0);
            else
                filename = sprintf('saved_data/RPCGM_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',rho,p,cbar,th,iter0);
            end
                data = load(filename);
                gap_track = data.gap_track;
                obj_track = data.obj_track;
                train_track = data.train_track;
                test_track = data.test_track;
                
                %                 catch
                %                     continue
                %                 end
                
                [iter0,th,rho,p,   train_track.P(end),    train_track.R(end),  train_track.f1(end),train_track.acc(end)]
                subplot(3,4,1)
                loglog(gap_track)
                hold on
                title('gap')
                
                subplot(3,4,2)
                plot(obj_track)
                hold on
                title('obj')
                
                subplot(3,4,3)
                semilogx(train_track.P)
                hold on
                title('train P')
                
                subplot(3,4,4)
                semilogx(train_track.R)
                hold on
                title('train R')
                
                subplot(3,4,5)
                semilogx(train_track.acc)
                hold on
                title('train acc')
                
                subplot(3,4,6)
                semilogx(train_track.f1)
                hold on
                title('train f1')
                
                
                subplot(3,4,7)
                semilogx(test_track.P)
                hold on
                title('test P')
                
                subplot(3,4,8)
                semilogx(test_track.R)
                hold on
                title('test R')
                
                subplot(3,4,9)
                semilogx(test_track.acc)
                hold on
                title('test acc')
                
                subplot(3,4,10)
                semilogx(test_track.f1)
                hold on
                title('test f1')
                
                drawnow
                
            end
        end
    end
% end

%% FOR PAPER


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
maxiter = 1000;

%%
figure(1)
clf


for p = [inf,  2,5]
    if p == 1.25
        iter0 = 1000;
        rho = .0001;
        filename = sprintf('saved_data/PCGM_rho%f_p%f_iter0%d.mat',rho,p,iter0);
    elseif p == 2
        iter0 = 1000;
        rho = 1e-5;
        filename = sprintf('saved_data/PCGM_rho%f_p%f_iter0%d.mat',rho,p,iter0);
    elseif p == 5
        rho = 1e-13;% = [1e-13,1e-10];
        iter0 = 1000;
        if rho < 1e-5
            filename = sprintf('saved_data/PCGM_rho%e_p%f_iter0%d.mat',rho,p,iter0)
        else
            filename = sprintf('saved_data/PCGM_rho%f_p%f_iter0%d.mat',rho,p,iter0);
        end
    elseif p == 10
        rho = [1e-26];
        iter0 = 1000;
        if rho < 1e-5
            filename = sprintf('saved_data/PCGM_rho%e_p%f_iter0%d.mat',rho,p,iter0)
        else
            filename = sprintf('saved_data/PCGM_rho%f_p%f_iter0%d.mat',rho,p,iter0);
        end
    elseif isinf(p)
        alpha = 100;
        filename = sprintf('saved_data/CGM_alpha%f.mat',alpha);
    end
    
    
    data = load(filename);
    gap_track = data.gap_track;
    obj_track = data.obj_track;
    train_track = data.train_track;
    test_track = data.test_track;
    
    
    figure(1)
    subplot(1,4,1)
    loglog(gap_track,'linewidth',2)
    hold on
%     title('gap')
    
    subplot(1,4,2)
    semilogx(obj_track,'linewidth',2)
    hold on
%     title('obj')
    ylim([0,2])
    
%     
%     subplot(3,2,3)
%     semilogx(train_track.acc)
%     hold on
%     title('train acc')
%     
%     figure(2)
    subplot(1,4,3)
    plot(train_track.f1,'linewidth',2)
    hold on
%     title('train f1')
    ylim([0,0.6])
    
%     plot(test_track.acc)
%     hold on
%     title('test acc')
    
    subplot(1,4,4)
    plot(test_track.f1,'linewidth',2)
    hold on
%     title('test f1')
    
    drawnow
end


iter0 = 1000;
% figure(3)
% clf




for p =[2,5]%,10]
    if p == 2
        rho = 1e-5;
        thvec = [1/2];
    elseif p == 5
        thvec = [1/2]
        rho = 1e-13;
    elseif p == 10
        %             rhovec = [1e-10];%
    end
    
    for th = thvec
        
     
        
        
        cbar = 1-th;
        
        %                 try
        
        
        if rho < 1e-5
            filename = sprintf('saved_data/RPCGM_rho%e_p%f_cbar%f_theta%f_iter0%d.mat',rho,p,cbar,th,iter0);
        else
            filename = sprintf('saved_data/RPCGM_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',rho,p,cbar,th,iter0);
        end
        
        
    data = load(filename);
    gap_track = data.gap_track;
    obj_track = data.obj_track;
    train_track = data.train_track;
    test_track = data.test_track;
    
    %                 catch
    %                     continue
    %                 end
    
    [iter0,th,rho,p,   train_track.P(end),    train_track.R(end),  train_track.f1(end),train_track.acc(end)]
    figure(1)
    subplot(1,4,1)
    loglog(gap_track,'linewidth',2)
    hold on
    ylabel('gap')
    xlabel('iterations')
    
    subplot(1,4,2)
    plot(obj_track,'linewidth',2)
    hold on
    ylabel('obj')
    xlabel('iterations')
    drawnow
    
    
%     subplot(3,2,3)
%     semilogx(train_track.acc)
%     hold on
%     title('train acc')
%     
%     figure(2)
    subplot(1,4,3)
    plot(train_track.f1,'linewidth',2)
    hold on
    ylabel('train F1')
    ylim([0,0.6]) 
    xlabel('iterations')
%     legend('CGM','PCGM-1','PCGM-2','RPCGM-1','RPCGM-2','location','southeast')
    
%     figure(2)
%     subplot(1,2,1)
%     plot(test_track.acc)
%     hold on
%     title('test acc')
    
    subplot(1,4,4)
    plot(test_track.f1,'linewidth',2)
    hold on
    ylabel('test F1')
    ylim([0,0.6])
    xlabel('iterations')
    legend('CGM','PCGM-1','PCGM-2','RPCGM-1','RPCGM-2','location','southeast')
    drawnow
    
    
    end
end

