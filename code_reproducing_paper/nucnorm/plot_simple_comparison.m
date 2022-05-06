clc
clear

%%

co = get(gca,'colororder');
fignum = 1;
ntrials = 1;
for rank_percent = .1%[0.1,0.5]
    
    %     figure(fignum)
    %     clf
    %     fignum = fignum + 1;
    %     '----'
    for noise_variance = 1%[0,.01,.1,1]
        for m = 100%[10,25,50,100,250,500,1000]
            
            r = round(rank_percent*m);
            
            if rank_percent == 0.1
                alphavec = [10000];
            elseif rank_percent == 0.5
                alphavec = 10000;
            end
            figure(1)
            clf
            
            for alpha = alphavec
                for trial = 1:ntrials
                    filename = sprintf('saved_data/CGM_m%d_rank%f_noise%f_trial%d_alpha%f.mat',m,rank_percent,noise_variance,trial,alpha);
                    
                    data = load(filename);
                    X = data.X;
                    gap_track = data.gap_track;
                    obj_track = data.obj_track;
                    rec_err = data.rec_err;
                    sing_rec_err = data.sing_rec_err;
                    
                    subplot(1,4,1)
                        loglog(gap_track,'k','linewidth',2)
                        xlim([0,200])
                        hold on
                        subplot(1,4,2)
                        loglog(obj_track,'k','linewidth',2)
%                         xlim([0,200])
                        hold on
                     
                        subplot(1,4,3)
                        loglog(rec_err,'k','linewidth',2)
%                         xlim([0,200])
                        hold on
                        
                        subplot(1,4,4)
                        
                        sx = svd(X);
                        sx = sx / max(sx);
                        sx(1:5)
                        plot(sx,'k','linewidth',2)
                        hold on
                        
                    drawnow
                end
            end
            
            
            
            
            %%
            
            %                 figure(2)
            %                 clf
            co_idx = 1;
            for iter0 = [100]%,10,100,1000]%[0,10,100,1000]
                for p = [2,5,10]
                    %                     rhovec = [0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000];
                    if p == 1.1
                        rhovec = [1e-10,1e-5,1];
                    elseif p == 1.25
                        rhovec = [1];
                    elseif p == 2
                        rhovec = [1e-5];
                    elseif p == 5
                        rhovec = [1e-14];
                    elseif p == 10
                        rhovec = [1e-30];
                    end
%                         rhovec =  [1e-35,1e-30,1e-25,1e-20,1e-15,1e-10];
                    %                         rhovec = [1e-15,1e-10,1e-5,1,1e5,1e10];
                    for rho = rhovec
                        for trial = 1:ntrials
                            if rho < 1e-5 || rho > 1e5
                                filename = sprintf('saved_data/PCGM_m%d_rank%f_noise%f_trial%d_rho%e_p%f_iter0%d.mat',m,r, noise_variance, trial,rho,p,iter0);
                            else
                                filename = sprintf('saved_data/PCGM_m%d_rank%f_noise%f_trial%d_rho%f_p%f_iter0%d.mat',m,r, noise_variance, trial,rho,p,iter0);
                                
                            end
                            data = load(filename);
                            X = data.X;
                            gap_track = data.gap_track;
                            rec_err = data.rec_err;
                            obj_track = data.obj_track;
                        end
                        
                        
                        subplot(1,4,1)
                        loglog(gap_track,'color',co(co_idx,:),'linewidth',2)
                        xlim([0,200])
                        hold on
                        subplot(1,4,2)
                        loglog(obj_track,'color',co(co_idx,:),'linewidth',2)
%                         xlim([0,200])
                        hold on
                     
                        subplot(1,4,3)
                        loglog(rec_err,'color',co(co_idx,:),'linewidth',2)
%                         xlim([0,200])
                        hold on
                        
                        subplot(1,4,4)
                        
                        sx = svd(X);
                        sx = sx / max(sx);
                        sx(1:5)
                        plot(sx,'color',co(co_idx,:),'linewidth',2)
                        hold on
                        
                        % ylim([1e-5,10])
                        drawnow
                        co_idx = co_idx+ 1;
                        
                    end
                    
                end
            end
            
            %%
            
            
            %              figure(3)
            %                 clf
            
            co_idx2 = 1;
            for iter0 = [100]%,10,100,1000]%,10,100,1000]
                
                
                
                for theta = [1/4,1/2,3/4]%[1/8,1/4,1/2,3/4]
                    for p =[2]%,5,10]
                        
                        if p == 1.25
                            rhovec = [1];
                        elseif p == 2
                            rhovec = [1e-5];
                        elseif p == 5
                            rhovec = [1e-10];
                        elseif p == 10
                            rhovec = [1e-15];
                        end
                            rhovec = [1];
                           
                        for rho = rhovec
                            cbar = 1-theta;
                            
                            for trial = 1%:ntrials
                                
                                if rho < 1e-5 || rho > 1e5
                                    filename = sprintf('saved_data/RPCGM_m%d_rank%f_noise%f_trial%d_rho%e_p%f_cbar%f_theta%f_iter0%d.mat',m,r, noise_variance,trial,rho,p,cbar,theta,iter0);
                                else
                                    filename = sprintf('saved_data/RPCGM_m%d_rank%f_noise%f_trial%d_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,r, noise_variance,trial,rho,p,cbar,theta,iter0);
                                    
                                end
                                data = load(filename);
                                
                                
                                X = data.X;
                                gap_track = data.gap_track
                                rec_err = data.rec_err;
                                sing_rec_err = data.sing_rec_err;
                                obj_track = data.obj_track;
                                
                            end
                            
                            
                            
                            subplot(1,4,1)
                            loglog(gap_track,'color',co(co_idx2,:),'linewidth',2,'linestyle','--')
%                             xlim([0,200])
                            hold on
                            subplot(1,4,2)
                            loglog(obj_track,'color',co(co_idx2,:),'linewidth',2,'linestyle','--')
                            %                         xlim([0,200])
                            hold on
                            
                            subplot(1,4,3)
                            loglog(rec_err,'color',co(co_idx2,:),'linewidth',2,'linestyle','--')
                            %                         xlim([0,200])
                            hold on
                            
                            subplot(1,4,4)
                            
                            sx = svd(X);
                            sx = sx / max(sx);
                            sx(1:5)
                            plot(sx,'color',co(co_idx2,:),'linewidth',2,'linestyle','--')
                            hold on
                            
                            legend('CGM','PCGM-1','PCGM-2','PCGM-3','RPCGM-1','RPCGM-2','RPCGM-3')
                            drawnow
                            
                            
                            
                            co_idx2 = co_idx2+ 1;
                            
                            
                            
                            
                            
                        end
                    end
                end
            end
        end
    end
end
