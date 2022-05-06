clc
clear

%%

% for trial = 1:10
%     trial
%     for rank_percent = [0.1,0.5]
%         for noise_variance = [0,.01,.1,1]
%             for m = [10,25,50,100,250,500,1000]
%                 filename = sprintf('saved_data/problem_m%d_rank%f_noise%f_trial%d.mat',m,rank_percent,noise_variance,trial);
%
%
%                 try
%                     data = load(filename);
%                     r = data.r;
%                     U = data.U;
%                     V = data.V;
%                     Y = data.Y;
%                 catch
%                     'data forming'
%                     r = round(rank_percent*m);
%                     U = randn(m,r);
%                     V = randn(m,r);
%                     Y = U*V' + randn(m)*noise_variance;
%
%                     save(filename,'r','U','V','Y')
%                 end
%
%
%
%             end
%         end
%     end
% end

%%

maxiter = 1000;

for trial = 1%$:10
    trial
    for rank_percent = 0.1%[0.1,0.5]
        for noise_variance = 1%[0,.01,.1,1]
            for m = 100%[10,25,50,100,250,500,1000]
                [rank_percent, noise_variance, m]
                loadfilename = sprintf('saved_data/problem_m%d_rank%f_noise%f_trial%d.mat',m,rank_percent,noise_variance,trial);
                data = load(loadfilename);
                
                r = data.r
                U = data.U;
                V = data.V;
                Y = data.Y;
                
                smoothobj = @(X)( sum(sum((X-Y).^2))/r/2/m );
                smoothgrad = @(X)((X-Y)/m/r);
                
                
                %
                for alpha = [.001,.01,.1,1,10,100,1000,10000,1e6]
                    filename = sprintf('saved_data/CGM_m%d_rank%f_noise%f_trial%d_alpha%f.mat',m,rank_percent,noise_variance,trial,alpha);
                    try
                        
                        data = load(filename);
                        X = data.X;
                        gap_track = data.gap_track;
                        obj_track = data.obj_track;
                        rec_err = data.rec_err;
                        
                    catch
                        [X,rec_err, sing_rec_err,gap_track,obj_track] = CGM(U,V,maxiter, alpha, smoothgrad, smoothobj, 'nucnorm');
                        save(filename,'X','gap_track','obj_track','rec_err','sing_rec_err')
                    end
                    
                    
                end
                
                        
                %%
%                 for iter0 = [100]%[0,10,100,1000]%[0,10,100,1000]
%                     for p = [10]%[1.1,1.25,2,5,10]
%                         rhovec =  [1e-30,1e-35,1e-25,1e-20,1e-15,1e-10];
%                     %         
%                         for rho = rhovec% [1e-15,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1,1e5,1e10]
%                             
%                             [rho,p,iter0]
%                             if rho < 1e-5 || rho > 1e5
%                                 filename = sprintf('saved_data/PCGM_m%d_rank%f_noise%f_trial%d_rho%e_p%f_iter0%d.mat',m,r, noise_variance, trial,rho,p,iter0);
%                             else
%                                 filename = sprintf('saved_data/PCGM_m%d_rank%f_noise%f_trial%d_rho%f_p%f_iter0%d.mat',m,r, noise_variance, trial,rho,p,iter0);
%                                 
%                             end
%                             
%                             try
%                                 
%                                 data = load(filename);
%                                 X = data.X;
%                                 gap_track = data.gap_track;
%                                 rec_err = data.left_rec_err;
%                                 sing_rec_err = data.sing_rec_err;
%                                 obj_track = data.obj_track;
%                             catch
%                                 
%                                 [X,rec_err, sing_rec_err,gap_track,obj_track] = PCGM(U,V,maxiter, rho, smoothgrad, smoothobj, 'nucnorm',p,iter0);
%                                 
%                                 save(filename,'X','gap_track','obj_track','rec_err','sing_rec_err')
%                             end
%                             
%                             
%                         end
%                         
%                     end
%                 end
                
                
                %%
                
                
                
                
                for iter0 = [100]%[0,10,100,1000]%[0,10,100,1000]
                    for theta = [1/4,1/2,3/4]
                        for p = 2%1.25,2,5,10]
                            p
                            if p == 1.25
                                rhovec = [1];
                            elseif p == 2
                                rhovec = [1e-5];
                            elseif p == 5
                                rhovec = [1e-10];
                            elseif p == 10
                                rhovec = [1e-15];
                            end
                            
                                             
                        
                        
                            %                             rhovec = [1e-15,1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,.1,1,10,100,1000,1e5];
                            
                            rhovec = [1e-15,1e-10,1e-5,1];
                            rhovec = [1e-6,1e-5,1e-4];
%                             rhovec = 1e10;
                            rhovec = 10;
                            for rho = rhovec%[.0001,.001,.01,.1,1]
                                cbar = 1-theta;
                                
                                
                                
                                if rho < 1e-5 || rho > 1e5
                                    filename = sprintf('saved_data/RPCGM_m%d_rank%f_noise%f_trial%d_rho%e_p%f_cbar%f_theta%f_iter0%d.mat',m,r, noise_variance,trial,rho,p,cbar,theta,iter0);
                                else
                                    filename = sprintf('saved_data/RPCGM_m%d_rank%f_noise%f_trial%d_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',m,r, noise_variance,trial,rho,p,cbar,theta,iter0);
                                    
                                end
                                try
                                    sadf
                                    data = load(filename);
                                    X = data.X;
                                    gap_track = data.gap_track;
                                    rec_err = data.left_rec_err;
                                    sing_rec_err = data.sing_rec_err;
                                    obj_track = data.obj_track;
                                catch
                                    
                                    [X,rec_err, sing_rec_err,gap_track,obj_track] = ...
                                        RPCGM(U,V,maxiter, rho,  theta, cbar, smoothgrad, smoothobj, 'nucnorm',p,iter0);
                                    
                                    save(filename,'X','gap_track','obj_track','rec_err','sing_rec_err')
                                end
                                gap_track
                                filename
                                
                             asdf
                            end
                        end
                    end
                end
                
            end
        end
    end
end









