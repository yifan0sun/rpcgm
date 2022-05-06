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
for alpha = 100%[.01,.1,1,10,100,1000]
    alpha
    filename = sprintf('saved_data/CGM_alpha%f.mat',alpha);
    try
        data = load(filename);
        gap_track = data.gap_track;
        obj_track = data.obj_track;
        train_track = data.train_track;
        test_track = data.test_track;
    catch
        
        
        [theta,obj_track, train_track, test_track, gap_track] ...
            = CGM(n,maxiter, alpha, get_trainobj, smoothgrad,  '1norm');
        
        
        save(filename,'theta','obj_track', 'gap_track','train_track','test_track')
    end
    
end

%%


            
for iter0 = [1000]%[0,10,100,1000]
%     for p = [1.25,2,5,10]
    for p = [1.25,2,5,10]
        for rho =  [ 1e-26,1e-13,1e-10,1e-5,0.0001]
            [iter0,p,rho]
            if rho < 1e-5
            filename = sprintf('saved_data/PCGM_rho%e_p%f_iter0%d.mat',rho,p,iter0)
            else
            filename = sprintf('saved_data/PCGM_rho%f_p%f_iter0%d.mat',rho,p,iter0);
            end
            try
                data = load(filename);
                gap_track = data.gap_track;
                obj_track = data.obj_track;
                train_track = data.train_track;
                test_track = data.test_track;
            catch
                [theta,obj_track, train_track, test_track, gap_track] = PCGM(n,maxiter, rho, get_trainobj, smoothgrad,  '1norm',p,iter0);
                save(filename,'theta','obj_track', 'gap_track','train_track','test_track')
            end
            
        end
    end
end

%%
for iter0 = [0,10,1000]
    for th = [1/8,1/4,1/2,3/4]
        for rho =  [1e-26,1e-14,1e-13,1e-12,1e-11,1e-5,1e-4,1e-3]
            for p = [1.25,2,5]
                [iter0,th,rho,p]
                cbar = 1-th;
                if rho < 1e-5
                    filename = sprintf('saved_data/RPCGM_rho%e_p%f_cbar%f_theta%f_iter0%d.mat',rho,p,cbar,th,iter0);
                else
                    filename = sprintf('saved_data/RPCGM_rho%f_p%f_cbar%f_theta%f_iter0%d.mat',rho,p,cbar,th,iter0);
                end
                filename
                try
                    data = load(filename);
                    gap_track = data.gap_track;
                    obj_track = data.obj_track;
                    train_track = data.train_track;
                    test_track = data.test_track;
                catch
                    [theta,obj_track, train_track, test_track, gap_track] = RPCGM(n,maxiter, rho, th, cbar, get_trainobj, smoothgrad,  '1norm',p,iter0);
                    
                    save(filename,'theta','obj_track', 'gap_track','train_track','test_track')
                    
                end
                
            end
        end
    end
end









