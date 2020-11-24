classdef setup
methods(Static)
    
function par = baseline()
            
    par = struct();
       
    % a. discrete states and choices   
    par.Nd = 2; 
    par.Nz = 2;
    par.Nk = 3;
    
    % b. grids
    par.NM = 300;
    par.NX = 200;
    par.NP = 150;
    par.NN = 100;
    par.NA = 100;
    
        % step
        par.NM_step = 80;
        par.NX_step = 60; 
        par.NP_step = 40;
        par.NN_step = 30;
        par.NA_step = 30;

        % min/max
        par.Pmin = 2e-1;
        par.Pmax = 6.0;
        par.mmax = 10.0;
        par.xmax = 10.0;
        par.nmax = 10.0;
        par.amax = 10.0;
                
        % curvature
        par.phi_P = 1.25;
        par.phi_m = 1.25;
        par.phi_n = 1.25;
        par.phi_x = 1.25;           
        par.phi_a = 1.25;
    
    % c. demographics, t in (1,2,...,T)
    par.T  = 60; % die at end-of-period
    par.TR = 35; % retire at end-of-period
    par.age_min = 24;

    % d. preferences
    par.rho  = 2.000;
    par.beta = 0.935;
    par.zeta = 1.000;
            
    par.epstein_zin = 1;
    par.sigma       = 2/3;

    % e. income process
    par.kappa = 0.90;
    
    % f. shocks    
    par.Npsi = 6;
    par.Nxi  = 6;
    
        par.alpha = 1.00;
        par.alpha_tilde = 1.00;
        
    % g. assets
    par.R = 1.02;
    par.Rb = 1.057;
    par.Rneg = 1.078;
   
    par.omega_work = 0.25;
    
    par.adjcost_share = 0.02;
    par.adjcostinf    = 0.00; % prob. of infinity adjustment costs
    
    % h. inheritance
    par.age_diff   = 30;
    par.death_mean = 77;
    par.death_std  = 9;
    
    par.inh_A_raw = 0.93;
    par.inh_R     = 1.00;
    par.inh_base_age = 20;
           
    % i. simulate
    
        % i. time and number of households
        par.T0sim      = 0;
        par.Tendsim    = par.TR;
        par.Nsim       = 50000;
        par.post_Nsim  = 100000;
        
        % ii. settings
        par.forget_inh    = 0;   % run simulation without inheritance expectations
        par.do_forget_inh = nan;
        par.save_to_disc  = 0;   % save simulation results in log_sim.txt    
        
        % iii. initial distribution
        par.P_lag_ini_sigma = 0.40;
        par.N_ini_sigma     = 0.80;
        par.N_ini_pos       = 0.30;
        par.N_add           = 0.0;
        
        % iv. MPC
        par.MPC_add_share = 0.01;
        par.cutoff_share  = 0.02;
        
        % v. IRF
        par.sim_IRF = 0;
        par.P_IRF = [];
        par.A_IRF = [];
        par.B_IRF = [];
    
    % j. technical
    par.tmin      = 0; % solve until tmin
    par.vfi_pure  = 0; % 0 = NGEM, 1 = VFI
    
    par.no_multistart = 0; % multistart or not in pure VFI
    par.no_interp_vec = 0; % vectorized interpolation i NEGM
    
    % k. estimation    
    par.moms.age  = 25:59;  
    par.moms.vars = {};
    par.target_IRF = 0;
    par.tol_x = 1e-3;
    par.tol_fun = 1e-4;
    par.maxiter = 200;
    par.iter_rough = 50;
    
        % heterogeneity
        par.het = 0;
        par.hetvar    = '';    
        par.het_est   = nan;
        par.het_min   = nan;
        par.het_max   = nan;    
        par.delta_est = nan;
        
        % post-estimation
        par.do_post_figs = 1;
        par.do_post_estimation_do_forget_inh = 1;
        par.do_post_save_to_disc = 0;
        
    par = setup.loaddata(par);
    
end

function par = loaddata(par)
    
    % a. from income variances
    tbl = readtable('data\income_stds.csv');
    par.sigma_xi = tbl.sigma_xi;
    par.sigma_psi = tbl.sigma_psi;
    
    % b. from all.csv
    tbl = readtable('data\all.csv');
    
    par.P0           = tbl.x_perminc_mean(1);
    par.moms.P_mean  = tbl.x_perminc_mean/par.P0;
    par.moms.P_p25   = tbl.x_perminc_p25/par.P0;
    par.moms.P_p50   = tbl.x_perminc_p50/par.P0;
    par.moms.P_p75   = tbl.x_perminc_p75/par.P0;    

    par.moms.a_mean = tbl.x_liqworth_mean;     
    par.moms.a_p25  = tbl.x_liqworth_p25;     
    par.moms.a_p50  = tbl.x_liqworth_p50; 
    par.moms.a_p75  = tbl.x_liqworth_p75; 
    
    par.moms.b_mean = tbl.x_illiqworth_mean;     
    par.moms.b_p25  = tbl.x_illiqworth_p25;     
    par.moms.b_p50  = tbl.x_illiqworth_p50; 
    par.moms.b_p75  = tbl.x_illiqworth_p75; 
    
    par.moms.ab_mean = tbl.x_networth_mean;     
    par.moms.ab_p25  = tbl.x_networth_p25;     
    par.moms.ab_p50  = tbl.x_networth_p50; 
    par.moms.ab_p75  = tbl.x_networth_p75; 
        
    par.moms.ab_iq  = par.moms.ab_p75-par.moms.ab_p25;
    
    par.moms.b_pos_p25  = tbl.x_illiqworth_cond_p25;
    par.moms.b_pos_p25  = tbl.x_illiqworth_cond_p25;    
    par.moms.b_pos_p75  = tbl.x_illiqworth_cond_p75;    
    
    par.moms.B_pos_mean  = tbl.x_pos_illiqworth;      
    par.moms.h_cond_mean = tbl.x_inheritance_mean;
            
    % c. from ageatdeath.csv    
    tbl = readtable('data/ageatdeath.csv');
    par.moms.H_pos_mean = [nan; tbl.count];
        
    % d. load IRF data
    tbl = readtable('data\empirical_networth.csv');
    par.IRF.ab.beta_emp = tbl.beta_networth(6:end);
    par.IRF.ab.w_emp = sum(tbl.se_networth(6:end))./tbl.se_networth(6:end);
    
    tbl = readtable('data\empirical_liqworth.csv');
    par.IRF.a.beta_emp = tbl.beta_liqworth(6:end);
    par.IRF.a.w_emp = sum(tbl.se_liqworth(6:end))./tbl.se_liqworth(6:end);
    
end
function par = adj_grid(par,fac)
    
    par.NM = par.NM + fac*par.NM_step;
    par.NP = par.NP + fac*par.NP_step;
    par.NX = par.NX + fac*par.NX_step;    
    par.NA = par.NA + fac*par.NA_step;      
    if par.NN > 1
        par.NN = par.NN + fac*par.NN_step;
    end

end

function par = baseline_bs()
    
    % a. basis
    par = setup.baseline();
    
    % b. changes      
        
        % i. grids
        par.Nz = 1;       
        
        par.NM = 600;
        par.NP = 150;
        par.NX = 1;
        par.NN = 1;
        par.NA = 150;

            % step
            par.NM_step = 100;
            par.NP_step = 40;
            par.NX_step = 0; 
            par.NN_step = 0;
            par.NA_step = 40;
                
            % cuv
            par.phi_P = 1.2;
            par.phi_m = 1.2;
            par.phi_n = 1.2;
            par.phi_x = 1.2;           
            par.phi_a = 1.2;

        % ii. preferences
        par.beta = 0.97;
        
        % iii. assets
        par.R  = 1.02;
        par.Rb = 1.00;
        par.adjcost_share = 0.0;
        
end
function par = test_bs_PIH()
    
    % a. basis
    par = setup.baseline_bs();
    
    % b. grids
    par.Nd = 1;

    % c. preferences
    par.epstein_zin = 0;
    par.sigma       = 2;
    par.rho         = 2;
    par.beta        = 0.97; 
    par.zeta        = 0.0;

    % d. assets
    par.R              = 1.0/par.beta;
    par.Rneg           = par.R;
    par.Rb             = 1.00;
    par.omega_work     = 2.0;

    % e. inheritance
    par.inh_A_raw = 0.0;
    par.inh_R     = par.R;
        
    % e. income shocks
    par.sigma_psi = 0.0;
    par.sigma_xi = 0.0;
    par.Npsi = 1;
    par.Nxi  = 1;

    par.P_lag_ini_sigma = 0.0;
    par.N_ini_sigma     = 0.0;
    par.N_ini_pos       = 0.0;        
    par.N_add           = 5;
        
end

end
end