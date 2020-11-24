classdef model
methods(Static)

%%%%%%%%%%
% 1. mex %
%%%%%%%%%%

function [] = mex_solve(threads,do_AVX515)
    
    str = sprintf('mex -largeArrayDims cfuncs/solve_model.cpp -DMAXTHREADS=%d',threads);
    if nargin == 1 || do_AVX515 == 0
        str = sprintf('%s %s',str,' COMPFLAGS=''$COMPFLAGS /o3 /openmp''');
    else
        str = sprintf('%s %s',str,' COMPFLAGS=''$COMPFLAGS /o3 /openmp /arch:CORE-AVX512  ''');
    end    
    str = sprintf('%s %s',str,' cfuncs\libnlopt-0.lib');
    
    eval(str);
    fprintf('\n');
    
end
function [] = mex_simulate(threads,do_AVX515)
    
    str = sprintf('mex -largeArrayDims cfuncs/simulate_model.cpp -DMAXTHREADS=%d',threads);
    if nargin == 1 || do_AVX515 == 0
        str = sprintf('%s %s',str,' COMPFLAGS=''$COMPFLAGS /o3 /openmp''');
    else
        str = sprintf('%s %s',str,' COMPFLAGS=''$COMPFLAGS /o3 /openmp /arch:CORE-AVX512  ''');
    end 
    eval(str);
    fprintf('\n');
        
end


%%%%%%%%%%%%
% 2. solve %
%%%%%%%%%%%%

function par = inheritance(par)
        
    % 1. full matrix of conditional probabilities
    par.inh_p_full = zeros(par.T,par.T);
    for t = 1:par.T
        
        for tplus = t:par.T
            
            parent_age = par.age_min + par.age_diff + tplus;
            
            die_now = normcdf(parent_age, par.death_mean, par.death_std);
            die_next = normcdf(parent_age+1, par.death_mean, par.death_std);
            
            % probability of dying within next year
            par.inh_p_full(t,tplus) = (die_next-die_now);
        
        end
        
        % scale probabilities to sum to one
        futsum = sum(par.inh_p_full(t,:));
        if futsum <= 1e-8
            par.inh_p_full(t,:) = 0.0;
        else
            par.inh_p_full(t,:) =  par.inh_p_full(t,:)/futsum;
        end
    end
    
    % 2. probability at each age
    par.inh_p = zeros(par.T,1);
    for t = 1:par.T-1
        if par.inh_p_full(t+1,t+1) == 0.0
            par.inh_p(t) = 1.0;
        else
            par.inh_p(t) = par.inh_p_full(t,t);
        end
    end
    par.inh_p(end) = 1.0;
    
end
function par = create_grids(par)
    
    assert((par.NN == 1 && par.Nz == 1) || (par.NN > 1 && par.Nz == 2));
    assert(par.rho > 0);
    assert(abs(par.sigma-1) > 1e-8);
    assert(par.alpha == 1 || par.alpha_tilde == 1);
    
        % if not epstein-zin then sigma=rho
        if par.epstein_zin == 0
            par.sigma = par.rho;
        end
    
        % clean up
        folder = ['figs_tabs\' par.prefix];
        if exist(folder,'dir') > 0
            delete(sprintf('%s/*.png',folder));       
        else
            mkdir(folder);    
        end
        
    % 1. income growth rate 
    par.Gamma = ones(par.TR,1);  
    par.G = ones(par.TR,1);    
    if isfield(par,'G_raw') == 1
        par.G(1:par.TR-1) = par.G_raw;
    else
        P =  par.moms.P_mean;  
        par.G(1:par.TR-1) = P(2:end)./P(1:end-1);           
    end
    par.G(par.TR) = 1.0;        
    
    for t = 2:par.TR
        par.Gamma(t) = prod(par.G(1:t-1));
    end
    
    % 2. adjcost and MPC
    par.avg_P   = mean(par.Gamma(1:par.TR));
    par.adjcost = par.adjcost_share*par.avg_P;
    par.MPC_add = par.MPC_add_share*par.avg_P;
    par.cutoff  = par.cutoff_share*par.avg_P;
    
    % 3. borrowing
    par.omega             = zeros(par.TR,1);
    par.omega(1:par.TR-1) = par.omega_work;
    
    % 4. grids
    par.grid_P = nan(par.NP,par.TR);        
    for t = 1:par.TR
        par.grid_P(:,t) = funs.nonlinspace(par.Pmin,par.Pmax,par.NP,par.phi_P)*par.Gamma(t);       
    end
    
    par.grid_M = nan(par.NM,par.NP,par.TR);
    par.grid_N = nan(par.NN,par.NP,par.TR);
    par.grid_X = nan(par.NX,par.NP,par.TR);
    par.grid_A = nan(par.NA,par.NP,par.TR);
    
    for t = 1:par.TR
        
        grid_m = funs.nonlinspace(-par.omega(t),par.mmax,par.NM,par.phi_m);
        grid_n = funs.nonlinspace(0,par.nmax,par.NN,par.phi_n);
        grid_x = funs.nonlinspace(0,par.xmax,par.NX,par.phi_x);        
        
        if par.vfi_pure == 0
            if par.omega(t) == 0
                grid_a = funs.nonlinspace(1e-8,par.amax,par.NA-1,par.phi_a);
            else
                NA_con = floor(par.NA*0.20);
                NA_rest = par.NA-1-NA_con;           
                grid_a_con = linspace(-par.omega(t)+1e-8,-1e-8,NA_con)';
                grid_a_rest = funs.nonlinspace(1e-8,par.amax,NA_rest,par.phi_a);
                grid_a = [grid_a_con; grid_a_rest];            
            end
        end
        
        for i_P = 1:par.NP               
            P = par.grid_P(i_P,t);
            par.grid_M(:,i_P,t) = grid_m*P;
            par.grid_N(:,i_P,t) = grid_n*P;   
            par.grid_X(:,i_P,t) = grid_x*P; 
            if par.vfi_pure == 0
                par.grid_A(:,i_P,t) = [-par.omega(t); grid_a]*P; 
            end
        end
        
    end
    
    % 5. shocks
    [xi, xi_w] = funs.GaussHermite_lognorm(par.sigma_xi,par.Nxi);    
    [psi, psi_w] = funs.GaussHermite_lognorm(par.sigma_psi,par.Npsi);
    
        % a. vectorized tensor products
        [psi, xi] = ndgrid(psi,xi);
        par.psi = psi(:);
        par.xi = xi(:);
    
        [psi_w, xi_w] = ndgrid(psi_w,xi_w);
        weights = psi_w.*xi_w;
        par.weights = weights(:);

        % b. maximum number of shocks
        par.Nshocks_max = numel(par.weights);
            
        % c. convert into matrices
        par.psi = [repmat(par.psi,[1 par.TR])];
        par.xi = [repmat(par.xi,[1 par.TR])];
        par.weights = [repmat(par.weights,[1 par.TR])];

        par.Nshocks =  int32([repmat(par.Nshocks_max,[1 par.TR])]);   
       
    % 6. inheritance
    par = model.inheritance(par);
    par.inh_fac = par.inh_R.^((1:par.T)-par.inh_base_age);
    par.inh_A   = par.inh_A_raw;
    
end
function [par,sol] = solve(par)
        
    sigma_psi = par.sigma_psi;
    sigma_xi = par.sigma_xi;       
    if par.alpha ~= 1
        par.sigma_psi = par.alpha*par.sigma_psi;     
        par.sigma_xi = par.alpha*par.sigma_xi;          
    elseif par.alpha_tilde ~= 1
        par.sigma_psi = par.alpha_tilde*par.sigma_psi;     
        par.sigma_xi = par.alpha_tilde*par.sigma_xi;                  
    end
        
    % 1. grids
    par = model.create_grids(par);

    % 2. solve

        t0 = tic;        
        
    sol = solve_model(par);        
        
        par.time_solve = toc(t0);
        
        % return shock stds to their original value
        par.sigma_psi = sigma_psi;     
        par.sigma_xi = sigma_xi; 
        
    % 3. saving
    sol.A_ast = cell(par.Nz,par.Nd,par.T);
    for t = 1:par.TR
    for d = 1:par.Nd
    for z = 1:par.Nz   

        if z == 2
            sol.A_ast{z,d,t} = nan(par.NX,par.NP);
            for i_P = 1:par.NP               
                sol.A_ast{z,d,t}(:,i_P) = par.grid_X(:,i_P,t)  - sol.C_ast{z,d,t}(:,i_P) - sol.B_ast{z,d,t}(:,i_P) - par.omega(t)*par.grid_P(i_P,t);
            end
        else
            sol.A_ast{z,d,t} = nan(par.NM,par.NN,par.NP);
            for i_P = 1:par.NP            
                sol.A_ast{z,d,t}(:,:,i_P) = par.grid_M(:,i_P,t) - sol.C_ast{z,d,t}(:,:,i_P);
            end
        end           

    end
    end
    end
    
end


%%%%%%%%%%%%%%%
% 3. simulate %
%%%%%%%%%%%%%%%

function par = setup_simulate(par)
    
    rng(2017)
    
        sigma_psi = par.sigma_psi;
        sigma_xi = par.sigma_xi;       
        if par.alpha ~= 1
            par.sigma_psi = par.alpha*par.sigma_psi;     
            par.sigma_xi = par.alpha*par.sigma_xi;                  
        end
    
    % 1. shocks - income
    par.psi_sim = ones(par.Nsim, par.TR);
    par.xi_sim  = ones(par.Nsim, par.TR);    
    par.xi_sim(:,1:par.TR) = exp(par.sigma_xi*randn(par.Nsim,par.TR) - 0.5*par.sigma_xi^2);        
    par.psi_sim(:,1:par.TR) = exp(par.sigma_psi*randn(par.Nsim,par.TR) - 0.5*par.sigma_psi^2);

    % 2. shocks -inheritance
    par.inh_pval = rand(par.Nsim,par.TR);
    
    % 3. shocks - A shocks
    par.A_shock = zeros(par.Nsim,par.TR);
    
    % 4. initial values
    par.P_ini_dist = exp(par.P_lag_ini_sigma*randn(par.Nsim,1) - 0.5*par.P_lag_ini_sigma^2);
    par.P_ini_dist = max(par.P_ini_dist,par.grid_P(1));
    par.N_ini_dist = exp(par.N_ini_sigma*randn(par.Nsim,1));
    par.N_ini_dist(rand(par.Nsim,1) > par.N_ini_pos) = 0; % N = 0
    par.N_ini_dist = par.N_ini_dist + par.N_add;
    
    % 5. fixed costs
    par.adjcost_pval = rand(par.Nsim,par.TR);
    
        % return shock stds to their original value
        par.sigma_psi = sigma_psi;     
        par.sigma_xi = sigma_xi;  
        
end
function par = select(par,ib,iu)
    
    % 1. shocks - income
    par.psi_sim = par.psi_sim(ib:iu,:);
    par.xi_sim = par.xi_sim(ib:iu,:);

    % 2. shocks -inheritance
    par.inh_pval = par.inh_pval(ib:iu,:);
    
    % 3. shocks - A shocks
    par.A_shock = par.A_shock(ib:iu,:);        
    
    % 4. initial values
    par.P_ini_dist = par.P_ini_dist(ib:iu);        
    par.N_ini_dist = par.N_ini_dist(ib:iu);
    
    % 5. fixed costs
    par.adjcost_pval = par.adjcost_pval(ib:iu,:);    
        
end
function [par,sim] = simulate(par,sol,do_setup)
           
    % 1. setup
    if do_setup == 1
        par = model.setup_simulate(par);
    end
    
    % 2. simulate
    
        t0 = tic;
            
    sim = simulate_model(par,sol);
    
        par.time_simulate = toc(t0);    
        
    sim = model.calc_additional_sim_var(par,sim);
    
    % 3. csv
    if par.save_to_disc > 0
        copyfile('sim_output.txt',sprintf('csv/%s.txt',par.prefix));
    end
    
    % 4. unexpected
    if par.do_forget_inh == 1
        par.forget_inh = 1;
        simulate_model(par,sol);          
        if par.save_to_disc > 0
            copyfile('sim_output.txt',sprintf('csv/%s_unexp.txt',par.prefix));
        end
        par.forget_inh = 0;
    end
    
    if par.save_to_disc > 0
        delete('sim_output.txt');
    end        
    
end
function sim = calc_additional_sim_var(par,sim)
   
    sim.a  = sim.A./(sim.P);
    sim.b  = sim.B./(sim.P);    
    sim.AB = sim.A + sim.B;
    sim.ab = sim.a + sim.b;
    sim.h  = sim.H./(sim.P);    

end


%%%%%%%%%%
% 4. ALL %
%%%%%%%%%%

function [par,sol,sim] = all(LOAD,setup_func,parnames,parvals,prefix)
    
    fprintf('model.all: %s',prefix);
    matfile = sprintf('csv/test_suite_%s.mat',prefix);
    if LOAD
        load(matfile);
        sol = struct();
        return;
    end
    
     % a. setup
    par = setup.(setup_func);
    par = estimate.updatepar(par,parnames,parvals);    
    par.prefix = prefix;    
    
    % c. solve   
    t0 = tic;
    [par, sol] = model.solve(par);
    fprintf(' (time = %4.2f secs)\n',toc(t0));
    % d. simulate
    if par.tmin == 0
        [par, sim] = model.simulate(par,sol,1);
    else
        sim = struct();
    end
    
    % c. save disc space
    par.psi_sim = [];
    par.xi_sim = [];
    par.inh_pval = [];
    par.A_shock = [];
    par.P_ini_dist = [];
    par.N_ini_dist = [];
    par.adjcost_pval = [];
    par.grid_P = [];
    par.grid_M = [];
    par.grid_N = [];
    par.grid_X = [];
    par.grid_A = [];
        
    if par.tmin == 0    
        sim_alt = struct();
        sim_alt.A = sim.A;
        sim_alt.B = sim.B;
        sim_alt.AB = sim.AB;
        sim_alt.C = sim.C;  
        sim_alt.P = sim.P;          
        sim = sim_alt;
    end
    
    % d. save
    save(matfile,'par','sim');

end


end
end