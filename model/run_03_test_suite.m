clc;
clear;
close all;
funs.layout();

%% MEX

LOAD = 0;

do_bs_PIH = 1;
do_bs_vfi = 1;
do_bs_ez = 1;
do_bs_zeta = 1;
do_bs_grid = 1;

do_grid = 1;
do_lambda = 1;
do_singles = 1;
do_vfi = 1;

%% 1. PIH with beta*R = flat consumption

if do_bs_PIH

    % a. standard
    setup_func = 'test_bs_PIH';
    parnames = {};
    parvals = [];
    prefix = 'test_bs_PIH';
    [par,sol,sim] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    % b. Epstein-Zin - high rho
    setup_func = 'test_bs_PIH';
    parnames = {'epstein_zin','rho'};
    parvals = [1,4];
    prefix = 'test_bs_PIH_ez_high_rho';
    [par_ez_highrho,sol_ez_highrho,sim_ez_highrho] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    % c. zeta = 1
    setup_func = 'test_bs_PIH';
    parnames = {'zeta'};
    parvals = [1];
    prefix = 'test_bs_PIH_high_zeta';
    [par_highzeta,sol_highzeta,sim_highzeta] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    % d. H > 0
    setup_func = 'test_bs_PIH';
    parnames = {'Nd','inh_A_raw'};
    parvals = [2,1];
    prefix = 'test_bs_PIH_high_zeta_inh';
    [par_inh,sol_inh,sim_inh] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    % e. figure
    pars = {par,par_ez_highrho,par_highzeta,par_inh};
    sims = {sim,sim_ez_highrho,sim_highzeta,sim_inh};
    names = {'CRRA','Epstein-Zin ($\rho = 4$)','CRRA ($\zeta = 1$)','CRRA ($h_{45} = 1, \eta = R$)'};
    figs.compare_lcp('compare','C','mean','$C_t$',pars,sims,names);
    figs.compare_lcp('compare','AB','mean','$A_t+B_t$',pars,sims,names);

end
clearvars -except do_* LOAD;

%% 2. Buffer-stock: VFI

if do_bs_vfi
    
    pars = cell(1,1);
    sims = cell(1,1);
    names = cell(1,1);
    
    zeta = 1.1;
        
    % a. standard
    j = 1;
    setup_func = 'baseline_bs';
    parnames = {'zeta'};
    parvals = [zeta];
    prefix = 'test_bs_vfi_negm';
    [pars{j},~,sims{j}] = model.all(LOAD,setup_func,parnames,parvals,prefix);     
    names{j} = 'baseline';    
    j = j+1;
        
    % b. vfi
    setup_func = 'baseline_bs';
    parnames = {'zeta','vfi_pure'};
    parvals = [zeta,1];
    prefix = 'test_bs_vfi_vfi';
    [pars{j},sols{j},sims{j}] = model.all(LOAD,setup_func,parnames,parvals,prefix);
    names{j} = 'VFI';  
    j = j+1;
    
    % e. figure
    figs.compare_lcp('compare','C','mean','$C_t$',pars,sims,names);
    figs.compare_lcp('compare','AB','mean','$A_t+B_t$',pars,sims,names);
    
end
clearvars -except do_* LOAD;

%% 3. Buffer-stock: Epstein-Zin -> CRRA 

if do_bs_ez
    
    Nsigma = 11;
    sigmas = linspace(1.5,2.5,Nsigma);
    
    pars = cell(Nsigma,1);
    sims = cell(Nsigma,1);

    % a. standard
    setup_func = 'baseline_bs';
    parnames = {'epstein_zin','rho','sigma'};
    parvals = [0,2,2];
    prefix = 'test_bs_ez';
    [par_base,~,sim_base] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    % b. Epstein-Zin
    for i = 1:Nsigma
        setup_func = 'baseline_bs';
        parnames = {'epstein_zin','rho','sigma'};
        parvals = [1,2,sigmas(i)];
        prefix = sprintf('test_bs_ez_%d',i);
        [pars{i},~,sims{i}] = model.all(LOAD,setup_func,parnames,parvals,prefix);
    end

    % c. figure
    name = 'CRRA ($\rho = \sigma = 2$)';
    name_alt = 'Epstein-Zin ($\rho = 2$)';
    figs.compare_AB(par_base.TR,'ez','$\sigma$','$A_t+B_t$',par_base,sim_base,sigmas,pars,sims,name,name_alt);
    
end
clearvars -except do_* LOAD;


%% 4. Buffer-stock: zeta -> 0

if do_bs_zeta
    
    Nzeta = 11;
    zetas = linspace(0.01,1.5,Nzeta);
    
    pars = cell(Nzeta,1);
    sims = cell(Nzeta,1);

    % a. standard
    setup_func = 'baseline_bs';
    parnames = {'zeta'};
    parvals = [0];
    prefix = 'test_bs_zeta';
    [par_base,~,sim_base] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    % b. various zeta
    for i = 1:Nzeta
        parnames = {'zeta'};
        parvals = [zetas(i)];
        prefix = sprintf('test_bs_zeta_%d',i);
        [pars{i},~,sims{i}] = model.all(LOAD,setup_func,parnames,parvals,prefix);
    end

    % c. figure
    name = '$\zeta = 0$';
    name_alt = 'different $\zeta$';
    for t = [par_base.TR]
        figs.compare_AB(t,'zetas','$\zeta$','$A_t+B_t$',par_base,sim_base,zetas,pars,sims,name,name_alt);
    end
    
end
clearvars -except do_* LOAD;

%% 5. Buffer-stock: grids

if do_bs_grid
    
    % a. standard
    setup_func = 'baseline_bs';
    prefix = 'test_bs_grid';
    [par_base,~,sim_base] = model.all(LOAD,setup_func,{},[],prefix);
    
    pars = cell(1,1);
    simss = cell(1,1);
    
    % b. setup
    vec = [-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3];
    
    basenames = {'NA','NM','NP'};
    basevals = [par_base.NA, par_base.NM, par_base.NP];
    stepsize = [par_base.NA_step, par_base.NM_step, par_base.NP_step];
    
    % c. run
    for i = 1:numel(vec)

        prefix = sprintf('test_bs_grid_%d',i);
        
        parnames = basenames;      
        parvals = basevals + vec(i)*stepsize;       
        
        [pars{i},~,sims{i}] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    end
    
    % d. figure
    name = 'baseline';
    name_alt = 'different grids';
    for t = [20 35]
       figs.compare_AB(t,'grid','$j$','$A_t+B_t$',par_base,sim_base,vec,pars,sims,name,name_alt);
    end
    
end
clearvars -except do_* LOAD;

%% 6. Two-asset: grids

if do_grid
        
    % a. standard
    setup_func = 'baseline';
    prefix = 'test_grid';
    [par_base,sol_base,sim_base] = model.all(LOAD,setup_func,{},[],prefix);
    
    pars = cell(1,1);
    sims = cell(1,1);
    
    % b. setup   
    vec = [-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
    
    basenames = {'NA','NM','NN','NP','NX'};
    basevals = [par_base.NA, par_base.NM, par_base.NN, par_base.NP, par_base.NX];
    stepsize = [par_base.NA_step, par_base.NM_step, par_base.NN_step, par_base.NP_step, par_base.NX_step];
    
    % c. run
    for i = 1:numel(vec)

        prefix = sprintf('test_grid_%d',i);
        
        parnames = basenames;   
        parvals = basevals + vec(i)*stepsize;  
        
        [pars{i},~,sims{i}] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    end
    
    % d. figure
    name = 'baseline';
    name_alt = 'different grids';
    for t = [20 par_base.TR]
       figs.compare_AB(t,'grid','$j$','$A_t+B_t$',par_base,sim_base,vec,pars,sims,name,name_alt);
    end
    
end
clearvars -except do_* LOAD;

%% 7. Two-asset: lambda -> 0

if do_lambda
    
    % a. standard    
    Nlambdas = 7;
    lambdas = linspace(0.0,0.02,Nlambdas);
    
    pars = cell(Nlambdas,1);
    sims = cell(Nlambdas,1);

    setup_func = 'baseline_bs';
    parnames = {'beta','R','Rb'};
    parvals = [0.935,1.057,1.00];
    prefix = 'test_lambda';
    [par_base,~,sim_base] = model.all(LOAD,setup_func,parnames,parvals,prefix);

    % b. various lambda
    for i = 1:Nlambdas
        setup_func = 'baseline';
        parnames = {'adjcost_share'};
        parvals = [lambdas(i)];
        prefix = sprintf('test_lambda_%d',i);
        [pars{i},~,sims{i}] = model.all(LOAD,setup_func,parnames,parvals,prefix);
    end

    % c. figure
    name = 'buffer-stock model ($R = 1.057, R_B = 1.00$)';
    name_alt = 'different $\lambda$';
    for t = [20 par_base.TR]
        figs.compare_AB(t,'lambas','$\lambda$','$A_t+B_t$',par_base,sim_base,lambdas,pars,sims,name,name_alt);
    end
    
end
clearvars -except do_* LOAD;

%% 8. Two-asset: various

if do_singles
    
    % a. setup
    N = 5;
    
    parname_vec = {'beta','rho','zeta','sigma','kappa','sigma_psi','sigma_xi'};
    latex_vec = {'$\beta$','$\rho$','$\zeta$','$\sigma$','$\kappa$','$\sigma_{\psi}$','$\sigma_{\xi}$'};
     
    vals = struct();
    vals.beta      = linspace(0.92,0.94,N);
    vals.rho       = linspace(1.5,4,N);    
    vals.zeta      = linspace(0,1.5,N);    
    vals.sigma     = linspace(0.4,2.0,N);
    vals.kappa     = linspace(0.5,1,N);
    vals.sigma_psi = linspace(0.01,0.14,N);
    vals.sigma_xi  = linspace(0.01,0.14,N);
        
    for j = 1:numel(parname_vec)
    
        pars = cell(N,1);
        sims = cell(N,1);
        
        parname_now = parname_vec{j};
        latex_now = latex_vec{j};        
        vals_now = vals.(parname_vec{j});

        % b. various
        for i = 1:numel(vals_now)
            setup_func = 'baseline';
            parnames = {parname_now};
            parvals = [vals_now(i)];
            prefix = sprintf('test_singles_%s_%d',parname_now,i);
            [pars{i},~,sims{i}] = model.all(LOAD,setup_func,parnames,parvals,prefix);
        end

        % c. figure
        par_base = pars{1};
        par_base.prefix = sprintf('test_singles_%s',parname_now);
        name = '';
        name_alt = '';
        for t = [20 par_base.TR]
            figs.compare_AB(t,parname_now,latex_now,'$A_t+B_t$',par_base,[],vals_now,pars,sims,name,name_alt);
        end
        
    end
    
end
clearvars -except do_* LOAD;

%% 9. Two-asset: VFI

if do_vfi
            
    pars = cell(1,1);
    sims = cell(1,1);
    names = cell(1,1);
            
    % a. standard
    j = 1;
    setup_func = 'baseline';
    parnames = {};
    parvals = [];
    prefix = 'test_vfi_negm';
    [pars{j},~,sims{j}] = model.all(LOAD,setup_func,parnames,parvals,prefix);     
    names{j} = 'baseline';    
    j = j+1;
        
    % b. vfi
    setup_func = 'baseline';
    parnames = {'vfi_pure','Nm','Nn'};
    parvals = [1,100,80];
    prefix = 'test_vfi_vfi';
    [pars{j},sols{j},sims{j}] = model.all(LOAD,setup_func,parnames,parvals,prefix);
    names{j} = 'VFI';  
    j = j + 1;
    
    % c. figure
    figs.compare_lcp('compare','C','mean','$C_t$',pars,sims,names);
    figs.compare_lcp('compare','A','mean','$A_t$',pars,sims,names);
    figs.compare_lcp('compare','B','mean','$B_t$',pars,sims,names);
    figs.compare_lcp('compare','AB','mean','$A_t+B_t$',pars,sims,names);
    
end
clearvars -except do_* LOAD;

%% clean up

rmdir('figs_tabs/test_bs*','s');
rmdir('figs_tabs/test_grid*','s');
rmdir('figs_tabs/test_lambda*','s');
rmdir('figs_tabs/test_singles*','s');
rmdir('figs_tabs/test_vfi*','s');
delete('csv/test_suite*');
fprintf('Done\n');