clc;
clear;
close all;
funs.layout();

%% 1. MEX

% threads = 56;
% model.mex_solve(threads,1);
% model.mex_simulate(threads,1);
% clc;

%% 2. specifications

% specs =
% { 
% 1. model,
% 2. name,
% 3. [epstein_zin,target_IRF,do_figures],
% 4. calibration parameters,
% 5. calibration values (including starting values)
% 6. estimation parameters
% 7. calibration parameters to loop over
% 8. matrix of calibration values
% }

DO_MAIN = 1;
DO_PREFS = 1;
DO_CRRA = 1;
DO_ALPHAS = 1;
DO_HET = 1;
DO_KV = 1;
DO_ROBUST = 1;

specs = {};

if DO_MAIN
    specs{end+1} = {'bs','ez_LCP_rho',[1,0,1],{},[],{'beta','zeta'},{'rho'},[1.5;2.0;4.0;6.0]};
    specs{end+1} = {'bs','ez_LCP_IRF',[1,1,1],{'beta','rho','zeta'},[0.94,6,1.6],{'beta','zeta','rho'},{},[]};
    specs{end+1} = {'bs','ez_LCP_IRF_per',[1,1,1],{'beta','rho','zeta','alpha_tilde'},[0.94,4,1.6,1.2],{'beta','zeta','alpha_tilde'},{},[]};
    specs{end+1} = {'bs','ez_LCP_IRF_var',[1,1,1],{'beta','rho','zeta','alpha'},[0.94,4,1.6,1.2],{'beta','zeta','alpha'},{},[]};
end
if DO_PREFS
    sigmas = [0.5 0.8 1.5];
    rhos = [1.5 2 4 6 8 10];
    [sigmas,rhos] = meshgrid(sigmas,rhos);
    specs{end+1} = {'bs','ez_LCP_sigma_rho',[1,0,0],{},[],{'beta','zeta'},{'sigma','rho'},[sigmas(:) rhos(:)]};
end
if DO_CRRA
    specs{end+1} = {'bs','CRRA_LCP_rho',[0,0,0],{},[],{'beta','zeta'},{'rho'},[1.5;2.0;4.0;6.0]};
    specs{end+1} = {'bs','CRRA_LCP_IRF_per',[0,1,0],{'beta','rho','zeta','alpha_tilde'},[0.80,4,5.0,1.4],{'beta','zeta','alpha_tilde'},{},[]};
    specs{end+1} = {'bs','CRRA_LCP_IRF_var',[0,1,0],{'beta','rho','zeta','alpha'},[0.80,4,5.0,1.4],{'beta','zeta','alpha'},{},[]};
end
if DO_ALPHAS
    specs{end+1} = {'bs','ez_LCP_IRF_per_rho',[1,1,0],{'alpha_tilde'},[1.4],{'beta','zeta','alpha_tilde'},{'rho'},[1.5;2.0;3.0]};
    specs{end+1} = {'bs','ez_LCP_IRF_var_rho',[1,1,0],{'alpha'},[1.4],{'beta','zeta','alpha'},{'rho'},[1.5;2.0;3.0]};
end
if DO_HET
    specs{end+1} = {'bs','het_LCP_rho',[1,0,1],{},[],{'het_est','delta_est','zeta'},{'rho'},[1.5;2.0;4.0;6.0]};                
    specs{end+1} = {'bs','het_LCP_IRF',[1,1,1],{'het_est','delta_est','rho','zeta'},[2.5,-1.0,6,1.6],{'het_est','delta_est','zeta','rho'},{},[]};
end
if DO_KV
    specs{end+1} = {'kv','ez_LCP_rho',[1,0,1],{},[],{'beta','zeta'},{'rho'},[1.5;2.0;4.0;6.0]};
    specs{end+1} = {'kv','ez_LCP_IRF',[1,1,1],{'beta','rho','zeta'},[0.90,8,1.5],{'beta','zeta','rho'},{},[]};   
end
if DO_ROBUST
    specs{end+1} = {'bs','robust_kappa',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'kappa'},[0.80;1.00],'\kappa'};
    specs{end+1} = {'bs','robust_sigma_psi',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'sigma_psi'},[0.08;0.10],'\sigma_\psi'};
    specs{end+1} = {'bs','robust_sigma_xi',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'sigma_xi'},[0.05;0.07],'\sigma_\xi'};   
    specs{end+1} = {'bs','robust_omega_work',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'omega_work'},[0.15;0.35],'\omega'};     
    specs{end+1} = {'bs','robust_R',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'R'},[1.01;1.03],'R'}; 
    specs{end+1} = {'bs','robust_Rneq',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'Rneq'},[1.06;1.10],'R_{-}'}; 
    specs{end+1} = {'bs','robust_death_mean',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'death_mean'},[70;85],'\mu_H'}; 
    specs{end+1} = {'bs','robust_death_std',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'death_std'},[6;12],'\sigma_H'}; 
    specs{end+1} = {'bs','robust_inh_A_raw',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'inh_A_raw'},[0.7;1.1],'h_{45}'}; 
    specs{end+1} = {'bs','robust_inh_R',[1,1,0],{'rho'},[6.0],{'beta','zeta','rho'},{'inh_R'},[0.99;1.01],'\eta'};
end

% robustness table
file = fopen(sprintf('%s/figs_tabs/table_lines/robustness.txt',pwd),'w');
for i = 1:numel(specs)
    if contains(specs{i}{2},'robust')
    for j = 1:size(specs{i}{8},1)
        name = sprintf('%s_%s_%d',specs{i}{1},specs{i}{2},j);
        val = specs{i}{8}(j);
        if ~any(specs{i}{8} ~= round(specs{i}{8}))
            fprintf(file,'$%s = %d $ & \\input{figs_tabs/table_lines/%s_noalphas.txt} \\\\ \n',specs{i}{9},val,name);
        else
            fprintf(file,'$%s = %3.2f$ & \\input{figs_tabs/table_lines/%s_noalphas.txt} \\\\ \n',specs{i}{9},val,name);    
        end        
    end    
    end
end
fclose(file);
        
%% 3. run

EST = 1;
DO_POST = 1;
DO_INFO = 1;
DO_DECOMPO = 1;

for i = 1:numel(specs)
	
    % 1. unpack
    spec = specs{i};
    modeltype = spec{1};
    name = spec{2};
    dos = spec{3};
    
        epstein_zin = dos(1);
        target_IRF = dos(2);
        do_post_save_to_disc = dos(3); 
        
    parnames = spec{4};
    parvals = spec{5};
    
    estvars = spec{6};
    
    parrnames = spec{7};
    parrvals = spec{8};
    
    if size(parrvals,1) > 1
        low = 1;
        high = size(parrvals,1);
    else
        low = 1;
        high = 1;
    end
    
for j = low:high
            
    % 2. filenames and print        
    if strcmp(name,'')
        model = modeltype;
    else
        model = sprintf('%s_%s',modeltype,name);
    end
    if size(parrvals,1) > 0
        model = sprintf('%s_%d',model,j);
    end
    matfile = sprintf('estimates/%s',model);  
    
        fprintf('\n\n');
        fprintf('********************************\n');
        fprintf('** %-26s **\n',model);
        fprintf('********************************\n\n');
        
    % 3. main
    if EST == 1 || (EST == 2 && exist(sprintf('%s.mat',matfile),'file') == 0)
        
        % a. baseline setup
        if strcmp(modeltype,'bs')
            par = setup.baseline_bs();   
            par.moms.vars = {{'ab',50}};                 
        elseif strcmp(modeltype,'kv')
            par = setup.baseline();
            par.moms.vars = {{'ab',50}};               
        else
            error('unknown setup');
        end
        par.prefix = model;
        
        % b. update
        [par] = estimate.updatepar(par,parnames,parvals);
            
            % i. do's
            par.epstein_zin = epstein_zin;
            par.target_IRF = target_IRF;
            par.do_post_save_to_disc = do_post_save_to_disc;

            % ii. heterogeneity
            if contains(name,'het')
                par.het = 5;
                par.hetvar = 'beta';    
                par.het_min = 0;
                par.het_max = 1;             
                par.moms.vars = {{'ab',50},{'ab','iq'}};
            end
            
            % iii. range
            for k = 1:size(parrvals,2)
                par.(parrnames{k}) = parrvals(j,k);
            end
        
        % c. estimate
        if par.het == 0 && size(parrvals,1) > 1
            [par.beta,par.zeta] = estimate.starting_values(par);
        elseif size(parrvals,1) > 1
            [par.het_est,par.delta_est,par.zeta] = estimate.starting_values_het(par);            
        end
        par.estvars = estvars;
        par = estimate.run(par,estvars,1);
        
        % d. save
        save(matfile,'par');
        
    else
        
        load(matfile,'par');
        par = setup.loaddata(par);
        
    end        
    par.estvars = estvars;
            
    % 4. postestimation
    if DO_POST && EST ~= 1
        par = estimate.post_estimation(par);
        save(matfile,'par');
    end      
    if DO_INFO && DO_POST ~= 1 && EST ~=1
        estimate.info(par);
    end
    
 	% 5. decomposition
    if DO_DECOMPO && par.adjcost_share == 0.0

        par_copy = par;
        
        % a. no retirement
        par = par_copy;
        par.zeta = 0;
        par.prefix = sprintf('%s_nr',par.prefix);
    
        fprintf('\n');
        fprintf(' ********************************\n');
        fprintf(' ** %-26s **\n',par.prefix);
        fprintf(' ********************************\n\n');

        par.do_post_figs = 0; 
        par.do_post_estimation_do_forget_inh = 0; 
        par = estimate.post_estimation(par);
        rmdir(sprintf('figs_tabs/%s',par.prefix));
        par_copy.mean_ab_nr = par.mean_ab;
        
        % b. no variance
        par = par_copy;
        
        par.sigma_psi = 0;
        par.sigma_xi = 0;        
        par.prefix = sprintf('%s_nv',par.prefix);
    
        fprintf('\n');
        fprintf('********************************\n');
        fprintf('** %-26s **\n',par.prefix);
        fprintf('********************************\n\n');

        par.do_post_figs = 0; 
        par.do_post_estimation_do_forget_inh = 0;  
        par = estimate.post_estimation(par);
        rmdir(sprintf('figs_tabs/%s',par.prefix));
        par_copy.mean_ab_nv = par.mean_ab;
        
        % c. no retirement, no variance 
        par = par_copy;
        
        par.sigma_psi = 0;
        par.sigma_xi = 0;      
        par.zeta = 0;
        par.prefix = sprintf('%s_nvr',par.prefix);
    
        fprintf('\n');
        fprintf('********************************\n');
        fprintf('** %-26s **\n',par.prefix);
        fprintf('********************************\n\n');

        par.do_post_figs = 0; 
        par.do_post_estimation_do_forget_inh = 0;   
        par = estimate.post_estimation(par);
        rmdir(sprintf('figs_tabs/%s',par.prefix));
        par_copy.mean_ab_nvr = par.mean_ab;
        
        % d. save
        par = par_copy;        
        save(matfile,'par');
        
    end
    
end
end