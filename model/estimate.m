classdef estimate
methods(Static)

%%%%%%%%%%%%%%%%
% 1. auxiliary %
%%%%%%%%%%%%%%%%

function dist = distance(par,sim,vars)
    
    dist = struct();
    dist.all = 0.0;
    for j = 1:numel(vars)            
                
        % 1. variable
        varnow = vars{j}{1};
        method = vars{j}{2};        
        if isnumeric(method)
            methodnow = sprintf('p%d',method);
        else
            methodnow = method;
        end
        dist.(varnow).(methodnow) = 0.0;
        
        % 2. age loop
        for age = par.moms.age     
            
            % a. time
            t = age-par.age_min;
                         
            % b. method
            if strcmp(method,'iq') == 1
                
                estmom = (prctile(sim.(varnow)(:,t),75)-prctile(sim.(varnow)(:,t),25));
                estmom_base = prctile(sim.(varnow)(:,1),75)-prctile(sim.(varnow)(:,1),25);
                estmom = estmom-estmom_base;

                datamom = par.moms.(sprintf('%s_iq',varnow))(t);
                datamom_base = par.moms.(sprintf('%s_iq',varnow))(1);
                datamom = datamom - datamom_base;

            elseif strcmp(method,'mean') == 1
            
                estmom = mean(sim.(varnow)(:,t));
                datamom = par.moms.(sprintf('%s_mean',varnow))(t);
            
            else
                
               estmom = prctile(sim.(varnow)(:,t),method);
               datamom = par.moms.(sprintf('%s_p%d',varnow,method))(t);
            
            end
            
            % c. sum
            w = 1/numel(par.moms.age);
            if strcmp(method,'iq') == 1
                w = 0.1*w;
            end
            dist.(varnow).(methodnow) = dist.(varnow).(methodnow) + w*(estmom - datamom)^2;
                
        end
        
        % 3. IRF
        if par.target_IRF
            dist.IRF = mean((par.IRF.ab.beta-par.IRF.ab.beta_emp).^2.*par.IRF.ab.w_emp);
            if isfield(par.IRF.a,'beta')
                dist.IRF_liq = mean((par.IRF.a.beta(2:end)-par.IRF.a.beta_emp(2:end)).^2.*par.IRF.a.w_emp(2:end));            
            end
        end
        
        % 4. total sum        
        dist.all = dist.all + dist.(varnow).(methodnow);
    
    end
    
    % 4. add IRF
    if par.target_IRF == 1
        dist.all = dist.all + dist.IRF;
    end
        
end
function par = het_vec(par)
     
    % 1. mean
    fac = exp(par.het_est)/(1+exp(par.het_est));
    par.(par.hetvar) = par.het_min + fac*(par.het_max - par.het_min);
    
    % 2. spread
    delta = exp(par.delta_est)/(1+exp(par.delta_est));                
    delta_max = delta*(par.het_max-par.(par.hetvar));
    delta_min = delta*(par.(par.hetvar)-par.het_min);
    par.delta = min(delta_max,delta_min);
    
    % 3. vector
    steps = linspace(-1,1,par.het);
    par.het_vec = par.(par.hetvar) + steps*par.delta;

end
function par = IRF_reg(par,sim,varname,print)
        
    % 1. age and time to inheritance
    age_vec = par.age_min+(1:par.TR);
    age = repmat(age_vec,[par.Nsim 1]);
    
    age_at_inheritance = (sim.h > 0)*age_vec';
    time_to_inheritance = age-age_at_inheritance;
    I = time_to_inheritance >= 0;
    time_to_inheritance(I) = time_to_inheritance(I) + 1;
    I = time_to_inheritance <= -5;
    time_to_inheritance(I) = -5;    
        
    % 2. vectorize
    y = sim.(varname)(:)+sim.h(:);
    age = age(:);
    time_to_inheritance = time_to_inheritance(:);
    
    % 3. select sample
    I = (age <= 54) & (time_to_inheritance <= 10);
    y = y(I);
    age = age(I);
    time_to_inheritance = time_to_inheritance(I);
    
    % 4. count number of obs.
    N = numel(y);
    
    % 5. build X
    X = nan(N,43);
    j = 1;
    
        % a. constant
        X(:,1) = 1;
        
        % b. age dummies
        for i = 26:54 % age = 25 is excluded
            j = j+1;
            X(:,j) = (age == i);
        end
        
        % c. time to inheritance dummies
        for i = [-5:-2 1:10]
            j = j+1;
            X(:,j) = (time_to_inheritance == i);
        end
        
        assert(j == size(X,2))
        
    % 6. weights
    w = nan(N,1);
    for i = 25:54
        for j = -5:10
            I = (age == i) & (time_to_inheritance == j);
            w(I) = N/sum(I);
        end
    end
    assert(sum(isnan(w)) == 0);
    sqrt_w = sqrt(w);
    
    % 7. transformation
    y_tilde = y.*sqrt_w;
    X_tilde = X.*repmat(sqrt_w,[1 size(X,2)]);
    
    % 8. regression
    beta = regress(y_tilde,X_tilde);
    
    if print
        
        %print
        j = 1;
        fprintf('              constant: %6.3f\n',beta(j));
        for i = 26:54
            j = j+1;
            fprintf('                age %2d: %6.3f\n',i,beta(j));
        end
        for i = [-5:-2 1:10]
            j = j+1;
            fprintf('time to inheritance %2d: %6.3f\n',i,beta(j));
        end
        
    end
    
    % 9. only 9 next-to-last elements
    beta = beta(end-9:end-1);
    if strcmp(varname,'ab')
        par.IRF.ab.norm = beta(1);
        par.IRF.ab.beta = beta(1:end)/beta(1);
    else
        par.IRF.a.beta = beta(1:end)/par.IRF.ab.norm;        
    end
    
end
function [par,dist] = solve_and_simulate(par,do_rough)
        
        global it;
        if do_rough && it < par.iter_rough
            par = setup.adj_grid(par,-2);
        end       
        
    % 1. solve and simulate   
    if par.het == 0
        
        [par,sol] = model.solve(par);
        [par,sim] = model.simulate(par,sol,1);
        
    else
           
        % a. setup
        par = estimate.het_vec(par);  
        par = model.create_grids(par);
        par = model.setup_simulate(par);
        Nsim_het = par.Nsim/par.het;
                       
        % b. loop
        as = cell(par.het,1);       
        bs = cell(par.het,1);         
        abs = cell(par.het,1);
        hs = cell(par.het,1);           
        for i = 1:par.het
            
            % i. parameters
            parh = par;            
            parh.(parh.hetvar) = parh.het_vec(i);              
            parh.Nsim = Nsim_het;
            ib = (i-1)*Nsim_het + 1;
            iu = i*Nsim_het;            
            
            % ii. solve
            [parh,sol] = model.solve(parh);
            
            % iii. simulate
            parh = model.select(parh,ib,iu);                      
            [parh,sim] = model.simulate(parh,sol,0);
            
            % iv. save
            abs{i} = sim.ab;
            as{i} = sim.a;                
            bs{i} = sim.b;                        
            hs{i} = sim.h;
            
        end
			
            par.time_solve = parh.time_solve*parh.het;
            par.time_simulate = parh.time_simulate*parh.het;
		
        clearvars parh sol sim;
        
        % c. combine
        sim.a = cell2mat(as);
        sim.h = cell2mat(hs);        
        sim.b = cell2mat(bs);        
        sim.ab = cell2mat(abs);

    end
    
    % 2. IRF
    if par.target_IRF == 1
        par = estimate.IRF_reg(par,sim,'ab',0);
    end
    
    % 3. distance
    dist = estimate.distance(par,sim,par.moms.vars);

end
function [par] = updatepar(par,parnames,parvals)

    for i = 1:numel(parnames)           
        parname = parnames{i};
        parval  = parvals(i);            
        par.(parname) = parval;            
    end      

end 
function [dist_all] = wrapper(x,parnames,par)
    
    global show_time it
   
    % 1. update
    par  = estimate.updatepar(par,parnames,x);    
    fprintf('%3d: x = [',it); it = it+1;
    for i = 1:numel(x)
        fprintf(' %6.4f',x(i));
    end
    
    % 2. distance
    [par,dist] = estimate.solve_and_simulate(par,1);
    if par.target_IRF == 1 && par.het == 0
        fprintf('], dist = %12.8f [LCP %12.8f, IRF: %12.8f]\n',dist.all,dist.ab.p50,dist.IRF);
    elseif par.target_IRF == 1
        fprintf('], dist = %12.8f [LCP %12.8f, IRF: %12.8f, IQ: %12.8f ]\n',dist.all,dist.ab.p50,dist.IRF,dist.ab.iq);
    elseif par.het > 1
        fprintf('], dist = %12.8f [LCP %12.8f, IQ: %12.8f ]\n',dist.all,dist.ab.p50,dist.ab.iq);        
    else
        fprintf('], dist = %12.8f\n',dist.all);        
    end
    dist_all = dist.all;
    
    % 3. print
    if show_time 
        show_time = false;
        fprintf('    time solve    = %5.2f secs\n',par.time_solve);
        fprintf('    time simulate = %5.2f secs\n',par.time_simulate);
    end
   
end


%%%%%%%%%%%
% 2. main %
%%%%%%%%%%%

function par = run(par,parnames,do_post)

    global show_time it
    show_time = true;
    it = 1;
        
    fprintf('estimate.run():\n')    
        
    % 1. objective function
    f = @(x) estimate.wrapper(x,parnames,par);

        options_fminsearch = optimset('Display','none',...
                   'TolX',par.tol_x,'TolFun',par.tol_fun,...
                   'MaxFunEvals', par.maxiter);
               
    % 2. initial guess
    initial_guess = nan(numel(parnames),1);
    for i = 1:numel(parnames)
        initial_guess(i) = par.(parnames{i});
    end
    
    % 3. estimate
    [x,par.fval,par.exitflag] = fminsearch(f,initial_guess,options_fminsearch);

        if it < par.iter_rough
            it = par.iter_rough;            
            [x,par.fval,par.exitflag] = fminsearch(f,x,options_fminsearch);
        end
        
    for i = 1:numel(x)
        fprintf('%s is estimated to be %g\n',parnames{i},x(i));
    end
    fprintf('\n')
    
    % 4. update par
    par = estimate.updatepar(par,parnames,x);
    
    % 5. post-estimation
    if do_post
        par = estimate.post_estimation(par);
    end
    
end
function [beta,zeta] = starting_values(par)
   
    global it;
    it = 0;
    
    fprintf('starting values:\n')
    
    % 1. less dense grid
    par = setup.adj_grid(par,-2);

    % 2. evaluate
    beta_zeta = [0.85 3.0; 0.90 2.0; 0.91 1.9; 0.92 1.8; 0.93 1.7; 0.94 1.6; 0.95 1.5; 0.96 1.4; 0.98 1.2];       
    if contains(par.prefix,'kv_')
        beta_zeta(:,1) = beta_zeta(:,1) - 0.03;
    end    
    dist_all = nan(size(beta_zeta ,1),1);             
    for i = 1:size(beta_zeta,1)
        par.beta = beta_zeta(i,1);
        par.zeta = beta_zeta(i,2);
        [~,dist] = estimate.solve_and_simulate(par,0);
        dist_all(i) = dist.all;
        fprintf('  beta = %6.4f, zeta = %6.4f: dist = %12.8f\n',par.beta,par.zeta,dist.all);
    end
    fprintf('\n');
    
    % 3. find minimum
    [~,imin] = min(dist_all);       
    beta = beta_zeta(imin,1);
    zeta = beta_zeta(imin,2);
        
end
function [het_est,delta_est,zeta] = starting_values_het(par)
   
    global it;
    it = 0;
    
    fprintf('starting values:\n')
    
    % 1. less dense grid
    par = setup.adj_grid(par,-2);

    % 2. evaluate
    delta_est = -1.0;
    het_est_zeta = [2.0 1.8; 2.5 1.6; 3.0 1.4; 3.5 1.2; 4.0 1.0];        
    dist_all = nan(size(het_est_zeta,1),1);             
    for i = 1:size(het_est_zeta,1)
        par.het_est = het_est_zeta(i,1);
        par.delta_est = delta_est;        
        par.zeta = het_est_zeta(i,2);
        [~,dist] = estimate.solve_and_simulate(par,0);
        dist_all(i) = dist.all;
        fprintf('  het_est = %6.4f, zeta = %6.4f: dist = %12.8f\n',par.het_est,par.zeta,dist.all);
    end
    fprintf('\n');
    
    % 3. find minimum
    [~,imin] = min(dist_all);       
    het_est = het_est_zeta(imin,1);
    zeta = het_est_zeta(imin,2);
            
end

%%%%%%%%%%%%%%%%%%%%%%
% 3. post-estimation %
%%%%%%%%%%%%%%%%%%%%%%

function par = post_estimation(par)
        
    % 1. update settings
    par.Nsim = par.post_Nsim;    
    par.do_forget_inh = par.do_post_estimation_do_forget_inh;

    % 2. solve and simulate
    par = model.create_grids(par);    
    if par.het == 0
        
        if par.do_post_save_to_disc == 1
            par.save_to_disc = par.Nsim;        
        end
        [par,sol] = model.solve(par);
        [par,sim] = model.simulate(par,sol,1);
        
    else
                
        % a. setup
        par = estimate.het_vec(par);    
        par = model.create_grids(par);        
        par = model.setup_simulate(par);
        Nsim_het = par.Nsim/par.het;
                                   
        % b. loop
        Ps = cell(par.het,1);
        As = cell(par.het,1);
        Bs = cell(par.het,1);
        Hs = cell(par.het,1);
        Cs = cell(par.het,1);
        ds = cell(par.het,1);
        zs = cell(par.het,1);
        Ns = cell(par.het,1);        
        MPCs = cell(par.het,1);                
        for i = 1:par.het

            % i. parameters
            parh = par;
            parh.prefix = sprintf('%s_%d',parh.prefix,i);            
            parh.(parh.hetvar) = parh.het_vec(i);
            parh.Nsim = Nsim_het;
            if par.do_post_save_to_disc == 1
                parh.save_to_disc = Nsim_het;            
            end
            ib = (i-1)*Nsim_het + 1;
            iu = i*Nsim_het;        
            
            % ii. solve
            [parh,sol] = model.solve(parh);
            rmdir(sprintf('figs_tabs\\%s',parh.prefix));
            
            % iii. simulate
            parh = model.select(parh,ib,iu);
            [parh,sim] = model.simulate(parh,sol,0);
            
            % iv. save
            Ps{i} = sim.P;
            As{i} = sim.A;
            Bs{i} = sim.B;
            Hs{i} = sim.H;
            Cs{i} = sim.C;
            ds{i} = double(sim.d);
            zs{i} = double(sim.z);
            Ns{i} = sim.N;
            MPCs{i} = sim.MPC;            
                        
        end
        
            par.time_solve = parh.time_solve*parh.het;
            par.time_simulate = parh.time_simulate*parh.het;
            
        clearvars parh sol sim;
                            
        % c. copy standard
        if par.do_post_save_to_disc == 1
            commandstr = sprintf('copy csv\\%s_1.txt',par.prefix);
            for i = 2:par.het
                commandstr = sprintf('%s+csv\\%s_%d.txt',commandstr,par.prefix,i);
            end
            commandstr = sprintf('%s csv\\%s.txt',commandstr,par.prefix);
            system(commandstr);
        end
        for i = 1:par.het
            delete(sprintf('csv\\%s_%d.txt',par.prefix,i));
        end        

        % d. copy unexp
        if par.do_post_save_to_disc == 1 && par.do_post_estimation_do_forget_inh
            commandstr = sprintf('copy csv\\%s_1_unexp.txt',par.prefix);
            for i = 2:par.het
                commandstr = sprintf('%s+csv\\%s_%d_unexp.txt',commandstr,par.prefix,i);
            end
            commandstr = sprintf('%s csv\\%s_unexp.txt',commandstr,par.prefix);
            system(commandstr);
        end            
        for i = 1:par.het
            delete(sprintf('csv\\%s_%d_unexp.txt',par.prefix,i));
        end  
        
        % e. combine
        sim.P = cell2mat(Ps);
        sim.A = cell2mat(As);
        sim.B = cell2mat(Bs);
        sim.H = cell2mat(Hs);
        sim.C = cell2mat(Cs);
        sim.d = cell2mat(ds);
        sim.z = cell2mat(zs);
        sim.N = cell2mat(Ns);        
        sim.MPC = cell2mat(MPCs);     
        
    end    
    
    % 3. more variables
    sim = model.calc_additional_sim_var(par,sim);
    
        % a. positive B and H and h_cond
        sim.B_pos = nan(1,par.TR);
        sim.H_pos = nan(1,par.TR);
        sim.h_cond = nan(1,par.TR);          
        for t = 1:par.TR
            I = sim.B(:,t) > 0;
            sim.B_pos(t) = mean(I);
            I = sim.H(:,t) > 0;
            sim.H_pos(t) = mean(I);    
            sim.h_cond(t) = mean(sim.H(I,t)./sim.P(I,t));
        end
                       
        % b. interquartile range
        sim.ab_iq = nan(1,par.TR);
        for t = 1:par.TR
            sim.ab_iq(t) = (prctile(sim.ab(:,t),75)-prctile(sim.ab(:,t),25));
        end
    
        % c. IRF
        par = estimate.IRF_reg(par,sim,'ab',0);
        if par.NN > 1
            par = estimate.IRF_reg(par,sim,'a',0);
        end
        
        % d. life-cycle profiles
        par.mean_a = mean(sim.a);
        par.median_a = median(sim.a);
        par.mean_ab = mean(sim.ab);
        par.median_ab = median(sim.ab);
        par.mean_MPC = mean(funs.vec(sim.MPC(:,5:par.TR)));
        par.median_MPC = median(funs.vec(sim.MPC(:,5:par.TR)));
        
    % 3. fit info    
    par.target_IRF = 1;
    par.dist = estimate.distance(par,sim,{{'a','mean'},{'a',25},{'a',50},{'a',75},...
                                          {'b','mean'},{'b',25},{'b',50},{'b',75},...
                                          {'ab','mean'},{'ab',25},{'ab',50},{'ab',75},...                                          
                                          {'ab','iq'}});
    par.target_IRF = 0;                                      
    par = estimate.info(par);
    
    % 4. figures
    if par.do_post_figs
        
        % a. fit
        estimate.fit_fig(par,sim,'P','mean','$P_t$');
        estimate.fit_fig(par,sim,'ab',50,'$a_t + b_t$');
        estimate.fit_fig(par,sim,'ab','iq','');    
        estimate.fit_fig(par,sim,'h_cond','mean','$h_t$');
        estimate.fit_fig(par,sim,'H_pos','mean','Parent die');        
        close all;
        
        % b. IRFS
        estimate.IRF(par,'ab');
        if par.NN > 1
            estimate.IRF(par,'a');
        end
        
        % c. inheritance     
        if par.het == 0
            figs.dist_inh(par)
        end
                
        % d. marginal prob.
        if par.het == 0
            sim.inh_p = par.inh_p(1:par.TR)';   
            sim.inh_p(sim.inh_p == 1) = nan;
            figs.lcp('inh_p','mean','probability',par,sim)
        end
        
        close all;
                
    end
        
    % 5. save space when saving 
    
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

end
function par = info(par)
    
    fprintf(' Parameters:\n')
        if par.het == 0
            fprintf(' %-15s: %4.3f\n','beta',par.beta);
        else
            fprintf(' %-15s: [%4.3f;%4.3f]\n','beta',par.het_vec(1),par.het_vec(end));
        end       
    pars = {'rho','sigma','zeta','alpha_tilde','alpha'};
    for i = 1:numel(pars)
        fprintf(' %-15s: %4.3f\n',pars{i},par.(pars{i}))
    end    
    fprintf('\n')
    fprintf(' Fit:\n')
    fprintf(' %-10s: %6.4f\n','LCP',par.dist.ab.p50);
    fprintf(' %-10s: %6.4f\n','IRF',par.dist.IRF);
    fprintf(' %-10s: %6.4f\n','IQ',par.dist.ab.iq);
    
end

function [] = fit_fig(par,sim,vary,method,vary_latex)

    if isnumeric(method)
        fig = figure('Name',sprintf('fit_%s_p%d',vary,method));
    else
        fig = figure('Name',sprintf('fit_%s_%s',vary,method));
    end
    
    hold('on') ;

    % 1. simulation
    x = (1:par.TR)+par.age_min;
    x_data = (1:1:par.TR)+par.age_min;
    
    if strcmp(method,'iq') == 1
        y_sim  = sim.(sprintf('%s_iq',vary));
        y_data = par.moms.(sprintf('%s_iq',vary));
    elseif strcmp(method,'mean') == 1
        y_sim  = mean(sim.(vary),1);  
        y_data = par.moms.(sprintf('%s_mean',vary));      
    else
        y_sim  = prctile(sim.(vary),method,1);  
        y_data = par.moms.(sprintf('%s_p%d',vary,method));
    end
    
    h = plot(x,y_sim(1:par.TR),'o','MarkerSize',3,'DisplayName','simulation');
    set(h, 'MarkerFaceColor', get(h, 'Color'));

    % 2. data
    h = plot(x_data,y_data,'s','MarkerSize',5','DisplayName','targets');
    set(h, 'MarkerFaceColor', get(h, 'Color'));

    % 3. layout
    if strcmp(vary,'a') == 1 
        ylim([-1 1.0]);        
    elseif strcmp(vary,'b') == 1 || strcmp(vary,'ab') == 1
        ylim([-1 5.0]);
    elseif strcmp(vary,'P') == 1
        ylim([0 3]);
    elseif strcmp(vary,'B_pos') == 1
        ylim([0 1]);             
    end    
    
    xlabel('age');
    ylabel(vary_latex);
    legend('Location','best');        
    box('on');
    grid on;

    funs.printfig(par,fig);

end
function [] = IRF(par,varname)
    
    fig = figure('Name',sprintf('IRF_%s',varname));
    hold('on') ;
        
        % data
        h = plot(1:9,par.IRF.(varname).beta_emp,'-o','MarkerSize',5,'DisplayName','data');
        set(h, 'MarkerFaceColor', get(h, 'Color'));
        
        % model
        h = plot(1:9,par.IRF.(varname).beta,'-o','MarkerSize',5,'DisplayName','model');
        set(h, 'MarkerFaceColor', get(h, 'Color'));
    
    ylim([0.0 1.2]);
    
    % layout
    legend('Location','northeast')
    xlabel('Years since inheritance');
    ylabel('Wealth change ratio');
    box('on');
    grid on;
    
    funs.printfig(par,fig);
    
end
function [] = overview_figs(models,prefix,do_IRF_liq)
    
    % a. setup
    fig_median = figure('Name',sprintf('%s_median',prefix));
    ax_median = axes;
    xlabel(ax_median,'age');
    ylabel(ax_median,'Years of permanent income');    
    grid(ax_median,'on');
    hold(ax_median,'on');

    fig_nr = figure('Name',sprintf('%s_nr',prefix));
    ax_nr = axes;
    xlabel(ax_nr,'age');
    ylabel(ax_nr,'Years of permanent income');
    grid(ax_nr,'on');
    hold(ax_nr,'on');
    
    fig_precau = figure('Name',sprintf('%s_precau',prefix));
    ax_precau = axes;
    xlabel(ax_precau,'age');
    ylabel(ax_precau,'Years of permanent income');
    grid(ax_precau,'on');
    hold(ax_precau,'on');

    fig_retire = figure('Name',sprintf('%s_retire',prefix));
    ax_retire = axes;
    xlabel(ax_retire,'age');
    ylabel(ax_retire,'Years of permanent income');
    grid(ax_retire,'on');
    hold(ax_retire,'on');
    
    fig_IRF = figure('Name',sprintf('%s_IRF',prefix));
    ax_IRF = axes;
    xlabel(ax_IRF,'Years since inheritance');
    ylabel(ax_IRF,'Wealth change ratio');
    grid(ax_IRF,'on');
    hold(ax_IRF,'on');
    ylim(ax_IRF,[0.0 1.2]);
    
    if do_IRF_liq
        fig_IRF_liq = figure('Name',sprintf('%s_IRF_liq',prefix));
        ax_IRF_liq = axes;
        xlabel(ax_IRF_liq,'Years since inheritance');
        ylabel(ax_IRF_liq,'Wealth change ratio');
        grid(ax_IRF_liq,'on');
        hold(ax_IRF_liq,'on');
        ylim(ax_IRF_liq,[0.0 1.2]);
    end
    
    % b. lines
    for i = 1:numel(models)

        name = models{i}{1};
        latexname = models{i}{2};
        load(sprintf('estimates/%s',name),'par');

        I = 5:par.TR;
        x = I+par.age_min;
        
        % add to median      
        y = par.median_ab;
        h = plot(ax_median,x,y(I),'-o','DisplayName',latexname);
        set(h, 'MarkerFaceColor', get(h, 'Color')); 
        if i == numel(models)
            y = par.moms.ab_p50;
            h = plot(ax_median,x,y(I),'-o','Color','black','DisplayName','Data');
            set(h, 'MarkerFaceColor', get(h, 'Color')); 
        end  
        
        if par.adjcost_share == 0

            % add to precau
            y = par.mean_ab_nr;
            h = plot(ax_nr,x,y(I),'-o','DisplayName',latexname);
            set(h, 'MarkerFaceColor', get(h, 'Color')); 

            % add to precau
            y = par.mean_ab_nr - par.mean_ab_nvr;
            h = plot(ax_precau,x,y(I),'-o','DisplayName',latexname);
            set(h, 'MarkerFaceColor', get(h, 'Color')); 

            % add to retire
            y = par.mean_ab_nv - par.mean_ab_nvr;
            h = plot(ax_retire,x,y(I),'-o','DisplayName',latexname);
            set(h, 'MarkerFaceColor', get(h, 'Color')); 
        
        end
        
        % add to IRF
        if i == numel(models)
            h = plot(ax_IRF,1:9,par.IRF.ab.beta_emp,'-o','Color','black','MarkerSize',5,'DisplayName','Data');
            set(h, 'MarkerFaceColor', get(h, 'Color'));            
        end
        h = plot(ax_IRF,1:9,par.IRF.ab.beta,'-o','MarkerSize',5,'DisplayName',latexname);
        set(h, 'MarkerFaceColor', get(h, 'Color'));
        
        % add to IRF_liq
        if do_IRF_liq
            h = plot(ax_IRF_liq,1:9,par.IRF.a.beta,'-o','MarkerSize',5,'DisplayName',latexname);
            set(h, 'MarkerFaceColor', get(h, 'Color'));        
            if i == numel(models)
                h = plot(ax_IRF_liq,1:9,par.IRF.a.beta_emp,'-o','Color','black','MarkerSize',5,'DisplayName','Data');
                set(h, 'MarkerFaceColor', get(h, 'Color'));            
            end            
        end
    
    end

    % c. legends
    legend(ax_median,'Location','northwest');
    
    % d. print
    par.prefix = '';
    funs.printfig(par,fig_median)
    if par.adjcost_share == 0
        funs.printfig(par,fig_nr)
        funs.printfig(par,fig_precau)
        funs.printfig(par,fig_retire)    
    end
    funs.printfig(par,fig_IRF)
    if do_IRF_liq
        funs.printfig(par,fig_IRF_liq)
    end
    
end
function [] = table_lines(models)
   
    for i = 1:numel(models)

        name = models{i}{1};
        load(sprintf('estimates/%s',name),'par');
        
        % a. strings
            
            est_str = '\textsuperscript{\dag}';
            target_str = '\textsuperscript{\ddag}';
            
            % i. beta
            if par.het == 0
                beta_str = sprintf('%4.3f',par.beta);
                if any(contains(par.estvars,'beta'))
                    beta_str = sprintf('%s%s',beta_str,est_str);
                end
            else
                beta_str = sprintf('[%4.3f;%4.3f]%s',par.het_vec(1),par.het_vec(end));
            end         
            
            % ii. rho
            rho_str = sprintf('%4.2f',par.rho);
            if any(contains(par.estvars,'rho'))
                rho_str = sprintf('%s%s',rho_str,est_str);
            end
        
            % iii. sigma
            sigma_str = sprintf('%4.2f',par.sigma);
            if any(contains(par.estvars,'sigma'))
                sigma_str = sprintf('%s%s',sigma_str,est_str);
            end

            % iv. zeta
            zeta_str = sprintf('%4.2f',par.zeta);
            if any(contains(par.estvars,'zeta'))
                zeta_str = sprintf('%s%s',zeta_str,est_str);
            end

            % v. alpha_tilde
            alpha_tilde_str = sprintf('%4.2f',par.alpha_tilde);
            if any(contains(par.estvars,'alpha_tilde'))
                alpha_tilde_str = sprintf('%s%s',alpha_tilde_str,est_str);
            end

            % vi. alpha
            alpha_str = sprintf('%4.2f',par.alpha);
            if any(contains(par.estvars,'alpha')) && any(contains(par.estvars,'alpha_tilde')) == 0
                alpha_str = sprintf('%s%s',alpha_str,est_str);
            end
            
            % vii. targets
            median_str = sprintf('%4.3f',par.dist.ab.p50);
            median_str = sprintf('%s%s',median_str,target_str);
            IRF_str = sprintf('%4.3f',par.dist.IRF);
            if contains(name,'_IRF') || contains(name,'robust_')
                IRF_str = sprintf('%s%s',IRF_str,target_str);
            end
            if isfield(par.dist,'IRF_liq')
                IRF_liq_str = sprintf('%4.3f',par.dist.IRF_liq);
            else
                IRF_liq_str = '-';
            end       
            iq_str = sprintf('%4.3f',par.dist.ab.iq);
            if par.het >0
                iq_str = sprintf('%s%s',iq_str,target_str);
            end
            
        % b. alphas
        file = fopen(sprintf('%s/figs_tabs/table_lines/%s_alphas.txt',pwd,name),'w');
        fprintf(file,'%s & %s & %s & %s & %s & %s & %s & %s & %4.2f \\\\',...
            beta_str,rho_str,sigma_str,zeta_str,...
            alpha_tilde_str,alpha_str,...
            median_str,IRF_str,par.median_MPC);
        fclose(file);
        
        % c. noalphas
        file = fopen(sprintf('%s/figs_tabs/table_lines/%s_noalphas.txt',pwd,name),'w');
        fprintf(file,'%s & %s & %s & %s & %s & %s & %4.2f \\\\',...
            beta_str,rho_str,sigma_str,zeta_str,...
            median_str,IRF_str,par.median_MPC);
        fclose(file);
        
        % d. irf_liq
        file = fopen(sprintf('%s/figs_tabs/table_lines/%s_iq_irf_liq.txt',pwd,name),'w');
        fprintf(file,'%s & %s & %s & %s & %s & %s & %s & %s & %4.2f \\\\',...
            beta_str,rho_str,sigma_str,zeta_str,...
            median_str,iq_str,IRF_str,IRF_liq_str,par.median_MPC);
        fclose(file);        
    
    end
    
end

end
end