classdef figs
methods(Static)

%%%%%%%%%%%%%%%%%%%%%%%%
% 1. simulation figure %
%%%%%%%%%%%%%%%%%%%%%%%%

function [] = dist_inh(par)
    
    Ts = 5:5:par.T-30;

    fig = figure('Name','inh_cycle');
    hold('on');
    c = hsv(numel(Ts));
    i = 1;
    for t = Ts
        bar(par.age_min + (1:par.T),par.inh_p_full(t,:),...
            'DisplayName',sprintf('age = %2d',par.age_min+t),...        
            'EdgeColor','none','FaceColor',c(i,:),'FaceAlpha',0.6);
        i = i + 1;
    end

    % layout
    xlim([30 70])  
    xlabel('age when parent dies');
    ylabel('probability');
    legend('Location','best');
    box('on');
    grid on;

    funs.printfig(par,fig);

end
function [] = lcp(vary,method,vary_latex,par,sim)
    
    if isnumeric(method)
        fig = figure('Name',sprintf('lcp_%s_p%d',vary,method));
    else
        fig = figure('Name',sprintf('lcp_%s_%s',vary,method));        
    end
    hold('on') ;

    x = (1:par.TR)+par.age_min;

    if strcmp(method,'iq') == 1
        y_sim  = sim.(sprintf('%s_iq',vary));
    elseif strcmp(method,'mean') == 1
        y_sim  = mean(sim.(vary),1);  
    else
        y_sim  = prctile(sim.(vary),method,1);  
    end
    
    h = plot(x,y_sim,'o','MarkerSize',3);
    set(h, 'MarkerFaceColor', get(h, 'Color'));
    ylim('auto');
    line([par.TR,par.TR]+par.age_min,[-1 6],'color','black')

    % layout
    if strcmp(vary,'A') == 1 || strcmp(vary,'B') == 1 || strcmp(vary,'AB') == 1
        ylim([-5 6])
    elseif strcmp(vary,'C') == 1 || strcmp(vary,'P') == 1 || strcmp(vary,'Y') == 1
        ylim([0.8 2.2]);
    elseif strcmp(vary,'H') == 1
        ylim([0 0.2]);  
    elseif strcmp(vary,'H_cond') == 1
        ylim([0 6]);          
    elseif strcmp(vary,'z') == 1 || strcmp(vary,'d') == 1 ||  strcmp(vary,'inh_p') == 1 ...
            || strcmp(vary,'age_parent_dies') == 1 || strcmp(vary,'MPC') == 1
        ylim([0 1]);        
    end
    xlabel('age');
    ylabel(vary_latex);
    box('on');
    grid on;
    
    funs.printfig(par,fig);
    
end

%%%%%%%%%%%%%%
% 2. compare %
%%%%%%%%%%%%%%

function [] = compare_lcp(prefix,vary,method,vary_latex,pars,sims,names)
    
    par = pars{1};
    
    if isnumeric(method)
        fig = figure('Name',sprintf('lcp_%s_%s_p%d',prefix,vary,method));
    else
        fig = figure('Name',sprintf('lcp_%s_%s_%s',prefix,vary,method));        
    end
    hold('on') ;

    x = (1:par.TR)+par.age_min;
        
    for i = 1:numel(pars)
        
        sim = sims{i};
        if strcmp(method,'iq') == 1
            y_sim  = sim.(sprintf('%s_iq',vary));
        elseif strcmp(method,'mean') == 1
            y_sim  = mean(sim.(vary),1);  
        else
            y_sim  = prctile(sim.(vary),method,1);  
        end

        if i == 1 
            h = plot(x,y_sim,'s','Color','black','MarkerSize',5,'DisplayName',names{i});
        else
            h = plot(x,y_sim,'--o','MarkerSize',3,'DisplayName',names{i});            
        end
        set(h, 'MarkerFaceColor', get(h, 'Color'));
    
    end
    ylim('auto');

    % layout
    if strcmp(vary,'A') == 1 || strcmp(vary,'B') == 1 || strcmp(vary,'AB') == 1
        ylim([-1 6])
    elseif strcmp(vary,'C') == 1 || strcmp(vary,'P') == 1 || strcmp(vary,'Y') == 1
        ylim([0.8 2.2]);
    elseif strcmp(vary,'H') == 1
        ylim([0 0.2]);  
    elseif strcmp(vary,'H_cond') == 1
        ylim([0 6]);          
    elseif strcmp(vary,'z') == 1 || strcmp(vary,'d') == 1 ||  strcmp(vary,'inh_p') == 1 ...
            || strcmp(vary,'age_parent_dies') == 1 || strcmp(vary,'MPC') == 1
        ylim([0 1]);        
    end
    xlabel('age');
    ylabel(vary_latex);
    box('on');
    grid on;
    legend('show','Location','best');
    
    funs.printfig_test_suite(par,fig);
    
end
function [] = compare_AB(t,prefix,varx_latex,vary_latex,par_base,sim_base,x,pars,sims,name,name_alt)
    
    par = par_base;
    
    fig = figure('Name',sprintf('compare_AB_t%d_%s',t,prefix));
    hold('on') ;

    % base
    if isstruct(sim_base)
        plot(x,mean(sim_base.AB(:,t))*ones(numel(x),1),'-','Color','black','DisplayName',name);
    end
        
    alt = nan(numel(pars),1);
    for i = 1:numel(pars)
        alt(i) = mean(sims{i}.AB(:,t));
    end

    h = plot(x,alt,'--o','MarkerSize',3,'DisplayName',name_alt);
    set(h, 'MarkerFaceColor', get(h, 'Color'));
       
    % layout
    xlabel(varx_latex);
    ylabel(vary_latex);
    box('on');
    grid on;
    if isstruct(sim_base)
        legend('show','Location','best');
    end
    
    funs.printfig_test_suite(par,fig);
    
end

end
end