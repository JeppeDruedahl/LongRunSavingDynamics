classdef funs
methods(Static)
    
function [] = layout()

    % set layout parameters
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');
    set(groot, 'defaultAxesFontSize', 12); 
    warning('off','all')

end  
function y = vec(x)
    y = x(:);
end
function [x, w] = GaussHermite(n)

    i   = 1:n-1;
    a   = sqrt(i/2);
    CM  = diag(a,1) + diag(a,-1);
    [V, L]   = eig(CM);
    [x, ind] = sort(diag(L));
    V       = V(:,ind)';
    w       = sqrt(pi) * V(:,1).^2;

end
function [x, w] = GaussHermite_lognorm(sigma,n)

    [x,w] = funs.GaussHermite(n);

    x = exp(x*sqrt(2)*sigma-0.5*sigma^2);
    w = w./sqrt(pi);

    % assert a mean of one
    assert(1-sum(w.*x) < 1e-8)

end
function x = nonlinspace(lo,hi,n,phi)
    % recursively constructs an unequally spaced grid.
    % phi > 1 -> more mass at the lower end of the grid.
    % lo can be a vector (x then becomes a matrix).

    x      = NaN(n,length(lo));
    x(1,:) = lo;
    for i = 2:n
        x(i,:) = x(i-1,:) + (hi-x(i-1,:))./((n-i+1)^phi);
    end

end
function [] = printfig(par,figin)

    fig = figure(figin);
    fig.PaperUnits = 'centimeters';   
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 16 12];
    fig.PaperSize = [16 12];

    folder = ['figs_tabs\' par.prefix];
    if exist(folder,'dir') == 0
        mkdir(folder);
    end
    filename = [folder '\' get(fig,'name') ''];
    print('-dpng',['' filename '.png']);

end
function [] = printfig_test_suite(par,figin)

    fig = figure(figin);
    fig.PaperUnits = 'centimeters';   
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 16 12];
    fig.PaperSize = [16 12];

    folder = ['figs_tabs\test_suite'];
    if exist(folder,'dir') == 0
        mkdir(folder);
    end
    filename = [folder '\' par.prefix '_' get(fig,'name') ''];
    print('-dpng',['' filename '.png']);

end

end
end

