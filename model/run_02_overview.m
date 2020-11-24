clc;
clear;
close all;
funs.layout();

%% 1. tables

models = {};

files = dir('estimates/*.mat');
for i = 1:numel(files)
    filename = files(i).name;
    name = strrep(filename,'.mat','');
    models{end+1} = {name};
end

estimate.table_lines(models);


%% 2. figures

close all;

% rho
models = {};
models{end+1} = {'bs_ez_LCP_rho_1','$\rho = 1.5$'};
models{end+1} = {'bs_ez_LCP_rho_2','$\rho = 2.0$'};
models{end+1} = {'bs_ez_LCP_rho_3','$\rho = 4.0$'};
models{end+1} = {'bs_ez_LCP_rho_4','$\rho = 6.0$'};
estimate.overview_figs(models,'bs_ez_LCP_rho',0);

% main
models = {};
models{end+1} = {'bs_ez_LCP_rho_2','LCP: $\rho = 2.0$'};
models{end+1} = {'bs_ez_LCP_IRF','LCP+IRF: Free $\rho$'};
models{end+1} = {'bs_ez_LCP_IRF_per','LCP+IRF: Free $\tilde{\alpha}$'};
models{end+1} = {'bs_ez_LCP_IRF_var','LCP+IRF: Free $\alpha$'};

estimate.overview_figs(models,'bs_main',0);

% het
models = {};
models{end+1} = {'bs_het_LCP_rho_2','LCP: $\rho = 2.0$'};
models{end+1} = {'bs_het_LCP_IRF','LCP+IRF: Free $\rho$'};

estimate.overview_figs(models,'bs_het',0);


% het
models = {};
models{end+1} = {'kv_ez_LCP_rho_2','LCP: $\rho = 2.0$'};
models{end+1} = {'kv_ez_LCP_IRF','LCP+IRF: Free $\rho$'};

estimate.overview_figs(models,'kv',1);
