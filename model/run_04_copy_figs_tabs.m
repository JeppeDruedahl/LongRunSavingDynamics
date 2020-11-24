% 1. destimation path
dst_path = 'output';

% 2. table_lines
dst_dir = sprintf('%s/table_lines',dst_path);
if exist(dst_dir,'dir') > 0
    rmdir(dst_dir,'s')
end
copyfile('figs_tabs/table_lines',dst_dir);

% 3. figures - bs_ez_LCP_IRF
dst_dir = sprintf('%s/bs_ez_LCP_IRF',dst_path);
if exist(dst_dir,'dir') > 0
    rmdir(dst_dir,'s')
end
copyfile('figs_tabs/bs_ez_LCP_IRF',dst_dir);

% 4. figures - test_suite
dst_dir = sprintf('%s/test_suite',dst_path);
if exist(dst_dir,'dir') > 0
    rmdir(dst_dir,'s')
end
copyfile('figs_tabs/test_suite',dst_dir);