clc; close all; clearvars;

% Load data
load('/home/mark/NMR Data/XmX processing/ACAC_325K/ACAC_peakB_offset_arbitrary_325K_FA23.mat','peaks','TR_vals');

% Constants 
R1 = 1/13;          % Longitudinal relaxation rate (s^-1)
FLIP = 23;          % Flip angle
skip_points = 0;    % points to omit from fitting
peak_choose = 'B';

% Acquisition times
Tacq = TR_vals(1+skip_points:end)*1e-3;  % seconds
M_0 = peaks(1+skip_points:end);          % experimental signal
err = 0.1*M_0;

% Parameter bounds
% [MA0, MB0, kex_AB, kex_BC, kex_AC, nuA, nuB, nuC, R2]
lb = [0.1, 0.5, 0.1, 0, 0, 100, 0, 0, 0.1];   
ub = [0.4, 0.99, 100, 0, 0, 1000, 800, 0, 5];
Nstart = 1000;  % number of random starting points within bounds
N_boot = 50; % number of bootstrap runs
n_grid = 50; % chi square map resolution

% Objective function with exact scale
function residuals = objective_function(params, Tacq, My_0, R1, FLIP, peak_choose)
    
pop_raw = [params(1), params(2)];
    kex = [params(3), params(4), params(5)];
    nu  = [params(6), params(7), params(8)];
    R2  = params(9);
    
    % Population constraint
    if nu(3) == 0 && kex(2) == 0 && kex(3) == 0 
        % 2-site regime 
        MC0 = 0;
        total_active = sum(pop_raw);
        MA0 = params(1)/total_active;
        MB0 = params(2)/total_active;
    else
        % 3-site regime
        MA0_raw = params(1);
        MB0_raw = params(2);
        MC0_raw = max(0, 1 - MA0_raw - MB0_raw);
        
        total = MA0_raw + MB0_raw + MC0_raw;
        MA0 = MA0_raw/total;
        MB0 = MB0_raw/total;
        MC0 = MC0_raw/total;
    end

    pop_constrained = [MA0, MB0, MC0];
    [M_A, M_B, M_C] = chem_exchange_sim(FLIP, Tacq, pop_constrained, nu, kex, R1, R2);

    switch upper(peak_choose)
        case 'A', M_model = M_A;
        case 'B', M_model = M_B;
        case 'C', M_model = M_C;
    end
    
    M_model = M_model(:); 
    My_0 = My_0(:);
    scale = (M_model'*My_0)/(M_model'*M_model);
    res_signal = (scale*M_model - My_0)/mean(abs(My_0));

    residuals = res_signal;
end

% lsqnonlin options 
options = optimoptions('lsqnonlin','Display','off','MaxIterations',1000,'TolFun',1e-9,'TolX',1e-9);

% Create a problem for MultiStart
x0 = lb + (ub-lb)/2;  % single initial guess in middle of bounds (MultiStart will generate others)
problem = createOptimProblem('lsqnonlin', ...
    'x0', x0, ...
    'objective', @(p)objective_function(p,Tacq,M_0,R1,FLIP,peak_choose), ...
    'lb', lb, 'ub', ub, ...
    'options', options);

npar = numel(lb);

% Generate Latin Hypercube samples in [lb,ub]
Xlhs = bsxfun(@plus,lb,lhsdesign(Nstart,npar).*(ub-lb));
startSet = CustomStartPointSet(Xlhs);
ms = MultiStart('Display','off','UseParallel',true);
rng default
[best_params,fval,exitflag,output,all_solutions] = run(ms, problem, Nstart);

% Extract parameter matrix from all solutions
param_matrix = vertcat(all_solutions.X);  % Each row = one local minimum
fvals = [all_solutions.Fval];             % Objective (resnorm) values


% Plot histograms for each parameter
param_names = {'MA0','MB0','kAB','kBC','kAC','nuA','nuB','nuC','R2'};

figure
for i = 1:size(param_matrix,2)
    subplot(3,3,i)
    histogram(param_matrix(:,i),15,'FaceColor',[0.2 0.2 0.8],'EdgeColor','k')
    xlabel(param_names{i},'FontSize',12)
    ylabel('Count','FontSize',12)
    title(sprintf('%s (mean=%.3f)',param_names{i},mean(param_matrix(:,i))))
    grid on
end

if best_params(8) == 0 && best_params(4) == 0 && best_params(5) == 0 
    % 2-site regime
    total_act = best_params(1) + best_params(2);
    MA_final = best_params(1)/total_act;
    MB_final = best_params(2)/total_act;
    MC_final = 0;
else
    % 3-site regime
    MA_raw = max(0, best_params(1));
    MB_raw = max(0, best_params(2));
    MC_raw = max(0, 1 - MA_raw - MB_raw);
    total = MA_raw + MB_raw + MC_raw;
    MA_final = MA_raw/total;
    MB_final = MB_raw/total;
    MC_final = MC_raw/total;
end

% Extract optimized parameters
[M_A_opt, M_B_opt, M_C_opt] = chem_exchange_sim(FLIP, Tacq, [MA_final,MB_final], ...
                                               [best_params(6),best_params(7),best_params(8)], ...
                                               [best_params(3),best_params(4),best_params(5)], R1, best_params(9));
switch upper(peak_choose)
    case 'A', M_opt = M_A_opt;
    case 'B', M_opt = M_B_opt;
    case 'C', M_opt = M_C_opt;
end

% Goodness of fit tests
scale = (M_opt(:)'*M_0(:))/(M_opt(:)'*M_opt(:));
M_opt = scale*M_opt;
residuals = M_opt - M_0;
RMSE = rmse(M_opt,M_0);

boot_params = zeros(N_boot, numel(best_params));

% Setup variables and normalization factor
res_vec = objective_function(best_params,Tacq,M_0,R1,FLIP,peak_choose);
M_model_best = M_opt; 
norm_fact = mean(abs(M_0));

% Define opt_fast for the bootstrap fits
opt_fast = optimoptions('lsqnonlin', 'Display', 'off', ...
    'MaxIterations', 100, 'TolFun', 1e-8);

fprintf('Running Bootstrap... ');
num_pts = numel(M_model_best); % Number of data points
parfor b = 1:N_boot
    % Only pull the first N elements of res_vec
    signal_res_only = res_vec(1:num_pts); 
    
    % Shuffle only those signal residuals
    shuffled_res = signal_res_only(randi(num_pts, [num_pts, 1]));
    
    % Now the dimensions match perfectly (N points + N points)
    noisy_M0 = M_model_best(:) + (shuffled_res(:)*norm_fact);
    
    % Fit the noisy data
    p_start = best_params.*(1 + 0.02*randn(size(best_params)));
    p_start = max(min(p_start, ub), lb);
    
    [p_boot, ~, ~] = lsqnonlin(@(p) objective_function(p,Tacq,noisy_M0,R1,FLIP,peak_choose), ...
                            p_start, lb, ub, opt_fast);
    boot_params(b, :) = p_boot;
end

% Final outputs
param_errors = std(boot_params);

% Clean up: set errors for inactive parameters to 0
inactive_mask = (ub - lb) < 1e-5; 
param_errors(inactive_mask) = 0;

fprintf('Done. RMSE = %.2f\n', RMSE);

R2_range = linspace(lb(9), ub(9), n_grid); 
k_indices = [3,4,5]; 
k_label_names = {'k_{AB}', 'k_{BC}', 'k_{AC}'};

all_maps = cell(1,3);
K_mesh_list = cell(1,3);
R_mesh_list = cell(1,3);

for s = 1:3
    k_idx = k_indices(s);
    k_range = linspace(lb(k_idx), ub(k_idx), n_grid);
    [K_mesh, R_mesh] = meshgrid(k_range, R2_range);
    
    K_mesh_list{s} = K_mesh;
    R_mesh_list{s} = R_mesh;
    
    % Flatten the grid for a single parfor pass
    K_flat = K_mesh(:);
    R_flat = R_mesh(:);
    chi2_flat = zeros(size(K_flat));
    
    fprintf('Calculating Map %d/3 (%s)... ', s, k_label_names{s});
    
    % parfor here runs on the entire 2500-point grid at once
    parfor idx = 1:numel(K_flat)
        p_temp = best_params;        
        p_temp(k_idx) = K_flat(idx); 
        p_temp(9) = R_flat(idx);     
        res = objective_function(p_temp,Tacq,M_0,R1,FLIP,peak_choose);
        chi2_flat(idx) = sum(res.^2); 
    end
    
    % Reshape back to the 2D grid
    all_maps{s} = log10(reshape(chi2_flat, n_grid, n_grid));
    fprintf('Done.\n');
end

% Determine global limits for colorbar
global_min = min(cellfun(@(x) min(x(:)), all_maps));
global_max = max(cellfun(@(x) max(x(:)), all_maps));

% Plotting loop 
figure('Name', 'Unified Error Surface Analysis', 'Color', 'w', 'Position', [50, 200, 1600, 500]) 
t = tiledlayout(1, 3, 'TileSpacing', 'Loose', 'Padding', 'Compact');

for s = 1:3
    nexttile
    k_idx = k_indices(s);
    
    contourf(K_mesh_list{s}, R_mesh_list{s}, all_maps{s}, 25, 'LineColor', 'none')
    hold on
    plot(best_params(k_idx), best_params(9), 'r*', 'MarkerSize', 12, 'LineWidth', 2)
    
    colormap(jet)
    clim([global_min, global_max])
    
    c = colorbar;
    c.Label.String = 'log_{10}(\chi^2)';
    
    xlabel(sprintf('%s (s^{-1})', k_label_names{s}), 'FontSize', 11)
    ylabel('R_2 (s^{-1})', 'FontSize', 11)
    title(sprintf('Error Surface: %s vs R_2', k_label_names{s}), 'FontSize', 13)
    grid on
end
 
% Identify active exchange rates and frequencies
active_k = find(ub(3:5) > lb(3:5)) + 2; % Indices 3, 4, or 5
active_nu = find(ub(6:8) > lb(6:8)) + 5; % Indices 6, 7, or 8

% Define map pairings dynamically
if numel(active_k) == 1 && active_k == 3 % 2-Site Case (AB only)
    % Plot k_AB vs nu_A and k_AB vs nu_B
    k_map_indices = [3, 3];
    nu_map_indices = [6, 7];
    map_labels = {'k_{AB} vs \nu_A', 'k_{AB} vs \nu_B'};
else % 3-Site Case 
    % Plot k_AB vs nu_A, k_BC vs nu_B, k_AC vs nu_C
    k_map_indices = [3, 4, 5];
    nu_map_indices = [6, 7, 8];
    map_labels = {'k_{AB} vs \nu_A', 'k_{BC} vs \nu_B', 'k_{AC} vs \nu_C'};
end

n_maps = numel(k_map_indices);
all_maps_nu = cell(1, n_maps);
K_mesh_list_nu = cell(1, n_maps);
NU_mesh_list = cell(1, n_maps);

% Dynamic data generation loop
for s = 1:n_maps
    k_idx = k_map_indices(s);
    nu_idx = nu_map_indices(s);
    
    k_range = linspace(lb(k_idx), ub(k_idx), n_grid);
    nu_range = linspace(lb(nu_idx), ub(nu_idx), n_grid);
    [K_mesh, NU_mesh] = meshgrid(k_range, nu_range);
    
    K_mesh_list_nu{s} = K_mesh;
    NU_mesh_list{s} = NU_mesh;
    
    K_flat = K_mesh(:);
    NU_flat = NU_mesh(:);
    chi2_flat = zeros(size(K_flat));
    
    fprintf('Calculating Dynamic Map %d/%d (%s)... ', s, n_maps, map_labels{s});
    
    parfor idx = 1:numel(K_flat)
        p_temp = best_params;        
        p_temp(k_idx) = K_flat(idx); 
        p_temp(nu_idx) = NU_flat(idx);     
        res = objective_function(p_temp,Tacq,M_0,R1,FLIP,peak_choose);
        chi2_flat(idx) = sum(res.^2); 
    end
    all_maps_nu{s} = log10(reshape(chi2_flat, n_grid, n_grid));
    fprintf('Done.\n');
end

% Determine global limits for colorbar
global_min_nu = min(cellfun(@(x) min(x(:)), all_maps_nu));
global_max_nu = max(cellfun(@(x) max(x(:)), all_maps_nu));

% Dynamic plotting loop 
all_names = cell(1, 9);
all_names{1} = 'M_{A0}';
all_names{2} = 'M_{B0}';
all_names{3} = 'k_{AB}';
all_names{4} = 'k_{BC}';
all_names{5} = 'k_{AC}';
all_names{6} = '\nu_A';
all_names{7} = '\nu_B';
all_names{8} = '\nu_C';
all_names{9} = 'R_2';

figure('Name', 'Dynamic Frequency vs Exchange Analysis', 'Color', 'w', 'Position', [50, 50, 500*n_maps, 500]); 
t2 = tiledlayout(1, n_maps, 'TileSpacing', 'Loose', 'Padding', 'Compact');

for s = 1:n_maps
    nexttile
    k_idx = k_map_indices(s);
    nu_idx = nu_map_indices(s);
    
    contourf(K_mesh_list_nu{s}, NU_mesh_list{s}, all_maps_nu{s}, 25, 'LineColor', 'none')
    hold on
    plot(best_params(k_idx), best_params(nu_idx), 'r*', 'MarkerSize', 12, 'LineWidth', 2)
    
    colormap(jet)

    % Set limits based on the specific map's data for better contrast
    clim([min(all_maps_nu{s}(:)), max(all_maps_nu{s}(:))]);

    c = colorbar;
    c.Label.String = 'log_{10}(\chi^2)';

    xlabel_str = sprintf('%s (s^{-1})', all_names{k_idx});
    ylabel_str = sprintf('%s (Hz)', all_names{nu_idx});

    xlabel(xlabel_str, 'FontSize', 12, 'Interpreter', 'tex');
    ylabel(ylabel_str, 'FontSize', 12, 'Interpreter', 'tex');

    title(map_labels{s}, 'FontSize', 13, 'Interpreter', 'tex')
    grid on
end

% Plot experimental and fitted data
title_str = { ...
    sprintf('Fit %s: M_{A0}=%.3f\\pm%.3f, M_{B0}=%.3f\\pm%.3f, M_{C0}=%.3f', ...
            peak_choose, MA_final, param_errors(1), MB_final, param_errors(2), MC_final), ...
    sprintf('k_{AB}=%.2f\\pm%.2f, k_{BC}=%.2f\\pm%.2f, k_{AC}=%.2f\\pm%.2f s^{-1}', ...
            best_params(3), param_errors(3), best_params(4), param_errors(4), best_params(5), param_errors(5)), ...
    sprintf('\\nu_{A}=%.1f\\pm%.1f, \\nu_{B}=%.1f\\pm%.1f, \\nu_{C}=%.1f\\pm%.1f Hz', ...
            best_params(6), param_errors(6), best_params(7), param_errors(7), best_params(8), param_errors(8)), ...
    sprintf('R_{2}=%.3f\\pm%.3f s^{-1}', best_params(9), param_errors(9)) ...
};

figure
plot(Tacq*1e3, M_0/max(max(M_0),max(M_opt)), 'b-','LineWidth',2,'DisplayName','Experimental')
hold on
plot(Tacq*1e3, M_opt/max(max(M_0),max(M_opt)), 'r--','LineWidth',2,'DisplayName','Fitted')
xlabel('Tacq [ms]')
ylabel('Signal Intensity');
legend
title(title_str,'Interpreter','tex')
set(gca,'FontSize',12)

% Plot residuals 
figure
errorbar(Tacq*1e3,residuals,err,'o','MarkerSize',6,'CapSize',8,'LineWidth',1)
hold on
yline(0,'--k','LineWidth',1)
xlabel('Tacq [ms]')
ylabel('Residuals')
title(sprintf('RMSE = %.e', RMSE))
set(gca,'FontSize',12)
