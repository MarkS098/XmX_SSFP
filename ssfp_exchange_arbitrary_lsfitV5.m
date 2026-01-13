clc; close all; clearvars;

% Load data
load('/home/mark/NMR Data/XmX processing/16-12-2025 acac/peakB_ao','peaks','TR_vals');

% Constants 
R1 = 1/13;          % Longitudinal relaxation rate (s^-1)
FLIP = 23;    % Flip angle (radians)
skip_points = 0;       % points to omit from fitting
peak_choose = 'B';

% Acquisition times
Tacq = TR_vals(1+skip_points:end)*1e-3;  % seconds
M_0 = peaks(1+skip_points:end);          % experimental signal
err = 0.1*abs(M_0);                     % uniform noise

% Parameter bounds
% [MA0, MB0, kex_AB, kex_BC, kex_AC, nuA, nuB, nuC, R2]
lb = [0.01, 0.01, 0.1, 0.1, 0.1, 0, 0, 0, 0.1];   
ub = [0.99, 0.99, 25, 25, 25, 1000, 1000, 1000, 10];
Nstart = 1000;  % number of random starting points within bounds

% Objective function with exact scale
function residuals = objective_function(params, Tacq, My_0, R1, FLIP, peak_choose)
    
    pop = [params(1), params(2)];
    kex = [params(3), params(4), params(5)];
    nu = [params(6), params(7), params(8)];
    R2 = params(9);
    
    if sum(pop) >= 1
        residuals = ones(size(My_0)) * 1e6;
        return;
    end
    
    [M_A, M_B, M_C] = chem_exchange_sim(FLIP, Tacq, pop, nu, kex, R1, R2);
    
    switch upper(peak_choose)
        case 'A', M_model = M_A;
        case 'B', M_model = M_B;
        case 'C', M_model = M_C;
    end
    
    % scale = (M_model(:)'*My_0(:))/(M_model(:)'*M_model(:));
    residuals = M_model(:) - My_0(:);
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
ms = MultiStart('Display','off');
rng default
[best_params,fval,exitflag,output,all_solutions] = run(ms, problem, Nstart);

% Extract parameter matrix from all solutions
param_matrix = vertcat(all_solutions.X);  % Each row = one local minimum
fvals = [all_solutions.Fval];             % Objective (resnorm) values

% Keep only finite, successful solutions
valid_idx = isfinite(fvals);
param_matrix = param_matrix(valid_idx,:);
fvals = fvals(valid_idx);

% Plot histograms for each parameter
param_names = {'MA0','MB0','kAB','kBC','kAC','nuA','nuB','nuC','R2'};

figure;
for i = 1:size(param_matrix,2)
    subplot(3,3,i)
    histogram(param_matrix(:,i),15,'FaceColor',[0.2 0.2 0.8],'EdgeColor','k');
    xlabel(param_names{i},'FontSize',12);
    ylabel('Count','FontSize',12);
    title(sprintf('%s (mean=%.3f)',param_names{i},mean(param_matrix(:,i))));
    grid on
end

% Extract optimized parameters
[M_A_opt, M_B_opt, M_C_opt] = chem_exchange_sim(FLIP, Tacq, [best_params(1),best_params(2)], ...
                                               [best_params(6),best_params(7),best_params(8)], ...
                                               [best_params(3),best_params(4),best_params(5)], R1, best_params(9));
switch upper(peak_choose)
    case 'A', M_opt = M_A_opt;
    case 'B', M_opt = M_B_opt;
    case 'C', M_opt = M_C_opt;
end

% Goodness of fit tests
residuals = M_opt - M_0;
RMSE = sqrt(mean(residuals.^2));

% Display results 
fprintf('MA0 = %.3f, MB0 = %.3f, MC0 = %.3f\n', best_params(1), best_params(2), 1-best_params(1)-best_params(2));
fprintf('kAB = %.3f, kBC = %.3f, kAC = %.3f\n', best_params(3), best_params(4), best_params(5));
fprintf('nuA = %.3f, nuB = %.3f, nuC = %.3f\n', best_params(6), best_params(7), best_params(8));
fprintf('R2 = %.3f\n', best_params(9));
fprintf('RMSE = %.5f\n', RMSE);

% Plot experimental and fitted data
title_str = {sprintf('Fit %s: MA0=%.2f, MB0=%.2f, MC0=%.2f', peak_choose, best_params(1), best_params(2), 1-best_params(1)-best_params(2)), ...
             sprintf('kAB=%.1f, kBC=%.1f, kAC=%.1f', best_params(3), best_params(4), best_params(5)), ...
             sprintf('nuA=%.1f, nuB=%.1f, nuC=%.1f, R2=%.2f', best_params(6), best_params(7), best_params(8), best_params(9))};

figure
plot(Tacq*1e3, M_0, 'b-','LineWidth',2,'DisplayName','Experimental')
hold on
plot(Tacq*1e3, M_opt, 'r--','LineWidth',2,'DisplayName','Fitted')
xlabel('Tacq [ms]')
ylabel('Signal Intensity'); ylim([0,0.17])
legend
title(title_str)
set(gca,'FontSize',12)

% Plot residuals 
figure
errorbar(Tacq*1e3,residuals,err,'o','MarkerSize',6,'CapSize',8,'LineWidth',1)
hold on
yline(0,'--k','LineWidth',1)
xlabel('Tacq [ms]')
ylabel('Residuals')
ylim([-0.05 0.05])
title(sprintf('RMSE = %.5f', RMSE))
set(gca,'FontSize',12)

