%% Q1_Newton_Method.m
% CO1: Formulate the nonlinear operating-temperature equation and solve
% using Newton's Method.
%
% Governing (steady-state energy balance) equation:
%   eta*G*A = U*A*(T - Ta) + eps*sigma*A*[(T+273.15)^4 - (Ta+273.15)^4]
%
% Rearranged as f(T) = 0:
%   f(T) = eta*G*A - U*A*(T - Ta) - eps*sigma*A*[(T+273.15)^4 - (Ta+273.15)^4]
%
% This script can be run standalone (defines its own parameters) or
% called from main.m (which passes shared parameters in, if present).

if ~exist('params', 'var')
    addpath(fullfile(pwd, 'functions'));
    params = struct();
    params.A     = 2;          % collector area, m^2
    params.G     = 700;        % solar irradiance, W/m^2
    params.eta   = 0.70;       % collector efficiency
    params.U     = 6;          % overall heat-loss coefficient, W/m^2K
    params.Ta    = 30;         % ambient temperature, deg C
    params.eps   = 0.90;       % emissivity
    params.sigma = 5.670374419e-8; % Stefan-Boltzmann constant, W/m^2K^4
end

fprintf('\n============================================================\n');
fprintf(' Q1 - CO1: NEWTON''S METHOD (Nonlinear Operating Temperature)\n');
fprintf('============================================================\n');

A = params.A; G = params.G; eta = params.eta; U = params.U;
Ta = params.Ta; eps = params.eps; sigma = params.sigma;

f  = @(T) eta*G*A - U*A*(T - Ta) - eps*sigma*A*((T + 273.15)^4 - (Ta + 273.15)^4);
fp = @(T) -U*A - eps*sigma*A*4*(T + 273.15)^3;

T0      = 60;      % initial guess, deg C
tol     = 1e-6;    % stopping tolerance, deg C
maxIter = 100;

[T_root, iterTable, iterCount, converged] = newton_temperature(f, fp, T0, tol, maxIter);

fprintf('\nNewton Iteration Table:\n');
fprintf(' Iter |     T_n (C)    |      f(T_n)     |     fprime(T_n)   |   T_n+1 (C)    |  |Error| (C)\n');
fprintf('------------------------------------------------------------------------------------------------\n');
for k = 1:iterCount
    fprintf(' %3d  | %14.8f | %15.6e | %17.6e | %14.8f | %10.3e\n', ...
        iterTable(k,1), iterTable(k,2), iterTable(k,3), iterTable(k,4), iterTable(k,5), iterTable(k,6));
end

residual = f(T_root);
fprintf('\nConverged operating temperature T = %.8f deg C\n', T_root);
fprintf('Number of iterations            = %d\n', iterCount);
fprintf('Residual f(T_root)               = %.6e\n', residual);
fprintf('Converged (tol met)              = %d\n', converged);

% --- Convergence / error plot ---
figNewton = figure('Name', 'Q1 Newton Convergence', 'Visible', 'off');
semilogy(iterTable(:,1), iterTable(:,6), '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on;
xlabel('Iteration number');
ylabel('|T_{n+1} - T_n|  (deg C, log scale)');
title('Q1: Newton''s Method Convergence (Error vs Iteration)');
saveas(figNewton, fullfile('figures', 'Q1_Newton_Convergence.png'));

% --- Save results for validation.m / results folder ---
Q1_results.T_root       = T_root;
Q1_results.iterTable    = iterTable;
Q1_results.iterCount    = iterCount;
Q1_results.residual     = residual;
Q1_results.converged    = converged;
save(fullfile('results', 'Q1_results.mat'), 'Q1_results');

fprintf('\n[Q1] Figure saved to figures/Q1_Newton_Convergence.png\n');
fprintf('[Q1] Results saved to results/Q1_results.mat\n');
