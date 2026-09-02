%% Q5_Finite_Difference.m
% CO5: Formulate the one-dimensional heat equation and implement an
% explicit finite-difference scheme. Plot temperature distributions and
% discuss stability.
%
% Governing PDE: dT/dt = alpha * d2T/dx2
% Explicit (FTCS) scheme: T_i^(n+1) = T_i^n + r*(T_(i+1)^n - 2T_i^n + T_(i-1)^n)
% r = alpha*dt/dx^2

if ~exist('params', 'var')
    addpath(fullfile(pwd, 'functions'));
end

fprintf('\n============================================================\n');
fprintf(' Q5 - CO5: EXPLICIT FINITE DIFFERENCE (1-D Heat Equation)\n');
fprintf('============================================================\n');

L         = 1;        % domain length, m (rod/slab representing collector plate direction)
Nx        = 11;       % number of spatial nodes
alpha     = 1e-5;     % thermal diffusivity, m^2/s
dx        = 0.1;      % spatial step, m
dt        = 400;      % time step, s
totalTime = 3600;     % total simulation time, s
T_left    = 80;       % Dirichlet BC at x = 0, deg C
T_right   = 40;       % Dirichlet BC at x = L, deg C
T_init    = 30;       % initial interior temperature, deg C

[x, tSteps, Tmatrix, r] = finite_difference_heat(L, Nx, alpha, dx, dt, totalTime, T_left, T_right, T_init);

dt_max = dx^2 / (2*alpha);

fprintf('\nGrid: Nx = %d nodes, dx = %.3f m, dt = %.1f s, total time = %.0f s\n', Nx, dx, dt, totalTime);
fprintf('Stability number r = alpha*dt/dx^2 = %.6f\n', r);
fprintf('Stability condition r <= 0.5      : %s\n', mat2str(r <= 0.5));
fprintf('Maximum stable time step dt_max    = %.4f s  (used dt = %.1f s)\n', dt_max, dt);

fprintf('\nFull Temperature Matrix (rows = time steps, cols = spatial nodes x = 0..%g m):\n', L);
fprintf('  t(s) \\ x(m) '); fprintf('%8.1f', x); fprintf('\n');
for n = 1:size(Tmatrix,1)
    fprintf('  %8.0f  ', tSteps(n));
    fprintf('%8.3f', Tmatrix(n,:));
    fprintf('\n');
end

% --- Extract requested snapshot times ---
targetTimes = [0 800 1600 2400 3200 3600];
fprintf('\nTemperature Distribution at Requested Snapshot Times:\n');
snapIdx = zeros(size(targetTimes));
for k = 1:numel(targetTimes)
    [~, idx] = min(abs(tSteps - targetTimes(k)));
    snapIdx(k) = idx;
    fprintf('t = %4d s : ', tSteps(idx));
    fprintf('%8.3f', Tmatrix(idx,:));
    fprintf('\n');
end

% --- Validation: finer grid (halved dt) comparison at final time ---
dt_fine = dt / 2;
[~, ~, Tmatrix_fine, r_fine] = finite_difference_heat(L, Nx, alpha, dx, dt_fine, totalTime, T_left, T_right, T_init);
diffFinal = abs(Tmatrix(end,:) - Tmatrix_fine(end,:));
fprintf('\n[Validation] Finer time-step check: dt_fine = %.1f s, r_fine = %.4f\n', dt_fine, r_fine);
fprintf('[Validation] Max |T(dt) - T(dt/2)| at final time = %.6f deg C\n', max(diffFinal));

% --- Plot temperature distributions ---
figFDM = figure('Name', 'Q5 FDM Temperature Distribution', 'Visible', 'off');
hold on;
colors = lines(numel(targetTimes));
for k = 1:numel(targetTimes)
    plot(x, Tmatrix(snapIdx(k), :), '-o', 'Color', colors(k,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('t = %d s', tSteps(snapIdx(k))));
end
grid on;
xlabel('Position x (m)');
ylabel('Temperature T (deg C)');
title('Q5: 1-D Transient Temperature Distribution (Explicit FDM)');
legend('show', 'Location', 'northeast');
saveas(figFDM, fullfile('figures', 'Q5_FDM_Temperature_Distribution.png'));

Q5_results.x        = x;
Q5_results.tSteps   = tSteps;
Q5_results.Tmatrix  = Tmatrix;
Q5_results.r        = r;
Q5_results.dt_max   = dt_max;
Q5_results.targetTimes = targetTimes;
Q5_results.snapIdx  = snapIdx;
Q5_results.diffFinal = diffFinal;
save(fullfile('results', 'Q5_results.mat'), 'Q5_results');

fprintf('\n[Q5] Figure saved to figures/Q5_FDM_Temperature_Distribution.png\n');
fprintf('[Q5] Results saved to results/Q5_results.mat\n');
