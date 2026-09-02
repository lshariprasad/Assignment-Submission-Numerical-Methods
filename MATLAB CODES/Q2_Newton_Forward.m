%% Q2_Newton_Forward.m
% CO2: Use equal-interval temperature observations to estimate a missing
% temperature using Newton's Forward Interpolation.

if ~exist('params', 'var')
    addpath(fullfile(pwd, 'functions'));
end

fprintf('\n============================================================\n');
fprintf(' Q2 - CO2: NEWTON FORWARD INTERPOLATION\n');
fprintf('============================================================\n');

t_data = [0 10 20 30 40];      % minutes
T_data = [42 47 51 54 56];     % deg C
h      = 10;                   % step size, minutes
t_query = 15;                  % point to estimate

[T_est, diffTable, p] = newton_forward_interpolation(t_data, T_data, h, t_query);

fprintf('\nForward Difference Table:\n');
fprintf('  t(min)   T       DeltaT   Delta2T   Delta3T   Delta4T\n');
n = numel(T_data);
for row = 1:n
    fprintf('  %5.1f  ', t_data(row));
    for col = 1:n
        val = diffTable(row, col);
        if isnan(val)
            fprintf('        ');
        else
            fprintf('%8.4f', val);
        end
    end
    fprintf('\n');
end

fprintf('\np = (t_query - t0)/h = (%g - %g)/%g = %.4f\n', t_query, t_data(1), h, p);
fprintf('Estimated T(%g min) = %.4f deg C\n', t_query, T_est);

% --- Validation: MATLAB built-in interp1 (spline / polynomial) reference ---
T_builtin = interp1(t_data, T_data, t_query, 'spline');
fprintf('\n[Validation] MATLAB built-in interp1 (spline) estimate = %.4f deg C\n', T_builtin);
fprintf('[Validation] Difference (manual vs built-in)            = %.6f deg C\n', abs(T_est - T_builtin));

% --- Interpolation graph ---
figInterp = figure('Name', 'Q2 Interpolation', 'Visible', 'off');
tt = linspace(min(t_data), max(t_data), 200);
TT = zeros(size(tt));
for i = 1:numel(tt)
    TT(i) = newton_forward_interpolation(t_data, T_data, h, tt(i));
end
plot(tt, TT, 'b-', 'LineWidth', 1.5); hold on;
plot(t_data, T_data, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
plot(t_query, T_est, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
grid on;
xlabel('Time t (min)');
ylabel('Temperature T (deg C)');
title('Q2: Newton Forward Interpolation Curve');
legend('Interpolating polynomial', 'Observed data', sprintf('Estimate at t=%g min', t_query), 'Location', 'southeast');
saveas(figInterp, fullfile('figures', 'Q2_Interpolation.png'));

Q2_results.t_data     = t_data;
Q2_results.T_data     = T_data;
Q2_results.diffTable  = diffTable;
Q2_results.p          = p;
Q2_results.t_query    = t_query;
Q2_results.T_est      = T_est;
Q2_results.T_builtin  = T_builtin;
save(fullfile('results', 'Q2_results.mat'), 'Q2_results');

fprintf('\n[Q2] Figure saved to figures/Q2_Interpolation.png\n');
fprintf('[Q2] Results saved to results/Q2_results.mat\n');
