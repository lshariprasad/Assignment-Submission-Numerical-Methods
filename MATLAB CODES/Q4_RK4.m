%% Q4_RK4.m
% CO4: Formulate collector temperature variation as a first-order ODE and
% solve using fourth-order Runge-Kutta. Plot temperature versus time.
%
% Governing ODE (lumped-capacitance energy balance):
%   m*cp*dT/dt = eta*G*A - U*A*(T - Ta)
%   =>  dT/dt = [eta*G*A - U*A*(T - Ta)] / (m*cp)

if ~exist('params', 'var')
    addpath(fullfile(pwd, 'functions'));
    params = struct();
    params.m   = 50;
    params.cp  = 4180;
    params.eta = 0.70;
    params.G   = 700;
    params.A   = 2;
    params.U   = 6;
    params.Ta  = 30;
end

fprintf('\n============================================================\n');
fprintf(' Q4 - CO4: FOURTH-ORDER RUNGE-KUTTA (RK4)\n');
fprintf('============================================================\n');

m = params.m; cp = params.cp; eta = params.eta;
G = params.G; A = params.A; U = params.U; Ta = params.Ta;

T0_rk4 = 35;      % initial temperature, deg C
h      = 600;     % step size, s
tfinal = 3600;     % final time, s

f_ode = @(t, T) (eta*G*A - U*A*(T - Ta)) / (m*cp);

[tOut, TOut, rkTable] = rk4_solver(f_ode, 0, T0_rk4, h, tfinal);

fprintf('\nRK4 Iteration Table:\n');
fprintf(' Step |  t_n (s)  |  T_n (C)   |     k1        |     k2        |     k3        |     k4        |  T_n+1 (C)\n');
fprintf('---------------------------------------------------------------------------------------------------------------\n');
for i = 1:size(rkTable, 1)
    fprintf(' %4d | %8.0f  | %9.5f  | %12.6e | %12.6e | %12.6e | %12.6e | %10.5f\n', ...
        rkTable(i,1), rkTable(i,2), rkTable(i,3), rkTable(i,4), rkTable(i,5), rkTable(i,6), rkTable(i,7), rkTable(i,8));
end

fprintf('\nFinal RK4 temperature at t = %d s : T = %.5f deg C\n', tfinal, TOut(end));

% --- Analytical reference solution (linear ODE, closed form) ---
Tinf = Ta + eta*G*A / (U*A);
tau  = m*cp / (U*A);
T_analytical = @(t) Tinf + (T0_rk4 - Tinf) .* exp(-t ./ tau);

Tanalytical_vec = T_analytical(tOut);
absErr = abs(TOut - Tanalytical_vec);
relErr = absErr ./ abs(Tanalytical_vec);
pctErr = relErr * 100;

fprintf('\nAnalytical model:  T_inf = %.6f deg C,  tau = %.4f s\n', Tinf, tau);
fprintf('\nRK4 vs Analytical Comparison:\n');
fprintf('   t(s)   |  RK4 T(C)  | Analytical T(C) | Abs Error | Rel Error   | %% Error\n');
for i = 1:numel(tOut)
    fprintf(' %7.0f  | %9.5f  | %14.5f  | %9.3e | %10.3e | %8.6f\n', ...
        tOut(i), TOut(i), Tanalytical_vec(i), absErr(i), relErr(i), pctErr(i));
end

% --- Plot: RK4 vs Analytical ---
figRK4 = figure('Name', 'Q4 RK4 vs Analytical', 'Visible', 'off');
plot(tOut, TOut, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
tt_fine = linspace(0, tfinal, 200);
plot(tt_fine, T_analytical(tt_fine), 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time t (s)');
ylabel('Collector Temperature T (deg C)');
title('Q4: RK4 Numerical Solution vs Analytical Solution');
legend('RK4 solution', 'Analytical solution', 'Location', 'southeast');
saveas(figRK4, fullfile('figures', 'Q4_RK4_vs_Analytical.png'));

% --- Plot: error vs time ---
figErr = figure('Name', 'Q4 RK4 Error', 'Visible', 'off');
plot(tOut, absErr, 'm-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
grid on;
xlabel('Time t (s)');
ylabel('Absolute Error |RK4 - Analytical| (deg C)');
title('Q4: RK4 Absolute Error vs Time');
saveas(figErr, fullfile('figures', 'Q4_RK4_Error.png'));

Q4_results.tOut   = tOut;
Q4_results.TOut   = TOut;
Q4_results.rkTable = rkTable;
Q4_results.Tinf   = Tinf;
Q4_results.tau    = tau;
Q4_results.Tanalytical = Tanalytical_vec;
Q4_results.absErr = absErr;
Q4_results.relErr = relErr;
Q4_results.pctErr = pctErr;
save(fullfile('results', 'Q4_results.mat'), 'Q4_results');

fprintf('\n[Q4] Figures saved to figures/Q4_RK4_vs_Analytical.png and figures/Q4_RK4_Error.png\n');
fprintf('[Q4] Results saved to results/Q4_results.mat\n');
