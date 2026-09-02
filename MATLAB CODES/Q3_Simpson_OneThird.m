%% Q3_Simpson_OneThird.m
% CO3: Calculate total heat supplied using Simpson's 1/3 Rule.

if ~exist('params', 'var')
    addpath(fullfile(pwd, 'functions'));
    params = struct();
    params.A   = 2;
    params.eta = 0.70;
end

fprintf('\n============================================================\n');
fprintf(' Q3 - CO3: SIMPSON''S 1/3 RULE (Total Heat Supplied)\n');
fprintf('============================================================\n');

A   = params.A;
eta = params.eta;

t_min = [0 10 20 30 40 50 60];             % minutes
G     = [600 700 800 900 850 750 650];     % W/m^2
h_min = 10;
h_s   = h_min * 60;                        % step size in seconds = 600 s
t_s   = t_min * 60;                        % time in seconds

q = eta .* A .* G;   % instantaneous heat rate, W

[Q_joules, tbl] = simpson_one_third(t_s, q, h_s);

fprintf('\nHeat-rate table:\n');
fprintf(' i   t(min)   t(s)     G(W/m^2)   q(t)=etaAG (W)   Weight   Weight*q\n');
for i = 1:size(tbl,1)
    fprintf(' %d   %5.0f   %6.0f   %8.1f   %14.2f   %6.0f   %10.2f\n', ...
        i-1, t_min(i), tbl(i,2), G(i), tbl(i,3), tbl(i,4), tbl(i,5));
end

Q_MJ  = Q_joules / 1e6;
Q_kWh = Q_joules / 3.6e6;

fprintf('\nn (intervals) = %d  (even -> valid for Simpson 1/3)\n', numel(q) - 1);
fprintf('h = %g s\n', h_s);
fprintf('Total heat supplied Q = %.2f J\n', Q_joules);
fprintf('                      = %.4f MJ\n', Q_MJ);
fprintf('                      = %.4f kWh\n', Q_kWh);

% --- Validation: MATLAB built-in trapz as an independent reference ---
Q_trapz = trapz(t_s, q);
fprintf('\n[Validation] MATLAB built-in trapz reference   = %.2f J (%.4f MJ)\n', Q_trapz, Q_trapz/1e6);
fprintf('[Validation] Simpson vs trapz absolute diff    = %.2f J\n', abs(Q_joules - Q_trapz));
fprintf('[Validation] Simpson vs trapz percent diff     = %.4f %%\n', 100*abs(Q_joules-Q_trapz)/Q_joules);

% --- Plot heat rate curve with area under curve context ---
figSimpson = figure('Name', 'Q3 Simpson Heat Rate', 'Visible', 'off');
plot(t_min, q, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on;
xlabel('Time (min)');
ylabel('Heat rate q(t) = \eta A G(t)  (W)');
title('Q3: Collector Heat-Rate Profile (Simpson''s 1/3 Rule Integration)');
saveas(figSimpson, fullfile('figures', 'Q3_HeatRate.png'));

Q3_results.t_min   = t_min;
Q3_results.G       = G;
Q3_results.q       = q;
Q3_results.table   = tbl;
Q3_results.Q_joules = Q_joules;
Q3_results.Q_MJ    = Q_MJ;
Q3_results.Q_kWh   = Q_kWh;
Q3_results.Q_trapz = Q_trapz;
save(fullfile('results', 'Q3_results.mat'), 'Q3_results');

fprintf('\n[Q3] Figure saved to figures/Q3_HeatRate.png\n');
fprintf('[Q3] Results saved to results/Q3_results.mat\n');
