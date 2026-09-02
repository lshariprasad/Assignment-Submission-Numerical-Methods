%% validation.m
% Q6 (part): Integrate and validate all five numerical methods.
% Assumes Q1-Q5 scripts have already been run in this session (main.m
% runs them first) so that Q1_results ... Q5_results structures exist.

fprintf('\n============================================================\n');
fprintf(' VALIDATION SUMMARY - ALL METHODS\n');
fprintf('============================================================\n');

validationTable = {};

% Q1 validation: residual should be ~0, and cross-check with fzero
f_check = @(T) params.eta*params.G*params.A - params.U*params.A*(T-params.Ta) ...
             - params.eps*params.sigma*params.A*((T+273.15)^4 - (params.Ta+273.15)^4);
T_fzero = fzero(f_check, 60);
fprintf('\n[Q1] Newton root                 = %.8f deg C\n', Q1_results.T_root);
fprintf('[Q1] MATLAB fzero() reference     = %.8f deg C\n', T_fzero);
fprintf('[Q1] Difference                   = %.3e deg C\n', abs(Q1_results.T_root - T_fzero));
fprintf('[Q1] Residual at root              = %.3e\n', Q1_results.residual);
validationTable(end+1,:) = {'Q1 Newton', Q1_results.T_root, T_fzero, abs(Q1_results.T_root - T_fzero)};

% Q2 validation: manual vs built-in interp1
fprintf('\n[Q2] Manual Newton-forward estimate T(15) = %.4f deg C\n', Q2_results.T_est);
fprintf('[Q2] MATLAB interp1 (spline) reference     = %.4f deg C\n', Q2_results.T_builtin);
fprintf('[Q2] Difference                            = %.4f deg C\n', abs(Q2_results.T_est - Q2_results.T_builtin));
validationTable(end+1,:) = {'Q2 Interp', Q2_results.T_est, Q2_results.T_builtin, abs(Q2_results.T_est - Q2_results.T_builtin)};

% Q3 validation: Simpson vs trapz
fprintf('\n[Q3] Simpson 1/3 integral   = %.2f J\n', Q3_results.Q_joules);
fprintf('[Q3] MATLAB trapz reference = %.2f J\n', Q3_results.Q_trapz);
fprintf('[Q3] Percent difference     = %.4f %%\n', 100*abs(Q3_results.Q_joules-Q3_results.Q_trapz)/Q3_results.Q_joules);
validationTable(end+1,:) = {'Q3 Simpson', Q3_results.Q_joules, Q3_results.Q_trapz, abs(Q3_results.Q_joules-Q3_results.Q_trapz)};

% Q4 validation: RK4 vs analytical (already computed in Q4 script)
fprintf('\n[Q4] RK4 final T(3600s)        = %.5f deg C\n', Q4_results.TOut(end));
fprintf('[Q4] Analytical final T(3600s) = %.5f deg C\n', Q4_results.Tanalytical(end));
fprintf('[Q4] Max absolute error        = %.3e deg C\n', max(Q4_results.absErr));
fprintf('[Q4] Max percent error         = %.6f %%\n', max(Q4_results.pctErr));
validationTable(end+1,:) = {'Q4 RK4', Q4_results.TOut(end), Q4_results.Tanalytical(end), max(Q4_results.absErr)};

% Q5 validation: stability + grid refinement check
fprintf('\n[Q5] Stability number r    = %.4f  (<= 0.5 required: %s)\n', Q5_results.r, mat2str(Q5_results.r <= 0.5));
fprintf('[Q5] dt_max                = %.4f s (used dt = 400 s)\n', Q5_results.dt_max);
fprintf('[Q5] Max diff vs finer dt/2 grid at final time = %.6f deg C\n', max(Q5_results.diffFinal));
validationTable(end+1,:) = {'Q5 FDM stability r', Q5_results.r, 0.5, NaN};

fprintf('\n============================================================\n');
fprintf(' Final Validation Table\n');
fprintf('============================================================\n');
fprintf('%-20s %16s %16s %16s\n', 'Method', 'Result', 'Reference', 'Abs Diff');
for i = 1:size(validationTable,1)
    row = validationTable(i,:);
    fprintf('%-20s %16.6f %16.6f', row{1}, row{2}, row{3});
    if isnan(row{4})
        fprintf('%16s\n', 'n/a');
    else
        fprintf('%16.3e\n', row{4});
    end
end

save(fullfile('results', 'validation_results.mat'), 'validationTable');
fprintf('\n[Validation] Results saved to results/validation_results.mat\n');
