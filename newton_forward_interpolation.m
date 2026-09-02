function [T_est, diffTable, p] = newton_forward_interpolation(t, T, h, t_query)
%NEWTON_FORWARD_INTERPOLATION  Estimate T at t_query using Newton's Forward
%   Interpolation formula (manual implementation, equal-interval data only).
%
%   [T_est, diffTable, p] = newton_forward_interpolation(t, T, h, t_query)
%
%   Inputs:
%       t       - vector of equally spaced independent-variable values (time, min)
%       T       - vector of dependent-variable values (temperature, deg C)
%       h       - step size (t(2) - t(1))
%       t_query - value of t at which to estimate T
%
%   Outputs:
%       T_est     - interpolated temperature at t_query
%       diffTable - forward difference table, padded with NaN (rows = data
%                   points, columns = difference order 0..n-1)
%       p         - the (t_query - t(1))/h value used in the formula

    n = numel(T);

    % --- Input validation ---
    if numel(t) ~= n
        error('newton_forward_interpolation:dimMismatch', 't and T must be the same length.');
    end
    spacing = diff(t);
    if any(abs(spacing - h) > 1e-9)
        error('newton_forward_interpolation:unequalSpacing', ...
            'Data must be equally spaced with step h = %.6f.', h);
    end

    % --- Build forward difference table ---
    diffTable = NaN(n, n);
    diffTable(:, 1) = T(:);
    for col = 2:n
        for row = 1:(n - col + 1)
            diffTable(row, col) = diffTable(row + 1, col - 1) - diffTable(row, col - 1);
        end
    end

    % --- Newton forward formula using the leading diagonal (row 1) ---
    p = (t_query - t(1)) / h;

    T_est = diffTable(1, 1);
    pProd = 1;
    for k = 1:(n - 1)
        pProd = pProd * (p - (k - 1));
        term = (pProd / factorial(k)) * diffTable(1, k + 1);
        T_est = T_est + term;
    end
end
