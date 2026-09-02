function [Q, table_out] = simpson_one_third(t, q, h)
%SIMPSON_ONE_THIRD  Numerically integrate q(t) using Simpson's 1/3 Rule
%   (manual implementation).
%
%   [Q, table_out] = simpson_one_third(t, q, h)
%
%   Inputs:
%       t - vector of equally spaced time values (seconds), length n+1
%       q - vector of integrand values q(t) (W), length n+1
%       h - step size (seconds)
%
%   Outputs:
%       Q         - integral (Joules) via Simpson's 1/3 rule
%       table_out - [index, t, q, weight, weight*q] per row

    n = numel(q) - 1; % number of intervals

    % --- Input validation: Simpson's 1/3 requires an even number of intervals ---
    if mod(n, 2) ~= 0
        error('simpson_one_third:oddIntervals', ...
            'Simpson''s 1/3 Rule requires an even number of intervals (got n = %d).', n);
    end
    if numel(t) ~= numel(q)
        error('simpson_one_third:dimMismatch', 't and q must be the same length.');
    end

    weights = ones(1, n + 1);
    weights(2:2:n) = 4;   % odd-indexed (1-based: 2,4,6,...) -> coefficient 4
    weights(3:2:n) = 2;   % even-indexed (1-based: 3,5,7,...) -> coefficient 2

    Q = (h / 3) * sum(weights .* q);

    table_out = [(1:n+1)', t(:), q(:), weights(:), (weights(:) .* q(:))];
end
