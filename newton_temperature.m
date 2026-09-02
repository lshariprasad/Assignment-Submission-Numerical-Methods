function [T_root, iterTable, iterCount, converged] = newton_temperature(f, fprime, T0, tol, maxIter)
%NEWTON_TEMPERATURE  Solve f(T) = 0 using Newton's Method (manual implementation).
%
%   [T_root, iterTable, iterCount, converged] = newton_temperature(f, fprime, T0, tol, maxIter)
%
%   Inputs:
%       f       - function handle, f(T)
%       fprime  - function handle, f'(T) (analytical derivative)
%       T0      - initial guess (deg C)
%       tol     - stopping tolerance on |T_(n+1) - T_n|
%       maxIter - maximum number of iterations allowed
%
%   Outputs:
%       T_root    - converged root (operating temperature, deg C)
%       iterTable - [iter, Tn, f(Tn), fprime(Tn), Tn+1, abs error] per row
%       iterCount - number of iterations actually performed
%       converged - logical flag, true if tolerance was met

    if nargin < 5
        maxIter = 100;
    end

    iterTable = zeros(maxIter, 6);
    T = T0;
    converged = false;

    for k = 1:maxIter
        fT  = f(T);
        fpT = fprime(T);

        % --- Input validation: derivative cannot be zero ---
        if abs(fpT) < 1e-12
            error('newton_temperature:zeroDerivative', ...
                'Newton derivative is (near) zero at T = %.6f. Method fails to converge.', T);
        end

        Tnext = T - fT / fpT;
        err   = abs(Tnext - T);

        iterTable(k, :) = [k, T, fT, fpT, Tnext, err];

        T = Tnext;

        if err < tol
            converged = true;
            iterCount = k;
            iterTable = iterTable(1:k, :);
            T_root = T;
            return;
        end
    end

    % Did not converge within maxIter
    iterCount = maxIter;
    iterTable = iterTable(1:maxIter, :);
    T_root = T;
    warning('newton_temperature:notConverged', ...
        'Newton''s method did not converge within %d iterations.', maxIter);
end
