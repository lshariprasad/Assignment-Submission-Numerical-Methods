function [tOut, TOut, rkTable] = rk4_solver(f, t0, T0, h, tfinal)
%RK4_SOLVER  Fourth-order Runge-Kutta ODE solver (manual implementation).
%
%   [tOut, TOut, rkTable] = rk4_solver(f, t0, T0, h, tfinal)
%
%   Inputs:
%       f      - function handle, dT/dt = f(t, T)
%       t0     - initial time (s)
%       T0     - initial temperature (deg C)
%       h      - step size (s)
%       tfinal - final time (s)
%
%   Outputs:
%       tOut    - vector of time points (s)
%       TOut    - vector of temperature values (deg C) at tOut
%       rkTable - [step, t_n, T_n, k1, k2, k3, k4, T_n+1] per row

    % --- Input validation ---
    if h <= 0
        error('rk4_solver:invalidStep', 'Step size h must be positive.');
    end
    if tfinal <= t0
        error('rk4_solver:invalidRange', 'tfinal must be greater than t0.');
    end

    nSteps = round((tfinal - t0) / h);
    tOut = zeros(nSteps + 1, 1);
    TOut = zeros(nSteps + 1, 1);
    rkTable = zeros(nSteps, 8);

    t = t0;
    T = T0;
    tOut(1) = t;
    TOut(1) = T;

    for i = 1:nSteps
        k1 = f(t, T);
        k2 = f(t + h/2, T + h*k1/2);
        k3 = f(t + h/2, T + h*k2/2);
        k4 = f(t + h,   T + h*k3);

        Tnext = T + (h/6) * (k1 + 2*k2 + 2*k3 + k4);

        rkTable(i, :) = [i, t, T, k1, k2, k3, k4, Tnext];

        t = t + h;
        T = Tnext;

        tOut(i + 1) = t;
        TOut(i + 1) = T;
    end
end
