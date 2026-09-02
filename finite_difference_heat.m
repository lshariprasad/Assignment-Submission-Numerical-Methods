function [x, tSteps, Tmatrix, r] = finite_difference_heat(L, Nx, alpha, dx, dt, totalTime, T_left, T_right, T_initial)
%FINITE_DIFFERENCE_HEAT  Explicit (forward-time, central-space) finite
%   difference solution of the 1-D transient heat equation
%   dT/dt = alpha * d2T/dx2 (manual implementation).
%
%   [x, tSteps, Tmatrix, r] = finite_difference_heat(L, Nx, alpha, dx, dt, ...
%                                   totalTime, T_left, T_right, T_initial)
%
%   Inputs:
%       L         - domain length (m)
%       Nx        - number of spatial nodes (including boundaries)
%       alpha     - thermal diffusivity (m^2/s)
%       dx        - spatial step (m)
%       dt        - time step (s)
%       totalTime - total simulation time (s)
%       T_left    - Dirichlet BC at x = 0 (deg C)
%       T_right   - Dirichlet BC at x = L (deg C)
%       T_initial - initial interior temperature (deg C)
%
%   Outputs:
%       x       - spatial node coordinates (m), 1 x Nx
%       tSteps  - time values at each stored step (s)
%       Tmatrix - temperature matrix, rows = time steps, cols = spatial nodes
%       r       - stability (Fourier) number r = alpha*dt/dx^2

    % --- Input validation ---
    if abs((Nx - 1) * dx - L) > 1e-9
        error('finite_difference_heat:gridMismatch', ...
            '(Nx-1)*dx must equal L. Got (Nx-1)*dx = %.6f, L = %.6f.', (Nx-1)*dx, L);
    end

    r = alpha * dt / dx^2;
    if r > 0.5
        error('finite_difference_heat:unstable', ...
            'Stability violated: r = %.6f > 0.5. Reduce dt or increase dx.', r);
    end

    x = linspace(0, L, Nx);
    nSteps = round(totalTime / dt);

    T = T_initial * ones(1, Nx);
    T(1)   = T_left;
    T(end) = T_right;

    Tmatrix = zeros(nSteps + 1, Nx);
    tSteps  = zeros(nSteps + 1, 1);
    Tmatrix(1, :) = T;
    tSteps(1) = 0;

    for n = 1:nSteps
        Tnew = T;
        for i = 2:(Nx - 1)
            Tnew(i) = T(i) + r * (T(i+1) - 2*T(i) + T(i-1));
        end
        Tnew(1)   = T_left;   % Dirichlet BC maintained every step
        Tnew(end) = T_right;

        T = Tnew;
        Tmatrix(n + 1, :) = T;
        tSteps(n + 1) = n * dt;
    end
end
