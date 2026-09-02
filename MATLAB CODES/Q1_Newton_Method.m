%% Q1_NEWTON_METHOD
% Newton-Raphson Method
% Steady-state operating temperature of a solar collector

clear;
clc;
close all;

%% Create folders
if ~exist('figures','dir')
    mkdir('figures');
end

%% Engineering Parameters

A     = 2;                 % Collector area (m^2)
G     = 700;               % Solar irradiance (W/m^2)
eta   = 0.70;              % Collector efficiency
U     = 6;                 % Heat-loss coefficient (W/m^2 K)
Ta    = 30;                % Ambient temperature (deg C)
eps_r = 0.90;              % Emissivity
sigma = 5.670374419e-8;    % Stefan-Boltzmann constant

%% Newton-Raphson Parameters

T0  = 60;                  % Initial guess (deg C)
tol = 1e-6;                % Stopping tolerance

%% Call Newton-Raphson Function

[T_final, iterTable] = newton_temperature( ...
    T0, eta, G, A, U, Ta, eps_r, sigma, tol);

%% Display Results

fprintf('\n');
fprintf('============================================================\n');
fprintf(' Q1: NEWTON-RAPHSON METHOD\n');
fprintf(' Steady-State Operating Temperature of Solar Collector\n');
fprintf('============================================================\n\n');

fprintf('%4s %14s %16s %14s %14s %14s\n', ...
    'Iter', 'Tn (C)', 'f(Tn)', 'f''(Tn)', 'Tn+1 (C)', '|Error| (C)');

fprintf('--------------------------------------------------------------------------\n');

for k = 1:size(iterTable,1)

    err = abs(iterTable(k,5) - iterTable(k,2));

    fprintf('%4d %14.8f %16.6e %14.5f %14.8f %14.6e\n', ...
        iterTable(k,1), ...
        iterTable(k,2), ...
        iterTable(k,3), ...
        iterTable(k,4), ...
        iterTable(k,5), ...
        err);
end

fprintf('--------------------------------------------------------------------------\n');

fprintf('\nConverged steady-state operating temperature:\n');
fprintf('T = %.7f degC\n', T_final);

fprintf('\nNumber of iterations = %d\n', size(iterTable,1));

%% Calculate Errors

errs = abs(iterTable(:,5) - iterTable(:,2));

%% Plot Newton-Raphson Convergence

figure('Name','Q1 Newton Convergence');

semilogy(iterTable(:,1), errs, '-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6);

xlabel('Iteration Number');
ylabel('|Error| (^oC), log scale');

title('Newton-Raphson Method Convergence');

grid on;

%% Save Figure

saveas(gcf, ...
    fullfile('figures','Q1_Newton_Convergence.png'));

fprintf('\nConvergence graph saved in the "figures" folder.\n');


%% ============================================================
% LOCAL FUNCTION: NEWTON-RAPHSON METHOD
% =============================================================

function [T, iterTable] = newton_temperature( ...
    T0, eta, G, A, U, Ta, eps_r, sigma, tol)

    % Maximum number of iterations
    maxIter = 100;

    % Initial temperature
    T = T0;

    % Initialize iteration table
    iterTable = [];

    %% Newton-Raphson Iteration

    for k = 1:maxIter

        % Convert Celsius to Kelvin
        Tk = T + 273.15;
        Tak = Ta + 273.15;

        % ----------------------------------------------------
        % Energy balance equation
        %
        % Solar energy absorbed:
        %       eta*G*A
        %
        % Heat loss:
        %       U*A*(T-Ta)
        %
        % Radiation loss:
        %       eps*sigma*A*(Tk^4 - Tak^4)
        %
        % f(T) = absorbed energy - heat losses
        % ----------------------------------------------------

        f = eta*G*A ...
            - U*A*(T-Ta) ...
            - eps_r*sigma*A*(Tk^4 - Tak^4);

        % ----------------------------------------------------
        % Derivative of f(T)
        % ----------------------------------------------------

        fprime = -U*A ...
            - 4*eps_r*sigma*A*Tk^3;

        % ----------------------------------------------------
        % Newton-Raphson formula
        %
        % T(n+1) = T(n) - f(Tn)/f'(Tn)
        % ----------------------------------------------------

        Tnew = T - f/fprime;

        % Store iteration information
        iterTable(k,:) = [k, T, f, fprime, Tnew];

        % Calculate error
        error = abs(Tnew - T);

        % Check convergence
        if error < tol
            T = Tnew;
            return;
        end

        % Update temperature
        T = Tnew;

    end

    % If maximum iterations reached
    warning('Newton-Raphson method did not converge within %d iterations.', ...
        maxIter);

end
