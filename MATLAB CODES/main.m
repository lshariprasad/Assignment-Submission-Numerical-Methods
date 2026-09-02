%% main.m
% =========================================================================
%  UBA10 - NUMERICAL METHODS
%  Assignment: Numerical Simulation of a Solar Water-Heating System
%  CO1-CO5 | SDG 4, SDG 7, SDG 9
%
%  Master driver script. Running this file executes the complete project:
%    Q1 - Newton's Method            (CO1)
%    Q2 - Newton Forward Interpolation (CO2)
%    Q3 - Simpson's 1/3 Rule          (CO3)
%    Q4 - Fourth-order Runge-Kutta    (CO4)
%    Q5 - Explicit Finite Difference  (CO5)
%    Validation of all methods
%    Engineering analysis / CO summary
% =========================================================================

clear; clc; close all;

fprintf('#########################################################\n');
fprintf('#  SOLAR WATER-HEATING SYSTEM - NUMERICAL METHODS STUDY  #\n');
fprintf('#  UBA10 Numerical Methods | CO1-CO5 | SDG 4, 7, 9       #\n');
fprintf('#########################################################\n');

% --- Ensure output folders exist ---
if ~exist('figures', 'dir'); mkdir('figures'); end
if ~exist('results', 'dir'); mkdir('results'); end
addpath(fullfile(pwd, 'functions'));

% =========================================================================
% Representative Engineering Input Data and Assumptions
% (No faculty-supplied numerical data sheet was provided for this
%  assignment. All numerical values below are clearly labelled
%  representative engineering assumptions used for demonstration.)
% =========================================================================
params = struct();

% Q1 - Newton's Method parameters
params.A     = 2;                    % collector area, m^2
params.G     = 700;                  % solar irradiance, W/m^2
params.eta   = 0.70;                 % collector efficiency
params.U     = 6;                    % overall heat-loss coefficient, W/m^2K
params.Ta    = 30;                   % ambient temperature, deg C
params.eps   = 0.90;                 % emissivity
params.sigma = 5.670374419e-8;       % Stefan-Boltzmann constant, W/m^2K^4

% Q4 - RK4 lumped-capacitance parameters (reuse eta,G,A,U,Ta above)
params.m  = 50;                      % water mass, kg
params.cp = 4180;                    % specific heat of water, J/kgK

fprintf('\nShared parameters loaded:\n');
disp(params);

% =========================================================================
% Execute Q1 - Q5
% =========================================================================
Q1_Newton_Method;
Q2_Newton_Forward;
Q3_Simpson_OneThird;
Q4_RK4;
Q5_Finite_Difference;

% =========================================================================
% Validation
% =========================================================================
validation;

% =========================================================================
% Engineering Analysis and Final CO Summary
% =========================================================================
engineering_analysis;

% =========================================================================
% Final overall summary + save consolidated results
% =========================================================================
fprintf('\n#########################################################\n');
fprintf('#  PROJECT EXECUTION COMPLETE                             #\n');
fprintf('#########################################################\n');
fprintf('\nAll figures saved in:  %s\n', fullfile(pwd, 'figures'));
fprintf('All .mat results saved in: %s\n', fullfile(pwd, 'results'));

all_results = struct('Q1', Q1_results, 'Q2', Q2_results, 'Q3', Q3_results, ...
                      'Q4', Q4_results, 'Q5', Q5_results, 'params', params);
save(fullfile('results', 'all_results_summary.mat'), 'all_results');
fprintf('Consolidated summary saved to results/all_results_summary.mat\n');

fprintf('\nDone. See README.md for a description of every file and\n');
fprintf('expected output.\n');
