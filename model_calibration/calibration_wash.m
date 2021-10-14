clear all
close all

% HAITI -- parameter estimation

% Problem specific parameter settings
MCMCPar.n = 4;                         % Dimension of the problem (number of parameters to be estimated)
MCMCPar.ndraw = 15000000;                % Maximum number of function evaluations
MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates

% Recommended parameter settings
MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
MCMCPar.nCR = 3;                        % Number of crossover values used
MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
MCMCPar.steps = 10;                     % Number of steps before calculating convergence diagnostics

% -----------------------------------------------------------------------------------------------------------------------
Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
% -----------------------------------------------------------------------------------------------------------------------

% --------------------------------------- Added for reduced sample storage ----------------------------------------------
Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
Extra.T = 1000;                         % Every Tth sample is collected
% -----------------------------------------------------------------------------------------------------------------------

% What type of initial sampling
Extra.InitPopulation = 'LHS_BASED';

% Give the parameter ranges (minimum and maximum values)
    
ParRange.minn = [0  0    1    1]; 
ParRange.maxn = [5  5    4*365  4*365]; 
    
% Define the boundary handling
Extra.BoundHandling = 'Reflect';
% Save in memory or not
Extra.save_in_memory = 'Yes';


% Define data
cases_week = csvread('data/cases.csv',1,1)';    %ski the first 2 weeks %from 20 Oct 2010 to 31 Aug 2011 --> the last week is full, so actually to 3 Sep 2011

Measurement.MeasData = cases_week(:,1:350); 

Measurement.Sigma = []; 
Measurement.N = length(Measurement.MeasData(:));
% Define modelName
ModelName = 'm5c';
% Define likelihood function -- Sum of Squared Error
option = 3;
edit=1;

% Run the distributed DREAM algorithm with sampling from past
[Sequences,Reduced_Seq,X,Z,output] = dream_zs(MCMCPar,ParRange,Measurement,ModelName,Extra,option,edit);
