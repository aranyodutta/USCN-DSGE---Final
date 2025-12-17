%% run_all.m
% Add Dynare to MATLAB path
addpath('C:\dynare\6.2\matlab')

% (Optional but recommended) Initialize Dynare paths
dynare_config;

% Clear previous Dynare variables
clearvars -except M_ oo_ options_

% Run the Dynare model
dynare final_model.mod
