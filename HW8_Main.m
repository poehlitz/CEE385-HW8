close all, clear all;

Intensity = load('Intensity.mat');

% Fragility parameters for post Nor
median = .05; 
sigma_ln = .31;

% Create Fragility
EDP = linspace(.001, .12, 1000);
P_C_EDP = normcdf((log(EDP) - log(median))/sigma_ln);

% Plot the fragility to check
figure
plot(EDP, P_C_EDP)
