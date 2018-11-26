clear all, close all
%hazard = xlsread('HW6_Hazard.xlsx','A10:B1809')
load('HW8_hazardcurve.mat')

% Mean Annual Frequency for MCE Event
MAF = -log(1 - 0.02) / 50;
% Extrapolate to Intensity Measure at MCE Level
IM_MCE = interp1(hazard(:,2),hazard(:,1),MAF);