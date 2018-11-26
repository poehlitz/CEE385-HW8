close all, clear all;
%% Hazard Curve Read in 
%hazard = xlsread('HW6_Hazard.xlsx','A10:B1809')
load('HW8_hazardcurve.mat')

% Mean Annual Frequency for MCE Event
MAF = -log(1 - 0.02) / 50;
% Extrapolate to Intensity Measure at MCE Level
IM_MCE = interp1(hazard(:,2),hazard(:,1),MAF);
IM_DBE = 2/3*IM_MCE;

%% Fragility Read in 
Intensity = load('Intensity.mat');
% Stripes
Stripe = [.0918, .198, .305, .411, .518, .625, .731, .838];

% Fragility parameters for post Northridge
med_fragility = .05; 
sigma_ln_fragility = .31;

% Create Fragility
EDP = linspace(.001, .12, 1000);
P_C_EDP = normcdf((log(EDP) - log(med_fragility))/sigma_ln_fragility);

% Plot the fragility to check
figure
plot(EDP, P_C_EDP)
title('1A Fragility')

%% A1
% Fields variable for ability to loop over all fields in structure
fields = fieldnames(Intensity.Intensity);
Int = Intensity.Intensity;

% Loop over all IL's (fields) and find the median at each story
for i=1:numel(fields)
    data = Int.(fields{i});
    for j = 1:size(data,1)
        med.(fields{i})(j,1) = median(data(j,:));
        sigmaln_t.(fields{i})(j,1) = std(log(data(j,:)));
    end
end

% 2nd story median, sigma for MCE event
index = 4; % For the second story
vector = zeros(1, length(Stripe));

for i = 1:numel(fields)
    m_temp(i) = med.(fields{i})(index);
    s_temp(i) = sigmaln_t.(fields{i})(index);  
end
    
MCE_median = interp1(Stripe, m_temp, IM_MCE); 
MCE_sig = interp1(Stripe, s_temp, IM_MCE);

PDF_MCE = (1./(EDP.*MCE_sig.*sqrt(2*pi))).*exp(-((log(EDP)-log(MCE_median)).^2)./(2*(MCE_sig^2)));

figure 
plot(EDP, PDF_MCE)
title('EDP/PDF Check MCE')

P_fracture_MCE = trapz(EDP, PDF_MCE.*P_C_EDP);
disp('Fracture Probability of Connection under MCE: ')
disp(P_fracture_MCE)

%% A2
DBE_median = interp1(Stripe, m_temp, IM_DBE); 
DBE_sig = interp1(Stripe, s_temp, IM_DBE);

PDF_DBE = (1./(EDP.*DBE_sig.*sqrt(2*pi))).*exp(-((log(EDP)-log(DBE_median)).^2)./(2*(DBE_sig^2)));
figure 
plot(EDP, PDF_DBE)
title('EDP/PDF Check DBE')

P_fracture_DBE = trapz(EDP, PDF_DBE.*P_C_EDP);

%% A3 
% Fragility parameters for post Northridge
med_fragility = .0181; 
sigma_ln_fragility = .44;

% Create Fragility
EDP = linspace(.001, .12, 1000);
P_C_EDP_prenorth = normcdf((log(EDP) - log(med_fragility))/sigma_ln_fragility);

% Plot the fragility to check
figure
plot(EDP, P_C_EDP)

P_fracture_MCE = trapz(EDP, PDF_MCE.*P_C_EDP_prenorth);

disp('Fracture Probability of Connection under MCE: ')
disp(P_fracture_MCE)

P_fracture_DBE = trapz(EDP, PDF_DBE.*P_C_EDP_prenorth);

disp('Fracture Probability of Connection under DBE: ')
disp(P_fracture_DBE)

%% B1




    



