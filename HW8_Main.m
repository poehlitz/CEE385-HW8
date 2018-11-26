close all, clear all;
%% Hazard Curve Read in 
%hazard = xlsread('HW6_Hazard.xlsx','A10:B1809')
load('HW8_hazardcurve.mat')

% Mean Annual Frequency for MCE Event
MAF = -log(1 - 0.02) / 50;
% Extrapolate to Intensity Measure at MCE Level
IM_MCE = interp1(hazard(:,2),hazard(:,1),MAF);
IM_DBE = 2/3*IM_MCE;
Sa = hazard(:,1);

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

for i = 1:numel(fields)
    m_temp(i) = med.(fields{i})(index);
    s_temp(i) = sigmaln_t.(fields{i})(index);  
end
    
MCE_median = interp1(Stripe, m_temp, IM_MCE); 
MCE_sig = interp1(Stripe, s_temp, IM_MCE);

PDF_MCE = (1./(EDP.*MCE_sig.*sqrt(2*pi))).*exp(-((log(EDP)-log(MCE_median)).^2)./(2*(MCE_sig^2)));

% Collapse Fragility Calculation
handles.numberCollapse= [0,0,0,0,1,2,4,3];
n = 7;
output = CollapseFragility(Stripe, n, handles, Sa);
param_collapse = output.parameters;

figure 
plot(EDP, PDF_MCE)
title('EDP/PDF Check MCE')

%Probability of Fracture Given No Collapse
P_fracture_MCE_NC = trapz(EDP, PDF_MCE.*P_C_EDP);
P_fracture_C = 1; %Assume if building collapses the connnection fractures

%Calculate Probability of Collapse Using Fragility Function
P_C_MCE = normcdf((log(IM_MCE)-log(param_collapse(1)))/param_collapse(2));
P_NC_MCE = 1 - P_C_MCE;

%Probability of Fracture Given MCE
P_fracture_MCE = P_fracture_MCE_NC*P_NC_MCE + P_fracture_C*P_C_MCE;

disp('Fracture Probability of Pre-North Connection under MCE: ')
disp(P_fracture_MCE)

%% A2
DBE_median = interp1(Stripe, m_temp, IM_DBE); 
DBE_sig = interp1(Stripe, s_temp, IM_DBE);

PDF_DBE = (1./(EDP.*DBE_sig.*sqrt(2*pi))).*exp(-((log(EDP)-log(DBE_median)).^2)./(2*(DBE_sig^2)));
figure 
plot(EDP, PDF_DBE)
title('EDP/PDF Check DBE')

%Probabilty of Fracture Given No Collapse
P_fracture_DBE_NC = trapz(EDP, PDF_DBE.*P_C_EDP);

%Probability of Collapse from Fragility Function
P_C_DBE = normcdf((log(IM_DBE)-log(param_collapse(1)))/param_collapse(2));
P_NC_DBE = 1 - P_C_DBE;

%Probability of Fracture Given DBE
P_fracture_DBE = P_fracture_DBE_NC*P_NC_DBE + P_fracture_C*P_C_DBE;

disp('Fracture Probability of Pre-North Connection under DBE: ')
disp(P_fracture_DBE)


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

P_fracture_MCE_NC = trapz(EDP, PDF_MCE.*P_C_EDP_prenorth);

%Probability of Fracture Given MCE
P_fracture_MCE = P_fracture_MCE_NC*P_NC_MCE + P_fracture_C*P_C_MCE;

disp('Fracture Probability of Connection under MCE: ')
disp(P_fracture_MCE)

P_fracture_DBE_NC = trapz(EDP, PDF_DBE.*P_C_EDP_prenorth);

%Probability of Fracture Given DBE
P_fracture_DBE = P_fracture_DBE_NC*P_NC_DBE + P_fracture_C*P_C_DBE;

disp('Fracture Probability of Connection under DBE: ')
disp(P_fracture_DBE)

%% B1
onesided = [.03, .04, .05; .3, .3, .3];
DS = zeros(length(EDP), size(onesided, 2)+1);

for i = 1:size(onesided, 2)
    DS(:, i) = normcdf((log(EDP)-log(onesided(1,i)))/(onesided(2, i)));
end

PDS = zeros(length(EDP), size(onesided, 2)+1);
PDS(:,1) = 1-DS(:,1);

for i = 1:size(onesided, 2)
    PDS(:,i+1) = DS(:,i)-DS(:,i+1);
end

figure
hold on
for i = 1:4
    plot(EDP, PDS(:,i))
end
legend('No Damage', 'DS1','DS2','DS3')
title('Probability of Being in certain Damage State')

EL_DS_1S = [0, 17400, 29300, 29300];
EL_DS_2S = [0, 30000, 61000, 52300];

% Expected Loss Given EDP for 1 sided connection or 2 sided conenction
EL_EDP_1S = PDS(:,2)*EL_DS_1S(2) + PDS(:,3)*EL_DS_1S(3) + PDS(:,4)*EL_DS_1S(4);
EL_EDP_2S = PDS(:,2)*EL_DS_2S(2) + PDS(:,3)*EL_DS_2S(3) + PDS(:,4)*EL_DS_2S(4);

numcxn_story = [2 4;
                2 4;
                4 4;
                4 4;
                4 4];

EL_EDP_Story = EL_EDP_1S'.*numcxn_story(:,1) + EL_EDP_2S'.*numcxn_story(:,2);

% Do stuff 
numfloors = 5;
PDF_MCE_story = zeros(numfloors, length(EDP));
for j = 1:numfloors
    for i = 1:numel(fields)
        m_temp(i) = med.(fields{i})(j);
        s_temp(i) = sigmaln_t.(fields{i})(j);  
    end
    
MCE_median = interp1(Stripe, m_temp, IM_MCE); 
MCE_sig = interp1(Stripe, s_temp, IM_MCE);

PDF_MCE_story(j, :) = (1./(EDP.*MCE_sig.*sqrt(2*pi))).*exp(-((log(EDP)-log(MCE_median)).^2)./(2*(MCE_sig^2)));
end

figure
hold on
for i = 1:numfloors
    plot(EDP, PDF_MCE_story(i,:))
end
legend('Story 5', 'Story 4', 'Story 3', 'Story 2', 'Story 1')

LStory_IM = zeros(1,numfloors);
for i = 1:numfloors
    LStory_IM(i) = trapz(EDP, PDF_MCE_story(i,:).*EL_EDP_Story(i,:));
end

L_IM_NC = sum(LStory_IM);

L_C_story = 45000*numcxn_story(:,1) + 70000*numcxn_story(:,2);
L_C = sum(L_C_story);

P_C_MCE = normcdf((log(IM_MCE)-log(param_collapse(1)))/param_collapse(2));
P_NC_MCE = 1 - P_C_MCE;

L_MCE = L_IM_NC*P_NC_MCE + L_C*P_C_MCE;
disp('Expected Loss MCE')
disp(L_MCE)






