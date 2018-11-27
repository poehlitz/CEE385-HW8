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
        med.(fields{i})(j,1) = geomean(data(j,:));
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
%Enter Fragility Functions for One Sided Connections
%Each Column is a damage state
onesided = [.03, .04, .05; .3, .3, .3];
DS = zeros(length(EDP), size(onesided, 2)+1);

%Probability of Exceeding Damage State
for i = 1:size(onesided, 2)
    DS(:, i) = normcdf((log(EDP)-log(onesided(1,i)))/(onesided(2, i)));
end

%Probabiilty of Being Within Damage State
PDS = zeros(length(EDP), size(onesided, 2)+1);
PDS(:,1) = 1-DS(:,1);

for i = 1:size(onesided, 2)
    PDS(:,i+1) = DS(:,i)-DS(:,i+1);
end

%Plot Probability of Being Within Damage States
figure
hold on
for i = 1:4
    plot(EDP, PDS(:,i))
end
legend('No Damage', 'DS1','DS2','DS3')
title('Probability of Being in certain Damage State')

%Enter median of losses given a damage state
m_DS_1S = [0, 17400, 29300, 29300];
m_DS_2S = [0, 30000, 61000, 52300];

%Enter dispersion of expected losses given a damage state
s_DS_1S = [0, 0.35, 0.31, 0.31];
s_DS_2S = [0, 0.28, 0.28, 0.28];

%Calculate mean(expected) loss given a damage state
EL_DS_1S = m_DS_1S.*exp(0.5*s_DS_1S.^2);
EL_DS_2S = m_DS_2S.*exp(0.5*s_DS_2S.^2);

% Expected Loss Given EDP for 1 sided connection or 2 sided conenction
EL_EDP_1S = PDS(:,2)*EL_DS_1S(2) + PDS(:,3)*EL_DS_1S(3) + PDS(:,4)*EL_DS_1S(4);
EL_EDP_2S = PDS(:,2)*EL_DS_2S(2) + PDS(:,3)*EL_DS_2S(3) + PDS(:,4)*EL_DS_2S(4);

numcxn_story = [2 4;
                2 4;
                4 4;
                4 4;
                4 4];

EL_EDP_Story = EL_EDP_1S'.*numcxn_story(:,1) + EL_EDP_2S'.*numcxn_story(:,2);

% Find Median and Dispersion of EDP for each story and each stripe
numfloors = 5;
PDF_MCE_story = zeros(numfloors, length(EDP));
for j = 1:numfloors %loop over stories
    for i = 1:numel(fields) %loop over stripes
        m_story_stripe(j,i) = med.(fields{i})(j); % take median of values from each GM
        s_story_stripe(j,i) = sigmaln_t.(fields{i})(j); %standard deviations of values from each GM
    end

    %Linearly interpolate for value at exact intensity performing analysis at
    MCE_median = interp1(Stripe, m_story_stripe(j,:), IM_MCE);
    MCE_sig = interp1(Stripe, s_story_stripe(j,:), IM_MCE);

    %Calculation the probability of observing each EDP value at each story
    %given MCE
    PDF_MCE_story(j, :) = (1./(EDP.*MCE_sig.*sqrt(2*pi))).*exp(-((log(EDP)-log(MCE_median)).^2)./(2*(MCE_sig^2)));
end

% Plot EDP v. PDF(EDP given MCE)
figure
hold on
grid on
for i = 1:numfloors
    plot(EDP, PDF_MCE_story(i,:))
end
legend('Story 5', 'Story 4', 'Story 3', 'Story 2', 'Story 1')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on');


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

%% B2 - Function of Loss v. IM

% Interpolate Median Values from Stripes to all IM values
m_Story_IM = zeros(numfloors,size(hazard,1));
for i = 1:numfloors
    clear a b 
    a = interp1([0,Stripe], [0,m_story_stripe(i,:)], hazard(:,1));
    b = a(~isnan(a));
    %a(isnan(a)) = b(end);
    a(isnan(a)) = b(end)/Stripe(end)*hazard(isnan(a),1);
    m_Story_IM (i,:) = a; 
end

% Interpolate Dispersion Values from Stripes to all IM values
s_Story_IM = zeros(numfloors,size(hazard,1));
for i = 1:numfloors
    clear a b
    a = interp1([0,Stripe], [0,s_story_stripe(i,:)], hazard(:,1));
    b = a(~isnan(a));
    a(isnan(a)) = b(end)/Stripe(end)*hazard(isnan(a),1);
    s_Story_IM (i,:) = a; 
end

%Find Probability of Each EDP at Each IM Level Per Story
PDF_story = zeros(numfloors, length(EDP), size(hazard,1));
for i = 1:numfloors % Loop through every story
    for j= 1:length(EDP) %Loop through every EDP
    PDF_story(i, j, :) = (1./(EDP(j).*s_Story_IM(i,:).*sqrt(2*pi))).*exp(-((log(EDP(j))-log(m_Story_IM(i,:))).^2)./(2*(s_Story_IM(i,:).^2)));
    end
end

figure
plot(hazard(:,1),(1./(EDP(140).*s_Story_IM(1,:).*sqrt(2*pi))).*exp(-((log(EDP(140))-log(m_Story_IM(1,:))).^2)./(2*(s_Story_IM(1,:).^2))))
title('PDF for EDP = 0.02')
xlabel('Sa (g)')
ylabel('Probability')
grid on
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on');

for i=1:size(hazard,1)
    for j = 1:numfloors
    L_IM_story_NC (i,j) = trapz(EDP,PDF_story(j,:,i).*EL_EDP_Story(j,:));
    end
end

L_IM_NC = sum(L_IM_story_NC,2);

%Probability of Collapse at Each IM
P_C_IM = normcdf((log(hazard(:,1))-log(param_collapse(1)))/param_collapse(2));
P_NC_IM = 1 - P_C_IM ; 

%Expected Loss at Each IM
L_IM = L_C*P_C_IM + P_NC_IM.*L_IM_NC;

figure
plot(hazard(:,1),L_IM,hazard(:,1),P_NC_IM.*L_IM_NC,hazard(:,1),L_C*P_C_IM)
legend('Total','No Collapse','Collapse')
title('Variation of Loss Over Intensity Measures')
xlabel('Sa (g)')
ylabel('Cost ($)')
grid on
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on');

figure
plot(hazard(:,1),L_IM_NC)
title('Total Loss Given No Collapse')
xlabel('Sa (g)')
ylabel('Cost ($)')
grid on
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on');

% figure
% plot(Stripe,m_story_stripe)
% title('Median EDP Values at Stripes')
% xlabel('Sa (g)')
% ylabel('Median IDR')
% Grid on
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on');

figure
plot(hazard(:,1),m_Story_IM)
title('Median EDP Values')
xlabel('Sa (g)')
ylabel('Median IDR')
legend('Story 5','Story 4', 'Story 3', 'Story 2', 'Story 1')
grid on
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on');

figure
plot(hazard(:,1),s_Story_IM)
title('Standard Deviation Median EDP Values')
xlabel('Sa (g)')
ylabel('Median IDR')
legend('Story 5','Story 4', 'Story 3', 'Story 2', 'Story 1')
grid on
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on');

deriv = [diff(hazard(:,2))/(hazard(1,1)-hazard(2,1));0];
EAL = trapz(hazard(:,1),L_IM.*deriv);