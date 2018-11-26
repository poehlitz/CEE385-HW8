function [handles] = CollapseFragility(stripes,n,handles)

%Hello

%Fit of Collapse Fragility Function based on Maximum Likelihood
fun2 = @(v) maxLikelihood(handles.numberCollapse, n, stripes, v(1), v(2));
v_guess = [.8, .4];
ML_minimumParameters = fminsearch(fun2, v_guess);

ML_Sa = handles.hazardDerivative(1,:);
ML_P = normcdf((log(ML_Sa)-log(ML_minimumParameters(1)))/ML_minimumParameters(2));

figure
plot(stripes,handles.numberCollapse/n,'o',ML_Sa,ML_P,'k')
grid on
title('Collapse Fragility Function')
legend('Stripe Analysis Median Collapse', 'Max Likelihood Fragility Fit')
xlabel('Sa (g)')
ylabel('P[C]')
legend('Location','northwest')
xlim([0,2])
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on');

%Calculate PDF of Collapse
pdf_collapse = normcdf((log(handles.hazardDerivative(1,:))-log(ML_minimumParameters(1)))./ML_minimumParameters(2)).*handles.hazardDerivative(2,:);
%Integrate PDF of Collapse to get Mean Annual Frequency
MAF_c = trapz(handles.hazardDerivative(1,:), pdf_collapse);
%Probability of Collapse in 50 years
Prob_50 = 1 - exp(-MAF_c*50);

handles.MAF_c = MAF_c;
handles.Prob_50 = Prob_50;
