%% Test Example TE_time for the following parameters: b_value=400 s*mm^-2, gmMax=0.016 T/m, gamma=42.577e6 Hz/T
clear all

addpath(genpath('./SNR_CNR_Study'));

b_value = 400*1e6;
gmMax   = 0.016;
gamma   = 42.577e6;
TEtime_pred = load('./Example_Data/TE_time_TestResults.mat'); % in ms
TEtime_act  = TE_time(b_value,gmMax,gamma); % in ms
assertValue(TEtime_act,TEtime_pred.TE_VectorTest(1))
clear b_value gmMax gamma

%% Test Example TE_time for the following parameters: b_value=1 s*mm^-2,gmMax=0.26 T/m,gamma=42.577e6 Hz/T
b_value = 1*1e6;
gmMax   = 0.26;
gamma   = 42.577e6;

% TEtime_pred = ; % in ms
TEtime_act  = TE_time(b_value,gmMax,gamma); % in ms
assertValue(TEtime_act,TEtime_pred.TE_VectorTest(2))
clear b_value gmMax gamma

%% Test Example TE_time for the following parameters: b_value=800 s*mm^-2,gmMax=0.22 T/m,gamma=42.577e6 Hz/T
b_value = 800*1e6;
gmMax   = 0.22;
gamma   = 42.577e6;

% TEtime_pred = ; % in ms
TEtime_act  = TE_time(b_value,gmMax,gamma); % in ms
assertValue(TEtime_act,TEtime_pred.TE_VectorTest(3))

clear b_value gmMax gamma

%%

fprintf('All test for TE_timeTest run!')    
%%
function assertValue(actVal,expVal)
% Helper function to assert equality
% Takes two images and compares them 

tf = isequal(actVal,expVal);
assert(tf, 'TE values do not Match!');
end