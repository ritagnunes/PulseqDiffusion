%% Test Example aux_delta for the following parameters: b_value=10 s*mm^-2, gmMax=0.030 T/m, part_Tend=0.0458 s, gamma=42.577e6 Hz/T
clear all

addpath(genpath('./SNR_CNR_Study'));

b_value   = 10*1e6;
gmMax     = 0.030;
part_Tend = 0.0458;
gamma     = 42.577e6; 
aDelta_pred = load('./Example_Data/aux_delta_TestResults.mat');
aDelta_act  = aux_delta(b_value,gmMax,part_Tend,gamma);
assertValue(aDelta_act,aDelta_pred.aDelta_VectorTest(1))
clear b_value gmMax part_Tend gamma

%% Test Example aux_delta for the following parameters: b_value=1000 s*mm^-2, gmMax=0.014 T/m, part_Tend=0.0458 s, gamma=42.577e6 Hz/T

b_value   = 1000*1e6;
gmMax     = 0.014;
part_Tend = 0.0458;
gamma     = 42.577e6;
aDelta_pred = load('./Example_Data/aux_delta_TestResults.mat');
aDelta_act = aux_delta(b_value,gmMax,part_Tend,gamma);
assertValue(aDelta_act,aDelta_pred.aDelta_VectorTest(2))
clear b_value gmMax part_Tend gamma

%% Test Example aux_delta for the following parameters: b_value=200 s*mm^-2, gmMax=0.020 T/m, part_Tend=0.0458 s, gamma=42.577e6 Hz/T

b_value   = 200*1e6;
gmMax     = 0.020;
part_Tend = 0.0458;
gamma     = 42.577e6;
aDelta_pred = load('./Example_Data/aux_delta_TestResults.mat');
aDelta_act = aux_delta(b_value,gmMax,part_Tend,gamma);
assertValue(aDelta_act,aDelta_pred.aDelta_VectorTest(3))
clear b_value gmMax part_Tend gamma

fprintf('All test for aux_deltaTest run!')    
%%
function assertValue(actVal,expVal)
% Helper function to assert equality
% Takes two images and compares them 

tf = isequal(actVal,expVal);
assert(tf, 'TE values do not Match!');
end