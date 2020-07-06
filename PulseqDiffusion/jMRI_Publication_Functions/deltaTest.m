%% Test Example delta for the following parameters: aux_delta=aux_delta_TestResults(1), gmMax=0.030 T/m, part_Tend=0.0458 s, gamma=42.577e6 Hz/T
clear all

addpath(genpath('./SNR_CNR_Study'));

aDelta = load('./Example_Data/aux_delta_TestResults.mat'); % in s
aux_delta = aDelta.aDelta_VectorTest(1);
gmMax     = 0.030;
part_Tend = 0.0458;
gamma     = 42.577e6; 
Delta_pred = load('./Example_Data/delta_TestResults.mat'); % in s
Delta_act  = delta(gmMax,part_Tend,aux_delta,gamma);
assertValue(Delta_act,Delta_pred.Delta_VectorTest(1))
clear b_value gmMax part_Tend gamma

%% Test Example aux_delta for the following parameters: aux_delta=aux_delta_TestResults(2), gmMax=0.014 T/m, part_Tend=0.0458 s, gamma=42.577e6 Hz/T

aux_delta = aDelta.aDelta_VectorTest(2);
gmMax     = 0.014;
part_Tend = 0.0458;
gamma     = 42.577e6; 
Delta_act  = delta(gmMax,part_Tend,aux_delta,gamma);
assertValue(Delta_act,Delta_pred.Delta_VectorTest(2))
clear b_value gmMax part_Tend gamma

%% Test Example aux_delta for the following parameters: aux_delta=aux_delta_TestResults(3), gmMax=0.020 T/m, part_Tend=0.0458 s, gamma=42.577e6 Hz/T

aux_delta = aDelta.aDelta_VectorTest(3);
gmMax     = 0.020;
part_Tend = 0.0458;
gamma     = 42.577e6; 
Delta_act  = delta(gmMax,part_Tend,aux_delta,gamma);
assertValue(Delta_act,Delta_pred.Delta_VectorTest(3))
clear b_value gmMax part_Tend gamma

fprintf('All test for deltaTest run!')    
%%
function assertValue(actVal,expVal)
% Helper function to assert equality
% Takes two images and compares them 

tf = isequal(actVal,expVal);
assert(tf, 'TE values do not Match!');
end