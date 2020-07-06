%% Test Example Ernstangle_d for the following parameters: TR=17.3011 s, T1=0.9982 s
clear all

addpath(genpath('./SNR_CNR_Study'));

TR = 17.3011;
T1 = 0.9982; 
Eangle_pred = load('./Example_Data/Ernstangle_d_TestResults.mat');
Eangle_act  = Ernstangle_d(TR,T1); % in degrees
assertValue(Eangle_act,Eangle_pred.Eangle_VectorTest(1))
clear TR T1

%% Test Example Ernstangle_d for the following parameters: TR=1.7642 s, T1=0.9982 s
addpath(genpath('./SNR_CNR_Study'));

TR = 1.7642;
T1 = 0.9982; 
Eangle_act  = Ernstangle_d(TR,T1); % in degrees
assertValue(Eangle_act,Eangle_pred.Eangle_VectorTest(2))
clear TR T1

%% Test Example Ernstangle_d for the following parameters: TR=1.2858 s, T1=0.6806 s
addpath(genpath('./SNR_CNR_Study'));

TR = 1.2858;
T1 = 0.6806; 
Eangle_act  = Ernstangle_d(TR,T1); % in degrees
assertValue(Eangle_act,Eangle_pred.Eangle_VectorTest(3))
clear TR T1

fprintf('All test for Ernstangle_dTest run!')    
%%
function assertValue(actVal,expVal)
% Helper function to assert equality
% Takes two images and compares them 

tf = isequal(actVal,expVal);
assert(tf, 'TE values do not Match!');
end