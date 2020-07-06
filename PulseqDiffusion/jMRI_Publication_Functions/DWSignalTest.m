%% Test Example DWSignal for the following parameters: b_value=400 s*mm^-2, D=0.8*1e-3 mm^2 * s^-1
clear all

addpath(genpath('./SNR_CNR_Study'));

b_value = 400;
D       = 0.8*1e-3; % Diffusivity WM
DWSignal_pred = load('./Example_Data/DWSignal_TestResults.mat');
DWSignal_act  = DWSignal(b_value,D); % in ms
assertValue(DWSignal_act,DWSignal_pred.DWs_VectorTest(1))
clear b_value D

%% Test Example DWSignal for the following parameters: b_value=1 s*mm^-2, D=0.8*1e-3 mm^2 * s^-1
b_value = 1;
D       = 0.8*1e-3; % Diffusivity WM
DWSignal_act  = DWSignal(b_value,D); % in ms
assertValue(DWSignal_act,DWSignal_pred.DWs_VectorTest(2))
clear b_value D

%% Test Example DWSignal for the following parameters: b_value=1000 s*mm^-2, D=0.55*1e-3 mm^2 * s^-1
b_value = 1000;
D       = 0.55*1e-3; % Diffusivity LS
DWSignal_act  = DWSignal(b_value,D); % in ms
assertValue(DWSignal_act,DWSignal_pred.DWs_VectorTest(3))
clear b_value D

fprintf('All test for DWSignalTest run!')    
%%
function assertValue(actVal,expVal)
% Helper function to assert equality
% Takes two images and compares them 

tf = isequal(actVal,expVal);
assert(tf, 'TE values do not Match!');
end