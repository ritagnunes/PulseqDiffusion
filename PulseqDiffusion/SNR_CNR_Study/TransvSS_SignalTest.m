%% Test Example delta for the following parameters: FlipAngle=90 degrees, TRep=17.3011 s, TE=0.1739 s, T1=0.9982 s, T2=0.0727 s
clear all

addpath(genpath('./SNR_CNR_Study'));

FlipAngle = 90;
TRep      = 17.3011;
TE        = 0.1739;
T1        = 0.9982; 
T2        = 0.0727; 
TransvSS_signal_pred = load('./Example_Data/TransvSS_Signal_TestResults.mat'); % in s
TransvSS_signal_act  = TransvSS_Signal(FlipAngle,TRep,TE,T1,T1);
assertValue(TransvSS_signal_act,TransvSS_signal_pred.TransvSSsignal_VectorTest(1))
clear FlipAngle TRep TE T1 T2

%% Test Example delta for the following parameters: FlipAngle=81 degrees, TRep=3.2833 s, TE=0.0143 s, T1=0.1004 s, T2=0.0639 s

FlipAngle = 81;
TRep      = 3.2833;
TE        = 0.0143;
T1        = 0.1004; 
T2        = 0.0639; 
TransvSS_signal_pred = load('./Example_Data/TransvSS_Signal_TestResults.mat'); % in s
TransvSS_signal_act  = TransvSS_Signal(FlipAngle,TRep,TE,T1,T1);
assertValue(TransvSS_signal_act,TransvSS_signal_pred.TransvSSsignal_VectorTest(2))
clear FlipAngle TRep TE T1 T2


%% Test Example delta for the following parameters: FlipAngle=90 degrees, TRep=2.9760 s, TE=0.0132 s, T1=0.1004 s, T2=0.0639 s

FlipAngle = 90;
TRep      = 2.9760;
TE        = 0.0132;
T1        = 0.1004; 
T2        = 0.0639; 
TransvSS_signal_pred = load('./Example_Data/TransvSS_Signal_TestResults.mat'); % in s
TransvSS_signal_act  = TransvSS_Signal(FlipAngle,TRep,TE,T1,T1);
assertValue(TransvSS_signal_act,TransvSS_signal_pred.TransvSSsignal_VectorTest(3))
clear FlipAngle TRep TE T1 T2

fprintf('All test for TransvSS_SignalTest run!')    
%%
function assertValue(actVal,expVal)
% Helper function to assert equality
% Takes two images and compares them 

tf = isequal(actVal,expVal);
assert(tf, 'TransvSS_signal values do not Match!');
end
