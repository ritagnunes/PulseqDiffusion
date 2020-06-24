%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior T�cnico - Universidade de Lisboa
%
% name: LongSS_Signal
% type: Function
% description: Longitudinal Steady-State signal calculation given the flip angle in degrees
%

%%
function [SS_Signal] = LongSS_Signal(FlipAngle,TRep,TE,T1,T2)

SS_Signal = sind(FlipAngle).*exp(-TE/T2).*(1-exp(-TRep/T1))./(1-(exp(-TRep/T1)).*cosd(FlipAngle));

end