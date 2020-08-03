function [Transv_SS_Signal] = TransvSS_Signal(FlipAngle,TRep,TE,T1,T2)
% Transverse Steady-State signal calculation given the flip angle in degrees
% by T.T. Fernandes, September 2019 - LarSys - Instituto Superior Tecnico - Universidade de Lisboa
% usage:SS_Signal = LongSS_Signal(FlipAngle,TRep,TE,T1,T2)
%
% :parameters: FlipAngle: Flip angle in degrees
% :parameters: TRep: Time of readout in s
% :parameters: TE: Echo time time in seconds
% :parameters: T1: T1 time in seconds
% :parameters: T2: T2 time in seconds
%
% :returns: Transv_SS_Signal: Transverse Steady-State signal

Transv_SS_Signal = sind(FlipAngle).*exp(-TE/T2).*(1-exp(-TRep/T1))./(1-(exp(-TRep/T1)).*cosd(FlipAngle));

end
