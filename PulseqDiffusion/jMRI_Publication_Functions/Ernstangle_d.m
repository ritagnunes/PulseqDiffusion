function [ernstAngle] = Ernstangle_d(TRep,T1)
% It obtains the Ernst Angle Calculation in degrees
% by T.T. Fernandes, September 2019 - LarSys - Instituto Superior Tecnico - Universidade de Lisboa
% usage:ernstAngle = Ernstangle_d(TRep,T1)
%
% :parameters: TRep: Time of readout in s
% :parameters: T1: T1 time in seconds
%
% :returns: ernstAngle: ernst Angle in degrees

ernstAngle = acosd( exp(-TRep./T1) );

end
