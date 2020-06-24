%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior T�cnico - Universidade de Lisboa
%
% name: Ernstangle_d
% type: Function
% description: It obtains the Ernst Angle Calculation (degrees)
%

%%
function [ernstAngle] = Ernstangle_d(TRep,T1)

ernstAngle = acosd( exp(-TRep./T1) );

end