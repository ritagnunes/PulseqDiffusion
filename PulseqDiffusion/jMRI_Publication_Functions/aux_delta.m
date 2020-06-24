%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior T�cnico - Universidade de Lisboa
%
% name: aux_delta
% type: Function
% description: Auxiliar for calculation of delta for an EPI sequence
%

%%
function [aDelta] = aux_delta(b_value,gmMax,part_Tend,gamma)

aDelta = (  6*b_value*(gamma.*2*pi*gmMax).^4  - (gamma.*2*pi*gmMax).^6.*part_Tend.^3  +  2*sqrt(3)*sqrt(  3*b_value^2*(gamma.*2*pi*gmMax).^8 - b_value*(gamma.*2*pi*gmMax).^10.*part_Tend.^3  )  ).^(1/3);

end