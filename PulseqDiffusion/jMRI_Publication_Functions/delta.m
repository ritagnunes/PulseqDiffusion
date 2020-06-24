%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior T�cnico - Universidade de Lisboa
%
% name: delta
% type: Function
% description: delta for TE in (s) for an EPI sequence. It follows the demonstration provided in the 'Calc_delta_EPI.pdf' GITHUB folder 
%

%%
function [DELTA] = delta(gmMax,part_Tend,aux_delta,gamma)

DELTA = (-part_Tend)/2  - ((1-1*i*sqrt(3))*(gamma.*2*pi*gmMax.*part_Tend).^2) ./ (4*aux_delta)  -  ((1-1*i*sqrt(3))*aux_delta) ./ (4*(gamma.*2*pi*gmMax).^2);

end