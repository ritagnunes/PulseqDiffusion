%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior T�cnico - Universidade de Lisboa
%
% name: DWSignal
% type: Function
% description: Diffusion signal calculation
%

%%
function [DWs] = DWSignal(b_value,ADC)

DWs = exp(-b_value*ADC);

end