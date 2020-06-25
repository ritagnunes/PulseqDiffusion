function [DWs] = DWSignal(b_value,ADC)
% Diffusion signal calculation
% by T.T. Fernandes, September 2019 - LarSys - Instituto Superior Tecnico - Universidade de Lisboa
% usage:DWs = DWSignal(b_value,ADC)
%
% :parameters: b_value: value associated with the test in s*mm^-2
% :parameters: ADC: ADCs for specific tissue in mm^2 * s^-1
%
% :returns: DWs value

DWs = exp(-b_value*ADC);

end
