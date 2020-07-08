function [DELTA] = delta(gmMax,part_Tend,aux_delta,gamma)
% Delta for TE in (s) for an EPI sequence. It follows the demonstration provided in the 'Calc_delta_EPI.pdf' GITHUB folder .
% by T.T. Fernandes, September 2019 - LarSys - Instituto Superior Tecnico - Universidade de Lisboa
% usage: DELTA = delta(gmMax,part_Tend,aux_delta,gamma)
%
% :parameters: gmMax: Gradient maximum amplitude in T/m
% :parameters: part_Tend: Sequence readout time until k_space center in s
% :parameters: aux_delta: see 'aux_delta.m' function
% :parameters: gamma: gyromagnetic constant in Hz/T
%
% :returns: DELTA: Time of diffusion gradient in s

DELTA = (-part_Tend)/2  - ((1-1*i*sqrt(3))*(gamma.*2*pi*gmMax.*part_Tend).^2) ./ (4*aux_delta)  -  ((1-1*i*sqrt(3))*aux_delta) ./ (4*(gamma.*2*pi*gmMax).^2);

end
