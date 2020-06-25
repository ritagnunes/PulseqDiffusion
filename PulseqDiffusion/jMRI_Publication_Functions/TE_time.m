function [TE_spiral] = TE_time(b_value,gmMax,gamma)
% Echo Time (TE) (s) for a Spiral sequence
% by T.T. Fernandes, September 2019 - LarSys - Instituto Superior Tecnico - Universidade de Lisboa
% usage: TE_spiral = TE_time(b_value,gmMax,gamma)
%
% :parameters: b_value: b-value in s*mm^-2
% :parameters: gmMax: Gradient maximum amplitude in T/m
% :parameters: gamma: gyromagnetic constant in Hz/T
%
% :returns: TE_spiral: Longitudinal Steady-State signal

TE_spiral = ((12 .* b_value) ./ ((gamma.*2*pi).^2 * gmMax.^2)).^(1/3);

end
