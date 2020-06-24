%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior T�cnico - Universidade de Lisboa
%
% name: TE_time
% type: Function
% description: Echo Time (TE) (s) for a Spiral sequence
%

%%
function [TE_spiral] = TE_time(b_value,gmMax,gamma)

TE_spiral = ((12 .* b_value) ./ ((gamma.*2*pi).^2 * gmMax.^2)).^(1/3);

end