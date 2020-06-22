function  generate_eddyfiles(dimpars)
% Generates  acqp and index text files required for running eddy (FSL) to perform 
% eddy current distortion and motion correction. 
% Assuming for the moment that all the data was acquired using a single phase encode
% direction and that topup B0 estimation/distortion correction will not be carried out
%
% usage: generate_eddyfiles(dimpars)
%
% :parameters: dimpars: structure containg data dimensions
%                       dimpars.Nb0s: number of non-DWI volumes
%                       dimpars.Nb:   number of non-zero b-values
%                       dimpars.Ndir: number of sampled gradient directions
%                       dimpars.Nsl:  number of slices

Nvol = dimpars.Nb0s + dimpars.Nb*dimpars.Ndir;

index = ones(Nvol,1);

%Phase encoding along the y direction, as topup will not be applied,
%inputting an accurate echo spacing value is not necessary
acqp = [0 1 0 0.05];

save('index','index','-ASCII')
save('acqp','acqp','-ASCII')

end