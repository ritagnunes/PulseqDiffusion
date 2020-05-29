function  generate_eddyfiles(dimpars)
%  generate_eddyfiles(dimpars)

% produce acqp and index text files for running eddy to perform eddy current
% distortion and motion correction 
% assuming for the moment that all the data was acquired using a single phase encode
% direction and that topup B0 estimation/distortion correction will not be carried out

Nvol = dimpars.Nb0s + dimpars.Nb*dimpars.Ndir;

index = ones(Nvol,1);

%Phase encoding along the y direction, as topup will not be applied,
%inputting an accurate echo spacing value is not necessary
acqp = [0 1 0 0.05];

save('index','index','-ASCII')
save('acqp','acqp','-ASCII')

end