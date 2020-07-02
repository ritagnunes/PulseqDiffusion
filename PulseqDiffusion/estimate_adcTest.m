%% Test Example In Vivo data 12 directions, Nb 1, 20 slices
t=load('../Example_Data/kInVivo12dirs1b20slices.mat');
tt=load('../Example_Data/ImgInVivo12dirs1b20slices.mat');
ttt=load('../Example_Data/ADCInVivo12dirs1b20slices.mat');
ADCpred=ttt.ADC;
ADCact=estimate_adc(tt.im_sos, ones(size(tt.im_sos(:,:,:,1))), t.dim_struct, t.dim_struct.bmax);
assertValue(ADCact,ADCpred)

%% Test Example In Vivo data 3 directions, Nb 3, 20 slices
t=load('../Example_Data/kInVivo3dirs3b20slices.mat');
tt=load('../Example_Data/ImgInVivo3dirs3b20slices.mat');
ttt=load('../Example_Data/ADCInVivo3dirs3b20slices.mat');
ADCpred=ttt.ADC;
ADCact=estimate_adc(tt.im_sos, ones(size(tt.im_sos(:,:,:,1))), t.dim_struct, t.dim_struct.bmax);
assertValue(ADCact,ADCpred)

%% Test Example Phantom data 3 directions, Nb 5, 3 slices
t=load('../Example_Data/kPhantom3dirs5b3slices.mat');
tt=load('../Example_Data/ImgPhantom3dirs5b3slices.mat');
ttt=load('../Example_Data/ADCPhantom3dirs5b3slices.mat');
ADCpred=ttt.ADC;
ADCact=estimate_adc(tt.im_sos, ones(size(tt.im_sos(:,:,:,1))), t.dim_struct, t.dim_struct.bmax);
assertValue(ADCact,ADCpred)

function assertValue(actVal,expVal)
% Helper function to assert equality
% Takes two images and compares them 

tf = isequal(actVal,expVal);
assert(tf, 'Recon images do not match!');
end