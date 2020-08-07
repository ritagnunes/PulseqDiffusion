%% Test Example In Vivo data 12 directions, Nb 1, 20 slices
t=load('../Example_Data/kInVivo12dirs1b20slices.mat');
tmp=load('../Example_Data/ImgInVivo12dirs1b20slices.mat');
impred=tmp.im_sos;
imact=recon_imgs(t.k, t.dim_struct);  
assertWithAbsTol(imact,impred,'Reconstructed images do not match predictions: InVivo12dirs1b20slices')

%% Test Example In Vivo data 3 directions, Nb 3, 20 slices
t=load('../Example_Data/kInVivo3dirs3b20slices.mat');
tmp=load('../Example_Data/ImgInVivo3dirs3b20slices.mat');
impred=tmp.im_sos;
imact=recon_imgs(t.k, t.dim_struct);  
assertWithAbsTol(imact,impred,'Reconstructed images do not match predictions: InVivo3dirs3b20slices')

%% Test Example In Vivo data 12 directions, Nb 1, 20 slices
t=load('../Example_Data/kPhantom3dirs5b3slices.mat');
tmp=load('../Example_Data/ImgPhantom3dirs5b3slices.mat');
impred=tmp.im_sos;
imact=recon_imgs(t.k, t.dim_struct);
assertWithAbsTol(imact,impred,'Reconstructed images do not match predictions: Phantom3dirs5b3slices')
