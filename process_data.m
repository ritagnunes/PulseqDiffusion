%Script for reconstructing the images, pre-process them 
%and estimate Diffusion Tensor

kname='kInVivo20sl1b60dirs';
img_name='DWI_InVivo1b60dirs';

%diffusion acquisition parameters
Nb0s=1;
Nb=1;
bmax=1000;
Ndir=60;
Nsl=1;

%Partial Fourier factor and image spatial resolution
PFourier=0.75;
Resx=2.5;
Resy=2.5;
%slice thickness including gaps if they exist
Resz=2.5;%

%include eddy current and motion correction pre-processing
eddy=0;

%%
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% MATLAB Central File Exchange
addpath('nifti_tools')

load(kname);

dim_struct=struct('Nb0s',Nb0s,'Nb',Nb,'Ndir',Ndir,'Nsl',Nsl,'PFourier',PFourier);
%recon sum-of-squares image including ghost correction
im_sos = recon_imgs(k, dim_struct, 1);

ni=make_nii(im_sos,[Resx Resy Resz]);
save_nii(ni,strcat(img_name,'.nii'));

[bvals, bvecs] = generate_bvalsbvecs(bmax,Ndir,Nb0s,Nb,1);
save('bvals','bvals','-ASCII')
save('bvecs','bvals','-ASCII')

if eddy
    %perform motion and eddy-current distortion correction
    disp('Performing eddy current distortion correction...')
    strCmd=sprintf('dtifit -k %s -m mask.nii -o dti -b bvals -r bvecs',img_name);
    [st,res]=system(strCmd);
end

%perform brain extraction
if eddy
    strCmd=sprintf('bet -k %s_cor %s_brain -Z -m',img_name,img_name);
    [st,res]=system(strCmd);
else    
    strCmd=sprintf('bet -k %s %s_brain -Z -m',img_name,img_name);
    [st,res]=system(strCmd);
end

%estimate Diffusion Tensor and derived maps
if eddy
    strCmd=sprintf('dtifit -k %s_cor -m %s_brain -o %s_dti -b bvals -r bvecs',img_name,img_name,img_name);
    [st,res]=system(strCmd);
else    
    strCmd=sprintf('dtifit -k %s -m %s_brain -o %s_dti -b bvals -r bvecs',img_name,img_name,img_name);
    [st,res]=system(strCmd);
end





