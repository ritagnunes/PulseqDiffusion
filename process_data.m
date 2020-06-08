%Script for reconstructing the images, pre-process them 
%and estimate Diffusion Tensor or simple ADC maps

%file name for the input k-space data
kname='Example_Data/kInVivo12dirs1b20slices';

%name for the output reconstructed image
img_name='DWI_InVivo12dirs1b20slices';

% In-vivo or phantom?
invivo=1;

%load data
load(kname);

%to be modified unless dim_struct is
%provided in example data
if ~exist('dim_struct','var')
    %diffusion acquisition parameters
    Nb=1;%excluding non-DWI
    bmax=1000;
    Ndir=12;
    Nsl=20;   

    %Partial Fourier factor 
    PFourier=0.75;
else
    bmax=dim_struct.bmax;
    Ndir=dim_struct.Ndir;
    Nb=dim_struct.Nb;    
end

%Image spatial resolution
Resx=2.5;
Resy=2.5;
%slice thickness including gaps if they exist
Resz=2.5;%

% model 'dti' or 'adc'
model = 'dti'

%include eddy current and motion correction pre-processing
eddy=0;

%for consistency with the FSL reference frame, 
%we may need to multiply the y gradient component by -1
fsl_flag=1;

%%
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% MATLAB Central File Exchange
addpath('nifti_tools')
setenv('FSLOUTPUTTYPE','NIFTI')

[bvals, bvecs, Nb0s] = generate_bvalsbvecs(bmax,Ndir,Nb,fsl_flag);
save('bvals','bvals','-ASCII')
save('bvecs','bvecs','-ASCII')

if ~exist('dim_struct','var')
    dim_struct=struct('Nb',Nb,'Ndir',Ndir,'Nb0s',Nb0s,'Nsl',Nsl,'PFourier',PFourier);
end
%recon sum-of-squares image including ghost correction
im_sos = recon_imgs(k, dim_struct, 1);

ni=make_nii(im_sos,[Resx Resy Resz]);
save_nii(ni,strcat(img_name,'.nii'));

if invivo
    %perform brain extraction
    strCmd=sprintf('fslroi %s nonDWI 0 1',img_name);
    [st,res]=system(strCmd);
    strCmd='bet nonDWI nonDWI_brain -Z -m -f 0.4';
    [st,res]=system(strCmd);
else
    %define masking threshold
    midSl=round(Nsl/2);
    midImg=im_sos(:,:,midSl,1);
    thr = prctile(midImg(:),80);
    mask=zeros(size(im_sos(:,:,:,1)));
    mask(im_sos(:,:,:,1)>=thr)=1;    
end

if eddy
    %perform motion and eddy-current distortion correction
    disp('Performing eddy current distortion correction...')
    generate_eddyfiles(dim_struct);
    strCmd=sprintf('eddy_openmp --imain=%s --mask=nonDWI_brain_mask --bvals=bvals --bvecs=bvecs --out=%s_cor --acqp=acqp --index=index',img_name,img_name);
    [st,res]=system(strCmd);
end

switch lower(model)
    case 'dti'
        %estimate Diffusion Tensor and derived maps
        disp('estimating diffusion tensor...')
        if eddy
            strCmd=sprintf('dtifit -k %s_cor -m nonDWI_brain_mask -o %s_dti -b bvals -r bvecs',img_name,img_name);
            [st,res]=system(strCmd);
        else
            strCmd=sprintf('dtifit -k %s -m nonDWI_brain_mask -o %s_cor_dti -b bvals -r bvecs',img_name,img_name);
            [st,res]=system(strCmd);
        end
    case 'adc'
         %estimate ADC maps
        disp('estimating ADC map...')
        if invivo            
            tmp=load_nii('nonDWI_brain_mask.nii');
            mask=tmp.img;
        end
        adc=estimate_adc(im_sos,mask, dim_struct, bmax);
        ni=make_nii(adc,[Resx Resy Resz]);
        save_nii(ni,strcat(img_name,'_adc.nii'));
    otherwise
        disp('model should be equal to either dti or adc')
end

disp('Done!')


