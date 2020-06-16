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
    Nb=dim_struct.Nb;  
    bmax=dim_struct.bmax;
    Ndir=dim_struct.Ndir;
    Nsl=dim_struct.Nsl;    
end

%Image spatial resolution
Resx=2.5;
Resy=2.5;
%slice thickness including gaps if they exist
Resz=2.5;%

% model 'dti' or 'adc'
model = 'dti'

%include eddy current and motion correction pre-processing 
%only considered for invivo acquisitions
eddy=1;

%for consistency with the FSL reference frame, 
%we may need to multiply the y gradient component by -1
fsl_flag=1;

%%
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% MATLAB Central File Exchange
addpath('nifti_tools')
if isempty(which('make_nii'))
    disp('Nifti tools could not be found ')
    return
end

if isempty(getenv('FSLDIR'))
    disp('FSL installation not found')
    if or(or(invivo,eddy),strcmp(lower(model),'dti'))
        disp('FSL is required to run with these options')
        return
    end
end

%nifti_tools is not compatible with NIFTI_GZ
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
    if st
        disp('Error found when trying to extract non-DWI volume...')
        disp(res)
        return
    end
    strCmd='bet nonDWI nonDWI_brain -Z -m -f 0.4';
    [st,res]=system(strCmd);
    if st
        disp('Error found when trying to perform brain extraction...')
        disp(res)
        return
    end
else
    %define masking threshold
    midSl=round(Nsl/2);
    midImg=im_sos(:,:,midSl,1);
    thr = prctile(midImg(:),80);
    mask=zeros(size(im_sos(:,:,:,1)));
    mask(im_sos(:,:,:,1)>=thr)=1;
    nimsk=make_nii(mask,[Resx Resy Resz]);
    save_nii(nimsk,'nonDWI_phantom_mask.nii');
end

if and(eddy,not(invivo))
    disp('Warning: eddy current correction only performed on in vivo data')
    disp('Switching off eddy correction')
    eddy=0;
end

if eddy
    %perform motion and eddy-current distortion correction
    disp('Performing eddy current distortion correction...')
    generate_eddyfiles(dim_struct);
    strCmd=sprintf('eddy_openmp --imain=%s --mask=nonDWI_brain_mask --bvals=bvals --bvecs=bvecs --out=%s_cor --acqp=acqp --index=index',img_name,img_name);
    [st,res]=system(strCmd);
    if st
        disp('Error found when trying to run eddy...')
        disp(res)
        return
    end
end

if and(Ndir<6,strcmp(lower(model),'dti'))
    fprintf('Warning: Ndir=%d is lower than the required minimum of 6 for estimating a tensor\n',Ndir)
    disp('Switching to adc model')
    model='adc';
end

switch lower(model)
    case 'dti'
        %estimate Diffusion Tensor and derived maps
        disp('Estimating diffusion tensor...')
        if invivo
            msk_name='nonDWI_brain_mask';
        else
            msk_name='nonDWI_phantom_mask'; 
        end
        if eddy
            strCmd=sprintf('dtifit -k %s_cor -m %s -o %s_cor_dti -b bvals -r bvecs',img_name,msk_name,img_name);
            [st,res]=system(strCmd);
        else
            strCmd=sprintf('dtifit -k %s -m %s -o %s_dti -b bvals -r bvecs',img_name,msk_name,img_name);
            [st,res]=system(strCmd);
        end
        if st
            disp('Error found when trying to run dtifit...')
            disp(res)
            return
        end
    case 'adc'
         %estimate ADC maps
        disp('Estimating ADC map...')
        if invivo            
            tmp=load_nii('nonDWI_brain_mask.nii');
            mask=tmp.img;
        else
            tmp=load_nii('nonDWI_phantom_mask.nii');
            mask=tmp.img;
        end
        if eddy
            tmp=load_nii(sprintf('%s_cor.nii',img_name));
            im_cor=tmp.img;
            adc=estimate_adc(im_cor,mask, dim_struct, bmax);
        else
            adc=estimate_adc(im_sos,mask, dim_struct, bmax);
        end
        ni=make_nii(adc,[Resx Resy Resz]);
        save_nii(ni,strcat(img_name,'_adc.nii'));
    otherwise
        disp('Model should be equal to either dti or adc')
end

disp('Done!')


