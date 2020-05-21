function im_sos = recon_imgs(k, dimpars ,gc)
%im_sos = recon_imgs(k, dimpars )
%k: k-space data should already have been be loaded into Matlab 
%using a raw data reading function
%compatible with the used MRI scanner's raw data structure
%
%dimpars: structure containg data dimensions
% dimpars.Nb0s - number of non-DWI volumes
% dimpars.Nb - number of non-zero b-values
% dimpars.Ndir - number of sampled gradient directions
% dimpars.Nsl - number of slices
%multi-slice 2D data containing multiple volumes
% assuming input k-space size is N x Nc x Nimgs
% N - matrix size (assuming square field-of-view)
% Nc - number of coil channels
% where N and Nc are estimated from the data and
% Nvol - number of volumes (Nb0 + Nb*Ǹdir)
% Nimg - number of imgs = (1+Nvol) * Nsl
% extra volume correponds to calibration data for EPI ghost-correction
%
%gc: flag for ghost correction - On by default
N = size(k,1);
Nc = size(k,2);

Nvol = dimpars.Nb0 + dimpars.Nb*dimpars.Ndir;
Nimg = dimpars.Nsl*(Nvol+1);

%gc: flag for ghost correction - On by default
if nargin<3
    gc=1;
end

if ~isequal(Nimg,size(k,3)/N)
    disp('Inconsistent data size! Check provided dimensions.')
end


k=reshape(k,[N Nc N Nimg]);

%calibration data is not required after reconstruction
im=zeros([N Nc N Nimg]);

% reorder to obtain
% N x N x Nc x Nimg
im=permute(im,[1 3 2 4]);

%read calibration lines first and reorder (N x N x Nc)
kcalib=permute(k(:,:,:,1),[1 3 2]);

%EPI re-ordering: mirror even k-space lines
kcalib(:,2:2:end,:)=permute(k(end:-1:1,:,2:2:end,1),[1 3 2]);

%basic ghost correction is carried out by subtracting 
%the phase of projection data (pcor)
hybcalib=fftshift(fft(kcalib,[],1),1);
pcor=exp(-1i*angle(hybcalib));

for ii=2:Nimg
    for c=1:Nc
        ktmp=squeeze(k(:,c,:,ii));
        %EPI re-ordering: mirror even k-space lines
        ktmp(:,2:2:end)=squeeze(k(end:-1:1,c,2:2:end,ii));
        if gc
            hybtmp=fftshift(fft(ktmp,[],1),1);
            hybcor=hybtmp.*pcor(:,:,c);
            im(:,:,c,ii-1)=fftshift(fft(hybcor,[],2),2);
        else
            im(:,:,c,ii-1)=fftshift(fft2(ktmp));
        end
    end
end

%calculate sum-of-squares image, combining all coils
im_sos=squeeze(sum(im.*conj(im),3)).^0.5;

