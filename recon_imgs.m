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
% dimpars.PFourier - partial Fourier factor
%
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
PFourier = dimpars.PFourier;

%gc: flag for ghost correction - On by default
if nargin<3
    gc=1;
end

if ~isequal(Nimg,size(k,3)/N)
    disp('Inconsistent data size! Check provided dimensions.')
end

if or(PFourier < 0.6, PFourier >1)
    disp('The partial Fourier parameter has to been within the [0.6, 1.0] interval')
    return
end
    
if PFourier<1
    Ny=N*PFourier;
    pad=N-Ny;
    k=reshape(k,[N Nc Ny Nimg]);
else
    k=reshape(k,[N Nc N Nimg]);
end

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
            hybtmp=fftshift(fft(ifftshift(ktmp,1),[],1),1);
            hybcor=hybtmp.*pcor(:,:,c);
             %perform homodyne reconstruction
            if PFourier<1
                kcor=fftshift(ifft(ifftshift(hybcor,1),[],1),1);
                % zero-fill
                kn=zeros(N);
                kn(:,pad+1:end)=kcor;
                %keep only central portion for phase corrections
                kn(:,end:-1:end-pad+1)=0;
                %perform phase-correction
                img=fftshift(fft2(ifftshift(kn)));
                img=img.*exp(-1i*angle(img));
                % back to k-space
                kn=fftshift(ifft2(ifftshift(img)));
                kn(:,1:pad)=0;
                
                % apply Hermitian symmetry to fill in unacquired data
                kn(:,1:pad)=conj(kn(end:-1:1,end:-1:end-pad+1));
                im(:,:,c,ii-1)=fftshift(fft2(ifftshift(kn)));
            else
                im(:,:,c,ii-1)=fftshift(fft(ifftshift(hybcor,2),[],2),2);
            end
        else
            %perform homodyne reconstruction
            if PFourier<1
                % zero-fill
                kn=zeros(N);
                kn(:,pad+1:end)=ktmp;
                %keep only central portion for phase corrections
                kn(:,end:-1:end-pad+1)=0;
                %perform phase-correction
                img=fftshift(fft2(ifftshift(kn)));
                img=img.*exp(-1i*angle(img));
                % back to k-space
                kn=fftshift(ifft2(ifftshift(img)));
                kn(:,1:pad)=0;
                
                % apply Hermitian symmetry to fill in unacquired data
                kn(:,1:pad)=conj(kn(end:-1:1,end:-1:end-pad+1));
                im(:,:,c,ii-1)=fftshift(fft2(ifftshift(kn)));
            else
                im(:,:,c,ii-1)=fftshift(fft2(ifftshift(ktmp)));
            end
        end
    end
end
    
end
%calculate sum-of-squares image, combining all coils
im_sos=squeeze(sum(im.*conj(im),3)).^0.5;

