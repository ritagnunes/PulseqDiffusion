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
% assuming input k-space size is N x Nc x Nimg.N
% N - matrix size (assuming square field-of-view)
% Nc - number of coil channels
% where N and Nc are estimated from the data and
% Nvol - number of volumes (1 + Nb0s + Nb*Ǹdir)
% the extra volume correponds to calibration data for EPI ghost-correction
% Nimg - number of imgs + calibration data = Nvol * Nsl
%
%gc: flag for ghost correction - On by default
N = size(k,1);
Nc = size(k,2);

Nsl = dimpars.Nsl;
Nvol = dimpars.Nb0s + dimpars.Nb*dimpars.Ndir;
Nimg = Nsl*Nvol;
PFourier = dimpars.PFourier;

%gc: flag for ghost correction - On by default
if nargin<3
    gc=1;
end

if or(PFourier < 0.6, PFourier >1)
    disp('The partial Fourier parameter has to been within the [0.6, 1.0] interval')
    return
end

if ~isequal(Nimg+Nsl,size(k,3)/(N*PFourier))
    disp('Inconsistent data size! Check provided dimensions.')
    return
end

% extra volume (Nsl slices) correponds to calibration data for EPI ghost-correction
if PFourier<1
    Ny=N*PFourier;
    pad=N-Ny;
    k=reshape(k,[N Nc Ny Nsl Nvol+1]);
else
    k=reshape(k,[N Nc N Nsl Nvol+1]);
end

%calibration data is not required after reconstruction
im=zeros([N N Nc Nsl Nvol]);

%read calibration lines first and reorder (N x N x Nc x Nsl)
kcalib=permute(k(:,:,:,:,1),[1 3 2 4]);

%EPI re-ordering: mirror even k-space lines
kcalib(:,2:2:end,:,:)=permute(k(end:-1:1,:,2:2:end,:,1),[1 3 2 4]);

%basic ghost correction is carried out by subtracting
%the phase of projection data (pcor)
hybcalib=fftshift(fft(ifftshift(kcalib,1),[],1),1);
pcor=exp(-1i*angle(hybcalib));

for v=2:Nvol+1
    for s=1:Nsl
        for c=1:Nc
            ktmp=squeeze(k(:,c,:,s,v));
            %EPI re-ordering: mirror even k-space lines
            ktmp(:,2:2:end)=squeeze(k(end:-1:1,c,2:2:end,s,v));
            if gc
                hybtmp=fftshift(fft(ifftshift(ktmp,1),[],1),1);
                hybcor=hybtmp.*pcor(:,:,c,s);
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
                    im(:,:,c,s,v-1)=fftshift(fft2(ifftshift(kn)));
                else
                    im(:,:,c,s,v-1)=fftshift(fft(ifftshift(hybcor,2),[],2),2);
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
                    im(:,:,c,s,v-1)=fftshift(fft2(ifftshift(kn)));
                else
                    im(:,:,c,s,v-1)=fftshift(fft2(ifftshift(ktmp)));
                end
            end
        end
    end
end
%calculate sum-of-squares image, combining all coils
im_sos=squeeze(sum(im.*conj(im),3)).^0.5;
im_sos = reshape(im_sos,[N N Nsl Nvol]);

%reorder slices
im_sos = im_sos(:,:,end:-1:1,:);

end

