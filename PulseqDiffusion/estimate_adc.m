function ADC = estimate_adc(img, mask, dimpars, bmax)
% Estimates one ADC map per input direction
% usage: ADC=estimate_adc(img, mask, dimpars, bmax)
%
% :parameters: img:  input image (minimum 1 non-zero b-value, 1 direction)
% :parameters: mask: defines where to estimate the ADC maps 
%                    mask should have the same dimensions as img
% :parameters: dimpars: structure containg data dimensions
%                       dimpars.Nb0s: number of non-DWI volumes
%                       dimpars.Nb:   number of non-zero b-values
%                       dimpars.Ndir: number of sampled gradient directions
%                       dimpars.Nsl:  number of slices
% :parameters: bmax: maximum non-zero b-value in s/mm2 (optional - default 1000 s/mm2)
%                    when dimpars.Nb > 1, a linear spacing is assumed for the sampled b-values
%
% :returns: ADC map in mm2/s

if nargin<4 
    bmax=1000;
end

Nb=dimpars.Nb;
Nb0s=dimpars.Nb0s;
Ndir=dimpars.Ndir;
Nsl=dimpars.Nsl;

N=size(img,1);
N2=N^2;

b=linspace(0,bmax,Nb+1);
b=cat(2,b,zeros(1,Nb0s-1));

ADC=zeros(N2,Nsl,Ndir);
img=reshape(img,[N2,Nsl,Nb0s+Nb*Ndir]);

for s=1:Nsl
    for d=1:Ndir
        imgd=squeeze(img(:,s,[1:Nb0s,(Nb0s+d):Ndir:size(img,3)]));
        parfor ii=1:N2
            if mask(ii)>0           
                sig=imgd(ii,:);
                if size(sig,2)>2
                    [~,m,~]=regression(b,log(sig));
                    ADC(ii,s,d)=-m;
                else
                    ADC(ii,s,d)= (log(sig(1))-log(sig(2)))/bmax;
                end
            end
        end
    end
end
ADC=reshape(ADC,[N N Nsl Ndir]);

end

