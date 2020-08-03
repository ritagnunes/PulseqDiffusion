%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior Tecnico - Universidade de Lisboa
%
% name: SNR_CNR_Study_Simulation
% type: Script
% description:
% A tool for predicting the signal-to-noise ratio (SNR) observed for
% different brain tissues (gray matter, white matter, cerebrospinal 
% fluid - CSF) and the contrast-to-noise ratio (CNR) of an acute stroke
% lesion relative to those tissues. The implemented code predicts the 
% SNR and CNR per time unit for a spin-echo diffusion-weighted sequence, 
% taking into account the steady-state value for the transverse
% magnetization. For that purpose, we adapted an expression used to
% predict the SNR per tissue per time unit for a spoiled gradient echo
% sequence (Marques, JP et al, 2019). This simulation considers the possibility to use either
% EPI or spiral readouts and considers the impact of: B0, max gradient
% amplitude, spatial resolution and b-value (comparing the achieved SNR
% with that of a typical 1.5T clinical scanner).
% Inspired in the function SNRandCNRfieldDependenceJMRI.m 
% from Marques, JP 2019 - DOI: 10.1002/jmri.26637
%
% It loads two '.mat' files:
%   'Tend_vector_1shot.mat' - vector of values for time of readout for a
%      one shot spiral with different max. gradients (between [1 2:2:30]
%      in (T/m) ) and with different spatial resolution (between [2:0.5:6] 
%      in (mm) ).
%   'Tend_vector_EPI.mat'   - vector of values for time of readout for a
%      EPI sequence with different max. gradients (between [1 2:2:30]
%      in (T/m) ) and with different spatial resolution (between [2:0.5:6] 
%      in (mm) ).
%

%% Add path

addpath(genpath('./SNR_CNR_Study'));

clear all
% close all
clc
tic

%% 0 - Select values to visualize in the figures

% Values for plot
valuegmMax  = 0.016;   % value for plot - Max Gradient (T/m)
valueB      = 400;     % value for plot - b-value (s*mm^-2)
valueRes    = 2.5;     % value for plot - Resolution (mm)
valueB0     = 0.05;     % value for plot - B0 field (T)

% defaults
seq         = 'S';     % Choose sequence: 'S' - Spiral or 'E' - partial Echo Planar Imaging
percSNRinit = 1/40;    % value for colour bar relative to SNRinit = 40
sizeB0res   = 9;       % for resolution plots
logscale    = 0;       % Logorithmic scale = 1 or = 0 no log scale+
refScan     = 1;       % '1' - ref SNR & CNR values to EPI clinical Scanner (SNR_WM = 40 per time unit) or '0' - arbitary units


%% 1 - initializing variables

% ... initial variables ... 
gamma        = 42.577e6;     % gyromagnetic constant (Hz/T)
SNR_PowerLaw = 1.5;          % Power law for the SNR, see JP Marques, 2019

% ... range of b-value ... 
b               = 10;               % b-value (s*mm^-2)
b_vector        = [10 50:50:1000];  % b-value vector (s*mm^-2)
convertB        = b*1e6;            % b-value (s*m^-2)
convertB_vector = b_vector*1e6;     % b-value vector (s*m^-2)

% ... range of resolution ... 
resINI = 2;             % original resolution (mm)
resNEW = [2:0.5:6];     % resolution vector (mm)

% ... range of GMmax ... 
gmMaxVector  = [1 2:2:30]*1e-3;     % Maximum Amplitude (T/m)
gmMax        = 0.030;               % Maximum Amplitude (T/m)

% ... range of B0 fields for simulations ... 
B0       = [0.01 0.05 0.1 0.2 0.35 0.5 0.7 1 1.5];  % B0 (T)
refB0init  = 1.5;                                   % refB0 (T)

% ... Load Tend Vector ... 
if seq == 'S'       % load Tend_vector_Spiral
    load('./Example_Data/Tend_vector_1shot.mat') 
    plot_seq = 'Spiral';
elseif seq == 'E'   % load Tend_vector_EPI
    load('./Example_Data/Tend_vector_EPI.mat')
    plot_seq = 'EPI';
end

% ... ADCs for each tissue ... 
D_CSF = 3.2*1e-3;   % Diffusivity of CSF (mm^2 * s^-1) - Breger et. al, 1986 - DOI: 10.1002/mrm.1910030502 
D_GM  = 0.8*1e-3;   % Diffusivity of GM  (mm^2 * s^-1) - Breger et. al, 1986 - DOI: 10.1002/mrm.1910030502 
D_WM  = 0.7*1e-3;   % Diffusivity of WM  (mm^2 * s^-1) - Breger et. al, 1986 - DOI: 10.1002/mrm.1910030502 
D_LS  = 0.55*1e-3;  % Diffusivity of LS  (mm^2 * s^-1) - Knight et. al, 2019 - DOI: 10.3233/bsi-190185

% ... Counting the number of repetitions -N slices ... 
Z_size    = 120;    % dimention of the brain in zz axes (mm)

% ... reference values for SNR study ...
SNR_refinit = 40;                       % for gmmax = 0.030 T/m & Nslices = 11 & WM;
refgmMax    = find(gmMaxVector==gmMax); % gmMax      = 30 mT/m
refbValue   = find(b_vector==b);        % b_value    = 10 s*mm^-2
refRes      = find(resNEW==resINI);     % resolution = 2 mm
refB0idx    = find(B0==refB0init);      % ref B0     = 1.5T

% ... Plots - on (1) / off (0) ...
PlotIntermediate = 0;
PlolPowerLaw     = 0;
PlotTE_TR_Ttotal = 0;
PlotB0variation  = 0;
tests    = 0;          % Tests for each variable that contribute to SNR


%% 2 - initializing relaxation times field dependence
% Longitudinal and transverse relaxation rates using the models presented
% assuming LS occurs WM - Only ADC change
% in:
% - Rooney et al MRM 2007  - DOI: 10.1002/mrm.21122 ;
% - Pohmann et al MRM 2016 - DOI: 10.1002/mrm.25677 ;
% - Bottemley et al, Med Phys, 1984 - DOI:  10.1118/1.595535 ;

% R1 as a function of magnetic field according to Rooney et al, MRM, 2007
% and using a model suggested by Bottemley et al, Med Phys, 1984 - Time (s)
T1_GM = 0.00116 * (gamma*B0).^0.376;
T1_WM = 0.00071 * (gamma*B0).^0.382;
T1_LS = T1_WM;
T1_BL = 0.00335 * (gamma*B0).^0.340;
T1_CSF =1/0.231 * ones(size(B0));

% using Pohmann et al, MRM, 2016 - Time (s)
T2s_GM = 0.090 * exp(-0.142 *B0);
T2s_WM = 0.064 * exp(-0.132 *B0);
T2s_LS = T2s_WM;
T2s_CSF = 0.1*ones(size(B0));


%% 3 - Select ref SNR expressions for partial EPI and for a value of SNR=40

SNR_B0ref     = B0(refB0idx).^SNR_PowerLaw;

% TE reference with partial EPI - results from the derivation of the b-value formula
percTrout     = 0.285714;                                               % percentage of time to get to the centre k_space due to partial EPI
T_endRef      = load('./Example_Data/Tend_vector_EPI.mat');                            % load reference Time of read out for partial EPI
part_Tend_ref = percTrout*T_endRef.Tend_Vector(refgmMax,refRes)*1e-3;   % sequence readout time until k_space center in (s)
auxD_ref      = aux_delta(convertB,gmMax,part_Tend_ref,gamma);          % auxiliar
deltaEPI_ref  = delta(gmMax,part_Tend_ref,auxD_ref,gamma);              % delta in (s)
TE_B0ref      = 2*(part_Tend_ref+abs(deltaEPI_ref));                    % TE = 2 x (delta + sequence readout until k_space center) in (s)

TTotal_B0ref  = T_endRef.Tend_Vector(refgmMax,refRes)*1e-3*(1-percTrout) + TE_B0ref;  % TR per slice in (s)
TRo_B0ref     = T_endRef.Tend_Vector(refgmMax,refRes)*1e-3*(1-percTrout);             % Time of Read Out in (s)
Nslic_ref     = Z_size/resNEW(refRes);                                                % Number of slices in reference

clear part_Tend_ref auxD_ref deltaEPI_ref T_endRef

BW_B0ref = 1 ./ TRo_B0ref;                                          % Band-width for readout time at B0 ref
TR_B0ref = TTotal_B0ref * Nslic_ref;                                % TR at B0 ref
OptimumAlpha_GM_B0ref = Ernstangle_d(TR_B0ref,T1_GM(refB0idx));
OptimumAlpha_WM_B0ref = Ernstangle_d(TR_B0ref,T1_WM(refB0idx));
OptimumAlpha_LS_B0ref = Ernstangle_d(TR_B0ref,T1_LS(refB0idx));

% References values for each tissue (GM, WM, LS) and its respective SNR
% ... GM tissue ...
TransvSSs_GM_B0ref = TransvSS_Signal(OptimumAlpha_GM_B0ref,TR_B0ref,TE_B0ref,T1_GM(refB0idx),T2s_GM(refB0idx));
DWs_GM_B0ref       = DWSignal(b,D_GM);
S_GM_B0ref         = TransvSSs_GM_B0ref * DWs_GM_B0ref;
SNR_GM_ref         = SNR_B0ref * 1/sqrt(TR_B0ref) * 1/sqrt(BW_B0ref) * S_GM_B0ref * (resNEW(refRes)/resINI)^3;

% ... WM tissue ...
TransvSSs_WM_B0ref = TransvSS_Signal(OptimumAlpha_WM_B0ref,TR_B0ref,TE_B0ref,T1_WM(refB0idx),T2s_WM(refB0idx));
DWs_WM_B0ref       = DWSignal(b,D_WM);
S_WM_B0ref         = TransvSSs_WM_B0ref * DWs_WM_B0ref;
SNR_WM_ref         = SNR_B0ref * 1/sqrt(TR_B0ref) * 1/sqrt(BW_B0ref) * S_WM_B0ref * (resNEW(refRes)/resINI)^3;

% ... LS tissue ...
TransvSSs_LS_B0ref = TransvSS_Signal(OptimumAlpha_LS_B0ref,TR_B0ref,TE_B0ref,T1_LS(refB0idx),T2s_LS(refB0idx));
DWs_LS_B0ref       = DWSignal(b,D_LS);
S_LS_B0ref         = TransvSSs_LS_B0ref * DWs_LS_B0ref;
SNR_LS_ref         = SNR_B0ref * 1/sqrt(TR_B0ref) * 1/sqrt(BW_B0ref) * S_LS_B0ref * (resNEW(refRes)/resINI)^3;


%% 4 - variation of SNR for B0, gm, Nslices, res, b_value ...

clear TE TTotal TR BW SNR_B0_SNRinit SNR_SE_GM SNR_SE_WM SNR_SE_LS TRo CNR_SE_GM_LS ...
      SNR_SE_WM CNR_SE_WM_LS S_GM_SNRinit S_WM_SNRinit S_LS_SNRinit 
for bb = 1:size(B0,2)
    
    % ... part 1 - calculation SNR due to the field ...
    SNR_B0_SNRinit{bb} = B0(bb).^SNR_PowerLaw ;
    
    for bval=1:size(b_vector,2)
        for r=1:size(resNEW,2)
            % ... Auxiliar Calculus ...
            if (seq == 'S')       % calculate TE accordingly to Spiral
                TE{bval}       = TE_time(convertB_vector(bval),gmMaxVector,gamma);   % TE per gm in (s)
                TTotal{bval,r} = Tend_Vector(:,r)'*1e-3 + TE{bval};                  % TR per slice for each gm in (s)
                TRo{bval,r}    = Tend_Vector(:,r)'*1e-3;                             % Time of Readout in (s)
            elseif (seq == 'E')   % calculate TE accordingly to partial EPI
                part_Tend      = percTrout*Tend_Vector(:,r)'*1e-3;                   % sequence readout time until k_space center in (s)
                auxD           = aux_delta(convertB_vector(bval),gmMaxVector,part_Tend,gamma); % auxiliar
                deltaEPI       = delta(gmMaxVector,part_Tend,auxD,gamma);            % delta in (s)
                TE{bval}       = 2*(part_Tend+abs(deltaEPI));                        % TE = delta + sequence readout until k_space center in (s)
                TTotal{bval,r} = Tend_Vector(:,r)'*1e-3*(1-percTrout) + TE{bval};    % TR per slice for each gm in (s)
                TRo{bval,r}    = Tend_Vector(:,r)'*1e-3*(1-percTrout);               % Time of Readout in (s)
            end
            Nslices(r)     = Z_size/resNEW(r);
            
            % ... part 2 - calculation variables ...            
            TR{bval,r}           = TTotal{bval,r}' .* Nslices(r);           % TR for each slice in (s)
            TE{bval}             = ones(size(TR{bval,r})).*TE{bval}';       % TE per gm by size of TR in (s)
            TTotal{bval,r}       = ones(size(TR{bval,r})).*TTotal{bval,r}'; % TTotal per gm by size of TR in (s)
            
            OptimumAlphaGM_SNRinit{bb,bval,r}  = Ernstangle_d(TR{bval,r},T1_GM(bb));
            OptimumAlphaWM_SNRinit{bb,bval,r}  = Ernstangle_d(TR{bval,r},T1_WM(bb));
            OptimumAlphaLS_SNRinit{bb,bval,r}  = Ernstangle_d(TR{bval,r},T1_LS(bb));

            SNR_B0_SNRinit{bb}   = ones(size(TR{bval,r})).*SNR_B0_SNRinit{bb};
            BW{bval,r}           = 1 ./ TRo{bval,r};
            
            % ... part 3 - Contribution of noise due to Signal Sequence (GM, WM & Lesion) ...
            TransvSS_GM_SNRinit{bb,bval,r} = TransvSS_Signal(OptimumAlphaGM_SNRinit{bb,bval,r}, ...
                                           TR{bval,r},TE{bval},T1_GM(bb),T2s_GM(bb));
            TransvSS_WM_SNRinit{bb,bval,r} = TransvSS_Signal(OptimumAlphaWM_SNRinit{bb,bval,r}, ...
                                           TR{bval,r},TE{bval},T1_WM(bb),T2s_WM(bb));
            TransvSS_LS_SNRinit{bb,bval,r} = TransvSS_Signal(OptimumAlphaLS_SNRinit{bb,bval,r}, ...
                                           TR{bval,r},TE{bval},T1_LS(bb),T2s_LS(bb));
            
            % ... part 4 - Contribution of noise due to Diffusion ...
            DWs_GM_SNRinit = DWSignal(b_vector(bval),D_GM);
            DWs_WM_SNRinit = DWSignal(b_vector(bval),D_WM);
            DWs_LS_SNRinit = DWSignal(b_vector(bval),D_LS);
            
            % ... part 5 - Combining of signal and diffusion...
            S_GM_SNRinit{bb,bval,r} = TransvSS_GM_SNRinit{bb,bval,r} .* DWs_GM_SNRinit;
            S_WM_SNRinit{bb,bval,r} = TransvSS_WM_SNRinit{bb,bval,r} .* DWs_WM_SNRinit;
            S_LS_SNRinit{bb,bval,r} = TransvSS_LS_SNRinit{bb,bval,r} .* DWs_LS_SNRinit;
            
            % ... part 5.5 - SNR for time unit & referenced to EPI clinical scanner ...
            if refScan == 1
                refSNRinit = SNR_refinit./sqrt(TR_B0ref);
            elseif refScan == 0
                refSNRinit = 1;
            end
            
            % ... part 6 - SNR projection due to field, Sequence, Hardware and with the refSNR ...
            SNR_SE_GM(bb,bval,:,r) = ( (SNR_B0_SNRinit{bb} .* 1./sqrt(TR{bval,r}) ...
                                       .*1./sqrt(BW{bval,r}') .*S_GM_SNRinit{bb,bval,r} ...
                                       .*(resNEW(r)/resINI)^3) ./ SNR_WM_ref ...
                                     ) .* refSNRinit;
            SNR_SE_WM(bb,bval,:,r) = ( (SNR_B0_SNRinit{bb} .* 1./sqrt(TR{bval,r}) ...
                                       .*1./sqrt(BW{bval,r}') .*S_WM_SNRinit{bb,bval,r} ...
                                       .*(resNEW(r)/resINI)^3) ./ SNR_WM_ref ...
                                     ) .* refSNRinit;
            SNR_SE_LS(bb,bval,:,r) = ( (SNR_B0_SNRinit{bb} .* 1./sqrt(TR{bval,r}) ...
                                       .*1./sqrt(BW{bval,r}') .*S_LS_SNRinit{bb,bval,r} ...
                                       .*(resNEW(r)/resINI)^3) ./ SNR_WM_ref ...
                                     ).* refSNRinit;
            
            % ... part  7 - CNR between GM-LS & WM-LS ...
            CNR_SE_GM_LS(bb,bval,:,r) = abs(   SNR_B0_SNRinit{bb} .* 1./sqrt(TR{bval,r}) ...
                                          .*1./sqrt(BW{bval,r}') .*(resNEW(r)/resINI)^3 ...
                                          .* ...
                                            ( S_GM_SNRinit{bb,bval,r}./SNR_WM_ref - ...
                                              S_LS_SNRinit{bb,bval,r}./SNR_WM_ref )  );%.*refSNRinit ;   
            CNR_SE_WM_LS(bb,bval,:,r) = abs(   SNR_B0_SNRinit{bb} .* 1./sqrt(TR{bval,r}) ...
                                           .*1./sqrt(BW{bval,r}') .*(resNEW(r)/resINI)^3 ...
                                           .* ...
                                            ( S_WM_SNRinit{bb,bval,r}./SNR_WM_ref - ...
                                              S_LS_SNRinit{bb,bval,r}./SNR_WM_ref ));%.*refSNRinit;
        end
    end
end


%% pre_5 - cut the top right figure of the gm plots vs SNR & vs CNR
aux_Cut = [1 9; 2 11; 3 13; 4 15];
gm_idxMax = size(gmMaxVector,2);

for uu=1:size(aux_Cut,1)
    % ... cut for SNR ...
    SNR_SE_GM(uu,:,aux_Cut(uu,2):gm_idxMax,:) = NaN;
    SNR_SE_WM(uu,:,aux_Cut(uu,2):gm_idxMax,:) = NaN;
    SNR_SE_LS(uu,:,aux_Cut(uu,2):gm_idxMax,:) = NaN;
    % ... cut for CNR ...
    CNR_SE_GM_LS(uu,:,aux_Cut(uu,2):gm_idxMax,:) = NaN;
    CNR_SE_WM_LS(uu,:,aux_Cut(uu,2):gm_idxMax,:) = NaN;
end


%% 5 - SNR plots for a reference SNR for each Condition of B0
% Using the relationship SNR_new / SNR_ref = (B0new/B0ref)^b

plotgmMax  = find(round(gmMaxVector,3)==round(valuegmMax,3)); % gmMax   = 31 mT/m
plotRes    = find(resNEW==valueRes);        % res     = 4 mm between Slice
plotb      = find(b_vector==valueB);        % b-value = 10 s*mm^-2
plotB0     = find(B0==valueB0);             % B0 = 500mT

valueNslice      = Z_size/valueRes;                % Number of Slices
valueMaxCM       = SNR_refinit*percSNRinit;
valueMaxCMlog    = log(max(max(SNR_SE_WM(:,plotb,:,plotRes))));
valueMinCMlog    = log(min(min(SNR_SE_WM(:,plotb,:,plotRes))));
valueMaxContrast = 2;
valuelogMinContrast = log(min(min(CNR_SE_WM_LS(:,plotb,:,plotRes))));
valuelogMaxContrast = log(max(max(CNR_SE_WM_LS(:,plotb,:,plotRes))));

if logscale == 1    
    % ... 1 - Plot Matrices of SNR for GM Max...
    fSNR_gmMax = figure('name',['logScale - SNR (SE - ',plot_seq,') variation w/ Gm Max (Diffusion) - b_value=', ... 
                                num2str(valueB),' - Nslices=',num2str(valueNslice), ...
                                ' - res=',num2str(valueRes)], ...
                        'numbertitle','off','Color',[1 1 1]);
    set(fSNR_gmMax,'Position',[ 461   474   917   369],'Color',[1 1 1])
    subplot(131)
    imagesc(log(reshape(SNR_SE_GM(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2))))
    hold on
    grid on
    title('logscale - SNR for GM','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
%     xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([valueMinCMlog valueMaxCMlog])
    colorbar
    
    subplot(132)
    imagesc(log(reshape(SNR_SE_WM(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2))))
    hold on
    grid on
    title('logscale - SNR for WM','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([valueMinCMlog valueMaxCMlog])
    colorbar
    
    subplot(133)
    imagesc(log(reshape(SNR_SE_LS(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2))))
    hold on
    grid on
    title('logscale - SNR for LS','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
%     xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([valueMinCMlog valueMaxCMlog])
    colorbar
    
    % ... 2 - Plot Matrices of SNR for resolution ...
    fSNR_Slice = figure('name',['logScale - SNR (SE - ',plot_seq,') variation w/ resolution (Diffusion) - gm=', ... 
                        num2str(valuegmMax),' - Nslices=',num2str(valueNslice),   ...
                        ' - b-value=',num2str(valueB)],'numbertitle','off','Color',[1 1 1]);
    set(fSNR_Slice,'Position',[ 461   474   917   369],'Color',[1 1 1])
    
    subplot(131)
    imagesc(log(reshape(SNR_SE_GM(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2))))
    hold on
    grid on
    title('logscale - SNR for GM','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
%     xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([valueMinCMlog valueMaxCMlog])
    colorbar
    
    subplot(132)
    imagesc(log(reshape(SNR_SE_WM(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2))))
    hold on
    grid on
    title('logscale - SNR for WM','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([valueMinCMlog valueMaxCMlog])
    colorbar
    
    subplot(133)
    imagesc(log(reshape(SNR_SE_LS(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2))))
    hold on
    grid on
    title('logscale - SNR for LS','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
%     xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([valueMinCMlog valueMaxCMlog])
    colorbar
    
else
    % ... 1 - Plot Matrices of SNR for GM Max...
    fSNR_gmMax = figure('name',['SNR (SE - ',plot_seq,') variation w/ Gm Max (Diffusion) - b_value=', ... 
                                num2str(valueB),' - Nslices=',num2str(valueNslice), ...
                                ' - res=',num2str(valueRes)], ...
                        'numbertitle','off','Color',[1 1 1]);
    set(fSNR_gmMax,'Position',[ 461   474   917   369],'Color',[1 1 1])
    subplot(131)
    imagesc(reshape(SNR_SE_GM(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2)))
    hold on
    grid on
    title('SNR for GM','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
%     xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
%     caxis([0 valueMaxCM])
    colorbar
    
    subplot(132)
    imagesc(reshape(SNR_SE_WM(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2)))
    hold on
    grid on
    title('SNR for WM','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
%     caxis([0 valueMaxCM])
    colorbar
    
    subplot(133)
    imagesc(reshape(SNR_SE_LS(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2)))
    hold on
    grid on
    title('SNR for LS','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
%     xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
%     caxis([0 valueMaxCM])
    colorbar
    
    % ... 2 - Plot Matrices of SNR for resolution ...
    fSNR_Slice = figure('name',['SNR (SE - ',plot_seq,') variation w/ resolution (Diffusion) - gm=', ... 
                        num2str(valuegmMax),' - Nslices=',num2str(valueNslice),   ...
                        ' - b-value=',num2str(valueB)],'numbertitle','off','Color',[1 1 1]);
    set(fSNR_Slice,'Position',[ 461   474   917   369],'Color',[1 1 1])
    
    subplot(131)
    imagesc(abs(reshape(SNR_SE_GM(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2))))
    hold on
    grid on
    title('SNR for GM','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
%     xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
%     caxis([0 valueMaxCM])
    colorbar
    
    subplot(132)
    imagesc(reshape(SNR_SE_WM(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2)))
    hold on
    grid on
    title('SNR for WM','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
%     caxis([0 valueMaxCM])
    colorbar
    
    subplot(133)
    imagesc(reshape(SNR_SE_LS(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2)))
    hold on
    grid on
    title('SNR for LS','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
%     xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
%     caxis([0 valueMaxCM])
    colorbar
end

% ... 3 - Plot Values of SNR for b_values ...
fSNR_Slice = figure('name',['SNR (SE) matrix variation w/ b_values (Diffusion) - gm=', ...
    num2str(valuegmMax),' - B0=', num2str(valueB0),' - Nslices=',num2str(valueNslice), ...
    ' - res=',num2str(valueRes)],'numbertitle','off','Color',[1 1 1]);
set(fSNR_Slice,'Position',[ 461   474   917   369],'Color',[1 1 1])

% subplot(131)
plot(b_vector,SNR_SE_GM(plotB0,:,plotgmMax,plotRes),'Color','r')
ylabel([' SNR (au)'],'fontsize',20)
% xlabel(['b_{value}'],'fontsize',20)
title('SNR variation with B-value, for each tissue ','fontsize',20)
set(gca,'FontSize',15)
hold on
% subplot(132)
plot(b_vector,SNR_SE_WM(plotB0,:,plotgmMax,plotRes),'Color','b')
% ylabel([' SNR (au)'],'fontsize',20)
xlabel(['b_{value}'],'fontsize',20)
% title('SNR for WM','fontsize',20)
set(gca,'FontSize',15)
hold on
% subplot(133)
plot(b_vector,SNR_SE_LS(plotB0,:,plotgmMax,plotRes),'Color','g')
% ylabel([' SNR (au)'],'fontsize',20)
% xlabel(['b_{value}'],'fontsize',20)
% title('SNR for LS','fontsize',20)
set(gca,'FontSize',15)
hold on
grid on
lgd = legend('GM','WM','LS');
lgd.FontSize = 20;
lgd.LineWidth = 4;
    
    
%% 6 - auxiliar (Plots for CNR)
if logscale == 1
    % ... Plot Matrices of CNR for gmax & TE ...
    fCNR_gmMax_Contrast = figure('name',['logscale - CNR (SE - ',plot_seq,') variation w/ Gm Max (Diffusion) - b_value=', ...
        num2str(valueB),' - Nslices=',num2str(valueNslice), ...
        ' - res=',num2str(valueRes)], ...
        'numbertitle','off','Color',[1 1 1]);
    set(fCNR_gmMax_Contrast,'Position',[ 461   474   917   369],'Color',[1 1 1])
    
    subplot(121)
    imagesc(log(reshape(CNR_SE_GM_LS(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2))))
    hold on
    grid on
    title('logscale - CNR for GM LS contrast','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([valuelogMinContrast valuelogMaxContrast])
    colorbar
    
    subplot(122)
    imagesc(log(reshape(CNR_SE_WM_LS(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2))))
    hold on
    grid on
    title('logscale - CNR for WM LS contrast','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([valuelogMinContrast valuelogMaxContrast])
    colorbar  
    
    % ... Plot Matrices of CNR for resolution ...
    fCNR_Res_Contrast = figure('name',['logscale - CNR (SE - ',plot_seq,') variation w/ resolution (Diffusion) - gm=', ...
        num2str(valuegmMax),' - Nslices=',num2str(valueNslice),' - b-value=',num2str(valueB)], ...
        'numbertitle','off','Color',[1 1 1]);
    set(fCNR_Res_Contrast,'Position',[ 461   474   917   369],'Color',[1 1 1])
    
    subplot(121)
    imagesc(log(reshape(CNR_SE_GM_LS(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2))))
    hold on
    grid on
    title('logscale - CNR for GM LS contrast','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([valuelogMinContrast valuelogMaxContrast])
    colorbar
    
    subplot(122)
    imagesc(log(reshape(CNR_SE_WM_LS(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2))))
    hold on
    grid on
    title('logscale - CNR for WM LS contrast')
%     ylabel(['B0 (T)'])
    xlabel(['Resolution (mm)'])
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([valuelogMinContrast valuelogMaxContrast])   
    colorbar
    
else         % no logscale
    % ... Plot Matrices of CNR for gmax & TE ...
    fCNR_gmMax_Contrast = figure('name',['CNR (SE - ',plot_seq,') variation w/ Gm Max (Diffusion) - b_value=', ...
        num2str(valueB),' - Nslices=',num2str(valueNslice), ...
        ' - res=',num2str(valueRes)], ...
        'numbertitle','off','Color',[1 1 1]);
    set(fCNR_gmMax_Contrast,'Position',[ 461   474   917   369],'Color',[1 1 1])
    
    subplot(121)
    imagesc(reshape(CNR_SE_GM_LS(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2)))
    hold on
    grid on
    title('CNR for GM LS contrast','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([0 valueMaxContrast])
    colorbar
    
    subplot(122)
    imagesc(reshape(CNR_SE_WM_LS(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2)))
    hold on
    grid on
    title('CNR for WM LS contrast','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([0 valueMaxContrast])
    colorbar
    
    % ... Plot Matrices of CNR for resolution ...
    fCNR_Res_Contrast = figure('name',['CNR (SE - ',plot_seq,') variation w/ resolution (Diffusion) - gm=', ...
        num2str(valuegmMax),' - Nslices=',num2str(valueNslice),' - b-value=',num2str(valueB)], ...
        'numbertitle','off','Color',[1 1 1]);
    set(fCNR_Res_Contrast,'Position',[ 461   474   917   369],'Color',[1 1 1])
    
    subplot(121)
    imagesc(reshape(CNR_SE_GM_LS(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2)))
    hold on
    grid on
    title('CNR for GM LS contrast','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([0 valueMaxContrast])
    colorbar
    
    subplot(122)
    imagesc(reshape(CNR_SE_WM_LS(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2)))
    hold on
    grid on
    title('CNR for WM LS contrast','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([0 valueMaxContrast])
    colorbar
end

% ... Plot Values of CNR contrast for b-value ...
fCNR_bValue_Contrast = figure('name',['CNR (SE - ',plot_seq,') matrix variation w/ b_values (Diffusion) - gm=', ...
    num2str(valuegmMax),' - B0=', num2str(valueB0),' - Nslices=',num2str(valueNslice),' - res=',num2str(valueRes)], ...
    'numbertitle','off','Color',[1 1 1]);
set(fCNR_bValue_Contrast,'Position',[ 461   474   917   369],'Color',[1 1 1])

% subplot(121)
plot(b_vector,CNR_SE_GM_LS(plotB0,:,plotgmMax,plotRes),'Color','c')
ylabel([' CNR (au)'],'fontsize',20)
xlabel(['b_{value}'],'fontsize',20)
% title('CNR for GM LS contrast','fontsize',20)
set(gca,'FontSize',15)
hold on
% subplot(122)
plot(b_vector,CNR_SE_WM_LS(plotB0,:,plotgmMax,plotRes),'Color','b')
% ylabel([' CNR (au)'],'fontsize',20)
xlabel(['b_{value}'],'fontsize',20)
title('CNR variation with B-value, for each contrast','fontsize',20)
set(gca,'FontSize',15)
hold on
grid on
lgd = legend('GM-LS Contrast','WM-LS Contrast');
lgd.FontSize = 20;

%% 7 - Paper Figure 5: plot for SNR & CNR variation with b-value

if logscale == 1

    % Initialization of variables for the plots
    plotgmMax  = find(round(gmMaxVector,3)==round(valuegmMax,3));
    plotRes    = find(resNEW==valueRes);
    plotb      = find(b_vector==valueB);
    plotB0     = find(B0==valueB0);
    
    valueNslice      = Z_size/valueRes;                % Number of Slices
    valueMaxCM       = SNR_refinit*percSNRinit;
    valueMaxCMlog    = log(max(max(SNR_SE_WM(:,plotb,:,plotRes))));
    valueMinCMlog    = log(min(min(SNR_SE_WM(:,plotb,:,plotRes))));
    valueMaxContrast = 2;
    valuelogMinContrast = log(min(min(CNR_SE_WM_LS(:,plotb,:,plotRes))));
    valuelogMaxContrast = log(max(max(CNR_SE_WM_LS(:,plotb,:,plotRes))));
    
    figure()
    
    yyaxis left
    plot(b_vector,log(SNR_SE_GM(plotB0,:,plotgmMax,plotRes)),'Color','r')
    set(gca,'FontSize',22)
    hold on
    
    yyaxis left
    plot(b_vector,log(SNR_SE_WM(plotB0,:,plotgmMax,plotRes)),'-','Color','b')
    set(gca,'FontSize',22)
    hold on
    
    yyaxis left
    plot(b_vector,log(SNR_SE_LS(plotB0,:,plotgmMax,plotRes)),'-','Color','g')
    set(gca,'FontSize',22)
    hold on
    
    ax = gca; % Get handle to current axes.
    ax.YColor = 'k'; % black
    
    yyaxis right
    plot(b_vector,log(CNR_SE_GM_LS(plotB0,:,plotgmMax,plotRes)),'--','Color','c')
    set(gca,'FontSize',22)
    hold on
    
    yyaxis right
    plot(b_vector,log(CNR_SE_WM_LS(plotB0,:,plotgmMax,plotRes)),'--','Color','m')
    set(gca,'FontSize',22)
    hold on
    
    % title('SNR & CNR variation with B-value','fontsize',20)
    yyaxis left
    ylabel(['SNR log scale (a.u.)'],'fontsize',20)
    xlabel(['b_{value}'],'fontsize',20)
    
    yyaxis right
    ylabel(['CNR log scale (a.u.)'],'fontsize',20)
    
    hold on
    
    grid on
    lgd = legend('SNR GM','SNR WM','SNR LS','CNR GM-LS','CNR WM-LS');
    lgd.FontSize = 20;set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    lgd.LineWidth = 4;
    
    ax = gca; % Get handle to current axes.
    ax.YColor = 'k'; % Black
    ax.XColor = 'k'; % Black
 
elseif logscale == 0
    
    figure()
    
    yyaxis left
    plot(b_vector,SNR_SE_GM(plotB0,:,plotgmMax,plotRes),'Color','r')
    set(gca,'FontSize',22)
    hold on
    
    yyaxis left
    plot(b_vector,SNR_SE_WM(plotB0,:,plotgmMax,plotRes),'-','Color','b')
    set(gca,'FontSize',22)
    hold on
    
    yyaxis left
    plot(b_vector,SNR_SE_LS(plotB0,:,plotgmMax,plotRes),'-','Color','g')
    set(gca,'FontSize',22)
    hold on
    
    ax = gca; % Get handle to current axes.
    ax.YColor = 'k'; % black
    
    yyaxis right
    plot(b_vector,CNR_SE_GM_LS(plotB0,:,plotgmMax,plotRes),'--','Color','c')
    set(gca,'FontSize',22)
    hold on
    
    yyaxis right
    plot(b_vector,CNR_SE_WM_LS(plotB0,:,plotgmMax,plotRes),'--','Color','m')
    set(gca,'FontSize',22)
    hold on
    
    % title('SNR & CNR variation with B-value','fontsize',20)
    yyaxis left
    ylabel(['SNR scale (dB)'],'fontsize',20)
    xlabel(['b_{value}'],'fontsize',20)
    
    yyaxis right
    ylabel(['CNR scale (dB)'],'fontsize',20)
    
    hold on
    
    grid on
    lgd = legend('SNR GM','SNR WM','SNR LS','CNR GM-LS','CNR WM-LS');
    lgd.FontSize = 20;set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    lgd.LineWidth = 4;
    
    ax = gca; % Get handle to current axes.
    ax.YColor = 'k'; % Black
    ax.XColor = 'k'; % Black
    
end

%% 8 - Paper Figure 5: Plots WM - SNR & CNR (WM-LS) - gm, res

if logscale == 1
    % SNR WM - gm
    figure()
    subplot(121)
    imagesc(log(reshape(SNR_SE_WM(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2))))
    hold on
    grid on
%     title('logscale - SNR','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([valueMinCMlog valueMaxCMlog])
    colorbar
    
    hold on
    plot(plotgmMax,plotB0, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    
    % SNR WM - res
%     figure()
    subplot(122)
    imagesc(log(reshape(SNR_SE_WM(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2))))
    hold on
    grid on
%     title('logscale - SNR','fontsize',20)
%     ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['Resolution (mm)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([valueMinCMlog valueMaxCMlog])
    colorbar
    
    hold on
    plot(plotRes,plotB0, 'r+', 'MarkerSize', 20, 'LineWidth', 2);   
    
    % CNR WM_LS - gm
    figure()
    subplot(121)
    imagesc(log(reshape(CNR_SE_WM_LS(:,plotb,:,plotRes),size(B0,2),size(gmMaxVector,2))))
    hold on
    grid on
%     title('logscale - CNR for WM LS contrast','fontsize',20)
    ylabel(['B0 (T)'],'fontsize',20)
    xlabel(['GM max Value (mT/m)'],'fontsize',20)
    set(gca,'FontSize',15)
    yticks(1:1:size(B0,2)*2)
    xticks(1:1:size(gmMaxVector,2)*2)
    yticklabels(B0)
    xticklabels(gmMaxVector*10^3)
    caxis([valuelogMinContrast valuelogMaxContrast])
    colorbar
    hold on
    plot(plotgmMax,plotB0, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    
    % CNR WM_LS - res
%     figure()
    subplot(122)
    imagesc(log(reshape(CNR_SE_WM_LS(1:sizeB0res,plotb,plotgmMax,:),sizeB0res,size(resNEW,2))))
    hold on
    grid on
%     title('logscale - CNR for WM LS contrast')
    %     ylabel(['B0 (T)'])
    xlabel(['Resolution (mm)'])
    set(gca,'FontSize',15)
    yticks(1:1:sizeB0res*2)
    xticks(1:1:size(resNEW,2)*2)
    yticklabels(B0(1:sizeB0res))
    xticklabels(resNEW)
    caxis([valuelogMinContrast valuelogMaxContrast])
    colorbar
    hold on
    plot(plotRes,plotB0, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    
end
