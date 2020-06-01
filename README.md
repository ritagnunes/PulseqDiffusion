<p align="center">
<img src="Logo.png"/>
</p>

# PulseqDiffusion:  Diffusion-Weighted Echo Planar Imaging sequence using the Open Source PyPulseq

Diffusion-weighted Imaging (DWI) is an important sequence for many clinical applications such as stroke and tumor characterization [[1]](#references).  Though several vendor-neutral post-processing tools are available for diffusion, open-source acquisition implementations have not been provided for research purposes.

We propose the development of a tool `PulseqDiffusion` cross-vendor, open-source package of a multi-slice single-shot spin-echo echo-planar imaging (EPI)-based diffusion pulse sequence using `PyPulseq`[[2]](#references) - [[4]](#references), which can be extended to support multiple b-values and directions. We demonstrate this on (i) in vitro phantom to measure the Apparent Diffusion Coefficient  (ii) in vivo human brain data to obtain good quality Fractional Anisotropy (FA) and Mean Diffusivity maps. We provide basic reconstrucion software implemented in Matlab, for generating the images from the raw k-space data (main Matlab script: process_data.m). Example k-space data acquired at the University of Columbia is provided in the folder Example_Data (Phantom data: 3 directions, 5 b-values, 3 slices; In-Vivo data: 12 directions, 1 b-value, 20 slices and 3 directions, 3 b-values, 20 slices).

The images can then be processed utilizing freely available post-processing tools for generating quantitative diffusion maps. Example code using the image analysis software FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) [[5]](#references) is included in this package.

A tool for predicting the signal-to-noise ratio (SNR) observed for different brain tissues (gray matter, white matter, cerebrospinal fluid - CSF) and the contrast-to-noise ratio (CNR) of an acute stroke lesion relative to those tissues can be found in `jMRI_Publication_Functions`. The implemented code predicts the SNR and CNR per time unit for a spin-echo diffusion-weighted sequence, taking into account the steady-state value for the longitudinal magnetization. For that purpose, we adapted an expression used to predict the SNR per tissue per time unit for a spoiled gradient echo sequence [[6]](#references). This simulation considers the possibility to use either EPI or spiral readouts and considers the impact of: B0, max gradient amplitude, spatial resolution and b-value (comparing the achieved SNR with that of a typical 1.5T clinical scanner).

Software dependencies:
Matlab toolbox for writing and reading Nifti images: Nifti_tools [[7]](#references)
FSL - Fmrib's Software Library [[5]](#references)

---
## [Citations][scholar-citations]

1. Schellinger, P. D., R. N. Bryan, L. R. Caplan, J. A. Detre, R. R. Edelman, C. Jaigobin, C. S. Kidwell et al. "Evidence-based guideline: the role of diffusion and perfusion MRI for the diagnosis of acute ischemic stroke: report of the Therapeutics and Technology Assessment Subcommittee of the AmericanAcademy of Neurology." Neurology 75, no. (2010): 177-185.
2. Ravi, K., Geethanath, S., and John Thomas Vaughan Jr. PyPulseq: A Python Package for MRI Pulse Sequence Design. 10.21105/joss.01725
3. Ravi, K., Potdar, S., Poojar. P, Reddy, A., Kroboth, S., Nielsen, J., Zaitsev, M., Venkatesan, R & Geethanath, S. (2018). Pulseq-Graphical Programming Interface: Open source visual environment for prototyping pulse sequences and integrated magnetic resonance imaging algorithm development. Magnetic resonance imaging 52: 9-15.
4. Nunes, R.G., Ravi, K.S.,  Geethanath, S., Vaughan Jr, J.T. (2020). “Implementation of a Diffusion-Weighted Echo Planar Imaging sequence using the Open Source Hardware-Independent PyPulseq Tool”, 28th Annual Meeting of the International Society for Magnetic Resonance in Medicine, Sidney, Australia.
5. Jenkinson, M., Beckmann, C.F., Behrens, T.E., Woolrich, M.W., Smith, S.M. (2012). "FSL". NeuroImage, 62:782-90.
6. Fernandes, T.T., Gobul, M., Freitas, A.F., Geethanath, S., Nunes, R.G., (2020). "Diffusion-weighted imaging for acute stroke detection at low field MR - a feasibility study”, 28th Annual Meeting of the International Society for Magnetic Resonance in Medicine, Sidney, Australia.
7. https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image 
