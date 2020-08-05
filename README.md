<p align="center">
<img src="Logo.png"/>
</p>

[<img title="PyPulseq Badge" src="https://img.shields.io/badge/made%20using-pypulseq-brightgreen" width="100">](https://github.com/imr-framework/pypulseq)

# PulseqDiffusion:  Diffusion-Weighted Echo Planar Imaging sequence using the Open Source PyPulseq

Diffusion-weighted Imaging (DWI) is an important sequence for many clinical applications such as stroke and tumor characterization [[1]](#references).  Though several vendor-neutral post-processing tools are available for diffusion, open-source acquisition implementations have not been provided for research purposes.

We propose the development of a tool `PulseqDiffusion` cross-vendor, open-source package of a multi-slice single-shot spin-echo echo-planar imaging (EPI)-based diffusion pulse sequence using `PyPulseq` [[2]](#references) - [[4]](#references), which can be extended to support multiple b-values and directions. We demonstrate this on (i) in vitro phantom to measure the Apparent Diffusion Coefficient  (ii) in vivo human brain data to obtain good quality Fractional Anisotropy (FA) and Mean Diffusivity (MD) maps. We provide basic reconstrucion software implemented in Matlab, for generating the images from the raw k-space data. Example data acquired at the University of Columbia is provided in the folder `Example_Data` (Phantom data: 3 directions, 5 b-values, 3 slices; In-Vivo data: 12 directions, 1 b-value, 20 slices and 3 directions, 3 b-values, 20 slices).

The images can then be processed utilizing freely available post-processing tools for generating quantitative diffusion maps. Example code using the image analysis software FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) [[5]](#references) is included in this package.

A tool for predicting the signal-to-noise ratio (SNR) observed for different brain tissues (gray matter, white matter, cerebrospinal fluid - CSF) and the contrast-to-noise ratio (CNR) of an acute stroke lesion relative to those tissues can be found in `SNR_CNR_Study`. The implemented code predicts the SNR and CNR per time unit for a spin-echo diffusion-weighted sequence, taking into account the steady-state value for the longitudinal magnetization. For that purpose, we adapted an expression used to predict the SNR per tissue per time unit for a spoiled gradient echo sequence [[6]](#references). This simulation considers the possibility to use either EPI or spiral readouts and considers the impact of: B0, max gradient amplitude, spatial resolution and b-value (comparing the achieved SNR with that of a typical 1.5T clinical scanner).

## Installation and Dependencies
For sequence optimization depending on available hardware (`SNR_CNR_Study`): Matlab (www.mathworks.com)

For sequence development: PyPulseq [[2]](#references) \>=Python 3.6, virtual environment recommended:
```pip install pypulseq```

To perform basic image reconstruction and processing:
Matlab toolbox for writing and reading Nifti images: Nifti_tools downloadable from [[7]](#references) to be placed within a folder called 'nifti_tools'

FSL - Fmrib's Software Library [[5]](#references)
FSL installation instructions can be found at:
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation

Note: raw MRI data readers are vendor specific and hence cannot be provided.

## Usage and Examples
1- The user should start by running one of the three Python scripts provided to generate a diffusion-weighted EPI sequence:
`write_se_dw_epi_bValue.py`, `write_se_dw_rsepi_bValue.py`, `write_trse_dw_epi_bValue.py` 
The first two use single spin-echo preparation modules, while in the latter a doubly-refocused diffusion preparation module is provided for reducing eddy-current distortions.

Each of these scripts contains a section called "Acquisition Parameters" which can be modified as desired (e.g. k-space matrix size - Nx, Ny; slice_thickness).
Importantly, the system limits should be set according to the used MRI scanner.
Basic diffusion schemes are provided in the auxiliary file `diff_funcs.py` (3, 6, 12 and 60 directions), which can easily be extended by the user.

A .seq file will be written which can be run in different MRI platforms. In order to be able to do so, users will need to contact the developers of [Pulseq](https://github.com/pulseq/pulseq) [[8]](#references) to obtain the required vendor-specific software interpreter and follow their instructions for uploading the sequence onto the scanner.

2 - Basic reconstruction Matlab code is provided for reconstructing the images from the exported raw data main Matlab script: `process_data.m`.
By specifying the name of the k-space data file provided in the [Example_Data](https://github.com/ritagnunes/PulseqDiffusion/tree/master/Example_Data) folder, the user will be able to reconstruct the corresponding images and estimate either the full diffusion tensor (using FSL's ´dtifit´) or directional apparent diffusion coefficient maps (by setting the "model" parameter). Note that this script will need to be consistent with the selected acquisition parameters (specified spatial resolution and slice thickness).
The script enables to use the FSL "eddy" pre-processing tool by setting a parameter with the same name to 1.

3 - To visualize the obtained Nifti format images, as well as the results of the diffusion tensor fit (FA. MD and the principal diffusion eigenvector), the FSL image viewer "FSLeyes" can be used.

## Code testing
For checking that the provided software is properly installed, the following tests should be run:
1- In Matlab: ´runTests_ReconAndProcessing´;
2- Using Python: ´running test_suite.py´

---
## [Citations][scholar-citations]

1. Schellinger, P. D., R. N. Bryan, L. R. Caplan, J. A. Detre, R. R. Edelman, C. Jaigobin, C. S. Kidwell et al. "Evidence-based guideline: the role of diffusion and perfusion MRI for the diagnosis of acute ischemic stroke: report of the Therapeutics and Technology Assessment Subcommittee of the AmericanAcademy of Neurology." Neurology 75, no. (2010): 177-185.
2. Ravi, K., Geethanath, S., and John Thomas Vaughan Jr. PyPulseq: A Python Package for MRI Pulse Sequence Design. 10.21105/joss.01725
3. Ravi, K., Potdar, S., Poojar. P, Reddy, A., Kroboth, S., Nielsen, J., Zaitsev, M., Venkatesan, R & Geethanath, S. (2018). Pulseq-Graphical Programming Interface: Open source visual environment for prototyping pulse sequences and integrated magnetic resonance imaging algorithm development. Magnetic resonance imaging 52: 9-15.
4. Nunes, R.G., Ravi, K.S.,  Geethanath, S., Vaughan Jr, J.T. (2020). “Implementation of a Diffusion-Weighted Echo Planar Imaging sequence using the Open Source Hardware-Independent PyPulseq Tool”, 28th Annual Meeting of the International Society for Magnetic Resonance in Medicine, Sidney, Australia.
5. Jenkinson, M., Beckmann, C.F., Behrens, T.E., Woolrich, M.W., Smith, S.M. (2012). "FSL". NeuroImage, 62:782-90.
6. Fernandes, T.T., Golub, M., Freitas, A.C., Geethanath, S., Nunes, R.G., (2020). "Diffusion-weighted imaging for acute stroke detection at low field MR - a feasibility study”, 28th Annual Meeting of the International Society for Magnetic Resonance in Medicine, Sidney, Australia.
7. https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
8. Layton, K., Kroboth, S., Jia, F., Littin, S., Yu, H, Leupold, J., Nielsen, J., Stöcker, T., Zaitsev, M. (2017). "Pulseq: A rapid and hardware-independent pulse sequence prototyping framework" Magnetic Resonance in Medicine, 77(4):1544-1552.
