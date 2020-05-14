<p align="center">
<img src="Logo.png"/>
</p>

# PulseqDiffusion:  Diffusion-Weighted Echo Planar Imaging sequence using the Open Source PyPulseq

Diffusion-weighted Imaging (DWI) is an important sequence for many clinical applications such as stroke and tumor characterization [[1]](#references).  Though several vendor-neutral post-processing tools are available for diffusion, open-source acquisition implementations have not been provided for research purposes.

We propose the development of a tool `PulseqDiffusion` cross-vendor, open-source package of a multi-slice single-shot spin echo-planar imaging-based diffusion pulse sequence using `PyPulseq`[[2]](#references) - [[4]](#references), which can be extended to support multiple b values and directions. We demonstrate this on (i) in vitro phantom to measure the Apparent Diffusion Coefficient values  (ii) in vivo human brain data to obtained good quality Fractional Anisotropy (FA) and Mean Diffusivity maps. The data can be processed utilizing the freely available post-processing tools and generate quantitative diffusion maps.

A SNR and CNR tool can be found in `jMRI_Publication_Functions` for predict the noise impact per tissue per time unit for a spoiled gradient echo sequence according to a set of parameters [[5]](#references). This focus on EPI and Spiral sequences and, allows to understand the impact of: B0, max gradient, resolution and b-value in the noise (comparing it with typical SNR for a clinical scan at 1.5T).

---
## [Citations][scholar-citations]

1. Schellinger, P. D., R. N. Bryan, L. R. Caplan, J. A. Detre, R. R. Edelman, C. Jaigobin, C. S. Kidwell et al. "Evidence-based guideline: the role of diffusion and perfusion MRI for the diagnosis of acute ischemic stroke: report of the Therapeutics and Technology Assessment Subcommittee of the AmericanAcademy of Neurology." Neurology 75, no. (2010): 177-185.
2. Ravi, K., Geethanath, S., and John Thomas Vaughan Jr. PyPulseq: A Python Package for MRI Pulse Sequence Design. 10.21105/joss.01725
3. Ravi, K., Potdar, S., Poojar. P, Reddy, A., Kroboth, S., Nielsen, J., Zaitsev, M., Venkatesan, R & Geethanath, S. (2018). Pulseq-Graphical Programming Interface: Open source visual environment for prototyping pulse sequences and integrated magnetic resonance imaging algorithm development. Magnetic resonance imaging 52: 9-15.
4. Nunes, R.G., Ravi, K.S.,  Geethanath, S., Vaughan Jr, J.T. (2020). “Implementation of a Diffusion-Weighted Echo Planar Imaging sequence using the Open Source Hardware-Independent PyPulseq Tool”, 28th Annual Meeting of the International Society for Magnetic Resonance in Medicine, Sidney, Australia.
5. Fernandes, T.T., Gobul, M., Freitas, A.F., Geethanath, S., Nunes, R.G., (2020). "Diffusion-weighted imaging for acute stroke detection at low field MR - a feasibility study”, 28th Annual Meeting of the International Society for Magnetic Resonance in Medicine, Sidney, Australia.
