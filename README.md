<p align="center">
<img src="logo.png"/>
</p>

# PulseqDiffusion:  Diffusion-Weighted Echo Planar Imaging sequence using the Open Source PyPulseq

Diffusion-weighted Imaging (DWI) is an important sequence for many clinical applications such as stroke and tumor characterization [[1]](#references).  Though several vendor-neutral post-processing tools are available for diffusion, open-source acquisition implementations have not been provided for research purposes.

We propose the development of a tool `PulseqDiffusion` cross-vendor, open-source package of a multi-slice single-shot spin echo-planar imaging-based diffusion pulse sequence using `PyPulseq`, which can be extended to support multiple b values and directions. We demonstrate this on (i) in vitro phantom to measure the Apparent Diffusion Coefficient values  (ii) in vivo human brain data to obtained good quality Fractional Anisotropy (FA) and Mean Diffusivity maps. The data can be processed utilizing the freely available post-processing tools and generate quantitative diffusion maps.

---
## [Citations][scholar-citations]

1. Schellinger, P. D., R. N. Bryan, L. R. Caplan, J. A. Detre, R. R. Edelman, C. Jaigobin, C. S. Kidwell et al. "Evidence-based guideline: the role of diffusion and perfusion MRI for the diagnosis of acute ischemic stroke: report of the Therapeutics and Technology Assessment Subcommittee of the AmericanAcademy of Neurology." Neurology 75, no. (2010): 177-185.
2. Ravi, K., Potdar, S., Poojar. P, Reddy, A., Kroboth, S., Nielsen, J., Zaitsev, M., Venkatesan, R & Geethanath, S. (2018). Pulseq-Graphical ProgrammingInterface: Open source visual environment for prototyping pulse sequences and integrated magnetic resonance imaging algorithm development. Magneticresonance imaging 52: 9-15
3. Ravi, K., Geethanath, S., and John Thomas Vaughan Jr. PyPulseq: A Python Package for MRI Pulse Sequence Design. 10.21105/joss.01725
4. Layton, K., Kroboth, S., Jia, F., Littin, S., Yu, H., Leupold, J., Nielsen, J., Stöcker, T., Zaitsev, M. (2017) Pulseq: a rapid and hardware‐independent pulse sequence prototyping framework. Magnetic resonance in medicine77(4):1544-1552.
