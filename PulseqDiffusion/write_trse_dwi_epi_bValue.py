"""
Oscar Pena-Nogales(opennog(at)lpi(dot)tel(dot)uva(dot)es)
Rita G. Nunes
January 2020
Laboratorio de Procesado de Imagen, Universidad de Valladolid, Valladolid, Spain
LASEEB, Instituto Superior Técnico, Lisbon, Portugal

This function builds a twice refocused spin echo EPI diffusion weighting sequence optimizing the TE for a desired b-value.
All diffusion time is used to achieve larger b-values.
The b-value is computed according to the b-value integral definition.
"""

import math
import os

import matplotlib.pyplot as plt
import numpy as np
from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.calc_rf_center import calc_rf_center
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts

import PulseqDiffusion.diff_funcs as difunc

# Create new Sequence Object
seq = Sequence()

# Note in this code we do not work on integer times 'n'.

# Acquisition Parameters
# FOV and resolution
fov = 240e-3  # [m]
Nx = 64
Ny = 64
slice_thickness = 2.5e-3  # [m]
n_slices = 5

# Partial Fourier
pF = 1  # 0.75
Nyeff = int(pF * Ny)  # Number of Ny samples acquired
if pF is not 1:
    pF_str = "_" + str(pF) + "pF"
else:
    pF_str = ""

# b-value
bvalue = np.array([500])  # i.e., [100 200 500] in [s/mm2] - bvalue = 0 is always included.

# Spin-Echo parameters
TR = 5000e-3  # [s]

# Save seq
save_flag = False

# Fat saturation
fatsat_enable = 0
# pe_enable = 1
if fatsat_enable:
    fatsat_str = "_fatsat"
else:
    fatsat_str = ""

# Plot sequence and k-space trajectory
seqplot = 0
kplot = 0
waveplot = 0

# b-value parameters
nbvals = np.shape(bvalue)[0]
ndirs = 6

# Gradient Scaling
gscl = np.zeros(nbvals + 1)
gscl[1:] = np.sqrt(bvalue / np.max(bvalue))
gdir, nb0s = difunc.get_dirs(ndirs)

# Set system limits
system = Opts(max_grad=32, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', rf_ringdown_time=20e-6,
              rf_dead_time=100e-6, adc_dead_time=10e-6)

# Fat saturation
if fatsat_enable:
    b0 = 1.494
    sat_ppm = -3.45
    sat_freq = sat_ppm * 1e-6 * b0 * system.gamma
    rf_fs, _, _ = make_gauss_pulse(flip_angle=110 * math.pi / 180, system=system, duration=8e-3,
                                   bandwidth=abs(sat_freq),
                                   freq_offset=sat_freq)
    gz_fs = make_trapezoid(channel='z', system=system, delay=calc_duration(rf_fs), area=1 / 1e-4)

# Create 90 degree slice selection pulse and gradient
rf, gz, _ = make_sinc_pulse(flip_angle=math.pi / 2, system=system, duration=3e-3, slice_thickness=slice_thickness,
                            apodization=0.5, time_bw_product=4)

# Refocusing pulse with spoiling gradients.
# Note that the RF180 can be the same with the same phase (as long as it agrees with CPMG condition)
# Also, the spoiling gradients must be different to avoid un-wanted echos:
rf180, gz180, _ = make_sinc_pulse(flip_angle=math.pi, system=system, duration=5e-3, slice_thickness=slice_thickness,
                                  apodization=0.5, time_bw_product=4)
rf180.phase_offset = math.pi / 2
gz_spoil_1 = make_trapezoid(channel='z', system=system, area=6 / slice_thickness, duration=3e-3)
gz_spoil_2 = make_trapezoid(channel='z', system=system, area=6 * 2 / slice_thickness, duration=2 * 3e-3)

# Define other gradients and ADC events
delta_k = 1 / fov
k_width = Nx * delta_k
dwell_time = seq.grad_raster_time  # Full receiver bandwith
readout_time = Nx * dwell_time  # T_acq (acquisition time)
flat_time = math.ceil(readout_time / seq.grad_raster_time) * seq.grad_raster_time
gx = make_trapezoid(channel='x', system=system, amplitude=k_width / readout_time, flat_time=flat_time)
adc = make_adc(num_samples=Nx, duration=readout_time,
               delay=gx.rise_time + flat_time / 2 - (readout_time - dwell_time) / 2)

# Pre-phasing gradients
pre_time = 1e-3
gx_pre = make_trapezoid(channel='x', system=system, area=-gx.area / 2, duration=pre_time)
gz_reph = make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)
gy_pre = make_trapezoid(channel='y', system=system, area=-(Ny / 2 - 0.5 - (Ny - Nyeff)) * delta_k, duration=pre_time)

# Phase blip in shortest possible time
gy = make_trapezoid(channel='y', system=system, area=delta_k)
dur = math.ceil(calc_duration(gy) / seq.grad_raster_time) * seq.grad_raster_time

# Calculate some times constant throughout the process
# The time(gy) refers to the number of blips, thus we substract 0.5 since the number of lines is always even.
# The time(gx) refers to the time needed to read each line of the k-space. Thus, if Ny is even, it would take half of the lines plus another half.
duration_center = (
            calc_duration(gx) * (Ny / 2 + 0.5 - (Ny - Nyeff)) + dur * (Ny / 2 - 0.5 - (Ny - Nyeff)) + calc_duration(
        gx_pre, gy_pre))
rf_center_with_delay = rf.delay + calc_rf_center(rf)[0]
rf180_center_with_delay = rf180.delay + calc_rf_center(rf180)[0]

# Group variables
seq_sys_Dict = {"seq": seq,
                "system": system}

grads_times_Dict = {"rf180": rf180,
                    "rf180_center_with_delay": rf180_center_with_delay,
                    "rf_center_with_delay": rf_center_with_delay,
                    "gz_spoil_1": gz_spoil_1,
                    "gz_spoil_2": gz_spoil_2,
                    "gz": gz,
                    "duration_center": duration_center,
                    "pre_time": pre_time}

bvalue_Dict = {"bvalue": bvalue,
               "nbvals": nbvals,
               "gscl": gscl}

# Optimize TE for the desired b-value for the TRSE sequence
TE, bvalue_tmp, waveform, gscl_max, d1, d2, d3, d4 = difunc.opt_TE_bv_TRSE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)

# Crusher gradients
gx_crush = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
gz_crush = make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

# TR delay - Takes everything into account
# EPI reading time:
# Distance between the center of the RF90s must be TR
EPI_time = calc_duration(gx) * Nyeff + calc_duration(gy) * (Nyeff - 1)
tr_per_slice = TR / n_slices
if fatsat_enable:
    tr_delay = math.floor((tr_per_slice - (
                TE - duration_center + EPI_time + pre_time) - rf_center_with_delay - calc_duration(gx_crush, gz_crush) \
                           - calc_duration(rf_fs, gz_fs)) \
                          / seq.grad_raster_time) * seq.grad_raster_time
else:
    tr_delay = math.floor((tr_per_slice - (
                TE - duration_center + EPI_time + pre_time) - rf_center_with_delay - calc_duration(gx_crush, gz_crush)) \
                          / seq.grad_raster_time) * seq.grad_raster_time

# Check TR delay time
assert tr_delay > 0, "Such parameter configuration needs longer TR."

""" EPI calibration """
for s in range(n_slices):
    # Fat saturation
    if fatsat_enable:
        seq.add_block(rf_fs, gz_fs)

    # RF90
    rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
    seq.add_block(rf, gz)
    seq.add_block(gz_reph)

    # Delay for first RF180
    seq.add_block(make_delay(d1))

    # First RF180
    seq.add_block(gz_spoil_1)
    rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
    seq.add_block(rf180, gz180)
    seq.add_block(gz_spoil_1)

    # Delay for second RF180
    seq.add_block(make_delay(d2 + d3))

    # Second RF180
    seq.add_block(gz_spoil_2)
    rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
    seq.add_block(rf180, gz180)
    seq.add_block(gz_spoil_2)

    # Delay for EPI
    seq.add_block(make_delay(d4))

    # Locate k-space - Only on the frequency encoding direction.
    seq.add_block(gx_pre)

    for i in range(Nyeff):
        seq.add_block(gx, adc)  # Read one line of k-space
        if i is not Nyeff - 1:
            seq.add_block(make_delay(calc_duration(gy)))  # Imitate EPI reading but reading always the same line.
        gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient

    seq.add_block(gx_crush, gz_crush)

    # Wait TR
    if tr_delay > 0:
        seq.add_block(make_delay(tr_delay))

""" b-zero acquisition """
for d in range(nb0s):
    for s in range(n_slices):
        # Fat saturation
        if fatsat_enable:
            seq.add_block(rf_fs, gz_fs)

        # RF90
        rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf, gz)
        seq.add_block(gz_reph)

        # Delay for first RF180
        seq.add_block(make_delay(d1))

        # First RF180
        seq.add_block(gz_spoil_1)
        rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf180, gz180)
        seq.add_block(gz_spoil_1)

        # Delay for second RF180
        seq.add_block(make_delay(d2 + d3))

        # Second RF180
        seq.add_block(gz_spoil_2)
        rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf180, gz180)
        seq.add_block(gz_spoil_2)

        # Delay for EPI
        seq.add_block(make_delay(d4))

        # Locate k-space
        seq.add_block(gx_pre, gy_pre)

        for i in range(Nyeff):
            seq.add_block(gx, adc)  # Read one line of k-space
            if i is not Nyeff - 1:
                seq.add_block(gy)  # Phase blip
            gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient

        seq.add_block(gx_crush, gz_crush)

        # Wait TR
        if tr_delay > 0:
            seq.add_block(make_delay(tr_delay))

""" DWI acquisition """
for bv in range(1, nbvals + 1):
    for d in range(ndirs):
        for s in range(n_slices):
            # Fat saturation
            if fatsat_enable:
                seq.add_block(rf_fs, gz_fs)

            # RF90
            rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf, gz)
            seq.add_block(gz_reph)

            # Diffusion-weighting gradient d1
            gdiffx_d1 = make_trapezoid(channel='x', system=system,
                                       amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 0],
                                       duration=d1)
            gdiffy_d1 = make_trapezoid(channel='y', system=system,
                                       amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 1],
                                       duration=d1)
            gdiffz_d1 = make_trapezoid(channel='z', system=system,
                                       amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 2],
                                       duration=d1)

            seq.add_block(gdiffx_d1, gdiffy_d1, gdiffz_d1)

            # First RF180
            seq.add_block(gz_spoil_1)
            rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil_1)

            # Diffusion-weighting gradient d2
            gdiffx_d2 = make_trapezoid(channel='x', system=system,
                                       amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 0],
                                       duration=d2)
            gdiffy_d2 = make_trapezoid(channel='y', system=system,
                                       amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 1],
                                       duration=d2)
            gdiffz_d2 = make_trapezoid(channel='z', system=system,
                                       amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 2],
                                       duration=d2)

            seq.add_block(gdiffx_d2, gdiffy_d2, gdiffz_d2)

            # Diffusion-weighting gradient d3
            gdiffx_d3 = make_trapezoid(channel='x', system=system,
                                       amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 0],
                                       duration=d3)
            gdiffy_d3 = make_trapezoid(channel='y', system=system,
                                       amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 1],
                                       duration=d3)
            gdiffz_d3 = make_trapezoid(channel='z', system=system,
                                       amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 2],
                                       duration=d3)

            seq.add_block(gdiffx_d3, gdiffy_d3, gdiffz_d3)

            # Second RF180
            seq.add_block(gz_spoil_2)
            rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil_2)

            # Diffusion-weighting gradient d4
            gdiffx_d4 = make_trapezoid(channel='x', system=system,
                                       amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 0],
                                       duration=d4)
            gdiffy_d4 = make_trapezoid(channel='y', system=system,
                                       amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 1],
                                       duration=d4)
            gdiffz_d4 = make_trapezoid(channel='z', system=system,
                                       amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 2],
                                       duration=d4)

            seq.add_block(gdiffx_d4, gdiffy_d4, gdiffz_d4)

            # Locate k-space
            seq.add_block(gx_pre, gy_pre)

            for i in range(Nyeff):
                seq.add_block(gx, adc)  # Read one line of k-space
                if i is not Nyeff - 1:
                    seq.add_block(gy)  # Phase blip
                gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient

            seq.add_block(gx_crush, gz_crush)

            # Wait TR
            if tr_delay > 0:
                seq.add_block(make_delay(tr_delay))

if seqplot:
    seq.plot()

if kplot:
    ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

    time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
    plt.plot(time_axis, ktraj.T)
    plt.plot(t_adc, ktraj_adc[0, :], '.')
    plt.figure()
    plt.plot(ktraj[0, :], ktraj[1, :], 'b')
    plt.axis('equal')
    plt.plot(ktraj_adc[0, :], ktraj_adc[1, :], 'r.')
    plt.show()

if waveplot:
    plt.figure()
    plt.plot(waveform)

# Save the sequence
if save_flag:
    seqfname = "trse_dwi_" + str(n_slices) + "slices_" + str(bvalue) + "bvalues_" + str(ndirs) + "dirs_" + str(
        d4 * 1e3) + "d4_" + \
               str(round(TE * 1e3, 2)) + "TE_" + str(round(TR * 1e3)) + "TR_" + str(
        system.max_grad / system.gamma * 1e3) + "G_Max_" + \
               str(system.max_slew / system.gamma) + "SR_Max" + fatsat_str + pF_str
    os.mkdir("tests/" + seqfname)
    seq.write("tests/" + seqfname + "/" + seqfname + ".seq")
    print("Seq file saved --", seqfname)
