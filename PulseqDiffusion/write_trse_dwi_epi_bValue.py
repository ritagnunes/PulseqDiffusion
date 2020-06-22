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
import sys
import numpy as np
import os
import diff_funcs as difunc

from pypulseq.make_adc import make_adc
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts
from pypulseq.calc_duration import calc_duration
from pypulseq.calc_rf_center import calc_rf_center
from pypulseq.make_delay import make_delay
from pypulseq.Sequence.sequence import Sequence

import matplotlib.pyplot as plt

# Create new Sequence Object
seq = Sequence()
# This code is old, and we do not work with integer times 'n', it would be convenient to update it.

# Acquisition Parameters
# FOV and resolution
fov = 240e-3 # [m]
Nx = 64
Ny = 64
slice_thickness = 2.5e-3 # [m]
n_slices = 5

# Partial Fourier
pF = 1 # 0.75
Nyeff = int(pF*Ny) # Number of Ny samples acquired
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
gscl = np.zeros(nbvals+1)
gscl[1:] = np.sqrt(bvalue/np.max(bvalue))
gdir, nb0s = difunc.get_dirs(ndirs)

# Set system limits
system = Opts(max_grad=32, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', rf_ringdown_time=20e-6,
              rf_dead_time=100e-6, adc_dead_time=10e-6)

# Fat saturation
if fatsat_enable:
    b0 = 1.494
    sat_ppm = -3.45
    sat_freq = sat_ppm * 1e-6 * b0 * system.gamma
    rf_fs, _, _ = make_gauss_pulse(flip_angle=110 * math.pi / 180, system=system, duration=8e-3, bandwidth=abs(sat_freq),
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
rf180.phase_offset = math.pi/2
gz_spoil_1 = make_trapezoid(channel='z', system=system, area=6/slice_thickness, duration=3e-3)
gz_spoil_2 = make_trapezoid(channel='z', system=system, area=6*2/slice_thickness, duration=2*3e-3)

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
gy_pre = make_trapezoid(channel='y', system=system, area=-(Ny / 2 - 0.5 - (Ny - Nyeff)) * delta_k, duration=pre_time) # Es -0.5 y no +0.5 porque hay que pensar en areas, no en rayas!

# Phase blip in shortest possible time
gy = make_trapezoid(channel='y', system=system, area=delta_k)
dur = math.ceil(calc_duration(gy) / seq.grad_raster_time) * seq.grad_raster_time

"""" Obtain TE and diffusion-weighting gradient waveform """
# For S&T monopolar waveforms
# From an initial TE, check we satisfy all constraints -> otherwise increase TE.
# Once all constraints are okay -> check b-value, if it is lower than the target one -> increase TE
# Finally, with TRSE low b-values can not be acquired, thus proper scaling is needed.
# Looks time-inefficient but it is fast enough to make it user-friendly.
# TODO: Add ramps to the waveforms to compure the b-values more accurately.
# TODO: Re-scale the waveform to the exact b-value because increasing the TE might produce slightly higher ones.

# Calculate some times constant throughout the process
# The time(gy) refers to the number of blips, thus we substract 0.5 since the number of lines is always even.
# The time(gx) refers to the time needed to read each line of the k-space. Thus, if Ny is even, it would take half of the lines plus another half.
duration_center = (calc_duration(gx)*(Ny/2 + 0.5 - (Ny - Nyeff)) + dur * (Ny/2 - 0.5 - (Ny - Nyeff)) + calc_duration(gx_pre, gy_pre))
rf_center_with_delay = rf.delay + calc_rf_center(rf)[0]
rf180_center_with_delay = rf180.delay + calc_rf_center(rf180)[0]

# Find minimum TE considering the readout times and RF + spoil gradient durations
TE = 80e-3  # [s]
delay_tela = -1 # Large TE/2
while delay_tela <= 0:
    TE = TE + 0.02e-3  # [ms]
    delay_tela = math.ceil((TE / 2 + (- calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_2) - \
                       duration_center ) + (- calc_duration(gz) + rf_center_with_delay - \
                       pre_time - calc_duration(gz_spoil_1) - rf180_center_with_delay)) / seq.grad_raster_time) * seq.grad_raster_time

# Sice there are 3 equations and 4 parameters (TGReese2002 - https://doi.org/10.1002/mrm.10308)
# One parameter can be tuned, while the other 3 are then fixed. This is why we fix d4.
d4 = 3e-3  # [ms]

# Find minimum TE for the target d4
# Waveform Ramp time
gdiff_rt = math.ceil(system.max_grad / system.max_slew / seq.grad_raster_time) * seq.grad_raster_time
d1 = -1
while d1 <= 2 * gdiff_rt: # Include this condition to have trapepoids everywhere
    TE = TE + 2 * seq.grad_raster_time  # [ms]
    # Large TE/2
    delay_tela = math.ceil((TE / 2 + (- calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_2) - \
                          duration_center) + (- calc_duration(gz) + rf_center_with_delay - \
                        pre_time - calc_duration(gz_spoil_1) - rf180_center_with_delay)) / seq.grad_raster_time) * seq.grad_raster_time
    # Short TE/2 (Time between RF180s)
    delay_tes = math.ceil((TE / 2 - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_1) - \
                           calc_duration(gz_spoil_2) - rf180_center_with_delay) / seq.grad_raster_time) * seq.grad_raster_time

    d1 = math.ceil((delay_tela - d4) / seq.grad_raster_time) * seq.grad_raster_time

# Find minimum TE for the target b-value
bvalue_tmp = 0
# We initially substract 2 times the raster time because we don't know if will achieve the desired b-value from the beginning.
# Because we would not need to increase the TE.
# Note that for low b-values, the TE is the same for all of them due to the large number of elements of this sequence.
TE = TE - 2*seq.grad_raster_time
while bvalue_tmp < np.max(bvalue):
    # The following scalar multiplication is just to increase speed
    TE = TE + int(math.ceil(np.max(bvalue)/500)) * 4 * seq.grad_raster_time  # [ms]

    # Large TE/2
    delay_tela = math.ceil((TE / 2 + (- calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_2) - \
                          duration_center) + (- calc_duration(gz) + rf_center_with_delay - \
                        pre_time - calc_duration(gz_spoil_1) - rf180_center_with_delay)) / seq.grad_raster_time) * seq.grad_raster_time
    # Short TE/2 (Time between RF180s)
    delay_tes = math.ceil((TE / 2 - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_1) - \
                           calc_duration(gz_spoil_2) - rf180_center_with_delay) / seq.grad_raster_time) * seq.grad_raster_time

    d1 = math.ceil((delay_tela - d4) / seq.grad_raster_time) * seq.grad_raster_time
    d3 = math.ceil(((d1 + delay_tes - d4) / 2) / seq.grad_raster_time) * seq.grad_raster_time
    d2 = math.ceil((delay_tes - d3) / seq.grad_raster_time) * seq.grad_raster_time

    # Due to the complexity of the sequence and that I have not found its b-value, I implement it by its definition.
    # However, for simplicity I compute them as ideal rectangular pulses.

    # d1 start after the RF90
    n1 = int(math.ceil((calc_duration(gz) - rf_center_with_delay + pre_time + seq.grad_raster_time) / seq.grad_raster_time) * seq.grad_raster_time/seq.grad_raster_time)
    nd1 = int(math.ceil(d1 / seq.grad_raster_time) * seq.grad_raster_time/seq.grad_raster_time)

    # d2 start after the first RF180
    n2 = int(math.ceil(((n1 + nd1)*seq.grad_raster_time + 2 * calc_duration(gz_spoil_1) + calc_duration(rf180) + seq.grad_raster_time ) / seq.grad_raster_time) * seq.grad_raster_time/seq.grad_raster_time)
    nd2 = int(math.ceil(d2 / seq.grad_raster_time) * seq.grad_raster_time/seq.grad_raster_time)

    # d3 starts after d2
    n3 = int(math.ceil(((n2 + nd2)*seq.grad_raster_time + seq.grad_raster_time) / seq.grad_raster_time) * seq.grad_raster_time/seq.grad_raster_time)
    nd3 = int(math.ceil(d3 / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)

    # d4 starts after the second RF180
    n4 = int(math.ceil(((n2 + nd2 + nd3)*seq.grad_raster_time + 2 * calc_duration(gz_spoil_2) + calc_duration(rf180) + seq.grad_raster_time) / seq.grad_raster_time) * seq.grad_raster_time/seq.grad_raster_time)
    nd4 = int(math.ceil(d4 / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)

    # Due to the large TE of TRSE and the high resolution we might need to shrink the following array.
    # To speed up computation we calculate the b-value with short_f times lower of the time resolution.
    # Shortening factor:
    short_f = 5
    n = int(np.ceil(TE / seq.grad_raster_time)/short_f)
    waveform = np.zeros(n)

    # Compose waveform
    waveform[math.ceil(n1/short_f):math.floor((n1 + nd1 + 1)/short_f)] = system.max_grad
    waveform[math.ceil(n2/short_f):math.floor((n2 + nd2 + 1)/short_f)] = -system.max_grad
    waveform[math.ceil(n3/short_f):math.floor((n3 + nd3 + 1)/short_f)] = system.max_grad
    waveform[math.ceil(n4/short_f):math.floor((n4 + nd4 + 1)/short_f)] = -system.max_grad

    # Include ramps
    nrt = int(math.floor(math.floor(system.max_grad / system.max_slew / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time / short_f))
    ramp_values = np.linspace(1, nrt, nrt) * seq.grad_raster_time * system.max_slew * short_f
    # d1
    if n1 % short_f == 0:
        waveform[math.ceil(n1 / short_f):math.ceil(n1 / short_f) + nrt] = ramp_values
    else:
        waveform[math.ceil(n1/short_f):math.ceil((n1 + 1)/short_f) + nrt] = ramp_values

    waveform[math.ceil((n1 + nd1 - 1)/short_f - nrt):math.ceil((n1 + nd1 - 1) / short_f)] = np.flip(ramp_values)

    # d2
    if n2 % short_f == 0:
        waveform[math.ceil(n2 / short_f):math.ceil(n2 / short_f) + nrt] = -ramp_values
    else:
        waveform[math.ceil(n2 / short_f):math.ceil((n2 + 1) / short_f) + nrt] = -ramp_values

    waveform[math.ceil((n2 + nd2 - 1) / short_f - nrt):math.ceil((n2 + nd2 - 1) / short_f)] = -np.flip(ramp_values)

    # d3
    if n3 % short_f == 0:
        waveform[math.ceil(n3 / short_f):math.ceil(n3 / short_f) + nrt] = ramp_values
    else:
        waveform[math.ceil(n3 / short_f):math.ceil((n3 + 1) / short_f) + nrt] = ramp_values

    waveform[math.ceil((n3 + nd3 - 1) / short_f - nrt):math.ceil((n3 + nd3 - 1) / short_f)] = np.flip(ramp_values)

    # d4
    if n4 % short_f == 0:
        waveform[math.ceil(n4 / short_f):math.ceil(n4 / short_f) + nrt] = -ramp_values
    else:
        waveform[math.ceil(n4 / short_f):math.ceil((n4 + 1) / short_f) + nrt] = -ramp_values

    waveform[math.ceil((n4 + nd4 - 1) / short_f - nrt):math.ceil((n4 + nd4 - 1) / short_f)] = -np.flip(ramp_values)

    # Prepare Integral
    nRF180_1 = int(math.ceil((n2 * seq.grad_raster_time - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_1))/short_f/ seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
    nRF180_2 = int(math.ceil((n4 * seq.grad_raster_time - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_2))/short_f/ seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
    INV = np.ones(n)
    INV[nRF180_1:nRF180_2+1] = -1

    C = np.tril(np.ones([n, n]))
    C2 = np.matmul(np.transpose(C), C)

    F = waveform * INV * short_f * seq.grad_raster_time
    bv = (2*math.pi)**2 * np.matmul(np.matmul(np.transpose(F), C2), F) * short_f * seq.grad_raster_time
    bvalue_tmp = bv * 1e-6  # [s/mm2]

print("b-value low time resolution:", round(bvalue_tmp, 2), "s/mm2")
# To now the exact b-value - repeat the above with full time resolution
short_f = 1
n = int(np.ceil(TE / seq.grad_raster_time)/short_f)
waveform = np.zeros(n)

# Compose waveform
waveform[math.ceil(n1/short_f):math.floor((n1 + nd1 + 1)/short_f)] = system.max_grad
waveform[math.ceil(n2/short_f):math.floor((n2 + nd2 + 1)/short_f)] = -system.max_grad
waveform[math.ceil(n3/short_f):math.floor((n3 + nd3 + 1)/short_f)] = system.max_grad
waveform[math.ceil(n4/short_f):math.floor((n4 + nd4 + 1)/short_f)] = -system.max_grad

# Include ramps
nrt = int(math.floor(math.floor(system.max_grad / system.max_slew / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time / short_f))
ramp_values = np.linspace(0, nrt, nrt) * seq.grad_raster_time * system.max_slew * short_f
# d1
if n1 % short_f == 0:
    waveform[math.ceil(n1 / short_f):math.ceil(n1 / short_f) + nrt] = ramp_values
else:
    waveform[math.ceil(n1/short_f):math.ceil((n1 + 1)/short_f) + nrt] = ramp_values

if (n1 + nd1) % short_f == 0:
    waveform[math.ceil((n1 + nd1)/short_f - nrt) + 1:math.floor((n1 + nd1 + 1)/short_f)] = np.flip(ramp_values)
else:
    waveform[math.ceil((n1 + nd1) / short_f - nrt):math.floor((n1 + nd1 + 1) / short_f)] = np.flip(ramp_values)

# d2
if n2 % short_f == 0:
    waveform[math.ceil(n2 / short_f):math.ceil(n2 / short_f) + nrt] = -ramp_values
else:
    waveform[math.ceil(n2 / short_f):math.ceil((n2 + 1) / short_f) + nrt] = -ramp_values

if (n2 + nd2) % short_f == 0:
    waveform[math.ceil((n2 + nd2) / short_f - nrt) + 1:math.floor((n2 + nd2 + 1)/short_f)] = -np.flip(ramp_values)
else:
    waveform[math.ceil((n2 + nd2) / short_f - nrt):math.floor((n2 + nd2 + 1) / short_f)] = -np.flip(ramp_values)

# d3
if n3 % short_f == 0:
    waveform[math.ceil(n3 / short_f):math.ceil(n3 / short_f) + nrt] = ramp_values
else:
    waveform[math.ceil(n3 / short_f):math.ceil((n3 + 1) / short_f) + nrt] = ramp_values

if (n3 + nd3) % short_f == 0:
    waveform[math.ceil((n3 + nd3) / short_f - nrt) + 1:math.floor((n3 + nd3 + 1)/short_f)] = np.flip(ramp_values)
else:
    waveform[math.ceil((n3 + nd3) / short_f - nrt):math.floor((n3 + nd3 + 1) / short_f)] = np.flip(ramp_values)

# d4
if n4 % short_f == 0:
    waveform[math.ceil(n4 / short_f):math.ceil(n4 / short_f) + nrt] = -ramp_values
else:
    waveform[math.ceil(n4 / short_f):math.ceil((n4 + 1) / short_f) + nrt] = -ramp_values

if (n4 + nd4) % short_f == 0:
    waveform[math.ceil((n4 + nd4) / short_f - nrt) + 1:math.floor((n4 + nd4 + 1)/short_f)] = -np.flip(ramp_values)
else:
    waveform[math.ceil((n4 + nd4) / short_f - nrt):math.floor((n4 + nd4 + 1) / short_f)] = -np.flip(ramp_values)

# Prepare Integral
nRF180_1 = int(math.ceil((n2 * seq.grad_raster_time - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_1))/short_f/ seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
nRF180_2 = int(math.ceil((n4 * seq.grad_raster_time - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_2))/short_f/ seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
INV = np.ones(n)
INV[nRF180_1:nRF180_2+1] = -1

# These following lines might crush due to system overload - too large matrix
# In that case:
# 1) Reduce b-value to decrease TE
# 2) Omit the computation of the b-value with full time resolution.
# 3) Run code in a better computer
C = np.tril(np.ones([n, n]))
C2 = np.matmul(np.transpose(C), C)

F = waveform * INV * short_f * seq.grad_raster_time
bv = (2*math.pi)**2 * np.matmul(np.matmul(np.transpose(F), C2), F) * short_f * seq.grad_raster_time
bvalue_tmp = bv * 1e-6  # [s/mm2]
print("b-value high time resolution:", round(bvalue_tmp, 2), "s/mm2")

# Scale gradients amplitude accordingly.
if bvalue_tmp > np.max(bvalue):
    gscl_max = np.sqrt(np.max(bvalue) / bvalue_tmp)
else:
    gscl_max = 1

# Show final TE and b-values:
print("Final times:")
print("TE:", round(TE*1e3, 2), "ms")
for bv in range(1, nbvals+1):
    F = waveform * gscl_max * gscl[bv] * INV * short_f * seq.grad_raster_time
    bv = (2 * math.pi) ** 2 * np.matmul(np.matmul(np.transpose(F), C2), F) * short_f * seq.grad_raster_time * 1e-6  # [s/mm2]
    print(round(bv, 2), "s/mm2")

# Crusher gradients
gx_crush = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
gz_crush = make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

# TR delay - Takes everything into account
# EPI reading time:
# Distance between the center of the RF90s must be TR
EPI_time = calc_duration(gx) * Nyeff + calc_duration(gy) * (Nyeff - 1)
tr_per_slice = TR / n_slices
if fatsat_enable:
    tr_delay = math.floor((tr_per_slice - (TE - duration_center + EPI_time + pre_time) - rf_center_with_delay - calc_duration(gx_crush, gz_crush) \
                           - calc_duration(rf_fs, gz_fs)) \
                          / seq.grad_raster_time) * seq.grad_raster_time
else:
    tr_delay = math.floor((tr_per_slice - (TE - duration_center + EPI_time + pre_time) - rf_center_with_delay - calc_duration(gx_crush, gz_crush)) \
                          /seq.grad_raster_time)*seq.grad_raster_time

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
            seq.add_block(gx, adc)          # Read one line of k-space
            if i is not Nyeff-1:
                seq.add_block(gy)               # Phase blip
            gx.amplitude = -gx.amplitude    # Reverse polarity of read gradient

        seq.add_block(gx_crush, gz_crush)

        # Wait TR
        if tr_delay > 0:
            seq.add_block(make_delay(tr_delay))

""" DWI acquisition """
for bv in range(1, nbvals+1):
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
            gdiffx_d1 = make_trapezoid(channel='x', system=system, amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 0],
                                    duration=d1)
            gdiffy_d1 = make_trapezoid(channel='y', system=system, amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 1],
                                    duration=d1)
            gdiffz_d1 = make_trapezoid(channel='z', system=system, amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 2],
                                    duration=d1)

            seq.add_block(gdiffx_d1, gdiffy_d1, gdiffz_d1)

            # First RF180
            seq.add_block(gz_spoil_1)
            rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil_1)

            # Diffusion-weighting gradient d2
            gdiffx_d2 = make_trapezoid(channel='x', system=system, amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 0],
                                    duration=d2)
            gdiffy_d2 = make_trapezoid(channel='y', system=system, amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 1],
                                    duration=d2)
            gdiffz_d2 = make_trapezoid(channel='z', system=system, amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 2],
                                    duration=d2)

            seq.add_block(gdiffx_d2, gdiffy_d2, gdiffz_d2)

            # Diffusion-weighting gradient d3
            gdiffx_d3 = make_trapezoid(channel='x', system=system, amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 0],
                                       duration=d3)
            gdiffy_d3 = make_trapezoid(channel='y', system=system, amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 1],
                                       duration=d3)
            gdiffz_d3 = make_trapezoid(channel='z', system=system, amplitude=system.max_grad * gscl_max * gscl[bv] * gdir[d, 2],
                                       duration=d3)

            seq.add_block(gdiffx_d3, gdiffy_d3, gdiffz_d3)

            # Second RF180
            seq.add_block(gz_spoil_2)
            rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil_2)

            # Diffusion-weighting gradient d4
            gdiffx_d4 = make_trapezoid(channel='x', system=system, amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 0],
                                       duration=d4)
            gdiffy_d4 = make_trapezoid(channel='y', system=system, amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 1],
                                       duration=d4)
            gdiffz_d4 = make_trapezoid(channel='z', system=system, amplitude=-system.max_grad * gscl_max * gscl[bv] * gdir[d, 2],
                                       duration=d4)

            seq.add_block(gdiffx_d4, gdiffy_d4, gdiffz_d4)

            # Locate k-space
            seq.add_block(gx_pre, gy_pre)

            for i in range(Nyeff):
                seq.add_block(gx, adc)          # Read one line of k-space
                if i is not Nyeff-1:
                    seq.add_block(gy)               # Phase blip
                gx.amplitude = -gx.amplitude    # Reverse polarity of read gradient

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
    seqfname = "trse_dwi_" + str(n_slices) + "slices_" + str(bvalue) + "bvalues_" + str(ndirs) + "dirs_" + str(d4*1e3) + "d4_" + \
        str(round(TE * 1e3, 2)) + "TE_" + str(round(TR * 1e3)) + "TR_" + str(
        system.max_grad / system.gamma * 1e3) + "G_Max_" + \
               str(system.max_slew / system.gamma) + "SR_Max" + fatsat_str + pF_str
    os.mkdir("tests/" + seqfname)
    seq.write("tests/" + seqfname + "/" + seqfname + ".seq")
    print("Seq file saved --", seqfname)














