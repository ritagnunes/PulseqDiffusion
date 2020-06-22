"""
Oscar Pena-Nogales (opennog(at)lpi(dot)tel(dot)uva(dot)es)
Rita G. Nunes
January 2020
Laboratorio de Procesado de Imagen, Universidad de Valladolid, Valladolid, Spain
LASEEB, Instituto Superior Técnico, Lisbon, Portugal

This function builds a spin echo EPI with ramp sampling diffusion weighting sequence optimizing the TE for a desired b-value.
"""

import math
import sys
import numpy as np
import scipy.io
import os
import matplotlib.pyplot as plt

import diff_funcs as difunc

from pypulseq.make_adc import make_adc
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.make_arbitrary_grad import make_arbitrary_grad
from pypulseq.opts import Opts
from pypulseq.calc_duration import calc_duration
from pypulseq.calc_rf_center import calc_rf_center
from pypulseq.make_delay import make_delay
from pypulseq.Sequence.sequence import Sequence
from pypulseq.split_gradient_at import split_gradient_at
from pypulseq.align import align
from pypulseq.add_gradients import add_gradients

# Create new Sequence Object
seq = Sequence()
# Due to the floating point uncertainty (https://docs.python.org/3/tutorial/floatingpoint.html),
# sometimes I'm having trouble when '* seq.grad_raster_time',  it returns number like:
# 0.005560000000000001 which give me problems afterwards. However, such problem is solved if I used a manually inputed
# raster time (1/100000). To see this problem do: 556*1e-5 vs 556/1e5. Thus, whenever we have to multiply we divide
# by the equivalent (100000)

# I'm still getting some problems, so the best option to avoid them is to work with integers and only use the actual
# timings by the end. Thus making all intermediate operations (+ and -) with n.

i_raster_time = 100000
# I need to check the manually inputed inverse raster time and the actual value are the same.
assert 1/i_raster_time == seq.grad_raster_time, "Manualy inputed inverse raster time does not match the actual value."

# Acquisition Parameters
# FOV and resolution
fov = 240e-3  # [m]
Nx = 64
Ny = 64
slice_thickness = 2.5e-3  # [m]
n_slices = 5

# Partial Fourier
pF = 1
Nyeff = int(pF*Ny) # Number of Ny samples acquired
if pF is not 1:
    pF_str = "_" + str(pF) + "pF"
else:
    pF_str = ""

# b-value
bvalue = np.array([500])  # [s/mm2]

# Spin-Echo parameters
TR = 5000e-3  # [s]
n_TR = math.ceil(TR * i_raster_time)

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
kplot = 1

# b-value parameters
nbvals = np.shape(bvalue)[0]
ndirs = 3

# Gradient Scaling
gscl = np.zeros(nbvals+1)
gscl[1:] = np.sqrt(bvalue/np.max(bvalue))
gdir, nb0s = difunc.get_dirs(ndirs)

# Set system limits
system = Opts(max_grad=32, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s', rf_ringdown_time=20e-6,
              rf_dead_time=100e-6, adc_dead_time=10e-6, rf_raster_time=1e-6)

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

# Refocusing pulse with spoiling gradients
rf180, gz180, _ = make_sinc_pulse(flip_angle=math.pi, system=system, duration=5e-3, slice_thickness=slice_thickness,
                            apodization=0.5, time_bw_product=4)
rf180.phase_offset = math.pi/2
gz_spoil = make_trapezoid(channel='z', system=system, area=6/slice_thickness, duration=3e-3)

"""" Define other gradients and ADC events """
# Gx Without ramp sampling
delta_k = 1 / fov
k_width = Nx * delta_k
dwell_time = seq.grad_raster_time  # Full receiver bandwith
# For initialization purposes - Minimum readout time possible for standard EPI.
readout_time = Nx * dwell_time  # T_acq (acquisition time)

# Phase blip in shortest possible time
gy = make_trapezoid(channel='y', system=system, area=delta_k)
dur = math.ceil(calc_duration(gy) / seq.grad_raster_time) / i_raster_time

extra_area = dur / 2 * dur / 2 * system.max_slew

""" Obtain the minimum readout time with a fixed maximum gradient intensity and min dwell time"""
condition = False
while not condition:
    try:
        gx = make_trapezoid(channel='x', system=system, area=k_width + extra_area, duration=readout_time + dur)
        # We control the ruduction factor manually. I don't know if the readout gradients ever in practice achieve the
        # max amplitude of the gradients. Also, an increase on the amplitude of the gx produces a decrease of the
        # adc_dwell time. However, I don't know if there is a limit for that value either.
        if gx.amplitude < 0.5 * system.max_grad:
            condition = True
        else:
            readout_time = math.ceil(readout_time / seq.grad_raster_time + 1) / i_raster_time
    except:
        readout_time = math.ceil(readout_time / seq.grad_raster_time + 1) / i_raster_time

# Actual Area
# Para los desplazamientos por el k-spacio mientras se lee solo hay que considerar el area normal (sin extra),
# por ello, al haber calculado el gradiente con el area extra, esta hay que reducirla para muestrear en los sitios
# deseados.

# Lo unico, que ahora el area es ligeramente superior porque en el primer paso calculo el area extra con el
# system.max_slew y despues lo resto con el slew rate real (gx.amplitude / gx.rise_time), lo que deja un pelin mas de
# area
actual_area = gx.area - gx.amplitude / gx.rise_time * dur / 2 * dur / 2 / 2 - gx.amplitude / gx.fall_time * dur / 2 * dur / 2 / 2
gx.amplitude = gx.amplitude / actual_area * k_width
gx.area = gx.amplitude * (gx.flat_time + gx.rise_time / 2 + gx.fall_time / 2)
gx.flat_area = gx.amplitude * gx.flat_time

# Update number of samples and ADC sampling rate (increase receiver bandwidth) to maintain the FOV.
adc_samples = math.ceil(readout_time / dwell_time / 4) * 4  # Apparently, on Siemens the number of samples needs to be divisible by 4...
adc = make_adc(num_samples=adc_samples, system=system, dwell=dwell_time, delay=dur / 2)

# Center the sampling where needed - This is needed due to the new sampling rate
# To align odd and even readouts. However, on the real hardware this misalignment is
# much stronger due to the gradient delays -> thus we can not achieve perfect matching.
time_to_center = dwell_time * (adc_samples - 1) / 2
adc.delay = np.floor((gx.rise_time + gx.flat_time / 2 - time_to_center) * 1e6) / 1e6  # For some reason is has to be multiple of the rf_raster_time....

# Split Gy to put half on each consecutive Gx and produce a combined synthetic gradient
gy_parts = split_gradient_at(grad=gy, time_point=dur / 2, system=system)
gy_blipup, gy_blipdown, _ = align('right', gy_parts[0], 'left', gy_parts[1], gx)
# Incluir la siguiente linea se queremos que las muestras esten centradas y no queremos que las ultimas se monten con
# la activacion dl gy. - Aunque una mejor opcion es cambiar la forma del blip partido para que ocupe un instante menos
# de tiempo manteniendo el mismo area.
gy_blipup.delay = math.ceil(gy_blipup.delay / system.grad_raster_time - 1) / i_raster_time
gy_blipdownup = add_gradients([gy_blipdown, gy_blipup], system=system)
assert math.ceil(calc_duration(gy_blipdownup) / system.grad_raster_time) == math.ceil(calc_duration(gy_blipup) / system.grad_raster_time), \
    "Timing mismatch between gy_blipup and gy_blipdownup."
assert math.ceil(calc_duration(gy_blipdownup) / system.grad_raster_time) == math.ceil(calc_duration(gx) / system.grad_raster_time), \
    "Timing mismatch between gy_blipup and gy_blipdownup."

# Pre-phasing gradients
pre_time = 1e-3
n_pre_time = math.ceil(pre_time * i_raster_time)
gx_pre = make_trapezoid(channel='x', system=system, area=-gx.area / 2, duration=pre_time)
gz_reph = make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)
gy_pre = make_trapezoid(channel='y', system=system, area=-(Ny / 2 - 0.5 - (Ny - Nyeff)) * delta_k, duration=pre_time) # Es -0.5 y no +0.5 porque hay que pensar en areas, no en rayas!

"""" Obtain TE and diffusion-weighting gradient waveform """
# For S&T monopolar waveforms
# From an initial TE, check we satisfy all constraints -> otherwise increase TE.
# Once all constraints are okay -> check b-value, if it is lower than the target one -> increase TE
# Looks time-inefficient but it is fast enough to make it user-friendly.

# Calculate some times constant throughout the process
# We need to compute the exact time sequence. For the normal SE-MONO-EPI sequence micro second differences
# are not important, however, if we wanna import external gradients the allocated time for them needs to
# be the same, and thus exact timing is mandatory. With this in mind, we establish the following rounding rules:
# Duration of RFs + spoiling, and EPI time to the center of the k-space is always math.ceil().

# El time(gy) se refiere a los numeros de blips por lo que se resta 0.5 ya que se asume que el numero de lineas Ny es par.
# El time(gx) se refiere a lo que tarda en leer cada linea del k-space. Por lo que si tenemos Ny par, se tarda la mitad de lineas mas media (0.5)
n_duration_center = math.ceil((calc_duration(gx)*(Ny/2 + 0.5 - (Ny - Nyeff)) + calc_duration(gx_pre, gy_pre)) / seq.grad_raster_time)
rf_center_with_delay = rf.delay + calc_rf_center(rf)[0]

n_rf90r = math.ceil((calc_duration(gz) - rf_center_with_delay + pre_time) / seq.grad_raster_time)
n_rf180r = math.ceil((calc_duration(rf180) + 2 * calc_duration(gz_spoil)) / 2 / seq.grad_raster_time)
n_rf180l = math.floor((calc_duration(rf180) + 2 * calc_duration(gz_spoil)) / 2 / seq.grad_raster_time)

# Find minimum TE considering the readout times.
n_TE = math.ceil(40e-3 / seq.grad_raster_time)
n_delay_te2 = -1
while n_delay_te2 <= 0:
    n_TE = n_TE + 2

    n_tINV = math.floor(n_TE / 2)
    n_delay_te2 = n_tINV - n_rf180r - n_duration_center

# Find minimum TE for the target b-value
bvalue_tmp = 0
while bvalue_tmp < np.max(bvalue):
    n_TE = n_TE + 2

    n_tINV = math.floor(n_TE / 2)
    n_delay_te1 = n_tINV - n_rf90r - n_rf180l
    delay_te1 = n_delay_te1 / i_raster_time
    n_delay_te2 = n_tINV - n_rf180r - n_duration_center
    delay_te2 = n_delay_te2 / i_raster_time

    # Waveform Ramp time
    n_gdiff_rt = math.ceil(system.max_grad / system.max_slew / seq.grad_raster_time)

    # Select the shortest available time
    n_gdiff_delta = min(n_delay_te1, n_delay_te2)
    n_gdiff_Delta = n_delay_te1 + 2 * math.ceil(calc_duration(gz_spoil) / seq.grad_raster_time) + math.ceil(calc_duration(gz180) / seq.grad_raster_time)

    gdiff = make_trapezoid(channel='x', system=system, amplitude=system.max_grad, duration=n_gdiff_delta / i_raster_time)

    # delta only corresponds to the rectangle.
    n_gdiff_delta = n_gdiff_delta - 2 * n_gdiff_rt

    bv = difunc.calc_bval(system.max_grad, n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time, n_gdiff_rt / i_raster_time)
    bvalue_tmp = bv * 1e-6

# Show final TE and b-values:
print("TE:", round(n_TE / i_raster_time * 1e3, 2), "ms")
for bv in range(1, nbvals+1):
    print(round(difunc.calc_bval(system.max_grad * gscl[bv], n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time, n_gdiff_rt / i_raster_time) * 1e-6, 2), "s/mm2")

# Crusher gradients
gx_crush = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
gz_crush = make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

# TR delay - Takes everything into account
# EPI reading time:
# Distance between the center of the RF90s must be TR
n_EPI_time = int(math.ceil(calc_duration(gx) * i_raster_time) * Nyeff)
n_tr_per_slice = math.ceil(TR / n_slices * i_raster_time)
if fatsat_enable:
    n_tr_delay = n_tr_per_slice - (n_TE - n_duration_center + n_EPI_time + n_pre_time) \
                        - math.ceil(rf_center_with_delay * i_raster_time) \
                        - math.ceil(calc_duration(gx_crush, gz_crush) * i_raster_time) \
                        - math.ceil(calc_duration(rf_fs, gz_fs) * i_raster_time)
else:
    n_tr_delay = n_tr_per_slice - (n_TE - n_duration_center + n_EPI_time + n_pre_time) \
                    - math.ceil(rf_center_with_delay * i_raster_time) \
                    - math.ceil(calc_duration(gx_crush, gz_crush) * i_raster_time)
tr_delay = n_tr_delay / i_raster_time

# Check TR delay time
assert n_tr_delay > 0, "Such parameter configuration needs longer TR."

# Delay time
# Time between the gradient and the RF180. This time might be zero some times, although it is not normal.
n_gap_te1 = n_delay_te1 - n_delay_te2
gap_te1 = n_gap_te1 / i_raster_time

""" EPI calibration """
for s in range(n_slices):
    # Fat saturation
    if fatsat_enable:
        seq.add_block(rf_fs, gz_fs)

    # RF90
    rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
    seq.add_block(rf, gz)
    seq.add_block(gz_reph)

    # Delay for RF180
    seq.add_block(make_delay(delay_te1))

    # RF180
    seq.add_block(gz_spoil)
    rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
    seq.add_block(rf180, gz180)
    seq.add_block(gz_spoil)

    # Delay for EPI
    seq.add_block(make_delay(delay_te2))

    # Locate k-space - Only on the frequency encoding direction.
    seq.add_block(gx_pre)

    for i in range(Nyeff):
        seq.add_block(gx, adc)
        gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient

    seq.add_block(gx_crush, gz_crush)

    # Wait TR
    if tr_delay > 0:
        seq.add_block(make_delay(tr_delay))

# """ b-zero acquisition """
for d in range(nb0s):
    for s in range(n_slices):
        # Fat saturation
        if fatsat_enable:
            seq.add_block(rf_fs, gz_fs)

        # RF90
        rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf, gz)
        seq.add_block(gz_reph)

        # Delay for RF180
        seq.add_block(make_delay(delay_te1))

        # RF180
        seq.add_block(gz_spoil)
        rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf180, gz180)
        seq.add_block(gz_spoil)

        # Delay for EPI
        seq.add_block(make_delay(delay_te2))

        # Locate k-space
        seq.add_block(gx_pre, gy_pre)

        for i in range(Nyeff):
            if i == 0:
                seq.add_block(gx, gy_blipup, adc)  # Read first line of k-space
            elif i == Nyeff - 1:
                seq.add_block(gx, gy_blipdown, adc)  # Read last line of k-space
            else:
                seq.add_block(gx, gy_blipdownup, adc)  # Read one line of k-space
            gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient

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

            # Diffusion-weighting gradient
            gdiffx = make_trapezoid(channel='x', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 0],
                                    duration=calc_duration(gdiff))
            gdiffy = make_trapezoid(channel='y', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 1],
                                    duration=calc_duration(gdiff))
            gdiffz = make_trapezoid(channel='z', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 2],
                                    duration=calc_duration(gdiff))

            seq.add_block(gdiffx, gdiffy, gdiffz)

            # Delay for RF180
            seq.add_block(make_delay(gap_te1))

            # RF180
            seq.add_block(gz_spoil)
            rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil)

            # Diffusion-weighting gradient
            seq.add_block(gdiffx, gdiffy, gdiffz)

            # Locate k-space
            seq.add_block(gx_pre, gy_pre)

            for i in range(Nyeff):
                if i == 0:
                    seq.add_block(gx, gy_blipup, adc)  # Read first line of k-space
                elif i == Nyeff - 1:
                    seq.add_block(gx, gy_blipdown, adc)  # Read last line of k-space
                else:
                    seq.add_block(gx, gy_blipdownup, adc)  # Read one line of k-space
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
    plt.grid(True)
    plt.show()

# Save the sequence
if save_flag:
    print("Test Report: " + seq.test_report())
    TE = n_TE / i_raster_time
    TR = n_TR / i_raster_time
    seqfname = "se_dwi_rs_" + str(n_slices) + "slices_" + str(bvalue) + "bvalues_" + str(ndirs) + "dirs_" + str(
        round(TE * 1e3, 2)) + "TE_" + str(round(TR * 1e3)) + "TR_" + str(system.max_grad/system.gamma*1e3) + "G_Max_" +\
               str(system.max_slew/system.gamma) + "SR_Max" + fatsat_str + pF_str
    os.mkdir("tests/" + seqfname)
    seq.write("tests/" + seqfname + "/" + seqfname + ".seq")
    print("Seq file saved --", seqfname)











