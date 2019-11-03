import math

import numpy as np

import diff_funcs as difunc

from pypulseq.Sequence.sequence import Sequence
from pypulseq.make_adc import make_adc
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts
from pypulseq.calc_duration import calc_duration
from pypulseq.make_delay import make_delay

import matplotlib.pyplot as plt

seq = Sequence()
seqfname = 'dwepi_3slices_2bval_3dirsnofs_New.seq'
fov = 240e-3
Nx = 96
Ny = 96
slice_thickness = 2.5e-3
n_slices = 3

# Partial Fourier
pF = 0.75
Nyeff = int (pF*Ny)
te=112e-3
tr=5

fatsat_enable=0
pe_enable = 1

nbvals=2
ndirs=3
#gradient scaling
gscl=np.sqrt(np.linspace(0., 1., nbvals+1))
gdir, nb0s = difunc.get_dirs(ndirs)

tr_per_slice = tr / n_slices

system = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s', rf_ringdown_time=30e-6,
              rf_dead_time=100e-6)

if fatsat_enable:
    b0 = 2.89
    sat_ppm = -3.45
    sat_freq = sat_ppm * 1e-6 * b0 * system.gamma
    rf_fs, _, _ = make_gauss_pulse(flip_angle=110 * math.pi / 180, system=system, duration=8e-3, bandwidth=abs(sat_freq),
                               freq_offset=sat_freq)
    gz_fs = make_trapezoid(channel='z', system=system, delay=calc_duration(rf_fs), area=1 / 1e-4)

rf, gz, _ = make_sinc_pulse(flip_angle=math.pi / 2, system=system, duration=3e-3, slice_thickness=slice_thickness,
                                  apodization=0.5, time_bw_product=4)

rf180, gz180, _ =  make_sinc_pulse(flip_angle=math.pi, system=system, duration=5e-3, slice_thickness=slice_thickness,
                            apodization=0.5, time_bw_product=4)
rf180.phase_offset = math.pi/2

gz_spoil = make_trapezoid(channel='z', system=system, area=6/slice_thickness, duration=3e-3)

delta_k = 1 / fov
k_width = Nx * delta_k
dwell_time = seq.grad_raster_time
readout_time = Nx * dwell_time
flat_time = math.ceil(readout_time / seq.grad_raster_time) * seq.grad_raster_time
gx = make_trapezoid(channel='x', system=system, amplitude=k_width / readout_time, flat_time=flat_time)
adc = make_adc(num_samples=Nx, duration=readout_time,
               delay=gx.rise_time + flat_time / 2 - (readout_time - dwell_time) / 2)

pre_time = 1e-3
gx_pre = make_trapezoid(channel='x', system=system, area=-gx.area / 2, duration=pre_time)
gz_reph = make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)
gy_pre = make_trapezoid(channel='y', system=system, area=- (Ny-Nyeff) * delta_k, duration=pre_time)

dur = math.ceil(2 * math.sqrt(delta_k / system.max_slew) / seq.grad_raster_time) * seq.grad_raster_time
gy = make_trapezoid(channel='y', system=system, area=delta_k, duration=dur)

duration_center = (calc_duration(gx))*(Ny-Nyeff-0.5)

delay_te1 = math.ceil((te/2 - calc_duration(gz)/2 - + calc_duration(gz_reph) - calc_duration(gz_spoil) - calc_duration(gz180)/2)/seq.grad_raster_time)*seq.grad_raster_time
delay_te2 = math.ceil((te/2 - calc_duration(gz180)/2 - calc_duration(gz_spoil) - calc_duration(gx_pre,gy_pre) - duration_center)/seq.grad_raster_time)*seq.grad_raster_time

min_te = math.ceil(2*max(calc_duration(gz)/2 + calc_duration(gz_reph) + calc_duration(gz_spoil) + calc_duration(gz180)/2, calc_duration(gz180)/2 + calc_duration(gz_spoil) + calc_duration(gx_pre,gy_pre) + duration_center)/ seq.grad_raster_time) * seq.grad_raster_time

assert np.all(te >= min_te)

#adding diffusion-weighting - based on TE for now
gdiff_dur = min(delay_te1,delay_te2)

rt= math.ceil(system.max_grad/system.max_slew/seq.grad_raster_time)*seq.grad_raster_time

gdiff = make_trapezoid(channel='x', system=system, amplitude=system.max_grad, duration=gdiff_dur)
gdiff_dur = calc_duration(gdiff)

delay_te1 = math.ceil((delay_te1 - gdiff_dur)/seq.grad_raster_time)*seq.grad_raster_time
delay_te2 = math.ceil((delay_te2 - gdiff_dur)/seq.grad_raster_time)*seq.grad_raster_time

print(gdiff_dur),print(gdiff_dur + + 2*calc_duration(gz_spoil) + calc_duration(gz180))
bv=difunc.calc_bval(system.max_grad,gdiff_dur,gdiff_dur + 2*calc_duration(gz_spoil) + calc_duration(gz180))
bv_smm_2=bv*1e-6
print(bv_smm_2)

gx_crush = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
gz_crush = make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

if fatsat_enable:
    min_tr = n_slices*round((calc_duration(gz_fs)  + calc_duration(gz) + calc_duration(gz_reph) + delay_te1 + 2*gdiff_dur + 2*calc_duration(gz_spoil) + calc_duration(gz180) + delay_te2 + calc_duration(gx_pre,gy_pre)
        + calc_duration(gx)*Nyeff + calc_duration(gx_crush,gz_crush))/seq.grad_raster_time)*seq.grad_raster_time

    assert np.all(tr >= min_tr)

    tr_delay = math.ceil((tr_per_slice - (calc_duration(gz_fs) + calc_duration(gz) + calc_duration(gz_reph) + delay_te1 + 2*calc_duration(gz_spoil) + calc_duration(gz180) + delay_te2 + 2*gdiff_dur + calc_duration(gx_pre,gy_pre)
        + calc_duration(gx)*Nyeff + calc_duration(gx_crush,gz_crush)))/seq.grad_raster_time) * seq.grad_raster_time
else:
    min_tr = n_slices * round(( calc_duration(gz) + calc_duration(gz_reph) + delay_te1 + 2 * gdiff_dur + 2 * calc_duration(gz_spoil) + calc_duration(gz180) + delay_te2 + calc_duration(gx_pre, gy_pre)
                               + calc_duration(gx) * Nyeff + calc_duration(gx_crush, gz_crush)) / seq.grad_raster_time) * seq.grad_raster_time

    assert np.all(tr >= min_tr)

    tr_delay = math.ceil((tr_per_slice - (calc_duration(gz) + calc_duration(gz_reph) + delay_te1 + 2 * calc_duration(gz_spoil) + calc_duration(gz180) + delay_te2 + 2 * gdiff_dur + calc_duration(gx_pre, gy_pre)
                + calc_duration(gx) * Nyeff + calc_duration(gx_crush, gz_crush))) / seq.grad_raster_time) * seq.grad_raster_time


#EPI calibration
for s in range(n_slices):
    print(s - (n_slices - 1) / 2)
    if fatsat_enable:
        seq.add_block(rf_fs, gz_fs)
    rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
    seq.add_block(rf, gz)
    seq.add_block(gz_reph)

    seq.add_block(make_delay(gdiff_dur + delay_te1))

    seq.add_block(gz_spoil)
    rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
    seq.add_block(rf180, gz180)
    seq.add_block(gz_spoil)

    seq.add_block(make_delay(delay_te2 + gdiff_dur))

    seq.add_block(gx_pre)

    for i in range(1, Nyeff + 1):
        seq.add_block(gx, adc)
        gx.amplitude = -gx.amplitude

    seq.add_block(gx_crush, gz_crush)

    if tr_delay > 0:
        seq.add_block(make_delay(tr_delay))


#b-zero acquisition
for d in range(nb0s):
    for s in range(n_slices):
        print(s - (n_slices - 1) / 2)
        if fatsat_enable:
            seq.add_block(rf_fs, gz_fs)
        rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf, gz)
        seq.add_block(gz_reph)

        seq.add_block(make_delay(gdiff_dur + delay_te1))

        seq.add_block(gz_spoil)
        rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf180, gz180)
        seq.add_block(gz_spoil)

        seq.add_block(make_delay(delay_te2 + gdiff_dur))

        seq.add_block(gx_pre)

        for i in range(1, Nyeff + 1):
            seq.add_block(gx, adc)
            seq.add_block(gy)
            gx.amplitude = -gx.amplitude

        seq.add_block(gx_crush, gz_crush)

        if tr_delay > 0:
            seq.add_block(make_delay(tr_delay))

#dwi acquisitions
for bv in range(1,nbvals+1):
    for d in range(ndirs):
        for s in range(n_slices):
            if fatsat_enable:
                seq.add_block(rf_fs, gz_fs)
            rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf, gz)
            seq.add_block(gz_reph)

            gdiffx = make_trapezoid(channel='x', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 0],duration=gdiff_dur)
            gdiffy = make_trapezoid(channel='y', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 1],duration=gdiff_dur)
            gdiffz = make_trapezoid(channel='z', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 2],duration=gdiff_dur)

            seq.add_block(gdiffx, gdiffy, gdiffz)
            if delay_te1 > 0:
                seq.add_block(make_delay(delay_te1))

            seq.add_block(gz_spoil)
            rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil)

            seq.add_block(gdiffx, gdiffy, gdiffz)

            if delay_te2 > 0:
                seq.add_block(make_delay(delay_te2))

            seq.add_block(gx_pre, gy_pre)

            for i in range(1, Nyeff + 1):
                seq.add_block(gx, adc)
                seq.add_block(gy)
                gx.amplitude = -gx.amplitude

            seq.add_block(gx_crush, gz_crush)

            if tr_delay > 0:
                seq.add_block(make_delay(tr_delay))

seq.plot()

#ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

#time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
#plt.plot(time_axis, ktraj.T)
#plt.plot(t_adc, ktraj_adc[0, :], '.')
#plt.figure()
#plt.plot(ktraj[0, :], ktraj[1, :], 'b')
#plt.axis('equal')
#plt.plot(ktraj_adc[0, :], ktraj_adc[1, :], 'r.')
#plt.show()

seq.write(seqfname)
