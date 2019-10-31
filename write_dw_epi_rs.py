import math
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np

import diff_funcs as difunc

from pypulseq.Sequence.sequence import Sequence
from pypulseq.add_gradients import add_gradients
from pypulseq.align import align
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.make_delay import make_delay
from pypulseq.opts import Opts
from pypulseq.split_gradient_at import split_gradient_at

seq = Sequence()
seqfname = 'nondwepi_1slice.seq'
fov = 240e-3
Nx = 96
Ny = 96
slice_thickness = 2.5e-3
n_slices = 1

# Partial Fourier
pF = 0.75
Nyeff = int (pF*Ny)
te=84e-3
tr=1

pe_enable = 1

nbvals=0
ndirs=1
#gradient scaling
gscl=np.sqrt(np.linspace(0., 1., nbvals+1))
gdir, nb0s = difunc.get_dirs(ndirs)

#to avoid exceeding gradient limits
gdfact = 1.0

tr_per_slice = tr / n_slices

system = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s', rf_ringdown_time=30e-6,
              rf_dead_time=100e-6, adc_dead_time=20e-6)

b0 = 2.89
sat_ppm = -3.45
sat_freq = sat_ppm * 1e-6 * b0 * system.gamma
rf_fs, _, _ = make_gauss_pulse(flip_angle=110 * math.pi / 180, system=system, duration=8e-3, bandwidth=abs(sat_freq),
                               freq_offset=sat_freq)
gz_fs = make_trapezoid(channel='z', system=system, delay=calc_duration(rf_fs), area=1 / 1e-4)

rf, gz, gz_reph = make_sinc_pulse(flip_angle=math.pi / 2, system=system, duration=3e-3, slice_thickness=slice_thickness,
                                  apodization=0.5, time_bw_product=4)

rf180, gz180, _ =  make_sinc_pulse(flip_angle=math.pi, system=system, duration=5e-3, slice_thickness=slice_thickness,
                            apodization=0.5, time_bw_product=4)
rf180.phase_offset = math.pi/2

gz_spoil = make_trapezoid(channel='z', system=system, area=2*gz.area, duration=3e-3)

delta_k = 1 / fov
k_width = Nx * delta_k
readout_time = 1.05e-3

blip_dur = math.ceil(2 * math.sqrt(delta_k / system.max_slew) / seq.grad_raster_time / 2) * seq.grad_raster_time * 2
gy = make_trapezoid(channel='y', system=system, area=delta_k, duration=blip_dur)

extra_area = blip_dur / 2 * blip_dur / 2 * system.max_slew
gx = make_trapezoid(channel='x', system=system, area=k_width + extra_area, duration=readout_time + blip_dur)
actual_area = gx.area - gx.amplitude / gx.rise_time * blip_dur / 2 * blip_dur / 2 / 2 - gx.amplitude / gx.fall_time * blip_dur / 2 * blip_dur / 2 / 2
gx.amplitude = gx.amplitude / actual_area * k_width
gx.area = gx.amplitude * (gx.flat_time + gx.rise_time / 2 + gx.fall_time / 2)
gx.flat_area = gx.amplitude * gx.flat_time

adc_dwell_nyquist = delta_k / gx.amplitude
adc_dwell = math.floor(adc_dwell_nyquist * 1e7) * 1e-7
#both adc_dwell and grad_time must fall on the same grid
adc_dwell = math.floor(adc_dwell / seq.grad_raster_time)*seq.grad_raster_time
#making sure the number of samples is even
adc_samples = math.floor(readout_time / adc_dwell / 4) * 4

if adc_samples == 0:
    print('error! The readout time needs to be longer!')

adc = make_adc(num_samples=adc_samples, dwell=adc_dwell, delay=blip_dur / 2)

time_to_center = system.adc_dead_time * (adc_samples - 1) / 2
adc.delay = round((gx.rise_time + gx.flat_time / 2 - time_to_center) / seq.grad_raster_time) * seq.grad_raster_time

gy_parts = split_gradient_at(grad=gy, time_point=blip_dur / 2, system=system)
gy_blipup, gy_blipdown, _ = align('right', gy_parts[0], 'left', gy_parts[1], gx)
gy_blipdownup = add_gradients([gy_blipdown, gy_blipup], system=system)

gy_blipup.waveform = gy_blipup.waveform * pe_enable
gy_blipdown.waveform = gy_blipdown.waveform * pe_enable
gy_blipdownup.waveform = gy_blipdownup.waveform * pe_enable

gx_pre = make_trapezoid(channel='x', system=system, area=-gx.area / 2)
gy_pre = make_trapezoid(channel='y', system=system, area= -(1-pF)*Ny * delta_k)

gy_pre = make_trapezoid(channel='y', system=system, area=gy_pre.area, duration=calc_duration(gx_pre, gy_pre))
gy_pre.amplitude = gy_pre.amplitude * pe_enable

dur2center = (calc_duration(gx))*Ny*(1-pF)

delay_te1 = math.ceil((te/2 - calc_duration(gz)/2 - + calc_duration(gz_reph) - calc_duration(gz_spoil) - calc_duration(rf180)/2)/seq.grad_raster_time)*seq.grad_raster_time
delay_te2 = math.ceil((te/2 - calc_duration(rf180)/2 - calc_duration(gz_spoil) - calc_duration(gx_pre,gy_pre) - dur2center)/seq.grad_raster_time)*seq.grad_raster_time

min_te = math.ceil(2*max(calc_duration(gz)/2 + calc_duration(gz_reph) + calc_duration(gz_spoil) + calc_duration(rf180)/2, calc_duration(rf180)/2 + calc_duration(gz_spoil) + calc_duration(gx_pre,gy_pre) + dur2center)/ seq.grad_raster_time) * seq.grad_raster_time

assert np.all(te >= min_te)

#adding diffusion-weighting - based on TE for now
gdiff_dur = min(delay_te1,delay_te2)

gdiff = make_trapezoid(channel='x', system=system, amplitude=gdfact*system.max_grad, duration=gdiff_dur)

delay_te1 = math.ceil((delay_te1 - gdiff_dur)/seq.grad_raster_time)*seq.grad_raster_time
delay_te2 = math.ceil((delay_te2 - gdiff_dur)/seq.grad_raster_time)*seq.grad_raster_time

gdiff_dur = calc_duration(gdiff)

print(gdiff_dur),print(gdiff_dur +  2*calc_duration(gz_spoil) + calc_duration(rf180))
bv=difunc.calc_bval(gdfact*system.max_grad,gdiff_dur,gdiff_dur + 2*calc_duration(gz_spoil) + calc_duration(rf180))
bv_smm_2=bv*1e-6
print(bv_smm_2)

gx_crush = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
gz_crush = make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

min_tr = n_slices*round((calc_duration(rf_fs) + calc_duration(gz_fs) + calc_duration(rf) + calc_duration(gz_reph) + 2*calc_duration(gz_spoil) + calc_duration(rf180) + calc_duration(gx_pre,gy_pre)
        + calc_duration(gx)*Ny*pF + calc_duration(gx_crush,gz_crush))/seq.grad_raster_time)*seq.grad_raster_time

assert np.all(tr >= min_tr)

tr_delay = round((tr_per_slice - (calc_duration(rf_fs) + calc_duration(gz_fs) + calc_duration(rf) + calc_duration(gz_reph) + 2*calc_duration(gz_spoil) + calc_duration(rf180) + calc_duration(gx_pre,gy_pre)
        + calc_duration(gx)*Ny*pF + calc_duration(gx_crush,gz_crush)))/seq.grad_raster_time) * seq.grad_raster_time


#EPI calibration
for s in range(n_slices):
    print(s - (n_slices - 1) / 2)
    seq.add_block(rf_fs, gz_fs)
    rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
    seq.add_block(rf, gz)
    seq.add_block(gz_reph)

    seq.add_block(make_delay(gdiff_dur + delay_te1))

    seq.add_block(gz_spoil)
    seq.add_block(rf180, gz180)
    seq.add_block(gz_spoil)

    seq.add_block(make_delay(gdiff_dur))

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
        seq.add_block(rf_fs, gz_fs)
        rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
        seq.add_block(rf, gz)
        seq.add_block(gz_reph)

        seq.add_block(make_delay(gdiff_dur + delay_te1))

        seq.add_block(gz_spoil)
        seq.add_block(rf180, gz180)
        seq.add_block(gz_spoil)

        seq.add_block(make_delay(gdiff_dur))

        seq.add_block(make_delay(delay_te2 + gdiff_dur))

        seq.add_block(gx_pre,gy_pre)

        for i in range(1, Nyeff + 1):
            if i == 1:
                seq.add_block(gx, gy_blipup, adc)
            elif i == Nyeff:
                seq.add_block(gx, gy_blipdown, adc)
            else:
                seq.add_block(gx, gy_blipdownup, adc)
            gx.amplitude = -gx.amplitude

        seq.add_block(gx_crush, gz_crush)

        if tr_delay > 0:
            seq.add_block(make_delay(tr_delay))

#dwi acquisitions
for bv in range(1,nbvals+1):
    for d in range(ndirs):
        for s in range(n_slices):
            seq.add_block(rf_fs, gz_fs)
            rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf, gz)
            seq.add_block(gz_reph)

            if delay_te1 > 0:
                seq.add_block(make_delay(delay_te1))

            rt = math.ceil(system.max_grad * gscl[bv] * gdir[d, 0] / system.max_slew / seq.grad_raster_time) * seq.grad_raster_time
            ft = gdiff_dur - rt
            gdiffx = make_trapezoid(channel='x', system=system, amplitude=gdfact*system.max_grad * gscl[bv] * gdir[d, 0],duration=gdiff_dur)
                                   # area=(gdiff_dur - rt) * system.max_grad * gscl[bv] * gdir[d, 0], flat_time=ft)
            gdiffy = make_trapezoid(channel='y', system=system, amplitude=gdfact*system.max_grad * gscl[bv] * gdir[d, 1],duration=gdiff_dur)
                                   # area=(gdiff_dur - rt) * system.max_grad * gscl[bv] * gdir[d, 1], flat_time=ft)
            gdiffz = make_trapezoid(channel='z', system=system, amplitude=gdfact*system.max_grad * gscl[bv] * gdir[d, 2],duration=gdiff_dur)
                                    #area=(gdiff_dur - rt) * system.max_grad * gscl[bv] * gdir[d, 2], flat_time=ft)

            seq.add_block(gdiffx, gdiffy, gdiffz)

            seq.add_block(gz_spoil)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil)

            seq.add_block(gdiffx, gdiffy, gdiffz)

            if delay_te2 > 0:
                seq.add_block(make_delay(delay_te2))

            seq.add_block(gx_pre, gy_pre)

            for i in range(1, Nyeff + 1):
                if i == 1:
                    seq.add_block(gx, gy_blipup, adc)
                elif i == Nyeff:
                    seq.add_block(gx, gy_blipdown, adc)
                else:
                    seq.add_block(gx, gy_blipdownup, adc)
                gx.amplitude = -gx.amplitude

            seq.add_block(gx_crush, gz_crush)

            if tr_delay > 0:
                seq.add_block(make_delay(tr_delay))



ok, error_report = seq.check_timing()

if ok:
    print('Timing check passed succesfully.')
else:
    warn('Timing check failed! Error listing follows:\n')
    print(error_report)
    print('\n')

seq.plot()

#ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

#time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
#plt.plot(time_axis, ktraj.T)
#plt.plot(t_adc, ktraj_adc[0], '.')
#plt.figure()
#plt.plot(ktraj_adc[0], ktraj_adc[1], 'b')
#plt.axis('equal')
#plt.plot(ktraj_adc[0], ktraj_adc[1], 'r.')
#plt.show()

seq.set_definition('FOV', np.array([fov, fov, slice_thickness]) * 1e3)
seq.set_definition('Name', 'dw_epi')

seq.write(seqfname)
