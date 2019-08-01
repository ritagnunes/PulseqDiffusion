import math
import numpy as np

import diff_funcs as difunc

from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.make_delay import make_delay
from pypulseq.opts import Opts
import matplotlib.pyplot as plt

seq = Sequence()
fov = 220e-3
Nx = 64
Ny = 64
slice_thickness = 3e-3
n_slices = 1
te=100e-3
tr=200e-3

nbvals=1
ndirs=12
gscl=np.linspace(0., 1., nbvals+1)

system = Opts(max_grad=32, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s', rf_ringdown_time=30e-6,
              rf_dead_time=100e-6)

gdir = difunc.get_dirs(ndirs)

rf, gz, _ = make_sinc_pulse(flip_angle=math.pi / 2, system=system, duration=3e-3, slice_thickness=slice_thickness,
                            apodization=0.5, time_bw_product=4)

rf180, gz180, _ =  make_sinc_pulse(flip_angle=math.pi, system=system, duration=3e-3, slice_thickness=slice_thickness,
                            apodization=0.5, time_bw_product=4)
rf180.phase_offset = math.pi/2

gz_spoil = make_trapezoid(channel='z', system=system, area=2*gz.area, duration=3e-3)

delta_k = 1 / fov
k_width = Nx * delta_k
dwell_time = 4e-6
readout_time = Nx * dwell_time
flat_time = math.ceil(readout_time * 1e5) * 1e-5
gx = make_trapezoid(channel='x', system=system, amplitude=k_width / readout_time, flat_time=flat_time)
adc = make_adc(num_samples=Nx, duration=readout_time,
               delay=gx.rise_time + flat_time / 2 - (readout_time - dwell_time) / 2)

pre_time = 8e-4
gx_pre = make_trapezoid(channel='x', system=system, area=-gx.area / 2, duration=pre_time)
gz_reph = make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)
gy_pre = make_trapezoid(channel='y', system=system, area=-Ny / 2 * delta_k, duration=pre_time)

dur = math.ceil(2 * math.sqrt(delta_k / system.max_slew) / 10e-6) * 10e-6
gy = make_trapezoid(channel='y', system=system, area=delta_k, duration=dur)

duration_center = (calc_duration(gx)+calc_duration(gy))*Ny/2+0.5*calc_duration(gx)

delay_te1 = math.ceil((te/2 - calc_duration(gz)/2 - pre_time - calc_duration(gz_spoil) - calc_duration(rf180)/2)/seq.grad_raster_time)*seq.grad_raster_time
delay_te2 = math.ceil((te/2 - calc_duration(rf180)/2 - calc_duration(gz_spoil) - pre_time - duration_center)/seq.grad_raster_time)*seq.grad_raster_time

min_te = math.ceil(2*max(calc_duration(gz)/2 + pre_time + calc_duration(gz_spoil) + calc_duration(rf180)/2, calc_duration(rf180)/2 + calc_duration(gz_spoil) + pre_time + duration_center)/ seq.grad_raster_time) * seq.grad_raster_time

assert np.all(te >= min_te)

#adding diffusion-weighting - based on TE for now
gdiff_dur = min(delay_te1,delay_te2)

gdiff = make_trapezoid(channel='x', system=system, amplitude=system.max_grad, duration=gdiff_dur)

delay_te1 = math.ceil((delay_te1 - gdiff_dur)/seq.grad_raster_time)*seq.grad_raster_time
delay_te2 = math.ceil((delay_te2 - gdiff_dur)/seq.grad_raster_time)*seq.grad_raster_time

delta = round(gdiff.rise_time + gdiff.flatTime,5)
Delta = delta + delay_te1 + 2*calc_duration(gz_spoil) + calc_duration(rf180) + delay_te2

bmax = difunc.calc_bval(system.max_grad, delta, Delta, gdiff.rise_time)
bmax_smm2 = bmax*1e-6
print(bmax_smm2)

gx_crush = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
gz_crush = make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

delay_tr = math.ceil((tr/n_slices - (calc_duration(gz) + pre_time + delay_te1 + 2*calc_duration(gz_spoil) + calc_duration(rf180) + 
        delay_te2 + pre_time + (calc_duration(gx)+calc_duration(gy))*Ny))/seq.grad_raster_time)*seq.grad_raster_time

assert np.all(delay_tr >= calc_duration(gx_crush, gz_crush))


for bv in range(1,(nbvals+1)):
    for d in range(ndirs):    
        for s in range(n_slices):
            rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf, gz)
            seq.add_block(gz_reph)
            
            gdiffx=make_trapezoid(channel='x',system=system, amplitude=system.max_grad*gscl(bv)*gdir[d,1], duration=gdiff_dur);
            gdiffy=make_trapezoid(channel='y',system=system, amplitude=system.max_grad*gscl(bv)*gdir[d,2], duration=gdiff_dur);
            gdiffz=make_trapezoid(channel='z',system=system, amplitude=system.max_grad*gscl(bv)*gdir[d,3], duration=gdiff_dur);
            seq.addBlock(gdiffx,gdiffy,gdiffz);
            
            if delay_te1>0: 
                seq.add_block(make_delay(delay_te1))
            seq.add_block(gz_spoil)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil)
            if delay_te2>0:
                seq.add_block(make_delay(delay_te2))
            
            gdiffx=make_trapezoid(channel='x',system=system, amplitude=system.max_grad*gscl(bv)*gdir[d,1], duration=gdiff_dur);
            gdiffy=make_trapezoid(channel='y',system=system, amplitude=system.max_grad*gscl(bv)*gdir[d,2], duration=gdiff_dur);
            gdiffz=make_trapezoid(channel='z',system=system, amplitude=system.max_grad*gscl(bv)*gdir[d,3], duration=gdiff_dur);
            seq.addBlock(gdiffx,gdiffy,gdiffz);
            
            seq.add_block(gx_pre, gy_pre)
            for i in range(Ny):
                seq.add_block(gx, adc)
                seq.add_block(gy)
                gx.amplitude = -gx.amplitude
            
            seq.add_block(gx_crush,gz_crush)
            
            if delay_tr>0:
                seq.add_block(make_delay(delay_tr))

seq.plot()

ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
plt.plot(time_axis, ktraj.T)
plt.plot(t_adc, ktraj_adc[0, :], '.')
plt.figure()
plt.plot(ktraj[0, :], ktraj[1, :], 'b')
plt.axis('equal')
plt.plot(ktraj_adc[0, :], ktraj_adc[1, :], 'r.')
plt.show()

seq.write('dwepi_pypulseq.seq')