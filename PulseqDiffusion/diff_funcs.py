#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:09:08 2019

@authors: Rita G. Nunes, Oscar Pena-Nogales 
"""

import math

import numpy as np
from pypulseq.calc_duration import calc_duration
from pypulseq.make_trap_pulse import make_trapezoid


def get_dirs(ndirs):
    """
    Retrieves list of gradient directions and number of non\-DWI images.

    Parameters
    ----------
    ndirs : integer
            number of diffusion directions to sample
    Returns
    -------
    g : numpy.ndarray
	gradient components for each direction (gx,gy,gz)
    nb0s: integer
	  number of non\-DWI volumes to acquire
    """

    if ndirs == 3:
        g = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        nb0s = 1
    elif ndirs == 6:
        g = np.zeros((3, ndirs))
        # obtained using gen_scheme (mrtrix)
        g = np.array([[-0.283341, -0.893706, -0.347862],
                      [-0.434044, 0.799575, -0.415074],
                      [0.961905, 0.095774, 0.256058],
                      [-0.663896, 0.491506, 0.563616],
                      [-0.570757, -0.554998, 0.605156],
                      [-0.198848, -0.056534, -0.978399]])
        nb0s = 1
    elif ndirs == 12:
        g = np.zeros((3, ndirs))
        # obtained using gen_scheme (mrtrix)
        g = np.array([[0.648514, 0.375307, 0.662249],
                      [-0.560493, 0.824711, 0.075496],
                      [-0.591977, -0.304283, 0.746307],
                      [-0.084472, -0.976168, 0.199902],
                      [-0.149626, -0.006494, -0.988721],
                      [-0.988211, -0.056904, -0.142130],
                      [0.864451, -0.274379, -0.421237],
                      [-0.173549, -0.858586, -0.482401],
                      [-0.039288, 0.644885, -0.763269],
                      [0.729809, 0.673235, -0.118882],
                      [0.698325, -0.455759, 0.551929],
                      [-0.325340, 0.489774, 0.808873]])
        nb0s = 1
    elif ndirs == 60:
        g = np.zeros((3, ndirs))
        # obtained using gen_scheme (mrtrix)
        g = np.array([[-0.811556, 0.245996, -0.529964],
                      [-0.576784, -0.313126, 0.754502],
                      [-0.167946, -0.899364, -0.403655],
                      [0.755699, -0.512113, -0.408238],
                      [0.116846, 0.962654, -0.244221],
                      [0.495465, 0.208081, 0.843337],
                      [0.901459, 0.385831, -0.196230],
                      [-0.248754, 0.420519, -0.872516],
                      [-0.047525, 0.444671, 0.894432],
                      [- 0.508593, 0.857494, -0.077699],
                      [0.693558, 0.614042, 0.376737],
                      [0.990394, -0.134781, -0.030898],
                      [0.019140, -0.684235, 0.729010],
                      [0.385221, - 0.339346, - 0.858166],
                      [- 0.440289, - 0.853536, 0.278609],
                      [0.680515, - 0.559825, 0.472752],
                      [- 0.146029, 0.872237, 0.466774],
                      [0.317352, 0.195118, - 0.928018],
                      [- 0.796280, - 0.129004, - 0.591013],
                      [- 0.711299, 0.249255, 0.657211],
                      [- 0.838383, - 0.538321, - 0.085587],
                      [0.202544, - 0.966710, 0.156357],
                      [-0.296747, - 0.476761, - 0.827430],
                      [0.545225, 0.637023, - 0.544914],
                      [-0.887097, 0.451265, 0.097048],
                      [0.034752, - 0.124211, 0.991647],
                      [0.469222, - 0.766720, - 0.438145],
                      [-0.948457, - 0.088803, 0.304209],
                      [-0.354311, 0.664176, - 0.658281],
                      [-0.462117, - 0.833550, - 0.302724],
                      [0.949202, - 0.022353, 0.313871],
                      [0.791248, 0.065354, - 0.607994],
                      [-0.004026, 0.992213, 0.124488],
                      [- 0.357034, 0.107223, - 0.927917],
                      [0.414504, 0.685422, 0.598651],
                      [-0.331743, - 0.552720, 0.764492],
                      [-0.749058, 0.641807, - 0.164305],
                      [0.238666, - 0.655691, - 0.716315],
                      [0.619125, 0.784982, 0.022082],
                      [0.123966, - 0.872070, 0.473419],
                      [-0.185240, 0.122014, 0.975089],
                      [-0.980282, - 0.189151, - 0.057166],
                      [0.637873, - 0.084335, 0.765510],
                      [-0.668960, 0.723273, 0.171375],
                      [-0.775822, - 0.381543, 0.502519],
                      [-0.636044, - 0.425049, - 0.644035],
                      [0.229220, 0.809364, - 0.540729],
                      [-0.538340, 0.531550, 0.653946],
                      [0.906105, - 0.354129, 0.231445],
                      [-0.166743, - 0.191681, - 0.967189],
                      [0.324636, - 0.927784, - 0.183922],
                      [0.551291, - 0.398459, 0.733014],
                      [0.537753, - 0.032690, - 0.842468],
                      [-0.306182, - 0.951456, - 0.031381],
                      [0.875976, 0.329209, 0.352545],
                      [0.902989, - 0.218632, - 0.369879],
                      [-0.456427, 0.801551, - 0.386251],
                      [0.089001, 0.716134, 0.692265],
                      [-0.714965, - 0.648438, 0.261444],
                      [0.076308, 0.420804, - 0.903936]])
        nb0s = 3
    else:
        g = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        nb0s = 1
    return g, nb0s


def calc_bval(G, delta, Delta, gdiff_rt):
    """
    Calculates the achieved diffusion-weighting (b-value in s/mm2)

    Parameters
    ----------
    G : float
        amplitude of the diffusion gradient (Hz)
    delta : float
            duration of the diffusion gradients (s)
    Delta : float
            time between start of the first and second diffusion gradients (s)
    gdiff_rt : float
               diffusion gradient ramp time (s)
    Returns
    -------
    bval : float
	   b-value in s/mm2
    """

    bval = (2 * math.pi * G) ** 2 * (
            (Delta - delta / 3) * (delta ** 2) + (gdiff_rt ** 3) / 30 - delta * (gdiff_rt ** 2) / 6)

    return bval


def calc_exact_bval(waveform, INV, short_f, seq):
    """
    Calculates the achieved exact diffusion-weighting (b-value in s/mm2)

    Parameters
    ----------
    waveform : numpy.ndarray
               diffusion weighting gradient waveform

    INV : numpy.ndarray
          array indicating whether the elements of waveform is in phase or in opposite phase compared to the
          beginning of the sequence (phase +1, opposite phase -1).

    short_f : int
              shortening factor to speed up compuation of the b-value.
              higher values will short the time resolution of the waveform.

    seq : sequence
          sequence

    Returns
    -------
    bval : float
           b-value in s/mm2
    """

    # These following lines might crush due to system overload - too large matrix C.
    # Compared to function 'calc_bval' this computes the exact b-value depending on the diffusion weighting gradient waveform.
    # It's independent of the diffusion-weighting sequence as long as the waveform and the INV variables are in accordance
    # with the corresponding sequence.

    # In that case:
    # 1) Reduce b-value to decrease TE
    # 2) Omit the computation of the b-value with full time resolution.
    # 3) Run code in a better computer (larger RAM)
    n = len(waveform)

    C = np.tril(np.ones([n, n]))
    C2 = np.matmul(np.transpose(C), C)

    F = waveform * INV * short_f * seq.grad_raster_time
    bval = (2 * math.pi) ** 2 * np.matmul(np.matmul(np.transpose(F), C2),
                                          F) * short_f * seq.grad_raster_time * 1e-6  # [s/mm2]
    return bval


def opt_TE_bv_SE(bvalue_Dict, grads_times_Dict, seq_sys_Dict):
    """
    Obtain optimal TE and diffusion-weighting gradient waveforms for twice-refocued spin echo (TRSE) DW sequence.

    Parameters
    ----------
    bvalue_Dict : dictionary
                  b-value related parameters needed in the TE optimization.

    grads_times_Dict : dictionary
                       gradient times related parameters needed in the TE optimization.

    seq_sys_Dict : dictionary
                   sequence and system related parameters needed in the TE optimization.

    Returns
    -------
    n_TE : int
	   echo time in integers units.

    bval : float
           b-value in s/mm2

    gdiff : pypulseq gradient
	    diffusion-weighting gradient waveform.

    n_delay_te1 : int
	          delay between RF90 and RF180 in integers units.

    n_delay_te2 : int
	          delay between RF180 and EPI readout in integers units.
    """

    # Description
    # For S&T monopolar waveforms.
    # From an initial TE, check we satisfy all constraints -> otherwise increase TE.
    # Once all constraints are okay -> check b-value, if it is lower than the target one -> increase TE
    # Looks time-inefficient but it is fast enough to make it user-friendly.

    # We need to compute the exact time sequence. For the normal SE-MONO-EPI sequence micro second differences
    # are not important, however, if we wanna import external gradients the allocated time for them needs to
    # be the same, and thus exact timing is mandatory. With this in mind, we establish the following rounding rules:
    # Duration of RFs + spoiling, and EPI time to the center of the k-space is always math.ceil().

    bvalue = bvalue_Dict["bvalue"]  # Target b-value.
    nbvals = bvalue_Dict["nbvals"]  # number of b-values.
    gscl = bvalue_Dict["gscl"]  # gradient scaling.

    n_rf90r = grads_times_Dict["n_rf90r"]  # Half right duration of RF90 (integer).
    n_rf180r = grads_times_Dict["n_rf180r"]  # Half right duration of the RF180 (integer).
    n_rf180l = grads_times_Dict["n_rf180l"]  # Half left duration of the RF180 (integer).
    gz_spoil = grads_times_Dict["gz_spoil"]  # Spoil gradient to be placed around the RF180.
    gz180 = grads_times_Dict["gz180"]  # RF180 simultaneous gradient.
    n_duration_center = grads_times_Dict[
        "n_duration_center"]  # Time needed to achieve the center of the k-space (integer).

    seq = seq_sys_Dict["seq"]  # Sequence
    system = seq_sys_Dict["system"]  # System
    i_raster_time = seq_sys_Dict["i_raster_time"]  # Manually inputed inverse raster time.

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
        n_delay_te2 = n_tINV - n_rf180r - n_duration_center

        # Waveform Ramp time
        n_gdiff_rt = math.ceil(system.max_grad / system.max_slew / seq.grad_raster_time)

        # Select the shortest available time
        n_gdiff_delta = min(n_delay_te1, n_delay_te2)
        n_gdiff_Delta = n_delay_te1 + 2 * math.ceil(calc_duration(gz_spoil) / seq.grad_raster_time) + math.ceil(
            calc_duration(gz180) / seq.grad_raster_time)

        gdiff = make_trapezoid(channel='x', system=system, amplitude=system.max_grad,
                               duration=n_gdiff_delta / i_raster_time)

        # delta only corresponds to the rectangle.
        n_gdiff_delta = n_gdiff_delta - 2 * n_gdiff_rt

        bv = calc_bval(system.max_grad, n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time,
                       n_gdiff_rt / i_raster_time)
        bval = bv * 1e-6

    # Show final TE and b-values:
    print("TE:", round(n_TE / i_raster_time * 1e3, 2), "ms")
    for bv in range(1, nbvals + 1):
        print(round(
            calc_bval(system.max_grad * gscl[bv], n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time,
                      n_gdiff_rt / i_raster_time) * 1e-6, 2), "s/mm2")

    return n_TE, bval, gdiff, n_delay_te1, n_delay_te2


def opt_TE_bv_TRSE(bvalue_Dict, grads_times_Dict, seq_sys_Dict):
    """
    Obtain optimal TE and diffusion-weighting gradient waveforms for twice-refocued spin echo (TRSE) DW sequence.

    Parameters
    ----------
    bvalue_Dict : dictionary
                  b-value related parameters needed in the TE optimization.

    grads_times_Dict : dictionary
                       gradient times related parameters needed in the TE optimization.

    seq_sys_Dict : dictionary
                   sequence and system related parameters needed in the TE optimization.

    Returns
    -------
    TE : float
       	 echo time (s).

    bval : float
     	   b-value in s/mm2

    waveform : float
               diffusion-weighting gradient waveform.

    gscl_max : float
               gradient scale factor of the diffusion-weighting gradient waveform.

    d1 : float
         duration of the first diffusion-weighting gradient lobe between the RF90 and the first RF180 (s).

    d2 : float
         duration of the second diffusion-weighting gradient lobe after the first RF180 (s).

    d3 : float
         duration of the third diffusion-weighting gradient lobe before the second RF180 (s).

    d4 : float
         duration of the last diffusion-weighting gradient lobe after the second RF180 (s).
    """

    # Description
    # For TRSE waveforms.
    # From an initial TE, check we satisfy all constraints -> otherwise increase TE.
    # Once all constraints are okay -> check b-value, if it is lower than the target one -> increase TE
    # Finally, with TRSE low b-values can not be acquired, thus proper scaling is needed.
    # Looks time-inefficient but it is fast enough to make it user-friendly.

    bvalue = bvalue_Dict["bvalue"]  # Target b-value.
    nbvals = bvalue_Dict["nbvals"]  # number of b-values.
    gscl = bvalue_Dict["gscl"]  # gradient scaling.

    rf180 = grads_times_Dict["rf180"]  # RF180
    rf180_center_with_delay = grads_times_Dict["rf180_center_with_delay"]  # time to the center of the RF180
    rf_center_with_delay = grads_times_Dict["rf_center_with_delay"]  # time to the center of the RF90
    gz_spoil_1 = grads_times_Dict["gz_spoil_1"]  # Spoil gradient to be placed around the first RF180
    gz_spoil_2 = grads_times_Dict["gz_spoil_2"]  # Spoil gradient to be placed around the second RF180
    gz = grads_times_Dict["gz"]  # RF90 gradient
    duration_center = grads_times_Dict["duration_center"]  # Time needed to achieve the center of the k-space.
    pre_time = grads_times_Dict["pre_time"]  # (s)

    seq = seq_sys_Dict["seq"]  # Sequence
    system = seq_sys_Dict["system"]  # System

    # Find minimum TE considering the readout times and RF + spoil gradient durations
    TE = 80e-3  # [s] It's a reasonable estimate for TRSE.
    delay_tela = -1  # Large TE/2
    while delay_tela <= 0:
        TE = TE + 0.02e-3  # [ms]
        delay_tela = math.ceil((TE / 2 + (- calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_2) - \
                                          duration_center) + (- calc_duration(gz) + rf_center_with_delay - \
                                                              pre_time - calc_duration(
                    gz_spoil_1) - rf180_center_with_delay)) / seq.grad_raster_time) * seq.grad_raster_time

    # Sice there are 3 equations and 4 parameters (TGReese2002 - https://doi.org/10.1002/mrm.10308)
    # One parameter can be tuned, while the other 3 are then fixed. This is why we fix d4.
    d4 = 3e-3  # [ms]

    # Find minimum TE for the target d4
    # Waveform Ramp time
    gdiff_rt = math.ceil(system.max_grad / system.max_slew / seq.grad_raster_time) * seq.grad_raster_time
    d1 = -1
    while d1 <= 2 * gdiff_rt:  # Include this condition to have trapezoids everywhere
        TE = TE + 2 * seq.grad_raster_time  # [ms]
        # Large TE/2
        delay_tela = math.ceil((TE / 2 + (- calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_2) - \
                                          duration_center) + (- calc_duration(gz) + rf_center_with_delay - \
                                                              pre_time - calc_duration(
                    gz_spoil_1) - rf180_center_with_delay)) / seq.grad_raster_time) * seq.grad_raster_time

        # Short TE/2 (Time between RF180s)
        delay_tes = math.ceil((TE / 2 - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_1) - \
                               calc_duration(
                                   gz_spoil_2) - rf180_center_with_delay) / seq.grad_raster_time) * seq.grad_raster_time

        d1 = math.ceil((delay_tela - d4) / seq.grad_raster_time) * seq.grad_raster_time

    # Find minimum TE for the target b-value
    bvalue_tmp = 0

    # We initially substract 2 times the raster time because we don't know if it will achieve the desired b-value from the beginning.
    # Note that for low b-values, the TE is the same for all of them due to the large number of elements of this sequence.
    TE = TE - 2 * seq.grad_raster_time
    while bvalue_tmp < np.max(bvalue):
        # The following scalar multiplication is just to increase speed
        TE = TE + int(math.ceil(np.max(bvalue) / 500)) * 4 * seq.grad_raster_time  # [ms]

        # Large TE/2
        delay_tela = math.ceil((TE / 2 + (- calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_2) - \
                                          duration_center) + (- calc_duration(gz) + rf_center_with_delay - \
                                                              pre_time - calc_duration(
                    gz_spoil_1) - rf180_center_with_delay)) / seq.grad_raster_time) * seq.grad_raster_time

        # Short TE/2 (Time between RF180s)
        delay_tes = math.ceil((TE / 2 - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil_1) - \
                               calc_duration(
                                   gz_spoil_2) - rf180_center_with_delay) / seq.grad_raster_time) * seq.grad_raster_time

        d1 = math.ceil((delay_tela - d4) / seq.grad_raster_time) * seq.grad_raster_time
        d3 = math.ceil(((d1 + delay_tes - d4) / 2) / seq.grad_raster_time) * seq.grad_raster_time
        d2 = math.ceil((delay_tes - d3) / seq.grad_raster_time) * seq.grad_raster_time

        # Due to the complexity of the sequence and that I have not found its b-value, I implement it by its definition.
        # However, for simplicity I compute each gradient lobe as ideal rectangular pulses.

        # d1 start after the RF90
        n1 = int(math.ceil((calc_duration(
            gz) - rf_center_with_delay + pre_time + seq.grad_raster_time) / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
        nd1 = int(math.ceil(d1 / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)

        # d2 start after the first RF180
        n2 = int(math.ceil(((n1 + nd1) * seq.grad_raster_time + 2 * calc_duration(gz_spoil_1) + calc_duration(
            rf180) + seq.grad_raster_time) / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
        nd2 = int(math.ceil(d2 / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)

        # d3 starts after d2
        n3 = int(math.ceil(((
                                    n2 + nd2) * seq.grad_raster_time + seq.grad_raster_time) / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
        nd3 = int(math.ceil(d3 / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)

        # d4 starts after the second RF180
        n4 = int(math.ceil(((n2 + nd2 + nd3) * seq.grad_raster_time + 2 * calc_duration(gz_spoil_2) + calc_duration(
            rf180) + seq.grad_raster_time) / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
        nd4 = int(math.ceil(d4 / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)

        # Due to the large TE of TRSE and the high resolution we might need to shrink the following array.
        # To speed up computation we calculate the b-value with short_f times lower of the time resolution.

        # Shortening factor:
        short_f = 5
        n = int(np.ceil(TE / seq.grad_raster_time) / short_f)
        waveform = np.zeros(n)

        # Compose waveform
        waveform[math.ceil(n1 / short_f):math.floor((n1 + nd1 + 1) / short_f)] = system.max_grad
        waveform[math.ceil(n2 / short_f):math.floor((n2 + nd2 + 1) / short_f)] = -system.max_grad
        waveform[math.ceil(n3 / short_f):math.floor((n3 + nd3 + 1) / short_f)] = system.max_grad
        waveform[math.ceil(n4 / short_f):math.floor((n4 + nd4 + 1) / short_f)] = -system.max_grad

        # Include ramps
        nrt = int(math.floor(math.floor(
            system.max_grad / system.max_slew / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time / short_f))
        ramp_values = np.linspace(1, nrt, nrt) * seq.grad_raster_time * system.max_slew * short_f
        # d1
        if n1 % short_f == 0:
            waveform[math.ceil(n1 / short_f):math.ceil(n1 / short_f) + nrt] = ramp_values
        else:
            waveform[math.ceil(n1 / short_f):math.ceil((n1 + 1) / short_f) + nrt] = ramp_values

        waveform[math.ceil((n1 + nd1 - 1) / short_f - nrt):math.ceil((n1 + nd1 - 1) / short_f)] = np.flip(ramp_values)

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
        nRF180_1 = int(math.ceil((n2 * seq.grad_raster_time - calc_duration(
            rf180) + rf180_center_with_delay - calc_duration(
            gz_spoil_1)) / short_f / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
        nRF180_2 = int(math.ceil((n4 * seq.grad_raster_time - calc_duration(
            rf180) + rf180_center_with_delay - calc_duration(
            gz_spoil_2)) / short_f / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
        INV = np.ones(n)
        INV[nRF180_1:nRF180_2 + 1] = -1

        bvalue_tmp = calc_exact_bval(waveform, INV, short_f, seq)

    print("b-value low time resolution:", round(bvalue_tmp, 2), "s/mm2")

    # To know the exact b-value - repeat the above with full time resolution
    short_f = 1
    n = int(np.ceil(TE / seq.grad_raster_time) / short_f)
    waveform = np.zeros(n)

    # Compose waveform
    waveform[math.ceil(n1 / short_f):math.floor((n1 + nd1 + 1) / short_f)] = system.max_grad
    waveform[math.ceil(n2 / short_f):math.floor((n2 + nd2 + 1) / short_f)] = -system.max_grad
    waveform[math.ceil(n3 / short_f):math.floor((n3 + nd3 + 1) / short_f)] = system.max_grad
    waveform[math.ceil(n4 / short_f):math.floor((n4 + nd4 + 1) / short_f)] = -system.max_grad

    # Include ramps
    nrt = int(math.floor(math.floor(
        system.max_grad / system.max_slew / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time / short_f))
    ramp_values = np.linspace(0, nrt, nrt) * seq.grad_raster_time * system.max_slew * short_f
    # d1
    if n1 % short_f == 0:
        waveform[math.ceil(n1 / short_f):math.ceil(n1 / short_f) + nrt] = ramp_values
    else:
        waveform[math.ceil(n1 / short_f):math.ceil((n1 + 1) / short_f) + nrt] = ramp_values

    if (n1 + nd1) % short_f == 0:
        waveform[math.ceil((n1 + nd1) / short_f - nrt) + 1:math.floor((n1 + nd1 + 1) / short_f)] = np.flip(ramp_values)
    else:
        waveform[math.ceil((n1 + nd1) / short_f - nrt):math.floor((n1 + nd1 + 1) / short_f)] = np.flip(ramp_values)

    # d2
    if n2 % short_f == 0:
        waveform[math.ceil(n2 / short_f):math.ceil(n2 / short_f) + nrt] = -ramp_values
    else:
        waveform[math.ceil(n2 / short_f):math.ceil((n2 + 1) / short_f) + nrt] = -ramp_values

    if (n2 + nd2) % short_f == 0:
        waveform[math.ceil((n2 + nd2) / short_f - nrt) + 1:math.floor((n2 + nd2 + 1) / short_f)] = -np.flip(ramp_values)
    else:
        waveform[math.ceil((n2 + nd2) / short_f - nrt):math.floor((n2 + nd2 + 1) / short_f)] = -np.flip(ramp_values)

    # d3
    if n3 % short_f == 0:
        waveform[math.ceil(n3 / short_f):math.ceil(n3 / short_f) + nrt] = ramp_values
    else:
        waveform[math.ceil(n3 / short_f):math.ceil((n3 + 1) / short_f) + nrt] = ramp_values

    if (n3 + nd3) % short_f == 0:
        waveform[math.ceil((n3 + nd3) / short_f - nrt) + 1:math.floor((n3 + nd3 + 1) / short_f)] = np.flip(ramp_values)
    else:
        waveform[math.ceil((n3 + nd3) / short_f - nrt):math.floor((n3 + nd3 + 1) / short_f)] = np.flip(ramp_values)

    # d4
    if n4 % short_f == 0:
        waveform[math.ceil(n4 / short_f):math.ceil(n4 / short_f) + nrt] = -ramp_values
    else:
        waveform[math.ceil(n4 / short_f):math.ceil((n4 + 1) / short_f) + nrt] = -ramp_values

    if (n4 + nd4) % short_f == 0:
        waveform[math.ceil((n4 + nd4) / short_f - nrt) + 1:math.floor((n4 + nd4 + 1) / short_f)] = -np.flip(ramp_values)
    else:
        waveform[math.ceil((n4 + nd4) / short_f - nrt):math.floor((n4 + nd4 + 1) / short_f)] = -np.flip(ramp_values)

    # Prepare Integral
    nRF180_1 = int(math.ceil((
                                     n2 * seq.grad_raster_time - calc_duration(
                                 rf180) + rf180_center_with_delay - calc_duration(
                                 gz_spoil_1)) / short_f / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
    nRF180_2 = int(math.ceil((
                                     n4 * seq.grad_raster_time - calc_duration(
                                 rf180) + rf180_center_with_delay - calc_duration(
                                 gz_spoil_2)) / short_f / seq.grad_raster_time) * seq.grad_raster_time / seq.grad_raster_time)
    INV = np.ones(n)
    INV[nRF180_1:nRF180_2 + 1] = -1

    bval = calc_exact_bval(waveform, INV, short_f, seq)

    print("b-value high time resolution:", round(bval, 2), "s/mm2")

    # Scale gradients amplitude accordingly.
    if bval > np.max(bvalue):
        gscl_max = np.sqrt(np.max(bvalue) / bval)
    else:
        gscl_max = 1

    # Show final TE and b-values:
    print("Final times:")
    print("TE:", round(TE * 1e3, 2), "ms")
    for bv in range(1, nbvals + 1):
        bval_tmp = calc_exact_bval(waveform * gscl_max * gscl[bv], INV, short_f, seq)
        print(round(bval_tmp, 2), "s/mm2")

    return TE, bval, waveform, gscl_max, d1, d2, d3, d4
