# writeB0.m
#
# Create the file 'b0.seq', a 3D spoiled GRE B0 mapping sequence.
# Edit this script as desired.

import math
import matplotlib.pyplot as plt
import numpy as np

import pypulseq as mr
import sigpy.mri.rf as rf
import pp_visualization as vis



# RF/gradient delay (sec). 
# Conservative choice that should work across all GE scanners.
psd_rf_wait = 200e-6;  

# Set system limits
system = mr.Opts(max_grad=40, grad_unit='mT/m',
                 max_slew=150, slew_unit='T/m/s', 
                 rf_dead_time=100e-6,
                 rf_ringdown_time=60e-6 + psd_rf_wait, 
                 adc_dead_time=40e-6,
                 adc_raster_time = 2e-6, 
                 block_duration_raster = 10e-6, 
                 B0 = 3.0)


# Create a new sequence object
seq=mr.Sequence(system)


# Acquisition parameters
fov = [240e-3, 240e-3, 240e-3]   # FOV (m)
Nx = 60; Ny = Nx; Nz = 60      # Matrix size
dwell = 10e-6                  # ADC sample time (s)
alpha = 3                      # flip angle (degrees)
fatChemShift = 3.5e-6          # 3.5 ppm
fatOffresFreq = system.gamma*system.B0*fatChemShift  # Hz
TE = 1/fatOffresFreq*[1, 2]     # fat and water in phase for both echoes
TR = 6e-3*[1, 1]                # constant TR
nCyclesSpoil = 2               # number of spoiler cycles
Tpre = 0.5e-3                  # prephasing trapezoid duration
rfSpoilingInc = 117            # RF spoiling increment


# Create non-selective pulse
rf = mr.make_block_pulse(flip_angle = alpha/180*np.pi, system=system, duration = 0.2e-3)

# Define other gradients and ADC events
# Cut the redaout gradient into two parts for optimal spoiler timing
deltak = 1./fov;
Tread = Nx*dwell;

# maximum PE1 gradient, max positive amplitude
gyPre = mr.make_trapezoid(channel='y', system=system, area=Ny*deltak(2)/2, duration=Tpre)

# maximum PE2 gradient, max positive amplitude
gzPre = mr.make_trapezoid(channel='z', system=system, area= Nz*deltak(3)/2, duration=Tpre)

# readout trapezoid, temporary object
gxtmp = mr.make_trapezoid(channel='x', system=system, amplitude=Nx*deltak(1)/Tread, flat_time=Tread)

gxPre = mr.make_trapezoid(channel='x', system=system, area=-gxtmp.area/2, duration = Tpre)

adc = mr.make_adc(num_samples=Nx, system=system, duration=Tread, delay= gxtmp.rise_time)

# extend flat time so we can split at end of ADC dead time
gxtmp2 = mr.make_trapezoid(channel='x', system=system, amplitude=Nx*deltak(1)/Tread, flat_time=Tread + adc.dead_time)   

gx, _ = mr.split_gradient_at(gxtmp2, gxtmp2.rise_time + gxtmp2.flat_time) #????
#[gx, ~] = mr.split_gradient_at(gxtmp2, gxtmp2.rise_time + gxtmp2.flat_time) #????
#%gx = gxtmp;

gzSpoil = mr.make_trapezoid(channel='z', system=system, area=Nx*deltak(1)*nCyclesSpoil)
gxSpoil = mr.make_extended_trapezoid_area(channel='x', gxtmp.amplitude,0, gzSpoil.area,system)
#%gxSpoil = mr.makeTrapezoid('x', sys, ...
#%    'Area', Nx*deltak(1)*nCyclesSpoil);

# y/z PE steps
pe1Steps = ((0:Ny-1)-Ny/2)/Ny*2;
pe2Steps = ((0:Nz-1)-Nz/2)/Nz*2;

# Calculate timing
TEmin = rf.shape_dur/2 + rf.ringdown_time + mr.calc_duration(gxPre) + adc.delay + Nx/2*dwell
delayTE = ceil((TE-TEmin)/seq.grad_raster_time)*seq.grad_raster_time
TRmin = mr.calc_duration(rf) + delayTE + mr.calc_duration(gxPre) + mr.calc_duration(gx) + mr.calc_duration(gxSpoil);
delayTR = ceil((TR-TRmin)/seq.gradRaster_time)*seq.gradRaster_time;

#% Loop over phase encodes and define sequence blocks
#% iZ < 0: Dummy shots to reach steady state
#% iZ = 0: ADC is turned on and used for receive gain calibration on GE scanners
#% iZ > 0: Image acquisition

nDummyZLoops = 1
rf_phase = 0
rf_inc = 0

for iZ = -nDummyZLoops:Nz
    isDummyTR = iZ < 0;

    for iY = 1:Ny
        #% Turn on y and z prephasing lobes, except during dummy scans and
        #% receive gain calibration (auto prescan)
        yStep = (iZ > 0) * pe1Steps(iY);
        zStep = (iZ > 0) * pe2Steps(max(1,iZ));

        for c = 1:length(TE)

            #% RF spoiling
            rf.phase_offset = rf_phase/180*pi;
            adc.phase_offset = rf_phase/180*pi;
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
            
            #% Excitation
            #% Mark start of segment (block group) by adding label.
            #% Subsequent blocks in block group are NOT labelled.
            seq.add_block(rf, mr.make_label('SET', 'LIN', 2-isDummyTR));
            
            #% Encoding
            seq.add_block(mr.make_delay(delayTE(c)));
            seq.add_block(gxPre, ...
                mr.scale_grad(gyPre, yStep), ...
                mr.scale_grad(gzPre, zStep));
            if isDummyTR
                seq.add_block(gx);
            else
                seq.add_block(gx, adc);
            end

            #% rephasing/spoiling
            seq.add_block(gxSpoil, ...
                mr.scale_grad(gyPre, -yStep), ...
                mr.scale_grad(gzPre, -zStep));
            seq.add_block(mr.make_delay(delayTR(c)));
        end
    end
end

print('Sequence ready\n');

### check whether the timing of the sequence is correct
ok, error_report = seq.check_timing() 
if ok:
    print('Timing check passed successfully')
else:
    print('Timing check failed. Error listing follows:')
    [print(e) for e in error_report]


#% Output for execution
seq.set_definition('FOV', fov)
seq.set_definition('Name', 'b0')
seq.write('b0.seq')


# ======
# VISUALIZATION
# ======

seq.plot()
# Calculate trajectory
ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

# Plot k-spaces
time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
plt.figure()
plt.plot(time_axis, ktraj.T)  # Plot entire k-space trajectory
plt.plot(t_adc, ktraj_adc[0], '.')  # Plot sampling points on kx-axis
plt.figure()
plt.plot(ktraj[0], ktraj[1], 'b', ktraj_adc[0], ktraj_adc[1], 'r.')  # 2D plot
plt.axis('equal')
plt.show()

