% Low-resolution 3D B0 mapping sequence in Pulseq, optimized
% for identical execution on Siemens and GE scanners.
%
% This script creates the file 'b0.seq', that can be executed directly
% on Siemens MRI scanners using the Pulseq interpreter.
% The .seq file can also be converted to a .tar file that can be executed on GE
% scanners, see main.m

% Define experimental parameters.
% These are chosen such that the sequence can be executed identically
% (to us precision) on Siemens and GE systems.
sys=mr.opts('maxGrad', 30, 'gradUnit','mT/m', ...
            'maxSlew', 100, 'slewUnit', 'T/m/s', ...
            'rfDeadTime', 100e-6, ...
            'rfRingdownTime', 60e-6, ...
            'adcDeadTime', 40e-6, ...
            'adcRasterTime', 2e-6, ...
            'gradRasterTime', 10e-6, ...
            'blockDurationRaster', 10e-6, ...
            'B0', 3.0);
            %'riseTime', 400e-6, ...

fov = [240e-3 240e-3 240e-3];     % FOV (m)
Nx = 60; Ny = Nx; Nz = 20;        % Matrix size
dwell = 16e-6;                    % ADC sample time (s). For GE, must be multiple of 2us.
alpha = 5;                        % flip angle (degrees)

fatChemShift = 3.5*1e-6;          % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz
TE = 3e-3; % + [0 1/fatOffresFreq];                % fat and water in phase for both echoes
TR = 8e-3; %*[1 1];

nCyclesSpoil = 4;
Tpre = 1e-3;

alphaPulseDuration = 0.2e-3;

% Create a new sequence object
seq = mr.Sequence(sys);           

% Create non-selective pulse
[rf, rfDelay] = mr.makeBlockPulse(alpha/180*pi, sys, 'Duration', alphaPulseDuration);

% Define other gradients and ADC events
% We cut the RO gradient into two parts for the optimal spoiler timing
deltak = 1./fov;
Tread = Nx*dwell;

gxtmp = mr.makeTrapezoid('x', sys, ...  % temporary object
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread);
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gxtmp.area/2);
adc = mr.makeAdc(Nx, sys, ...
    'Duration', Tread,...
    'Delay', gxtmp.riseTime);
gx = mr.makeTrapezoid('x', sys, ...
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread + adc.deadTime); % extend flat time so we can split at end of ADC dead time
gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', Ny*deltak(2)/2);   % maximum PE1 gradient, positive amplitude
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', Nz*deltak(3)/2);   % maximum PE2 gradient, positive amplitude
gxSpoil = mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);
gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', gxSpoil.area);

gxSpoil = mr.makeExtendedTrapezoidArea('x', gx.amplitude, 0, gzSpoil.area, sys);
[gx, ~] = mr.splitGradientAt(gx, gx.riseTime + gx.flatTime);

% y/z PE steps
pe1Steps = ((0:Ny-1)-Ny/2)/Ny*2;
pe2Steps = ((0:Nz-1)-Nz/2)/Nz*2;

% Calculate timing
delayTE = ceil((TE - mr.calcDuration(rf) + mr.calcRfCenter(rf) + rf.delay - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = ceil((TR - mr.calcDuration(rf) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - mr.calcDuration(gxSpoil) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;

% Loop over phase encodes and define sequence blocks
% iZ < 0: Dummy shots to reach steady state
% iZ = 0: ADC is turned on and used to set receive gain calibration on GE scanners.
% iZ > 0: Actual image acquisition
for iZ = -1:Nz
    for ib = 1:40
        fprintf('\b');
    end
    if iZ > 0
        fprintf('Writing kz encode %d of %d', iZ, Nz);
    end
    if iZ < 0
        adcOn = false;
    else
        adcOn = true;
    end
    for iY = 1:Ny
        for c = 1:length(TE)
            % RF spoiling
            rf.phaseOffset = mod(117*(iY^2+iY+2)*pi/180, 2*pi);
            adc.phaseOffset = rf.phaseOffset;
            
            % Excitation
            seq.addBlock(rf, rfDelay);
            
            % Encoding
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gxPre, ...
                mr.scaleGrad(gyPre, pe1Steps(iY)), ...
                mr.scaleGrad(gzPre, pe2Steps(max(1,iZ))));
            if adcOn
                seq.addBlock(gx, adc);
            else
                seq.addBlock(gx);
            end
            seq.addBlock(mr.scaleGrad(gyPre, -pe1Steps(iY)), ...
                mr.scaleGrad(gzPre, -pe2Steps(max(1,iZ))), ...
                gxSpoil);
            seq.addBlock(mr.makeDelay(delayTR(c)), gzSpoil);
        end
    end
end

fprintf('\nSequence ready\n');

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% Visualise sequence and output for execution
Ndummy = length(TE)*2*Ny;
seq.plot('TimeRange',[Ndummy+1 Ndummy+10]*TR(1))

seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'b0');
seq.write('b0.seq', false);

% visualize the 3D k-space (only makes sense for low-res, otherwise one sees nothing)
if Nx<=32
    tic;
    [kfa,ta,kf]=seq.calculateKspacePP();
    toc
    figure;plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
end

return

%% create a smoothly rotating plot
if Nx<=16
    figure;plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
    kabsmax=max(abs(kf)')';
    kxyabsmax=max(kabsmax(1:2));
    kxyzabsmax=max(kabsmax);
    %axis([-kxyabsmax kxyabsmax -kxyabsmax kxyabsmax -kabsmax(3) kabsmax(3)])
    axis([-kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax])
    [caz,cel] = view;
    for caz_add=0:1:359 
        view(caz+caz_add,cel);
        drawnow;
    end
end

%% create a smoothly rotating plot (rotated to read along z)
if Nx<=16
    figure;plot3(kf(2,:),-kf(3,:),kf(1,:));
    hold on;plot3(kfa(2,:),-kfa(3,:),kfa(1,:),'r.');
    set(gca,'visible','off'); % hide axes
    set(gca, 'CameraViewAngle',get(gca, 'CameraViewAngle')); % freeze the view
    kabsmax=max(abs(kf)')';
    kxyabsmax=max(kabsmax(1:2));
    kxyzabsmax=max(kabsmax);
    %axis([-kxyabsmax kxyabsmax -kxyabsmax kxyabsmax -kabsmax(3) kabsmax(3)])
    %axis([-kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax])
    s1=1.2;    
    axis([ -kabsmax(2)*s1 kabsmax(2)*s1  -kabsmax(2)*s1 kabsmax(2)*s1 min(kf(1,:)) kabsmax(1)]);
    [caz,cel] = view;
    folder='kspace3d';
    mkdir(folder);
    for caz_add=0:1:359 
        view(caz+caz_add,cel);
        drawnow;
        print( '-r100', '-dpng', [folder '/frame_' num2str(caz_add,'%03d') '.png']);
        % use convert frame_???.png -gravity center -crop 300x300+0+0 +repage -delay 0.1 -loop 0 kspace_gre3d.gif
        % to create a GIF movie
    end
end
