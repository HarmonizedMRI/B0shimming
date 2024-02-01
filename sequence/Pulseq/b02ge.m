% b02ge.m

sysGE = toppe.systemspecs('maxGrad', 5, ...   % G/cm
    'maxSlew', 20, ...               % G/cm/ms
    'maxRF', 0.15, ...               % Gauss. Must be >= peak RF in sequence.
    'maxView', 120, ...              % Determines slice/view index in data file
    'rfDeadTime', 100, ...           % us
    'rfRingdownTime', 60, ...        % us
    'adcDeadTime', 40, ...           % us
    'psd_rf_wait', 148, ...          % RF/gradient delay (us)
    'psd_grd_wait', 156);            % ADC/gradient delay (us)

seq2ge('b0.seq', sysGE, 'b0.tar');
