% signal loss model is sinc(dx*g)
% dx = voxel size (cm)
% g = B0 gradient (Hz/cm)

% synthesize g
x = -10:0.05:10;  % fine grid
dx = x(2)-x(1);
f0 = x.^2;        % field map
g = diff(f0)/dx;  % cycles/cm (Hz/cm)

% plot resulting signal loss in voxel of size dx
dx = 0.2;   % voxel size (cm)
subplot(121); plot(x(2:end),g);
subplot(122); plot(x(2:end), sinc(dx*g))
