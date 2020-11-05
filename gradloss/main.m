x = -10:0.05:10;  % fine grid
f0 = x.^2;  % field map
g = diff(f0)/dx;  % cycles/cm (Hz/cm)

dx = 0.3;   % voxel size (cm)

subplot(121); plot(x(2:end),g);
subplot(122); plot(x(2:end), sinc(dx*g))
