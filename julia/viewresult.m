
load fa fa           % acquired, after shimming
load result f0 fp mask % see shim.jl

figure;
x = 11:50;
y = 11:50;
z = 12:41;
f0 = f0(x,y,z);
fp = fp(x,y,z);
fa = fa(x,y,z);
im(cat(1,f0, fp, fa, fp-fa), [-80 80])
colormap jet
h = colorbar;
h.TickLabels{end} = 'Hz';
title(sprintf('UM MR750 3T (BIRB) 11-Dec-2020\noriginal -- predicted -- acquired -- difference'))

% 
load result f0 fp mask % see shim.jl
load redhead_roi roi3d mask z

x = 1:60;
y = 1:60;

lims = [-200 200];

figure; im(f0, lims); colormap jet; colorbar; title('field map (baseline)');
print -dpng baseline.png

figure;
im(cat(1,f0(x,y,z), fplin(x,y,z), fphos(x,y,z)), lims)
colormap jet
h = colorbar;
h.TickLabels{end} = 'Hz';
t = sprintf('red head phantom\noriginal -- predicted (linear) -- predicted (hos) ');
rmsf0 = norm(f0orig(mask));
nrmslin = norm(fplin(mask))/rmsf0
nrmshos = norm(fphos(mask))/rmsf0
t = sprintf('%s\n nrms_{lin}, nrms_{hos} = %.2f, %.2f', t, nrmslin, nrmshos);
title(t)
print -dpng shimmed.png

figure; im(mask);
print -dpng mask.png


