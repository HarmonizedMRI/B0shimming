% UM data
load f0_redhead.mat
masklocal = mask;
load ~/tmp/fa_redhead.mat
maskglobal = mask;

x = 7:52;
y = 1:60;
z = 21:28;
f0 = f0(x,y,z);
fa = fa(x,y,z);
maskglobal = maskglobal(x,y,z);
masklocal = masklocal(x,y,z);

rmsf0 = norm(f0(masklocal));
nrmsfa = norm(fa(masklocal))/rmsf0;

f0(~maskglobal) = -inf;
fa(~maskglobal) = -inf;
fa(~masklocal) = -inf;
f0(:,:,end+1) = -inf;
fa(:,:,end+1) = -inf;
masklocal(:,:,end+1) = 0;
maskglobal(:,:,end+1) = 0;

figure;
im(cat(1,f0.*maskglobal, fa.*masklocal), [-100 100]); colormap jet; 
h = colorbar;
h.TickLabels{end} = 'Hz';
axis off
title ''
print -dpng redhead.png
t = sprintf('head phantom\nafter built-in shimming -- after localized shimming ');
t = sprintf('%s\n nrms = %.2f, %.2f', t, 1.0, nrmsfa);
title(t)


clear all
load f0_jar.mat
load fa_jar.mat

x = 13:50;
y = 13:50;
z = 18:3:42;
f0 = f0(x,y,z);
fa = fa(x,y,z);
mask = mask(x,y,z);

f0(~mask) = -inf;
fa(~mask) = -inf;
%f0(:,:,end+1) = -inf;
%fa(:,:,end+1) = -inf;
%mask(:,:,end+1) = 0;

figure;
im(cat(1,f0.*mask, fa.*mask), [-160 100]); colormap jet; 
h = colorbar;
h.TickLabels{end} = 'Hz';
axis off
title ''
print -dpng jar.png
t = sprintf('Agar jar\nafter built-in shimming -- after proposed shimming ');
title(t)


% MGH data
clear all
load f0_mgh     % f0, fov, mask
load fa_mgh     % fa, fov, mask
load result;    % fp
f0(~mask) = -inf;
fa(~mask) = -inf;
fp(~mask) = -inf;

[nx ny nz] = size(f0);
z = 1:4:(nz-10);
f0 = f0(:,:,z);
fa = fa(:,:,z);
fp = fp(:,:,z);
mask = mask(:,:,z);

figure;
im(cat(1,f0.*mask, fp.*mask, fa.*mask), [-20 20]); colormap jet; 
h = colorbar;
h.TickLabels{end} = 'Hz';
axis off
title ''
print -dpng mgh.png
t = sprintf('FBIRN phantom\nBaseline -- optimized shims (predicted) -- after manual shimming');
title(t)
rmsfp = norm(fp(mask))
rmsfa = norm(fa(mask))

