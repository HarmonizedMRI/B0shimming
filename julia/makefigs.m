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
f0(:,:,end+1) = -inf;
fa(:,:,end+1) = -inf;
masklocal(:,:,end+1) = 0;
maskglobal(:,:,end+1) = 0;

figure;
r1 =  max(abs(f0(maskglobal)));
r1 = 170;
load mask_readhead_local_continuous
roi = mask_readhead_local_continuous;
for iz = 1:size(roi,3)
	roiedge(:,:,iz) = edge(roi(:,:,iz));
end
im(cat(1,f0.*maskglobal, fa.*maskglobal), r1*[-1 1]); colormap jet; 
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
r2 = max(abs(f0(mask)));
r2 = 150;
f0(~mask) = -r2;
fa(~mask) = -r2;
im(cat(1,f0, fa), r2*[-1 1]); colormap jet; 
m = colormap;
m(1,:) = 0;   % black background
colormap(m);
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
z = 1:3:(nz-15);
f0 = f0(:,:,z);
fa = fa(:,:,z);
fp = fp(:,:,z);
mask = mask(:,:,z);

figure;
r3 = max(abs(f0(mask)));
%im('row', 3, cat(1,f0.*mask, fa.*mask, fp.*mask), r*[-1 1]); colormap jet; 
im('row', 3, cat(1,f0.*mask), r3*[-1 1]); colormap jet; 
h = colorbar; h.TickLabels{end} = 'Hz'; axis off; title '';
print -dpng mgh1.png
t = sprintf('FBIRN phantom\nBaseline -- after built-in shim -- predicted best shims');
title(t)
n = sum(mask(:));
[norm(fp(mask))/sqrt(n) norm(fa(mask))/sqrt(n)]


figure;
r4 = max(abs(fa(mask)));
r4 = 15;
im('row', 3, cat(1,fa.*mask, fp.*mask), r4*[-1 1]); colormap jet; 
h = colorbar; h.TickLabels{end} = 'Hz'; axis off; title '';
print -dpng mgh2.png

