
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

im(cat(1,f0.*maskglobal, fa.*masklocal), [-100 100]); colormap jet; 
h = colorbar;
h.TickLabels{end} = 'Hz';
t = sprintf('head phantom\nafter built-in shimming -- after localized shimming ');
t = sprintf('%s\n nrms = %.2f, %.2f', t, 1.0, nrmsfa);
disp(t)
print -dpng redhead.png


return;

mask = maskloc;

x = 1:60;
y = 1:60;

lims = [-200 200];


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


