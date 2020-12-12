load result f0 fp  % see shim.jl
load fa            % acquired, after shimming

x = 11:50;
y = 11:50;
z = 13:40;
f0 = f0(x,y,z);
fp = fp(x,y,z);
fa = fa(x,y,z);
im(cat(1,f0, fp, fa, fp-fa), [-80 80])
colormap jet
h = colorbar;
h.TickLabels{end} = 'field map (Hz)';
title('original -- predicted -- acquired -- difference')


