function f = evalspharm(x,y,z,l,m)
%
% Evaluate solution to Laplace's equation at position (x,y,z),
% for spherical harmonic degree l and order m.
%
% In:
% x/y/z    [N 1]     spatial locations
% l                  degree of Legendre polynomial
% m                  order (0 <= m <= l)
%
% Out:
% f        [N 1]     r^l * exp(1i*m*ph) * legendre(l, cos(th))   (complex)

if isstring(x)
	sub_test();
	return;
end

% Spherical coordinate system as defined in https://en.wikipedia.org/wiki/Spherical_harmonics
% and Romeo and Hoult MRM 1984
[ph,el,r] = cart2sph(x,y,z);  % ph = azimuth, el = elevation, r = radius
th = pi/2 - el; % polar angle

lp = legendre(l, cos(th));   % [l+1 N]

f = r.^l .* exp(1i*m*ph) .* lp(m+1,:)';

return


function sub_test

% define spatial locations 
n = 36;
r = linspace(-10,10,n);
[X,Y,Z] = meshgrid(r,r,r);
R = sqrt([X.^2 + Y.^2 + Z.^2]);
mask = R > 0;

f = 0*X;
f_mz = f;

lmax = 2;  % degree
nr = lmax+1;
nc = nr;
ii = 1;
for l = 0:lmax
	for m = 0:l
		fm = evalspharm(X(mask), Y(mask), Z(mask), l, m);
		f(mask) = fm;
		figure(1); subplot(nr,nc,l*nc+m+1); im(real(f)); title(sprintf('l,m = %d,%d', l, m'));
		colorbar; colormap jet;
		figure(2); subplot(nr,nc,l*nc+m+1); im(imag(f)); title(sprintf('l,m = %d,%d', l, -m'));
		colorbar; colormap jet;
		%figure(1); subplot(nr,nc,l*nc+m+1); im(cat(1,f(:,:,end/2), squeeze(f(:,end/2,:)), squeeze(f(end/2,:,:)))); 

		% compare with Maxim's Cartesian expressions
		fm_mz = poly_harmonic(l, m, X(mask), Y(mask), Z(mask));
		f_mz(mask) = fm_mz;
		figure(3); subplot(nr,nc,l*nc+m+1); im(f_mz); title(sprintf('mz, l,m = %d,%d', l, m'));
		colorbar; colormap jet;
		fm_mz = poly_harmonic(l, -m, X(mask), Y(mask), Z(mask));
		f_mz(mask) = fm_mz;
		figure(4); subplot(nr,nc,l*nc+m+1); im(f_mz); title(sprintf('mz, l,m = %d,%d', l, -m'));
		colorbar; colormap jet;
	end
end

return


