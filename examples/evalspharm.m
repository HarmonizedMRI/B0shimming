function f = evalbasis(x,y,z,l,m)
%
% In:
% x/y/z    [N 1]
% l                 degree
% m                 order
%
% Out:
% f 

if isstring(x)
	sub_test();
	return;
end

x = x(:)';
y = y(:)';
z = z(:)';

r = norm([x y z]);   % radius

% Follow spherical coordinate system as defined in https://en.wikipedia.org/wiki/Spherical_harmonics
[ph,el,r] = cart2sph(x,y,z);  % ph = azimuth, el = elevation, r = radius
th = el + pi/4; % polar angle

lp = legendre(l, cos(th));   % [l+1 N]

f = real(r.^l .* exp(1i*m*ph) .* lp(m+1,:));

return


function sub_test

n = 25;
r = linspace(-1,1,n);
[X,Y,Z] = meshgrid(r,r,r);
R = sqrt([X.^2 + Y.^2 + Z.^2]);
mask = R < 1;
N = n^3;

f = 0*X;

lmax = 1;  % degree
nr = lmax+1;
nc = nr;
ii = 1;
for l = 0:lmax
	for m = 0:l
		fm = evalbasis(X(mask), Y(mask), Z(mask), l, m);
		f(mask) = fm;
		subplot(nr,nc,ii); im(f); title(sprintf('l,m = %d,%d', l, m'));
		ii = ii + 1;
	end
end

return


