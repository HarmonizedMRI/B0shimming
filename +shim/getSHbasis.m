function H = getSHbasis(X,Y,Z,l)
% function H = getSHbasis(X,Y,Z,l)
%
% Evaluate spherical harmonic basis up to degree 'l'
%
% Inputs:
%  X/Y/Z    [N 1]     x/y/z coordinates at which to evaluate (discretize) the basis (cm)
%  l        int       SH degree. Default: 4.
%
% Output:
%   H       [N ...]    SH basis, including dc (B0 offset) term

if ~isvector(X) | ~isvector(Y) | ~isvector(Z)
	error('X, Y, Z must be vectors');
end
N = numel(X);
if numel(Y) ~= N | numel(Z) ~= N
	error('X, Y, Z must be the same length');
end
if nargin < 4
	l = 4;
end
if rem(l, 1) | l < 0
	error('Order must be non-negative integer');
end

% Spherical coordinate system as defined in https://en.wikipedia.org/wiki/Spherical_harmonics
% and Romeo and Hoult MRM 1984
[ph,el,r] = cart2sph(X,Y,Z);  % ph = azimuth, el = elevation, r = radius
th = pi/2 - el; % polar angle

% construct basis matrix H
H = zeros(size(X,1), sum(2*(0:l)+1));
ic = 1;
for l1 = 0:l
	lp = legendre(l1, cos(th));   % [l1+1 N]
	for m = 0:l1
		f = r.^l1 .* exp(1i*m*ph) .* lp(m+1,:)';
		H(:,ic) = real(f);
		ic = ic+1;
		if m > 0
			H(:,ic) = imag(f);
			ic = ic+1;
 		end
 	end
end

% Equivalent Cartesian expressions for l=2
%x = X(:); y = Y(:); z = Z(:);
%H = [ones(length(x),1) x y z z.^2-1/2*(x.^2+y.^2) x.*y z.*x x.^2-y.^2 z.*y]; 


