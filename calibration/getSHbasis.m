function H = getSHbasis(X,Y,Z,ord)
% function H = getSHbasis(X,Y,Z,ord)
%
% Get spherical harmonic basis of order 'ord'
%
% Inputs:
%  X/Y/Z    [N 1]    x/y/z coordinates at which to evaluate (discretize) the basis (cm)
%  ord      int      SH order (0, 1, or 2)
%
% Output:
%   H       [N 9]    SH basis, including dc (B0 offset) term

if ~isvector(X) | ~isvector(Y) | ~isvector(Z)
	error('X, Y, Z must be vectors');
end
N = numel(X);
if numel(Y) ~= N | numel(Z) ~= N
	error('X, Y, Z must be the same length');
end
if rem(ord, 1) | ord < 0
	error('Order must be non-zero integer');
end

x = X(:);
y = Y(:);
z = Z(:);

switch ord
	case 0
		H = [ones(N,1)];
	case 1
		H = [ones(N,1) x y z ];
	case 2
		H = [ones(N,1) x y z z.^2 x.*y z.*x x.^2-y.^2 z.*y];
	otherwise
		error(('Order > 2 not supported'));
end

