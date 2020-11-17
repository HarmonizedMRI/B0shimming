
using LinearAlgebra: norm, opnorm
using Flux
using MIRT: embed!
using Random: seed!
using LaTeXStrings
using Plots; default(markerstrokecolor=:auto)

"""
	function getSHbasis(x,y,z,l)

Evaluate spherical harmonic basis functions up to degree 'l' at locations [x y z]

Inputs:
  x/y/z    length-N vector     x/y/z coordinates at which to evaluate (discretize) the basis (cm)
  l        int                 SH degree. Default: 4.

Output:
  H       [N ...]    SH basis, including dc (B0 offset) term
"""
function getSHbasis(
	x::Vector{<:Real},
	y::Vector{<:Real},
	z::Vector{<:Real},
	l::Int = 4
)

% Spherical coordinate system as defined in https://en.wikipedia.org/wiki/Spherical_harmonics
% and Romeo and Hoult MRM 1984
[ph,el,r] = cart2sph(x,y,z);  % ph = azimuth, el = elevation, r = radius
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

