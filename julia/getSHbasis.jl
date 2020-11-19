using SphericalHarmonics
using CoordinateTransformations: SphericalFromCartesian

function getSHbasis(
	x::Vector{<:Real}, 
	y::Vector{<:Real}, 
	z::Vector{<:Real}, 
	L::Int64
	)

% Spherical coordinate system as defined in https://en.wikipedia.org/wiki/Spherical_harmonics
% and Romeo and Hoult MRM 1984
[ph,el,r] = cart2sph(X,Y,Z);  % ph = azimuth, el = elevation, r = radius
th = pi/2 - el; % polar angle

	c2s = SphericalFromCartesian()
	function sh(x, y, z, l, m)
		a = c2s([x,y,z])
		(r, ϕ, θ) = (a.r, a.θ, π/2 - a.ϕ)    # (radius, aziumuth, colatitude)
		Y = computeYlm(θ, ϕ; lmax=l)
		Ylm[(l,m)]
	end

% construct basis matrix H
H = zeros(size(X,1), sum(2*(0:l)+1));
ic = 1;
for l1 = 0:l
	lp = legendre(l1, cos(th));   % [l1+1 N]
	for m = 0:l1
		f = r.^l1 .* exp(1i*m*ph) .* lp(m+1,:)';
		f = 
		H(:,ic) = real(f);
		ic = ic+1;
		if m > 0
			H(:,ic) = imag(f);
			ic = ic+1;
 		end
 	end
end

