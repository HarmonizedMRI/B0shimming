using LinearAlgebra

function getSHbasis(X,Y,Z,l)

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

	return H

end
