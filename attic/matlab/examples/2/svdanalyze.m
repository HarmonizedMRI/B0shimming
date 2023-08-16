% load calibration data
load Ffull                             % [nx ny nz 8]  (does not include DC offset term)
[nx, ny, nz, nShim] = size(Ffull);
fov = 19.8*[1 1 1];   % cm

% generate mask 
mask = sum(abs(Ffull),4) > 10;         % [nx ny nz]     
N = sum(mask(:));

% mask and reshape
F = zeros(N, nShim);
for ii = 1:nShim
   tmp = Ffull(:,:,:,ii);
   F(:,ii) = tmp(mask);
end

% reduce dimensionality
F = F(1:5:end,:);

% get principal components
[U,S,V] = svd(F);
