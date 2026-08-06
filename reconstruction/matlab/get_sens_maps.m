% see main.m

% create low-resolution copy of d
% size(d) = [n_x 2*n_y n_z n_c]
d_te1 = d(:,1:2:end,:,:);
d_lo = d_te1(end/2+1-n/2:end/2+n/2, ...
             end/2+1-n/2:end/2+n/2, ...
             end/2+1-n/2:end/2+n/2, ...
             :);

% Reshape from [48 x 48 x 48 x 32] to [48 x 48 x 48 x 1 x 32]
%d_lo_5d = reshape(d_lo, [size(d_lo,1), size(d_lo,2), size(d_lo,3), 1, size(d_lo,4)]);

fprintf('Computing sensitivity maps (this may take a while)\n');
tic; sense = bart('ecalib -r 20', d_te1); toc;
smap = bart('slice 4 0', sense); % Equivalently: smap = squeeze(sense(:,:,:,:,1));

% Run ecalib
%sens = bart('ecalib -m 1', d_lo_5d);
