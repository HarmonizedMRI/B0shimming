load exampledata    % A f0 X Y Z mask

H = shim.getSHbasis(X(mask),Y(mask),Z(mask),l);   % [N 2l+1]
f0 = f0(mask);
W = diag_sp(ones(size(f0,1),1));
shat = -(W*H*A)\(W*f0);    % [9 1]. NB! May need to be rounded before applying settings on scanner.

f = f0 + H*A*shat;  % predicted fieldmap after applying shims

% put back into 3D matrix
f0 = embed(f0, mask);
f = embed(f, mask);

% compare baseline and predicted fieldmaps
subplot(121)
im(f0); 
title(sprintf('before shimming'));
h = colorbar; h.TickLabels{end} = 'Hz'; % h.Label.String = 'B0 field (Hz)';
subplot(122)
im(f,20*[-1 1]); 
title(sprintf('Residual'));
h = colorbar; h.TickLabels{end} = 'Hz';

