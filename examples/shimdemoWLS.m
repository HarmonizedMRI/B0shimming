% demoWLS.m
%
% Demonstrate the proposed B0 shimming workflow using a synthesized baseline fieldmap f0.
% 
% Goal is to calculate shim setting changes s that minimize the weighted RMS B0 inhomogeneity.
% Provided here as an example; other loss functions may be more useful depending on the application.

% createexampledata; 

load exampledata   % A f0 X Y Z mask

% get shim amplitudes 'shat' that minimize RMS fieldmap
f0 = f0(mask);   % [N 1]
N = length(f0);
W = diag_sp(ones(N,1));
shat = -(W*H*A)\(W*f0);    % [9 1]. NB! May need to be rounded before applying settings on scanner.

f = f0 + H*A*shat;  % predicted fieldmap after applying shims

f0 = embed(f0, mask);
f = embed(f, mask);

% compare baseline and predicted fieldmaps
im(cat(1, f0, f)); 
title(sprintf('fieldmaps: original (left in each panel),\nand after 2nd order shimming (right)'));
h = colorbar; h.TickLabels{end} = 'Hz'; % h.Label.String = 'B0 field (Hz)';



