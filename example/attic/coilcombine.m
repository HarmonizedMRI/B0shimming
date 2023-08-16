function x = coilcombine(xin, sens, mask)
% function imout = coilcombine(imsin, sens, mask)
%
% Divide each coil image by its sensitivity map and add across coils
%
% imsin  [nx ny (nz) nCoil] complex coil images
% sens   [nx ny (nz) nCoil] sensitivity map

% object mask (across coils)
sensmask = abs(sens) ~= 0;

% divide by coil sensitivity
x = xin./sens;
x(~sensmask) = 0;

% sum over coils and apply mask
x = sum(x, ndims(xin));
x(~mask) = 0;

% remove outliers
mx = mean(abs(x(mask)));
x(abs(x)>5*mx) = 0;

