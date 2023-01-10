% shim mask
dth = angle(x1./x2);
mask = abs(x1) > 0.05*max(abs(x1(:))) ;
mask = abs(dth)<1.0;
