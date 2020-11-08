function getexample()

	# voxel dimensions (cm)
	Δ  = [0.3, 0.3, 0.3];

	# See ../examples/shimdemoWLS.m
	A = [1.0000e+00   1.7664e-01  -1.2696e-01  -1.2506e-01  -3.3787e+00   2.6074e-02   2.7725e-02   1.1804e-02   1.8739e-02;
   5.6691e-18   9.5727e-01   3.5927e-05   6.7842e-04  -1.2248e-03  -5.0951e-03  -5.7149e-04   3.7393e-03  -2.4849e-04;
  -7.5061e-17   1.5019e-03   9.5910e-01  -1.5123e-04   3.3408e-02   3.5609e-03   1.2344e-03   5.1306e-03  -3.8807e-04;
  -1.0248e-16   4.2560e-03  -1.6435e-04   9.4806e-01  -1.6389e-02  -3.7812e-03   4.5479e-03  -3.4316e-03  -2.6878e-03;
   6.0715e-18  -1.7399e-04  -4.9292e-05   6.6762e-05   2.5122e-02  -2.0259e-05   7.9206e-05   2.0473e-05  -6.6071e-05;
  -1.2062e-18  -2.9468e-05  -1.9254e-04  -5.9752e-05   5.7926e-05   1.9376e-02  -1.3662e-05   6.0089e-05  -2.7733e-05;
  -1.2265e-18   2.4758e-04   6.6401e-05   2.5912e-05  -2.1676e-04  -3.6114e-05   2.2779e-02   7.6101e-05   1.8713e-05;
   3.3644e-18  -1.9634e-04  -4.1191e-05   1.4665e-05  -1.8696e-05  -3.9703e-05  -2.9870e-05   9.6190e-03   2.5017e-05;
   5.6107e-18   1.8954e-04   5.3357e-05  -8.4566e-05  -3.8129e-05   1.1254e-04   1.1934e-04   2.3482e-05   1.1446e-02]

	Δs = [0,10,10,10,300,200,300.,500,500.];  # change in shim settings

	te = 30e-3; # echo time (sec)

	# voxel locations
	fov = 20;   # cm
	dx = 0.2;   # cm
	x = (-fov/2+dx/2):dx:(fov-dx/2)
	y = -x
	z = x/2
	N = length(x)
	r = map((x,y,z) -> [x, y, z], x, y, z)  # r = a vector of vectors

	# B0 gradients 
	g = Vector{Vector{Float64}}(undef,N)
	for i in 1:N
		g[i] = [0.,0.,0.];
	end

end
