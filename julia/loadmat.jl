using MAT

# Load shim calibration data
matf = matread("shimcal.mat");
F = matf["F"]   # [nx ny nz 8] for 2nd order shim systems
S = matf["S"]   # [8 8] matrix with shim amplitudes (typically diagonal)
maskcal = matf["mask"];
Xcal = matf["X"];
Ycal = matf["Y"];
Zcal = matf["Z"];

# Load f0 
matf = matread("f0.mat");
f0 = matf["f0"];       # 3D fieldmap we wish to shim (Hz)
mask = matf["mask"];
X = matf["X"];
Y = matf["Y"];
Z = matf["Z"];


