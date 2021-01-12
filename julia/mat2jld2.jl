using MAT, JLD2

if false
# shim calibration data
matf = matread("CalibrationDataUM10Dec2020.mat");
matf = matread("CalibrationDataSiemensMGH12Dec2020.mat");
F = matf["F"];
S = matf["S"];
fov = matf["fov"];
mask = matf["mask"];
desc = matf["desc"];
@save "CalibrationDataSiemensMGH12Dec2020.jld2" F S fov mask desc
end

# fieldmap we wish to shim
# matf = matread("f0_redhead_localmask.mat");
# matf = matread("f0_redhead.mat");
matf = matread("f0_mgh.mat");
f0 = matf["f0"];
fov = matf["fov"];
mask = matf["mask"];
@save "f0_mgh.jld2" f0 fov mask
#@save "f0_redhead.jld2" f0 fov mask

# Psub1.mat
matf = matread("Psub1_localmask.mat")
f0 = matf["f0"];
fov = vec(matf["fov"]);
mask = matf["mask"];
@save "Psub1_localmask.jld2" f0 fov mask

# Psub1.mat
matf = matread("Psub1_z41_70.mat")
f0 = matf["f0"];
fov = vec(matf["fov"]);
mask = matf["mask"];
@save "Psub1_z41_70.jld2" f0 fov mask
