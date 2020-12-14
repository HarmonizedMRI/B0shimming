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
matf = matread("f0.mat");
f0 = matf["f0"];
fov = matf["fov"];
mask = matf["mask"];
@save "f0_redhead.jld2" f0 fov mask

