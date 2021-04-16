using MAT, JLD2

# Read .mat file and write to .jld2 file

matf = matread("f0_jar.mat");
f0 = matf["f0"];          # 3D fieldmap we wish to shim (Hz)
fov = matf["fov"];
mask = matf["mask"];
@save "f0_jar.jld2" f0 fov mask

