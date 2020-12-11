using MAT, JLD2

matf = matread("CalibrationDataUM10Dec2020.mat");
F = matf["F"];
S = matf["S"];
fov = matf["fov"];
desc = matf["desc"];

@save "CalibrationDataUM10Dec2020.jld2" F S fov desc
