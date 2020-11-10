# Get example data generated with createexampledata.m

using JLD, MAT

matf = matread("exampledata.mat");
A = matf["A"];
f0 = matf["f0"];
fov = vec(matf["fov"]);
mask = matf["mask"];

f0masked = f0[mask];
