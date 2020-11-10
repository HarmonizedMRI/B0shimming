# Get example data generated with createexampledata.m

using JLD, MAT

matf = matread("exampledata.mat");
A = matf["A"];
f0 = matf["f0"];
(X,Y,Z) = (matf["X"], matf["Y"], matf["Z"])  
mask = matf["mask"];

f0masked = f0[mask];
