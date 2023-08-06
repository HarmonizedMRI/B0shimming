using MAT, MIRTjim, ROMEO

# Load b0 map (radians) and magnitude image
matf = matread("th.mat")
th = matf["th"]    # wrapped 3D phase map (radians)
# mask = matf["mask"]
mag = matf["mag"]    

thuw = unwrap(th; mag=mag);

# Write to .mat file for viewing
matwrite("thuw.mat", Dict(
   "thuw" => thuw
))

