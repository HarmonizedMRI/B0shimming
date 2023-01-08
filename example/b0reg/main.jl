using MRIFieldmaps

# Load b0 map (radians) and magnitude image
matf = matread("images.mat")
images = matf["images"]   # [nx ny ... ncoil nechotime]
images = ComplexF32.(images) # 32-bit floats saves memory and thus time
matf = matread("echotime.mat")
tmp = matf["echotime"]
echotime = [tmp[1], tmp[2]]  # vector

b0, _, _ = b0map(images, echotime) # regularized fieldmap in Hz

# Write to .mat file for viewing
matwrite("b0.mat", Dict(
   "b0" => b0
))

