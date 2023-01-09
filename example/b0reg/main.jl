using MRIFieldmaps

# Load complex coil images 
matf = matread("../images.mat")
images = matf["images"]   # [nx ny ... ncoil nechotime]
images = ComplexF32.(images) # 32-bit floats saves memory and thus time

# spatial mask
matf = matread("../mask.mat")
mask = BitArray(matf["mask"]) 

# echo times (sec)
matf = matread("../echotime.mat")
tmp = matf["echotime"]
echotime = [tmp[1], tmp[2]]  # vector

# coil sensitivity maps (see getsens.m)
matf = matread("../sens.mat")
sens = matf["sens"] 
sens = ComplexF32.(sens) # 32-bit floats saves memory and thus time

# initial (unwrapped) field map
matf = matread("../b0init.mat")
b0init = matf["b0init"]
b0init = ComplexF32.(b0init) # 32-bit floats saves memory and thus time

# get regularized fieldmap in Hz
b0, _, _ = b0map(images, echotime; smap=sens, finit=b0init) # , mask=mask) 

# Write to .mat file for viewing
matwrite("b0.mat", Dict(
   "b0" => b0
))

