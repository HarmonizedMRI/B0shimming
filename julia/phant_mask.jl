
using MIRTjim: jim

function phant_mask(dat,thresh)

    img1f = dat
    thresh = thresh*maximum(img1f)
    threshHigh = maximum(img1f) ### Adjust for different images
    new_img = zeros(size(img1f))
for i in 1:size(img1f,1)
    for j in 1:size(img1f,2)
        for k in 1:size(img1f,3)
            if img1f[i,j,k,1] > thresh
                if img1f[i,j,k,1] < threshHigh
                new_img[i,j,k] = img1f[i,j,k,1]
                end
            end
        end
    end
end
mask = iszero.(iszero.(new_img))



return(mask)

end