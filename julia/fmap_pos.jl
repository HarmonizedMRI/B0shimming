using MIRTjim: jim, prompt; jim(:prompt, false)
using Unitful: Hz, s

function fmap_pos(fmap, mask, img)

    img1 = img[:,:,:,1]

    fmap1=permutedims(fmap,(2,3,1))
    mask1=permutedims(mask,(2,3,1))
    img2=permutedims(img1,(2,3,1))

    # fmap1=permutedims(fmap,(1,3,2))
    # mask1=permutedims(mask,(1,3,2))
    # img2=permutedims(img1,(1,3,2))

    fmapp = fmap1[end:-1:1,end:-1:1,:]
    maskp = mask1[end:-1:1,end:-1:1,:]
    img3 = img2[end:-1:1,end:-1:1,:]

    return fmapp, maskp, img3
end