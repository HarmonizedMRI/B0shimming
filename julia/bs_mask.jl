
using MIRTjim: jim

function bs_mask(dat,loc)
    s = size(dat)
    mask = zeros(s)

    sz = [6,6,12]

    for i = 1:2*sz[1]
        for j = 1:2*sz[2]
            for k = 1:2*sz[3]
                x = i+loc[1]-sz[1]
                y = j+loc[2]-sz[2]
                z = k+loc[3]-sz[3]
                r = sqrt((x-loc[1])^2+(y-loc[2])^2)
                if sz[1] >= r
                    if x > 0 && y > 0 && z>0
                        if x<s[1] && y<s[2] && z<s[3]
                        mask[x, y, z] = 1
                        end
                    end
                end
            end
        end
    end
    return Bool.(mask)
end