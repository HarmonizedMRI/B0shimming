
include("phant_mask.jl")

function mask_THR(img, TH)
    global mask_flag=0
    global mask=[]
    TH2=TH
    while mask_flag == 0
        global mask_flag, mask
        mask = phant_mask(img,TH2)
        pl = jim(mask, title = "Mask")
        display(pl)
        println("Is sufficient? [y/n]")
        mask_good = readline();
        if mask_good == "y"
            global mask_flag = 1
        else
            println("Enter threshold ratio? [0:1]")
            println("Current Threshold:", TH2)
            TH2 = parse(Float32,readline());
        end
    end
    mask = mask[:,:,:,1]
    return mask


end