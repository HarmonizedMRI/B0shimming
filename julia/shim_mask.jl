include("mask_THR.jl")
function shim_mask(img1, mtype,TH=0.4)
    mask=[]
    maskbs=[]
    maskbrain=[]
    mask_w=[]

    if mtype == 1 # Threshold mask
        
        mask=mask_THR(img1,TH)
    end


    if mtype == 2 # BET mask
        niii = NIVolume(img1)

        niwrite("image1_for_mask.nii", niii)
        
        cmd = `bet image1_for_mask.nii bet_term.nii -m -n -f 0.85`   # 0.65 for phantom   0.75 for brain
        unz = `gunzip bet_term_mask.nii.gz`
        run(cmd)
        run(unz)
        
        mask1 = niread("bet_term_mask.nii")
        mask1 = mask1.raw
        mask1 = iszero.(iszero.(mask1))
        maskTh = mask_THR(img1,TH)
        mask = mask1.*maskTh
        
    end



    if mtype == 3 # Cylindircal mask
        mask_flag=0
        pos = [30,30,30]
        global mask1 = zeros(size(img1[:,:,:,1]))
        while mask_flag == 0
            global TH, mask_flag, mask1
            global mask1 = bs_mask(img1,pos)
            mask1p,img1p = fmap_pos(mask1[:,:,:,1],mask1[:,:,:,1],img1)[2:3]
            pl1 = jim(cat(cat(4*mask1[:,:,10,1]+img1[:,:,10,1],4*mask1[:,:,15,1]+img1[:,:,15,1], 4*mask1[:,:,20,1]+img1[:,:,20,1];dims=1),cat(4*mask1p[:,:,27]+img1p[:,:,27],4*mask1p[:,:,31]+img1p[:,:,31], 4*mask1p[:,:,34]+img1p[:,:,34];dims=1);dims=2), title = "Mask")
            display(pl1)
            println("Is mask sufficient? [y/n]")
            mask_good = readline();
            if mask_good == "y"
                global mask_flag =1
            else
                println("Enter mask location? [x,y,z]")
                println("Current location:", pos)
                pos[1] = parse(Float32,readline());
                pos[2] = parse(Float32,readline());
                pos[3] = parse(Float32,readline());
            end
        end
        maskTh = mask_THR(img1,TH)
        mask = mask1.*maskTh
        mask = mask[:,:,:,1]
    end



    if mtype == 4  # Weighted cylindircal mask (combo 2 and 3)
        niii = NIVolume(img1)
        wt=17
        niwrite("image1_for_mask.nii", niii)
        
        cmd = `bet image1_for_mask.nii bet_term.nii -m -n -f 0.75`   # 0.65 for phantom   0.75 for brain
        unz = `gunzip bet_term_mask.nii.gz`
        run(cmd)
        run(unz)
        
        mask1 = niread("bet_term_mask.nii")
        mask1 = mask1.raw
        mask1 = iszero.(iszero.(mask1))
        maskTh = mask_THR(img1,TH)
        maskbrain = mask1.*maskTh
        pos = [30,30,30]
        global mask_flag=0
        while mask_flag == 0
            global TH, mask_flag, mask1
            mask1 = bs_mask(img1,pos)
            mask1p,img1p = fmap_pos(mask1[:,:,:,1],mask1[:,:,:,1],img1)[2:3]
            pl1 = jim(cat(cat(4*mask1[:,:,10,1]+img1[:,:,10,1],4*mask1[:,:,15,1]+img1[:,:,15,1], 4*mask1[:,:,20,1]+img1[:,:,20,1];dims=1),cat(4*mask1p[:,:,27]+img1p[:,:,27],4*mask1p[:,:,31]+img1p[:,:,31], 4*mask1p[:,:,34]+img1p[:,:,34];dims=1);dims=2), title = "Mask")
            display(pl1)
            println("Is mask sufficient? [y/n]")
            mask_good = readline();
            if mask_good == "y"
                global mask_flag =1
            else
                println("Enter mask location? [x,y,z]")
                println("Current location:", pos)
                pos[1] = parse(Float32,readline());
                pos[2] = parse(Float32,readline());
                pos[3] = parse(Float32,readline());
            end
        end
        maskTh1 = mask_THR(img1,TH)
        maskbs = mask1.*maskTh1
        mask_w = wt*maskbs+maskbrain
        mask_w[mask_w .> wt] .= wt
        mask_w = mask_w[:,:,:,1]/wt
        mask = maskbs+maskbrain
        mask = iszero.(iszero.(mask))
        mask = mask[:,:,:,1]
    end


    if mtype == 5  # Weighted saved bs mask
        reffile = matread("bsmask_matched1006.mat")
        refmask = reffile["bsmaskref"]
        refimg = reffile["imgref"]
        
        niii = NIVolume(img1)
        wt=17
        niwrite("image1_for_mask.nii", niii)
        
        cmd = `bet image1_for_mask.nii bet_term.nii -m -n -f 0.75`   # 0.65 for phantom   0.75 for brain
        unz = `gunzip bet_term_mask.nii.gz`
        run(cmd)
        run(unz)
        
        mask1 = niread("bet_term_mask.nii")
        mask1 = mask1.raw
        mask1 = iszero.(iszero.(mask1))
        maskTh = mask_THR(img1,TH)
        maskbrain = mask1.*maskTh
        pos = [30,30,30]


        global mask_flag=0
        while mask_flag == 0
            global TH, mask_flag, mask1
            mask1 = circshift(refmask,(30-pos[1],30-pos[2],30-pos[3]))
            mask1p,img1p = fmap_pos(mask1[:,:,:,1],mask1[:,:,:,1],img1)[2:3]
            pl1 = jim(cat(cat(4*mask1[:,:,10,1]+img1[:,:,10,1],4*mask1[:,:,15,1]+img1[:,:,15,1], 4*mask1[:,:,20,1]+img1[:,:,20,1];dims=1),cat(4*mask1p[:,:,27]+img1p[:,:,27],4*mask1p[:,:,31]+img1p[:,:,31], 4*mask1p[:,:,34]+img1p[:,:,34];dims=1);dims=2), title = "Mask")
            display(pl1)
            println("Is mask sufficient? [y/n]")
            mask_good = readline();
            if mask_good == "y"
                global mask_flag =1
            else
                println("Enter mask new location? [R/L,A/P,S/I]")
                println("Current Location:", pos)
                pos[1] = parse(Float32,readline());
                pos[2] = parse(Float32,readline());
                pos[3] = parse(Float32,readline());
            end
        end
        maskTh1 = mask_THR(img1,TH)
        maskbs = mask1.*maskTh1
        mask_w = wt*maskbs+maskbrain
        mask_w[mask_w .> wt] .= wt
        mask_w = mask_w[:,:,:,1]/wt
        mask = maskbs+maskbrain
        mask = iszero.(iszero.(mask))
        mask = mask[:,:,:,1]



    end

    return mask,mask_w,maskbs,maskbrain

end