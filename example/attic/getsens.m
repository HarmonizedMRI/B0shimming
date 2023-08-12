% Figure out correct image orientation for sens map obtained with BART
% We will do this by visually comparing against a rough sens map 
% obtained using root-sum-of-squares coil-combined image as reference
if false
    pfile = [datDir 'P,b0.7'];
    readoutFile = [datDir 'readout_b0.mod'];
    [ims, imsos] = toppe.utils.recon3dft(pfile, ...
        'echo', 1, ...
        'readoutFile', readoutFile, ...
        'alignWithUCS', true);  

    mask = imsos > max(0.1*imsos(:));

    sens_rough = bsxfun(@rdivide, ims, imsos);
end

% load sensitivity maps. Permute so it aligns with UCS.
load([datDir 'sens_bart']);
sens = sens_bart(:,:,:,:,1);
sens = flipdim(sens, 1);
sens = flipdim(sens, 2);
sens = flipdim(sens, 3);
save sens sens
