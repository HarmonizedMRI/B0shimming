system('rm b0init.mat echotime.mat images.mat f0.mat shimvol.mat mask.mat');
a = input('Delete Pfiles? (y/n) ', 's');
if strcmp(a, 'y')
    system('rm P,b0.7 P,b0,post.7');
    end
