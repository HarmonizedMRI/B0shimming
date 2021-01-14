
load result_p2_l6.mat
gp_p2_l6 = gp .* mask;
fp_p2_l6 = fp .* mask;

load result_p10_l6.mat
gp_p10_l6 = gp .* mask;
fp_p10_l6 = fp .* mask;

g0 = g0 .* mask;

im(cat(1,g0, gp_p2_l6, gp_p10_l6), [-80 80]); colormap jet;

