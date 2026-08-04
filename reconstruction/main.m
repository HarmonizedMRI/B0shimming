
% TE difference between the two images
sys = mr.opts('B0', 3.0);
fatChemShift = 3.5e-6;          % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz
TE = 1/fatOffresFreq*[1 2];     % fat and water in phase for both echoes
n_TE = length(TE);                 
assert(n_TE == 2, 'acquisition must have 2 echoes');
delta_TE = diff(TE);        % sec 

% acquisition matrix size
n_x = 100; n_y = n_x; n_z = n_x;
n_z_dummy = 1;

FOV = [24 24 24]*1e-2;     % m

% load GE ScanArchive file and check size
%d_cell = pge2.utils.loaddata('data.h5');   % cell array
n_shots = length(d_cell);
assert(n_shots == n_TE * n_y * (n_z_dummy + n_z), ...
    'The number of TRs in data file does not match stated acquisition parameters');
assert(n_x == size(d_cell{1}, 1), 'b0.N(1) does not match the data file');

n_c = size(d_cell{1}, 2);  % number of receive coils

% sort data into an array
d = zeros(n_x, n_c, n_TE*n_y, n_z);

shot = n_z_dummy * n_y + 1;   % first shot after dummy shots
for z = 1 : n_z
    for y = 1 : n_TE*n_y
        d(:, :, y, z) = d_cell{shot};
        shot = shot + 1;
    end
end

d = permute(d, [1 3 4 2]);

