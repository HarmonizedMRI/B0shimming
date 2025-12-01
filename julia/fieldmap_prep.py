import numpy as np
#import matplotlib.pyplot as plt
#import nibabel as nib
import scipy.io as sio
import sys
#import loaddata as ld
from mapvbvd import mapVBVD as mapVBVD

if len(sys.argv) > 1:
    print(f"The script name is: {sys.argv[0]}")
    print("File Path: ",sys.argv[1])
else:
    print("No arguments provided.")

def loaddata_siemens(data_path):
    """
    Load data-file obtained with the b0.seq.

    Parameters:
    data_path (str): .dat-file name

    Returns:
    np.ndarray: Complex coil data for the 2 echoes [nx, ny, nz, ncoils, 2]
    """

    # Load data from .dat-file using pyMapVBVD
    twixObj = mapVBVD(data_path,quiet=True)
    twixObj = twixObj[1]
    twixObj.image.squeeze = True
    data_unsorted = twixObj.image['']
    twixObj.image.flagRemoveOS = False
    data_unsorted = twixObj.image.unsorted()

    # Rearrange the dimensions [nfid, nview, nslice, ncoil]
    #print('data_unsorted.shape',data_unsorted.shape)
    din = np.transpose(data_unsorted, (0, 2, 1))  # [nfid, nview, nslice, ncoil]

    # Remove dummy shots
    nz_dummy = 1
    nx = din.shape[0]
    ny = nx  # twix.hdr.Meas.ImageColumns, twix.hdr.Meas.ImageLines
    nz = nx  # Adjust accordingly
    n_coils = din.shape[2]
    #print('n_coils',n_coils)

    n_te = 2  # Number of echo times (interleaved)
    din = din[:, (nz_dummy * ny * n_te):, :]  # [nx, ny * nz * n_te, ncoils]

    d = np.zeros((nx, ny * nz, n_coils, n_te), dtype=np.complex_)
    for i_te in range(n_te):
        d[:, :, :, i_te] = din[:, i_te::n_te, :]

    d = d.reshape((nx, ny, nz, n_coils, n_te))
    d = np.transpose(d, (0, 2, 1, 3, 4))  #to stay conistent with matlab implementation 

    return d


# Echo times
echo_times = [2.22, 4.45]  # NW this should be read from the data or the pulseq file

# define Data Path

data_file_path = sys.argv[1]

# Load the data

data = loaddata_siemens(data_file_path)

# Number of coils
n_coil = data.shape[3]
im_t1 = data[..., 0]
im_t2 = data[..., 1]

# Perform FFT operations
im_te1 = np.zeros_like(im_t1, dtype=np.complex_)
im_te2 = np.zeros_like(im_t2, dtype=np.complex_)

for ic in range(n_coil):
    im_te1[..., ic] = np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(im_t1[..., ic])))
    im_te2[..., ic] = np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(im_t2[..., ic])))

# Compute the images
img1 = np.sqrt(np.sum(np.abs(im_te1), axis=3))
img2 = np.sqrt(np.sum(np.abs(im_te2), axis=3))

# Show the images
#print('img1.shape', img1.shape)
#plt.figure()
#plt.imshow(img1[:, :, 30])
#plt.title('Sliced Data')
#plt.show()

# Save the NIfTI file
#nib.save(nib.Nifti1Image(img1, np.eye(4)), 'sball.nii')

# Save the .mat file
#sio.savemat(data_path + 'sball.mat', {'im_te1': im_te1, 'im_te2': im_te2, 'echotimes': echo_times})
sio.savemat('/home/jpothoof/Julia_Files/sball.mat', {'im_te1': im_te1, 'im_te2': im_te2, 'echotimes': echo_times})

print('fieldmap prep done!')
