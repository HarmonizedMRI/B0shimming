import numpy as np
from mapvbvd import mapVBVD as mapVBVD

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


def loaddata_ge(data_path):
    """
    Load data-file obtained with the b0.seq.

    Parameters:
    data_path (str): .dat-file name

    Returns:
    np.ndarray: Complex coil data for the 2 echoes [nx, ny, nz, ncoils, 2]
    """
    #???

    return d
