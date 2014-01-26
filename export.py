from glob import glob

import numpy as np
import h5py
import pyfits
import healpy as hp

from planck.LFI import LFI
from testenv.todtools import check_bit
from planck.pointing import DiskPointing

lfi = LFI()

base_folder = "/project/projectdirs/planck/data/mission/lfi_ops_dx10/"
freq = 70

for od in range(91, 1604+1):
    matching_files = glob(base_folder + "%04d/*%03d*-R-*.fits" % (od, freq)) 
    assert len(matching_files) == 1, "More than one file for OD %d, freq %d" % (od, freq)
    fits_file = pyfits.open(matching_files[0])
    pointing_fits_file = pyfits.open(matching_files[0])
    pnt = DiskPointing(od, freq)

    fields = ['od', 'ring', 'glon', 'glat', 'psi', 'healpix_2048', 'tsky', 'dipole', 'utc']

    from dipole import Dipole
    dip = Dipole(obt=fits_file["OBT"].data["OBT"]/2**16, type="total") 

    output_folder = "out/"
    output_filename = "%d_%04d.h5" % (freq, od)

    with h5py.File(output_folder + output_filename, 'w') as output_h5:
        for ch in lfi.f[freq].ch[:1]:
            print ch
            good_data = np.logical_or(~check_bit(fits_file["OBT"].data["FLAG"], 0), fits_file[ch.tag].data["FLAG"]==0) 
            n_good_data = good_data.sum()

            data=np.zeros(n_good_data,dtype={
                    'names':fields,
                    'formats':['int64','int64','float64','float64','float64','int64','float64','float64','int64']
                                            })

            # timing
            # TODO use UTC instead of OBT
            data["od"] = od
            data["utc"] = fits_file["OBT"].data["OBT"][good_data]/2**16
            data["ring"] = 0

            # pointing
            theta, phi, psi = pnt.get_3ang(ch)
            data['glon'] =      np.degrees(phi[good_data])
            data['glat'] = 90 - np.degrees(theta[good_data])
            data['psi']  =      np.degrees(psi[good_data])
            data['healpix_2048'] = hp.ang2pix(2048, theta[good_data], phi[good_data])

            # dipole
            data['dipole'] = dip.get_4piconv_dx10(ch, theta, phi, psi)[good_data]

            # tsky (LFI data are already dipoleremoved
            # TODO load and remove Madam baselines
            data['tsky'] = fits_file[ch.tag].data[ch.tag][good_data]
            output_h5.create_dataset(ch.tag, data=data, compression='lzf')
