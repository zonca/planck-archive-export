from glob import glob

import numpy as np
import h5py
import pyfits
import healpy as hp

from planck.LFI import LFI
from testenv.todtools import check_bit
from planck.pointing import DiskPointing
from planck.metadata import obt2utc

class MadamBaselines(object):
    """MadamBaselines class

    Loads madam fits baselines and removes them from data"""

    def __init__(self, filename):
        self.filename = filename
        with pyfits.open(self.filename) as f:
            self.baselines_obt = np.array(f[1].data["time"]).flatten()

    def get_baselines(self, chtag):
        with pyfits.open(self.filename) as f:
            # need to remove white spaces from name
            detector_index = [cht.strip() for cht in np.array(f[2].data["detector"])].index(chtag)
            return np.array(f[2].data["baseline"][detector_index]).flatten()

    def baseline_remove(self, obt, data, chtag):
            return data - np.interp(obt, self.baselines_obt, self.get_baselines(chtag))


lfi = LFI()

base_folder = "/project/projectdirs/planck/data/mission/lfi_ops_dx11_delta/"
baselines_file ="/global/project/projectdirs/planck/data/mission/baselines/lfi/dx11_delta/base_dx11_delta_%03d_full.fits"
freq = 70

for od in range(700, 1604+1):
    print "*" * 30 + str(od) + "*" * 30
    matching_files = glob(base_folder + "%04d/?%03d-*-R-*.fits" % (od, freq)) 
    assert len(matching_files) == 1, "More than one file for OD %d, freq %d, %s" % (od, freq, str(matching_files))
    fits_file = pyfits.open(matching_files[0])
    pointing_fits_file = pyfits.open(matching_files[0])
    pnt = DiskPointing(od, freq)

    fields = ['od', 'ring', 'glon', 'glat', 'psi', 'healpix_2048', 'tsky', 'dipole', 'utc']

    from dipole import Dipole
    dip = Dipole(obt=fits_file["OBT"].data["OBT"]/2**16, type="total") 

    madam_baselines = MadamBaselines(baselines_file % freq)

    for ch in lfi.f[freq].ch[:1]:

        output_folder = "out/"
        output_filename = "%03d-%s_%04d.hdf5" % (freq, ch, od)
    
        with h5py.File(output_folder + output_filename, 'w') as output_h5:
            print ch
            good_data = np.logical_or(~check_bit(fits_file["OBT"].data["FLAG"], 0), fits_file[ch.tag].data["FLAG"]==0) 
            n_good_data = good_data.sum()
    
            data=np.zeros(n_good_data,dtype={
                    'names':fields,
                    'formats':['int64','int64','float64','float64','float64','int64','float64','float64','int64']
                                            })
    
            # timing
            data["od"] = od
            obt = fits_file["OBT"].data["OBT"][good_data]/2**16
            data["utc"] = obt2utc(obt)
            data["ring"] = 0
    
            # pointing
            theta, phi, psi = pnt.get_3ang(ch)
            data['glon'] =      np.degrees(phi[good_data])
            data['glat'] = 90 - np.degrees(theta[good_data])
            data['psi']  =      np.degrees(psi[good_data])
            data['healpix_1024'] = hp.ang2pix(1024, theta[good_data], phi[good_data])
    
            # dipole
            data['dipole'] = dip.get_4piconv_dx10(ch, theta, phi, psi)[good_data]
    
            # tsky (LFI data are already dipoleremoved
            data['tsky'] = madam_baselines.baseline_remove(obt, fits_file[ch.tag].data[ch.tag][good_data], ch.tag)
            output_h5.create_dataset("%03d-%s" % (freq, ch), data=data, compression='lzf')
    
            # attributes
            output_h5.attrs['PROCVER']  = 'DX11D'
            output_h5.attrs['TELESCOP'] = 'PLANCK'
            output_h5.attrs['INSTRUME'] = 'LFI'
            output_h5.attrs['DETNAM']   = ch
            output_h5.attrs['FREQ']     = str(freq).zfill(3)
            output_h5.attrs['OBJECT']   = 'MISSION OPERATIONAL DAY '+str(od)
            output_h5.attrs['FILENAME'] = output_filename
            output_h5.attrs['TIMEZERO'] = '1958-01-01z00:00'
