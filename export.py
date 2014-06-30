import sys

def export_od(od, freq):
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

    print "*" * 30 + str(od) + "*" * 30
    matching_files = glob(base_folder + "%04d/?%03d-*-R-*.fits" % (od, freq)) 
    assert len(matching_files) == 1, "More than one file for OD %d, freq %d, %s" % (od, freq, str(matching_files))
    fits_file = pyfits.open(matching_files[0])
    pnt = DiskPointing(od, freq)

    #fields = ['od', 'ring', 'glon', 'glat', 'psi', 'healpix_1024', 'tsky', 'dipole', 'utc']
    fields = ['glon', 'glat', 'psi', 'healpix_2048', 'tsky', 'utc', 'sso']

    from dipole import Dipole
    dip = Dipole(obt=fits_file["OBT"].data["OBT"]/2**16, type="total") 

    madam_baselines = MadamBaselines(baselines_file % freq)

    for ch in lfi.f[freq].ch:

        output_folder = "out/"
        output_filename = "%03d-%s_%04d.hdf5" % (freq, ch, od)
    
        with h5py.File(output_folder + output_filename, 'w') as output_h5:
            print ch
            good_data = np.logical_not(
                check_bit(fits_file[ch.tag].data['FLAG'], 0) | \
                check_bit(fits_file[ch.tag].data['FLAG'], 2) | \
                check_bit(fits_file[ch.tag].data['FLAG'], 4) | \
                check_bit(fits_file["OBT"].data['FLAG'], 0)  | \
                check_bit(fits_file["OBT"].data['FLAG'], 2)  | \
                check_bit(fits_file["OBT"].data['FLAG'], 4)
            )
            n_good_data = good_data.sum()
    
            data=np.zeros(n_good_data,dtype={
                    'names':fields,
                    #'formats':['int64','int64','float64','float64','float64','int64','float64','float64','int64']
                    'formats':['float32','float32','float32','int32','float32','float64','u1']
                                            })
    
            print "timing"
            #data["od"] = od
            obt = fits_file["OBT"].data["OBT"][good_data]/2**16
            data["utc"] = float(obt2utc(obt))
            #data["ring"] = 0
    
            print "pointing"
            ch_pnt = ch if ch.arm == "M" else ch.pair
            theta, phi, psi = pnt.get_3ang(ch_pnt)
            data['glon'] =      np.degrees(phi[good_data])
            data['glat'] = 90 - np.degrees(theta[good_data])
            data['psi']  =      np.degrees(psi[good_data])
            if ch.arm == "S":
                # S is using the same pointing as M, so also the same psi,
                # to recover the S psi, we remove the M PSI_POL angle (~90 deg) which
                # was added by the pointing code
                data['psi'] -= np.radians(ch_pnt.get_instrument_db_field("PSI_POL"))
                data['psi'][data['psi'] < - np.pi] += 2*np.pi

            data['healpix_2048'] = hp.ang2pix(2048, theta[good_data], phi[good_data])
    
            #print "dipole"
            #data['dipole'] = dip.get_4piconv_dx10(ch, theta, phi, psi)[good_data]
    
            print "dipole-removed baseline-removed data"
            data['tsky'] = madam_baselines.baseline_remove(obt, fits_file[ch.tag].data[ch.tag][good_data], ch.tag)
            output_h5.create_dataset("%03d-%s" % (freq, ch), data=data, compression='lzf')

            # solar system object flag
            data['sso'] = check_bit(fits_file["OBT"].data['FLAG'], 3)  | \
                          check_bit(fits_file[ch.tag].data["FLAG"], 3)
    
            print "attributes"
            output_h5.attrs['PROCVER']  = 'DX11D'
            output_h5.attrs['TELESCOP'] = 'PLANCK'
            output_h5.attrs['INSTRUME'] = 'LFI'
            output_h5.attrs['DETNAM']   = ch.tag
            output_h5.attrs['FREQ']     = "%03d" % freq
            output_h5.attrs['OBJECT']   = 'MISSION OPERATIONAL DAY '+str(od)
            output_h5.attrs['FILENAME'] = output_filename
            output_h5.attrs['TIMEZERO'] = '1958-01-01z00:00'

    fits_file.close()
    print "/" * 30 + str(od) + "/" * 30

if __name__ == "__main__":
    # Run as script with od and frequency arguments
    export_od(int(sys.argv[1]), int(sys.argv[2]))
