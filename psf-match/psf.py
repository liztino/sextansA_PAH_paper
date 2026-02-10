import webbpsf
import numpy as np
import matplotlib.pyplot as plt
import glob

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate

savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/psfs/'
nircam_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/nircam/'
miri_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/miri/'

m = webbpsf.MIRI()
nc = webbpsf.NIRCam()

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

filts = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W' ]
# nircam_filts = ['F300M', 'F335M', 'F360M']
filts = ['F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W']

# do just the SW nircam filts
filts = ['F115W', 'F150W', 'F200W']

# for nc_filt in nc_filts:
#     savename = f'{nc_filt}_webbpsfv121.fits'
#     nc.filter = nc_filt
#     nc.calc_psf(outfile = savedir + savename, oversample=4,  fov_pixels = 1024)
    
# # size Aniano+ 2011 uses to optimize FFT algorithms
# sample = 3
# psf_size = 3645
# pixels = 3645/3

# changing the optimized values for the SW nircam 
sample = 3
psf_size = 1125
pixels = psf_size/sample
    
for filt in filts:
    savename = f'{filt}_webbpsfv121_norot.fits'
    
    if filt in miri_filts:
        infile = glob.glob(miri_dir + f'*{filt}_skysub_i2d.fits')[0]
        
    elif filt in nircam_filts:
        infile = glob.glob(nircam_dir + '*{:s}_*i2d.fits'.format(filt.lower()))[0]

    sim = webbpsf.setup_sim_to_match_file(infile)
    sim.options['parity'] = 'odd'
    
    sim.calc_psf(outfile = savedir + savename, fov_pixels = pixels, oversample=sample)
    
    # load in header of image
    head = fits.open(infile)['SCI'].header
    
    orig_psf_hdu = fits.open(savedir + savename)[0]
    orig_psf_data = orig_psf_hdu.data
    orig_psf_header = orig_psf_hdu.header

    # calculate angle to rotate, will depend on the initial rotation of PSF
    ang = head['PA_APER']
    print(ang)

    # rotate with scipy tools
    rot_psf = rotate(orig_psf_data, -1 * ang, reshape = False, prefilter = False)

    # save the new rotated PSF
    orig_psf_header['PA_APER'] = ang

    rotname = f'{filt}_webbpsfv130_rot.fits'
    fits.writeto(savedir + rotname, rot_psf, orig_psf_header, overwrite = True )
    

# # test PSF rotation before applying it to all PSFs
# miri = webbpsf.MIRI()
# miri.filter =  'F770W'
# miri.calc_psf(savedir + "F770W_test.fits")         # you can also write the output directly to disk if you prefer.

# orig_psf_hdu = fits.open(savedir + "F770W_test.fits")[0]
# orig_psf_data = orig_psf_hdu.data
# orig_psf_header = orig_psf_hdu.header

# ang1 = orig_psf_header['DET_V3']
# ang2 = 306.8926240348045

# delta_ang = ang2 - ang1

# rot_psf = rotate(orig_psf_data, -1 * delta_ang, reshape = False, prefilter = False)

# orig_psf_header['PA_APER'] = ang2

# # fits.writeto(savedir + "F770W_test_rot.fits", rot_psf, orig_psf_header, overwrite = True )
 
