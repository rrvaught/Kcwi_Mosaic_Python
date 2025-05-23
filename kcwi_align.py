import numpy as np

from kcwi_repro import repro


def align_sf(list_of_images,directory):
    from kcwi_helper_functions import make_head_2d
    from kcwi_helper_functions import WcsUpdate
    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.io.fits as fits
    from brutus import calc_offsets

    anchor = list_of_images[0]
    anchdata = np.nansum(fits.open(anchor)[0].data[1000:,:,:], axis=0)
    plt.imshow(anchdata)
    plt.show()
    anchhdr = make_head_2d(fits.open(anchor)[0].header)
    for image in list_of_images[1:]:
        data = np.nansum(fits.open(image)[0].data[1000:,:,:], axis=0)
        datahdr = make_head_2d(fits.open(image)[0].header)
        sx, sy, sr= calc_offsets(anchdata, anchhdr, data, datahdr, 1.5)
        WcsUpdate(image, sx, sy, sr, directory)
    return

def align_ri(list_of_images,ref_image,filter, directory):
    from astropy.nddata.utils import Cutout2D
    from astropy.wcs import WCS
    from scipy import interpolate
    from astropy import units as u
    from astropy.wcs.utils import skycoord_to_pixel
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import ICRS, Galactic, FK4, FK5
    from kcwi_trim import make_wave
    from kcwi_helper_functions import make_head_2d
    from kcwi_helper_functions import WcsUpdate
    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.io.fits as fits
    from brutus import calc_offsets
    from reproject.mosaicking import find_optimal_celestial_wcs
    from reproject import reproject_interp

    # ensure ref image is north aligned
    wcs_out, shape_out = find_optimal_celestial_wcs((fits.open(ref_image)[0].data,fits.open(ref_image)[0].header))
    anchor, footprint = reproject_interp((fits.open(ref_image)[0].data,fits.open(ref_image)[0].header),
                                           wcs_out, shape_out=shape_out)

    anchorhdr = wcs_out.to_header()

    filter_name = '/Users/rrickardsvaught/PycharmProjects/Kcwi_Mosaic_Python/Data/transmission_curves/'+filter
    filt = np.genfromtxt(filter_name, dtype=[('wave', float), ('flux', float)])
    x = filt['wave']
    xliml, xlimh=x[0], x[-1]
    y = filt['flux']
    f = interpolate.interp1d(x, y)

    for image in list_of_images:
        # Read in image #
        data = fits.open(image)[0].data
        datahdr = make_head_2d(fits.open(image)[0].header)

        wave = make_wave(fits.open(image)[0].header)
        RA, DEC = datahdr['CRVAL1'], datahdr['CRVAL2']

        # convolve image with transmission curve #
        wl, yy, xx = np.shape(data)
        filt_image = np.zeros((yy, xx))
        wgood = np.logical_and(wave >  xliml, wave <  xlimh)

        it = np.nditer(filt_image, flags=['multi_index'])
        while not it.finished:
            i, j = it.multi_index[0], it.multi_index[1]

            cvolve = f(wave[wgood]) * data[wgood, i, j]
            filt_image[i, j] = np.nansum(cvolve) / np.nansum(f(wave[wgood]))
            it.iternext()

        wcs = WCS(anchorhdr)

        # Create SDSS cutout
        c = SkyCoord(ra=RA * u.deg, dec=DEC * u.deg, frame=FK5)
        (x, y) = c.to_pixel(wcs, origin=0, mode='all')
        xx, yy = int(x), int(y)
        cutout = Cutout2D(anchor, (xx, yy), (800, 800), wcs=wcs)

        # Write out coutout.wcs to new header
        n_anchhead = cutout.wcs.to_header()
        n_anchhead['NAXIS'] = 2
        n_anchhead['NAXIS1'] = np.shape(cutout.data)[1]
        n_anchhead['NAXIS2'] = np.shape(cutout.data)[0]
        print(n_anchhead)

        sx, sy, sr = calc_offsets(cutout.data, n_anchhead, filt_image, datahdr,8)
        WcsUpdate(image, sx, sy, sr, directory)

    return


def align(align_type,directory,ref_image=None, filter=None):
    import glob
    list_of_images = glob.glob(directory + '/working_directory/scubes/projected/*.fits')
    if align_type=='single-field':
        # align by first image #
        align_sf(list_of_images,directory)

    if align_type=='reference-image':
        align_ri(list_of_images,ref_image,filter, directory)


    return

