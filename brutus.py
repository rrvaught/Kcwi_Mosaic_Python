
def calc_offsets(aimg, ahdr, img, hdr, box_length):
    """
    :param aimg: anchor image
    :param ahdr: anchor header
    :param img: image to align
    :param hdr: header for image to align
    box_length: size of box in arcsecs
    :return:
    """
    from reproject import reproject_interp
    from astropy.utils.data import get_pkg_data_filename
    from astropy.wcs import WCS
    import astropy.io.fits as fits
    import numpy as np
    from astropy.wcs.utils import proj_plane_pixel_scales
    from scipy.ndimage.interpolation import shift
    from scipy.ndimage.interpolation import rotate

    def my_callback(params):
        print("Trying parameters " + str(params))
        return

    def trans_image(array, x_y_trans):
        from scipy.ndimage.interpolation import shift
        x_y_trans = np.array(x_y_trans)
        trans_slice = shift(array, [x_y_trans[1], x_y_trans[0]], cval=-10,
                        order=1)
        rot_slice = rotate(trans_slice, x_y_trans[2], reshape=False)
        return rot_slice

    def correl_mismatch(slice0, slice1):
        #     """ Negative correlation between the two images, flattened to 1D """
        w = np.where(slice0 > -9)
        slice0 = slice0[w] - np.median(slice0[w])
        slice1 = slice1[w] - np.median(slice1[w])
        return -np.corrcoef(slice0.ravel(), slice1.ravel())[0, 1]

    def fancy_cost_at_xy(x_y_trans):
        unshifted = trans_image(array, x_y_trans)
        return correl_mismatch(unshifted, cim)

    def rot_image(unshifted, rotations):
        rotations = np.array(rotations)
        rotations = np.squeeze(rotations)
        rot_slice = rotate(unshifted, rotations, reshape=False)
        return rot_slice

    awcs, wcs= WCS(ahdr), WCS(hdr)

    array, footprint = reproject_interp((img,hdr), awcs)

    w_nans =np.isnan(array)
    array[w_nans] = -10

    w_nans = np.isnan(aimg)
    aimg[w_nans] = -10

    # Initilaize grid search parameters #
    degperpix = abs(proj_plane_pixel_scales(awcs)[0])
    arcperpix = degperpix * 3600

    attempt=0

    while attempt < 2:
        if attempt==0:
            npix = round(box_length / arcperpix)
            stepr = 0.5
            step = 2.0
            gx = -11.0
            gy = 3.0
            tryx = np.arange(-npix + gx, gx + npix + 1, step)
            tryy = np.arange(-npix + gy, gy + npix + 1, step)
            tryr = np.arange(-1, 1 + stepr, stepr)

            match = 0.0
            bestx, besty, bestr, test = 0, 0, 0, 0

            i = 0

            import itertools
            for x in itertools.product(tryx, tryy, tryr):
                im = rotate(shift(array, [x[1], x[0]], order=0, cval=-10), x[2], reshape=False)
                test = abs(correl_mismatch(im, aimg))

                if np.isfinite(test) == False:
                    test = 0
                while i == 0:
                    print('Initial Correlation = ', test)
                    match = test
                    i = 1
                if test > match:
                    bestx = x[0]
                    besty = x[1]
                    bestr = x[2]
                    print('First pass parameters updated to:')
                    print(' bestx =',bestx,' besty =', besty, ' rot =',bestr, ' Corr =',abs(test))
                    match = abs(test)
                else:
                    continue
            attempt += 1

        else:
            npix = round(0.4/ arcperpix)
            stepr = 0.5
            step = 1.0
            gx = bestx
            gy = besty
            tryx = np.arange(-npix + gx, gx + npix + 1, step)
            tryy = np.arange(-npix + gy, gy + npix + 1, step)
            tryr = np.arange(-1, 1 + stepr, stepr)

            bestx, besty, bestr, test = 0, 0, 0, match

            i = 0

            import itertools
            for x in itertools.product(tryx, tryy, tryr):
                im = rotate(shift(array, [x[1], x[0]], order=0, cval=-10), x[2], reshape=False)
                test = abs(correl_mismatch(im, aimg))

                if np.isfinite(test) == False:
                    test = 0
                while i == 0:
                    print('Initial Correlation = ', test)
                    match = test
                    i = 1
                if test > match:
                    bestx = x[0]
                    besty = x[1]
                    bestr = x[2]
                    print('Updating best shift parameters to:')
                    print(' bestx =', bestx, ' besty =', besty, ' rot =', bestr, ' Corr =', abs(test))
                    match = abs(test)
                else:
                    continue
            attempt += 1

    shift_x = degperpix * bestx
    shift_y = degperpix * besty
    print('Using shifts')
    print(bestx, besty, bestr, abs(match))
    return  shift_x, shift_y, bestr