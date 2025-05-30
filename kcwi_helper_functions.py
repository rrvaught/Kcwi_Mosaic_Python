def make_dir(path):
    import os
    try:
        os.mkdir(path)
    except OSError:
        print('Directory {0} already created'.format(path))
    return

def make_directories(working_dir):
    make_dir(working_dir + '/working_directory/scubes/raw')
    make_dir(working_dir + '/working_directory/scubes/projected')
    make_dir(working_dir + '/working_directory/scubes/final')
    make_dir(working_dir + '/working_directory/scubes/aligned')
    make_dir(working_dir + '/working_directory/scubes/final_project')

    make_dir(working_dir + '/working_directory/vcubes/raw')
    make_dir(working_dir + '/working_directory/vcubes/projected')
    make_dir(working_dir + '/working_directory/vcubes/final')
    make_dir(working_dir + '/working_directory/vcubes/aligned')
    make_dir(working_dir + '/working_directory/vcubes/final_project')

    return

def unpack_raw_cubes(list_of_files,rawdir):
    import numpy as np
    import astropy.io.fits as fits
    import glob
    from kcwi_trim import path_leaf

    print(list_of_files)
    make_dir(rawdir+'/working_directory')

    output_dir = rawdir+'/working_directory/scubes'
    output_mask = rawdir+'/working_directory/mcubes'
    output_var = rawdir+'/working_directory/vcubes'

    make_dir(output_dir)
    make_dir(output_mask)
    make_dir(output_var)

    for name in list_of_files:
        filename = path_leaf(name)[0]
        # Open the multi extension fits out put by python DRP
        DATA = fits.open(name)

        # if science
        try:
            sname = filename[:-11] + 'icubes.fits'
            hdu = fits.PrimaryHDU(data=DATA[0].data, header=DATA[0].header)
            hdul = fits.HDUList([hdu])
            hdul.writeto(output_dir + '/' + sname, overwrite=True)
        except:
            print('Science cubes stored already')

        # if mask
        try:
            mname = filename[:-11] + 'mcubes.fits'
            hdu = fits.PrimaryHDU(data=DATA[1].data, header=DATA[0].header)
            hdul = fits.HDUList([hdu])
            hdul.writeto(output_mask + '/' + mname, overwrite=True)
        except:
            print('Mask cubes stored already')

        # if var
        try:
            vname = filename[:-11] + 'vcubes.fits'
            hdu = fits.PrimaryHDU(data=DATA[2].data * DATA[2].data, header=DATA[0].header)
            hdul = fits.HDUList([hdu])
            hdul[0].header['UTYPE'] = 'Variance'
            hdul.writeto(output_var + '/' + vname, overwrite=True)
        except:
            print('Variance cubes stored already')
    return

def clean_overlap(directory):
    import glob
    import astropy.io.fits as fits
    import numpy as np
    import matplotlib.pyplot as plt
    from reproject import reproject_interp
    from astropy.stats import sigma_clip
    from kcwi_trim import make_wave
    print('Correcting Overlap')
    fnames=glob.glob(directory+'/*fits')

    anchor=fits.open(fnames[0])[0].data
    anchor_header=fits.open(fnames[0])[0].header

    wave=make_wave(anchor_header)

    for name in fnames[1:]:
        data=fits.open(name)[0].data
        header=fits.open(name)[0].header

        # Reproject image to correct onto anchor wcs/grid
        uncorrected_image, footprint = reproject_interp((data,header),  anchor_header)

        diff=np.nanmedian(anchor-uncorrected_image, axis=(1,2))
        plt.plot(wave,diff)
        plt.show()

        out_fits = fits.HDUList([fits.PrimaryHDU(data=data+np.nanmedian(diff), header=header)])
        out_fits.writeto(name, overwrite=True)
    return


def join_fits(science, variance, outname):
    import numpy as np
    import astropy.io.fits as fits
    sci = fits.open(science)
    var = fits.open(variance)

    hdu1 = fits.PrimaryHDU(data=sci[0].data, header=sci[0].header)
    hdu2 = fits.ImageHDU(data=var[0].data, header=var[0].header)

    w = np.abs(sci[0].data / np.sqrt(var[0].data)) < 0.01

    pix_mask = np.zeros(np.shape(sci[0].data))
    pix_mask[w] = 1

    hdu3 = fits.ImageHDU(data=pix_mask, header=sci[0].header)
    sci.close()
    var.close()

    new_hdu = fits.HDUList([hdu1, hdu2, hdu3])
    new_hdu.writeto(outname, overwrite=True)
    return print('The multi-cube fits file has been created.')

def cosmic_ray_reject(pairs, raw_dir,out_dir):
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.table as table
    import astropy.io.fits as fits
    plt.rcParams["figure.figsize"] = (10, 10)
    import scipy.stats
    files = np.loadtxt(pairs, dtype=str)
    for i in range(len(files)):
        f1 = raw_dir + files[i][0]
        f2 = raw_dir + files[i][1]

        f1_out = out_dir + files[i][0]
        f2_out = out_dir + files[i][1]

        print('working on ', files[i][0], ' and ', files[i][1])
        f1_hud = fits.open(f1)
        f2_hud = fits.open(f2)

        ratio_map = np.log10(f1_hud[0].data / f2_hud[0].data)
        ysize, xsize = np.shape(ratio_map)
        mask = np.zeros((ysize, xsize))
        xbox_size = int(xsize / 12)
        ybox_size = int(ysize / 12)
        print(xbox_size, xsize, ybox_size, ysize)
        i = ybox_size
        while i <= ysize:
            j = xbox_size
            while j <= xsize:
                segment = ratio_map[i - ybox_size:i + ybox_size, j - xbox_size:j + xbox_size]

                med = np.nanmedian(segment)
                std = np.nanstd(segment)
                # print(med, std)

                wmin = np.where(segment < med - 1.5 * std)
                wmax = np.where(segment > med + 1.5 * std)
                mask[i - ybox_size:i + ybox_size, j - xbox_size:j + xbox_size][wmin] = -1
                mask[i - ybox_size:i + ybox_size, j - xbox_size:j + xbox_size][wmax] = 1
                j += xbox_size
            i += ybox_size

        top_mask = mask == 1
        bot_mask = mask == -1

        f1_hud[0].data[top_mask] = f2_hud[0].data[top_mask]
        f2_hud[0].data[bot_mask] = f1_hud[0].data[bot_mask]
        # plt.imshow(f1_hud[0].data)
        # plt.show()

        f1_hud.writeto(f1_out, overwrite=True)
        f2_hud.writeto(f2_out, overwrite=True)

    return

def make_head_2d(head):
    nhead=head.copy()
    try:
        nhead['NAXIS'] = 2
        nhead.remove('NAXIS3')
        nhead.remove('CDELT3')
        nhead.remove('CRVAL3')
        nhead.remove('CRPIX3')
        nhead.remove('CD3_3')
        return nhead
    except:
        return nhead


def WcsUpdate(f1, sx, sy, sr, raw_dir):
    import astropy.io.fits as fits
    import numpy
    def pathleaf(path):
        import ntpath
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)
    hdu = fits.open(f1)
    name = pathleaf(f1)
    vhdu = fits.open(raw_dir+'/working_directory/vcubes/projected/'+ name[:-17]+'vcubes.vtrims.fits')

    hdu[0].header['CRVAL1'] = hdu[0].header['CRVAL1'] - sx / numpy.cos(numpy.radians(hdu[0].header['CRVAL2']))
    hdu[0].header['CRVAL2'] = hdu[0].header['CRVAL2'] + sy
    hdu[0].header['CROTA2'] = hdu[0].header['CROTA2'] - sr

    vhdu[0].header['CRVAL1'] = hdu[0].header['CRVAL1']
    vhdu[0].header['CRVAL2'] = hdu[0].header['CRVAL2']
    vhdu[0].header['CROTA2'] = hdu[0].header['CROTA2']

    try:
        hdu['dec_off'] = sy
        hdu['ra_off'] = sx
        hdu['rot_off'] = sr
    except:
        print('No Offsets Stored')
    out_dir=raw_dir+'/working_directory/scubes/aligned'
    vout_dir = raw_dir + '/working_directory/vcubes/aligned'
    new_name = out_dir + '/' + name[:-5] + '.fits'
    vnew_name = vout_dir + '/' + name[:-17] + 'vcubes.fits'

    hdu.writeto(new_name, overwrite=True, output_verify='ignore')
    hdu.close()

    vhdu.writeto(vnew_name, overwrite=True, output_verify='ignore')
    vhdu.close()
    print('The file has been aligned and saved as, ' + new_name)
    return

def check_cdelta(old_cube, new_cube):
    import astropy.io.fits as fits
    new=fits.open(new_cube)
    new = fits.open(old_cube)
	new[0].header['CRVAL3']=old[0].header['CRVAL3']
	new[0].header['CDELT3']=old[0].header['CDELT3']
	new[0].header['CRPIX3']=1.0
	new.writeto(new_cube, overwrite=True)
	new.close()
    return