from array import array
def path_leaf(path):
    import ntpath
    head, tail = ntpath.split(path)

    return tail, head

def ster2pix(fname,which):
    import astropy.io.fits as fits
    hdu=fits.open(fname)
    header=hdu[0].header
    try:
        sqrarc_pix = header['PXSCL'] * header['SLSCL'] * 3600.0 ** 2
    except:
        sqrarc_pix = header['CDELT1'] * header['CDELT1'] * 3600.0 ** 2
    ster_to_sqrarc=4.255e10
    convfac=sqrarc_pix / ster_to_sqrarc

    if which=='sqrarc2ster':
        hdu[0].data=hdu[0].data / convfac
    elif which=='ster2sqrarc':
        hdu[0].data = hdu[0].data * convfac
    hdu.writeto(fname,overwrite=True)
    hdu.close()
    return convfac

def make_wave(header, z=0):
    import numpy as np
    """
    :param header:
    :param z:
    :return: restframe wavelength array
    """

    import numpy
    wcen=header['CRVAL3']
    length=header['NAXIS3']

    try:
        delta=header['CDELT3']
    except:
        delta=header['CD3_3']

    warray=np.arange(wcen,wcen+length*delta,delta)/(1+z)
    return warray


def trim(fname,is_var=False):
    """
    :param fname:
    :param is_var:
    :return:
    """
    import astropy.io.fits as fits
    import numpy as np

    hdu=fits.open(fname)
    header=hdu[0].header
    new_header=header.copy()
    data = hdu[0].data
    header.rename_keyword('CD3_3','CDELT3', force = False)
    # Montage likes CDELT3 Keyword
    wave=make_wave(header)

    # Trim Wavelength
    wlow, whigh = header['WAVGOOD0'],header['WAVGOOD1']
    gwave=np.logical_and(wave>=wlow,wave<=whigh)

    # Trim Data
    pady, padx=header['DARPADY']+3,header['DARPADX']
    new_data=data[gwave,pady:-pady,padx:-padx]

    new_header['NAXIS2'],new_header['NAXIS1']=np.shape(new_data)[1:]
    new_header['CRVAL3']=wave[gwave][0]
    new_header['CDELT3'] = header['CDELT3']

    tag='.trims.fits' if not is_var else '.vtrims.fits'
    out_fits=fits.HDUList([fits.PrimaryHDU(data=new_data,header=new_header)])

    file_name, head= path_leaf(fname)

    out_fits.writeto(head+'/raw/'+file_name[:-5]+tag,overwrite=True)
    print('Save Filed as '+file_name[:-5]+tag)

    return
# # Convert to SB units to insure flux conservation
    #sqrarc_pix=header['PXSCL']*header['SLSCL']*3600.0**2
    #ster_to_sqrarc=4.255e10