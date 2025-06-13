def path_leaf(path):
    import ntpath
    head, tail = ntpath.split(path)
    return tail , head

def repro(directory,out_directory):
    from MontagePy.main import mImgtbl
    from MontagePy.main import mMakeHdr
    from MontagePy.main import mProjectCube
    from kcwi_trim import ster2pix
    import glob
    import os
    from kcwi_helper_functions import check_cdelta

    # Read in Images
    tname = glob.glob(directory+'*.fits')
    # Write Montage header file
    ITable = mImgtbl(directory, out_directory+'first.tbl')
    print("mImgtbl (raw):    " + str( ITable ), flush=True)

    rtn = mMakeHdr(out_directory+'first.tbl', out_directory+'first.hdr')
    with open(out_directory+'first.hdr', 'r') as fin:
        print(fin.read(), end='')

    # Project KCWI image to Square Pixels via Montage
    for name in tname:
        print('Projecting {0}'.format(name))
        _=ster2pix(name,'sqrarc2ster')

        file_name, _ = path_leaf(name)

        rtn = mProjectCube(name,
                       out_directory+file_name,
                       out_directory+'first.hdr')#,fluxScale=scale)
        # ensure crpix3 and cdelt3 have been copied

        _ = ster2pix(out_directory+file_name, 'ster2sqrarc')
        check_cdelta(name, out_directory+file_name)
    # clean area files
    files = glob.glob(out_directory+'*_area.fits')
    for f in files:
        os.remove(f)
    return