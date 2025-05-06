def path_leaf(path):
    import ntpath
    head, tail = ntpath.split(path)
    return tail , head

def coadd(directory,out_directory,obj_name,is_var=False):
    from MontagePy.main import mImgtbl
    from MontagePy.main import mMakeHdr
    from MontagePy.main import mAddCube
    from kcwi_trim import ster2pix
    import glob

    # Read in Images
    # Write Montage header file
    ITable = mImgtbl(directory, directory+'add.tbl')
    print("mImgtbl (raw):    " + str( ITable ), flush=True)

    rtn = mMakeHdr(directory+'add.tbl', directory+'add.hdr')
    with open(directory+'add.hdr', 'r') as fin:
        print(fin.read(), end='')



    rtn = mAddCube(directory,
                       directory+'add.tbl',
                       directory+'first.hdr',
                       out_directory+obj_name+'_mosaic.fits',coadd=0)
    print(rtn)
    return