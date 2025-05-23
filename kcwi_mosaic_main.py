import os
import sys
import shutil

from MontagePy.main    import *
from MontagePy.archive import *
from kcwi_trim import *
from kcwi_coadd import coadd
from kcwi_align import align
from kcwi_repro import repro
from kcwi_helper_functions import *
import glob

# Input Information
obj_name='j0823_p2806'
rawdir='/Users/rrickardsvaught/Library/CloudStorage/Box-Box/JWST MIRI projects/Keck/220127/J0823_p2806'
side='blue' # or 'red'
fnames=glob.glob(rawdir+'/*icubes.fits')

# Indicate which steps to run
do_trim=1  # 1=on, 0 is off
do_repro=1 # 1=on, 0 is off
do_align=1 # 1=on, 0 is off
align_type='reference-image' # options single-field, pairs, reference-image
ref_image=rawdir+'/J0823_p2806_g.fits'
filter='SLOAN_SDSS.g.dat'

clean_dir=0 # Delete all intermediate files

####### Begin Main Section of Code #########

# housekeeping - Unpack the DRP output multi-extension cubes and save in working
#directories: scubes, vcubes and mcubes.
unpack_raw_cubes(fnames,rawdir)

# make required folders
make_directories(rawdir)

##### Trim trash pixels and bad wavelengths from cubes ####
if do_trim:
    # grab science frames
    print('Working on trim')
    snames=glob.glob(rawdir+'/working_directory/scubes/*icubes.fits')
    print(snames)
    for name in snames:
        print(name)
        trim(name)

    vnames = glob.glob(rawdir+'/working_directory/vcubes/*vcubes.fits')
    for name in vnames:
        trim(name, is_var=True)

elif not do_trim:
    print('Skipping trim step')

#### Reproject onto square pixels #####
if do_repro:
    print('Working on 1st Pass Reprojection')
    directory=rawdir + '/working_directory/scubes/raw/'
    out_directory = rawdir + '/working_directory/scubes/projected/'
    repro(directory,out_directory)
    directory = rawdir + '/working_directory/vcubes/raw/'
    out_directory = rawdir + '/working_directory/vcubes/projected/'
    repro(directory, out_directory)
else:
    print('Skipping 1st Pass Reprojection')

if do_align:
    align(align_type,rawdir,ref_image=ref_image, filter=filter)

    # Need to redo the projection with the aligned cubes #
    directory = rawdir + '/working_directory/scubes/aligned/'
    out_directory = rawdir + '/working_directory/scubes/final_project/'
    repro(directory, out_directory)

    directory = rawdir + '/working_directory/vcubes/aligned/'
    out_directory = rawdir + '/working_directory/vcubes/final_project/'
    repro(directory, out_directory)
    clean_overlap(out_directory)

    #### Co-Add Cubes ####
    directory = rawdir + '/working_directory/scubes/final_project/'
    out_directory = rawdir + '/working_directory/scubes/final/'
    coadd(directory, out_directory, obj_name=obj_name)

    directory = rawdir + '/working_directory/vcubes/final_project/'
    out_directory = rawdir + '/working_directory/vcubes/final/'
    coadd(directory, out_directory, obj_name=obj_name)

else:
    #### Co-Add Cubes ####
    directory = rawdir + '/working_directory/scubes/projected/'
    out_directory = rawdir + '/working_directory/scubes/final/'
    coadd(directory, out_directory, obj_name=obj_name)

    directory = rawdir + '/working_directory/vcubes/projected/'
    out_directory = rawdir + '/working_directory/vcubes/final/'
    coadd(directory, out_directory, obj_name=obj_name)


sin=rawdir + '/working_directory/scubes/final/'+obj_name+'_mosaic.fits'
vin=rawdir + '/working_directory/vcubes/final/'+obj_name+'_mosaic.fits'
outname=rawdir + '/'+obj_name+'.fits'
join_fits(sin, vin,outname)

### Clean directories ###
if clean_dir:
    os.system('rm -r '+rawdir+'/working_directory/')