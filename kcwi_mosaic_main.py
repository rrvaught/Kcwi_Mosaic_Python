import os
import sys
import shutil

from MontagePy.main    import *
from MontagePy.archive import *
from kcwi_trim import *
from kcwi_coadd import coadd
from kcwi_repro import repro
from kcwi_helper_functions import *
import glob

# Input Information
obj_name='leop_rm8895'
rawdir='/Users/rrickardsvaught/Desktop/Data/KCWI/leop/rm8895'
side='blue' # or 'red'
fnames=glob.glob(rawdir+'/*icubes.fits')

# Indicate which steps to run
do_trim=1  # 1=on, 0 is off
do_repro=1 # 1=on, 0 is off
align=0 # 1=on, 0 is off ## Not Implemented yet ##
clean_dir=1 # Delete all intermediate files

# Begin Main Section of Code

# housekeeping - Unpack the DRP output multi-extension cubes and save in working
#directories: scubes, vcubes and mcubes.
unpack_raw_cubes(fnames,rawdir)

# make required folders
make_dir(rawdir+'/working_directory/scubes/raw')
make_dir(rawdir+'/working_directory/scubes/projected')
make_dir(rawdir+'/working_directory/scubes/final')

make_dir(rawdir+'/working_directory/vcubes/raw')
make_dir(rawdir+'/working_directory/vcubes/projected')
make_dir(rawdir+'/working_directory/vcubes/final')


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
    clean_overlap(out_directory)

    directory = rawdir + '/working_directory/vcubes/raw/'
    out_directory = rawdir + '/working_directory/vcubes/projected/'
    repro(directory, out_directory)
else:
    print('Skipping 1st Pass Reprojection')



#### Co-Add Cubes ####
directory=rawdir + '/working_directory/scubes/projected/'
out_directory = rawdir + '/working_directory/scubes/final/'
coadd(directory,out_directory,obj_name=obj_name)

directory=rawdir + '/working_directory/vcubes/projected/'
out_directory = rawdir + '/working_directory/vcubes/final/'
coadd(directory,out_directory,obj_name=obj_name)

sin=rawdir + '/working_directory/scubes/final/'+obj_name+'_mosaic.fits'
vin=rawdir + '/working_directory/vcubes/final/'+obj_name+'_mosaic.fits'
outname=rawdir + '/'+obj_name+'.fits'
join_fits(sin, vin,outname)

### Clean directories ###
if clean_dir:
    os.system('rm -r '+rawdir+'/working_directory/')