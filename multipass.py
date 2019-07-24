from astropy.table import Table
import numpy as np
import subprocess
import glob
import os


from desisim.quickcat import get_median_obsconditions
import yaml
from collections import Counter
from pkg_resources import resource_filename
from time import asctime
import desitarget.mtl 

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, vstack
import sys
import scipy.special as sp
import desisim
from desisim.targets import get_simtype

import astropy.constants
c = astropy.constants.c.to('km/s').value

from desitarget.targetmask import desi_mask, bgs_mask, mws_mask

from desiutil.log import get_logger
log = get_logger()


def quickcat(tilefiles, targets, truth, zcat=None, obsconditions=None, perfect=False):
    """
    Generates quick output zcatalog
    Args:
        tilefiles : list of fiberassign tile files that were observed
        targets : astropy Table of targets
        truth : astropy Table of input truth with columns TARGETID, TRUEZ, and TRUETYPE
        zcat (optional): input zcatalog Table from previous observations
        obsconditions (optional): Table or ndarray with observing conditions from surveysim
        perfect (optional): if True, treat spectro pipeline as perfect with input=output,
            otherwise add noise and zwarn!=0 flags
    Returns:
        zcatalog astropy Table based upon input truth, plus ZERR, ZWARN,
        NUMOBS, and TYPE columns
    """
    #- convert to Table for easier manipulation
    if not isinstance(truth, Table):
        truth = Table(truth)

    #- Count how many times each target was observed for this set of tiles
    log.info('{} QC Reading {} tiles'.format(asctime(), len(tilefiles)))
    nobs = Counter()
    targets_in_tile = {}
    tileids = list()
    for infile in tilefiles:
        
        fibassign, header = fits.getdata(infile, 'FASSIGN', header=True)
 
        # hack needed here rnc 7/26/18
        if 'TILEID' in header:
            tileidnew = header['TILEID']
        else:
            fnew=infile.split('/')[-1]
            tileidnew=fnew.split("_")[-1]
            tileidnew=int(tileidnew[:-5])
            log.error('TILEID missing from {} header'.format(fnew))
            log.error('{} -> TILEID {}'.format(infile, tileidnew))
            
        tileids.append(tileidnew)

        ii = (fibassign['TARGETID'] != -1)  #- targets with assignments
        nobs.update(fibassign['TARGETID'][ii])
        targets_in_tile[tileidnew] = fibassign['TARGETID'][ii]


    #- Trim truth down to just ones that have already been observed
    log.info('{} QC Trimming truth to just observed targets'.format(asctime()))
    obs_targetids = np.array(list(nobs.keys()))
    iiobs = np.in1d(truth['TARGETID'], obs_targetids)
    truth = truth[iiobs]
    targets = targets[iiobs]

    #- Construct initial new z catalog
    log.info('{} QC Constructing new redshift catalog'.format(asctime()))
    newzcat = Table()
    newzcat['TARGETID'] = truth['TARGETID']
    if 'BRICKNAME' in truth.dtype.names:
        newzcat['BRICKNAME'] = truth['BRICKNAME']
    else:
        newzcat['BRICKNAME'] = np.zeros(len(truth), dtype=(str, 8))

    #- Copy TRUETYPE -> SPECTYPE so that we can change without altering original
    newzcat['SPECTYPE'] = truth['TRUESPECTYPE'].copy()

    #- Add ZERR and ZWARN
    log.info('{} QC Adding ZERR and ZWARN'.format(asctime()))
    nz = len(newzcat)
    newzcat['Z'] = truth['TRUEZ'].copy()
    newzcat['ZERR'] = np.zeros(nz, dtype=np.float32)
    newzcat['ZWARN'] = np.zeros(nz, dtype=np.int32)
   
    #- Add numobs column
    log.info('{} QC Adding NUMOBS column'.format(asctime()))
    newzcat.add_column(Column(name='NUMOBS', length=nz, dtype=np.int32))
    for i in range(nz):
        newzcat['NUMOBS'][i] = nobs[newzcat['TARGETID'][i]]

    #- Merge previous zcat with newzcat
    log.info('{} QC Merging previous zcat'.format(asctime()))
    if zcat is not None:
        #- don't modify original
        #- Note: this uses copy on write for the columns to be memory
        #- efficient while still letting us modify a column if needed
        zcat = zcat.copy()

        #- targets that are in both zcat and newzcat
        repeats = np.in1d(zcat['TARGETID'], newzcat['TARGETID'])

        #- update numobs in both zcat and newzcat
        ii = np.in1d(newzcat['TARGETID'], zcat['TARGETID'][repeats])
        orig_numobs = zcat['NUMOBS'][repeats].copy()
        new_numobs = newzcat['NUMOBS'][ii].copy()
        zcat['NUMOBS'][repeats] += new_numobs
        newzcat['NUMOBS'][ii] += orig_numobs

        #- replace only repeats that had ZWARN flags in original zcat
        #- replace in new
        replace = repeats & (zcat['ZWARN'] != 0)
        jj = np.in1d(newzcat['TARGETID'], zcat['TARGETID'][replace])
        zcat[replace] = newzcat[jj]

        #- trim newzcat to ones that shouldn't override original zcat
        discard = np.in1d(newzcat['TARGETID'], zcat['TARGETID'])
        newzcat = newzcat[~discard]

        #- Should be non-overlapping now
        assert np.all(np.in1d(zcat['TARGETID'], newzcat['TARGETID']) == False)

        #- merge them
        newzcat = vstack([zcat, newzcat])

    #- check for duplicates
    targetids = newzcat['TARGETID']
    assert len(np.unique(targetids)) == len(targetids)

    #- Metadata for header
    newzcat.meta['EXTNAME'] = 'ZCATALOG'

    log.info('{} QC done'.format(asctime()))
    return newzcat


tilefile = '1_percent_20_pass.fits'
tiles = Table.read(tilefile)

passes = list(set(tiles['PASS']))
passes.sort()
print('passes:', passes)

#split tileids into passes
tiles_pass = dict()

for i_pass in passes:
    pass_tilefile = 'pass_{:02d}_{}'.format(i_pass, tilefile)
    print(pass_tilefile)
    tiles[tiles['PASS']==i_pass].write(pass_tilefile, overwrite=True)
    tiles_pass[i_pass] = np.array(tiles['TILEID'][tiles['PASS']==i_pass])

#print(tiles_pass)
truth = Table.read("truth.fits")

for i_pass in passes:
    pass_tilefile = 'pass_{:02d}_{}'.format(i_pass, tilefile)
    input_targets_file = 'mtl_targets_{:02d}.fits'.format(i_pass)
    output_targets_file = 'mtl_targets_{:02d}.fits'.format(i_pass+1)
    fba_output_dir = 'fiberassign_pass_{:02d}'.format(i_pass)
    old_zcat_file = 'zcat_{:02}.fits'.format(i_pass-1)
    new_zcat_file = 'zcat_{:02}.fits'.format(i_pass)
    
    p = subprocess.call(["fba_run", "--targets", input_targets_file,'--sky', 'sky.fits',
                             '--footprint', pass_tilefile,  
                             '--dir',  fba_output_dir, '--overwrite'])
    
    fba_files = glob.glob(os.path.join(fba_output_dir,"fiberassign*.fits"))
    print(i_pass, fba_files)
    
    targets = Table.read(input_targets_file)
    if i_pass==0:
        zcat = quickcat(fba_files, targets, truth)
    else:
        old_zcat = Table.read(old_zcat_file)
        zcat = quickcat(fba_files, targets, truth, zcat=old_zcat, perfect=True)        
    
    zcat.write(new_zcat_file, overwrite=True)
    mtl = desitarget.mtl.make_mtl(targets, zcat)
    mtl.write(output_targets_file, overwrite=True)
        
    