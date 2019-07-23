from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import fitsio
import desimodel.io
import desitarget.mtl
from desitarget.targetmask import desi_mask, obsconditions

data = fitsio.FITS('targets.fits', 'r')
target_data = data[1].read(columns=['TARGETID', 'DESI_TARGET', 'MWS_TARGET', 'BGS_TARGET', 'SUBPRIORITY', 'NUMOBS_INIT', 'PRIORITY_INIT', 'RA', 'DEC'])
data.close()

mtl = desitarget.mtl.make_mtl(target_data)
mtl_file = 'mtl_targets.fits'
mtl.write(mtl_file, overwrite=True)
