from astropy.table import Table
import numpy as np
import subprocess


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

for i_pass in passes:
    pass_tilefile = 'pass_{:02d}_{}'.format(i_pass, tilefile)
    targets_file = 'mtl_targets.fits'
    fba_output_dir = 'fiberassign_pass_{:02d}'.format(i_pass)
    
    p = subprocess.call(["fba_run", "--targets", targets_file,'--sky', 'sky.fits',
                             '--footprint', pass_tilefile,  
                             '--dir',  fba_output_dir, '--overwrite'])