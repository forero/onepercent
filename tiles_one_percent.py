from astropy.table import Table, hstack
import numpy as np
import matplotlib.pyplot as plt
import desimodel.io
from astropy.time import Time

tiles = desimodel.io.load_tiles()

kk = (tiles['RA'] > 240) & (tiles['RA']<260) & (tiles['DEC']>20.0) & (tiles['DEC']<30) 
kk &= tiles['PASS']==0
print(np.count_nonzero(kk))
print(set(tiles['PASS'][kk]))
pattern = tiles[kk].copy()

delta_ra = 0.3
delta_dec = 0.01
l = list()
for i_pass in range(20):
    tmp_pattern = pattern.copy()
    tmp_pattern['RA'] = pattern['RA'] + i_pass*delta_ra
    tmp_pattern['DEC'] = pattern['DEC'] + i_pass*delta_dec
    tmp_pattern['PASS'] = i_pass
    l.append(tmp_pattern.copy())
#print(l)
full_pattern = np.concatenate(l)
n_tiles = len(full_pattern)
full_pattern['TILEID'] = np.arange(n_tiles, dtype='int')
ii_pass = full_pattern['PASS']>=10
full_pattern['PROGRAM'][ii_pass] = 'BRIGHT'
full_pattern['OBSCONDITIONS'][ii_pass] = 4

Table(full_pattern).write('1_percent_20_pass.fits', overwrite=True)
