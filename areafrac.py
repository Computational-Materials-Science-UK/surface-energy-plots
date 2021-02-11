from wulffpack import SingleCrystal
from ase.build import bulk
from ase.io import write
import numpy as np
import scipy.io as spio
import pandas as pd

mat = spio.loadmat('data.mat',squeeze_me=True)
matrix = mat['allSEs']

xx_list = []
yy_list = []
zz_list = []



for i in range(33366):
	x = matrix[i,3]
	y = matrix[i,4]
	z = matrix[i,5]
	

	# Specify parameters
	prim = bulk('W',
	            crystalstructure='bcc',
	            a=3.172)
	surface_energies = {(1, 0, 0): x,
	                    (1, 1, 0): y,
	                    (2, 1, 1): z}
	
	goodshape = SingleCrystal(surface_energies)
	#goodshape.view()
	#write('single_crystal.xyz', goodshape.atoms)
		
	                 
	#print('Fraction of {100} facets:')
	for name, particle in zip(['goodshape'], [goodshape]):
	    fraction = particle.facet_fractions.get((1, 0, 0), 0.)
	    xx = fraction
	    xx_list.append(xx)
	    #print('{}: {:.4f}'.format(name, xx))
	#print()
	
	#print('Fraction of {110} facets:')
	for name, particle in zip(['goodshape'], [goodshape]):
	    fraction = particle.facet_fractions.get((1, 1, 0), 0.)
	    yy = fraction
	    yy_list.append(yy)
	    #print('{}: {:.4f}'.format(name, yy))
	#print()
	
	#print('Fraction of {211} facets:')
	for name, particle in zip(['goodshape'], [goodshape]):
	    fraction = particle.facet_fractions.get((2, 1, 1), 0.)
	    zz = fraction
	    zz_list.append(zz)
	    #print('{}: {:.4f}'.format(name, zz))
	#print()
	print(i)
	
df = pd.DataFrame(matrix)
df['xx'] = xx_list
df['yy'] = yy_list
df['zz'] = zz_list
#print(df)
df.to_csv(r'/Users/mujanseif/Documents/Research/cathodes/W-surfaces/thermal-properties/surface-energy-plots/alldata.csv',\
index = False, header=False)
	
	
	