import numpy as np
import yt 
from yt.units import kpc, pc
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser();


parser.add_argument('--name', type = str, nargs = 1, help = 'base name of the files');

parser.add_argument('--nmin', type = int, nargs = 1, help = 'starting number for plot file');
parser.add_argument('--nmax', type = int, nargs = 1, help = 'last number for the plot file'); 

#parser.add_argument('--field', type = str, nargs = 1, help = 'name of the field to plot.');

args = parser.parse_args();

nmi = args.nmin[0];
nma = args.nmax[0];


name = args.name[0]

fig, ax = plt.subplots();
plt.style.use('default');

for i in range(nmi, nma):
    if (i<10):
        ds= yt.load('{}000{}'.format(name, i));
    elif(i<100): 
        ds= yt.load('{}00{}'.format(name, i));
    else:
        ds= yt.load('{}0{}'.format(name, i));

    ad = ds.all_data();
    disc = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 2*kpc, 200*pc);

    radius_1 = disc['radius'].v
    v_phi_1 = disc['velocity_cylindrical_theta']
        
    ax.scatter(radius_1/(3.086e+21), v_phi_1/(1e5), s = 0.1, label = '{} Myr'.format(i));


plt.legend();
ax.set_xlabel('R (kpc)');
ax.set_ylabel('V$_{\phi}$ (km/s)');
plt.savefig('velocity_cylindrical_scatterPlot');
