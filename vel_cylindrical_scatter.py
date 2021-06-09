import numpy as np
import yt 
from yt.units import kpc, pc
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser();


parser.add_argument('--name', type = str, nargs = 1, help = 'file to load');
#parser.add_argument('--field', type = str, nargs = 1, help = 'name of the field to plot.');

args = parser.parse_args();


ds = yt.load(args.name[0]);


ad = ds.all_data();
disc = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 5*kpc, 100*pc);

radius = disc['radius'].v
v_phi = disc['velocity_cylindrical_theta']



fig, ax = plt.subplots();
plt.style.use('default');
ax.scatter(radius/(3.086e+21), v_phi/(1e5), s = 0.1);
ax.set_xlabel('R (kpc)');
ax.set_ylabel('V$_{\phi}$ (km/s)');
plt.savefig('velocity_cylindrical_scatterPlot');
