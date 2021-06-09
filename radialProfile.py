import numpy as np
import yt
from yt.units import kpc
import argparse
parser = argparse.ArgumentParser();

parser.add_argument('--prefix', type = str, nargs = 1, help = 'name of the plot file');
parser.add_argument('--nmin', type = int, nargs = 1, help = 'name of the plot file');
parser.add_argument('--nmax', type = int, nargs = 1, help = 'name of the plot file');

args = parser.parse_args();
prefix = args.prefix[0];
nmin = args.nmin[0];
nmax = args.nmax[0];

i = nmin; 

while(i<=nmax):
	if (i<10):
		plotFile = "{}_000{}".format(prefix, i)
	elif(i<100):
		plotFile = "{}_00{}".format(prefix, i);
	else:
		plotFile = "{}_0{}".format(prefix, i);

	ds = yt.load(plotFile); 
	ad = ds.all_data();
	disc = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 10*kpc, 3.5*kpc);
	plot = yt.ProfilePlot(disc, "radius", ["density"], accumulation = False, weight_field = None, x_log= False);
	plot.set_unit("radius", "kpc")
	plot.set_ylim('density', 1e-24, 1e-18);
	plot.save();
	i = i+1;

