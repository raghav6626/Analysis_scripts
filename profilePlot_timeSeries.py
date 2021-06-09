import numpy as np
import yt
import yt.units as u
from yt.units import kpc, pc
import argparse
parser = argparse.ArgumentParser();

parser.add_argument('--name', type = str, nargs = 1, help = 'name of the plot file');
parser.add_argument('--field', type = str, nargs = 1, help = 'name of the field to plot');
parser.add_argument('--bin', type = str, nargs = 1, help = 'name of the field to be binned in (the x-axis)');

args = parser.parse_args();
#plotFile = args.name[0];

#this will be a glob pattern now.

es = yt.load('Disc_hdf5_plt_cnt_000?');

profiles = [];
labels = [];

x_field = "{}".format(args.bin[0])
y_field = "{}".format(args.field[0])


for ds in es:
	ad = ds.all_data();
	
	disc = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 5*kpc, 200*pc);

	
	profiles.append(yt.create_profile(disc, x_field,[y_field] , accumulation = False,  logs= {x_field: False, y_field: True} ));	
	
	
	labels.append('t = {}'.format(ds.current_time.v/(365*24*3600*1e6)));

plot = yt.ProfilePlot.from_profiles(profiles, labels = labels);

plot.set_unit("radius", "kpc")
plot.save();

