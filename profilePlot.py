import numpy as np
import yt
from yt.units import kpc
import argparse
parser = argparse.ArgumentParser();

parser.add_argument('--name', type = int, nargs = 1, help = 'name of the plot file');
parser.add_argument('--field', type = str, nargs = 1, help = 'name of the field to plot');
parser.add_argument('--bin', type = str, nargs = 1, help = 'name of the field to be binned in (the x-axis)');

args = parser.parse_args();
plotFile = args.name[0];


ds = yt.load(plotFile);
ad = ds.all_data();
disc = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 10*kpc, 3*kpc);

plot = yt.ProfilePlot(disc, "{}".format(args.bin[0]), ["{}".format(args.field[0])], accumulation = False, weight_field = None, x_log= False);
plot.set_unit("radius", "kpc")
plot.save();

