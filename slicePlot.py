import numpy as np
import yt 
from yt.units import kpc
import argparse

parser = argparse.ArgumentParser();


#parser.add_argument('--name', type = str, nargs = 1, help = 'file to load');
parser.add_argument('--field', type = str, nargs = 1, help = 'name of the field to plot.');

parser.add_argument('--nmin', type = int, nargs = 1, help = 'starting number for plot file');
parser.add_argument('--nmax', type = int, nargs = 1, help = 'last number for the plot file'); 
parser.add_argument('--width', type = int, nargs = 1, help = 'The width of the plot window')




args = parser.parse_args();


nmi = args.nmin[0];
nma = args.nmax[0];


for i in range(nmi, nma):
    if (i<10):
        ds= yt.load('Disc_hdf5_plt_cnt_000{}'.format(i));
    elif(i<100): 
        ds= yt.load('Disc_hdf5_plt_cnt_00{}'.format(i));
    else:
        ds= yt.load('Disc_hdf5_plt_cnt_0{}'.format(i));


    slc = yt.SlicePlot(ds, 'z', fields = args.field[0], width = (args.width[0], 'kpc'));
    slc.annotate_quiver('velocity_x', 'velocity_y');
    slc.annotate_timestamp();

    slc.save()
