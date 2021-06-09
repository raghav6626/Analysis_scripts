#do this with 10 cores. 
import numpy as np
import matplotlib.pyplot as plt
import yt
from mpi4py import MPI 
import argparse
import math

#comm = MPI.COMM_WORLD
#size = comm.Get_size()
#my_rank = comm.Get_rank();
#comm.Barrier()

parser = argparse.ArgumentParser();

parser.add_argument('--nmin', type = int, nargs = 1, help = 'starting number for plot file');
parser.add_argument('--nmax', type = int, nargs = 1, help = 'last number for the plot file'); 

parser.add_argument('--grid', type = bool, nargs = 1, help = 'to overplot grid or not'); 
parser.add_argument('--name', type = str, nargs = 1, help = 'a suffix to add to the plotFile name'); 

args = parser.parse_args();

nmi = args.nmin[0];
nma = args.nmax[0];

if(args.grid!=None):
    grid = args.grid[0];
    print(grid)
else:
    grid = False;
for i in range(nmi, nma):
    if (i<10):
        ds= yt.load('Disc_hdf5_plt_cnt_000{}'.format(i));
    elif(i<100): 
        ds= yt.load('Disc_hdf5_plt_cnt_00{}'.format(i));
    else:
        ds= yt.load('Disc_hdf5_plt_cnt_0{}'.format(i));
    prj = yt.ProjectionPlot(ds, 'z', 'density', width = (20, 'kpc'));
#    prj.set_zlim('density', 1e-6, 1e-1)
#    if(args.particles[0]):
#        prj.annotate_particles( (3, 'kpc'))
    prj.annotate_timestamp();
    if(grid == True):
        print('annotating grids...')
        prj.annotate_grids()
    if(args.name!=None):
        prj.save('{}_{}'.format(args.name[0], i));        
    else:
        prj.save();

