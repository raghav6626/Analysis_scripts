#do this with 10 cores. 
import numpy as np
import matplotlib.pyplot as plt
import yt
from mpi4py import MPI 
import argparse
import math
from yt.units import kpc, pc





def plotParticle(ds, plt_number, p_id_check):
    
    ad = ds.all_data();
    
    disc = ds.disk(ds.domain_center, [0.0, 0.0, 1], 13*kpc, 5000*pc);

    p_id = disc["particle_index"];

    x = disc['all', 'particle_position'][:,0].v/(one_parsec*1000)
    y = disc['all', 'particle_position'][:,1].v/(one_parsec*1000)
    z = disc['all', 'particle_position'][:,2].v/(one_parsec*1000)    

    
    particle_density = (disc['all', 'particle_pden']*avogadro_no)/mean_mass;

  
    indices = [];
    for i in range(len(particle_density)):
        if(p_id[i] in p_id_check):
            indices .append(i)
  

    #rounding off time to two decimal places...    
    time = math.trunc((ds.current_time.v*100)/3.154e13)
    time = time/100    
    
    colorbar_label = 'Number Density [$cm^{-3}$]';

    fig, ax = plt.subplots();
    ax.set_xlabel('x (kpc)');
    ax.set_ylabel('y (kpc)');
    ax.minorticks_on();
    
    im = ax.scatter(x[indices], y[indices], c = particle_density[indices], s = 0.01);
    ax.text(np.amin(x),np.amin(y), "T = {} Myrs".format(time), size = 10);
    ax.text(np.amin(x),np.amin(y)+ 1.5, "|z| < {} kpc".format(args.cutz[0]), size = 10);

    clb = fig.colorbar(im);
    clb.set_label(colorbar_label);
    if(args.threshDensity!=None):
        plt.savefig('particlePlot_{}_threshold_{}.png'.format(plt_number, args.threshDensity[0]));
    else:
        plt.savefig('particlePlot_{}.png'.format(plt_number));
    
    return;

#===============================
#Important Constants
#===============================

one_parsec = 3.086e+18 #in cms
avogadro_no = 6.022E23; 
mean_mass = 1.3;

#=============================================
#Reading the file arguments
#=============================================
parser = argparse.ArgumentParser();
#parser.add_argument('--nmin', type = int, nargs = 1, help = 'starting number for plot file');
parser.add_argument('--plt', type = int, nargs = 1, help = 'number for the plot file'); 
parser.add_argument('--cutr', type = float, nargs = 1, help = 'the minumum cylindrical radius from the center'); 
parser.add_argument('--cutz', type = float, nargs = 1, help = 'the z distance from the plane of the galaxy'); 
parser.add_argument('--threshDensity', type = float, nargs = 1, help = 'threshold (minimum) number density of the particles to plot. If not given, all the particles are plotted.'); 

args = parser.parse_args();

thresholdDensity = args.thresholdDensity[0];



#==============================================
i = args.plt_number[0]

if (i<10):
    ds= yt.load('Disc_hdf5_part_000{}'.format(i));
elif(i<100): 
    ds= yt.load('Disc_hdf5_part_00{}'.format(i));
else:
    ds= yt.load('Disc_hdf5_part_0{}'.format(i));

particle_density = (disc['all', 'particle_pden']*avogadro_no)/mean_mass;
indices_density = np.where(particle_density>= thresholdDensity)[0];

p_id = disc["particle_index"]
p_id = p_id[indices_density];

nmi = plt_number - 10 
nma = plt_number;

#================================================
for i in range(nmi, nma): 
    if (i<10):
        ds= yt.load('Disc_hdf5_part_000{}'.format(i));
    elif(i<100): 
        ds= yt.load('Disc_hdf5_part_00{}'.format(i));
    else:
        ds= yt.load('Disc_hdf5_part_0{}'.format(i));
    plotParticle(ds, i, particle_indices);




