#This is used to pick up a random cloud and then plot particle plots
#starting 100 Myrs before the current time (given in the arguments)
#The cloud particles are highlighted in black 

#Important argument -> Linking length
#Works in three major steps
#1. Making a KDTree for spatial awareness 
#2. Running a query ball that gets all the points that lie close to each point
#3. My FOF where I make disjoint sets
#Result disjointSets[i]: containing the indices of the points of the original data that lie 
#in the ith disjointSet. 
import numpy as np 
import scipy.spatial as sp
import os
home = os.environ["HOME"]
import time
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mtlb
from scipy.spatial import KDTree
import yt
from yt.units import kpc, pc
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import h5py 
import gc
import math 




one_kpc = 3.08567E+21
one_parsec = one_kpc*1e-3
avogadro_no = 6.023e23
mean_mass = 1.3;
# scaleHeight = 100*pc;


#=================================================
#Arguments for the fof algorithm 
#=================================================

import argparse
parser = argparse.ArgumentParser();


parser.add_argument('--plt', type = int, nargs = 1, required = True , help = 'the plot file numbers to read, separated by space');
parser.add_argument('--densityCut', type = float, nargs = 1, required = True , help = 'the threshold density of the particles to plot (in cm^-3)');
parser.add_argument('--ll', type = float, nargs = 1, required = True , help = 'the linkingLength for binning in constant phi bins');
#parser.add_argument('--weighted', type = bool, nargs = 1, required = False , help = 'Whether to have mass-weighted particle density values or not.');
#parser.add_argument('--radius', type = float, nargs = 1, required = True, help = 'The radius at which to take the fourier transform (in kpc)');

args = parser.parse_args();

densityCutOff = args.densityCut[0] ; #particles/cm3
linkingLength = args.ll[0]*one_kpc; #
plt_number = args.plt[0];

#====================================================


#particlePositions[i] in the form [px[i], py[i], pz[i]]
#linkinLength should be in same units as positions
#cellThreshold == The minimum size of the disjoint set
#returns[i] : Indices of all the particles that lie in one cloud/cluster
def getDisjointSets(*, particlePositions, particleDensities, tags,densityCut, linkingLength, cellThreshold = 10):
    
    indices = np.where(particleDensities>=densityCut)[0];
    points_cut = particlePositions[indices];
    

    tree = KDTree(points_cut)

#Find all the points whose distance from each other is at most r
    clusters = tree.query_ball_tree(tree, linkingLength);


    visited = np.zeros(len(clusters));
    disjointSets = [];
    for i in range(len(clusters)):
        clustertoAppend = clusters[i];
        if(visited[i] == 0):
            for x in clustertoAppend:
                for y in clusters[x]:
                    if(clustertoAppend.count(y) == 0):
                        clustertoAppend.append(y);
    
                visited[x] = 1;
            
            if(len(clustertoAppend) > cellThreshold): 
                tags_to_append = tags[indices[clustertoAppend]] #converting back to original indices 
                disjointSets.append(tags_to_append); 

    return disjointSets;


def plotCloud(set_, ax):
    toPlot = set_;
    ax.scatter(x[toPlot]/one_kpc,y[toPlot]/one_kpc , c = particle_density[toPlot], cmap = cmap, norm = normalizer,s= 1, edgecolor = 'face');    
    return ;

def get_particleID(tags, tag_to_track):

    pid = [];

    for i in range(len(tag_to_track)):
        pid.append(np.where(tag_to_track[i] == tags)[0][0]);

    pid = np.asarray(pid);
    return pid;

#Plot all the particles in a scatter plot
def plotParticle(fig, ax, ds, tag_to_track):

    ad = ds.all_data();
    
    disc = ds.disk(ds.domain_center, [0.0, 0.0, 1], 13*kpc, 5000*pc);

    x = disc['all', 'particle_position'][:,0].v/(one_parsec*1000)
    y = disc['all', 'particle_position'][:,1].v/(one_parsec*1000)
    # z = disc['all', 'particle_position'][:,2].v/(one_parsec*1000)   
    tag = ad['particle_tag'].v;

    #matching the tags with the particle_IDs of the file we just read.
    p_id = get_particleID(tag, tag_to_track);
    
    print(p_id);
    particle_density = (disc['all', 'particle_pden']*avogadro_no)/mean_mass;

    #rounding off time to two decimal places... 
    time = math.trunc((ds.current_time.v*100)/3.154e13)
    time = time/100 

    colorbar_label = 'Number Density [$cm^{-3}$]';
    
    ax.set_xlabel('x (kpc)');
    ax.set_ylabel('y (kpc)');
    ax.minorticks_on();
    
    im = ax.scatter(x, y, c = particle_density, s = 0.01);
    ax.scatter(x[p_id], y[p_id], c= 'black', s= 0.05);
    ax.text(np.amin(x),np.amin(y), "T = {} Myrs".format(time), size = 10);

    clb = fig.colorbar(im);
    clb.set_label(colorbar_label);

    return;



ds = yt.load('Disc_hdf5_part_0{}'.format(plt_number));
ad = ds.all_data();
density = ad['particle_pden'];
particle_density = (density.v*avogadro_no)/mean_mass;
x = ad['particle_position_x'].v;
y = ad['particle_position_y'].v;
z = ad['particle_position_z'].v;
tag = ad['particle_tag'].v;
particlePositions = np.transpose([x,y,z]);

# ds = yt.load('Disc_hdf5_plt_cnt_0{}'.format(plt_number));
# ad = ds.all_data();
# density = ad['density'];
# particle_density = (density.v*avogadro_no)/mean_mass;
# x = ad['x'].v;
# y = ad['y'].v;
# z = ad['z'].v;
# particlePositions = np.transpose([x,y,z]);

#====================================================
#making a radial cut for tracking the feathers 
#===================================================
# r = np.sqrt(x**2 + y**2);
# indices = np.where((r>5.5*one_kpc) & (r < 8*one_kpc))[0];

# x = x[indices];
# y = y[indices];
# z = z[indices]
# particlePositions = particlePositions[indices];
# particle_density = particle_density[indices];
# tag = tag[indices]



print('calculating the sets');
clouds_tags = getDisjointSets(particlePositions=particlePositions, particleDensities = particle_density, tags = tag,densityCut=densityCutOff, linkingLength = linkingLength, cellThreshold = 100);



print('Number of clouds found - ', len(clouds_tags));

#======================================
#1. Take all the IDs of the particles in one cloud
#2. Go to ~ 100 Myrs before. 
#3. Just plot the particle IDs in each plot
#4. Save and exit
#Since you are opening up a lot of particle files, you also have to 
#free up some memory each time I think
#======================================

print('making animation ')

tag_to_track = clouds_tags[0];

nmin = plt_number - 80 ;
nmax = plt_number;

for i in range(nmin, nmax): 
    if (i<10):
        ds= yt.load('Disc_hdf5_part_000{}'.format(i));
    elif(i<100): 
        ds= yt.load('Disc_hdf5_part_00{}'.format(i));
    else:
        ds= yt.load('Disc_hdf5_part_0{}'.format(i));
    fig, ax = plt.subplots();
    plotParticle(fig, ax, ds, tag_to_track);
    plt.savefig('animation_particles_feather_{}.png'.format(i));

    # #collecting memory from the previously opened plot files
    # del ds 
    # gc.collect();
