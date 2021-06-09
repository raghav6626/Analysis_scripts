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
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import h5py 


#=================================================
#Arguments for the fof algorithm 
#=================================================


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


parser.add_argument('--fileName', type = str, nargs = 1, required = True , help = 'the plot file name (assuming the data file is in the same directory)');
parser.add_argument('--densityCut', type = float, nargs = 1, required = True , help = 'the threshold density of the particles (in code units)');
parser.add_argument('--ll', type = float, nargs = 1, required = True , help = 'the linkingLength for the friends of friends algo (in code units)');
parser.add_argument('--cellThreshold', type = int, nargs = 1, required = True , help = 'the minimum number of particles/cells you want a clump to have');



args = parser.parse_args();

densityCutOff = args.densityCut[0];
linkingLength = args.ll[0]; #
fileName = args.fileName[0];
cellThreshold = args.cellThreshold[0];
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
    pid_sets = [];
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
                pid_sets.append(indices[clustertoAppend]);
                disjointSets.append(tags_to_append); 

    return disjointSets, pid_sets;


def plotCloud(set_, ax):
    toPlot = set_;
    ax.scatter(x[toPlot],y[toPlot] , c = particle_density[toPlot], cmap = cmap, norm = normalizer,s= 0.5, edgecolor = 'face');    
    return ;

dirr = './'

ds = yt.load('{}{}'.format(dirr, fileName));

ad = ds.all_data();
density = ad['particle_pden'];
#particle_density = (density.v*avogadro_no)/mean_mass;
particle_density = density; 
x = ad['particle_position_x'].v;
y = ad['particle_position_y'].v;
z = ad['particle_position_z'].v;
particlePositions = np.transpose([x,y,z]);
tag = ad['particle_tag'].v;


(clouds_tags, cloud_p_ids) = getDisjointSets(particlePositions=particlePositions, particleDensities = particle_density, tags = tag,densityCut=densityCutOff, linkingLength = linkingLength, cellThreshold = cellThreshold);

print('Number of clouds/clumps found - ', len(clouds_tags));


#============================================
#Writing the particle_IDs and properties in a file
#============================================

print('Writing the particle IDs to the file - clouds_{}'.format(fileName));

hf = h5py.File('clouds_{}'.format(fileName), 'w');

for i in range(len(clouds_tags)):
	g = hf.create_group('cloud_{}'.format(i))
	g.create_dataset('cloud_ID', data = i);
	g.create_dataset('particle_tags', data = clouds_tags[i]);
    # g.create_dataset('particle_IDs', data = cloud_p_ids[i]);
	g.create_dataset('particle_density', data = particle_density[cloud_p_ids[i]]);


hf.close();
#We need a common colorbar for all the clouds. So we normalize 

#=================================================
#Plotting in a lame plot
#================================================

print('Plotting a scatter plot of the particles')

# vmax = np.amax(particle_density);
vmax = np.amax(particle_density);

vmin = densityCutOff;

cmap=cm.get_cmap('viridis')
normalizer=Normalize(vmin, vmax);

im=cm.ScalarMappable(norm = normalizer, cmap = cmap)


fig, ax = plt.subplots(figsize=(12, 10))
plt.style.use('classic');
for set_ in cloud_p_ids:
    plotCloud(set_, ax);

plt.xlim(-11, 11)
plt.ylim(-11, 11);

clb = fig.colorbar(im);
clb.ax.minorticks_on();
clb.set_label('Number Density (cm$^{-3}$)', fontsize = 17)
#plt.colorbar(im);
plt.savefig('fof_clouds_{}.png'.format(fileName));