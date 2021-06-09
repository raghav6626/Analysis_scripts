import numpy as np
import matplotlib.pyplot as plt
import sys
import statistics
import argparse
from scipy import constants as sc
import yt
from mpi4py import MPI
# yt.toggle_interactivity();
from numpy import linalg as LA
import time
import math

m_0 = 1.989e+33 #in grams
one_parsec = 3.086e+18 #in cms


parser = argparse.ArgumentParser();

parser.add_argument('--plt_number', type = int, nargs = 1, help = 'starting number for plot file');


args = parser.parse_args();
plt_number = args.plt_number[0];

if (plt_number<10):
	ds= yt.load('Disc_hdf5_plt_cnt_000{}'.format(plt_number));
elif(plt_number<100): 
	ds= yt.load('Disc_hdf5_plt_cnt_00{}'.format(plt_number));
else:
	ds= yt.load('Disc_hdf5_plt_cnt_0{}'.format(plt_number));



ad = ds.all_data();


import random

level = ad['index', 'grid_level'];
radius = ad['index', 'radius'];


x = np.amin(level.v);
index_3 = np.where(level == x);
index_5 = np.where(level == x+1);
index_7 = np.where(level == x+2);
index_8 = np.where(level == x+3);
index_9 = np.where(level == x+4);
print(np.amin(level), np.amax(level));

index_max = np.where(level == np.amax(level.v));
print('Number of blocks at the highest refinement level are = ',len(index_max[0]));

x3 = []
x5 = []
x7 = []
x8 = []
x9 = []
# print(len(index_3[0]))

#generate random 100000 positions out of the array indices available...
for i in range(1, 100000):
#     x3[i] = int(random.randint(0, len(index_3[0])));
    x3.append(random.randint(0, len(index_3[0]) -1));
    if(len(index_5[0])>0):
    	x5.append(random.randint(0, len(index_5[0]) -1));
    if(len(index_7[0])>0):
        x7.append(random.randint(0, len(index_7[0]) -1));
    if(len(index_8[0])>0):    
        x8.append(random.randint(0, len(index_8[0]) -1));
    if(len(index_9[0])>0):    
        x9.append(random.randint(0, len(index_9[0]) -1));

indices_three_final = index_3[0][x3];
if(len(index_5[0])>0):
    indices_five_final = index_5[0][x5];
if(len(index_7[0])>0):
    indices_seven_final = index_7[0][x7];
if(len(index_8[0])>0):  
    indices_eight_final = index_8[0][x8];
if(len(index_9[0])>0):  
    indices_nine_final = index_9[0][x9];
# print(level[indices_three
#             _final])
# print(level[indices_five_final])

#print(level[np.amin(indices_seven_final)], level[np.amax(indices_seven_final)])
    
plt.figure()
plt.xlabel('Radius(kpc)');
plt.ylabel('PDF');

x = np.amin(level.v);

plt.hist(radius[indices_three_final].v/(1000*one_parsec), label = '{}'.format(x), density = True, alpha = 0.5);
if(len(index_5[0])>0):
	plt.hist(radius[indices_five_final].v/(1000*one_parsec), label = '{}'.format(x+1), density = True, alpha = 0.5);
if(len(index_7[0])>0):
	plt.hist(radius[indices_seven_final].v/(1000*one_parsec), label = '{}'.format(x+2), density = True, alpha = 0.5);
if(len(index_8[0])>0): 
	plt.hist(radius[indices_eight_final].v/(1000*one_parsec), label = '{}'.format(x+3), density = True, alpha = 0.5);
if(len(index_9[0])>0): 
	plt.hist(radius[indices_nine_final].v/(1000*one_parsec), label = '{}'.format(x+4), density = True, alpha = 0.5);
plt.legend();

plt.savefig('refine_info_{}.png'.format(plt_number))


