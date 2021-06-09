from mpi4py import MPI
import time
import yt 
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import math
from numpy import linalg as LA



def get_center_of_mass(x, m):
    
    cx = np.multiply(x, m);
    cx = np.sum(cx);
    M = np.sum(m) 

    return cx/M;


def get_distance_in_plane(leaf_clump):

    cx = get_center_of_mass(leaf_clump['grid', 'x'], leaf_clump['grid', 'cell_mass']);
    cy = get_center_of_mass(leaf_clump['grid', 'y'], leaf_clump['grid', 'cell_mass']);


    return np.sqrt(cx**(2) + cy**(2));

#=======================================
#Takes in a lead clumps and returns the COM velocity
#========================================

def get_center_of_mass_velocity(leaf_clump):
    
    mass = leaf_clump['grid', 'cell_mass'].v
    v = np.zeros((3, len(mass)));
#     print(v[0])
    v[0] = np.multiply(leaf_clump['grid', 'velx'].v, mass)
    v[1] = np.multiply(leaf_clump['grid', 'vely'].v, mass)
    v[2] = np.multiply(leaf_clump['grid', 'velz'].v, mass)
    v = np.transpose(v);
#     print(v);
    return np.sum(v,0)/np.sum(mass);

def get_com(r, m):
    r = np.transpose(r);
    
    rcmx = get_center_of_mass(r[0], m)
    rcmy = get_center_of_mass(r[1], m)
    rcmz = get_center_of_mass(r[2], m)

    rcm = [rcmx, rcmy, rcmz];
    rcm = np.array(rcm);
    
    return rcm;

def get_angular_momentum(leaf_clump, v_com):

    m = leaf_clump['grid', 'cell_mass'];
    r = np.zeros((3, len(m)));
    v = np.zeros((3, len(m)));

    v[0] = leaf_clump['grid', 'velx'].v;
    v[1] = leaf_clump['grid', 'vely'].v;
    v[2] = leaf_clump['grid', 'velz'].v;
    v = np.transpose(v);
    
    r[0] = leaf_clump['grid', 'x'].v;
    r[1] = leaf_clump['grid', 'y'].v;
    r[2] = leaf_clump['grid', 'z'].v;
    r = np.transpose(r);
    
    r_com = get_com(r, m); 
    
    r_rel = r-r_com;    

    v_rel = v - v_com;

    r_cross_v =  np.cross(r_rel, v_rel);
    
    r_cross_v = np.transpose(r_cross_v);
    mvr = np.zeros((3, len(m)));
    mvr[0] = np.multiply( r_cross_v[0], m);
    mvr[1] = np.multiply( r_cross_v[1], m);
    mvr[2] = np.multiply( r_cross_v[2], m);
    mvr = np.transpose(mvr);
    
    
    return np.sum(mvr,0);


def get_moment_of_Inertia(leaf_clump):
    m = leaf_clump['grid', 'cell_mass'];
    
    r = np.zeros((3, len(m)));
    r[0] = leaf_clump['grid', 'x'].v;
    r[1] = leaf_clump['grid', 'y'].v;
    r[2] = leaf_clump['grid', 'z'].v;
    r = np.transpose(r);
    
    r_com = get_com(r, m);

    r_rel = r-r_com;

    
    r_square = LA.norm(r_rel, axis =1) **2;
    r_square = np.multiply(m, r_square)
    
    return (np.sum(r_square).v, r_rel);


#takes in a single leaf_clump
#returns the velocity dispersion of a single clump
def velocity_dispersion(leaf_clump):

    v_cm = get_center_of_mass_velocity(leaf_clump); #a vector
    st = time.time();
    L = get_angular_momentum(leaf_clump, v_cm);	#a vector again
    print('angular momentum takes ', time.time() - st);
    st = time.time();
    (I, r_rel) = get_moment_of_Inertia(leaf_clump);	#a scalar (not realy but okay)
 
    omega = L/I;

    v_rot = np.cross(omega, r_rel); #v = omega x r_rel 

    
    #initializing the v 
    v = np.zeros((3, len(leaf_clump['grid', 'cell_mass'])));
    v[0] = leaf_clump['grid', 'velx'].v;
    v[1] = leaf_clump['grid', 'vely'].v;
    v[2] = leaf_clump['grid', 'velz'].v;
    v = np.transpose(v);
    
    v = v - v_cm - v_rot; # cm/sec  
    v = v*10**(-5)        #km/sec #these are n vectors right now. 

    v_squared = LA.norm(v, axis =1) **2;

    v_rms = np.sqrt(np.average(v_squared)); #v_rms = sqrt( mean (v^(2)))

#     print('moment of Inertia takes', time.time() - st);
#     print('Angular momenta is, ', L);
#     print('Moment of Inertia is', I)
#     print('vcm is', v_cm)
#     print('Omega is ', omega);
#     print('v rot is ', v_rot);
#     print('v squared', v_squared);
#     print('velocity dispersion is , ', v_rms);

    return (v_rms);



def velocity_dispersion_histogram(leaf_clumps, plot = False):
    v_d = np.zeros(len(leaf_clumps));
    

    for i in range(len(v_d)):
        v_d[i] = velocity_dispersion(leaf_clumps[i]);
    
    if(plot):
        plt.figure();
        plt.style.use('classic');
        plt.xlabel('Velocity Dispersion (km/sec)');
        plt.ylabel('Frequency');
        plt.hist(v_d, edgecolor = 'black');
        plt.savefig('{}{}_{}_velocity_dispersion_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));
        
    return v_d;

def select_clump(clump):
    print('selecting clumps')
    if(get_distance_in_plane(clump)>inner_disc_radius*1000*one_parsec and len(clump['grid', 'cell_mass'])>cell_threshold):
        return True;

    return False;


m_0 = 1.989e+33 #in grams
one_parsec = 3.086e+18 #in cms
inner_disc_radius = 2 #in kpc
cell_threshold = 50 #The minimum number of cells that should be there in the clump.



#=======================================
#Reading the parsed arguments.
#=======================================

parser = argparse.ArgumentParser();

 #density_parameter in particles/cm^(3)


parser.add_argument('-i', type = str, nargs = 1, required = True , help = 'the input hdf5 file directory');#the input hdf5 file

parser.add_argument('-o', type = str, nargs = 1, help = 'output directory name. If not provided, the current directory is taken');#the input hdf5 file

parser.add_argument('--cell_threshold', type = int, nargs = 1, help = 'The min. number of cells in a clump...')

parser.add_argument('--clumpName', required = True ,type = str, nargs = 1, help =
'The name of the clump to read.')


args = parser.parse_args();


if(args.cell_threshold != None):
	cell_threshold = args.cell_threshold[0];


if(args.o == None):
	output_dir = './';
else:
	output_dir = args.o[0];



clumpName = args.clumpName[0];




#========================================
# Loading the clumps and the plot file 
#========================================

#print(plotFileName);
print(clumpName);

#ds = yt.load(plotFileName); 
#data = ds.all_data();
ds_clumps = yt.load(args.i[0]+clumpName);
leaf_clumps = ds_clumps.leaves


#=============================================
#Calculating the v_d using MPI parallel
#=============================================

# get MPI ranks
comm = MPI.COMM_WORLD
size = comm.Get_size()
my_rank = comm.Get_rank();

#comm.Barrier synchronises all the threads/processes. 
#All the threads have to reach this point of execution, only then the program will move forward
comm.Barrier()

n = len(leaf_clumps);

count = np.zeros(size)

#dividing the work into chunks
my_start = my_rank*(n//size);     #// is integer division, that is, divide folowed by floor
my_end = (my_rank+1)*(n//size) -1;

#giving the rest over to the last chunk 
if(my_rank == size - 1):
    my_end = n-1;



if(my_rank == 0):
	print('Each processor has about {} clumps to analyze'.format(my_end-my_start));
	print('Total number of clumps are - ', len(leaf_clumps));

f = open('{}velocity_dispersion_parallel_{}.txt'.format(output_dir, clumpName), 'a')
	

i = my_start;


while(i<=my_end):
	to_select = select_clump(leaf_clumps[i]);
	
	if(to_select):
		count[my_rank] = count[my_rank]+1
		#Updating the value of count across all the processors using bcast.
		print('Procesor {} calculating for the clump number {} that has {} cells'.format(my_rank, i, len(leaf_clumps[i]['grid','cell_mass'])));
		
		st = time.time();
		v_d  = velocity_dispersion(leaf_clumps[i]);
		f.write('size ={} \t v_d = {} \t ID: {} \n'.format(len(leaf_clumps[i]['grid', 'cell_mass']), v_d, i));
		print('Time taken for clump number {} is {} secs'.format(i, time.time() - st));		
			
	else:
		print('Clump number {} is not selected'.format(i));

	i = i+1;


comm.Barrier();
print('Total number of clumps selected are - {}'.count);
























