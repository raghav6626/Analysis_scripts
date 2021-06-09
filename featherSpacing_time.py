import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sc
import yt
from mpi4py import MPI
from yt.units import kpc, pc
from scipy.fftpack import fft, fftfreq
import math
import argparse
from scipy.signal import find_peaks

parser = argparse.ArgumentParser();


def binInRadius(rBin, delta_r, radius):
    
#     dr = np.subtract(radius. rBin);
    dr = radius - rBin;
    
    indices = np.where(np.absolute(dr)<delta_r)
    
    return indices[0];

def binInZ(delta_z, z):
    
    indices = np.where(np.absolute(z)<delta_z);

    return indices[0];

def binDensity(densityCutOff, numberDensity):
    
    indices = np.where(numberDensity>=densityCutOff);
    
    return indices[0];



#takes in the sortedTheta, theta_spacing, quntity to bun
#returns the center of the theta_bins and the quantity in each bin.
def createBins(dtheta, theta_array, quantity, mass):
    theta_left = theta_array[0];
    binned_quantity = [];
    theta_bins_middle = [];
    while theta_left<=theta_array[-1]:
        theta_right = theta_left+dtheta;
        indices = np.where((theta_array>=theta_left) & (theta_array<theta_right))[0];
        if(len(indices)==0):
            binned_quantity.append(0);
        else:
            binned_quantity.append(np.average(quantity[indices], weights = mass[indices]));
        theta_bins_middle.append((theta_left + theta_right)/2);
        theta_left = theta_right
        
    return theta_bins_middle,binned_quantity;


#returns the indices of the edge of the spiral arms in the sorted thetaBins array
def get_spiralEdges(thetaBins, particleDensity_binned):
    index = np.where(particleDensity_binned == 0)[0];
    count = 1;
    countArray = [];
    for i in range(len(index) -1):
        if(index[i+1] - index[i] == 1):
            count = count + 1;
        else:
            countArray.append(count);
            count = 1;
    countArray.append(count);
#    countArray = np.asarray(countArray)    
    print('Count array is ', countArray);
    print('Index array is ', index);
#1 Finding the largest and the second largest counts, I have also added an exception for when both are equal. 

    largest_id = np.where(countArray == np.amax(countArray))[0];
    if(len(largest_id)>1):
        l = largest_id[0];
        slargest_id = [largest_id[1]];
        largest_id = [];
        largest_id.append(l);

    else:
        x = countArray.copy();
        countArray.remove(np.amax(countArray));
        secondLargest = np.amax(countArray);
    
        slargest_id = np.where(x == secondLargest)[0];

        countArray = x

    largest_count = countArray[largest_id[0]];
    slargest_count = countArray[slargest_id[0]];   
    print('largest and second larges counts are - ', largest_count, slargest_count);
    
#Now, we just go back to our original index array and get back the values 
#The left spiral arm's, right edge is taken
#The right spiral arm's left edge is taken

#a. Taking the left spiral arm
    if(largest_id[0]<slargest_id[0]):
        sp_left_id = largest_id[0];
        sp_right_id = slargest_id[0]
    else:
        sp_left_id = slargest_id[0];
        sp_right_id = largest_id[0];
    

    print('leftEdge Id, rightEdge Id,', sp_left_id, sp_right_id);
#right edge of the left spiral arm
    count = 0;
    if(sp_left_id>0):
        for i in range(sp_left_id+1):
            count = count + countArray[i];
    else:
        count = countArray[sp_left_id];
    
# print(index[count-1]);
    theta_left_edge_index = index[count-1];
    print('left edge theta = ', thetaBins[theta_left_edge_index], particleDensity_binned[theta_left_edge_index])

#left edge of the right spiral arm 
    count = 0;

    for i in range(sp_right_id):
        count = count + countArray[i];

# print(index[count])
    theta_right_edge_index = index[count];
    print('right edge theta = ', thetaBins[theta_right_edge_index], particleDensity_binned[theta_right_edge_index])
    
    return theta_left_edge_index, theta_right_edge_index;

#==========================================
#Getting the spacings
#We remove the largest two distances because they 
#correspond to the spiral arm locations.
#===============================================
def getSpacing(thetaBins, peaks):

    delta_theta = [];

    for i in range(len(peaks)-1):
        delta_theta.append(thetaBins[peaks[i+1]] - thetaBins[peaks[i]]);
   
    delta_theta.remove(np.amax(delta_theta));
    delta_theta.remove(np.amax(delta_theta));
    delta_theta = np.asarray(delta_theta);
    
    delta_x = delta_theta*rCutOff.v; #in kpc

    return  delta_x, delta_theta;


#===========================================
def plotFourier(ad):

#Gathering the data into numpy arrays

    density = ad['density'];
    particle_density = (density.v*avogadro_no)/mean_mass;
    x = ad['x'];
    y = ad['y'];
    z = ad['z'];
    r = np.sqrt(x**2 + y**2);
    mass = np.zeros(len(ad['cell_mass']));
    mass[:,] = 1;   



#Placing the various cuts - radius, z, min_density

    indices_gas = binInRadius(rCutOff, 100*pc, r)



    x = x[indices_gas];
    y = y[indices_gas];
    z = z[indices_gas];
    particle_density = particle_density[indices_gas];
    r = r[indices_gas]


    indices_gas_z = binInZ(scaleHeight, z)
    indices_gas = indices_gas_z;

    x = x[indices_gas];
    y = y[indices_gas];
    z = z[indices_gas];
    r = r[indices_gas];
    particle_density = particle_density[indices_gas];


#density cutOff of one-particle/cm^3??
    indices_gas_z_density = binDensity(densityCutOff, particle_density);
    indices_gas = indices_gas_z_density;
    x = x[indices_gas];
    y = y[indices_gas];
    z = z[indices_gas];
    r = r[indices_gas];
    particle_density = particle_density[indices_gas];



    theta = np.arctan2(y.v,x.v);#returns the angle back in radian


#==================================================================
#Getting the fourier transform of the density
#1. Creating bins of equal spacing
#===================================================================


#sorting the theta and the particle density arrays
    sorted_indices = np.argsort(theta);
    theta_sorted = theta[sorted_indices];
    particle_density_sorted = particle_density[sorted_indices];
    mass_sorted = mass[sorted_indices]

#bins of constant size in theta
    (thetaBins, particleDensity_binned) = createBins(spacing, theta_sorted, particle_density_sorted, mass_sorted);

    thetaBins = np.asarray(thetaBins);
    particleDensity_binned = np.asarray(particleDensity_binned);
    print('calculating peaks')   
    peaks, _ = find_peaks(particleDensity_binned, height=1)

    delta_x, delta_theta = getSpacing(thetaBins, peaks);

    return (np.percentile(delta_x, 16), np.percentile(delta_x, 50), np.percentile(delta_x, 64));    

#getting the position of the spiral arms and binning again 
#uncomment this if you want to use the binning 
'''
(theta_left_edge, theta_right_edge) = get_spiralEdges(thetaBins, particleDensity_binned)
thetaBins = thetaBins[theta_left_edge:theta_right_edge];
particleDensity_binned = particleDensity_binned[theta_left_edge:theta_right_edge]; 
'''


#===========================================
#For the spacings 
#============================================
def plotSpacings(dirr, label, ax, x_offset, dx_offset):
    percentiles = [];
    time_array = [];
    for i in args.plt:
        if (i<10):
            ds= yt.load('{}Disc_hdf5_plt_cnt_000{}'.format(dirr, i));
        elif(i<100): 
            ds= yt.load('{}Disc_hdf5_plt_cnt_00{}'.format(dirr, i));
        else:
            ds= yt.load('{}Disc_hdf5_plt_cnt_0{}'.format(dirr, i))
        ad = ds.all_data();
        time_array.append(ds.current_time.v/one_myr + x_offset);
        (sixteenth, fiftieth, sixty_fourth) = plotFourier(ad);
        percentiles.append([sixteenth, fiftieth, sixty_fourth]);
    
    percentiles = np.asarray(percentiles); 
    ax.scatter(time_array, percentiles[:,1], label = label, s = 5);
    ax.errorbar(time_array, percentiles[:,1], yerr = [ percentiles[:,2] - percentiles[:,1], percentiles[:,1] - percentiles[:,0] ], capsize=3.0, elinewidth = 0.5 )
    x_offset = x_offset + dx_offset;
    return;

#=======================================
#Defining constants
#======================================

one_kpc = 3.08567E+21
avogadro_no = 6.023e23
mean_mass = 1.3;
one_myr = 1e6 * (3600*24*365)
#========================================
#arguments 
#=======================================


scaleHeight = 100*pc;


parser.add_argument('--plt', type = int, nargs = '+', required = True , help = 'the plot file numbers to read, separated by space');
parser.add_argument('--densityCut', type = float, nargs = 1, required = True , help = 'the threshold density of the particles to plot (in cm^-3)');
parser.add_argument('--dphi', type = float, nargs = 1, required = True , help = 'the dphi for binning in constant phi bins');
parser.add_argument('--weighted', type = bool, nargs = 1, required = False , help = 'Whether to have mass-weighted particle density values or not.');
parser.add_argument('--radius', type = float, nargs = 1, required = True, help = 'The radius at which to take the fourier transform (in kpc)');

args = parser.parse_args();

directory = './'

densityCutOff = args.densityCut[0];  #particles/cm3

spacing = args.dphi[0]; #radians

rCutOff = (args.radius[0])*kpc #in kpc

omega = (2e7)/(rCutOff.v * one_kpc); #rad/sec


#so that the points are offset from each otherr!
x_offset = -1;
dx_offset = 1;

# directories = [ '/scratch/ek9/ra8040/flash_sims/sims/hollow_n2/sp_08/v_rms_10/highRes/', '/scratch/ek9/ra8040/flash_sims/sims/hollow_n2/sp_12/highRes/', '/scratch/ek9/ra8040/flash_sims/sims/mhd_runs/beta_100/sp_08/', '/scratch/ek9/ra8040/flash_sims/sims/mhd_runs/beta_100/sp_12/highRes/' ];

directories = ['/scratch/ek9/ra8040/flash_sims/sims/highRes_test/100Myr/', '/scratch/ek9/ra8040/flash_sims/sims/medium_fromStart/', '/scratch/ek9/ra8040/flash_sims/sims/hollow_n2/sp_08/v_rms_10/highRes/'];

# labels = ['F_25_hydro','F_35_hydro','F_25_beta_100','F_35_beta_100'];
labels = ['F_25_hydro_high_res', 'F_25_hydro_medium_res', 'F_25_hydro_low_res']
#=================================================================================

ft_array = [];

fig, ax = plt.subplots();

time_array = [];
percentiles = []

for dirr, label in zip(directories, labels):
    print(dirr, label);
    plotSpacings(dirr, label, ax, x_offset, dx_offset);



ax.set_ylabel('Spacing (kpc)');
ax.set_xlabel('Time (Myrs)');
ax.minorticks_on();
ax.set_ylim(0, )
plt.legend();
plt.savefig('spacing_{}_{}_radius_{}_multiRun.png'.format(args.plt[0], args.plt[-1], math.trunc(10*args.radius[0])));