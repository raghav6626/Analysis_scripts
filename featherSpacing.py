import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sc
import yt
from mpi4py import MPI
from yt.units import kpc, pc
from scipy.fftpack import fft, fftfreq
import math
import argparse

parser = argparse.ArgumentParser();


#=======================================
#Defining constants
#======================================

one_kpc = 3.08567E+21
avogadro_no = 6.023e23
mean_mass = 1.3;

#========================================
#arguments 
#=======================================


scaleHeight = 100*pc;

parser.add_argument('--plt_number', type = int, nargs = 1, required = True , help = 'the plot file number');
parser.add_argument('--densityCut', type = float, nargs = 1, required = True , help = 'the threshold density of the particles to plot (in cm^-3)');
parser.add_argument('--dphi', type = float, nargs = 1, required = True , help = 'the dphi for binning in constant phi bins');
parser.add_argument('--weighted', type = bool, nargs = 1, required = False , help = 'Whether to have mass-weighted particle density values or not.');
parser.add_argument('--radius', type = float, nargs = 1, required = True, help = 'The radius at which to take the fourier transform (in kpc)');

args = parser.parse_args();

directory = './'

densityCutOff = args.densityCut[0];  #particles/cm3

spacing = args.dphi[0]; #radians


rCutOff = (args.radius[0])*kpc #in kpc
#=================================================================================
if(args.plt_number[0]>=100):
    fileName = directory + 'Disc_hdf5_plt_cnt_0{}'.format(args.plt_number[0]);
elif(args.plt_number[0]>=10):
    fileName = directory + 'Disc_hdf5_plt_cnt_00{}'.format(args.plt_number[0]);
else:
    fileName = directory + 'Disc_hdf5_plt_cnt_000{}'.format(args.plt_number[0]);

print('Reading file ', fileName)
ds = yt.load(fileName);
ad = ds.all_data();

if(args.weighted):

    
    if(args.weighted[0] == True):
        print('Weighting the bins by cell_mass');
        mass = ad['cell_mass'];


else:
    mass = np.zeros(len(ad['cell_mass']));
    mass[:,] = 1;       

density = ad['density'];
particle_density = (density.v*avogadro_no)/mean_mass;
x = ad['x'];
y = ad['y'];
z = ad['z'];
r = np.sqrt(x**2 + y**2);

#=================================================================================


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



#==========================================
#Getting the 
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

#sorting the theta and the particle density arrays
sorted_indices = np.argsort(theta);
theta_sorted = theta[sorted_indices];
particle_density_sorted = particle_density[sorted_indices];
mass_sorted = mass[sorted_indices]

#creating bins of constant size in the theta direction for the fourier transform
(thetaBins, particleDensity_binned) = createBins(spacing, theta_sorted, particle_density_sorted, mass_sorted);




particleDensity_binned = np.asarray(particleDensity_binned);
thetaBins = np.asarray(thetaBins);



from scipy.signal import find_peaks

peaks, _ = find_peaks(particleDensity_binned, height=1)
print(peaks)

#checking if the density is smooth enough visually 
plt.figure( figsize = (12,8));
plt.plot(thetaBins, particleDensity_binned);
print('Plotting the individual peaks now')
plt.plot(thetaBins[peaks], particleDensity_binned[peaks], 'x');
plt.xlabel('$\phi$ (radians, r = 7kpc)');
plt.ylabel('Number Density (cm$^{-3}$)');
plt.savefig('number_density_theta_{}_densityCut_{}_radius_{}.png'.format(args.plt_number[0], args.densityCut[0], math.trunc(10*args.radius[0])));


delta_x, delta_theta = getSpacing(thetaBins, peaks);


fig, ax = plt.subplots();
ax.hist(delta_x, density = True);
ax.set_xlabel('Spacing (kpc)');
ax.set_ylabel('P (spacing)');
ax.minorticks_on();

plt.savefig('spacing_{}_radius_{}_withoutCut.png'.format(args.plt_number[0],math.trunc(10*args.radius[0])));





#======================================================
#getting the position of the spiral arms and binning again
(theta_left_edge, theta_right_edge) = get_spiralEdges(thetaBins, particleDensity_binned)
thetaBins = thetaBins[theta_left_edge:theta_right_edge];
particleDensity_binned = particleDensity_binned[theta_left_edge:theta_right_edge];    

peaks, _ = find_peaks(particleDensity_binned, height=1)

#finding the mean spacing of the peaks and their spreads (maybe)
#saving the distribution of the deltaas 

delta_x, delta_theta = getSpacing(thetaBins, peaks);

fig, ax = plt.subplots();
ax.hist(delta_x, density = True);
ax.set_xlabel('Spacing (kpc)');
ax.set_ylabel('P (spacing)');
ax.minorticks_on();

plt.savefig('spacing_{}_radius_{}.png'.format(args.plt_number[0],math.trunc(10*args.radius[0])));




#===================================================================
#2. Finally, taking the fourier transform of the thing. 
#===================================================================
#since we are taking the region between the spiral arms. It only makes sense to include the 
#k values that have length scales less than our sample theta values...
#we make an explicit cut before plotting - ask Robi if this is okay?
'''



N = len(thetaBins);
T = thetaBins[2] - thetaBins[1];
x = thetaBins;
y = particleDensity_binned - np.average(particleDensity_binned);

yf = fft(y);
xf = fftfreq(N, T)[:N//2]
yf = yf[:N//2];


#taking all the values xf>0
ind = np.where(xf>0.1)[0];
xf = xf[ind];
yf = yf[ind];


#slicing so that the periodicity is less than the sample size.
lmbda = 1/xf;
ind = np.where(lmbda<np.abs(thetaBins[-1] - thetaBins[0]))[0];
lmbda = lmbda[ind];
print(np.amax(lmbda), np.amin(lmbda));
print(thetaBins[-1] - thetaBins[0]);
yf = yf[ind]
x = lmbda*rCutOff.v;


#finally plotting

fig, ax = plt.subplots();
ax.plot(x, 2.0/N * np.abs(yf));
ax.minorticks_on();
ax.set_xscale('log')
#ax.set_xlabel('$ k \:(1/\lambda)$')
ax.set_xlabel('x (kpc)')
ax.set_ylabel('$\Sigma e^{-2 \pi j k x} \mathrm{n(x)}$ ')
ax.set_xlim(np.amax(x), np.amin(x));
plt.savefig('fourier_density_{}_densityCut_{}_radius_{}_distance'.format(args.plt_number[0], math.trunc(args.densityCut[0]), math.trunc(10*args.radius[0]) ));

#theta_graph
x = x/rCutOff.v


fig, ax = plt.subplots();
ax.plot(x, 2.0/N * np.abs(yf));
ax.minorticks_on();
#ax.set_xlabel('$ k \:(1/\lambda)$')
ax.set_xlabel('theta (rad)')
ax.set_ylabel('$\Sigma e^{-2 \pi j k x} \mathrm{n(x)}$ ')
ax.set_xlim(np.amax(x), np.amin(x));
plt.savefig('fourier_density_{}_densityCut_{}_radius_{}_theta'.format(args.plt_number[0], math.trunc(args.densityCut[0]), math.trunc(10*args.radius[0]) ));

'''

