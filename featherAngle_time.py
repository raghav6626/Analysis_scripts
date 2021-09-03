import h5py
import yt
import numpy as np 
import matplotlib.pyplot as plt 
import math
import random
import time
from mpi4py import MPI
from scipy.signal import find_peaks
import pickle
import h5py as h5

#Transforming to (R,Theta) from (x,y)
one_kpc = 3.08567E+21
avogadro_no = 6.023e23;
mean_mass = 1.3
density_cut = 10 #particle/cm3
one_myr = 60*60*24*365*1e6;

def get_radial_coordinates(x,y):
	r = np.sqrt(x**2 + y**2);
	theta = np.arctan2(y,x)
	return r, theta;



def getX_limits_new(edge1, edge2, yMax, pitch, y):
	pitch = pitch*np.pi/180;
	m = np.tan(pitch);
	
	
	x1 = (y-edge1[1])/m + edge1[0];
	x2 = x1 + (edge2[0] - edge1[0]);
	
	
	return (x1, x2);

#make an array y from (yMin, yMax)
#another array x from [(x1,x2)] for every value of y in (yMin, yMax)
#take a point from the grid
#check if it is in the parallelogram
#append the index if it os
def getPoints_new(edge1, edge2, yMax, pitch, y_grid, x_grid):
	
	indices = np.where( (y_grid>= edge1[1]) & (y_grid <= yMax) )[0]
   
	y_grid = y_grid[indices];
	x_grid = x_grid[indices];
	
	
	
	indices_parallelogram = []
	#for evert x and y, you have to check whether they lie within the parallelogram
	for i in range(len(y_grid)):
		x1, x2 = getX_limits_new(edge1, edge2, yMax, pitch, y_grid[i]);
		if(x_grid[i] > x1 and x_grid[i] < x2):
			indices_parallelogram.append(i);
	
	return indices[indices_parallelogram];



def getFilter_radius(pitch):
	
	P_degree = pitch #degrees
	P = P_degree *(np.pi/180);
	c = 0; #The centre (0,0)
	ml = -1/np.tan(P);	
	
#Step 1: Find the radius of the circle needed for the range of cosine
	r = (np.pi/2)*np.abs((ml/np.sqrt(ml**2 + 1)));
	
	return r;

def getF(pitch, x_grid_norm, y_grid_norm, center):
	P_degree = pitch #degrees
	P = P_degree *(np.pi/180);
	c = 0; #The centre (0,0)
	ml = -1/np.tan(P); 
	F = np.cos((1/ml) * (x_grid_norm - center[0]) + (y_grid_norm - center[1]) + c);
	return F;

#===============================================================
#We make a filter/kernel in unitless coordinates
#Then we scale the filter with the diameter of the circular region
#Then, for each point in the circular window -> we tranform the coordinates to the unitless ones and
#calculate the value of F for each pixel in the circular window.
#=================================================================
def getFilter_wrapper(theta_grid_window, lnR_grid_window, pitch, center, diameter):

	r_filter = getFilter_radius(pitch);
	ratio = (2*r_filter)/diameter
	F_array = np.zeros(len(theta_grid_window));	
	
   
	for i in range(len(theta_grid_window)):
		F_array[i] = getF(pitch, theta_grid_window[i]*ratio, lnR_grid_window[i]*ratio, center*ratio); 
	
	return F_array;


#=================================================================================
#create a circular window of a given diameter dW and a centre - C(lnR_c, theta_c)
#takes in - 
#y_grid: lnR grid
#x_grid: theta in radians grid
#returns: the indices of the cirlce centered at c, with a diameter = diameter. 
#===============================================================================
def getCircular_window(y_grid, x_grid, c, diameter):
	
	r = np.sqrt( (x_grid - c[0])**2 + (y_grid-c[1])**2);
	
	indices = np.where(r<=diameter/2)[0];
		
	return indices;

#=============================================================================
#we take in the 
#1. Filter Function
#2. the image - density values in a circle with a physical diameter
#we just do a multiplication for getting a cross-correlation
#============================================================================
def crossCorrelate(f, I):
	f = np.asarray(f);
	I = np.asarray(I);
	return (f*I);


def plotWindow_common(indices_circularWindow, theta_grid, lnR_grid, I):
	xLabel = 'theta'
	yLabel = 'lnR'
	#we just plot and check the filter and the data that we have
	fig, ax = plt.subplots();
	fig.set_size_inches(10,8);
	plt.style.use('classic')

	im = ax.scatter(theta_grid * 180/np.pi, lnR_grid, c = np.log10(I), edgecolor = 'face', cmap = 'viridis');
	ax.scatter((theta_grid[indices_circularWindow]) * 180/np.pi, lnR_grid[indices_circularWindow], c = 'black', alpha = 0.7, s = 0.05, edgecolor = 'face', cmap = 'viridis');   
	clb = plt.colorbar(im);
	ax.set_xlabel(xLabel);
	ax.set_ylim(lnr_min, lnr_max);
	ax.set_xlim(-180, 180)
	ax.set_ylabel(yLabel)
#	 plt.savefig('circularWindow_lnR_theta')  


	return;

def plotWindow_individual(indices_circularWindow, theta_grid, lnR_grid, I, diameter, cloud_id):
	xLabel = 'theta'
	yLabel = 'lnR'
	#we just plot and check the filter and the data that we have
	fig, ax = plt.subplots();
	fig.set_size_inches(10,8);
	plt.style.use('classic')

	im = ax.scatter(theta_grid[indices_circularWindow] * 180/np.pi, lnR_grid[indices_circularWindow], c = np.log10(I[indices_circularWindow]), edgecolor = 'face', cmap = 'viridis');
#	 ax.scatter((theta_grid[indices_circularWindow]) * 180/np.pi, lnR_grid[indices_circularWindow], c = 'black', alpha = 0.7, s = 0.05,edgecolor = 'face', cmap = 'viridis');   
	clb = plt.colorbar(im);
	ax.set_xlabel(xLabel);
#	 ax.set_ylim(lnr_min, lnr_max);
#	 ax.set_xlim(-180, 180)
	ax.set_ylabel(yLabel)	
	plt.savefig('window_{}_{}'.format(math.trunc(diameter*100), cloud_id));
	plt.close()
	
	
	return

#=======================================================
#returns the CCsum for a particular (pitch, diameter, center)
#Assuming you have a uniform grid. 
# I = density or the IR intensity of an image. 
# (x_grid, y_grid):  the (theta (radians), lnR)
# centre = center of the window.
# diameter: the diameter in lnR units
# pitch: pitch of the filter in degrees
#output: see eq (5) Puerari et. al. 2014.
#Modified - Now we divide by the perfect correlation for normalisation.
#==========================================================
def crossCorrelation_sum(I, x_grid, y_grid, center, diameter, pitch, toPlot, cloud_id):

	# Step 1 : Make a circular window in the data
	indices_circularWindow = getCircular_window(y_grid, x_grid, center, diameter);

	if(toPlot):
		plotWindow_individual(indices_circularWindow, x_grid,  y_grid, I, diameter, cloud_id);
	
	I = I[indices_circularWindow];
	
	#Step 2: Make a filter of the same resolution as that of the density window above...
	#we normalise our filter so that sum(F) = 1
	F = getFilter_wrapper(np.asarray(x_grid[indices_circularWindow]),np.asarray(y_grid[indices_circularWindow]), pitch, center, diameter);
	if(np.amin(F)<0):
		print('the filter is negative at :{}, having pitch = {}, with a value: {} '.format(center, pitch, np.amin(F)));
		plt.figure();
		im = plt.scatter(theta_grid[indices_circularWindow], lnR_grid[indices_circularWindow], c = F, cmap = 'viridis', edgecolor = 'face');
		plt.colorbar(im);
		
	F = F/np.sum(F);

	
	#Step 3: Cross-correlate. We normalize the cc with the np.sum(I)
	#As an additional step, we divide by the perfect cross-correlation at 
	#the point - when we cross correlate the filter with itself.
	CC_perfect = crossCorrelate(F, F);
	CC_perfect = np.sum(CC_perfect);


	#calculates the CC of the F with the I
	CC = crossCorrelate(F, I);
	CC = np.sum(CC)/np.sum(I);

	return CC/CC_perfect



def plot(indices_parallelo, indices_window):

	xLabel = 'theta (degree)'
	yLabel = 'lnR'
	fig, ax = plt.subplots();
	fig.set_size_inches(10,8);
	plt.style.use('classic')
# ax.plot(pitch_array, ccSum);
	ax.scatter((theta_grid[indices_parallelo])*180/np.pi, lnR_grid[indices_parallelo], c = np.log10(density[indices_parallelo]), edgecolor = 'face', cmap = 'viridis');
	ax.scatter((theta_grid[indices_window])*180/np.pi, lnR_grid[indices_window], s = 0.1, alpha = 1);
	ax.minorticks_on();
	ax.set_xlabel(xLabel);
	ax.set_ylabel(yLabel)
	return;


#==========================================
#The wrapper function for CCsum for a pitch_array, centre, diameter
#=========================================

def getCCSum_lambda(density, theta_grid, lnR_grid, centre, diameter, pitch_array, cloud_id):
	toPlot = False;		
	ccSum = np.zeros(len(pitch_array));
	
	#for each pitch angle, we get one number 
	for i in range(len(pitch_array)):
		ccSum[i]= crossCorrelation_sum(density, theta_grid, lnR_grid, centre, diameter, pitch_array[i], toPlot, cloud_id)
		toPlot = False;
 
	#an array of size == len(pitch angle);
	if(np.amin(ccSum)<0):
		print('Found a negative ccSum at: ', centre);
	return ccSum;



def plotCCSum(ccSum, pitch_array, cloud_id):
	fig, ax = plt.subplots();
	plt.style.use('default')
	plt.xlabel('Degrees');
	plt.ylabel('Correlation Sum (normalised)');
	ax.set_xlim(-90, 90)
	ax.plot(pitch_array,ccSum);
	plt.savefig('ccSum_onePoint_{}'.format(cloud_id));

	plt.close();

	return

#===========================================
#We calculate the height of the peak in case our data 
#is none other than a filter itself
#==============================================

def crossCorrelation_perfectHeight(x_grid, y_grid, center, diameter, pitch_array):
	
	#Step 1: Get a random filter of the same res as that of the run
	indices_circularWindow = getCircular_window(y_grid, x_grid, center, diameter);

	I = getFilter_wrapper(np.asarray(x_grid[indices_circularWindow]), np.asarray(y_grid[indices_circularWindow]), 10, center, diameter);
	I = I/np.sum(I);

	#Step 2 : Make a correlation
	ccSum_array = [];
	for i in range(len(pitch_array)):
		F = getFilter_wrapper(np.asarray(x_grid[indices_circularWindow]),np.asarray(y_grid[indices_circularWindow]), pitch_array[i], center, diameter);
		F = F/np.sum(F);

		ccSum_perfect = crossCorrelate(F, F);
		ccSum_perfect = np.sum(ccSum_perfect);

		ccSum = crossCorrelate(F, I)
		ccSum = np.sum(ccSum)/np.sum(I)

		ccSum_array.append(ccSum/ccSum_perfect);


	# plt.figure();
	# plt.style.use('classic');
	# plt.plot(pitch_array, ccSum_array);
	# plt.savefig('ccSum_perfect.png'.format(time.time()))


	return np.amax(ccSum_array) - np.amin(ccSum_array);




#===============================================================================
#A special version of the functions above that only calculate the sum 
# of the given centers
#Returns the peaks of respective clouds, with the respective cloud_ids of the 
#corresponding clouds too
#==============================================================================
def convolve_centers(centers, diameter, theta_grid, lnR_grid, density, cloud_ids):
	ccSum_array = [];
	peak_array = []
	useful_cloud_ids = [];
	height_array = [];
	max_array = [];

	pitch_array = np.linspace(0.1, 89.9);
	i = 0
	start_time = time.time();

	#for each center, we find a peak
	while i<len(centers):
		# print(i);
		centre = centers[i];


		#Step 1: Get the ccSum/ccSum perfect 

		ccSum = getCCSum_lambda(density, theta_grid, lnR_grid, centre, diameter, pitch_array, cloud_ids[i])
		ccSum_array.append(ccSum);


		#Step 2: Get the perfect height for this point and diameter 

		perfect_height = crossCorrelation_perfectHeight(theta_grid, lnR_grid, centre, diameter, pitch_array)

		#Step 3: Find the peaks and then save them with the relative heights of the peaks
	
		peaks, _ = find_peaks(ccSum);
		
		if(len(peaks) == 0):
			print('no peak found in cloudID number ', cloud_ids[i])
		
		if(len(peaks)>1):
			print('more than one peak (exactly {}) found in cloudID number '.format(len(peaks)), cloud_ids[i]);


		for j in range(len(peaks)):

			#saving the cloud_id with the peak found
			useful_cloud_ids.append(cloud_ids[i]);

			peak_array.append(pitch_array[peaks[j]]);

			height_array.append(100*(ccSum[peaks[j]] - np.amin(ccSum))/perfect_height);

			max_array.append(ccSum[peaks[j]]);

		# plotCCSum(ccSum, pitch_array, cloud_ids[i]);
		i = i+1;

	#we sum all the correlations found for each (point, pitch angle);
	return np.sum(ccSum_array, 0), peak_array, height_array, max_array, useful_cloud_ids;



#Takes: Name of a specific cloud property 
#Returns: An array propert[i] == property of the ith cloud
def getCenters(hf, propertyName = 'center'):
	centers = [];
	centers_xy = [];
	cloud_ids = [];

	print('Number of clouds are = ', len(hf.keys()));
	for i in range(len(hf.keys())):
		n = hf.get('cloud_{}'.format(i));
		p = n.get('com');
		centers_xy.append(np.asarray([p[0], p[1]]));
		# print(p)
		r, theta = get_radial_coordinates(p[0], p[1]);
		if(r>=6*one_kpc):
			centers.append(np.asarray([theta, np.log(r/one_kpc)]));
			cloud_ids.append(n.get('cloud_ID')[()]);
			
	return centers, cloud_ids, centers_xy;

#Takes: Name of a specific cloud property 
#Returns: An array propert[i] == property of the ith cloud
def getCellIDs(hf, cloud_ids):
	c_ids = [];
	print('number of useful clouds - ', len(cloud_ids))
	for i in cloud_ids:
		n = hf.get('cloud_{}'.format(i));
		c_ids.append(n.get('cell_IDs')[()]);
	
	return c_ids;

#=================================================================================

def wrapperFeatherAngle(plt_number, ds):

	proj = ds.proj('density', 'z')
	width = (26, 'kpc') # we want a 26 kpc view
	res = [1000, 1000] # create an image with 1000x1000 pixels
	frb = proj.to_frb(width, res)

	ds_frb = frb.export_dataset(fields=["density"]);
	ad = ds_frb.all_data();

	x = ad['x'];
	y = ad['y'];
	density = ad['density'];

	#==========================================================
	#defining the r region we are interested in
	#==========================================================
	r_min = 5*one_kpc;
	r_max = 12*one_kpc



	(r_grid, theta) = get_radial_coordinates(x,y);

	indices = np.where((r_grid>r_min) & (r_grid<r_max))[0]


	#===========================================================
	#Placing the cut uisng the rmin and rmax and making the (lnR, theta) grids
	#=============================================================

	r_grid = r_grid[indices];
	theta = theta[indices];
	density = density[indices];
	x = x[indices];
	y = y[indices];
	lnR_grid = np.log(r_grid/one_kpc);
	theta_grid = theta;


	#============================================================
	#Getting the centers from the cloud file to work on 
	#============================================================
	print('reading file - ')
	hf = h5.File('{}clouds_gas_{}'.format(dirr,plt_number), 'r')
	centers, cloud_ids, centers_xy = getCenters(hf, 'com');

	centers_lnR_theta = np.asarray(centers);
	print('Length of centers and cloud_IDs :', len(centers), len(cloud_ids));

	diameter = 0.10 #in units of lnR
	
	pitch_array = np.linspace(0.1, 89.9);


	ccSum, peaks, heights, max_values, useful_cloud_ids = convolve_centers(centers_lnR_theta, diameter, theta_grid, lnR_grid, density, cloud_ids);
	
	print('Number of clouds with a definite peak are - ', len(useful_cloud_ids), len(peaks));

	ccSum_normal = [];
	for i in range(len(ccSum)):
		ccSum_normal.append(ccSum[i] - np.amin(ccSum[i]));



	return peaks, heights, max_values, ccSum_normal, diameter, useful_cloud_ids



import argparse
parser = argparse.ArgumentParser();


parser.add_argument('--nmin', type = int, nargs = 1, required = True , help = 'the plot file numbers to read, separated by space');
parser.add_argument('--nmax', type = int, nargs = 1, required = True , help = 'the plot file numbers to read, separated by space');
parser.add_argument('--dir', type = str, nargs = 1, required = False, help = 'The directory of the plot files')
parser.add_argument('--saveName', type = str, nargs= 1, required = True, help = 'The name of the files to save the data in')



args = parser.parse_args();


dirr = args.dir[0]

#===========================================
#Splitting the work between different processors
#=========================================


comm = MPI.COMM_WORLD
size = comm.Get_size()
my_rank = comm.Get_rank();

	#comm.Barrier synchronises all the threads/processes. 
	#All the threads have to reach this point of execution, only then the program will move forward
comm.Barrier()

n = args.nmax[0] - args.nmin[0]

count = np.zeros(size)

	#dividing the work into chunks
my_start = my_rank*(n//size) + args.nmin[0];	 #// is integer division. I divide the total tasks by number of threads. my_rank will give unique my_starts
my_end = (my_rank+1)*(n//size) -1 + args.nmin[0];



	#giving the rest over to the last chunk 
if(my_rank == size - 1):
	my_end = n-1 + args.nmin[0];

print('My start and end is - ', my_start, my_end);

if(my_rank == 0):
	print('Each processor has about {} files to work on'.format(my_end-my_start+1));
		

i = my_start

while(i<=my_end):
	if (i<10):
		ds= yt.load('{}Disc_hdf5_plt_cnt_000{}'.format(dirr, i));
	elif(i<100): 
		ds= yt.load('{}Disc_hdf5_plt_cnt_00{}'.format(dirr, i));
	else:
		ds= yt.load('{}Disc_hdf5_plt_cnt_0{}'.format(dirr, i))

	peaks, heights, max_values, ccSum, diameter, cloud_ids = wrapperFeatherAngle(i, ds);

	time = ds.current_time/one_myr;
	#========================================================================
	#We get one file at each time. 
	#This format is easier because this is how I read the cloud files
	#========================================================================

	print('dumping things the file - {}_peaks_{}.p'.format(args.saveName[0], i))



	hf = h5py.File( "{}_peaks_{}".format(args.saveName[0], i), "w" );
	g = hf.create_group('angle');
	g.create_dataset('cloud_ids', data = cloud_ids);
	g.create_dataset('peaks', data = peaks);
	g.create_dataset('heights', data = heights);
	g.create_dataset('max_values', data = max_values);
	g.create_dataset('diameter', data = diameter);
	g.create_dataset('time', data = time);

	# print('appending cloud IDs', cloud_ids)
	print('done - ',i);
	i = i+1
