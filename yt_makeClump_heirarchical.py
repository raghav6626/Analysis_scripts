import yt
from yt.data_objects.level_sets.api import *
from yt.units import dimensions
import numpy as np
import sys
import argparse
import math 

yt.toggle_interactivity();
#takes the data set - ds and saves a slice plot of the density
def _gravitational_energy(field, data):
	return data['density'] * data['gravitational_potential'];

def _virial_parameter(field, data):
	return np.abs(data['gravitational_energy']) - data['kinetic_energy'] - data['magnetic_energy'];

def add_fields():	

	yt.add_field(('gas', 'gravitational_energy'), function = _gravitational_energy, units = 'dyne/cm**2');
	yt.add_field(('gas', 'virial_parameter'), function = _virial_parameter, units= "auto", dimensions = dimensions.density);

	return;


def make_clumps(*,ds, min_density, max_density, step_size, saveName):

	data_obj = ds.all_data(); #creating a data object for the clump
	master_clump = Clump(data_obj, ('gas', 'density'));	
	
	print('Master clump formed');

	stepsize = 2;
	find_clumps(master_clump, min_density, max_density, step_size);
	leaf_clumps = master_clump.leaves;

	print('Clumps found');
	fn = master_clump.save_as_dataset(saveName, fields=["density", "cell_mass","x", "y", "z", "gravitational_potential", "velx", "vely", "velz", "temperature"])

	return;

def print_some_random_info():
	print(ds.field_list);
	print(ds.field_info['gas', 'kinetic_energy'].get_source());
	print(ds.field_info['gas', 'magnetic_energy'].get_source());
	print(ds.field_info['gas', 'gravitational_potential'].get_units());
	print(ds.field_info['gas', 'density'].get_units());
	print(ds.field_info['gas', 'gravitational_energy'].get_units());
	print(ds.field_info['gas', 'plasma_beta'].get_source());
	return;


def return_density(particle_density):
	mean_mass = 1.3;
	avogadro_no = 6.022E23;

	return (mean_mass*particle_density)/(avogadro_no);

#takes in the data = ds.all_data() 
#returns the mean density in particles per cm^(-3)
def get_meanDensity(data):
	avogadro_no = 6.022E23; 
	mean_mass = 1.3;

	density = data['density'];

	average_cm = np.average(density)
	average_particle = (average_cm*avogadro_no)/mean_mass
	return (average_cm, average_particle);


#==========================================
#Parsing the arguments
#==========================================


parser = argparse.ArgumentParser();

parser.add_argument('--min_dens', type = float, nargs = 1, help  = 'minimum density of clump (in particles/cm^3)');  #overdensity parameter

parser.add_argument('--max_dens', type = float, nargs = 1, help  = 'maximum density of clump (in particles/cm^3)');  #overdensity parameter

parser.add_argument('--step', type = float, nargs = 1, help  = 'step size from min to max');  #overdensity parameter

parser.add_argument('-i', type = str, nargs = 1, required = True, help = 'input hdf5 file directory');# input hdf5 file directory

parser.add_argument('--plt_number', type = int, nargs = 1, required = True, help = 'plt file number');# input hdf5 file directory

parser.add_argument('-n', '--name',type = str, nargs = 1, required = False, help = 'the name of the suffix to be used - hydro, beta01, beta025 etc.');#the name of the suffix to be used - hydro, beta01, beta025 etc.

parser.add_argument('-o', type = str, nargs = 1, help = 'output directory name. If not provided, the current directory is taken');

args = parser.parse_args();

if(args.o == None):
	out_dir = './'
else:
	out_dir = args.o[0];




if(args.plt_number[0]>=100):
	fileName = '{}Disc_hdf5_plt_cnt_0{}'.format(args.i[0], args.plt_number[0]);
else:
	fileName = '{}Disc_hdf5_plt_cnt_00{}'.format(args.i[0], args.plt_number[0]);

if(args.min_dens != None):
	min_density = args.min_dens[0];
	min_density_cm = return_density(args.min_dens[0]);

if(args.max_dens != None):
	max_density = args.max_dens[0];
	max_density_cm = return_density(args.max_dens[0]);

if(args.step != None):
	step = args.step[0]

clumpName = "clump_min_{}_max_{}_{}_plt_{}_heirarchical.h5".format(min_density, max_density, args.name[0], args.plt_number[0]);

#======================================
#Reading the plot file
#======================================


ds = yt.load(fileName); 
data = ds.all_data();


#=====================================
#Making and saving the clumps
#====================================

print("Making the clump file with the fileName {} \n Min density = {} \n Max Density = {} \n Step Size = {} ".format(out_dir+clumpName, min_density, max_density, step));

make_clumps(ds = ds, min_density = min_density_cm, max_density = max_density_cm, step_size = step, saveName = out_dir+clumpName);




