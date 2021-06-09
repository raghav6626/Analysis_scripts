#File names will be {--name}_<plottyype>_{--ovd}
import yt 
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import math
from scipy import constants as sc


def get_center_of_mass(x, m):

	cx = 0;
	mx = 0
	for i in range(len(x)):
		cx = cx + x[i]*m[i] ;
		mx = mx + m[i];

	return cx/mx;



def get_cumulative_freq(mass_of_clumps):
	delta = 0.5; #dex. for forming different bins
	mass_min = np.amin(mass_of_clumps);
	mass_array = []
	n = [];

	while(mass_min<np.amax(mass_of_clumps)):
		count = 0;
		for i in range(len(mass_of_clumps)):
			if(mass_of_clumps[i]>mass_min):
				count = count + 1;

		mass_array.append(mass_min);
		n.append(count);
		mass_min = mass_min*10**(delta);
	print(n);
	print(mass_array);
	mass_array = np.asarray(mass_array);

	return (n, mass_array);

#returns a cumulative_mass histogram of each clump that was found
def mass_histogram(leaf_clumps):
	mass = [];

	for i in range(len(leaf_clumps)):
		mass.append(np.sum(leaf_clumps[i]['grid', 'cell_mass']));

	mass = np.asarray(mass);

	(n, mass_array) = get_cumulative_freq(mass);

	mass_array = mass_array/m_0 
	
	plt.figure();
	plt.xlabel('$ \log (Mass /(M_{\odot})) $')
	plt.ylabel('$\log(N(M>M))$');
	plt.style.use('classic');	
	# plt.hist(np.log10(mass_array), edgecolor = 'black', cumulative = True);
	plt.plot(np.log10(mass_array), np.log10(n), '--o');
	plt.savefig('{}{}_{}_cumulative_mass_histogram_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));

	mass = mass/m_0; #converting into solar mass units.

	plt.figure();
	# plt.xscale('log')
	plt.xlabel('$ \log (Mass (M_{\	odot})) $')
	plt.style.use('classic');
	plt.hist(np.log10(mass), edgecolor = 'black');
	plt.savefig('{}{}_{}_mass_histogram_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));


	return;
#average density in particles/cm^(3)
#should I also do median density?
def density_histogram(leaf_clumps):
	density = [];

	for i in range(len(leaf_clumps)):
		density.append(np.average(leaf_clumps[i]['grid','density']));
	
	density = np.asarray(density); #this is g/cm^(3)

	avogadro_no = 6.022E23; 
	mean_mass = 1.3;

	density_particle = (density*avogadro_no)/mean_mass;
	
	plt.figure();
	plt.xlabel(' $\log(Density(particles/cm^(-3))$');
	plt.ylabel('Frequency');
	plt.style.use('classic');
	plt.hist(np.log10(density_particle), edgecolor = 'black');
	plt.savefig('{}{}_{}_density_histogram_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));

	return;

def temperature_histogram(leaf_clumps):

	temperature = [];

	for i in range(len(leaf_clumps)):
		temperature.append(np.average(leaf_clumps[i]['grid','temperature']));
	
	temperature = np.asarray(temperature); #this is g/cm^(3)

	plt.figure();
	plt.xlabel('Temperature (K)');
	plt.ylabel('Frequency');
	plt.style.use('classic');
	plt.hist(temperature, edgecolor = 'black');
	plt.savefig('{}{}_{}_temperature_histogram_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));



	return;
#gives a histogram of the number of cells in leaf clumps. 
#and also the radius of the clumps  = (V)**(1/3)
def get_size(leaf_clumps):

	size = [];

	for i in range(len(leaf_clumps)):
		size.append(len(leaf_clumps[i]['grid', 'cell_mass']));
	
	size = np.asarray(size);
	plt.figure();
	# plt.xscale('log')
	plt.xlabel('log(Cells)')
	plt.style.use('classic');
	plt.hist(np.log10(size), edgecolor = 'black');
	plt.savefig('{}{}_{}_size_histogram_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));

	radius = 19.53125*(size)**(1/3); #in pc 

	plt.figure();
	# plt.xscale('log')
	plt.xlabel(' log(Radius (pc))')
	plt.style.use('classic');
	plt.hist(np.log10(radius), edgecolor = 'black');
	plt.savefig('{}{}_{}_radius_histogram_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));	

	# print(size);
	return;

#returns the virial parameter of a particular leaf_clump
def get_virial(leaf_clump):

	v_com = get_COM_velocity(leaf_clump);
	vx = leaf_clump['velx'] - v_com[0];
	vy = leaf_clump['vely'] - v_com[1];
	vz = leaf_clump['velz'] - v_com[2];

	kinectic_energy = 0.5*leaf_clump['cell_mass'] * (vx**(2) + vy**(2) + vz**(2));
	gravitational_energy = leaf_clump['cell_mass'] * leaf_clump['gravitational_potential'];

	return;

def virial_histogram(leaf_clumps):
	virial = np.zeros(len(leaf_clumps));

	for i in range(len(leaf_clumps)):
		virial[i] = get_virial(leaf_clumps[i]);


	plt.figure();
	plt.style.use('classic');
	plt.hist(mass, edgecolor = 'black');
	plt.savefig('virial_histogram.pdf');
	return;

#coloring the leaf_clumps' center of mass, based on the number of cells in it (volume)
def plot_clumps_only(leaf_clumps):

	# points = np.zeros(len(leaf_clumps));
	size = np.zeros(len(leaf_clumps));
	cx = np.zeros(len(leaf_clumps));
	cy = np.zeros(len(leaf_clumps));
	cz = np.zeros(len(leaf_clumps));

	for i in range(len(leaf_clumps)):
		size[i] = len(leaf_clumps[i]['grid', 'cell_mass']);
		cx[i] = get_center_of_mass(leaf_clumps[i]['grid', 'x']);
		cy[i] = get_center_of_mass(leaf_clumps[i]['grid', 'y']);
		cz[i] = get_center_of_mass(leaf_clumps[i]['grid', 'z']);

	plt.figure();
	plt.style.use('classic');
	plt.scatter(cx/one_parsec, cy/one_parsec, c = size);
	plt.colorbar();	
	plt.savefig('{}_clumps_with_color_{}.pdf'.format(args.name[0],density_parameter));

	return size, cx, cy, cz;

#make a scatter plot on top of the projection plot. You will have to be careful about the units na.
#I have failed to add the colorbar in this one. I still do not know how to proceed. 
def projection_with_colored_clumps(leaf_clumps):
	width_param = 20
	#(size, cx, cy, cz) = plot_clumps_only(leaf_clumps);

	prj = yt.ProjectionPlot(ds, 'z', 'density', width = (width_param, 'kpc'));
	#prj.set_cmap(field = 'density', cmap ='Greys');
	
#	ax = prj.annotate_marker([[cx], [cy], [cz]],coord_system='data',plot_args={'c':size, 's':250}, marker='.');
	for i in range(len(leaf_clumps)):
		ax = prj.annotate_marker([[leaf_clumps[i]['grid', 'x']], [leaf_clumps[i]['grid', 'y']], [leaf_clumps[i]['grid', 'z']]], coord_system='data' ,plot_args={'color':'red','s':10}, marker='.')	
	# fig = prj.plots['density'].figure;	
	# ax =  prj.plots['density'].axes ;
	# # cax = fig.add_axes([0.5, 1, 0.5, 0.6]);
	# im = ax.scatter(cx, cy, c = size, s = 200);
	# fig.colorbar(im, orientation = 'horizontal');
	prj.set_zlim('density', 1e-7, 10**(-0.8));
	prj.annotate_timestamp();
	prj.save('{}{}_{}projection_with_colored_clumps_selected_{}'.format(output_dir, args.name[0], args.year[0], density_parameter));
	return;

#I do not think I should include the z_coordinate? Just to be relevant with the observations ??
#but since it is negligble in comparison to the x and y, it won't hurt right?
#takes a simple_clump and returns the distance from the centre. 
def distance_from_centre(leaf_clump):
	cx = get_center_of_mass(leaf_clump['grid', 'x'], leaf_clump['grid', 'cell_mass']);
	cy = get_center_of_mass(leaf_clump['grid', 'y'], leaf_clump['grid', 'cell_mass']);
	cz = get_center_of_mass(leaf_clump['grid', 'z'], leaf_clump['grid', 'cell_mass']);

	return (np.sqrt(cx**(2) + cy**(2) +  cz**(2)));


#tries to find a correlation between the r and the mass of the clumps found.
def distance_mass_plot(leaf_clumps, inner_disc_radius):
	mass = [];
	r_c = [];
	size = [];

	for i in range(len(leaf_clumps)):
		if(get_distance_in_plane(leaf_clumps[i])>1000*inner_disc_radius*one_parsec):
			mass.append(np.sum(leaf_clumps[i]['grid', 'cell_mass'])/m_0);
			r_c.append(distance_from_centre(leaf_clumps[i]));
			size.append(len(leaf_clumps[i]['grid', 'cell_mass']));
	
	mass = np.asarray(mass);
	r_c = np.asarray(r_c);
	size = np.asarray(size);
	plt.figure();
	plt.xlabel('Distance from center (kpc)');
	plt.ylabel('$ \log (Mass(M_{\odot}))$');
	plt.xlim(0, 9)
	plt.scatter(r_c/(one_parsec*1000), np.log10(mass), c = size);
	plt.colorbar(label = 'Number of Cells');
	plt.savefig('{}_clump_correlations_{}.pdf'.format(args.name[0],density_parameter));
	return;

def get_distance_in_plane(leaf_clump):
	
	cx = get_center_of_mass(leaf_clump['grid', 'x'], leaf_clump['grid', 'cell_mass']);
	cy = get_center_of_mass(leaf_clump['grid', 'y'], leaf_clump['grid', 'cell_mass']);
	return np.sqrt(cx**(2) + cy**(2));

#takes in the data = ds.all_data() 
#returns the mean density in particles per cm^(-3)
def get_meanDensity(data):
	avogadro_no = 6.022E23; 
	mean_mass = 1.3;

	density = data['density'];

	average_cm = np.average(density)
	average_particle = (average_cm*avogadro_no)/mean_mass
	return (average_cm, average_particle);

def get_largest_clump_index(leaf_clumps):
	largest_index = 0
	for i in range(len(leaf_clumps)):
		if(len(leaf_clumps[i]['grid', 'density']) > len(leaf_clumps[largest_index]['grid', 'cell_mass'])):
			largest_index = i;

	return largest_index;

def zoom_in_plot(ds, leaf_clumps):
	index_of_largest = get_largest_clump_index(leaf_clumps);

	cx = get_center_of_mass(leaf_clumps[index_of_largest]['grid', 'x'], leaf_clumps[index_of_largest]['grid', 'cell_mass']);
	cy = get_center_of_mass(leaf_clumps[index_of_largest]['grid', 'y'], leaf_clumps[index_of_largest]['grid', 'cell_mass']);
	cz = get_center_of_mass(leaf_clumps[index_of_largest]['grid', 'z'], leaf_clumps[index_of_largest]['grid', 'cell_mass'])
	
	prject_plot = yt.ProjectionPlot(ds, 'z', 'density', width = (20, 'kpc'));
	prject_plot.annotate_marker((cx, cy , cz),  coord_system='data', plot_args={'color':'red', 's':50}, marker='x');	
	prject_plot.set_zlim('density', 1e-7, 10**(-0.8));
	prject_plot.save('{}{}_{}_zoom_out_{}'.format(output_dir, args.name[0], args.year[0], density_parameter));


	n = len(leaf_clumps[index_of_largest]['grid', 'cell_mass']);
	width_param = 4*(n)**(1/3)/50     				# n x 20/1000 kpc
	print('Number of cells in the larges clump : ', n);
	print('Width param in KPC', width_param);
	print('Center of mass ', cx, cy);

	prj = yt.ProjectionPlot(ds, 'z', 'density', width=(width_param,'kpc'));

	prj.set_center((cx, cy));
	prj.save('{}{}_{}_zoom_in_{}'.format(output_dir, args.name[0], args.year[0], density_parameter));
	return

#takes in the particular leaf clump to plot.
#then just a simple histrogram of the clump?
def density_dispersion(leaf_clump):

	density = leaf_clump['grid', 'density'];

	avogadro_no = 6.022E23; 
	mean_mass = 1.3;

	density_particle = (density*avogadro_no)/mean_mass;

	plt.figure();
	plt.xlabel(' $\log(Density(particles/cm^{-3})$');
	plt.ylabel('Frequency');
	plt.style.use('classic');
	plt.hist(np.log10(density_particle), edgecolor = 'black');
	plt.savefig('{}{}_{}_density_dispersion_particle_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));


	plt.figure();
	plt.xlabel(' $\log(Density(gm/cm^{-3})$');
	plt.ylabel('Frequency');
	plt.style.use('classic');
	plt.hist(np.log10(density), edgecolor = 'black');
	plt.savefig('{}{}_{}_density_dispersion_gm_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));

	return;

def get_center_of_mass_velocity(leaf_clump):
	mass = leaf_clump['grid', 'cell_mass'];
	vcm = [0,0,0];
	total_mass = 0;
	print(leaf_clump['grid', 'velx'][0]);

	for i in range(len(mass)):
		v = [leaf_clump['grid', 'velx'][i], leaf_clump['grid', 'vely'][i], leaf_clump['grid', 'velz'][i]] #ith cell vel
		vcm = np.add(vcm ,np.multiply(v,mass[i])); #vcm = sum(mi*vi)
		total_mass = total_mass + mass[i]; #m = sum(m);
#		print(vcm);

	return vcm/total_mass; #vcm =  sum(mi*vi)/sum(m);

def get_com(r, m):
	rcm = [0, 0, 0];
	m_tot = 0
	for i in range(len(m)):
		rcm = np.add(rcm,np.multiply(r[i],m[i])) ; #vector equation.
		m_tot = m_tot + m[i]; #scalar equaton

	return np.divide(rcm, m_tot);

#Problem: I get Vcom in units of cm/sec but not the individual arrays when I copy them
#Sol: Use the array.v to separate the values from the units that have been used.
def get_angular_momentum(leaf_clump, v_com):

    m = leaf_clump['grid', 'cell_mass'];
    r = np.zeros((len(m), 3));
    v = np.zeros((len(m), 3));

    for i in range(len(m)):
        r[i] = [leaf_clump['grid', 'x'][i], leaf_clump['grid', 'y'][i], leaf_clump['grid', 'z'][i]];
        v[i] = [leaf_clump['grid', 'velx'][i], leaf_clump['grid', 'vely'][i], leaf_clump['grid', 'velz'][i]];


    r_com = get_com(r, m); 

    r_rel = np.subtract(r, r_com.v);    
#     print('The average distance from the center of mass (in pc): ', r_rel/one_parsec);
    v_rel = np.subtract(v, v_com.v);

    L = [0,0,0];

    for i in range(len(m)):
        L = L +  np.multiply(m[i].v,np.cross(r_rel[i], v_rel[i])); #L = m(r x v) wrt to the center of mass of the clump.
    
    return L;

def get_moment_of_Inertia(leaf_clump):
	m = leaf_clump['grid', 'cell_mass'].v;
	r = np.zeros((len(m), 3));

	for i in range(len(m)):
		r[i] = [leaf_clump['grid', 'x'][i].v, leaf_clump['grid', 'y'][i].v, leaf_clump['grid', 'z'][i].v];	

	r_com = get_com(r, m);

	r_rel = np.subtract(r, r_com);

	I = 0;
	for i in range(len(m)):
		I = I + m[i] * (r_rel[i][0] **2 + r_rel[i][1]**2 + r_rel[i][2] ** 2);


	return (I, r_rel);

#takes in a single leaf_clump
#returns the velocity dispersion of a single clump
def velocity_dispersion(leaf_clump):

    v_cm = get_center_of_mass_velocity(leaf_clump); #a vector
    L = get_angular_momentum(leaf_clump, v_cm);	#a vector again
    (I, r_rel) = get_moment_of_Inertia(leaf_clump);	#a scalar (not realy but okay)

    omega = L/I;
    
    print(omega);
    
    v_rot = np.cross(omega, r_rel); #v = omega x r_rel 
    
    #initializing the v 
    for i in range(len(leaf_clump['grid', 'cell_mass'])): 
        v = [leaf_clump['grid', 'velx'][i].v, leaf_clump['grid', 'vely'][i].v, leaf_clump['grid', 'velz'][i].v] #ith cell velocity
    
    v = v - v_cm.v - v_rot; # cm/sec  
    v = v*10**(-5)        #km/sec #these are n vectors right now. 
    
    v_squared = np.zeros(len(v));
    
    for i in range(len(v)): 
        v_squared[i] = v[i][0]**(2) + v[i][1]**(2) + v[i][2]**(2); #now I have the v**(2)

    v_rms = np.sqrt(np.average(v_squared)); #v_rms = sqrt( mean (v^(2)))
    
    return (v_rms);

#takes in the velocity dispersion and saves the figure.
def velocity_dispersion_histogram(v_d):
    v_d_plot = []
    for i in range(len(v_d)):
        if(v_d[i]>45):
            print('v_d to cut ', v_d[i])
            i = i+1;
        else:
            v_d_plot.append(v_d[i]);

    plt.figure();
    plt.style.use('classic');
    plt.xlabel('Velocity Dispersion (km/sec)');
    plt.ylabel('Frequency');
    plt.hist(v_d_plot, edgecolor = 'black');
    plt.savefig('{}{}_{}_velocity_dispersion_{}.pdf'.format(output_dir, args.name[0], args    .year[0], density_parameter));

    return v_d;


#takes n the set of leaf clumps
#returns the density dispersion histogram of all the clumps taken together.
def density_dispersion(leaf_clumps):

    density_disp = np.zeros(len(leaf_clumps));
    
    for i in range(len(density_disp)):
        density_disp[i] = np.std(leaf_clumps[i]['grid', 'density'])/np.average(leaf_clumps[i]['grid', 'density']); #the standard deviation of the density of each clump.

    avogadro_no = 6.022E23; 
    mean_mass = 1.3;

    #density_disp_particles = (density_disp*avogadro_no)/mean_mass;

    plt.figure();
    plt.xlabel('$\log(Desnity_dipsersion/mean)$');
    plt.ylabel('Frequency');
    plt.style.use('classic');
    plt.hist(np.log10(density_disp), edgecolor = 'black');
    plt.savefig('{}{}_{}_density_dispersion_mean_log_particle_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));

    return;

def get_kinetic_energy_of_cloud(leaf_clump, v_d): #in CGS units
    m = leaf_clump['grid', 'cell_mass'];
    m_cloud = np.sum(m);
    v_d = v_d*10**(5);           #from km/sec to cm/sec
    return (m_cloud*(v_d)**(2))/2

def get_gravitational_energy(leaf_clump): #in CGS units
    size = len(leaf_clump['grid', 'cell_mass']);
    radius = 19.5*(size)**(1/3) * 3.08567758128E+18; #from pc to cms
    m_cloud = np.sum(leaf_clump['grid', 'cell_mass']); #in gms
    
    return (sc.G*10**(3)*m_cloud**(2))/(radius);

def get_rearranged_velocity_dispersion(number_of_clumps):
#     f = open('/home/picklerick/Downloads/Study/Hamburg/Work/clumps/velocity_dispersion_parallel.txt')
    f = open('velocity_dispersion_parallel.txt');
    
    v_d_arranged = np.zeros(number_of_clumps);
    
    x = f.readline();
    while(x):
        x = x.split();
        ids =  int(x[-1]);
        v_d_arranged[ids] = float(x[-3]);
        x = f.readline();

    return v_d_arranged;

#takes in the v_d and leaf_clumps to plot a hist of the virial parameter
def virial_histogram(leaf_clumps):
    v_d = get_rearranged_velocity_dispersion(len(leaf_clumps)); #this is arranged according to the clump number. 
    alpha_vir = np.zeros(len(leaf_clumps));
    alpha_vir_collapse = [];
    for i in range(len(leaf_clumps)):
        if(v_d[i]>45):
            i = i+1;
        else:
            KE = get_kinetic_energy_of_cloud(leaf_clumps[i], v_d[i]); #now this i corresponds to the same index for the velocity as well as for the leaf.
            gravitational_energy = get_gravitational_energy(leaf_clumps[i]);
            alpha_vir[i] = KE/gravitational_energy;
            if(alpha_vir[i]<=1):
                alpha_vir_collapse.append(alpha_vir[i]);

    alpha_vir_collapse = np.asarray(alpha_vir_collapse);
    print('number of collapsing clouds are : ', len(alpha_vir_collapse));
    plt.figure();
    plt.style.use('classic');
    plt.xlabel(r'$\alpha _{vir}$');
    plt.ylabel('$Frequency$');
    plt.hist(np.abs(alpha_vir), edgecolor = 'black');
    plt.savefig('{}{}_{}_virial_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));
    return;


def read_velocity_dispersion():
    v_d = np.genfromtxt('velocity_dispersion_parallel.txt');
    v_d = v_d[:,4];
    return v_d

def get_sizes_of_velocity_dispersion():
    f = open('velocity_dispersion_parallel.txt');
#     f = open('/home/picklerick/Downloads/Study/Hamburg/velocity_dispersion.txt');

    size = [];
    x = f.readline();
    while(x):
        x = x.split();
    
        y = x[1].split('=');
    
        size.append(int(y[1]));
        x = f.readline();
    
    size = np.asarray(size);
    return size;

def plot_v_d_size(v_d, size):
    
    v_d_to_plot = [];
    size_to_plot = [];
    for i in range(len(v_d)):
        if (i == 58 or i == 105 or i == 123 or i == 135):
            i = i+1;
        else:
            v_d_to_plot.append(v_d[i]);
            size_to_plot.append(size[i]);
    
    plt.figure();
    plt.ylabel(r'$\sigma (v)$');
    plt.xlabel('Radius (pc)');
    plt.style.use('classic');
    plt.scatter(19.53125*(size)**(1/3), v_d, label = '{}'.format(args.name[0]));
    plt.legend(loc = 'upper right');
    plt.savefig('{}{}_{}_velocity_disperson_with_size_{}.pdf'.format(output_dir, args.name[0], args.year[0], density_parameter));
    return;


parser = argparse.ArgumentParser();

parser.add_argument('-d', type = float, nargs = 1, help = 'overdensity_parameter'); #density_parameter in particles/cm^(3)
parser.add_argument('-ovd', type = float, nargs = 1, help = 'density_parameter');
parser.add_argument('-i', type = str, nargs = 1, required = True , help = 'the input hdf5 file directory');#the input hdf5 file
parser.add_argument('--year', type = int, nargs = 1, required = True , help = 'the year (has to be even)');
parser.add_argument('-n', '--name',type = str, nargs = 1, required = False, help = 'the name of the suffix to be used - hydro, beta01, beta025 etc.');#the name of the suffix to be used - hydro, beta01, beta025 etc.
parser.add_argument('-o', type = str, nargs = 1, help = 'output directory name. If not provided, the current directory is taken');#the input hdf5 file

parser.add_argument('--clumpName', type = str, nargs = 1, help = 'name of the clump being read...')



args = parser.parse_args();

if(args.d == None):
	density_parameter = args.ovd[0];
else:
	density_parameter = args.d[0];

if(args.o == None):
	output_dir = './';
else:
	output_dir = args.o[0];

if(args.clumpName==None):
    clumpName = "{}_{}_clumps_{}.h5".format(args.name[0], args.year[0], density_parameter);
else:
    clumpName = args.clumpName[0]
# clumpName = '../clump_threshold_density.h5';
if(args.year[0]/2>=100):
	fileName = '{}Disc_hdf5_plt_cnt_0{}'.format(args.i[0], math.trunc(args.year[0]/2));
else:
	fileName = '{}Disc_hdf5_plt_cnt_00{}'.format(args.i[0], math.trunc(args.year[0]/2));

# fileName = '{}Disc_hdf5_plt_cnt_0{}'.format(args.i[0], math.trunc(args.year[0]/2));
print(fileName);
ds = yt.load(fileName); 
data = ds.all_data();
ds_clumps = yt.load(clumpName);


m_0 = 1.989e+33 #in grams
one_parsec = 3.086e+18 #in cms
inner_disc_radius = 2 #in kpc
cell_threshold = 50 #The minimum number of cells that should be there in the clump.

leaf_clumps = ds_clumps.leaves

#selecting the leaf clumps so that they have 
#1. Enough number of cells>cell_threshold
#2. Do not lie in the inner disc

selected_leaf_clumps = [];

for i in range(len(leaf_clumps)):
	if(get_distance_in_plane(leaf_clumps[i])>inner_disc_radius*1000*one_parsec and len(leaf_clumps[i]['grid', 'cell_mass'])>cell_threshold):
		selected_leaf_clumps.append(leaf_clumps[i]);
		#print('selected IDs are : ', i);

print('Total Number of leaf clumps - {}'.format(len(leaf_clumps)));
print('Selected ones -', len(selected_leaf_clumps));

v_d = get_rearranged_velocity_dispersion(len(selected_leaf_clumps))
velocity_dispersion_histogram(v_d);
#size = get_sizes_of_velocity_dispersion();
#plot_v_d_size(v_d, size);


virial_histogram(selected_leaf_clumps);

	
#density_dispersion(selected_leaf_clumps);

# index = get_largest_clump_index(selected_leaf_clumps);

# density_dispersion(leaf_clumps[index]);
# mass_histogram(Selectedd_leaf_clumps);
# # #distance_mass_plot(leaf_clumps, inner_disc_radius);
# density_histogram(selected_leaf_clumps);
# get_size(selected_leaf_clumps);
#projection_with_colored_clumps(selected_leaf_clumps);
# # zoom_in_plot(ds, selected_leaf_clumps)
# temperature_histogram(selected_leaf_clumps);

#projection plot without the clumps
