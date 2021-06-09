import numpy as np
import matplotlib.pyplot as plt
import yt
from scipy import constants as sc
import argparse
import math
import sys

def get_center_of_mass(x, m):
    
    cx = np.multiply(x, m);
    cx = np.sum(cx);
    M = np.sum(m) 

    return cx/M;


def get_distance_in_plane(leaf_clump):

    cx = get_center_of_mass(leaf_clump['grid', 'x'], leaf_clump['grid', 'cell_mass']);
    cy = get_center_of_mass(leaf_clump['grid', 'y'], leaf_clump['grid', 'cell_mass']);


    return np.sqrt(cx**(2) + cy**(2));

def plot_scatter(*, quantity_x, quantity_y, xlabel, ylabel, name, legend):

    fig, ax = plt.subplots();
    plt.style.use('default');
    ax.tick_params(direction = 'in');
    plt.xlabel(xlabel);
    plt.ylabel(ylabel)
    plt.scatter(quantity_x, quantity_y, color = '#25D1E8' , alpha = 0.8 , label = legend);
    plt.legend();    
    plt.savefig('{}_comparison'.format(name));
    return;

#for the analytical part.
def plot_scatter_larson(*, quantity_x, quantity_y, larson_y,xlabel, ylabel, name, legend):

    fig, ax = plt.subplots();
    plt.style.use('default');
    ax.tick_params(direction = 'in');
    plt.xlabel(xlabel);
    plt.ylabel(ylabel)
    plt.scatter(quantity_x, quantity_y, color = '#25D1E8' , alpha = 0.8 , label = legend);
    plt.scatter(quantity_x, larson_y, color = '#FF404B', label = r'$\sigma_{v} = 0.9 L^{0.56}$');
    plt.legend();    
    plt.savefig('{}_larson'.format(name));
    return;

#plots a normalised historgram of the two quantities
def plot_compare_hist(*, hydro, beta01, xlabel, name, label1, label2):
    fig, ax = plt.subplots();
    plt.style.use('default');
    ax.tick_params(direction = 'in');
    plt.xlabel(xlabel);
    plt.ylabel('PDF')
    plt.hist(hydro, edgecolor = 'black', color = '#59FF83' , alpha = 0.5, label = label1, density = True);
    plt.hist(beta01, edgecolor = 'black', color = '#FF404B', alpha = 0.5, label = label2, density = True);
    plt.legend();
    plt.savefig('{}_comparison_{}_{}'.format(name, label1, label2));
    return;


def plot_compare_lines(*, freq_hydro, quantity_hydro, freq_beta01, quantity_beta01, xlabel, ylabel,name = 'mass'):

    fig, ax = plt.subplots();
    plt.style.use('default');
    ax.tick_params(direction = 'in');
    plt.xlabel(xlabel);
    plt.ylabel(ylabel)
    plt.plot(quantity_hydro, freq_hydro, '--*', color = '#59FF83' , alpha = 0.8, label = 'Hydro (180 Myrs)' );
    plt.plot(quantity_beta01, freq_beta01,  '--*', color = '#FF404B', alpha = 0.8, label = r'$\beta = 1$ (280 Myrs)');
    plt.legend();    
    plt.savefig('{}_comparison'.format(name));
    return;

#takes n the set of leaf clumps
#then just a simple histrogram of the clump?
def get_density_dispersion(leaf_clumps):

    density_disp = np.zeros(len(leaf_clumps));
    
    for i in range(len(density_disp)):
        density_disp[i] = np.std(leaf_clumps[i]['grid', 'density']); #the standard deviation of the density of each clump.

    avogadro_no = 6.022E23; 
    mean_mass = 1.3;

    density_disp_particles = (density_disp*avogadro_no)/mean_mass;
    
    return np.log10(density_disp_particles);

def trim_leaf_clumps(leaf_clumps):
    selected_leaf_clumps = [];

    for i in range(len(leaf_clumps)):
        if(get_distance_in_plane(leaf_clumps[i])>inner_disc_radius*1000*one_parsec and len(leaf_clumps[i]['grid', 'cell_mass'])>cell_threshold):
            selected_leaf_clumps.append(leaf_clumps[i]);

    print('Total Number of leaf clumps - {}'.format(len(leaf_clumps)));
    print('Selected ones -', len(selected_leaf_clumps));    
    
    return selected_leaf_clumps;

def get_rearranged_velocity_dispersion(dest, n):
    f = open(dest);
    
    v_d_arranged = np.zeros(n)
    
    x = f.readline();
    while(x):
        x = x.split();
        ids =  int(x[-1]);
        v_d_arranged[ids] = float(x[-3]);
        x = f.readline();

    return v_d_arranged;

def get_radius(leaf_clump):
    size = [];

    for i in range(len(leaf_clump)):
        size.append(len(leaf_clump[i]['grid', 'cell_mass']));
    
    size = np.asarray(size);
    radius = 19.53125*(size)**(1/3); #in pc 
    
    return radius;


def get_radius(leaf_clump):
    size = [];

    for i in range(len(leaf_clump)):
        size.append(len(leaf_clump[i]['grid', 'cell_mass']));
    
    size = np.asarray(size);
    radius = 19.53125*(size)**(1/3); #in pc 
    
    return radius;

#returns radius of a given single clump in pcs.
def get_radius_single(single_clump):
    size = len(single_clump['grid', 'cell_mass']);
    return (19.53125*(size)**(1/3));



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
    
    
#v_d is arranged according to the clump number.
#we put a cut_off in here: v_d<30 km/sec
def get_virial_parameter(leaf_clumps, v_d):
    alpha_vir = [];
    density = [];
    
    avogadro_no = 6.022E23; 
    mean_mass = 1.3;
    alpha_vir_collapse = []; 
    
    for i in range(len(leaf_clumps)): 
        if(v_d[i]<30):  #this cuts-off the unphysical clumps that we see in the hydro case.
            KE = get_kinetic_energy_of_cloud(leaf_clumps[i], v_d[i]); #now this i corresponds to the same index for the velocity as well as for the leaf.
            gravitational_energy = get_gravitational_energy(leaf_clumps[i]);
            alpha_vir.append(2*KE/gravitational_energy);
            
            d_cm = np.average(leaf_clumps[i]['grid', 'density']);
            d_p = (d_cm*avogadro_no)/mean_mass;
            density.append(d_p);
            
            if(2*KE/gravitational_energy<=1):
                alpha_vir_collapse.append(2*KE/gravitational_energy);

    alpha_vir_collapse = np.asarray(alpha_vir_collapse);
    print('the fraction of collapsing clouds are : ', len(alpha_vir_collapse)/len(alpha_vir)); 
    alpha_vir = np.asarray(alpha_vir);
    density = np.asarray(density);
    print(len(alpha_vir));
    print(len(density));
    
    return (alpha_vir, density);    

def get_cumulative_freq(mass_of_clumps):
    delta = 0.50; #dex. for forming different bins
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
#   print(n);
#   print(mass_array);
    mass_array = np.asarray(mass_array);
    n = np.asarray(n);
    
    n = n/len(mass_of_clumps);
    
    return (n, mass_array);

#returns (n, mass_array): where n = frequency in the marr_array[i] bin 
def get_cumulative_mass_freq(leaf_clumps):
    mass = [];

    for i in range(len(leaf_clumps)):
        mass.append(np.sum(leaf_clumps[i]['grid', 'cell_mass']));

    mass = np.asarray(mass);

    (n, mass_array) = get_cumulative_freq(mass);

    mass_array = mass_array/m_0 
    
    return (n, np.log10(mass_array));
    
#r in pc
def get_larson_radius_relaton(r):
    v_d = 0.9*(r)**(0.56)
    return(v_d);

#m in m_0
def get_larson_gravity(m):
    v_d = 0.42*m**(0.20);
    return(v_d);

def get_surface_density(leaf_clumps):
    surface_density = [];
    
    for i in range(len(leaf_clumps)):
        mass = np.sum(leaf_clumps[i]['grid', 'cell_mass']); #total mass in gms
        mass = mass/m_0; # in m_0 now
        
        r_c = get_radius_single(leaf_clumps[i]); #in pc already
        
        surface_density.append(mass/(4*np.pi*r_c**(2))); #in M_0/pc^(2)
    
    surface_density = np.asarray(surface_density);
    
    return surface_density;


def projection_with_colored_clumps(ds, leaf_clumps, name, plt_number,save = False):
    width_param = 20
    prj = yt.ProjectionPlot(ds, 'z', 'density', width = (width_param, 'kpc'));
    for i in range(len(leaf_clumps)):
        ax = prj.annotate_marker([[leaf_clumps[i]['grid', 'x']], [leaf_clumps[i]['grid', 'y']], [leaf_clumps[i]['grid', 'z']]], coord_system='data' ,plot_args={'color':'red','s':10}, marker='.')	

    prj.annotate_timestamp();
    if(save):
        prj.save('{}{}_plt_{}projection_with_colored_clumps_selected'.format(name, args.ovd[0], plt_number));
    else:
        prj.show();
    return;

def return_mass(leaf_clumps):
	mass = [];

	for i in range(len(leaf_clumps)):
		mass.append(np.sum(leaf_clumps[i]['grid', 'cell_mass']));

	mass = np.asarray(mass);	

	mass = mass/m_0 
	
	return mass;
#=============================================
#Defining constants 
#===============================================

m_0 = 1.989e+33 #in grams
one_parsec = 3.086e+18 #in cms
inner_disc_radius = 2 #in kpc
cell_threshold = 50 #The minimum number of cells that should be there in the clump.


#=========================================================
#Reading the command line arguments 
#========================================================

parser = argparse.ArgumentParser();
parser.add_argument('-d1', type = str, nargs = 1, help = 'The directory of the first clump');

parser.add_argument('-d2', type = str, nargs = 1, help = 'The directory of the second clump');

parser.add_argument('--plt_number_1', type = int, nargs = 1, help = 'The plotFile number of the first clump')

parser.add_argument('--plt_number_2', type = int, nargs = 1, help = 'The plotFile number of the second clump')

parser.add_argument('--cell_threshold', type = int, nargs = 1, help = 'The min. number of cells in a clump...')

parser.add_argument('-o', type = str, nargs = 1, help = 'output directory name. If not provided, the current directory is taken');

parser.add_argument('--name1',type = str, nargs = 1, required = False, help = 'the name of the suffix to be used - hydro, beta01, beta025 etc.');

parser.add_argument('--name2',type = str, nargs = 1, required = False, help = 'the name of the suffix to be used - hydro, beta01, beta025 etc.');

parser.add_argument('-ovd', type = float, nargs = 1, help = 'density_parameter');

parser.add_argument('--clumpName1', type = str, nargs = 1, help = 'name of the clump being read...')

parser.add_argument('--clumpName2', type = str, nargs = 1, help = 'name of the clump being read...')


args = parser.parse_args();

if(args.clumpName1==None):
    clumpFilename_one = None
else:
    clumpFilename_one = args.d1[0]+args.clumpName1[0]

if(args.clumpName2==None):
    clumpFilename_two = None
else:
    clumpFilename_two = args.d2[0]+args.clumpName2[0]

 
    



if(args.cell_threshold != None):
	cell_threshold = args.cell_threshold[0];


#=============================================
#Reading the clumps and the plot files
#==============================================

if(args.plt_number_1[0]<100):
	plotFileName_one = '{}Disc_hdf5_plt_cnt_00{}'.format(args.d1[0], args.plt_number_1[0])
	if(args.plt_number_2[0]<100):
		plotFileName_two = '{}Disc_hdf5_plt_cnt_00{}'.format(args.d2[0], args.plt_number_2[0])
	else:
		plotFileName_two = '{}Disc_hdf5_plt_cnt_0{}'.format(args.d2[0], args.plt_number_2[0])

else:
	plotFileName_one = '{}Disc_hdf5_plt_cnt_0{}'.format(args.d1[0], args.plt_number_1[0])
	if(args.plt_number_2[0]>100):
		plotFileName_two = '{}Disc_hdf5_plt_cnt_0{}'.format(args.d2[0], args.plt_number_2[0])
	else:
		plotFileName_two = '{}Disc_hdf5_plt_cnt_00{}'.format(args.d2[0], args.plt_number_2[0]) 
		

if(clumpFilename_one== None):
    clumpFilename_one = args.d1[0] +'clump_threshold_density_{}_{}_{}.h5'.format(args.ovd[0], args.name1[0], args.plt_number_1[0]);

if(clumpFilename_two == None):
    clumpFilename_two = args.d2[0] +'clump_threshold_density_{}_{}_{}.h5'.format(args.ovd[0], args.name2[0], args.plt_number_2[0]);


plot1 = yt.load(plotFileName_one);
plot2 = yt.load(plotFileName_two)

clump_hydro = yt.load(clumpFilename_one);
clump_beta01 = yt.load(clumpFilename_two);


leaf_clumps_hydro = clump_hydro.leaves;
leaf_clumps_beta01 = clump_beta01.leaves;
#========================================
#Plotting the clumps on top of a projectionPlot 
#=======================================

#projection_with_colored_clumps(plot1, leaf_clumps_hydro, args.name1[0], args.plt_number_1[0], True);

#projection_with_colored_clumps(plot2, leaf_clumps_beta01, args.name2[0], args.plt_number_2[0] , True);

#=========================================================
#selecting the leaf clumps so that they have 
#===========================================================
#1. Enough number of cells>cell_threshold
#2. Do not lie in the inner disc
#3. For the hydro case, we also trim it further by selecting clumps with v_d < 30 km/sec

selected_leaf_hydro = trim_leaf_clumps(leaf_clumps_hydro);
selected_leaf_beta01 = trim_leaf_clumps(leaf_clumps_beta01);
#v_d_hydro = get_rearranged_velocity_dispersion('./hydro/velocity_dispersion_parallel.txt', len(selected_leaf_hydro));

#selected_leaf_hydro_v_d = [];

#for i in range(len(v_d_hydro)):
#    if(v_d_hydro[i]<30):
#        selected_leaf_hydro_v_d.append(selected_leaf_hydro[i]);


###############################
#Density Dispersion historgrams
'''
density_disp_hydro = get_density_dispersion(selected_leaf_hydro);
density_dsip_beta01 = get_density_dispersion(selected_leaf_beta01);

density_disperson_label = r'$\log(\sigma (\rho))$ $particles/cm^{3}$' ;

plot_compare_hist(hydro= density_disp_hydro, beta01=density_dsip_beta01, xlabel=density_disperson_label, name ='density', label1 = args.name1[0]+'{}'.format(args.plt_number_1), label2 = args.name2[0]+'{}'.format(args.plt_number_2));
'''
################################
'''
#Velocity dispersion Histograms. 
#CAUTION - The cell threshold used in the vd calculation and here should be the same, else all the clump's velocity dispersions would not be read. 

v_d_hydro = get_rearranged_velocity_dispersion('./data/velocity_dispersion_parallel_{}.txt'.format(args.clumpName1[0]), len(leaf_clumps_hydro));

print(v_d_hydro);

v_d_beta =  get_rearranged_velocity_dispersion('./data/velocity_dispersion_parallel_{}.txt'.format(args.clumpName2[0]), len(leaf_clumps_beta01));


v_d_hydro_plot = [];
v_d_beta_plot = [];

radius_hydro = [];
radius_beta01 = [];

for i in range(len(v_d_hydro)):
    if(v_d_hydro[i] > 0): #checking against the seleted clumps. Only the selected clumps will have a non-zero vd
        v_d_hydro_plot.append(v_d_hydro[i]);  #adding the vd of the selected clump for plotting
        radius_hydro.append(get_radius_single(leaf_clumps_hydro[i])); #adding the radius of the same clump (using ID)
        
for i in range(len(v_d_beta)):
    if(v_d_beta[i]>0):
        v_d_beta_plot.append(v_d_beta[i]);
        radius_beta01.append(get_radius_single(leaf_clumps_beta01[i]));
        
v_d_hydro_plot = np.asarray(v_d_hydro_plot);
v_d_beta_plot = np.asarray(v_d_beta_plot);
radius_hydro = np.asarray(radius_hydro);
radius_beta01 = np.asarray(radius_beta01);

larson_v_d_hydro = get_larson_radius_relaton(radius_hydro);
larson_v_d_beta = get_larson_radius_relaton(radius_beta01);

#v_d histo gram
velocity_dispersion_label = r'$\sigma (v)$ (km/s)';
plot_compare_hist(hydro=v_d_hydro_plot, beta01=v_d_beta_plot, xlabel= velocity_dispersion_label, name ='velocity_dispersion', label1 = '{}_{}'.format(args.name1[0], args.plt_number_1[0]), label2 = '{}_{}'.format(args.name2[0], args.plt_number_2[0]));


#scatter plot
plot_scatter_larson(quantity_x= radius_hydro, quantity_y = v_d_hydro_plot , larson_y= larson_v_d_hydro, xlabel = 'R (pc)', ylabel = r'$\sigma (v)$ (km/s)', legend = '{}_{}'.format(args.
name1[0], args.plt_number_1[0]), name = '{}_{}'.format(args.name1[0],args.plt_number_1[0]));

plot_scatter_larson(quantity_x= radius_beta01, quantity_y = v_d_beta_plot , larson_y= larson_v_d_beta, xlabel = 'R (pc)', ylabel = r'$\sigma (v)$ (km/s)', 
                    legend = '{}_{}'.format(args.name2[0], args.plt_number_2[0]), name ="{}_{}".format(args.name2[0],args.plt_number_2[0]));

'''
#########################
#Radius histograms
hydro_radius = get_radius(selected_leaf_hydro);
beta01_radius = get_radius(selected_leaf_beta01);

radius_label = 'Radius (pc)'

plot_compare_hist(hydro=hydro_radius, beta01=beta01_radius, xlabel= radius_label,name= 'radius', label1 = args.name1[0]+'{}'.format(args.plt_number_1), label2 = args.name2[0]+'{}'.format(args.plt_number_2));


#=============================
#The virial parameter:
#=============================
'''
(alpha_vir_hydro, density_hydro) = get_virial_parameter(selected_leaf_hydro, v_d_hydro);
(alpha_vir_beta01, density_beta01) = get_virial_parameter(selected_leaf_beta01, v_d_beta);
print(density_hydro);
plot_scatter(quantity_x= density_hydro, quantity_y=alpha_vir_hydro , ylabel = r'$\alpha_{vir}$', xlabel= 'Density (particles/$cm^{3}$)', name = 'hydro_virial_density' , legend = 'Hydro (180 Myrs)');
plot_scatter(quantity_x= density_beta01, quantity_y=alpha_vir_beta01 , ylabel = r'$\alpha_{vir}$', xlabel= 'Density (particles/$cm^{3}$)', name = 'beta_virial_density' , legend = r'$\beta = 1$ (280 Myrs)');

# alpha_vir_hydro_cut = [];
# alpha_vir_beta01_cut = [];

# for i in range(len(alpha_vir_beta01)):
#     if(alpha_vir_beta01[i]<10):
#         alpha_vir_beta01_cut.append(alpha_vir_beta01[i]);
        
# for i in range(len(alpha_vir_hydro)):
#     if(alpha_vir_hydro[i]<10):
#         alpha_vir_hydro_cut.append(alpha_vir_hydro[i]);

alpha_label = r'$\log(\alpha_{vir})$'
plot_compare_hist(hydro = np.log10(alpha_vir_hydro), beta01= np.log10(alpha_vir_beta01), xlabel=alpha_label, name = 'virial');
'''
############################################################
#cumulative mass distributions: 
#########################################################
'''
(n_hydro, mass_hydro) = get_cumulative_mass_freq(selected_leaf_hydro_v_d);
(n_beta01, mass_beta01) = get_cumulative_mass_freq(selected_leaf_beta01);

xlabel = '$ \log (Mass /(M_{\odot})) $'
ylabel = '$N(M>M)$'
plot_compare_lines(quantity_hydro=mass_hydro, freq_hydro= n_hydro, quantity_beta01= mass_beta01, freq_beta01= n_beta01, xlabel = xlabel, ylabel = ylabel);
'''
###########################################
#surface density plot
###########################################
'''
surface_density_hydro = get_surface_density(selected_leaf_hydro);
surface_density_beta01= get_surface_density(selected_leaf_beta01);

xlabel = r'$\log(\Sigma_{gas}) (M_{\odot} pc^{-2})$';
plot_compare_hist(hydro= np.log10(surface_density_hydro), beta01= np.log10(surface_density_beta01), xlabel = xlabel, name = 'surface_density',label1 = args.name1[0]+'{}'.format(args.plt_number_1), label2 = args.name2[0]+'{}'.format(args.plt_number_2));
'''
#=====================================
#Mass Historgram
#=====================================

mass_one = return_mass(selected_leaf_hydro)
mass_two = return_mass(selected_leaf_beta01)

xlabel = r'$\log$ Mass (M_${\odot}$)';
plot_compare_hist(hydro = np.log10(mass_one), beta01 = np.log10(mass_two), xlabel = xlabel, name = 'mass_hist',label1 = args.name1[0]+'{}'.format(args.plt_number_1),  label2 = args.name2[0]+'{}'.format(args.plt_number_2));



