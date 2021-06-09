import numpy as np
import matplotlib.pyplot as plt
import sys
# sys.path.insert(1,'/home/picklerick/Downloads/AusThesis/Libraries/slug2');
import statistics
# from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:95% !important; }</style>"))
from scipy import constants as sc
G = sc.G*10**(3) #in CGS
print(G)
import yt
import matplotlib as mpl
import argparse
import math

#========================================
#The parameters of the spiral potential
#=======================================


#==================================
#Fixed Parameters
#====================================
G = sc.G*10**(3) #in CGS
mean_mass = 1.3;
avogadro_no = 6.022E23;

H = 0.18; #kpc
H_cgs = 0.18*3.086e+21; #in cms 
Rs = 7; #kpc
r_0 = 8; #kpc

v_0_kpc_sec = 200 * 3.24078e-17; #kpc/s 
v_0_kms = 200 #km/s
v_0_cgs = 2e7 #cm/s

q = 0.7; #kpc
R_c = 0.5; #kpc


#=================================
#Control Parameters
#================================
alpha = (20*np.pi)/180;
t = 0 #in years.
N = 2;
omega = 1.816*10**(-8) #rad/year coRot = 11kpc
rho_0 = 1.5 ; #particles/cm^{3}
rho_0 = (mean_mass*rho_0)/(avogadro_no); #gms/cm^(3)



constant = 4*np.pi*G*rho_0*H_cgs; #this is in CGS 




##=====================================================
##The functions
#=====================================================
def C(n):
    if(n == 1):
        return (8/3*np.pi);
    elif(n == 2):
        return 1/2;
    else:
        return 8/(15*np.pi);
    return;


#r is in kpc
#K is returned in cgs 
def K(r, n):
    r = r*3.086e+21;
    return (N*n)/(r*np.sin(alpha)); 

#r is in KPC
#D is returned in cgs
def D(r, n):
    num = 1 + K(r, n) * H_cgs + 0.3*(K(r, n) * H_cgs)**(2);
    den = 1 + 0.3*(K(r, n) * H_cgs);
    return num/den;

#r0, r are in kpc
def gamma(r, theta, t):
#     print(omega*t);
    return N*(theta - omega*t - np.log(r/r_0)/(np.tan(alpha)));

#Theta component of the potental.
#r is in kpc
#phi_theta is returned in cgs
def phi_theta(r, theta, n, t ):
    return (C(n)/(K(r, n)*D(r, n))) * np.cos(n*gamma(r, theta, t));

#r here is in kpc
#constant is in cgs
#Rs in kpc
def potential_spiral(r, theta, z, constant,t = 0):
    phi_sp_cgs = - (constant)*np.exp(-(r - r_0)/Rs) *(phi_theta(r, theta, 1, t) + phi_theta(r, theta, 2, t) + phi_theta(r, theta, 3, t));
    
#     phi_sp_km = phi_sp_cgs*1e-10
    return phi_sp_cgs;

def potential_disc(r, z):
#     phi_disc = (v_0_kms**(2)/2)*(np.log((R_c**(2) + r**(2) + (z/q)**(2))/R_c**(2)));
    phi_disc = (v_0_cgs**(2)/2)*(np.log((R_c**(2) + r**(2) + (z/q)**(2))/R_c**(2)));
    
    return phi_disc;


def potential_spiral_cartesian(x, y, t = 0,z = 0):
#     print(x, y)
    r = np.sqrt(x**2+y**2);
    theta = np.arctan2(y,x);#returns the angle back in radians
    print('calculating for t = ', t);
    phi_sp = potential_spiral(r, theta, z, constant, t);
#     print(phi_sp);
    return phi_sp;

def make_spiral(X,Y,Z,t_myrs):
    
    fig, ax = plt.subplots();

    ax.set_ylabel('y(kpc)');
    ax.set_xlabel('x(kpc)');
    fig.set_size_inches(10, 10);
# im  = ax.imshow(Z, extent = (-width_kpc/2, width_kpc/2, -width_kpc/2, width_kpc/2),norm = mpl.colors.SymLogNorm(linthresh = 0.1))
    cs = ax.contour(X, Y, Z, 3, linewidths = [0,0,2]);
#     cs = ax.contour(X, Y, Z);
    ax.clabel(cs, inline = 1, fontsize = 10);
#     CB = fig.colorbar(cs)
#     plt.savefig('/home/hslxrsrv3/bay5403/isolated_galaxies/hydro_turbulent/spiral_potential_minimas_{}.png'.format(t_myrs))   
#     plt.close();
    return;


#====================================================
#Program Main
#======================================================

parser = argparse.ArgumentParser();


parser.add_argument('-nmin', type = int, nargs = 1, help = 'starting number for plot file');
parser.add_argument('-nmax', type = int, nargs = 1, help = 'last number for the plot file'); 

parser.add_argument('-N', type = int, nargs = 1, help = 'number of spiral arms '); 

args = parser.parse_args();
 

N = args.N[0]

i = args.nmin[0];

width_kpc = 26


x = np.linspace(-width_kpc/2, width_kpc/2);
y = np.linspace(-width_kpc/2, width_kpc/2);
X,Y = np.meshgrid(x, y);

while(i<args.nmax[0]):
    if (i<10):
        ds= yt.load('Disc_hdf5_plt_cnt_000{}'.format(i));
    elif(i<100): 
        ds= yt.load('Disc_hdf5_plt_cnt_00{}'.format(i));
    else:
        ds= yt.load('Disc_hdf5_plt_cnt_0{}'.format(i));
#Creating an FRB
    res = [1000, 1000];
    width = (width_kpc, 'kpc')
    # proj = ds.slice('z',1);
    proj = ds.proj('density', 2);
    frb = proj.to_frb(width, res);
    
    t_secs = ds.current_time.v;
    t_years = t_secs/(3.154e7);

    Z = potential_spiral_cartesian(X,Y,t_years);
    

    fig, ax = plt.subplots();
    plt.style.use('classic');
    ax.tick_params(direction = 'in');
    fig.set_size_inches(12, 8);
    x_axis_label = 'x (kpc)'
    y_axis_label = 'y (kpc)'
    colobar_label = r'Projected Density $\left[\frac{\rm g}{\rm cm^{2}}\right]$'
    ax.set_xlabel(x_axis_label);
    ax.set_ylabel(y_axis_label);

    im = plt.imshow(np.array(frb['density']), cmap = 'viridis', extent = (-width_kpc/2, width_kpc/2, -width_kpc/2, width_kpc/2),norm = mpl.colors.LogNorm());

    # cs = ax.contour(X, -Y, Z, 3, linewidths = [0,2,0]);
    cs = ax.contour(X, -Y, Z, 3);
    ax.clabel(cs, inline = 1, fontsize = 10); #you can comment this out once you are sure that you are plotting the minima

    time = math.trunc((ds.current_time.v*100)/3.154e13)
    time = time/100
    
    ax.text(np.amin(x) + 1,np.amin(y) + 1, "T = {} Myrs".format(time), size = 10, color = 'white');    

    clb = fig.colorbar(im);
    clb.set_label(colobar_label);
    plt.savefig('contour_overPlot_{}'.format(i));
    plt.show()
    print('Saving the plot with the name - contour_overPlot_{}'.format(i));
    i = i+1;










