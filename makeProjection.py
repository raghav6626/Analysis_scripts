import numpy as np
import matplotlib.pyplot as plt
import yt
import matplotlib as mpl
import sys
import argparse

one_parsec = 3.08567e18 #in cms
#mpl.use('tkagg');
#=======================================
#Reading the parsed arguments.
#=======================================

parser = argparse.ArgumentParser();

 #density_parameter in particles/cm^(3)
parser.add_argument('--width', type = float, nargs = 1, help = 'width of the plane in KPC');

parser.add_argument('-i', type = str, nargs = 1, help = 'the input hdf5 file directory');#the input hdf5 file

parser.add_argument('--plt_number', type = int, nargs = 1, required = True , help = 'the plot file number');

parser.add_argument('-o', type = str, nargs = 1, help = 'output directory name. If not provided, the current directory is taken');#the input hdf5 file

parser.add_argument('--cmap', type = int, nargs = 1, help = 'The colorMap to use in the colorbar')

parser.add_argument('--axis', type = str, nargs = 1, help = 'x, y or z')

parser.add_argument('--quantity', type = int, nargs = 1, help = 'Density, magnetig_field etc...Default - density')

parser.add_argument('--clumpName', type = str, nargs = 1, help = 'name of the clump being read...')

parser.add_argument('--name', type = str, nargs = 1, help = 'name used to build the output fileName e.g. hydro_vrms_10 etc.')

parser.add_argument('--particle', type = bool, nargs= 1, help = 'whether to annotate particles or not')

args = parser.parse_args();


if(args.width!= None):
    width = args.width[0];
else:
    width = 20

if(args.i!= None):
    directory = args.i[0];
else:
    directory = './'

if(args.plt_number[0]>=100):
    fileName = directory + 'Disc_hdf5_plt_cnt_0{}'.format(args.plt_number[0]);
elif(args.plt_number[0]>=10):
    fileName = directory + 'Disc_hdf5_plt_cnt_00{}'.format(args.plt_number[0]);
else:
    fileName = directory + 'Disc_hdf5_plt_cnt_000{}'.format(args.plt_number[0]);
if(args.o!= None):
    saveName = args.o[0] + 'Projection_{}_{}.pdf'.format(args.name[0],args.plt_number[0]);
elif (args.name != None):
    saveName = './Projection_{}_{}.png'.format(args.name[0], args.plt_number[0]);
else:
    saveName = './Projection_{}_{}.png'.format(args.plt_number[0]);



if(args.cmap!= None):
    cmap = args.cmap[0]
else:
    cmap = 'viridis'
    
if(args.axis!= None):
    axis_to_use = args.axis[0]
else:
    axis_to_use = 2;
    
if(args.quantity!= None):
    quantity = args.quantity[0]
else:
    quantity = 'density';

if(args.clumpName==None):
    clumpName = None
else:
    clumpName = args.clumpName[0]


#Creating an FRB
print('Reading file ', fileName)
ds = yt.load(fileName);
res = [1000, 1000];
width_kpc = (width, 'kpc');
proj = ds.proj(quantity, axis_to_use);
frb = proj.to_frb(width_kpc, res);


fig, ax = plt.subplots();
plt.style.use('classic');
ax.minorticks_on();
# fig.set_size_inches(12, 8);
x_axis_label = 'x (kpc)'
y_axis_label = 'y (kpc)'
colobar_label = r'Projected Density $\left[\frac{\rm g}{\rm cm^{2}}\right]$'

# ,norm = mpl.colors.LogNorm()

ax.set_xlabel(x_axis_label);
ax.set_ylabel(y_axis_label);
im = plt.imshow(np.array(frb[quantity]), cmap = cmap, extent = (-width/2, width/2, -width/2, width/2),norm = mpl.colors.LogNorm());


#Annotating the leaf_clumps onto the projected graph.
if(clumpName != None):
    clump = yt.load(directory+clumpName);
    leaf_clumps = clump.leaves; 
    print('Reading clump ', clumpName);
    
    for i in range(len(leaf_clumps)):
        plt.scatter(leaf_clumps[i]['grid', 'x']/(1000*one_parsec), -leaf_clumps[i]['grid','y']/(1000*one_parsec), c = '#F54336', s = 2, edgecolors = 'face');

clb = fig.colorbar(im);
clb.set_label(colobar_label);
fig.show()


fig.savefig(saveName);
print('saving the plot  - ', saveName);




