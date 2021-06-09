import numpy as np
import matplotlib.pyplot as plt
import yt
from yt.units import kpc, pc
import argparse
parser = argparse.ArgumentParser();

parser.add_argument('--name', type = str, nargs = 1, help = 'name of the plot file');


args = parser.parse_args();
plotFile = args.name[0];


ds = yt.load(plotFile);
ad = ds.all_data();
disc = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 5*kpc, 200*pc);


profile = disc.profile([('index','radius')], [('gas','density')]);
profile1 = disc.profile([('index','radius')], [('gas','pressure_gradient_magnitude')]);


density = profile['gas','density'].v;
pressure_grad = profile1['gas','pressure_gradient_magnitude'].v;

radius = profile.x/3.086e+21;
f = pressure_grad/density;



fig, ax = plt.subplots();
ax.set_xlabel('r (kpc)');
ax.set_ylabel('fr$_{pressure}$ (dyne gm$^{-1}$)');
ax.minorticks_on();
ax.plot(np.log10(f), radius);
plt.savefig('f_pressure.png');






