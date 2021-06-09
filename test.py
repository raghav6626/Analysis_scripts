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

parser.add_argument('--plt', type = int, nargs = '+', required = True , help = 'the plot file numbers to read');

args = parser.parse_args();

nos = args.plt;

for i in args.plt:
	print(i);
