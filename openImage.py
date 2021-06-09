from PIL import Image
import argparse

parser = argparse.ArgumentParser();


parser.add_argument('--name', type = str, nargs = 1, help = 'name of the image to open and display');

args = parser.parse_args();

img = Image.open(args.name[0]);

print('Calling the display, takes a few seconds.')

img.show();
