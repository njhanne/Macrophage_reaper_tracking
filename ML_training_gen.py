import random
from math import ceil
from pathlib import Path
from tifffile import imread, imwrite

import tkinter as tk
from tkinter import filedialog

import sys
sys.path.append("./DirFileHelpers")
from find_all_files import find_all_filepaths


### MAIN ###
training_images_count = 10
channel = 3
xrange = (int((3884 / 2) - 512), int((3884 / 2) + 512))
yrange = (int((3892 / 2) - 512), int((3892 / 2) + 512))

# First find all needed directories and load them all up
root = tk.Tk()
root.withdraw()

image_directory = filedialog.askdirectory(title='Select directory of your images to be segmented')
train_directory = filedialog.askdirectory(title='Select directory to save training set')

image_dirs, image_paths = find_all_filepaths(Path(image_directory), '.tif')
number_per_stack = ceil(training_images_count / len(image_paths))
for image_path in image_paths:
  stack = imread(image_path)
  slices = random.sample(range(len(stack)), number_per_stack)
  print('hey')
  for slice in slices:
    output_name = Path(Path(image_path).stem).stem + '_sl' + str(slice) + '.tif'
    output_path = Path(train_directory) / output_name
    if xrange is not None and yrange is not None:
      imwrite(output_path, stack[slice,channel-1,xrange[0]:xrange[1],yrange[0]:yrange[1]])
    else:
      imwrite(output_path, stack[slice,channel-1,:,:])
