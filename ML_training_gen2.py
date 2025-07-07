import random
from math import ceil
from pathlib import Path
from tifffile import imread, imwrite

import tkinter as tk
from tkinter import filedialog

import sys
sys.path.append("./DirFileHelpers")
from find_all_files import find_all_filepaths
import re


### MAIN ###
xrange = (int((3884 / 2) - 512), int((3884 / 2) + 512))
yrange = (int((3892 / 2) - 512), int((3892 / 2) + 512))

# First find all needed directories and load them all up
root = tk.Tk()
root.withdraw()

image_directory = filedialog.askdirectory(title='Select directory of your images to be segmented')
train_directory = filedialog.askdirectory(title='Select directory to save training set')

image_dirs, image_paths = find_all_filepaths(Path(image_directory), '.tif')
old_train_img_dirs, old_train_img_paths = find_all_filepaths(Path(train_directory), '.tif')

for image_path in image_paths:
  print('processing image: ', image_path)

  matches = [old_p.parts[-1] for old_p in old_train_img_paths if old_p.parts[-1].startswith(image_path.parts[-1][:-4])]
  if len(matches) != 0:
    slices = [int( re.findall( r'(?<=_sl)(.*)(?=\.)', pn)[0]) for pn  in matches]
    stack = imread(image_path)

    for slice in slices:
      output_name = Path(Path(image_path).stem).stem + '_sl' + str(slice) + '.tif'
      output_path = Path(train_directory) / output_name
      if xrange is not None and yrange is not None:
        imwrite(output_path, stack[slice,:,xrange[0]:xrange[1],yrange[0]:yrange[1]])
      else:
        imwrite(output_path, stack[slice,:,:,:])
