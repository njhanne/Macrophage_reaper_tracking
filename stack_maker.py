import cv2
from tifffile import imwrite
import numpy as np

import tkinter as tk
from tkinter import filedialog

import re
from pathlib import Path

def find_all_filepaths(directory, extension, file_string=''):
  # this iterates through the directory and subdirectories and gets all the filepaths that end in your extension
  # '.extension' is how they should be passed
  # https://stackoverflow.com/a/59803793
  subfolders, filepaths = [], []

  for f in directory.glob('*'):
    if f.is_dir():
      subfolders.append(f)
    if f.is_file():
      if f.suffix in extension:
        if re.search(file_string, f.stem):
          filepaths.append(f)

  for dir in list(subfolders):
    sf, f = find_all_filepaths(dir, extension, file_string)
    subfolders.extend(sf)
    filepaths.extend(f)
  return subfolders, filepaths


### MAIN ###
# First ask for image directory then iteratively find all the images
root = tk.Tk()
root.withdraw()
input_directory = filedialog.askdirectory(title='Select input directory')
output_directory = filedialog.askdirectory(title='Select output directory')
image_dirs, image_paths = find_all_filepaths(Path(input_directory), '.tif')
print('hiya')


# loop through directories to get all the image trios
for dir in image_dirs:
  # get all the images for this directory
  dir_images = [image_path for image_path in image_paths if Path(image_path).parent == dir]
  stack_size = len(dir_images) # counts number of rgb images that will makeup stack. Should always be 3 for this project

  # create the stack
  if len(dir_images) != 0:
    image = cv2.imread(str(dir_images[0]))  # need this first one to get the xy size (they should all be same size)
    stack = np.empty([image.shape[0], image.shape[1], stack_size])  # create the stack (empty array)
    for idx, img_path in enumerate(dir_images):
      img = cv2.imread(str(img_path))  # read rgb image
      stack[:, :, idx] = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # convert to 8bit and add to stack

    new_filename = str(dir_images[0].stem) + '_stack'
    savename = (Path(output_directory) / (new_filename + '.tif')).resolve()
    # OME-TIFF should be TZCYX (frustrating) so we need to transpose
    stack = np.moveaxis(stack, -1, 0)  # moves the last (-1) dimension to the first dimension (0) so instead of yxz it will be cyx
    imwrite(savename, stack.astype('uint8'))