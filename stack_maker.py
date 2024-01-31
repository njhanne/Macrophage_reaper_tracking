import cv2
from tifffile import imwrite
import numpy as np

import tkinter as tk
from tkinter import filedialog

import re
from pathlib import Path

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths



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