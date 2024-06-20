import random
from math import ceil
from pathlib import Path
import numpy as np
# import pandas as pd
# import cv2
import nd2 # allows opening of Nikon nd2 images
# from nis2pyr.convertor import convert_nd2_to_pyramidal_ome_tiff # for converting Nikon image files
from tifffile import imread, imwrite

import tkinter as tk
from tkinter import filedialog

import sys
sys.path.append("./DirFileHelpers")
from find_all_files import find_all_filepaths


# 8bit conversion helpers
# https://stackoverflow.com/questions/14464449/using-numpy-to-efficiently-convert-16-bit-image-data-to-8-bit-for-display-with
def display(image, display_min, display_max):  # copied from Bi Rico
  # Here I set copy=True in order to ensure the original image is not
  # modified. If you don't mind modifying the original image, you can
  # set copy=False or skip this step.
  image = np.array(image, copy=True)
  image.clip(display_min, display_max, out=image)
  image -= display_min
  np.floor_divide(image, (display_max - display_min + 1) / 256,
                  out=image, casting='unsafe')
  return image.astype(np.uint8)


def lut_display(image, lut, display_min, display_max):
  lut = display(lut, display_min, display_max)
  return np.take(lut, image)


def find_display_range(image):
  # https://stackoverflow.com/questions/60059626/converting-3-channel16-bit-image-to-8-bit-while-preserving-color
  lo, hi = np.percentile(image, (0.1, 99.9)).astype('int')
  return lo, hi


def adaptive_8bit(image, lut):
  min_range, max_range = find_display_range(image)
  return lut_display(image, lut, min_range, max_range)


# End 8bit conversion helpers

### MAIN ###
# specify settings
# whether to convert nd2 to tiff
convert_nd2 = True
convert_8bit = True
adaptive_thresh = False

training_images_count = 100
channel = 3

# First find all needed directories and load them all up
root = tk.Tk()
root.withdraw()

image_directory = filedialog.askdirectory(title='Select directory of your images to be segmented')
train_directory = filedialog.askdirectory(title='Select directory to save training set')


if convert_nd2:
  image_dirs, image_paths = find_all_filepaths(Path(image_directory), '.nd2')
  convert_directory = filedialog.askdirectory(title='Select directory of your tiffs to be saved')
  for image_path in image_paths:
    image = nd2.imread(image_path)
    if convert_8bit:
      lut = np.arange(2 ** 16, dtype='uint16')
      if adaptive_thresh:
        for c in range(image.shape[1]):
          image[:,c,:,:] = adaptive_8bit(image[:,c,:,:], lut)
        image = image.astype(np.uint8)
      else:
        image = lut_display(image, lut, 0, 0xffff)
        # image = (image / 256).astype('uint8')
        # image = cv2.convertScaleAbs(image, alpha=1)
    nd2_name = Path(image_path).stem
    output_name = nd2_name + '_24hr.ome.tiff'
    # OME-TIFF should be TZCYX (frustrating) I think these nd2 are already like that
    output_path = Path(convert_directory) / output_name
    imwrite(output_path, image, imagej=True, metadata = {'axes': 'TCYX'})
    # convert_nd2_to_pyramidal_ome_tiff(image, output_path, max_levels = 1)
#
# image_dirs, image_paths = find_all_filepaths(image_directory, '.tiff')
# number_per_stack = ceil(training_images_count / len(image_paths))
# for image_path in image_paths:
#   stack = imread(image_path)
#   slices = random.sample(range(len(stack)), number_per_stack)
#   for slice in slices:
#     output_name = Path(Path(image_path).stem).stem + '_sl' + str(slice) + '.tif'
#     output_path = Path(train_directory) / output_name
#     imwrite(output_path, stack[slice,channel-1,:,:])
