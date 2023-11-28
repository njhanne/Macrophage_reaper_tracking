# image stuff
import cv2
from tifffile import imwrite, imread

# array / dataframe stuff
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import re

# path stuff
from pathlib import Path
import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


def jump_distance_calculations(jump_info):
  jump_info['delta_x'] = jump_info['x1'] - jump_info['x0']
  jump_info['delta_y'] = jump_info['y1'] - jump_info['y0']
  return jump_info


def initialize_image(x_size, y_size, t_size, jump_info):
  # need to include 0 into the cumsum (lol) or it won't calculate correctly
  y_range = max(np.append(jump_info['cumulative_y'].values, 0)) - min(np.append(jump_info['cumulative_y'].values, 0))
  x_range = max(np.append(jump_info['cumulative_x'].values, 0)) - min(np.append(jump_info['cumulative_x'].values, 0))
  # later, instead of zeros, I think it may make more sense to have it be 'average' pixel intensity background
  return np.zeros((t_size, y_size+y_range, x_size+x_range))


def jump_image(og_img, new_image, jump_info, x_size, y_size, t_size):
  x0 = abs(min(np.append(jump_info['cumulative_x'].values, 0)))
  y0 = abs(min(np.append(jump_info['cumulative_y'].values, 0)))
  for t in range(t_size):
    match_t_row = jump_info.loc[jump_info['t'] == t+1]
    if not match_t_row.empty:
      x0 = x0 + match_t_row['delta_x'].values[0]
      y0 = y0 + match_t_row['delta_y'].values[0]
    new_image[t, y0:y0+y_size, x0:x0+x_size] = og_img[t,:,:]
  return new_image


def bound_new_image(jump_info, og_img):
  x_size = og_img.shape[2]
  y_size = og_img.shape[1]
  t_size = og_img.shape[0]
  jump_info['cumulative_x'] = jump_info['delta_x'].cumsum()
  jump_info['cumulative_y'] = jump_info['delta_y'].cumsum()
  new_image = initialize_image(x_size, y_size, t_size, jump_info)
  new_image = jump_image(og_img, new_image, jump_info, x_size, y_size, t_size)
  return new_image


### Main ###
# We have time-series fluorescent stacks taken of in vitro culture of MSC and macrophages
# I think sometimes someone comes in and switches media on the cells which makes the video jump
# These jumps prevent celltracker from keeping track of cell movement as they move too far for it to register
# This code should create a larger image that moves the image to prevent jumps
# Unfortunately this requires a human to manually note when the jump occurs and how far they move
# We will read in a csv of the jump info and use it to create a new image stack

data_dir = (Path.cwd() / 'data').resolve()
process_dir = (data_dir / '8bit-tiffs' / 'test').resolve()
image_dirs, image_paths = find_all_filepaths(Path(process_dir), '.tif')

image_jump_info_csv = Path(data_dir / 'image_jump.csv')
jump_info = pd.read_csv(image_jump_info_csv)
jump_info = jump_distance_calculations(jump_info)
print('hey')


# loop through images
for image_path in image_paths:
  original_image = imread(str(image_path))

  # match video to dataframe info
  # finds string match between beginning of string and '_crop'. The group tells it to  only get the first match, not the whole string
  img_jump_info = jump_info.loc[jump_info['file name'].str.startswith(re.search('(.*?)_crop', image_path.stem).group(1))]
  new_img = bound_new_image(img_jump_info, original_image)

  plt.imshow(new_img[50,:,:])
  plt.show()


  # # create the stack
  # if len(dir_images) != 0:
  #   image = cv2.imread(str(dir_images[0]))  # need this first one to get the xy size (they should all be same size)
  #   stack = np.empty([image.shape[0], image.shape[1], stack_size])  # create the stack (empty array)
  #   for idx, img_path in enumerate(dir_images):
  #     img = cv2.imread(str(img_path))  # read rgb image
  #     stack[:, :, idx] = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # convert to 8bit and add to stack
  #
  #   new_filename = str(dir_images[0].stem) + '_stack'
  #   savename = (Path(output_directory) / (new_filename + '.tif')).resolve()
  #   # OME-TIFF should be TZCYX (frustrating) so we need to transpose
  #   stack = np.moveaxis(stack, -1, 0)  # moves the last (-1) dimension to the first dimension (0) so instead of yxz it will be cyx
  #   imwrite(savename, stack.astype('uint8'))