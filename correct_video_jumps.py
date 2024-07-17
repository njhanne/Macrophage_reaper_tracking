# image stuff
# import cv2
from tifffile import imwrite, imread

# array / dataframe stuff
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt

import re

# path stuff
from pathlib import Path
import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


def jump_distance_calculations(jump_info):
  jump_info['delta_x'] = jump_info['x0'] - jump_info['x1']
  jump_info['delta_y'] = jump_info['y0'] - jump_info['y1']
  return jump_info


def initialize_image(x_size, y_size, t_size, c_size, jump_info):
  # need to include 0 into the cumsum (lol) or it won't calculate correctly
  y_range = max(np.append(jump_info['cumulative_y'].values, 0)) - min(np.append(jump_info['cumulative_y'].values, 0))
  x_range = max(np.append(jump_info['cumulative_x'].values, 0)) - min(np.append(jump_info['cumulative_x'].values, 0))
  # later, instead of zeros, I think it may make more sense to have it be 'average' pixel intensity background
  # if you don't specify 8bit it chooses 64 bits LMAO
  return np.zeros((t_size, c_size, y_size+y_range, x_size+x_range), dtype='uint8')


def jump_image(og_img, new_image, jump_info, x_size, y_size, t_size):
  x0 = abs(min(np.append(jump_info['cumulative_x'].values, 0)))
  y0 = abs(min(np.append(jump_info['cumulative_y'].values, 0)))
  for t in range(t_size):
    match_t_row = jump_info.loc[jump_info['t'] == t+1]
    if not match_t_row.empty:
      x0 = x0 + match_t_row['dx'].values[0]
      y0 = y0 + match_t_row['dy'].values[0]
    new_image[t, :, y0:y0+y_size, x0:x0+x_size] = og_img[t,:,:,:]
  return new_image


def bound_new_image(jump_info, og_img):
  t_size = og_img.shape[0]
  c_size = og_img.shape[1]
  y_size = og_img.shape[2]
  x_size = og_img.shape[3]
  jump_info.loc[:,'cumulative_x'] = jump_info.loc[:,'dx'].cumsum()
  jump_info.loc[:,'cumulative_y'] = jump_info.loc[:,'dy'].cumsum()
  new_image = initialize_image(x_size, y_size, t_size, c_size, jump_info)
  new_image = jump_image(og_img, new_image, jump_info, x_size, y_size, t_size)
  return new_image


### Main ###
# We have time-series fluorescent stacks taken of in vitro culture of MSC and macrophages
# I think the camera moves to other wells or something which makes the video jump
# These jumps prevent celltracker from keeping track of cell movement as they move too far for it to register
# This code should create a larger image that moves the image to prevent jumps
# Unfortunately this requires a human to manually note when the jump occurs and how far they move
# We will read in a csv of the jump info and use it to create a new image stack

# so far it looks like nearly all the jumps happen on the same timepoint and similar distances
# I may just use the same jump info for all vids

data_dir = (Path.cwd() / 'data').resolve()
data_dir = Path("D:/UCSF/macrophage_video_analysis/")
process_dir = (data_dir / 'processed' / 'BG_corrected' / 'test').resolve()
output_dir = ( data_dir / 'processed' / 'stabilized_tiffs' / 'test').resolve()
image_dirs, image_paths = find_all_filepaths(Path(process_dir), '.tiff')

image_info_csv = Path(data_dir / 'image_log.csv')
image_info = pd.read_csv(image_info_csv)

image_jump_info_csv = Path(data_dir / 'image_jump.csv')
jump_info = pd.read_csv(image_jump_info_csv)

# jump_info = jump_distance_calculations(jump_info) this is old way, new way may not need this anymore
# print('hey')


# loop through images
for image_path in image_paths:
  img_info = image_info.loc[image_info['new_filename'].str.fullmatch(re.search('(.*?)\_BG', image_path.stem).group(1))]
  # if img_info['jump_correction'].values[0] not in ('d'):
  print(img_info['new_filename'].values[0])
  original_image = imread(str(image_path)) #TCYX

  # match video to dataframe info
  # finds string match between beginning of string and '_crop'. The group tells it to  only get the first match, not the whole string
  # img_jump_info = image_info.loc[image_info['new_filename'].str.startswith(re.search('(.*?)\.ome', image_path.stem).group(1))]
  img_jump_info = jump_info.loc[jump_info['file name'] == img_info['jump_correction'].values[0]].copy() # copy is needed to avoid annoying setwithcopy warning
  new_img = bound_new_image(img_jump_info, original_image)

  # plt.imshow(new_img[50,:,:])
  # plt.show()
  savename = (Path(output_dir) / (img_info['new_filename'].values[0] + '_corrected.tiff')).resolve()
  imwrite(savename, new_img.astype('uint8'), imagej=True, metadata={'axes': 'TCYX'})

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