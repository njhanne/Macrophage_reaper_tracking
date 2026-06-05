from cellpose import models, io
from pathlib import Path
import time

import numpy as np

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### main ###
# load in all the needed files and paths
# data_dir = (Path.cwd() / 'data' / 'processed').resolve()
data_dir = (Path.cwd() / 'data' / '20260304' / 'Timelapse_1').resolve()

# get all the stacks we created from rgb_to_stack.py
image_directory = (data_dir / 'BG_combined').resolve()
image_dirs, images = find_all_filepaths(image_directory, '.tif')
# find the cellpose models that Charlie helped train
model_path = (data_dir / 'train' / 'models' / 'NH_LB_LC5_HL_3chan').resolve()
# where to save the cellpose output
output_directory = (data_dir / 'cellpose_output').resolve()

### IMPORTANT! ###
# Change this to False if you did not install the cuda pytorch version when installing cellpose!
cuda = True
io.logger_setup()
model = models.CellposeModel(gpu=cuda, pretrained_model=str(model_path))

for image_path in images:
  filename = Path(image_path).name
  print('analyzing image: ', filename)
  image = io.imread(str(image_path))
  save_filename = filename[:-4]
  print('image loaded')

  big_masks = np.empty((image.shape[0], image.shape[2], image.shape[3]), dtype=np.uint16) # should be big enough
  # list of 3 arrays [z,x,y,c]
  #                  [2,z,x,y]
  #                  [z,x,y]
  # big_flows = [np.empty((image.shape[0], image.shape[2], image.shape[3], image.shape[1]), dtype=np.uint16),
  #              np.empty((2, image.shape[0], image.shape[2], image.shape[3]), dtype=np.uint16),
  #              np.empty((image.shape[0], image.shape[2], image.shape[3]), dtype=np.uint16)]
  for i in range(image.shape[0]):
    print(str(i+1) + ' of ' + str(image.shape[0]+1))
    masks, flows, styles = model.eval(image[i,:,:,:], channel_axis=0,
                                        normalize={'lowhigh': None, 'percentile': [1.0, 99.0], 'normalize': True,
                                                   'sharpen_radius': 0.0, 'smooth_radius': 0.0,
                                                   'tile_norm_blocksize': 0.0, 'invert': False})
    big_masks[i] = masks
    # big_flows[0][:half_img_k,:,:,:] = flows[0]
    # big_flows[1][:, :half_img_k, :, :] = flows[1]
    # big_flows[2][:half_img_k, :, :] = flows[2]
    i += 1
  print('saving')
  # flows is not optional in this function but it doesn't do anything with it for saving the mask image so we don't
  # need to figure out how to compile it. Can pass it a list I guess
  io.save_masks(image, big_masks, [1,2,3], save_filename, tif=True, savedir=str(output_directory), save_txt=False)
