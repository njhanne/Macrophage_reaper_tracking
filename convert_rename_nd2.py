from pathlib import Path
import numpy as np
import nd2 # allows opening of Nikon nd2 images
# from nis2pyr.convertor import convert_nd2_to_pyramidal_ome_tiff # for converting Nikon image files
from tifffile import imread, imwrite
import pandas as pd

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
adaptive_thresh = False
# rearrange = None
rearrange = [1,2,0] # this is a little confusing, you put what channel from nd2 to be in 0,1,2 position

## set directories for loading and saving
data_dir = Path("D:/UCSF/macrophage_video_analysis/")
nd2_dir = (data_dir / 'Exp3').resolve()
output_dir = (data_dir / 'processed' / '8bit_tiffs').resolve()

## Get sample info
sample_info = pd.read_csv((data_dir / 'image_log.csv').resolve())

## find the nd2 files
nd2_dirs, nd2_paths = find_all_filepaths(nd2_dir, '.nd2')

lut = np.arange(2 ** 16, dtype='uint16')
i = 0
for nd2_path in nd2_paths:
  print(str(i) + ' of ' + str(len(nd2_paths)))
  image = nd2.imread(nd2_path)
  if adaptive_thresh:
    for c in range(image.shape[1]):
      image[:,c,:,:] = adaptive_8bit(image[:,c,:,:], lut)
    image = image.astype(np.uint8)
  else:
    image = lut_display(image, lut, 0, 0xffff)
    # I think these two lines are unneeded
    # image = (image / 256).astype('uint8')
    # image = cv2.convertScaleAbs(image, alpha=1)
  if rearrange is not None:
    # axes are t,c,y,x
    image = image[np.ix_(np.arange(image.shape[0]),rearrange,np.arange(image.shape[2]),np.arange(image.shape[3]))]
  this_sample_info = sample_info.loc[(sample_info['original_filename'] == Path(nd2_path).stem + '.nd2') &
                                     (sample_info['old_subfolder'] == Path(nd2_path).parent.name)]
  output_name = this_sample_info['new_filename'].values[0] + '.ome.tiff'
  print('saving ' + output_name)
  # OME-TIFF should be TZCYX (frustrating) I think these nd2 are already like that
  output_path = Path(output_dir) / output_name
  imwrite(output_path, image, imagej=True, metadata = {'axes': 'TCYX'})
  i += 1
  # convert_nd2_to_pyramidal_ome_tiff(image, output_path, max_levels = 1)