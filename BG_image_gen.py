from pathlib import Path
import numpy as np
# from nis2pyr.convertor import convert_nd2_to_pyramidal_ome_tiff # for converting Nikon image files
from tifffile import imread, imwrite
import pandas as pd

import sys
sys.path.append("./DirFileHelpers")
from find_all_files import find_all_filepaths



### MAIN ###
## set directories for loading and saving
data_dir = Path("D:/UCSF/macrophage_video_analysis/")
input_dir = (data_dir / 'processed' / 'BG_corrected' / 'test2').resolve()
output_dir = (data_dir / 'processed' / 'BG_corrected' / 'test2').resolve()


## find the tiffs
img_dirs, img_paths = find_all_filepaths(input_dir, '.tif')

i = 0
for img_path in img_paths:
  print(str(i+1) + ' of ' + str(len(img_paths)))
  image = imread(img_path)
  if i == 0:
    avg = np.zeros((len(img_paths),image.shape[2], image.shape[3]))
  avg[i,:,:] = np.mean(image[:,2,:,:], axis=0)
  i += 1

# average all the averages
avg = np.mean(avg, axis=0).astype(np.uint8)
print('saving bg image')
# OME-TIFF should be TZCYX (frustrating) I think these nd2 are already like that
output_path = Path(output_dir) / 'bg_img.tif'
imwrite(output_path, avg)

# apply the average to make a bg correction
# for img_path in img_paths:
#   print(str(i+1) + ' of ' + str(len(img_paths)))
#   image = imread(img_path)
#   if i == 0:
#     avg = np.zeros((len(img_paths),image.shape[2], image.shape[3]))
#   avg[i,:,:] = np.mean(image[:,2,:,:], axis=0)
#   i += 1
