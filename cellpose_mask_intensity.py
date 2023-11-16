import cv2
import numpy as np
import pandas as pd

from pathlib import Path
import re

import tkinter as tk
from tkinter import filedialog


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


def get_avg_intensity(image, channel, nuc_mask, nuclei):
  masked_mask = np.where(nuc_mask == nuclei, 1, 0) # this creates a boolean mask of 1 where nuclei id = nuclei and 0 everywhere else
  masked_int = np.mean(image[:,:,channel], where=np.array(masked_mask, dtype=bool)) # gets avg only where the mask is, only in selected channel
  return masked_int


### MAIN ###
# First find all needed directories and load them all up
root = tk.Tk()
root.withdraw()

mask_directory = filedialog.askdirectory(title='Select cellpose output masks directory')
image_directory = filedialog.askdirectory(title='Select cellpose input stack directory')

image_dirs, image_paths = find_all_filepaths(Path(image_directory), '.tif')
mask_dirs, cp_mask_paths = find_all_filepaths(Path(mask_directory), '.png')
mask_paths = [mp for mp in cp_mask_paths if Path(mp).name.endswith('_masks.png')]

# results_df = initialize_results_dataframe(sample_info, len(stain_mask_paths.count))
overall_results = []

for mask_path in mask_paths:
  print(str('analyzing ' + Path(mask_path).stem))
  print((Path(mask_path).stem)[:-9])
  img_path = [ip for ip in image_paths if Path(ip).name.startswith(Path(mask_path).stem[:-9])]
  img = cv2.imread(str(img_path[0]), -1)
  mask = cv2.imread(str(mask_path), -1)

  # fill out a little dataframe of the name of the image and the number of nuclei found
  result = []
  result.append(str(Path(mask_path).stem))
  result.append(np.amax(mask)) #gets count of nuclei found by cellpose


  # a bit more complicated, make a csv of the name of each nuclei and the average intensity of TUNEL
  ind_results = []
  ### very important you need to change this number based on which channel is TUNEL!
  ### also very important - opencv opens image as BGR instead of RGB, so in this case the blue channel is 0
  channel = 0
  for i in range(1, np.amax(mask)):
    if i % 100 == 0:
      print('nuclei ' + str(i) + ' of ' + str(np.amax(mask)-1))
    ind_result = []
    ind_result.append(str('nuclei_' + str(i)))
    ind_result.append(get_avg_intensity(img, channel, mask, i))
    ind_results.append(ind_result)
  ind_results_df = pd.DataFrame(ind_results, columns = ['nuclei', 'avg_TUNEL'])
  ind_results_df['rel_TUNEL'] = ind_results_df['avg_TUNEL'] / ind_results_df['avg_TUNEL'].mean()
  ind_results_df = ind_results_df.sort_values(by=['avg_TUNEL'], ascending=False)
  ind_results_df.to_csv((Path(mask_directory) / str(Path(mask_path).stem + '_nuc_TUNEL_intensity.csv')).resolve())

  dead_count = len(ind_results_df[ind_results_df['rel_TUNEL'] > 5])
  result.append(dead_count)
  overall_results.append(result)

results_df = pd.DataFrame(overall_results, columns = ['filename', 'nuclei_count', 'TUNEL_pos'])
results_df.to_csv((Path(mask_directory) / 'nuclei_count.csv').resolve())