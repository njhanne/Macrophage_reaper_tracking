# image stuff
from tifffile import imwrite, imread

# array / dataframe stuff
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from trackmatexml import TrackmateXML # I didn't make this!
# https://github.com/rharkes/pyTrackMateXML
# requires lxml and (from pip) version-parser

from collections import defaultdict

# path stuff
from pathlib import Path
import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths

### XML Helpers
def get_all_tracks(tmxml, id_range=[], duplicate_split=True, break_split=True):
  all_tracks = defaultdict(list)
  if len(id_range) != 0:
    for i in id_range:
      all_tracks[i].extend(tmxml.analysetrackid(i, duplicate_split, break_split))
      # all_tracks[i].append()
  else:
    for i in range(len(tmxml.tracknames)):
      all_tracks[i].extend(tmxml.analysetrackid(i, duplicate_split, break_split))
  return all_tracks


# End XML Helpers


### Main ###
# The time-series images are named 24 and 48 hr. The number in the image name correspond to the same sample/well
# in the 24hr and 48hr images. So we want to continue the same tracks without having to combine them into one big image
# file or re-run all the tracking and detecting. What I will try to do here is just take the first and last slice
# from the 48 and 24 hr detected label mask files and do a simple track between them. There are several issues:
# - The spot_ids will not be the same as the previous spot ids, so we will need to link them together a different way
# - The images will not actually overlap because of the 'correct_video_jumps' script. We will need to 'uncrop' them so they
#   more closely match one another
# - Matching can maybe be done using the centroids?

# load the images and trackmate info
data_dir = (Path.cwd() / 'data').resolve()
process_dir = (data_dir / 'processed' / 'label_images').resolve()
processed_dir = (data_dir / 'processed' / 'stabilized_tiffs').resolve()
results_dir = (data_dir / 'results' / 'tracks_csv').resolve()


image_dirs, image_paths = find_all_filepaths(Path(process_dir), '.tiff')
link_xml_dirs, link_xml_paths = find_all_filepaths(process_dir, '.xml')
full_xml_dirs, full_xml_paths = find_all_filepaths(processed_dir, '.xml')
results_csv_dirs, results_csv_paths = find_all_filepaths(results_dir, '.csv')

# get all the sample info csvs
image_info_csv = Path(data_dir / 'image_log.csv')
image_info = pd.read_csv(image_info_csv)

image_jump_info_csv = Path(data_dir / 'image_jump.csv')
jump_info = pd.read_csv(image_jump_info_csv)

sample_info = pd.read_csv((data_dir / 'image_log.csv').resolve())

for sample_name in sample_info['new_filename'].unique():
  this_sample_info = sample_info.loc[sample_info['new_filename'] == sample_name]
  if this_sample_info['time'].values[0] == 24: # we don't want to run it on the pairs twice, only look at the 24hr ones
    paired_sample = sample_info.loc[sample_info['new_filename_timeless'] == this_sample_info['new_filename_timeless'].values[0]]
    if len(paired_sample) > 1: # only look at ones that actually have pairs
      paired_sample = paired_sample[paired_sample['time'] == 48]

      xml_path = [xp for xp in link_xml_paths if xp.parts[-1].startswith(this_sample_info['new_filename_timeless'].values[0]+'_link')]
      lmi_path = [lp for lp in image_paths if lp.parts[-1].startswith(this_sample_info['new_filename_timeless'].values[0]+'_link')]

      if len(lmi_path) == 0: # create the link mask file
        print('creating link mask file, you will need to run trackmate script in imagej to create the links!')
        lbl_24 = [lp for lp in image_paths if lp.parts[-1].startswith('LblImg_' + this_sample_info['new_filename'].values[0]+'_')]
        lbl_48 = [lp for lp in image_paths if lp.parts[-1].startswith('LblImg_' + paired_sample['new_filename'].values[0]+'_')]

        combined_mask = np.zeros((2, 2048, 2048), dtype='float32')
        temp_img = imread(lbl_24[0], key=177) # luckily they are all the same size
        combined_mask[0, :, :] = temp_img[0:2048, 0:2048]
        temp_img = imread(lbl_48[0], key=0)
        combined_mask[1, :, :] = temp_img[temp_img.shape[0]-2048:temp_img.shape[0], 0:2048]

        savename = (Path(process_dir) / (this_sample_info['new_filename_timeless'].values[0]+'_link.tiff')).resolve()
        imwrite(savename, combined_mask.astype('float32'), imagej=True, metadata={'axes': 'TYX'}, )

      if len(xml_path) != 0 and len(lmi_path) != 0: # link two existing csv with linking xml
        print('Analyzing ' + sample_info['new_filename_timeless'].values[0])
        tmxml_path = [tm for tm in full_xml_paths if tm.parts[-1].startswith(this_sample_info['new_filename'].values[0]+'_')]
        tmxml_24 = TrackmateXML()
        tmxml_24.loadfile(tmxml_path[0])  # import TrackMate XML

        tmxml_path = [tm for tm in full_xml_paths if tm.parts[-1].startswith(paired_sample['new_filename'].values[0]+'_')]
        tmxml_48 = TrackmateXML()
        tmxml_48.loadfile(tmxml_path[0])




