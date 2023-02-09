import numpy as np
from skimage.graph import pixel_graph
import pandas as pd
import matplotlib.pyplot as plt

from tifffile import imread
# from pytrackmate import trackmate_peak_import # if this isn't working check the version vs github version...
import sys
sys.path.append("./DirFileHelpers")
from trackmatexml import TrackmateXML # I didn't make this!
# requires lxml and (from pip) version-parser

import tkinter as tk
from tkinter import filedialog
from collections import defaultdict

### Neighbor Finding Helpers
def find_neighbors(label_mask, connectivity = 2):
  # https://stackoverflow.com/questions/72452267/finding-identity-of-touching-labels-objects-masks-in-images-using-python
  all_pairs = defaultdict(dict)
  if len(label_mask.shape) == 2:
    all_pairs[0] = find_touching_cells(label_mask, connectivity = 2)
  else:
    for slice in range(len(label_mask)):
      all_pairs[slice] = find_touching_cells(label_mask[slice, :, :], connectivity = 2)
  return all_pairs


def find_touching_cells(mask, connectivity = 2):
  # connectivity count diagonals in 2D
  g, nodes = pixel_graph(mask, mask=mask.astype(bool), connectivity=connectivity)

  g.eliminate_zeros()
  g = g.tocoo()
  center_coords = nodes[g.row]
  neighbor_coords = nodes[g.col]
  center_values = mask.ravel()[center_coords]
  neighbor_values = mask.ravel()[neighbor_coords]

  touching_masks = np.column_stack((center_values, neighbor_values))
  touching_masks = set(map(tuple, touching_masks))  # sets don't allow duplicates

  # this is ugly but it works
  pairs = defaultdict(list)
  for pair in touching_masks:
    if pair[0] != pair[1]:
      pairs[pair[0]].append(pair[1])
  return pairs


# End Neighbor Helpers


def get_all_tracks(tmxml, id_range=[], duplicate_split=False, break_split=False):
  all_tracks = defaultdict(list)
  if len(id_range) != 0:
    for i in id_range:
      all_tracks[i].append(tmxml.analysetrackid(i, duplicate_split, break_split))
      all_tracks[i].append()
  else:
    for i in range(len(tmxml.tracknames)):
      all_tracks[i].append(tmxml.analysetrackid(i, duplicate_split, break_split))
  return all_tracks

# image_path = filedialog.askopenfilename(title='Select directory of your images to be segmented')
#
# label_mask = imread(image_path)
# touching_cells = find_neighbors(label_mask, 2)

xml_track_path = filedialog.askopenfilename(title='Select directory of your images to be segmented')
tmxml = TrackmateXML()
tmxml.loadfile(xml_track_path)

print(f"the tracknames:{tmxml.tracknames}")
print(" ")
print(f"the spotheader:{tmxml.spotheader}")

# trackname = 'Track_0'
# tracks = tmxml.analysetrack(trackname, duplicate_split=False, break_split=True)

trackid_range = []
# tracks = tmxml.analysetrackid(trackid, duplicate_split=False, break_split=True)
all_tracks = get_all_tracks(tmxml, trackid_range, duplicate_split=False, break_split=True)

propertyname = 'MEAN_INTENSITY_CH1'

intensities = []
frames = []
for tracks in all_tracks:
  intensities[tracks] = [tmxml.getproperty(track['spotids'], propertyname) for track in all_tracks[tracks]]
  frames[tracks] = [tmxml.getproperty(track['spotids'], 'FRAME') for track in all_tracks[tracks]]

for track in range(len(all_tracks)):
  for i in range(len(frames[track])):
    plt.plot(frames[track][i], intensities[track][i]) #, label='Cell ' + str(tracks[i]['cell']) + ' ; Child of '+str(tracks[i]['parent']))
plt.xlabel('frame')
plt.ylabel(propertyname)
plt.legend()
# plt.title(trackid)



fig, axs = plt.subplots(len(frames),1)
fig.set_size_inches(5,15)
for i in range(len(frames)):
    axs[i].plot(frames[i], intensities[i], label='Cell ' + str(tracks[i]['cell']) + ' ; Child of '+str(tracks[i]['parent']))
    axs[i].set_xlabel('frame')
    axs[i].set_ylabel(propertyname)
    axs[i].legend()
    # axs[i].set_title(trackname)



print('hello')