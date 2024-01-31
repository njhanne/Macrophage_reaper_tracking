import numpy as np
from skimage.graph import pixel_graph
from skimage import segmentation
import pandas as pd
import matplotlib.pyplot as plt

from tifffile import imread
# from pytrackmate import trackmate_peak_import # if this isn't working check the version vs github version...
import sys
sys.path.append("./DirFileHelpers")
from find_all_files import find_all_filepaths
from trackmatexml import TrackmateXML # I didn't make this!
# https://github.com/rharkes/pyTrackMateXML
# requires lxml and (from pip) version-parser

import tkinter as tk
from tkinter import filedialog
from collections import defaultdict

### Neighbor Finding Helpers
def find_neighbors(label_mask, connectivity = 2):
  # https://stackoverflow.com/questions/72452267/finding-identity-of-touching-labels-objects-masks-in-images-using-python
  all_pairs = defaultdict(dict)
  if len(label_mask.shape) == 2:
    # I think I can expand the radius of the mask (i.e. what is considered 'touching') but making a 2D bool structure
    # for connectivity. i.e. skimage.morphology.disk(4) would be a disk of radius 4
    # no, this is not how this function works or what it does. It only finds edges that touch, not overlap
    # I should use the overlapping code from the GM130 work.
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
      pairs[pair[0]-1].append(pair[1]-1)
  return pairs


def find_distant_neighbors(label_mask, target_cells, distance):
  # may need a + or - here as the labelimage and xml are 1 off
  original_target_mask = np.isin(label_mask, target_cells) * label_mask
  # 0.65 um xy pixel
  distance_pixel = distance / 0.65

  # this works in 3D but I don't think we want it to since it is time not z, should loop over each t slice
  expanded_mask = np.zeros_like(original_target_mask)
  for frame in range(label_mask.shape[0]):
    expanded_mask[frame,:,:] = segmentation.expand_labels(original_target_mask[frame,:,:], distance_pixel)

  overlapping_stain = np.logical_and(label_mask, expanded_mask)
  overlap_cells = np.unique(label_mask[overlapping_stain])
  matches = []
  for cell in overlap_cells:
    matched_cell = np.unique(expanded_mask[np.where(label_mask == cell)])
# #   if len(matched_stain) > 1: # this prevents an error  where the entire stain overlaps with nuclei
# #     matched_stain = matched_stain[1:]
# #   if len(matched_stain) > 1: # this chooses the better match
# #     dist = distance.cdist([nuc_centroids[nuc-1]]*scaling, stain_centroids[matched_stain - 1]*scaling) # calculate distances
# #     matched_stain = matched_stain[dist.argsort()][0] # pick the match with the shortest distance
#   matched_stain = matched_stain[0] # its a list of lists, get rid of the layering
#
#   # nuc_mask = np.where(nuc_img[int(nuc_centroids[nuc][2]), :, :] == (nuc), 1, 0)
#   # stain_mask = np.where(stain_img[int(stain_centroids[matched_stain][2]), :, :] == (matched_stain), 1, 0)
#   # test_helper(nuc_mask, stain_mask, nuc_centroids, stain_centroids, [nuc, matched_stain])
#
#   matches.append([nuc, matched_stain])
# matches = np.array(matches) # matches are in image index, not centroid index
# # Find all the duplicates and delete them
# dups, dup_id, dup_count = np.unique(matches[:,1], return_inverse=True, return_counts=True)
# count_mask = dup_count > 1
# duplicates = dups[count_mask]
# to_del = np.argwhere(np.isin(matches[:,1], duplicates))
# first_matches = np.delete(matches, to_del, axis=0)



# End neighbor helpers

# xml Helpers
def get_all_tracks(tmxml, id_range=[], duplicate_split=False, break_split=False):
  all_tracks = defaultdict(list)
  if len(id_range) != 0:
    for i in id_range:
      all_tracks[i].extend(tmxml.analysetrackid(i, duplicate_split, break_split))
      # all_tracks[i].append()
  else:
    for i in range(len(tmxml.tracknames)):
      all_tracks[i].extend(tmxml.analysetrackid(i, duplicate_split, break_split))
  return all_tracks


def find_cells_greater_than_prop(tmxml, cell_propertyname, limit):
  cells  = []
  property_col_i = tmxml.spotheader.index(cell_propertyname)
  cell_property = tmxml.spots[:,property_col_i]
  cell_list = (i for i, x in enumerate(cell_property) if x >= limit)
  for i in cell_list:
    # gets the label id instead of the row number (i)
    cells.append(tmxml.spots[i,0])
  return cells


# End xml helpers

# Analysis Helpers
def get_frame_stats(tmxml, macrophages, apoptotic_cells, touching_cells, n_frames):
  afs = defaultdict(dict)
  frame_c = tmxml.spotheader.index('FRAME')
  for frame in range(n_frames):
    # get the index of all cells in this frame
    frame_cells_i = np.where(tmxml.spots[:,frame_c] == frame)[0]
    #  get the label of all cells in the frame
    afs[frame]['cells_label'] = tmxml.spots[frame_cells_i,0]
    # which of these are also known macrophages
    afs[frame]['macrophages'] = set(afs[frame]['cells_label']) & set(macrophages) # finds matches
    # which of these are known apoptotic
    afs[frame]['apoptotic'] = set(afs[frame]['cells_label']) & set(apoptotic_cells) # finds matches
    # from the list of cells that are touching, which are macrophages, i.e. which macrophages are touching cells
    afs[frame]['mac_touch'] = set(touching_cells[frame]) & set(afs[frame]['macrophages'])
    # from the touching list, which are dying
    afs[frame]['dead_touch'] = set(touching_cells[frame]) & set(afs[frame]['apoptotic'])
    # of the macrophages touching cells, which are touching apoptic cells?
    cells_macrophages_touch = [touching_cells[frame][x] for x in afs[frame]['mac_touch']]
    cells_macrophages_touch = [cell for cells in cells_macrophages_touch for cell in cells] #flatten list
    afs[frame]['mac_dead_touch'] = set(cells_macrophages_touch) & set(afs[frame]['apoptotic'])
  return(afs)


# End analysis helpers

### MAIN ###
# import label mask
image_path = filedialog.askopenfilename(title='Select segmented labelmask')
label_mask = imread(image_path)

# import xml
xml_track_path = filedialog.askopenfilename(title='Select trackmate xml')
tmxml = TrackmateXML()
tmxml.loadfile(xml_track_path)

# find touching cells
touching_cells = find_neighbors(label_mask, 2)

# find red cells
# contrast in channel 1 seems to work pretty well. Should be greater than 0.05
macrophages = find_cells_greater_than_prop(tmxml, 'MEAN_INTENSITY_CH1', 50)
#
# # find green cells
# # this one is much trickier due to the bright middle for channel 2...
# # try contrast greater SNR over 0.4
apoptotic_cells = find_cells_greater_than_prop(tmxml, 'MAX_INTENSITY_CH2', 255)


trackid_range = []
# tracks = tmxml.analysetrackid(trackid, duplicate_split=False, break_split=True)
all_tracks = get_all_tracks(tmxml, trackid_range, duplicate_split=False, break_split=True)

propertyname = 'MEAN_INTENSITY_CH1'
# property_LUT =
# frames_LUT =

# intensities = []
# frames = []
# for tracks in all_tracks:
#   intensities.append(tmxml.getproperty(tracks['spotids'], propertyname))
  # frames[tracks] = [tmxml.getproperty(track['spotids'], 'FRAME') for track in all_tracks[tracks]]


all_frame_stats = get_frame_stats(tmxml, macrophages, apoptotic_cells, touching_cells, len(label_mask))
print('test')
for frame in range(len(all_frame_stats)):
  chondro_reaper_rate = (len(all_frame_stats[frame]['dead_touch']) - len(all_frame_stats[frame]['mac_dead_touch'])) / \
                            (len(all_frame_stats[frame]['cells_label']) - len(all_frame_stats[frame]['apoptotic']) - len(all_frame_stats[frame]['macrophages']))
  macro_reaper_rate = (len(all_frame_stats[frame]['mac_dead_touch'])) / (len(all_frame_stats[frame]['macrophages']))
  print('\n frame ', frame)
  print('chondro death touch: ', chondro_reaper_rate * 100)
  print('macro death touch: ', macro_reaper_rate * 100)
  if chondro_reaper_rate == 0:
    print('macro effectiveness: inf')
  else:
    print('macro effectiveness: ', macro_reaper_rate / chondro_reaper_rate)


for track in all_tracks:
  for sub_track in all_tracks[track]:
    plt.plot(tmxml.getproperty(sub_track.spotids, 'FRAME'), tmxml.getproperty(sub_track.spotids, propertyname)) #, label='Cell ' + str(tracks[i]['cell']) + ' ; Child of '+str(tracks[i]['parent']))
plt.xlabel('frame')
plt.ylabel(propertyname)
plt.legend()
# plt.title(trackid)







print('hello')