import numpy as np
from skimage.graph import pixel_graph
from skimage import segmentation
import pandas as pd
import matplotlib.pyplot as plt
import pickle

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
     pairs[pair[0]].append(pair[1])
  return pairs


def find_distant_neighbors(label_mask, target_cells, distance):
  # may need a + or - here as the labelimage and xml are 1 off
  original_target_mask = np.zeros_like(label_mask)
  # for frame in range(label_mask.shape[0]):
    # original_target_mask[frame,:,:] = np.isin(label_mask[frame,:,:],  [corrections[frame][x] for x in target_cells[frame]]) * label_mask[frame,:,:]
  # original_target_mask = np.isin(label_mask, target_cells) * label_mask
  # 0.65 um xy pixel
  distance_pixel = distance / 0.65

  # this works in 3D but I don't think we want it to since it is time not z, should loop over each t slice
  matches = defaultdict(dict)
  for frame in range(label_mask.shape[0]):
    frame_macs = np.unique(original_target_mask[frame,:,:])
    expanded_mask = segmentation.expand_labels(original_target_mask[frame,:,:], distance_pixel)
    # overlapping_stain = np.logical_and(label_mask[frame,:,:], expanded_mask)
    # overlap_cells = np.unique(label_mask[frame,:,:][overlapping_stain])
    for mac in frame_macs:
      if mac != 0:
        matched_cells = np.unique(label_mask[frame,:,:][np.where(expanded_mask == mac)])
        matched_cells = matched_cells[(matched_cells != 0) & (matched_cells != mac)]
        if matched_cells.size != 0:
          matches[frame][mac] = matched_cells
  return matches


# End neighbor helpers

# xml Helpers
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


def create_cell_tracks(tmxml, all_frame_stats):
  tracks = get_all_tracks(tmxml)
  cell_id = 0
  all_cells = []
  for track in tracks:
    for cell in tracks[track]:
      cell_track = {}
      cell_track['cell_id'] = cell_id + cell.cell
      cell_track['spot_ids'] = list(cell.spotids)
      if cell.parent != 0:
        cell_track['parent_id'] = cell_id + cell.parent
      else:
        cell_track['parent_id'] = None
      all_cells.append(cell_track)
    cell_id += len(tracks[track])
  return pd.DataFrame(all_cells)



def find_cells_greater_than_prop(tmxml, cell_propertyname, limit):
  cells  = defaultdict(list)
  property_col = tmxml.spotheader.index(cell_propertyname)
  frame_c = tmxml.spotheader.index('FRAME')
  cell_property = tmxml.spots[:,property_col]
  cell_list = (i for i, x in enumerate(cell_property) if x >= limit)
  for i in cell_list:
    # gets the label id instead of the row number (i)
    cells[tmxml.spots[i,frame_c]].append(tmxml.spots[i,0])
  return cells


def get_label_correction(label_mask, tmxml):
  corrections = defaultdict(dict)
  rev_corrections = defaultdict(dict)

  cell_id_c = tmxml.spotheader.index('ID')
  positionx_c = tmxml.spotheader.index('POSITION_X')
  positiony_c = tmxml.spotheader.index('POSITION_Y')
  frame_c = tmxml.spotheader.index('FRAME')

  first_cell_id = min(tmxml.spots[:,cell_id_c])
  first_cell_i = np.where(tmxml.spots[:,cell_id_c] == first_cell_id)[0][0]

  label_start = label_mask[int(tmxml.spots[first_cell_i,frame_c]), int(tmxml.spots[first_cell_i,positiony_c]), int(tmxml.spots[first_cell_i,positionx_c])]
  label_delta = label_start - first_cell_id #is it always 16bit? probably how it's coded in trackmate...
  turnover_id = 65536 - label_start + first_cell_id # when it turns over we will lose a cell because it's mask is saved as 0, which is background

  for frame in range(label_mask.shape[0]):
    for cell_i in np.where(tmxml.spots[:,frame_c] == frame)[0]:
      corresponding_label = (tmxml.spots[cell_i,cell_id_c] + label_delta) % 65536
      if corresponding_label == 0:
        corresponding_label = 65537
      corrections[frame][tmxml.spots[cell_i,cell_id_c]] = corresponding_label
    rev_corrections[frame] = {v: k for k, v in corrections[frame].items()}
  return corrections, rev_corrections


# End xml helpers

# Analysis Helpers
def get_frame_stats(tmxml, macrophages, apoptotic_cells, touching_cells, neighbors, n_frames):
  afs = defaultdict(dict)
  frame_c = tmxml.spotheader.index('FRAME')
  for frame in range(n_frames):
    # get the index of all cells in this frame
    frame_cells_i = np.where(tmxml.spots[:,frame_c] == frame)[0]
    #  get the label of all cells in the frame
    afs[frame]['cells_label'] = tmxml.spots[frame_cells_i,0]
    # which of these are also known macrophages
    afs[frame]['macrophages'] = set(afs[frame]['cells_label']) & set(macrophages[frame]) # finds matches
    # which of these are known apoptotic
    afs[frame]['apoptotic'] = set(afs[frame]['cells_label']) & set(apoptotic_cells[frame]) # finds matches
    # from the list of cells that are touching, which are macrophages, i.e. which macrophages are touching cells
    afs[frame]['mac_touch'] = set(touching_cells[frame]) & set(afs[frame]['macrophages'])
    # from the list of cells that are neighbors, which are macrophages, i.e. which macrophages are near cells
    afs[frame]['mac_neighbor'] = list(neighbors[frame].keys())
    # from the touching list, which are dying
    afs[frame]['dead_touch'] = set(touching_cells[frame]) & set(afs[frame]['apoptotic'])
    # from the neighbor list, which are dying
    afs[frame]['dead_neighbor'] = set(touching_cells[frame]) & set(afs[frame]['apoptotic'])
    # of the macrophages touching cells, which are touching apoptic cells?
    cells_macrophages_touch = [touching_cells[frame][x] for x in afs[frame]['mac_touch']]
    cells_macrophages_touch = [cell for cells in cells_macrophages_touch for cell in cells] #flatten list
    afs[frame]['mac_dead_touch'] = set(cells_macrophages_touch) & set(afs[frame]['apoptotic'])
    # of the macrophages neighboring cells, which are neighboring apoptic cells?
    cells_macrophages_neighbor = [neighbors[frame][x] for x in afs[frame]['mac_neighbor']]
    cells_macrophages_neighbor = [cell for cells in cells_macrophages_neighbor for cell in cells] #flatten list
    afs[frame]['mac_dead_neighbor'] = set(cells_macrophages_neighbor) & set(afs[frame]['apoptotic'])

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

# I think this is fixed now
# fix the label image (hopefully this gets fixed soon!)
# correction_dictionary, rev_corrections = get_label_correction(label_mask, tmxml)

# find touching cells
touching_cells = find_neighbors(label_mask, 2)

# find red cells
# contrast in channel 1 seems to work pretty well. Should be greater than 0.05
macrophages = find_cells_greater_than_prop(tmxml, 'SNR_CH1', .1)
#
# # find green cells
# # this one is much trickier due to the bright middle for channel 2..
# # try contrast greater SNR over 0.4
apoptotic_cells1 = find_cells_greater_than_prop(tmxml, 'SNR_CH2', .3)
apoptotic_cells2 = find_cells_greater_than_prop(tmxml, 'MEAN_INTENSITY_CH2', 30)

apoptotic_cells = []
for frame in range(len(apoptotic_cells1)):
  apoptotic_cells.append(set(apoptotic_cells1[frame]).intersection(apoptotic_cells2[frame]))

# tracks have two columns, the first is the 'source' and the second is 'target',
# which I assume means spot index at t[x] and spot index at t[x+1]

distant_neighbors = find_distant_neighbors(label_mask, macrophages, 20)
all_frame_stats = get_frame_stats(tmxml, macrophages, apoptotic_cells, touching_cells, distant_neighbors, len(label_mask))

cell_tracks = create_cell_tracks(tmxml, all_frame_stats)

with open('test_48_24.pickle', 'wb') as handle:
  pickle.dump(all_frame_stats, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('test.pickle', 'rb') as handle:
  y = pickle.load(handle)

print('test')
# for frame in range(len(all_frame_stats)):
#   chondro_reaper_rate = (len(all_frame_stats[frame]['dead_touch']) - len(all_frame_stats[frame]['mac_dead_touch'])) / \
#                             (len(all_frame_stats[frame]['cells_label']) - len(all_frame_stats[frame]['apoptotic']) - len(all_frame_stats[frame]['macrophages']))
#   macro_reaper_rate = (len(all_frame_stats[frame]['mac_dead_touch'])) / (len(all_frame_stats[frame]['macrophages']))
#   macro_reaper_rate2 = (len(all_frame_stats[frame]['mac_dead_neighbor'])) / (len(all_frame_stats[frame]['macrophages']))
#   print('\n frame ', frame)
#   print('chondro death touch: ', chondro_reaper_rate * 100)
#   print('macro death touch: ', macro_reaper_rate * 100)
#   print('macro death neighbor: ', macro_reaper_rate2 * 100)
#   if chondro_reaper_rate == 0:
#     print('macro effectiveness: inf')
#   else:
#     print('macro effectiveness: ', macro_reaper_rate / chondro_reaper_rate)
#     print('macro reach eff: ', macro_reaper_rate2 / chondro_reaper_rate)

print('hello')