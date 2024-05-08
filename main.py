import numpy as np
import numpy.ma
from skimage.graph import pixel_graph
from skimage import segmentation
from scipy.spatial import distance

import pandas as pd
import matplotlib
matplotlib.use("qt5agg")
import matplotlib.pyplot as plt
# import pickle

from tifffile import imread
# from pytrackmate import trackmate_peak_import # if this isn't working check the version vs github version...

import sys
sys.path.append("./DirFileHelpers")
from find_all_files import find_all_filepaths
from trackmatexml import TrackmateXML # I didn't make this!
# https://github.com/rharkes/pyTrackMateXML
# requires lxml and (from pip) version-parser

from pathlib import Path
import tkinter as tk
from tkinter import filedialog
from collections import defaultdict
import time

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
  # https://scikit-image.org/docs/stable/api/skimage.graph.html#skimage.graph.pixel_graph

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
      # label mask and xml spot id are 1 off from each other
     pairs[int(pair[0] -1)].append(int(pair[1] - 1))
  return pairs


def find_distant_neighbors(label_mask, distance, target=None):
  # 0.65 um xy pixel
  distance_pixel = distance / 0.65

  # this works in 3D but I don't think we want it to since it is time not z, should loop over each t slice
  matches = defaultdict(dict)
  for frame in range(label_mask.shape[0]):
    if target != None:
      frame_cells = np.array(target[frame]) + 1 # add 1 for the labelimage xml offset

    else:
      frame_cells = np.unique(label_mask[frame,:,:])
      frame_cells = frame_cells[~np.isin(frame_cells, [0])] # gets rid of background

    for cell in frame_cells:
      target_mask = np.where(np.isin(label_mask[frame, :, :], cell), label_mask[frame, :, :], 0)
      expanded_mask = segmentation.expand_labels(target_mask, distance_pixel)
      matched_cells = np.unique(np.ma.masked_where(expanded_mask != cell, label_mask[frame,:,:])).data
      if matched_cells.size > 3: # this is so dumb but I think needed
        matched_cells = matched_cells[(matched_cells != 0) & (matched_cells != cell)]
        matches[frame][cell - 1] = matched_cells - 1
  return matches


def find_euclidian_neighbors(tmxml, touching_cells, thresh_distance=50):
  # 0.65 um xy pixel
  distance_pixel = thresh_distance / 0.65

  matches = defaultdict(dict)
  # this doesn't maintain order and won't work
  # [x_col, y_col, frame_col] = np.where(np.isin(tmxml.spotheader, ['POSITION_X', 'POSITION_Y', 'FRAME']))[0]
  # this does :/
  frame_col = tmxml.spotheader.index('FRAME')
  x_col = tmxml.spotheader.index('POSITION_X')
  y_col = tmxml.spotheader.index('POSITION_Y')
  slice_indices = [0, x_col, y_col] # id is always 1st col
  for frame in range(int(max(tmxml.spots[:,frame_col]))):
    centroid_list = tmxml.spots[tmxml.spots[:,frame_col] == frame][:,slice_indices]
    # calcs distance between all centroid pairs!
    distance_matrix = distance.cdist(centroid_list[:,[1,2]], centroid_list[:,[1,2]])
    # get the ones below distance
    matches_list = np.argwhere((distance_matrix <= distance_pixel) & (distance_matrix > 0))
    # match the row/column to the spot id
    matches_list = np.take(centroid_list[:,0], matches_list) # another banger function! love this!!
    # find neighbors that are already in touching list and drop them
    touching_list = np.asarray([(k, i) for k, v in touching_cells[frame].items() for i in v])
    idx, = np.where((matches_list == touching_list[:,None]).all(axis=-1).any(0)) # https://stackoverflow.com/a/75677754
    matches_list = np.delete(matches_list, idx, axis=0)
    # group them up
    matches_list = pd.DataFrame(matches_list).groupby([0]).agg(lambda x: list(x))
    matches[frame] = matches_list.to_dict()[1]
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
  child_track_alt = get_all_tracks(tmxml, [], False, True)
  cell_id = 0
  all_cells = []
  for track in tracks: # loop through all the tracks
    for cell in tracks[track]: # each track can have multiple 'splits' which are still called tracks...
      cell_track = {}
      cell_track['original_track_id'] = track # for making LUT later
      cell_track['cell_id'] = cell_id + cell.cell # make a new unique track/cell id for all the cells in the video
      cell_track['spot_ids'] = list(cell.spotids) # get all the spots in the track
      cell_track['frames'] = list(tmxml.getproperty(cell_track['spot_ids'], 'FRAME')) # get the frame # of each spot
      if cell.parent != 0:
        cell_track['parent_id'] = cell_id + cell.parent
        cell_track['independent_spot_ids'] = list(child_track_alt[track][cell.cell - 1].spotids)
        cell_track['independent_frames'] = list(tmxml.getproperty(cell_track['independent_spot_ids'], 'FRAME')) # get the frame # of each spot
      else:
        cell_track['parent_id'] = None
        cell_track['independent_spot_ids'] = None
        cell_track['independent_frames'] = None
      all_cells.append(cell_track)
    cell_id += len(tracks[track])
  return pd.DataFrame(all_cells) # very convenient we can convert list/dict to dataframe, good for csv export


def make_cell_spot_LUT(cell_tracks, tmxml):
  # create a dictionary of spot-id 'keys' for cell_track id 'values'
  LUT = defaultdict(dict)
  # you aren't really supposed to iterate on pandas objects but this is so much easier with the 'list of keys' issue
  for i, cell in cell_tracks[cell_tracks['parent_id'].isna()].iterrows():
    tmp_id = cell.cell_id
    for spot_id in cell.spot_ids:
      LUT[spot_id] = tmp_id
  # now for the children ones
  # setdefault command will only update if the key hasn't been used, which preserves the 'parent' ones instead
  # combined with the duplicate_split argument for get_all_tracks function this should make it where parent cells are the
  # 'default' cell and children are only considered when 'independent'. I think this is what we want for now.
  for i, cell in cell_tracks[cell_tracks['parent_id'].notna()].iterrows():
    tmp_id = cell.cell_id
    for spot_id in cell.spot_ids:
      LUT.setdefault(spot_id,  tmp_id)
  return LUT


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


# End xml helpers

# Analysis Helpers
def get_frame_stats(tmxml, macrophages, apoptotic_cells, touching_cells, n_frames, neighbors=None):
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

    if neighbors is not None:
      afs[frame]['mac_touch'] = set(touching_cells[frame]) & set(afs[frame]['macrophages'])
      # from the list of cells that are neighbors, which are macrophages, i.e. which macrophages are near cells
      afs[frame]['mac_neighbor'] = set(neighbors[frame]) & set(macrophages[frame])
      # from the touching list, which are dying
      afs[frame]['dead_touch'] = set(touching_cells[frame]) & set(afs[frame]['apoptotic'])
      # from the neighbor list, which are dying
      afs[frame]['dead_mac_neighbor'] = set(neighbors[frame]) & set(afs[frame]['apoptotic'])
      # of the macrophages touching cells, which are touching apoptic cells?
      cells_macrophages_touch = [touching_cells[frame][x] for x in afs[frame]['mac_touch']]
      cells_macrophages_touch = [cell for cells in cells_macrophages_touch for cell in cells] #flatten list
      afs[frame]['mac_dead_touch'] = set(cells_macrophages_touch) & set(afs[frame]['apoptotic'])
      # of the macrophages neighboring cells, which are neighboring apoptic cells?
      cells_macrophages_neighbor = [neighbors[frame][x] for x in afs[frame]['mac_neighbor']]
      cells_macrophages_neighbor = [cell for cells in cells_macrophages_neighbor for cell in cells] #flatten list
      afs[frame]['mac_dead_neighbor'] = set(cells_macrophages_neighbor) & set(afs[frame]['apoptotic'])
  return(afs)


def temp_cell_tracks_builder(cell_dict, LUT, type, other_dict=None):
  if type == 'apoptotic':
    tall_dict = [(i, k) for i,j in enumerate(list(cell_dict.values())) for k in j]
    col_names = ["cell_id", 'apoptotic_frames', 'apoptotic_spots']
  elif type == 'macrophages':
    tall_dict = [(i, k) for i,j in enumerate(list(cell_dict.values())) for k in j]
    col_names = ["cell_id", 'macrophage_frames', 'macrophage_spots']
  elif type == 'touching':
    tall_dict = [(i, k, l) for i, j in enumerate(list(cell_dict.values())) for k, l in j.items()]
    col_names = ["cell_id", 'touching_frames', 'touching_cells']
  elif type == 'neighboring':
    tall_dict = [(i, k, l) for i, j in enumerate(list(cell_dict.values())) for k, l in j.items()]
    col_names = ["cell_id", 'neighbor_frames', 'neighbor_cells']
    touch_tall_dict = [(i, k, l) for i, j in enumerate(list(other_dict.values())) for k, l in j.items()]

  temp_df = pd.DataFrame(tall_dict)  # convert to df
  if type == 'touching' or type == 'neighboring':
    # column '2' contains lists of all touching cells. This fn will 'flatten' and replicate the other columns.
    temp_df = temp_df.explode(2)  # this is amazing! I am so happy this function exists
    if type == 'neighboring':
      # this will remove all the neighbors that are actually touching (duplicate from other list)
      touch_temp_df = pd.DataFrame(touch_tall_dict)  # convert to df
      touch_temp_df = touch_temp_df.explode(2)
      exclude_list = temp_df.merge(touch_temp_df, on=[0,1,2], how='outer', indicator=True)
      temp_df = touch_temp_df.loc[exclude_list._merge == 'left_only']
    temp_df[1] = temp_df[1].map(LUT)
    temp_df[1] = pd.to_numeric(temp_df[1], downcast='integer')
    temp_df[2] = temp_df[2].map(LUT)
    temp_df[2] = pd.to_numeric(temp_df[2], downcast='integer')
    temp_df = temp_df.dropna()
    temp_df = temp_df.groupby([1]).agg(lambda x: list(x)) # group them up in lists by cell_id
    temp_df = temp_df.rename_axis("cell_id").reset_index() # get rid of 'rownames'
    temp_df.columns = col_names # set new colnames
  else:
    temp_df['cell_id'] = temp_df[1].map(LUT) #convert spots to track-cells
    temp_df['cell_id'] = pd.to_numeric(temp_df['cell_id'], downcast='integer') # get rid of 'series' datatype
    temp_df = temp_df.dropna()
    temp_df = temp_df.groupby('cell_id').agg(lambda x: list(x)) # group them up in lists by cell_id
    temp_df = temp_df.rename_axis("cell_id").reset_index() # get rid of 'rownames'
    temp_df.columns = col_names # set new colnames
  return temp_df


def get_cell_stats(tmxml, cell_tracks, cell_spot_LUT, macrophages, apoptotic_cells, touching_cells, n_frames, all_frame_stats, distant_neighbors=None):
  # apoptotic spots & frames
  temp_cell_tracks = temp_cell_tracks_builder(apoptotic_cells, cell_spot_LUT, 'apoptotic')
  cell_tracks = cell_tracks.merge(temp_cell_tracks, on='cell_id', how='outer') # outer merge into cell_track df

  # macrophage spots & frames
  temp_cell_tracks = temp_cell_tracks_builder(macrophages, cell_spot_LUT, 'macrophages')
  cell_tracks = cell_tracks.merge(temp_cell_tracks, on='cell_id', how='outer')

  # touching cells
  temp_cell_tracks = temp_cell_tracks_builder(touching_cells, cell_spot_LUT, 'touching')
  cell_tracks = cell_tracks.merge(temp_cell_tracks, on='cell_id', how='outer')

  # neighboring cells (this doesn't seem to work right now)
  if distant_neighbors is not None:
    temp_cell_tracks = temp_cell_tracks_builder(distant_neighbors, cell_spot_LUT, 'neighboring', touching_cells)
    cell_tracks = cell_tracks.merge(temp_cell_tracks, on='cell_id', how='outer') # outer merge into cell_track df

  return cell_tracks

# End analysis helpers


### MAIN ###
## First find all needed directories and load them all up
data_dir = (Path.cwd() / 'data').resolve()
results_dir = (data_dir / 'results' / 'tracks_csv').resolve()
xml_directory = (data_dir / 'processed' / 'stabilized_tiffs').resolve()
lblimgs_directory = (data_dir / 'processed' / 'label_images').resolve()

## Get sample info
# I use a csv file with info on all the images to save/load sample info
sample_info = pd.read_csv((data_dir / 'image_log.csv').resolve())
# set default values where they may not be set
sample_info['mac_filt'] = sample_info['mac_filt'].fillna('MEAN_INTENSITY_CH1')
sample_info['mac_filt_val'] = sample_info['mac_filt_val'].fillna(50.0)
sample_info['ap_filt1'] = sample_info['ap_filt1'].fillna('MEAN_INTENSITY_CH2')
sample_info['ap_filt1_val'] = sample_info['ap_filt1_val'].fillna(25.0)

## find all the Trackmate xml and label mask files
xml_dirs, xml_paths = find_all_filepaths(xml_directory, '.xml')
lmi_dirs, lmi_paths = find_all_filepaths(lblimgs_directory, '.tiff')

## actually loop through and process the trackmate results
for sample_name in sample_info['new_filename'].unique():
  this_sample_info = sample_info.loc[sample_info['new_filename'] == sample_name]
  xml_path = [xp for xp in xml_paths if xp.parts[-1].startswith(sample_name+'_')]
  lmi_path = [lp for lp in lmi_paths if lp.parts[-1].startswith('LblImg_'+sample_name+'_')]

  if len(xml_path) != 0 and len(lmi_path) != 0:
    print('Analyzing ' + sample_name)
    label_mask = imread(lmi_path[0]) # import label mask

    tmxml = TrackmateXML()
    tmxml.loadfile(xml_path[0]) # import TrackMate XML

    ## find touching cells
    print('Finding touching cells...')
    touching_cells = find_neighbors(label_mask, 2)

    ## find macrophages
    print('Finding macrophages...')
    macrophages = find_cells_greater_than_prop(tmxml, this_sample_info['mac_filt'].values[0], this_sample_info['mac_filt_val'].values[0])

    ## find apoptotic cells
    print('Finding apoptotic cells...')
    apoptotic_cells = find_cells_greater_than_prop(tmxml, this_sample_info['ap_filt1'].values[0], this_sample_info['ap_filt1_val'].values[0])
    # apoptotic_cells2 = find_cells_greater_than_prop(tmxml, 'MEAN_INTENSITY_CH2', 30)
    #
    # apoptotic_cells = []
    # for frame in range(len(apoptotic_cells1)):
    #   apoptotic_cells.append(set(apoptotic_cells1[frame]).intersection(apoptotic_cells2[frame]))

    ## Find neighboring cells
    print('Finding neighboring cells...')
    # This method is the most accurate but it is VERY slow since it has to separately look at each labelled cell and expand it
    # distant_neighbors = find_distant_neighbors(label_mask, 30)

    # This method just uses the centroid position to calculate euclidian distance between all the centroids, then matches up
    # the ones below a set threshold. It should be exponentially faster than the more accurate method
    distant_neighbors = find_euclidian_neighbors(tmxml, touching_cells, 50)

    ## Build and save results
    print('Building csv...')
    all_frame_stats = get_frame_stats(tmxml, macrophages, apoptotic_cells, touching_cells, len(label_mask), distant_neighbors)

    cell_tracks = create_cell_tracks(tmxml, all_frame_stats)
    cell_spot_LUT = make_cell_spot_LUT(cell_tracks, tmxml)
    all_cell_stats = get_cell_stats(tmxml, cell_tracks, cell_spot_LUT, macrophages, apoptotic_cells, touching_cells,  len(label_mask), all_frame_stats, distant_neighbors)

    print('Saving results...')
    save_name = (results_dir / (sample_name + '.csv')).resolve()
    all_cell_stats.to_csv(save_name)