import numpy as np
import numpy.ma
from skimage.graph import pixel_graph
from skimage import segmentation
from scipy.spatial import distance

import pandas as pd
import matplotlib
import math
matplotlib.use("qt5agg")
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

from pathlib import Path
from collections import defaultdict

import time
from joblib import Parallel, delayed


### Neighbor Finding Helpers
def find_neighbors(label_mask, connectivity = 2, frames_to_add=0, debug = False):
  # https://stackoverflow.com/questions/72452267/finding-identity-of-touching-labels-objects-masks-in-images-using-python
  all_pairs = defaultdict(dict)
  if len(label_mask.shape) == 2:
    # I think I can expand the radius of the mask (i.e. what is considered 'touching') but making a 2D bool structure
    # for connectivity. i.e. skimage.morphology.disk(4) would be a disk of radius 4
    # no, this is not how this function works or what it does. It only finds edges that touch, not overlap
    # I should use the overlapping code from the GM130 work.
    all_pairs[0] = find_touching_cells(label_mask, connectivity = 2)
  else:
    if debug:
      for slice in range(len(label_mask)):
        all_pairs[slice + frames_to_add] = find_touching_cells(label_mask[slice, :, :], connectivity = connectivity)
    else:
      all_pairs_temp = Parallel(n_jobs=6)(delayed(find_touching_cells)(label_mask[slice, :, :], connectivity = connectivity) for slice in range(0, len(label_mask)))
      all_pairs = {k+frames_to_add: all_pairs_temp[k] for k in range(len(all_pairs_temp))}
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


def find_euclidian_neighbors(tmxml, touching_cells, thresh_distance=50, frames_to_add=0):
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
    touching_list = np.asarray([(k, i) for k, v in touching_cells[frame + frames_to_add].items() for i in v])
    idx, = np.where((matches_list == touching_list[:,None]).all(axis=-1).any(0)) # https://stackoverflow.com/a/75677754
    matches_list = np.delete(matches_list, idx, axis=0)
    # group them up
    matches_list = pd.DataFrame(matches_list).groupby([0]).agg(lambda x: list(x))
    matches[frame + frames_to_add] = matches_list.to_dict()[1]
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


def create_cell_tracks(tmxml, frames_to_add=0):
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
      cell_track['frames'] = [x+frames_to_add for x in cell_track['frames']]
      if cell.parent != 0:
        cell_track['parent_id'] = cell_id + cell.parent
        cell_track['independent_spot_ids'] = list(child_track_alt[track][cell.cell - 1].spotids)
        cell_track['independent_frames'] = list(tmxml.getproperty(cell_track['independent_spot_ids'], 'FRAME')) # get the frame # of each spot
        cell_track['independent_frames'] = [x + frames_to_add for x in cell_track['independent_frames']]
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


def find_cells_greater_than_prop(tmxml, cell_propertyname, limit, frames_to_add=0):
  cells  = defaultdict(list)
  property_col = tmxml.spotheader.index(cell_propertyname)
  frame_c = tmxml.spotheader.index('FRAME')
  cell_property = tmxml.spots[:,property_col]
  cell_list = (i for i, x in enumerate(cell_property) if x >= limit)
  for i in cell_list:
    # gets the label id instead of the row number (i)
    cells[tmxml.spots[i,frame_c]+frames_to_add].append(tmxml.spots[i,0])
  return cells


# End xml helpers

# Analysis Helpers
def get_cellularity(tmxml, label_mask):
  # count all non-zero (not background) divided by 2048 x 2048 x slices
  # do non-zero because with the video correction we've introduced more background than there are in reality...
  cell_area = np.count_nonzero(label_mask) / (2048 * 2048 * label_mask.shape[0])
  # count number of spots (cells) / slices
  cell_count = tmxml.spots.shape[0] / label_mask.shape[0]
  return [cell_area, cell_count]


def get_frame_stats(tmxml, macrophages, apoptotic_cells, touching_cells, n_frames, neighbors=None, frames_to_add=0):
  afs = defaultdict(dict)
  frame_c = tmxml.spotheader.index('FRAME')
  for frame in range(n_frames):
    frame_i = frame
    frame = frame + frames_to_add
    # get the index of all cells in this frame
    frame_cells_i = np.where(tmxml.spots[:,frame_c] == frame_i)[0]
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
    tall_dict = [(k, i) for k,j in list(cell_dict.items()) for i in j]
    col_names = ["cell_id", 'apoptotic_frames', 'apoptotic_spots']
  elif type == 'macrophages':
    tall_dict = [(k, i) for k,j in list(cell_dict.items()) for i in j]
    col_names = ["cell_id", 'macrophage_frames', 'macrophage_spots']
  elif type == 'touching':
    tall_dict = [(i, k, l) for i, j in list(cell_dict.items()) for k, l in j.items()]
    col_names = ["cell_id", 'touching_frames', 'touching_cells']
  elif type == 'neighboring':
    tall_dict = [(i, k, l) for i, j in list(cell_dict.items()) for k, l in j.items()]
    col_names = ["cell_id", 'neighbor_frames', 'neighbor_cells']
    # touch_tall_dict = [(i, k, l) for i, j in list(other_dict.items()) for k, l in j.items()]

  temp_df = pd.DataFrame(tall_dict)  # convert to df
  if type == 'touching' or type == 'neighboring':
    # column '2' contains lists of all touching cells. This fn will 'flatten' and replicate the other columns.
    temp_df = temp_df.explode(2)  # this is amazing! I am so happy this function exists
    # if type == 'neighboring':
    # I've moved all this functionality into the finding neighbors code. I like the way this was done with the 'merge'
    # so I'm going to leave it in for future reference
    #   # this will remove all the neighbors that are actually touching (duplicate from other list)
    #   touch_temp_df = pd.DataFrame(touch_tall_dict)  # convert to df
    #   touch_temp_df = touch_temp_df.explode(2)
    #   exclude_list = temp_df.merge(touch_temp_df, on=[0,1,2], how='outer', indicator=True)
    #   # 'left only' is only present in temp_df, not in touch_temp_df
    #   temp_df = touch_temp_df.loc[exclude_list._merge == 'left_only']
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


def get_cell_stats(cell_tracks, cell_spot_LUT, macrophages, apoptotic_cells, touching_cells, distant_neighbors=None):
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


def analyze_xml(tmxml, label_mask, sample_info, this_sample_info, frames_to_add=0, cell_tracks_to_add=0, csv_links=None, cell_stats_24=None):
  # this is so that using debug mode doesn't crash pycharm. It doesn't work with the parallelized processing, unforunately
  debug = False

  link_results=False
  sample_name = this_sample_info['new_filename'].values[0]

  ## get cellularity
  info_rindex = sample_info.loc[sample_info['new_filename'] == sample_name].index.values[0]
  sample_info.iloc[info_rindex,-2:] = get_cellularity(tmxml, label_mask)

  # get cell track lookup table, correct if 48 hr sample
  cell_tracks = create_cell_tracks(tmxml, frames_to_add)
  cell_spot_LUT = make_cell_spot_LUT(cell_tracks, tmxml)

  # if it's the 48hr one, create the links b/w 24 and 48 hr videos for the tracks instead of spots
  if csv_links is not None and cell_stats_24 is not None:
    csv_links = link_video_tracks(cell_stats_24, cell_tracks, csv_links)
    # all these csv_link numbers are row numbers, not the track_id
    # luckily right now they are all in order so adding 1 will correct them
    csv_links_track_id = np.array(csv_links) + 1


  # here use the adjusted tracks to update the LUT, should make later steps much easier if it gets sorted out early!
  if cell_tracks_to_add != 0 and csv_links_track_id is not None:
    cell_spot_LUT_array = np.array((list(cell_spot_LUT.keys())[0:-1],list(cell_spot_LUT.values())[0:-1])).T
    # cell_spot_LUT_array = np.concatenate((cell_spot_LUT_array, cell_spot_LUT_array[:,1]), axis=1)
    # The above does not work - index slice notation returns a 1D array that can't be concatenated (DUMB!)
    # here we create a copy of the track index
    cell_spot_LUT_array = np.concatenate((cell_spot_LUT_array, cell_spot_LUT_array[:,-1:], cell_spot_LUT_array[:,-1:]), axis=1)
    # and overwrite all the ones that are in the 24hr links
    mask = np.isin(cell_spot_LUT_array[:, 2], csv_links_track_id[:,1])
    masked_mapped_vals = np.array(pd.DataFrame(cell_spot_LUT_array[:,2])[0].map({ k: v for v, k in csv_links_track_id}))
    cell_spot_LUT_array[:, 2] = np.where(mask, masked_mapped_vals, cell_spot_LUT_array[:, 2])
    # then apply the increased track number
    cell_spot_LUT_array[:, 1] = cell_spot_LUT_array[:, 1] + cell_tracks_to_add
    # and merge
    cell_spot_LUT_array[:, 1] = np.where(mask, cell_spot_LUT_array[:, 2], cell_spot_LUT_array[:, 1])
    # use this array to make two LUTs, one for the old to new tracks, one for the spots to new tracks
    cell_spot_LUT = {k: v for k, v in cell_spot_LUT_array[:,0:2]}
    new_track_LUT = {k: v for v, k in cell_spot_LUT_array[:,[1,3]]}

    # Now these adjusted LUT need to be applied to the cell tracks that have already been evaluated
    cell_tracks['cell_id'] = cell_tracks['cell_id'].map(new_track_LUT)
    cell_tracks['parent_id'] = cell_tracks['parent_id'].map(new_track_LUT)
    link_results = True

  ## find touching cells
  print('Finding touching cells...')
  tic = time.time()
  touching_cells = find_neighbors(label_mask, 2, frames_to_add, debug)
  toc = time.time()
  print(toc-tic)

  ## find macrophages
  print('Finding macrophages...')
  macrophages = find_cells_greater_than_prop(tmxml, this_sample_info['mac_filt'].values[0],
                                             this_sample_info['mac_filt_val'].values[0], frames_to_add)

  ## find apoptotic cells
  print('Finding apoptotic cells...')
  apoptotic_cells = find_cells_greater_than_prop(tmxml, this_sample_info['ap_filt1'].values[0],
                                                 this_sample_info['ap_filt1_val'].values[0], frames_to_add)
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
  distant_neighbors = find_euclidian_neighbors(tmxml, touching_cells, 50, frames_to_add)

  ## Build and save results
  print('Building csv...')
  # all_frame_stats = get_frame_stats(tmxml, macrophages, apoptotic_cells, touching_cells, len(label_mask),
  #                                   distant_neighbors, frames_to_add)


  all_cell_stats = get_cell_stats(cell_tracks, cell_spot_LUT, macrophages, apoptotic_cells, touching_cells, distant_neighbors)

  if link_results:
    compiled_cell_stats = link_cell_stats(cell_stats_24, all_cell_stats, csv_links)
  else:
    compiled_cell_stats = None

  return(all_cell_stats, sample_info, compiled_cell_stats)


def link_video_tracks(cell_stats_24, cell_stats_48, csv_links):
  # first figure out which tracks link from 24 to 48 hours
  final_24_spots = pd.DataFrame([t[-1] for t in cell_stats_24['spot_ids']])
  final_24_spots['track_index_24'] = final_24_spots.index
  first_48_spots = pd.DataFrame([t[0] for t in cell_stats_48['spot_ids']])
  first_48_spots['track_index_48'] = first_48_spots.index
  first_48_spots = first_48_spots[first_48_spots[0].duplicated() == False] # we only want the non-duplicated ones (not children)
  csv_links = pd.merge(csv_links, final_24_spots, 'inner', left_on='24hr_spots', right_on=[0])
  csv_links = pd.merge(csv_links, first_48_spots, 'inner', left_on='48hr_spots', right_on=[0])
  csv_links = csv_links[['track_index_24', 'track_index_48']]
  return csv_links


def link_cell_stats(cell_stats_24, cell_stats_48, csv_links):
  # These first two lines convert all the 'nan' values to None, which makes the logic of concatenating rows easier later
  # it doesn't get 'parent_id', but we aren't trying to combine that row!
  cell_stats_combined = cell_stats_24.where(cell_stats_24.notnull(), None)
  cell_stats_48 = cell_stats_48.where(cell_stats_48.notnull(), None)

  col_names_to_copy = ['spot_ids', 'frames', 'independent_spot_ids', 'independent_frames', 'apoptotic_frames',
                       'apoptotic_spots', 'macrophage_frames', 'macrophage_spots', 'touching_frames', 'touching_cells',
                       'neighbor_frames', 'neighbor_cells']
  for index, link_track in csv_links.iterrows():
    combined_row = cell_stats_combined.iloc[link_track.iloc[0],:]
    # new_row = cell_stats_48.iloc[link_track.iloc[1],:]
    for col in col_names_to_copy:
      og = combined_row.loc[col]
      new = cell_stats_48.loc[link_track.iloc[1],col]
      if og is None and new is None:
        combined = None
      elif og is None:
        combined = new
      elif new is None:
        combined = og
      else:
        combined = og + new
      combined_row.loc[col] = combined
    cell_stats_combined.loc[link_track.iloc[0]] = combined_row
  # append the 48hr stats to the 24, but remove the rows we copied in the awful loop above
  cell_stats_combined = pd.concat([cell_stats_24, cell_stats_48.drop(csv_links['track_index_48'])], axis=0)
  return cell_stats_combined

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
csv_dirs, csv_paths = find_all_filepaths(xml_directory, '.csv')

## actually loop through and process the trackmate results
for sample_name in sample_info['new_filename'].unique():
  pair_video = False # boolean to prevent 48 hr analysis if not loaded properly or if it doesn't exist
  loaded = False # boolean to prevent 24 hr analysis if not loaded properly

  this_sample_info = sample_info.loc[sample_info['new_filename'] == sample_name]
  # we don't want to run it on the pairs twice, only look at the 24hr ones, also skip ones we want to skip
  if this_sample_info['time'].values[0] == 24 and this_sample_info['skip'].values[0] != True:
    # get the 48 hour counterpart
    paired_sample_info = sample_info.loc[sample_info['new_filename_timeless'] == this_sample_info['new_filename_timeless'].values[0]]
    if len(paired_sample_info) > 1: # only look at ones that actually have pairs
      paired_sample_info = paired_sample_info[paired_sample_info['time'] == 48]
      pair_video = True
      paired_name = paired_sample_info['new_filename'].values[0]

    # load the xml files. We will do this first to make sure they load correctly. If not we won't run analysis...
    tmxml_path = [tm for tm in xml_paths if tm.parts[-1].startswith(this_sample_info['new_filename'].values[0] + '_')]
    tmxml_24 = TrackmateXML()
    try:
      tmxml_24.loadfile(tmxml_path[0])  # import TrackMate XML
      loaded = True
    except:
      loaded = False
      print('could not load xml for ' + this_sample_info['new_filename'].values[0])
    if loaded and pair_video:
      tmxml_path = [tm for tm in xml_paths if tm.parts[-1].startswith(paired_name + '_')]
      tmxml_48 = TrackmateXML()
      try:
        tmxml_48.loadfile(tmxml_path[0])
      except:
        pair_video = False
        print('could not load xml for ' + paired_name)

    # run the analysis
    if loaded:
      lmi_path = [lp for lp in lmi_paths if lp.parts[-1].startswith('LblImg_' + sample_name + '_')]
      if len(lmi_path) != 0:
        print('Analyzing ' + sample_name)
        label_mask = imread(lmi_path[0])  # import label mask
        cell_stats_24, sample_info, _ = analyze_xml(tmxml_24, label_mask, sample_info, this_sample_info)
        print('Saving 24hr results...')
        save_name = (results_dir / (sample_name + '.csv')).resolve()
        cell_stats_24.to_csv(save_name)

      if pair_video:
        lmi_path = [lp for lp in lmi_paths if lp.parts[-1].startswith('LblImg_' + paired_name + '_')]
        if len(lmi_path) != 0:
          print('Analyzing ' + paired_name)
          csv_link_file = [csv for csv in csv_paths if csv.parts[-1].startswith(this_sample_info['new_filename_timeless'].values[0] + '_links')]
          csv_links = pd.read_csv(csv_link_file[0])

          cell_tracks_to_add = max(cell_stats_24['cell_id'].to_numpy()) + 1
          frames_to_add = max([t[-1] for t in cell_stats_24['frames']]) + 1

          label_mask = imread(lmi_path[0])  # import label mask
          cell_stats_48, sample_info, combined_results = analyze_xml(tmxml_48, label_mask, sample_info, paired_sample_info, frames_to_add, cell_tracks_to_add, csv_links, cell_stats_24)
          print('Saving 48hr results...')
          save_name = (results_dir / (paired_name + '.csv')).resolve()
          cell_stats_48.to_csv(save_name)

          print('Saving combined results...')

          save_name = (results_dir / (this_sample_info['new_filename_timeless'].values[0] + '_combined.csv')).resolve()
          combined_results.to_csv(save_name)
sample_info.to_csv((data_dir / 'image_log_out.csv').resolve())