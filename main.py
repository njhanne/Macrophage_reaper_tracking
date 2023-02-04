import numpy as np
from skimage.graph import pixel_graph
import pandas as pd
import matplotlib.pyplot as plt

from tifffile import imread
from pytrackmate import trackmate_peak_import # if this isn't working check the version vs github version...

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

# image_path = filedialog.askopenfilename(title='Select directory of your images to be segmented')
#
# label_mask = imread(image_path)
# touching_cells = find_neighbors(label_mask, 2)

xml_track_path = filedialog.askopenfilename(title='Select directory of your images to be segmented')

spots = trackmate_peak_import(xml_track_path) # unfortunately I think this is deprecated
spots.head()


print('hello')