from cellpose import models, io
from pathlib import Path
import pandas as pd

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### main ###
info_csv_path = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/GM130_image_log.csv')
info_csv_path = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/GM130_image_log_test.csv')
sample_info = pd.read_csv(info_csv_path)

image_directory = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/new_images/mix')
image_dirs, images = find_all_filepaths(image_directory, '.tif')

# select 'GM130' or 'DAPI'
stain_type = 'GM130'

if stain_type == 'GM130':
  diam = 15
  model_path = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/cellpose_output/train/GM130/models/golgi_lab40X_GM5"
  output_directory = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/cellpose_output/golgi')
  type = 'golgi'
elif stain_type == 'DAPI':
  diam = 30
  model_path = "C:/Users/njhan/PycharmProjects/Golgi_angles/data/training/DAPI/models/Charlie_2D_HumanLoop10"
  output_directory = Path('C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/cellpose_output/nuclei')
  type = 'nuclei'
else:
  diam = None
  model_type = 'cyto2'

# DEFINE CELLPOSE MODEL
model = models.CellposeModel(gpu=True, pretrained_model=model_path)

chan = [0,0]

# or in a loop
for image in images:
  filename = Path(image).name
  sample = sample_info.loc[sample_info['old_filename'] == filename]
  if sample.size != 0:
    if sample['channel'].values[0] == type:
      print('analyzing sample', filename)
      anis = sample['z_space'].values[0]
      img = io.imread(image)
      masks, flows, styles = model.eval(img, diameter=diam, channels=chan, do_3D=True, anisotropy=anis)

      # save results so you can load in gui
      # io.masks_flows_to_seg(img, masks, flows, diam, filename, chan)
      print('saving sample', filename)
      # save results as tiff
      outname = filename + '_cp.tiff'
      io.save_masks(img, masks, flows, filename, png=False, tif=True, savedir=str(output_directory))
      print('\n')