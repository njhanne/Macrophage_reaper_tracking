from cellpose import models, io

from cellpose.contrib.distributed_segmentation import distributed_eval
from cellpose.contrib.distributed_segmentation import numpy_array_to_zarr
import numpy as np

from pathlib import Path
import time

import sys
sys.path.append("./DirFileHelpers")
from DirFileHelpers.find_all_files import find_all_filepaths


### main ###
# load in all the needed files and paths
data_dir = (Path.cwd() / 'data' / 'processed').resolve()
# get all the stacks we created from rgb_to_stack.py
zarr_directory = (data_dir / 'BG_corrected' / 'zarr').resolve()
image_directory = (data_dir / 'BG_corrected' / 'run').resolve()
image_dirs, images = find_all_filepaths(image_directory, '.tif')
# find the cellpose models that Charlie helped train
model_path = (Path.cwd().parent.parent / '.cellpose' / 'models' / 'NH_LB_LC4_HL1_3chan').resolve()
# where to save the cellpose output
output_directory = (data_dir / 'cellpose_output').resolve()

### IMPORTANT! ###
# Change this to False if you did not install the cuda pytorch version when installing cellpose!
cuda = True

# parameterize cellpose however you like
# cellpose needs c z y x
model_kwargs = {'gpu': cuda, 'pretrained_model': str(model_path)}  # can also use 'pretrained_model'
eval_kwargs = {'stitch_threshold': 1.0,
               'z_axis': 0,
               'channel_axis': 1,
               'normalize': {'lowhigh': None, 'percentile': [1.0, 99.0],
                             'normalize': True, 'norm3D': True, 'sharpen_radius': 0.0,
                             'smooth_radius': 0.0, 'tile_norm_blocksize': 0.0,
                             'tile_norm_smooth3D': 0.0, 'invert': False},
               }




# define compute resources for local workstation
cluster_kwargs = {
  'n_workers': 1,  # if you only have 1 gpu, then 1 worker is the right choice
  'ncpus': 8,
  'memory_limit': '64GB',
  'threads_per_worker': 1,
}


def stack_channels(image, crop):
    return np.stack((data_zarr[crop], image), axis=1)  # second channel is also a zarr array
preprocessing_steps = [(stack_channels, {}), ]

# for image_path in images:
    # filename = Path(image_path).name
filename = Path(images[0]).name
print('analyzing image: ', filename)
image = io.imread(str(images[0]))
chunk_dims = (1, image.shape[1], image.shape[2], image.shape[3])
data_zarr = numpy_array_to_zarr(str(zarr_directory), image, chunks=chunk_dims)
del image
print('image loaded')

chunk_dims = (1, image.shape[1], image.shape[2], image.shape[3])


# io.logger_setup()

print('running')
tic = time.time()
# masks, flows, styles = model.eval(image, channel_axis=1, z_axis=0, stitch_threshold=1.0, normalize={'lowhigh': None, 'percentile': [1.0, 99.0], 'normalize': True, 'norm3D': True, 'sharpen_radius': 0.0, 'smooth_radius': 0.0, 'tile_norm_blocksize': 0.0, 'tile_norm_smooth3D': 0.0, 'invert': False})

# run segmentation
# outputs:
#     segments: zarr array containing labels
#     boxes: list of bounding boxes around all labels (very useful for navigating big data)
segments, boxes = distributed_eval(
  input_zarr=data_zarr[:,2,:,:],
  write_path=str(output_directory),
  blocksize=chunk_dims,
  model_kwargs=model_kwargs,
  eval_kwargs=eval_kwargs,
  cluster_kwargs=cluster_kwargs,
)
toc = time.time()
print((toc-tic)/60)


print('saving')
# io.save_masks(image, masks, flows, filename, tif=True,  savedir=str(output_directory), save_txt=False)
tec=time.time()
print((tec - toc) / 60)