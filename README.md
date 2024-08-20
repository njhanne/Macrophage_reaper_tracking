# Macrophage_reaper_tracking
The purpose of this project is to track the movement and touching of cells *in vitro*, to test the hypothesis that macrophages
can induce apoptosis in other cells by touching them. I call this process 'reaping' because the macrophages, like the grim
reaper, touch to induce death.

Macrophages and chondrocytes were grown either in monoculture or together in coculture. The cultures were imaged every
8 min for 48 hours. The time-series images have three channels: 1) phase-contrast image that captures all cells; 2) a red 
fluor channel for the macrophages which express membrane bound tdTomato; and 3) a green fluor for
[NucView](https://biotium.com/technology/nucview-caspase-3-substrates/) which only fluoresces when cells undergo apoptosis. 

These scripts:
1) convert the image files to the appropriate format
2) segment the cells
3) track the cells' movement
4) track the cells' touching and neighbors
5) assign cell type based on fluorescence intensity
6) perform statistical analysis

These are still being developed so the code may not be too clean.

## Overall strategy & notes
- The tracking is done with [Trackmate](https://imagej.net/plugins/trackmate/), a plugin for imageJ. It works 
well! It is mainly limited by the quality of segmentation of cells.
- Segmentation is done in [Cellpose](https://github.com/MouseLand/cellpose). We have a pretty good training set for Cellpose, and segmentation works quite well 
for sparse cells but the high confluence plates are more difficult. 
  - At this point, with plenty of 'in the loop'
  hand training, only the most confluent plate doesn't segment well. The rest are fine! 
- Running cellpose through trackmate is slow. Cellpose in general is slow for these huge files. Luckily I  can set it and
forget it O/N and it does fine.
  - I saw a talk by Valentina Pedoia where she recommended using the ['Yolo'](https://github.com/ultralytics/ultralytics) tracker as it is much faster than cellpose.
  It is the kind of tracking used for car vision, it draws rectangles around objects. Unfortunately, we need the segmentations 
  because that's how we get intensity of our fluors... So this won't work for this project.
    - There is probably a way to export the rectangles to trackmate for segmentation, which could be faster... 
    but not faster than development time.
  - Running cellpose outside of Trackmate isn't ideal. It treats the time slices as z slices (3D segmentation). 
  You can set 3D to false in cellpose, but it really struggles to open these files. It's easier to run it in Trackmate, 
  and not much slower.
  - UCSF is getting a supercompute cluster with GPUs, but not too soon.

3. Run cellpose and Trackmate. Trackmate fortunately runs cellpose on t-stacks automatically
4. The trackmate data is input back to python where we can analyze and summarize the cell activities. At some future point
it may (unfortunately) be a good idea to switch the data into some kind of sql database since the ids and tracks are a bit
convoluted to deal with. Let's see how far we can get without adding that complexity.
   - The labelimages are automatically saved as 16bit which 'only' allows for 65,000 segmentated objects. This problem has
   been resolved with an update.
5. Python script that creates a 2 frame 'video' of the last 24hr and first 48hr frames. This needs to be run through 
Trackmate to create links between the two video files.
6. R script to do the remaining analysis.

## Work order

### 1. Formatting and cleaning images
#### Convert nd2 to 8bit tiff
Python script 'convert_rename_nd2.py'. This script converts the Nikon proprietary image format into 8bit tiffs. It also
rearranges the order of the channels so that they are macrophage(red):apoptosis(green):phase-contrast(grey). 
  - The .nd2 files are 16bit, but this is overkill for this segmentation. Also it makes the files twice as big :)
  - Setting the order of the channels needs to be changed manually in the script with the 'rearrange' variable. The images 
  were taken in 3 batches and each batch has a different order of channels.
  - There is a imagej script called 'convert_8bit.ijm' that will convert 16bit to 8bit if you want, too.

#### Background correction for NucView channel
ImageJ script 'BG_correct_TUNEL.ijm' runs the background correction command on the apoptosis fluor. This channel is signal
to noise ratio is awful and needs to be corrected for good assignment of 'apoptotic'.
- This command needs to run before anti-cropping (next step) otherwise there are edge effects near background to true 
black borders. It does a rolling average, which will make a halo effect on the artificial anti-cropped edge.

#### Anti-crop the time series images to correct for camera movement
Videos need to be anti-cropped to correct for jumps in the video. This is done semi-manually by noting when
and how severely the videos jump. This is done with 'correct_video_jumps.py'

### 2&3.Segment and track the cells w/ Cellpose and Trackmate
As discussed above, these steps are done by Cellpose and Trackmate, respectively. Since Cellpose runs inside of Trackmate
I lump the two steps together here. All of the segmentation and tracking is performed on the phase-contrast channel since 
it has all of the different cells visible.

#### Training Cellpose (optional)
I've already completed this step so it doesn't need to be re-done, unless you want to further refine the model.
A previous postdoc (Leila Bagheri) and a dental student (Karan Patel) spent a lot of time hand-tracing cells to train a model.
During that time a new version of Cellpose was released (v3 I think?) that allows for 'human in the loop' training where you
can manually adjust segmentations and use them to re-train the model before segmenting another image. It works wonderfully.
I used the hand drawn traces that Leila and Karan had previously generated to refine the live cell (LC) model that comes with 
Cellpose, then did a few iterations of 'in the loop' to get it working even better. All of these steps are performed inside
the Cellpose GUI.
- The'ML_training_gen.py' script will automatically crop out a specified number of random timepoints from a directory of images.
These can be used to get a good spread of training data to refine the model.

#### Running Cellpose and Trackmate
Trackmate is run through ImageJ, but it only will process one image at a time. I made a script 'run_trackmate.py' which runs
inside ImageJ to process a full directory of images. I usually run them in batches since the file size is too large to
have them all available at the same time.

- The script has a number of settings that also control Cellpose segmentation. The one that makes the biggest difference is
cell diameter and is a good place to start troubleshooting.
- The script also has the settings that control Trackmate tracking algorithm. These settings are pretty well documented on their
[website](https://imagej.net/media/plugins/trackmate/trackmate-manual.pdf).
- Trackmate assigns each Cellpose segmented object at each time point as a 'spot'. It then uses the tracking algorithm 
to assign spots to a 'track' if the spots are near one another between time points (it's more complicated than that, but
for explanatory sake...). The settings in the script that control the tracking algorithm relate to this process of forming links.

This Trackmate script will output two files – an xml file with all the Trackmate data on it (detailed below), and a label
image file where pixel values correspond to cell (spot) id.
- The xml file has statistics for every single spot (area, position, circularity, rel. velocity, ch1 brightness, and more)
and also a track_id if the spot has been linked to other spots across time.
- The label image file is a greyscale image of each segmentation with the pixel value corresponding to the assigned spot
id by Trackmate. Funny aside: this labelimage used to be 16bits which only allows for values up to 65,536. Many of these
videos have more 'spots' than that, and the spot_id would just roll over to 0 and restart, which would screw everything up.
Luckily the Trackmate folks fixed it by outputting to 32bit which allows for over 2billion spots! Unfortunately the file size is larger, though.

### 4. Tracking cell touching

### 5. Assign cell type

### 5.5 Housekeeping

### 6. Statistical analysis