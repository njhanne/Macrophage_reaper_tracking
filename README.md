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
- Currently the trackmate xml is converted into pandas which is kind of a pain to work with. The pandas file is then 
saved as a csv for import into R, which introduces even bigger proglems. R has a big problem with empty frames, lists inside
of cells, and many more. Basically pandas and R are not supposed to work with this kind of data so htere are a lot of hacky
workarounds in the code that obfuscate what I am doing and make complex what should be simple.
  - A workaround would be to not condense all the data but then the R steps will take a VERY long time and the file sizes
  will be large. I think this is a bad idea.
  - Another potential workaround is to use SQL instead of a condensed df with a foreign key linking spots to tracks. 
  This is not ideal either since many of the Python and R steps read and write in the dataframes, which would be dangerous
  on a database. We don't necessarily want these changes to be permanent.
  - Maybe xarray? I've used it before for n-dimensional image data but I'm not sure that fits here.

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
to noise ratio is awful (see left image below) and needs to be corrected for good assignment of 'apoptotic'.
- This command needs to run before anti-cropping (next step) otherwise there are edge effects near background to true 
black borders. It does a rolling average, which will make a halo effect on the artificial anti-cropped edge.

![NucView_BG_correction](/Readme_images/background_subtraction.png)

#### Anti-crop the time series images to correct for camera movement
Videos need to be anti-cropped to correct for jumps in the video. I am not sure what causes these jumps, but it totally 
messes up the tracking if you don't correct them. Most of the images have very similar jump patterns, so maybe it is someone
walking into the microscope room and slamming the door?

This is done semi-manually by noting when
and how severely the videos jump. I open each image in ImageJ and scrub through until I notice a jump. I note down the frame
number and how many pixels it jumps. The 'correct_video_jumps.py' reads the csv of corrections and saves a new tiff with 
black space that keeps the image stabilized. Kind of hard to explain, but a picture helps: 

![anti_cropping_example](/Readme_images/anti_crop_example.png)

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
Luckily the Trackmate folks fixed it by outputting to 32bit which allows for over 2 billion spots! Unfortunately the file size is larger, though.

#### Housekeeping
Each well is imaged for 48 hours, but the raw data is saved in two 24 hr long segments. This is useful for keeping the file
size down, but means we have to have a way to link the two image files together. The way I get around this, without creating
a 48hr long image file, is by taking the last 24hr frame and the first 48hr frame and putting them together as a 2 frame image.
This image is then re-run through Trackmate with the 'link_videos_trackmate.py' ImageJ script. Then the linking is completed
with the 'link_videos.py' script (in Python, not ImageJ).

### 4. Tracking cell touching
The trackmate xml data is imported to Python where we can analyze and summarize cell touching. The next few steps are all
performed in 'main.py' The xml is loaded using a plugin called 'trackmatexml'. None of that code is mine, but it is not
available on pip or conda so I had to just include the code directly in this project. Trackmatexml parses the xml data into
a python dictionary.

#### Touching
The corresponding labelimage is loaded with the xml and used for determining touching and nearest neighbors. The touching 
analysis is performed using [scikit-image pixel graph function](https://scikit-image.org/docs/stable/api/skimage.graph.html#skimage.graph.pixel_graph). 
It is kind of hard to explain with words, but basically it makes an [adjacency matrix](https://en.wikipedia.org/wiki/Adjacency_matrix) for all the pixel values. Since the pixel
values are equivalent to id, this can be used to determine which cells are touching.

#### Neighboring
Determining which cells are 'neighbors' (near each other but not touching) is more difficult. I have two methods in the script that each have pros and cons:

a) **Distance from centroid** just looks at the x,y coordinates of the centroid of each cell to compute the euclidean distance
between all other centroids in the time frame. If the distance is less than some chosen cutoff, they are near one another.
- This method is very fast.
- It does not take into account the actual edges of the cell. If a cell is shaped like a long sliver there could be parts
of the cell body near/touching a neighbor cell that is quite far from the centroid of the cell.

b) **Expanding the mask** grows each cell mask by a chosen number of pixels and re-runs the cell touching algorithm. This 
method needs to run for every single spot separately!
- Extraordinarily slow. Perfectly accurate to strangely shaped cells. Macrophages tend to be oddly shaped...
- It can be sped up by cropping the image around the cell of interest, but it is still super slow

![Neighbor_type_comparison](/Readme_images/distance_calc_rough.png)

Because of the time it takes to run the more accurate method I usually just use the euclidian distance. Because of the inaccuracy
of the Euclidean distance, I don't use the 'neighbor' feature in any further analysis. The 8min time between each frame 
seems to be sufficient to capture cells touching so I am not sure if the 'neighbor' info is that important, anyway.

### 5. Assign cell type
Assigning 'macrophage', 'not macrophage' and 'apoptotic' is performed for each spot based on the statitistics calculated 
by Trackmate. I set the desired feature and threshold for each image by hand, since they all have different brightnesses.
These feature names and thresholds are stored in the sample info csv file ('image_log.csv'), along with other important information about the 
image and settings to use. Generally I use the standard deviation of fluor intensity, mean fluor intensity, or contrast of fluor intensity.

![Trackmate_thresholding](/Readme_images/Trackmate_thresh.png)

The Trackmate files can be opened in Trackmate which allows you to play around with the thresholds to get an idea of what 
works best for each image. In the above triptych the left shows all segmented cells in phase contrast, the middle shows 
macrophages in red, and the right shows apoptotic cells in yellow. These images are from the Trackmate GUI.

#### Saving data
Finally the script concatenates all the xml data into lists which are saved into a pandas dataframe, and saved as a csv.
The xml data has a row for every single spot. Here we have a row for every single track, and a bunch of lists that contain
info for each spot saved into each cell in the row. It's ugly and pandas is not meant to be used this way, but ultimately 
it works. 

As an example the xml may look like this:
|spot_id|frame|track_id|macrophage_bool|many_more_cols....|
|----|----|----|----|----|
|0 |20 |1| 0| ...|
|1 |21 |1|1| ...|
|2 |22 |1|1| ...|
|3 |23 |1|0| ...|

but the pandas df will look like this:
|spot_id|frame|track_id|macrophage_bool|etc....|
|----|----|----|----|----|
|(0,1,2,3) |(20,21,22,23) |1| (0,1,1,0)| ...|

Before saving the script also calculates confluency two ways: as the avg % area of the frame that is spots, and as the
avg number of cells per frame. This is saved in 'image_log_out.csv'

After the 24hr and 48hr images are run together, the script will also automatically link the two pandas results together 
and save a combined csv with the data for the full 48 hr run. Further analysis is now completed in R.

### 6. Statistical analysis (and more housekeeping)
All the R analysis and housekeeping is done with the 'apoptosis_video_analysis.R' script.
The individual csv files are loaded individually and then combined into a large dataframe. Since the previous steps were all
performed in batches, the dataframes were also created in batches. The beginning steps of the R code load these batches
individually and combine them to df_all.

After they've all been loaded the dataframe is culled of parents. In the video files, many of the cells, especially the 
macrophages pile on top of one another. Cellpose has a very hard time differentiating these individual cells when they are 
on top of one another and when they seperate Trackmate treats these events as mitosis. The way I deal with this is with the
assumption that there is no mitosis in these samples. Obviously this is not true, but it is minimal compared to the amount
of times cells pass over one another. In the script, all parent cells are deleted and their parameters are copied over to their children.
The image below I hope makes more clear what this part of the code is doing.

<img src="/Readme_images/parents_culling.png" width=50% height=50%>

As a consequence of 'parent culling' we also need to handle children inheriting cell assignment (macrophage, not_macrophage, apoptotic).

<img src="/Readme_images/child_features.png" width=50% height=50%>

#### Re-assign cell features
From here we can perform the actual analysis. I change the assignment of cell features from spot to track. As shown with the 
'macrophage_bool' in the example table above, the features are saved for each spot. In R I set rules that reassign these 
features to the entire track based on conditions. I say a cell track is a macrophage if 50% of the spots are macrophage_bool
positive, for example. This is needed because the fluors are not always above threshold, and the apoptotic cells only fluoresce
for a bit before they ball up and float off the plate.

#### Compiling, graphing, stats
Here the df is concatenated even further and then changed from wide to tall for easier graphing and analysis. Theres a lot
here and a lot of it is under active development so I'll leave this as is to be flushed out further when it's more complete.
