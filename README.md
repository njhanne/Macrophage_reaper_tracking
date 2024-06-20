# CellTracker_LB

### General strategy & notes
- The actual tracking is done with 'trackmate', a plugin for imageJ. It works quite well! It is mainly limited by the 
quality of our segmentation of cells.
- Segmentation is done in Cellpose. We have a pretty good training set for Cellpose, and segmentation works quite well 
for sparse cells but the high confluence plates are difficult (probably impossible).
  - The segmentation
- Running cellpose through trackmate is very slow. Cellpose in general is very slow for these huge files.
  - I saw a talk by Valentina Pedoia where she recommended using the 'Yolo' tracker as it is much faster than cellpose.
  It is the kind of tracking used for car vision, it draws rectangles around objects. For this project it would probably
  be sufficient to find the cells since we aren't currently using the segmentations, but it does not output in a format
  that trackmate uses. I'm sure I could figure out a way to get them to work together but that time would probably
  be greater than the time to just run Cellpose ha! I could just convert the rectangles to a label mask or something...
  but again... time.
    - Actually we do need the segmentations because that's how we get intensity of our fluors...
  - Running cellpose separately isn't ideal. It treates the time slices as z slices (3D segmentation). You can set 3D to false
  in cellpose, but it really struggles to open these files. It's easier to run it in Trackmate, and not much slower.
  - UCSF is getting a supercompute cluster with GPUs, but not too soon.

1. ImageJ script (BG_correct_TUNEL.ijm) that does background correction on the apoptosis fluor. Needs to be run before anti-cropping (step 2) 
otherwise there are edge effects near background to true black borders. 
2. Videos need to be anti-cropped to correct for jumps in the video. Currently this is done semi-manually by noting when
and how severely the videos jump. This is done with 'correct_video_jumps.py'
3. Run cellpose and Trackmate. Trackmate fortunately runs cellpose on t-stacks automatically
4. The trackmate data is input back to python where we can analyze and summarize the cell activities. At some future point
it may (unfortunately) be a good idea to switch the data into some kind of sql database since the ids and tracks are a bit
convoluted to deal with. Let's see how far we can get without adding that complexity.
   - The labelimages are automatically saved as 16bit which 'only' allows for 65,000 segmentated objects. This problem has
   been resolved with an update.
5. Python script that creates a 2 frame 'video' of the last 24hr and first 48hr frames. This needs to be run through 
Trackmate to create links between the two video files.
6. R script to do the remaining analysis.