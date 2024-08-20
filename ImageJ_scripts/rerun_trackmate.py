# https://imagej.net/plugins/trackmate/scripting/trackmate-detectors-trackers-keys
# https://github.com/trackmate-sc/TrackMate/blob/master/scripts/ExampleScript_1.py
# https://github.com/trackmate-sc/TrackMate/blob/master/scripts/TrackMate-Cellpose-script.py
# https://github.com/trackmate-sc/TrackMate/blob/master/scripts/ExampleScript_ReadTrackMateFile.py

import sys
from datetime import datetime as dt
import re
import os


from ij import IJ
from ij import WindowManager

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter

from fiji.plugin.trackmate.util import LogRecorder
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.io import TmXmlReader
from fiji.plugin.trackmate.util import TMUtils
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.action import LabelImgExporter

from java.io import File

# We have to do the following to avoid errors with UTF8 chars generated in
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')


# ------------------------------------------------------
# 	EDIT FILE PATHS BELOW.
# ------------------------------------------------------

def find_all_filepaths(directory, extension):
  filepaths = []

  for subdir, dirs, files in os.walk(directory):
    for file in files:
      # print os.path.join(subdir, file)
      filepath = os.path.normpath(subdir + os.sep + file)

      if filepath.endswith(extension):
        filepaths.append(filepath)
        print(filepath)
  return(filepaths)

# Shall we display the results each time?
show_output = False

# Channel to process?
channel_to_process = 3

# Image files to analyse.
directory_to_process = os.path.normpath("C:/Users/njhan/Box/macrophage_coculture/processed/stabilized_tiffs")
file_paths = find_all_filepaths(directory_to_process, '.xml')
print(file_paths)


# -------------------------------------------------------
# 	ACTUAL CODE.
# ------------------------------------------------------
def run(xml_file):
  # We have to feed a logger to the reader.
  logger = Logger.IJ_LOGGER

  # Logger -> content will be saved in the XML file.
  logger = LogRecorder(Logger.VOID_LOGGER)

  logger.log('TrackMate-Cellpose analysis script\n')

  dt_string = dt.now().strftime("%d/%m/%Y %H:%M:%S")
  logger.log(dt_string + '\n\n')

  # -------------------
  # Instantiate reader
  # -------------------
  xml_file = File(xml_file)
  reader = TmXmlReader(xml_file)
  if not reader.isReadingOk():
    sys.exit(reader.getErrorMessage())

  # -----------------
  # Get a full model
  # -----------------

  # This will return a fully working model, with everything
  # stored in the file. Missing fields (e.g. tracks) will be
  # null or None in python
  model = reader.getModel()
  # model is a fiji.plugin.trackmate.Model

  # ----------------
  # Display results
  # ----------------

  # We can now plainly display the model. It will be shown on an
  # empty image with default magnification because we do not
  # specify an image to display it. We will see later how to
  # retrieve the image on which the data was generated.

  # A selection.
  sm = SelectionModel(model)

  # Read the display settings that was saved in the file.
  ds = reader.getDisplaySettings()

  # The viewer.
  # displayer = HyperStackDisplayer(model, sm, ds)
  # displayer.render()

  # ---------------------------------------------
  # Get only part of the data stored in the file
  # ---------------------------------------------

  # You might want to access only separate parts of the
  # model.

  # spots = model.getSpots()
  # spots is a fiji.plugin.trackmate.SpotCollection

  # logger.log(str(spots))

  # If you want to get the tracks, it is a bit trickier.
  # Internally, the tracks are stored as a huge mathematical
  # simple graph, which is what you retrieve from the file.
  # There are methods to rebuild the actual tracks, taking
  # into account for everything, but frankly, if you want to
  # do that it is simpler to go through the model:

  # trackIDs = model.getTrackModel().trackIDs(True)  # only filtered out ones
  # for id in trackIDs:
  #   logger.log(str(id) + ' - ' + str(model.getTrackModel().trackEdges(id)))

  # ---------------------------------------
  # Building a settings object from a file
  # ---------------------------------------

  # The settings object can be  used to regenerate the full
  # model from scratch. It contains all the parameters that
  # were configured when the GUI was used and the source image.

  # First load the image. TrackMate will retrieve the image *PATH*
  # from the XML file and try to reopen it. It works properly only
  # if the image was saved as an ImageJ TIFF.
  imp = reader.readImage()
  # Reading the image does not display it.

  # Now we can read the settings objects. Notice the method
  # asks for the image. This is the settings object will be built
  # with a link to the source image.
  settings = reader.readSettings(imp)

  # Let's print the settings object.
  logger.log(str('\n\nSETTINGS:'))
  logger.log(str(settings))

  # Display the image now.
  # imp.show()

  # With this, we can overlay the model and the source image:
  # displayer = HyperStackDisplayer(model, sm, imp, ds)
  # displayer.render()


  # here is the old code to run trackmate
  cal = imp.getCalibration()



  # ------------------------
  # Prepare settings object
  # ------------------------

  settings = Settings(imp)
  # setup = settings.toStringImageInfo()

  # Configure Cellpose default detector.

  settings.detectorFactory = CellposeDetectorFactory()

  settings.detectorSettings['TARGET_CHANNEL'] = 3
  settings.detectorSettings['OPTIONAL_CHANNEL_2'] = 0
  settings.detectorSettings['CELLPOSE_PYTHON_FILEPATH'] = "C:/Users/njhan/anaconda3/envs/Branches/python.exe"
  settings.detectorSettings['CELLPOSE_MODEL'] = PretrainedModel.CUSTOM
  settings.detectorSettings['CELLPOSE_MODEL_FILEPATH'] = "C:/Users/njhan/Box/macrophage_coculture/Cellpose_training/models/NH_LB_LC1_HL2"
  settings.detectorSettings['CELL_DIAMETER'] = 42.0
  settings.detectorSettings['USE_GPU'] = True
  settings.detectorSettings['SIMPLIFY_CONTOURS'] = False

  # Configure tracker
  settings.trackerFactory = SparseLAPTrackerFactory()
  settings.trackerSettings = settings.trackerFactory.getDefaultSettings()

  settings.trackerSettings['LINKING_MAX_DISTANCE'] = 100.0

  settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 50.0
  settings.trackerSettings['MAX_FRAME_GAP'] = 5

  settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
  settings.trackerSettings['SPLITTING_MAX_DISTANCE'] = 50.0

  settings.initialSpotFilterValue = -1.

  # Analyzers
  settings.addAllAnalyzers()

  # Add some filters for tracks/spots

  # filter on track duration = keep tracks > 75% of total duration
  # duration_threshold = 75
  # maxduration = (duration_threshold / 100.0) * (imp.getNFrames() * cal.frameInterval)
  # filter1_track = FeatureFilter('TRACK_DURATION', maxduration, True)
  # settings.addTrackFilter(filter1_track)

  # filter on spot = keep spots having radius > 1.6 um, and circularity > 0.7
  # filter1_spot = FeatureFilter('RADIUS', 1.6, True)
  # filter2_spot = FeatureFilter('CIRCULARITY', 0.7, True)
  # settings.addSpotFilter(filter1_spot)
  # settings.addSpotFilter(filter2_spot)

  print
  "Spot filters added = ", settings.getSpotFilters()
  print
  "Track filters added = ", settings.getTrackFilters(), "\n"

  # -------------------
  # Instantiate plugin
  # -------------------

  trackmate = TrackMate(model, settings)
  trackmate.computeSpotFeatures(True)
  trackmate.computeTrackFeatures(True)
  trackmate.getModel().setLogger(logger)

  # --------
  # Process
  # --------

  ok = trackmate.checkInput()
  if not ok:
    print(str(trackmate.getErrorMessage()))
    return

  # ok = trackmate.process()
  ok = trackmate.execTracking()
  if not ok:
    print(str(trackmate.getErrorMessage()))
    return

  # ----------------
  # Save results
  # ----------------

  saveFile = TMUtils.proposeTrackMateSaveFile(settings, logger)

  writer = TmXmlWriter(saveFile, logger)
  writer.appendLog(logger.toString())
  writer.appendModel(trackmate.getModel())
  writer.appendSettings(trackmate.getSettings())
  writer.writeToFile();
  print("Results saved to: " + saveFile.toString() + '\n');

  # lblImg = LabelImgExporter()
  # lblImg.createLabelImagePlus(trackmate, False, False, True, logger).show()
  # c = WindowManager.getCurrentImage()
  # IJ.save(c, "C:/Users/njhan/Box/macrophage_coculture/processed/label_images/" + c.getTitle())
  # c.close()

  # ----------------
  # Display results
  # ----------------

  if show_output:
    model = trackmate.getModel()
    selectionModel = SelectionModel(model)
    # ds = DisplaySettings()
    # ds = DisplaySettingsIO.readUserDefault()
    # ds.spotDisplayedAsRoi = True
    # displayer = HyperStackDisplayer(model, selectionModel, imp, ds)
    # displayer.render()
    # displayer.refresh()

  # capture overlay - RGB file
  # image = trackmate.getSettings().imp
  # capture = CaptureOverlayAction.capture(image, -1, imp.getNFrames(), logger)
  # capture.setTitle("TracksOverlay")
  # capture.show()


# ------------------------------------------------------

for file_path in file_paths:
  dt_string = dt.now().strftime("%d/%m/%Y %H:%M:%S")
  print('\nRunning analysis on %s - %s' % (str(file_path), dt_string))

  run(file_path)

print('Finished!')