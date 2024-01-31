# https://github.com/trackmate-sc/TrackMate/blob/master/scripts/ExampleScript_1.py
# https://github.com/trackmate-sc/TrackMate/blob/master/scripts/TrackMate-Cellpose-script.py

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
from fiji.plugin.trackmate.util import TMUtils
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.action import LabelImgExporter

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
directory_to_process = os.path.normpath("C:/Users/njhan/Box/Collaboration/macrophage_video")
file_paths = find_all_filepaths(directory_to_process, '.tif')


# ------------------------------------------------------
# 	ACTUAL CODE.
# ------------------------------------------------------
def run(image_file):
  # Open image.
  imp = IJ.openImage(image_file)
  cal = imp.getCalibration()

  # Logger -> content will be saved in the XML file.
  logger = LogRecorder(Logger.VOID_LOGGER)

  logger.log('TrackMate-Cellpose analysis script\n')

  dt_string = dt.now().strftime("%d/%m/%Y %H:%M:%S")
  logger.log(dt_string + '\n\n')

  # ------------------------
  # Prepare settings object
  # ------------------------

  settings = Settings(imp)
  setup = settings.toStringImageInfo()

  # Configure Cellpose default detector.

  settings.detectorFactory = CellposeDetectorFactory()

  settings.detectorSettings['TARGET_CHANNEL'] = 3
  settings.detectorSettings['OPTIONAL_CHANNEL_2'] = 0
  settings.detectorSettings['CELLPOSE_PYTHON_FILEPATH'] = "C:/Users/njhan/anaconda3/envs/Branches/python.exe"
  settings.detectorSettings['CELLPOSE_MODEL'] = PretrainedModel.CUSTOM
  settings.detectorSettings['CELLPOSE_MODEL_FILEPATH'] = "C:\Users\njhan\Box\Collaboration\macrophage_video\NH_LB_LC1_HL1"
  settings.detectorSettings['CELL_DIAMETER'] = 30.0
  settings.detectorSettings['USE_GPU'] = True
  settings.detectorSettings['SIMPLIFY_CONTOURS'] = False

  # Configure tracker
  settings.trackerFactory = SparseLAPTrackerFactory()
  settings.trackerSettings = settings.trackerFactory.getDefaultSettings()

  settings.trackerSettings['LINKING_MAX_DISTANCE'] = 70.0
  settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 10.0
  settings.trackerSettings['MAX_FRAME_GAP'] = 3
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

  trackmate = TrackMate(settings)
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

  ok = trackmate.process()
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

  lblImg = LabelImgExporter()
  lblImg.createLabelImagePlus(trackmate, False, False, True, logger).show()
  c = WindowManager.getCurrentImage()
  IJ.save(c, "C:/Users/njhan/Box/Collaboration/macrophage_video/" + c.getTitle())
  c.close()

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