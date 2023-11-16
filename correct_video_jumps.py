import cv2
from tifffile import imwrite
import numpy as np

import sys
sys.path.append("./DirFileHelpers")

### Main ###
# We have time-series fluorescent stacks taken of in vitro culture of MSC and macrophages
# I think sometimes someone comes in and switches media on the cells which makes the video jump
# These jumps prevent celltracker from keeping track of cell movement as they move too far for it to register
# This code should create a larger image that moves the image to prevent jumps
# Unfortunately this requires a human to manually note when the jump occurs and how far they move
# We will read in a csv of the jump info and use it to create a new image stack

