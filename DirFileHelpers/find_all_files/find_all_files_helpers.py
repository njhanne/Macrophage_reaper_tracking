import os
import re

def find_all_filepaths(directory, extension, file_string=''):
  # this iterates through the directory and subdirectories and gets all the filepaths that end in your extension
  # '.extension' is how they should be passed
  # https://stackoverflow.com/a/59803793
  subfolders, filepaths = [], []

  for f in os.scandir(directory):
    if f.is_dir():
      subfolders.append(f.path)
    if f.is_file():
      if os.path.splitext(f.name)[1].lower() in extension:
        if re.search(file_string, os.path.splitext(f.name)[0]):
          filepaths.append(f.path)

  for dir in list(subfolders):
    sf, f = find_all_filepaths(dir, extension, file_string)
    subfolders.extend(sf)
    filepaths.extend(f)
  return subfolders, filepaths