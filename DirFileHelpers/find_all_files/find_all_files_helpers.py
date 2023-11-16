from pathlib import Path
import re

def find_all_filepaths(directory, extension, file_string=''):
  # this iterates through the directory and subdirectories and gets all the filepaths that end in your extension
  # '.extension' is how they should be passed
  # https://stackoverflow.com/a/59803793
  subfolders, filepaths = [], []

  # for f in os.scandir(directory):
  for f in directory.glob('*'):
    if f.is_dir():
      subfolders.append(f)
    if f.is_file():
      if f.suffix in extension:
        if re.search(file_string, f.stem):
          filepaths.append(f)

  for dir in list(subfolders):
    sf, f = find_all_filepaths(dir, extension, file_string)
    subfolders.extend(sf)
    filepaths.extend(f)
  return subfolders, filepaths