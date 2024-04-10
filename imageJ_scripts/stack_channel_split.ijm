// https://gist.github.com/romainGuiet/cf42f3b1d31222a76d602dfe2f028894
dir = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/new_images/double_raw/";
output_dir = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/GM130/new_images/double/";
fileList = getFileList(dir);

//activate batch mode
setBatchMode(true);

// LOOP to process the list of files
for (i = 0; i < lengthOf(fileList); i++) {
	// define the "path" 
	// by concatenation of dir and the i element of the array fileList
	current_imagePath = dir+fileList[i];
	// check that the currentFile is not a directory
	if (!File.isDirectory(current_imagePath)){
		// open the image and split
		open(current_imagePath);
		// get some info about the image
		getDimensions(width, height, channels, slices, frames);
		// if it's a multi channel image , or a RGB
		if ((channels > 1)||( bitDepth() == 24 )) run("Split Channels");

		// now we save all the generated images as tif in the output_dir
		ch_nbr = nImages ; 
		suffix_array = newArray("_nuc", "_gm130");
		for ( c = 1 ; c <= ch_nbr ; c++){
			selectImage(c);
			currentImage_name = getTitle();
			print(currentImage_name);
			currentImage_name = substring(currentImage_name,3,lengthOf(currentImage_name));
			currentImage_name = substring(currentImage_name,0,lengthOf(currentImage_name)-4);
			print(currentImage_name);
			saveAs("tiff", output_dir + currentImage_name + suffix_array [c-1]); // -1 because array start at 0
		}
		// make sure to close every images befores opening the next one
		run("Close All");
	}
}
setBatchMode(false);
