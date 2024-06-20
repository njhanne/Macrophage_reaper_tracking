// https://gist.github.com/romainGuiet/cf42f3b1d31222a76d602dfe2f028894
dir = "D:/UCSF/macrophage_video_analysis/processed/8bit-tiffs/";
fileList = getFileList(dir);

//activate batch mode
setBatchMode(true);

// LOOP to process the list of files
for (i = 0; i < lengthOf(fileList); i++) {
	// define the "path" 
	// by concatenation of dir and the i element of the array fileList
	current_imagePath = dir+fileList[i];
	// check that the currentFile is not a directory
	if (!File.isDirectory(current_imagePath)) {
		if (endsWith(current_imagePath, ".tiff")) {
            // open the image and split
            run("Bio-Formats Macro Extensions");
            Ext.openImagePlus(current_imagePath);
            //open(current_imagePath);
            // get some info about the image
            img_name = getTitle();
            print(img_name);

            setOption("ScaleConversions", true);
            run("8-bit");

            saveAs("tiff", current_imagePath);
        }
        // make sure to close every images before opening the next one
        run("Close All");
	}
}
print('Done!');
setBatchMode(false);
