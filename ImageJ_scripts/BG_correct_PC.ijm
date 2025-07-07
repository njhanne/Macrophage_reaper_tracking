// https://gist.github.com/romainGuiet/cf42f3b1d31222a76d602dfe2f028894
//dir = "D:/UCSF/macrophage_video_analysis/processed/BG_corrected/test2/";
//output_dir = "D:/UCSF/macrophage_video_analysis/processed/BG_corrected/test2/";

dir = "E:/Nicholas/processed/8bit_tiffs/run/";
output_dir = "E:/Nicholas/processed/BG_corrected/";

fileList = getFileList(dir);
print(lengthOf(fileList));

//activate batch mode
setBatchMode(true);

// LOOP to process the list of files
for (i = 0; i < lengthOf(fileList); i++) {
	// define the "path" 
	// by concatenation of dir and the i element of the array fileList
	current_imagePath = dir+fileList[i];
	// check that the currentFile is not a directory
	if (!File.isDirectory(current_imagePath)) {
		if (endsWith(current_imagePath, ".tif")) {
            // open the image and split
            run("Bio-Formats Macro Extensions");
            Ext.openImagePlus(current_imagePath);
            //open(current_imagePath);
            // get some info about the image
            img_name = getTitle();
            print(img_name);
            currentImage_name = substring(img_name,0,lengthOf(img_name)-4);


            run("Split Channels");

			selectImage("C3-"+img_name);
            run("Subtract Background...", "rolling=15 sliding stack");
            
            // save near end of timestack image
            run("Make Subset...", "slices=355");
            selectImage("Substack (355)");
			run("Enhance Contrast", "saturated=0.35");
			run("RGB Color");
			saveAs("tiff", currentImage_name+"_endNucView");
            

            selectImage("C1-"+img_name);
            run("Enhance Contrast...", "saturated=1 normalize process_all");

            command_str = "c1=C2-" + img_name + " c2=C3-" + img_name + " c3=C1-" + img_name + " create";
            run("Merge Channels...", command_str);


            currentImage_name = currentImage_name+"_BG";
            print(currentImage_name);
            saveAs("tiff", currentImage_name);
        }
        // make sure to close every images before opening the next one
        run("Close All");
	}
}
print('Done!');
setBatchMode(false);
