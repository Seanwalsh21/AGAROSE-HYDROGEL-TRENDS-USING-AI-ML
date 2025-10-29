macro "Cryo-SEM Pore Analysis (Prompt for pixels/µm)" {
    // Ask user to select folder containing images
    dir = getDirectory("Choose a folder with Cryo-SEM images");
    list = getFileList(dir);

    // Set up output folder for results
    outputDir = dir + "Results/";
    File.makeDirectory(outputDir);

    // Start fresh summary file
    summaryFile = outputDir + "summary.csv";
    File.delete(summaryFile); // just in case it already exists
    File.append("Filename,Object,Area,CentroidX,CentroidY,Perimeter,Circularity,AR,Roundness,Feret,Major,Minor,Angle\n", summaryFile);

    setBatchMode(true); // run faster, don’t update GUI constantly

    // Ask for scale — user enters pixels per micron
    pixelsPerMicron = getNumber("Enter pixels per µm (e.g. 107):", 107.0);
    pixelSize = 1 / pixelsPerMicron;

    // Loop through each image in the folder
    for (i = 0; i < list.length; i++) {
        filename = list[i];

        // Accept only image files
        if (endsWith(filename, ".tif") || endsWith(filename, ".tiff") || endsWith(filename, ".jpg") || endsWith(filename, ".png")) {
            
            open(dir + filename);
            run("Set Scale...", "distance=1 known=" + pixelSize + " pixel=1 unit=µm global");

            // Preprocessing
            run("8-bit"); // ensure consistent bit depth
            run("Enhance Contrast", "saturated=0.35"); // stretch histogram a bit

            // Threshold and segment
            run("Auto Local Threshold", "method=Phansalkar radius=15 parameter_1=0 parameter_2=0 white");
            run("Invert"); // make pores white on black
            setOption("BlackBackground", false);
            run("Make Binary");
            run("Watershed"); // separates connected pores

            // Duplicate original for annotated overlay
            selectWindow(filename);
            run("Duplicate...", "title=Overlay_copy");

            // Run measurements on overlay copy
            selectWindow("Overlay_copy");
            run("Set Measurements...", "area centroid perimeter shape feret's fit redirect=None decimal=3");
            run("Analyze Particles...", "size=0.5-Infinity show=Overlay display clear include add");

            // Save overlay image and per-image results
            saveAs("PNG", outputDir + replace(filename, ".tif", "") + "_annotated_overlay.png");
            saveAs("Results", outputDir + replace(filename, ".tif", "") + "_results.csv");

            // Append each pore’s data to master summary CSV
            selectWindow("Results");
            count = nResults;
            for (j = 0; j < count; j++) {
                row = filename + "," + (j+1);
                row += "," + getResult("Area", j);
                row += "," + getResult("X", j);
                row += "," + getResult("Y", j);
                row += "," + getResult("Perim.", j);
                row += "," + getResult("Circ.", j);
                row += "," + getResult("AR", j);
                row += "," + getResult("Round", j);
                row += "," + getResult("Feret", j);
                row += "," + getResult("Major", j);
                row += "," + getResult("Minor", j);
                row += "," + getResult("Angle", j);
                File.append(row + "\n", summaryFile);
            }

            close("Results");
            close(); // close image
        }
    }

    setBatchMode(false);
    print("Results: " + outputDir);
}
