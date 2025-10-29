macro "Cryo-SEM Pore Analysis (pixel size prompt)" {

    // User inputs
    // Folder with images you want to analyse
    dir = "C:/Users/walsh/Downloads/SAMJ_Analysis_Results (1)/Cryo/Cryo-SEM 10000x/";
    // Put one or more .tif filenames here
    list = newArray("CRYO-SEM pp x10000a.tif");



    // Ask once for pixel size. This is just written to metadata; it does NOT change area filtering.
    pxSize = getNumber("Pixel size (nm/pixel). You can leave 0 if unknown:", 0.0);

    // Main loop over each image in list 
    for (f = 0; f < list.length; f++) {

        // timing (just to keep track of how long each image takes)
        startTime = getTime();

        // open image
        openPath = dir + list[f];
        open(openPath);
        originalname = getTitle();
        saveDir = dir;

        // base name without extension
        name = originalname;
        name = replace(name, ".tif", "");
        name = replace(name, ".tiff", "");

        // record run metadata so we can reconstruct later
        runMeta =
            "Image: " + openPath + "\n" +
            "Start(ms): " + startTime + "\n" +
            "PixelSize_nm_per_px: " + pxSize + "\n";

        outputDir = saveDir;
        File.saveString(runMeta, outputDir + name + "_run_metadata.txt");

        // keep a measuring copy for final measurements
        run("Duplicate...", "title=Measuring_Copy");
        MeasureCopy = getTitle();
        selectWindow(originalname);

        
        // Seed generation (pre-SAMJ)
        
        Dialog.create("Seed settings (values with * are important)");
        Dialog.addCheckbox("Generate seeds from threshold", true);
        Dialog.addNumber("(*) Threshold min (0-255):", 200);
        Dialog.addNumber("(*) Threshold max (0-255):", 255);
        Dialog.addNumber("Enhance Contrast (%% saturated):", 20);
        Dialog.addNumber("Ignore ROIs touching edge (<px):", 10);
        Dialog.show();

        useSeeds     = Dialog.getCheckbox();
        thrMin       = Dialog.getNumber();
        thrMax       = Dialog.getNumber();
        satPct       = Dialog.getNumber();
        edgeMarginPx = Dialog.getNumber();

        // Clean ROI Manager before we start
        if (roiManager("count") > 0) roiManager("Reset");

        // Build SAMJ prompt points by thresholding / particles
        if (useSeeds) {

            // Work on a temporary duplicate so we don't mangle the raw micrograph
            selectWindow(originalname);
            if (isOpen("Workable_Copy")) { selectWindow("Workable_Copy"); run("Close"); }
            run("Duplicate...", "title=Workable_Copy");
            selectWindow("Workable_Copy");
            run("Select None");

            // prep: 8-bit, invert contrast, boost
            run("8-bit");
            run("Invert");
            run("Enhance Contrast...", "saturated=" + satPct);

            // hard threshold -> binary mask
            setThreshold(thrMin, thrMax);
            run("Convert to Mask");

            // First pass particle detection (used here partly to knock out the dominant background blob)
            roiManager("Reset");
            run("Set Measurements...", "area redirect=None decimal=3");
            run("Analyze Particles...", "size=1-Infinity add show=Nothing");

            nC = roiManager("count");
            if (nC > 0) {
                maxArea = -1;
                maxIdx  = -1;
                for (i = 0; i < nC; i++) {
                    roiManager("Select", i);
                    getStatistics(area, mean, min, max, std);
                    if (area > maxArea) { maxArea = area; maxIdx = i; }
                }
                if (maxIdx >= 0) {
                    // kill the biggest ROI (usually background sheet / ice / etc.)
                    // Nota: “Apagar el ROI más grande pintándolo de negro (fondo).”
                    roiManager("Select", maxIdx);
                    run("Restore Selection");
                    run("Invert"); // fill that selection with black
                    run("Select None");

                    // also invert on the original, to be safe
                    selectWindow(originalname);
                    roiManager("Select", maxIdx);
                    run("Invert");
                    run("Select None");
                }
            }

            // Re-run particles to get actual pore-ish ROIs
            if (roiManager("count") > 0) roiManager("Reset");
            roiManager("Show All");
            selectWindow("Workable_Copy");
            run("Set Measurements...", "area redirect=None decimal=3");
            run("Analyze Particles...", "size=1-Infinity pixel clear add");

            // toss ROIs that touch the border (partial pores at edge)
            imgW = getWidth();
            imgH = getHeight();
            nC = roiManager("count");
            toDel = newArray(nC); k = 0;
            for (i = 0; i < nC; i++) {
                roiManager("Select", i);
                getSelectionBounds(x, y, w, h);
                del = 0;
                if (x < edgeMarginPx) del = 1;
                if (y < edgeMarginPx) del = 1;
                if ((x + w) > (imgW - edgeMarginPx)) del = 1;
                if ((y + h) > (imgH - edgeMarginPx)) del = 1;
                if (del == 1) { toDel[k] = i; k = k + 1; }
            }
            for (j = k - 1; j >= 0; j--) {
                idx = toDel[j];
                if (idx >= 0) { roiManager("Select", idx); roiManager("Delete"); }
            }

            // convert each ROI to its bounding-box centre point.
            nC = roiManager("count");
            if (nC < 1) exit("No ROIs to convert into prompts.");

            xarr = newArray(nC);
            yarr = newArray(nC);

            for (i = 0; i < nC; i++) {
                roiManager("Select", i);
                getSelectionBounds(x, y, w, h);
                xarr[i] = x + w/2.0;
                yarr[i] = y + h/2.0;

                // clamp to image
                if (xarr[i] < 0) xarr[i] = 0;
                if (yarr[i] < 0) yarr[i] = 0;
                if (xarr[i] > imgW-1) xarr[i] = imgW-1;
                if (yarr[i] > imgH-1) yarr[i] = imgH-1;
            }

            // Put those points back on the ORIGINAL image. SAMJ will read them.
            selectWindow(originalname);
            run("Invert"); // just to make pores bright for SAMJ prompting
            makeSelection("point", xarr, yarr);

            // Replace ROI Manager contents with this one multi-point selection
            roiManager("Reset");
            roiManager("Add");
            roiManager("Show All");

            // Supplemental seeds (peak finding away from existing seeds) 
            // Rationale: sometimes thresholding under-samples distinct pores,
            // so we add extra points at strong local maxima that are not too close
            // to what we already marked.

            supp_satPct        = 20;   // contrast boost before maxima finding
            supp_sigma         = 1.0;  // small blur to calm noise
            supp_prominence    = 28;   // higher = pick only strong maxima
            seedMinDistPx      = 9;    // don't add a seed if it's basically on top of an existing one
            seedEdgeMarginPx   = 10;   // avoid extrema right on the image edge

            // grab the current multi-point (these are our "existing seeds")
            selectWindow(originalname);
            roiManager("Select", 0);
            run("Restore Selection");
            getSelectionCoordinates(px, py);
            nSeeds = px.length;

            // helper image for maxima search
            if (isOpen("SeedHelper")) { selectWindow("SeedHelper"); run("Close"); }
            run("Duplicate...", "title=SeedHelper");
            selectWindow("SeedHelper");
            run("8-bit");
            run("Enhance Contrast...", "saturated=" + supp_satPct);
            if (supp_sigma > 0) run("Gaussian Blur...", "sigma=" + supp_sigma);

            // find maxima, but ignore a safety margin at the border
            margin = seedEdgeMarginPx;
            makeRectangle(margin, margin, getWidth()-2*margin, getHeight()-2*margin);
            run("Find Maxima...", "prominence=" + supp_prominence + " output=[Point Selection]");
            getSelectionCoordinates(cx, cy);
            nCand = cx.length;

            // keep only maxima that are not too close to existing seeds
            keepFlag = newArray(nCand);
            m = 0;
            for (i = 0; i < nCand; i++) {
                keep = 1;
                for (j = 0; j < nSeeds; j++) {
                    dx = cx[i] - px[j];
                    dy = cy[i] - py[j];
                    if (dx*dx + dy*dy < seedMinDistPx*seedMinDistPx) { keep = 0; break; }
                }
                if (keep == 1) {
                    keepFlag[i] = 1;
                    m = m + 1;
                } else {
                    keepFlag[i] = 0;
                }
            }

            // collect the "kept" supplemental seeds
            newX = newArray(m);
            newY = newArray(m);
            idx = 0;
            for (i = 0; i < nCand; i++) {
                if (keepFlag[i] == 1) {
                    newX[idx] = cx[i];
                    newY[idx] = cy[i];
                    idx = idx + 1;
                }
            }

            // merge original seeds + supplemental seeds
            total = nSeeds + m;
            combinedX = newArray(total);
            combinedY = newArray(total);

            for (i = 0; i < nSeeds; i++) {
                combinedX[i] = px[i];
                combinedY[i] = py[i];
            }
            for (i = 0; i < m; i++) {
                combinedX[nSeeds + i] = newX[i];
                combinedY[nSeeds + i] = newY[i];
            }

            // overwrite ROI Manager with merged seeds
            selectWindow(originalname);
            makeSelection("point", combinedX, combinedY);

            roiManager("Reset");
            roiManager("Add");
            roiManager("Show All");

            // cleanup temp helper
            if (isOpen("SeedHelper")) { selectWindow("SeedHelper"); run("Close"); }
            if (isOpen("Workable_Copy")) { selectWindow("Workable_Copy"); run("Close"); }
        }

       
        // SAMJ step (manual)
        // NOTE FOR USER:
        // 1. Plugins > SAMJ > SAMJ Annotator
        // 2. Click "Go!" to encode this image.
        // 3. Use "Preset prompts: ROI Manager" to import the points we just made.
        // 4. Models: I usually use "SAM2 Small" (GPU) or "EfficientViT-SAM-L2" (CPU).
        // 5. Turn OFF "only return largest ROI".
        // 6. Send the resulting ROIs back to the standard ROI Manager.
        waitForUser("SAMJ step",
            "Run SAMJ now (see notes in code). When ROIs are in ROI Manager, click OK.");

        // sanity check
        if (roiManager("count") < 1) exit("No ROIs came back from SAMJ.");

        
        // Filter ROIs by area and edge BEFORE merging them
        
        Dialog.create("Filter ROIs before merge");
        Dialog.addNumber("Min area (px^2):", 0.01);
        Dialog.addNumber("Max area (px^2):", 800);
        Dialog.addNumber("Ignore ROIs touching edge (<px):", 10);
        Dialog.show();

        AminPx        = Dialog.getNumber();
        AmaxPx        = Dialog.getNumber();
        edgeMarginPost = Dialog.getNumber();

        imgW = getWidth();
        imgH = getHeight();

        n = roiManager("count");
        markDelete = newArray(n);

        for (i = 0; i < n; i++) {
            roiManager("Select", i);
            getStatistics(area, mean, min, max, std);
            areaPx = area;

            getSelectionBounds(x, y, w, h);

            touchesEdge =
                (x < edgeMarginPost) ||
                (y < edgeMarginPost) ||
                ((x + w) > (imgW - edgeMarginPost)) ||
                ((y + h) > (imgH - edgeMarginPost));

            del = 0;
            if (areaPx < AminPx) del = 1;
            if (areaPx > AmaxPx) del = 1;
            if (touchesEdge)     del = 1;

            markDelete[i] = del;
        }

        // delete from bottom up so indices don't shift
        for (j = n - 1; j >= 0; j--) {
            if (markDelete[j] == 1) {
                roiManager("Select", j);
                roiManager("Delete");
            }
        }

        if (roiManager("count") < 1) exit("Everything got filtered out pre-merge.");

        
        // Merge overlapping ROIs


        // Paint all current ROIs into a fresh binary mask
        selectWindow(originalname);
        run("Select None");
        getDimensions(w, h, c, z, t);

        tempTitle = "Roi_Mask_Temp";
        if (isOpen(tempTitle)) { selectWindow(tempTitle); run("Close"); }
        run("New...", "name=" + tempTitle + " type=8-bit fill=Black width=" + w + " height=" + h + " slices=1");
        selectWindow(tempTitle);
        run("8-bit");
        run("Select None");

        // fill ROIs onto mask as white blobs
        selectWindow(originalname);
        roiManager("Deselect");
        selectWindow(tempTitle);
        roiManager("Fill");

        // binarise to be sure it's 0/255
        run("Convert to Mask");

        // make sure SAM Roi Manager isn't shadowing anything
        if (isOpen("SAM Roi Manager")) { selectWindow("SAM Roi Manager"); run("Close"); }
        run("ROI Manager...");

        // Re-extract connected components from that mask = merged ROIs
        selectWindow(tempTitle);

        // Manual cleanup step:
        // You'll get tiled view, draw rectangles over junk regions and press Ctrl+F
        // (that's "Fill" to black them out). Then click OK.
        run("Tile");
        roiManager("Show All without labels");
        roiManager("Show None");
        setTool("rectangle");
        waitForUser(
            "Manual cleanup",
            "Remove obvious noise:\n" +
            "1. Draw a box around junk.\n" +
            "2. Press Ctrl+F to fill black.\n" +
            "Repeat as needed, then click OK."
        );

        run("Select None");
        run("Set Measurements...", "area redirect=None decimal=3");
        run("Analyze Particles...", "size=1-Infinity add clear show=Nothing redirect=None");

        // close temp mask window
        close(tempTitle);

        // Filter ROIs again after merge
        imgW = getWidth();
        imgH = getHeight();

        n2 = roiManager("count");
        if (n2 < 1) exit("No ROIs after merge.");

        markDelete = newArray(n2);
        for (i = 0; i < n2; i++) {
            roiManager("Select", i);
            getStatistics(area, mean, min, max, std);
            areaPx = area;

            getSelectionBounds(x, y, w, h);

            touchesEdge =
                (x < edgeMarginPost) ||
                (y < edgeMarginPost) ||
                ((x + w) > (imgW - edgeMarginPost)) ||
                ((y + h) > (imgH - edgeMarginPost));

            del = 0;
            if (areaPx < AminPx) del = 1;
            if (areaPx > AmaxPx) del = 1;
            if (touchesEdge)     del = 1;

            markDelete[i] = del;
        }

        for (j = n2 - 1; j >= 0; j--) {
            if (markDelete[j] == 1) {
                roiManager("Select", j);
                roiManager("Delete");
            }
        }

        if (roiManager("count") < 1) exit("Nothing survived post-merge filter.");

        
        // Final measurements (Feret, perimeter, etc.)
        selectWindow(MeasureCopy);
        roiManager("Show All");
        run("Set Measurements...", "area centroid perimeter feret's shape redirect=None decimal=3");

        nFinal = roiManager("count");
        for (i = 0; i < nFinal; i++) {
            roiManager("Select", i);
            run("Measure");
        }

        
        // Save outputs
        

        // 1) Save the Results table as CSV alongside the input image
        if (nResults() > 0) {
            csvPathLocal = saveDir + name + "_Pore_Measurements.csv";
            saveAs("Results", csvPathLocal);
        } else {
            print("Warning: Results table is empty, nothing saved.");
        }

        // Save final pore mask that matches what we just measured
        if (roiManager("count") > 0) {

            // combine the final ROIs -> mask
            selectWindow(originalname);
            roiManager("Deselect");
            roiManager("Combine");
            run("Create Mask"); // opens new mask window
            rename("Final_Pore_Mask");

            selectWindow("Final_Pore_Mask");
            run("8-bit");
            run("Convert to Mask");

            maskPathLocal = saveDir + name + "_Final_Pore_Mask.tif";
            saveAs("Tiff", maskPathLocal);
            print("Saved CSV + mask to " + saveDir);

            // also store mask copy in the GitHub folder
            saveAs("Tiff", outputDir + name + "_Final_Pore_Mask.tif");

            // close mask windows
            close(name + "_Final_Pore_Mask");
            close("Final_Pore_Mask");

        } else {
            print("Note: no ROIs in manager at save step, so no final mask was written.");
        }

        // Timing log
        endTime   = getTime();
        elapsedSec = (endTime - startTime) / 1000.0;
        elapsedStr = d2s(elapsedSec, 3);

        print("Done with", name, "| ROIs:", roiManager("count"), "| ~", elapsedStr, "s");

        txtPathLocal = saveDir + name + "_elapsed_seconds.txt";
        logLine = "analysis_time_s=" + elapsedStr;
        File.saveString(logLine, txtPathLocal);
        File.saveString(logLine, outputDir + name + "_elapsed_seconds.txt");


        // append timing info to the run metadata file
        File.append("\nAnalysisTime_s=" + elapsedStr + "\n", outputDir + name + "_run_metadata.txt");

        // clean up windows before next image
        if (isOpen(MeasureCopy))   { selectWindow(MeasureCopy); run("Close"); }
        if (isOpen(originalname))  { selectWindow(originalname); run("Close"); }
    }
}
