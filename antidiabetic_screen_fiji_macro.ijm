// ============================================================
// PLATE A — FULL CZI SCREEN ANALYSIS
//
// One CZI file = one field of view
// Example filename: b2(2)_40x.czi
//
// Expected structure:
// - 4 channels
// - 5 Z-slices
//
// Channel order:
// 1 = pRab10
// 2 = MAP2
// 3 = Total Rab10
// 4 = DAPI
//
// Output:
// One row per image containing:
// Well, Field, Area, fluorescence IntDen values and
// pRab10 / Total Rab10 ratio
// ============================================================


// Choose folder containing all Plate A CZI files
inputDir = getDirectory("Choose the folder containing Plate A CZI files");

files = getFileList(inputDir);


// Create storage arrays
wellValues = newArray(files.length);
fieldValues = newArray(files.length);
areaValues = newArray(files.length);

pRab10Values = newArray(files.length);
totalRab10Values = newArray(files.length);
ratioValues = newArray(files.length);

map2Values = newArray(files.length);
dapiValues = newArray(files.length);


// Measure Area and Integrated Density
run(
    "Set Measurements...",
    "area integrated redirect=None decimal=3"
);


// Clear any old Results table
if (isOpen("Results")) {
    selectWindow("Results");
    run("Clear Results");
}


processedFiles = 0;
skippedFiles = 0;


// Batch mode ON for the full plate
setBatchMode(true);


for (i = 0; i < files.length; i++) {

    fileName = files[i];
    lowerName = toLowerCase(fileName);

    // Process only CZI files
    if (endsWith(lowerName, ".czi")) {

        fullPath = inputDir + fileName;

        // Remove .czi extension
        baseName = substring(
            fileName,
            0,
            lengthOf(fileName) - 4
        );


        // Extract well and field from b2(2)_40x
        openBracket = indexOf(baseName, "(");
        closeBracket = indexOf(baseName, ")");

        if (openBracket >= 0 && closeBracket > openBracket) {

            wellName = substring(
                baseName,
                0,
                openBracket
            );

            fieldNumber = substring(
                baseName,
                openBracket + 1,
                closeBracket
            );

        } else {

            wellName = baseName;
            fieldNumber = "NA";
        }


        print("Processing: " + fileName);


        // Open the CZI file
        open(fullPath);

        originalID = getImageID();


        // Read image dimensions
        getDimensions(
            imgWidth,
            imgHeight,
            channelCount,
            zCount,
            frameCount
        );


        // Check expected image structure
        if (channelCount != 4 || zCount != 5) {

            print(
                "SKIPPED: " + fileName +
                " has " + channelCount +
                " channels and " + zCount +
                " Z-slices."
            );

            close();
            skippedFiles++;
            continue;
        }


        // Create maximum-intensity projection
        run("Z Project...", "projection=[Max Intensity]");

        projectionID = getImageID();


        // Close original Z-stack
        selectImage(originalID);
        close();


        // Select maximum projection
        selectImage(projectionID);


        // Clear temporary measurements
        if (isOpen("Results")) {
            selectWindow("Results");
            run("Clear Results");
            selectImage(projectionID);
        }


        // Channel 1: pRab10
        Stack.setChannel(1);
        run("Select All");
        run("Measure");

        measuredArea = getResult("Area", 0);
        pRab10IntDen = getResult("IntDen", 0);


        // Channel 2: MAP2
        Stack.setChannel(2);
        run("Select All");
        run("Measure");

        map2IntDen = getResult("IntDen", 1);


        // Channel 3: Total Rab10
        Stack.setChannel(3);
        run("Select All");
        run("Measure");

        totalRab10IntDen = getResult("IntDen", 2);


        // Channel 4: DAPI
        Stack.setChannel(4);
        run("Select All");
        run("Measure");

        dapiIntDen = getResult("IntDen", 3);


        // Calculate pRab10 / Total Rab10 ratio
        if (totalRab10IntDen > 0) {
            pRab10Ratio = pRab10IntDen / totalRab10IntDen;
        } else {
            pRab10Ratio = NaN;
        }


        // Store image-level results
        wellValues[processedFiles] = wellName;
        fieldValues[processedFiles] = fieldNumber;
        areaValues[processedFiles] = measuredArea;

        pRab10Values[processedFiles] = pRab10IntDen;
        totalRab10Values[processedFiles] = totalRab10IntDen;
        ratioValues[processedFiles] = pRab10Ratio;

        map2Values[processedFiles] = map2IntDen;
        dapiValues[processedFiles] = dapiIntDen;


        processedFiles++;


        // Close maximum projection
        selectImage(projectionID);
        close();
    }
}


setBatchMode(false);


// Remove temporary channel-level measurements
if (isOpen("Results")) {
    selectWindow("Results");
    run("Clear Results");
}


// Build final clean table
for (r = 0; r < processedFiles; r++) {

    setResult("Well", r, wellValues[r]);
    setResult("Field", r, fieldValues[r]);
    setResult("Area", r, areaValues[r]);

    setResult(
        "pRab10_IntDen",
        r,
        pRab10Values[r]
    );

    setResult(
        "Total_Rab10_IntDen",
        r,
        totalRab10Values[r]
    );

    setResult(
        "pRab10_Total_Ratio",
        r,
        ratioValues[r]
    );

    setResult(
        "MAP2_IntDen",
        r,
        map2Values[r]
    );

    setResult(
        "DAPI_IntDen",
        r,
        dapiValues[r]
    );
}


updateResults();


if (processedFiles == 0) {

    showMessage(
        "No results",
        "No CZI files were processed.\n\nCheck the Log window."
    );

} else {

    // Save automatically inside the selected folder
    outputPath = inputDir + "Plate_A_image_level_measurements.csv";

    saveAs("Results", outputPath);

    showMessage(
        "Plate A finished",
        "Files processed: " + processedFiles +
        "\nFiles skipped: " + skippedFiles +
        "\n\nResults saved as:\n" +
        outputPath
    );
}
