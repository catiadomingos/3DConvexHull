// Macro to calculate 3D Convex Hull in Fiji (ImageJ)

// Retrieve the image path from the arguments
imagePath = getArgument();

// Open the binarized image
open(imagePath);

// Get the sample name from the image title
sampleName = getTitle();

// Measure all the 3D Convex Hull statistics
run("Measure All...");

// Check if the results window has data
if (nResults() == 0) {
    print("Error: No measurements found. Ensure that 'Measure All' command is working.");
    exit();
}

// Specify the output file path
outputFilePath = "D:/PhD_Catia Domingos/3. Programs/Matlab/2. Analysis/2. ExM analysis/Analysis Convex Hull/ConvexHull_Results_Raw.csv";
print("Output File Path: " + outputFilePath);

// Load existing results if the file exists
if (File.exists(outputFilePath)) {
    existingResults = File.openAsString(outputFilePath);
    print("Existing results loaded.");
} else {
    existingResults = "";
    print("No existing results found.");
}

// Save the new results to a temporary file
tempFilePath = "temp_results.csv";
saveAs("Results", tempFilePath);
print("Temporary results saved to: " + tempFilePath);

// Read the new results
newResults = File.openAsString(tempFilePath);
print("New results read from temp file. Length: " + newResults.length());

// Output the first few lines of new results for debugging
print("New results snippet:\n" + newResults);

// Prepare results for appending
newResultsLines = split(newResults, "\n");
headers = "";
newResultsNoHeaders = "";
dataStart = false;

// Check and process headers
if (newResultsLines.length > 1) {
    if (newResultsLines[0].contains("Object Voxels")) {
        headers = newResultsLines[0]; // Capture headers
        headers += ",Concave Volume Ratio"; // Append new column header
        dataStart = true;
        for (i = 1; i < newResultsLines.length; i++) {
            newResultsNoHeaders += newResultsLines[i] + "\n";
        }
    } else {
        newResultsNoHeaders = newResults;
    }
} else {
    newResultsNoHeaders = newResults;
}

// Prepare final results
newResultsLines = split(newResultsNoHeaders, "\n");
updatedResults = "";
for (i = 0; i < newResultsLines.length; i++) {
    if (newResultsLines[i] != "") {
        columns = split(newResultsLines[i], ","); // Use comma for CSV
        if (columns.length > 1) {
            // Calculate Concave Volume Ratio
            convexVolume = parseFloat(columns[4]);
            objectVolume = parseFloat(columns[3]);
            if (objectVolume != 0) {
                concaveVolumeRatio = (convexVolume - objectVolume) / objectVolume;
            } else {
                concaveVolumeRatio = 0; // Handle division by zero
            }
            // Append calculated value to the row
            columns[0] = sampleName;
            updatedResults += columns[0];
            for (j = 1; j < columns.length; j++) {
                updatedResults += "," + columns[j]; // Use comma for CSV
            }
            updatedResults += "," + concaveVolumeRatio; // Append the new column value
            updatedResults += "\n";
        }
    }
}

// Combine existing and new results without repeating headers
if (existingResults == "") {
    combinedResults = headers + "\n" + updatedResults;
} else {
    // Check if existing results end with a newline, to avoid duplicating headers
    if (existingResults.endsWith("\n")) {
        // Remove header from new results if they are the same as the existing ones
        if (headers.length() > 0 && existingResults.contains(headers)) {
            updatedResults = updatedResults.replace(headers + "\n", "");
        }
        combinedResults = existingResults + updatedResults;
    } else {
        combinedResults = existingResults + "\n" + updatedResults;
    }
}

// Save the combined results back to the original file
File.saveString(combinedResults, outputFilePath);

// Remove the temporary file
File.delete(tempFilePath);
