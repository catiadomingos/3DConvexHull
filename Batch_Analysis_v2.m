%% Batch analysis Convex Hull 3D
% Structure of interest: post-synaptic Shank2 cluster labeling (object) 
% Obtains ROI of post-synaptic structure (original cropped image)
% Binarizes isolated post-synaptic structure
% Calculates object and convex hull voxels, object and convex hull volumes and concave volume ratio 
% Plots: 
% - original cropped image,
% - binarized image,
% - 3D graphic of object,
% - 3D graphic of convex hull,
% - 3D graphic superimposition of object and convex hull 
% Exports previous representations + concave volume ratio, Shank2 volume (Fluorescence intensity), 
% EGFP intensity plotted in function of the distance to Shank2 center of mass
% and Shank2 center of mass (visualization) for each ROI per slide in a PPTX presentation.

%% Clear workspace and command window
clc;
clear;

%% Input parameters

% Astrocyte 2 - Re-analysis parameters
% xyp = 0.11686;
% zp = 0.4256000;
xfactor = 4.6288;
xyp = 0.11686 / xfactor;
zp = 0.4256000 / xfactor;

% Path to your Fiji installation and macro script, with the image path as an argument
fijiPath = 'D:\PhD_Catia Domingos\3. Programs\Fiji.app\ImageJ-win64.exe';  % Adjust for your operating system    
macroScript = 'D:\PhD_Catia Domingos\3. Programs\Fiji.app\macros\3D ConvexHull.ijm'; 

%% Batch Analysis

% Ask user to select the folder containing the images
imageFolder = uigetdir('D:\PhD_Catia Domingos\4. Experiments\7.  Expansion Microscopy + IHC (i)\2. Experiments (i)\1. ExM tripartite synapse final (i)\20240529 Analysis Convex Hull\1. Analysis', 'Select the folder containing the images to analyze');

%imageFolder = uigetdir('D:\', 'Select the folder containing the images to analyze');
if imageFolder == 0
    error('No folder selected. Please select a folder to proceed.');
end

% Define the raw results file path based on the selected image folder
rawResultsFilePath = fullfile(imageFolder, 'ConvexHull_Results_Raw.csv');

% Ask user to select the CSV file with cropping positions
[csvFile, csvPath] = uigetfile('*.csv', 'Select the CSV file with cropping positions', imageFolder);
if isequal(csvFile, 0)
    error('No CSV file selected. Please select a CSV file to proceed.');
end

% Read data from the selected CSV file (cropping positions for ROI)
csvFullPath = fullfile(csvPath, csvFile);
dataTable = readtable(csvFullPath, 'Delimiter', ',');

% Extract data from table (file name=fileTifList, x, y, and z cropping indices)
fileTifList = dataTable.FileTif;
x_refList = dataTable.x_ref;
y_refList = dataTable.y_ref;
z_refList = dataTable.z_ref;

% Calculate number of rows displaying in pptx (fixed for easiness of comparedness)
% Initialize an empty array to store intervals
intervals = zeros(1, length(z_refList));
% Loop through each cell in z_refList
for i = 1:length(z_refList)
    currentList = z_refList{i};           % Access the current cell's contents
    interval = currentList(end) - currentList(1) + 1;  % Calculate interval
    intervals(i) = interval;              % Store the interval
end
% Calculate the maximum interval
numCols = max(intervals);
%disp(['The maximum interval is: ', num2str(numCols)]);

% Convert x, y, and z cropping indices from string arrays to numeric
for file = 1:numel(fileTifList)
    x_refList{file} = eval(x_refList{file});
    y_refList{file} = eval(y_refList{file});
    z_refList{file} = eval(z_refList{file});
end

% Parameters
numSamples = numel(fileTifList);

% Calculate the PPCM (pixels per centimeter) from the micron resolution
micronsPerCm = 10000;  % 10,000 microns per centimeter
xyResolutionMicrons = 0.11686;
% Convert microns to PPCM
resolutionPPCM = micronsPerCm / xyResolutionMicrons;


%% Create a PowerPoint presentation
% Define Powerpoint presentation`s name based on images folder name 
% Extract the name of the selected folder
[~, folderName, ~] = fileparts(imageFolder);

% Construct the PowerPoint file name
pptFileName = ['ConvexHull_Analysis_', folderName, '.pptx'];

% Create the PowerPoint presentation with the constructed file name
ppt = mlreportgen.ppt.Presentation(pptFileName);
open(ppt);

%% Loop through each sample, process images and generate plots
for i = 1:numSamples
    % Process and binarize image
    % Debugging check
    if exist(fullfile(imageFolder, fileTifList{i}), 'file') ~= 2
        error('Image does not exist: %s', fileTifList{i});
    end
    
    try
        Image = fullfile(imageFolder, fileTifList{i});
        InfoImage = imfinfo(Image);
        mImage = InfoImage(1).Width;      % number of pixels x axis
        nImage = InfoImage(1).Height;     % number of pixels in y axis
        NumberImages = length(InfoImage); % number of images in z axis
    catch ME
        error('Error reading file: %s', ME.message);
    end
    
    % Cropping image (create an empty array and extract the x, y and z ref only)  
    Intensity_array = zeros(nImage, mImage, NumberImages);
    for m = 1:NumberImages
        Intensity_array(:, :, m) = imread(Image, 'Index', m);
    end
        
    % Isolate channel 2 
    Intensity_array_ch2 = Intensity_array(:, :, 2:3:NumberImages);
    
    % Debugging checks
    if isempty(Intensity_array_ch2)
        error('Intensity_array_ch2 is empty for file: %s', fileTifList{i});
    end
     
    % Ensure that indices are valid
    yIndex = x_refList{i} + 1;
    xIndex = y_refList{i} + 1;
    zIndex = z_refList{i};
    
    if any(xIndex < 1 | xIndex > size(Intensity_array_ch2, 1)) || ...
       any(yIndex < 1 | yIndex > size(Intensity_array_ch2, 2)) || ...
       any(zIndex < 1 | zIndex > size(Intensity_array_ch2, 3))
        error('Reference indices out of bounds for file: %s', fileTifList{i});
    end
    
    Shank2 = Intensity_array_ch2(xIndex, yIndex, zIndex);
       
    % Debugging check
    if isempty(Shank2)
        error('Shank2 is empty for file: %s', fileTifList{i});
    end

    %% Save the cropped image as Ch2_Shank2_"original file name"
    % File name: Concatenate the prefix with the original filename
    prefix = 'Ch2_Shank2_';  
    newFileName = [prefix, fileTifList{i}];
    
    % Open the Tiff object
    try
        t = Tiff(newFileName, 'w');
    
        % Get the size of the array
        [X, Y, Z] = size(Shank2);
               
        % Loop over the third dimension (e.g., image slices)        
        for sliceIndex = 1:Z    %length(z_refList{i})
            zIdx = z_refList{i}(sliceIndex);
            
            % Set Tiff tags for each slice
            setTag(t, 'ImageLength', X);
            setTag(t, 'ImageWidth', Y);
            setTag(t, 'Photometric', Tiff.Photometric.MinIsBlack);
            setTag(t, 'BitsPerSample', 8);
            setTag(t, 'SamplesPerPixel', 1);
            setTag(t, 'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
            setTag(t, 'Compression', Tiff.Compression.None);
            setTag(t, 'Software', 'MATLAB');
        
            % Add additional tags for original size
            setTag(t, 'XResolution', resolutionPPCM);
            setTag(t, 'YResolution', resolutionPPCM);
            setTag(t, 'ResolutionUnit',  Tiff.ResolutionUnit.Centimeter); 
                     
            % Write the image data
            %write(t, uint8(Shank2(:, :, zIndex)));
            t.write(uint8(Shank2(:, :, sliceIndex)));
                
            % If not the last slice, create a new directory for the next image
            if sliceIndex < Z    %length(z_refList{i})
               t.writeDirectory();
                
            %ifzIndex ~= Z
            %writeDirectory(t);
            end
        end
        
        % Close the Tiff object
        t.close;

    catch ME
          if exist('t', 'var') && isa(t, 'Tiff')
             t.close();
         end
         rethrow(ME);
    end
   

    %% Binarizes slices of a multi-page TIFF image and saves as a new TIFF file
    
    % Specify the filename of the multi-page TIFF image
    filename = newFileName;
    
    % Get information about the TIFF file
    info = imfinfo(filename);
    
    % Number of pages (slices) in the TIFF file
    numSlices = numel(info);
    
    % Create a new filename for the binarized image
    % Concatenate the prefix with the original filename
    outputFilename = ['Binarized_', newFileName];
    
    % Loop over each slice
    for k = 1:numSlices
        % Read the k-th slice
        slice = imread(filename, 'Index', k);
        
        % Convert to grayscale if the slice is in color
        if size(slice, 3) == 3
            slice = rgb2gray(slice);
        end
    
    % Ensure the image data type is one of the supported types
    slice = im2double(slice); % Convert to double for consistency
    
    % Binarize the slice
    binarySlice = imbinarize(slice);
    
    % Remove small noise and keep only the largest connected component in a binarized image
    % Identify connected components
    cc = bwconncomp(binarySlice);
    
    % Measure the size of each component
    numPixels = cellfun(@numel, cc.PixelIdxList);
    
    % Find the largest component
    [~, largestComponentIdx] = max(numPixels);
    
    % Create a new binary image with only the largest component
    largestComponentImage = false(size(binarySlice));
    largestComponentImage(cc.PixelIdxList{largestComponentIdx}) = true;
    
    % Convert the binary slice to uint8 format
    binarySlice = uint8(largestComponentImage) * 255; % Scale to 0-255 range
    
    % Open the TIFF file for writing
    if k == 1
        % For the first slice, create a new file
        t = Tiff(outputFilename, 'w');
    else
        % For subsequent slices, append to the existing file
        t = Tiff(outputFilename, 'a');
    end
    
    % Write the image data to the TIFF file
    tagstruct.ImageLength = size(binarySlice, 1);
    tagstruct.ImageWidth = size(binarySlice, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.Software = 'MATLAB';
    t.setTag(tagstruct);
    
    % Add additional tags for original size
    setTag(t, 'XResolution', resolutionPPCM);
    setTag(t, 'YResolution', resolutionPPCM);
    setTag(t, 'ResolutionUnit',  Tiff.ResolutionUnit.Centimeter); 

    % Write the image data
    t.write(binarySlice);
    
    % Close the current TIFF directory
    t.close();
    end       

     
    %% Calculates 3D Convex hull properties in Fiji: 
    % Define the image path and command to run Fiji
    imagePath = fullfile(imageFolder, ['Binarized_Ch2_Shank2_', fileTifList{i}]);   

    % Command to run Fiji and execute the macro
    command = sprintf('"%s" --headless -macro "%s" "%s"', fijiPath, macroScript, imagePath);
       
    try
        % Execute the command
        [status, cmdout] = system(command);
        
        if status == 0
            % disp('Fiji macro executed successfully.');
            % disp(cmdout);                         
        else
            error('Error executing Fiji macro. Status: %d. Output: %s', cmdout);
        end
    
    catch ME
        disp(['Error: ', ME.message]);
    end    

    
    %% Creating a figure for pptx: Method validation
    % with process of image binarization and convex hull representation
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Set up the figure with an appropriate layout
    % Set numCols to be either numSlices or at least 3
%     if numSlices < 3
%         numRows = 3;
%     else
%         numRows = numSlices;
%     end
%    numCols = 4;
    numRows = 4;
      
    % Set spacing values
    gap = [0.01,0.001]; % Example [gap_h, gap_w] between subplots
    xSpacing = 0.03; % Example horizontal margin
    ySpacing = 0; % Example vertical margin
    
    % Plot cropped slices
    for p = 1:numSlices
        croppedSlice = imread(filename, 'Index', p);
        threshold = 50;
        adjustedCroppedSlice = imadjust(croppedSlice, [0 threshold/255], []);
        binarizedSlice = imread(outputFilename, 'Index', p);
        
        % Plot the slices in a subplot (cropped)
        index = 1 + (p-1) * numCols;
        subtightplot(numRows, numCols, index, gap, xSpacing, ySpacing);
        imshow(croppedSlice)      
        % Add title only for the top row
        hold on;
        if index < (numCols + 1)
           title('Cropped Slices', 'FontSize', 16, 'FontWeight', 'bold');
        end
        hold off;
        
               
        % Plot the slices in a subplot (adjusted)
        index = 2 + (p-1) * numCols;
        subtightplot(numRows, numCols, index, gap, xSpacing, ySpacing);
        imshow(adjustedCroppedSlice)
        % Add title only for the top row
        hold on;
        if index < (numCols + 1)
           title('Threshold adjusted', 'FontSize', 16, 'FontWeight', 'bold');
        end
        hold off;
        
                
        % Plot the slice in a subplot (binarized)
        index = 3 + (p-1) * numCols;
        subtightplot(numRows, numCols, index, gap, xSpacing, ySpacing);
        imshow(binarizedSlice);
        % Add title only for the top row
        hold on;
        if index < (numCols + 1)
           title('Binarized', 'FontSize', 16, 'FontWeight', 'bold');
        end
        hold off;
        
        %title(['Binarized Slice ', num2str(p)]);
               
    end
              
    %% Plots representations of convex hull volume and object volume
    
    % Get information about the TIFF file
    InfoImage = imfinfo(outputFilename);
    numImages = length(InfoImage); % Number of slices in the TIFF file
    
    % Read the first image to get its dimensions
    firstImage = imread(outputFilename, 'Index', 1);
    [height, width] = size(firstImage);
    
    % Preallocate a 3D array for the images
    image3D = zeros(height, width, numImages);
    
    % Read each slice and store it in the 3D array
    for n = 1:numImages
        image3D(:, :, n) = imread(outputFilename, 'Index', n);
    end
    % Define spacing values
%     xSpacing = 0.05; % Horizontal space between subplots
%     ySpacing = 0.05; % Vertical space between subplots

    % Display the binary image volume, convex hull, and combined plot
    % Plot binary image volume
    index = numCols*1;
    %gap = [0.01, 0.1];
    subtightplot(numRows, numCols, index, gap, xSpacing, ySpacing);    
    [x, y, z] = meshgrid(1:width, 1:height, 1:numImages);
    % Convert voxel coordinates to physical dimensions
    x = x * xyp;
    y = y * xyp;
    z = z * zp;
    
    p1 = patch(isosurface(x, y, z, image3D, 0.5));
    isonormals(x, y, z, image3D, p1);
    p1.FaceColor = 'red';
    p1.EdgeColor = 'none';
    daspect([1, 1, 1]); % Ensure equal aspect ratio
    view(3); 
    camlight; 
    lighting phong;
    title('Binary Image Volume','FontSize', 16);
    xlabel('X (µm)');
    ylabel('Y (µm)');
    zlabel('Z (µm)');
    
    hold on;
    grid on; % Turns on the grid lines

    % Calculate convex hull on scaled coordinates
    [rows, cols, slices] = ind2sub(size(image3D), find(image3D));
    scaled_points = [cols * xyp, rows * xyp, slices * zp]; % Adjusted x and y for correct orientation
    [K, convex_hull_volume] = convhulln(scaled_points);             

    % Plot convex hull
    index = numCols*2;
    subtightplot(numRows, numCols, index, gap, xSpacing, ySpacing);
    scatter3(scaled_points(:,1), scaled_points(:,2), scaled_points(:,3), 'b.');
    hold on;
    trisurf(K, scaled_points(:,1), scaled_points(:,2), scaled_points(:,3), 'FaceColor', 'cyan', 'FaceAlpha', 0.3);
    xlabel('X (µm)');
    ylabel('Y (µm)');
    zlabel('Z (µm)');
    title('3D Convex Hull','FontSize', 14);
    daspect([1, 1, 1]); % Ensure equal aspect ratio
    view(3); 
    camlight; 
    lighting phong;
    hold off;
    
    % Combined plot with both binary image volume and convex hull
    index = numCols*3;
    subtightplot(numRows, numCols, index, gap, xSpacing, ySpacing);
    p1 = patch(isosurface(x, y, z, image3D, 0.5));
    isonormals(x, y, z, image3D, p1);
    p1.FaceColor = 'red';
    p1.EdgeColor = 'none';
    daspect([1, 1, 1]); % Ensure equal aspect ratio
    view(3); 
    camlight; 
    lighting phong;
    %title('Superimposition');
    xlabel('X (µm)');
    ylabel('Y (µm)');
    zlabel('Z (µm)');
    
    hold on;
    grid on; % Turns on the grid lines
    scatter3(scaled_points(:,1), scaled_points(:,2), scaled_points(:,3), 'b.');
    trisurf(K, scaled_points(:,1), scaled_points(:,2), scaled_points(:,3), 'FaceColor', 'cyan', 'FaceAlpha', 0.3);
    hold off;
    
    % Add Name + Concave volume ratio
    % Check if the file exists
    if exist(rawResultsFilePath, 'file') ~= 2
        error('The ConvexHull_Results_Raw.csv file does not exist.');
    end
    % Read the data from the CSV file into a table
    convexHullData = readtable(rawResultsFilePath);
    % Read the last row of the CSV file
    lastRow = convexHullData(end, :);
    % Extract text from column 10
    column10Text = lastRow{1, 10};   % Assuming column 10 is of numeric type
    ConcaveVolumeRatio = sprintf('%.3f', column10Text);
    % Prepare the text content
    textboxText1 = fileTifList{i};
    textboxText2 = sprintf('Concave Volume Ratio: %s');
    textboxText3 = sprintf(ConcaveVolumeRatio);
    % Plot the text string
    index = numCols*4;
    subtightplot(numRows, numCols, index, gap, xSpacing, ySpacing);
    % Remove the axes (make it invisible)
    axis off;
    % Add your text to the figure
    text(0.2, 0.6, textboxText1, 'FontSize', 20, 'HorizontalAlignment', 'left','Color', 'red');
    text(0.2, 0.4, textboxText2, 'FontSize', 20, 'HorizontalAlignment', 'left','Color', 'red');
    text(0.2, 0.2, textboxText3, 'FontSize', 35, 'HorizontalAlignment', 'left','Color', 'red');
    % Save the figure as an image
    imgFile = fullfile(tempdir, ['plot_' num2str(i) '.png']);
    saveas(gcf, imgFile);
    close(gcf);

    %% Laying out the slide content
    % Add a slide with a blank layout to the PowerPoint presentation
    slide = add(ppt, 'Title and Content');           

    % Set the title of the slide
    titleText = ['Sample: ', fileTifList{i}];
    titlePlaceholder = find(slide, 'Title');
    replace(titlePlaceholder, titleText);

    % Add the image to the slide
    contentPlaceholder = find(slide, 'Content');
    img = mlreportgen.ppt.Picture(imgFile);
    replace(contentPlaceholder, img);             
end

%% Save and close the PowerPoint presentation
close(ppt);

%% Save Convex Hull measurements as Excel file

%% Convert CSV to Excel

% Define the Excel file path based on the selected image folder
rawExcelFilePath = fullfile(imageFolder, 'ConvexHull_Results_Processed.xlsx');

% Read the CSV file into a table
dataTable = readtable(rawResultsFilePath);
% Write the table to an Excel file
writetable(dataTable, rawExcelFilePath);

%% Process and format data
%% Correct headers for rawResults
% Load data from the Excel file
filename = 'ConvexHull_Results_Processed.xlsx'; 
dataTable = readtable(filename);

% Modify column name
dataTable.Properties.VariableNames{'Var1'} = 'Sample';

% Save the modified table 
writetable(dataTable, filename);

%% Auto-fit columns and align text in Excel
% Start Excel
excel = actxserver('Excel.Application');
excel.Visible = true;  % Set to true for debugging, false for silent operation

% Open the workbook
workbook = excel.Workbooks.Open(fullfile(pwd, filename));

% Access the first sheet
sheet = workbook.Sheets.Item(1);

% Auto-fit columns
sheet.Columns.AutoFit();

% Center align text in all cells
range = sheet.UsedRange;
range.HorizontalAlignment = -4108;  % -4108 corresponds to 'xlCenter'

% Explicitly save the workbook
workbook.Save();

% Close the workbook and quit Excel
workbook.Close();
excel.Quit();
delete(excel);