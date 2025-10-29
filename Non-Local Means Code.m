% Change the directory to Github
cd('C:\Users\walsh\Documents\GitHub\AGAROSE-HYDROGEL-TRENDS-USING-AI-ML')

% Labels of the file name used

concentrationLabel = '1%';             % e.g. '1%', '5%', etc.
regionLabel        = 'Region1';        % e.g. 'Region1', 'Region2', etc.

rawImageName       = 'Cropped 1%.tif';  % the original AFM image
filteredImageName  = ['NLM_' concentrationLabel '_' regionLabel '.tif'];
comparisonImageName = ['Comparison_' concentrationLabel '_' regionLabel '.tif'];

% Locate the Main Folder containing the images
imageFolder = fullfile('C:', 'Users', 'walsh', 'Desktop', ...
    'MSc Thesis', 'Dataset', 'AFM', concentrationLabel, ['Images (' concentrationLabel ')']);

% Build the full path to the input image
inputPath = fullfile(imageFolder, rawImageName);

% For convenience, also keep track of the current working directory.
% We will save copies in both places.
currentFolder = pwd;

% Create an "outputs" folder in the current directory (repo) if it does not exist
outputFolder = fullfile(currentFolder, 'outputs');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Read the original image from disk
I = imread(inputPath);

% Convert to double [0,1] for processing
I = im2double(I);

% Apply NLM denoising.
% This can be adjusted based on 'DegreeOfSmoothing'
denoised = imnlmfilt(I, 'DegreeOfSmoothing', 10);

%Compare both images visually
comparisonFig = figure;

% Show original (left) and filtered (right)
imshowpair(I, denoised, 'montage');

% Title
title(['Original vs NLM Filtered: ' concentrationLabel ' ' regionLabel], ...
    'Interpreter', 'none');

% Remove axes ticks etc. for a cleaner export
axis off;

% Hide toolbar so it doesn't appear in export
set(comparisonFig, 'Toolbar', 'none');
set(comparisonFig, 'menubar', 'none');

% Saves both the NLM image and a comparison between that and the original
imwrite(denoised, fullfile(imageFolder, filteredImageName));
imwrite(denoised, fullfile(currentFolder, filteredImageName));
imwrite(denoised, fullfile(outputFolder, filteredImageName));

exportgraphics(comparisonFig, fullfile(imageFolder, comparisonImageName));
exportgraphics(comparisonFig, fullfile(currentFolder, comparisonImageName));
exportgraphics(comparisonFig, fullfile(outputFolder, comparisonImageName));

% Message to explain the status after running the code

disp([' - Filtered image:      ' fullfile(imageFolder, filteredImageName)]);
disp([' - Filtered image copy: ' fullfile(currentFolder, filteredImageName)]);
disp([' - Filtered image repo: ' fullfile(outputFolder, filteredImageName)]);
disp([' - Comparison image:    ' fullfile(imageFolder, comparisonImageName)]);
disp([' - Comparison copy:     ' fullfile(currentFolder, comparisonImageName)]);
disp([' - Comparison repo:     ' fullfile(outputFolder, comparisonImageName)]);
