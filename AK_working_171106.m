%% 171106 AK working through 2P analysis code
% From Harvey Lab & Laura Driscoll

% Data hosted on Neurobiology server
% https://help.ubuntu.com/community/MountWindowsSharesPermanently
% Backed up on Orch?
mousePath = '/home/aaron/NeuroShare/Lee_Lab/Aaron/data/scanimage/AK013';
sessionList = dir([mousePath '/AK*']);
session = 18;% can loop over session here

%% Produce motion corrected movies
% Needs acq2P package https://github.com/HarveyLab/Acquisition2P_class
% and Harvey lab helper functions https://github.com/HarveyLab/helperFunctions
sessionPath=fullfile(mousePath,sessionList(session).name);
output_dir = fullfile(mousePath, '/sessions/',sessionList(session).name);
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
cd(sessionPath);
% Will ask you to select which video files to process
% Uses Acquision2P class
acq = Acquisition2P([],@SC2Pinit); 

acq.motionRefChannel = 2; % red channel 
acq.motionCorrect; % Produce motion corrected 

%load([sessionList(session).name '_001']);
%eval(['acq = ' sessionList(session).name '_001;'])
%acq.metaDataSI = metaDataSI;
channelNum = 1; % green?

sizeBit = acq.correctedMovies.slice(1).channel(1).size;
sizeArray = repmat(sizeBit,[size(acq.Movies,2) 1]);

% Here we assume there are 4 planes in the image
for sliceNum = 1:4
acq.correctedMovies.slice(sliceNum).channel(1).size = sizeArray;
acq.indexMovie(sliceNum,channelNum,output_dir);
acq.calcPxCov([],[],[],sliceNum,channelNum,output_dir);
save(fullfile(output_dir,sessionList(session).name),'acq','-v7.3')
end
% clear acq
% end

%% Compare motion-corrected and raw movies
movNum = 1;
castType = 'single';
sliceNum = 1;
channelNum = 1;%acq.motionRefChannel;
mov = readCor(acq,movNum,castType,sliceNum,channelNum);
rawMov = readRaw(acq,movNum,castType);
implay(cat(2,mov,rawMov(:,:,sliceNum:8:end))/1e3,30),

%% Precalculations for ROI selection and trace extraction
% First we will calculate pixel-pixel covariances. We want to do this for
% the GCaMP signal for some slice, so be sure the following variables match
% your settings (all arguments here are optional, but provided for clarity).

sliceNum = 1; %Choose a slice to analyze
channelNum = 1; %Choose the GCaMP channel
movNums = []; %this will default to all movie files
radiusPxCov = 11; %default, may need zoom level adjustment
temporalBin = 8; %default (tested w/ 15-30hz imaging), may need adjustment based on frame rate
writeDir = []; %empty defaults to the directory the object is saved in (the 'defaultDir')

% Now call the function:
acq.calcPxCov(movNums,radiusPxCov,temporalBin,sliceNum,channelNum,writeDir);

% Now we will create the binary file 'indexed movie' for the same
% slice and channel. 

acq.indexMovie(sliceNum,channelNum,writeDir);

% These functions do not automatically save the Acquisition object to disk.
% I suggest you overwrite the old object with the new one as you progress with each stage, so
% that you don't have to manually enter in filenames into fields of the
% 'old' object if you load the object again in the future. The
% acquisition2P class has a save method, with optional arguments, but
% default behavior is to overwrite the acqName file in defaultDir 

%help Acquisition2P.save
acq.save;

%% make sessionList_all variable for each mouse and save 
sessionList_all = cell(size(sessionList,1),1);
for s = 1:size(sessionList,1)
    sessionList_all{s,1} = sessionList(s,1).name;
end

%% ROI selection

% Now we can use these files to select ROIs. We can call the selectROIs
% method without input arguments, but again for clarity I will define some
% inputs here

% Use the built-in function meanRef to get a mean reference image, then
% process its square root with adaptive histogram equalization. Since this
% empirically produces nice looking images, this is actually the default
% behavior if no reference image is passed to the function call, so the
% code below is redundant but provided for demonstration

img = acq.meanRef;
img(img<0) = 0;
img(isnan(img)) = 0;
img = sqrt(img);
img = adapthisteq(img/max(img(:)));

% An alternative is to use an 'activity overview image', which has been
% precalculated in the calcPxCov call. This image highlights pixels which
% share strong correlations with neighboring pixels, and can be used
% independently or shared with an anatomical image, e.g.

actImg = acq.roiInfo.slice(sliceNum).covFile.activityImg;
% img = img/2 + actImg/2;

% Note that the reference image is only used for display purposes, and has no impact
% on the segmentation algorithm itself.

% Now start the ROI selection GUI. This tool is complex enough to have its
% own tutorial, located in the same folder as this file. Again, all
% arguments are optional, provided here just for clarity.
smoothWindow = 15; % Gaussian window with std = smoothWin/5, for displaying traces
excludeFrames = []; %List of frames that need to be excluded, e.g. if they contain artifacts

sliceNum = 1;
channeNum = 1; % green

acq.selectROIs(img,sliceNum,channelNum,smoothWindow,excludeFrames);

% Take the time to read through 'Using the ROI selection tool', and then
% select your cells. Once you've selected ROIs, be sure to save the acquisition again
acq.save;
