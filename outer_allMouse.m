% ATK 180709
% outer wrapper for 2P analysis 
% Github Repo: https://github.com/lauradriscoll/pre_process_shared
% AKdev branch

% Requires 

% (1) Harvey Lab Acq2P class
% https://github.com/HarveyLab/Acquisition2P_class
% (2) Harvey Lab Helper Functions
% https://github.com/HarveyLab/helperFunctions

mouseName = 'AK001';
mousePath = ['Z:\Lee_Lab\Aaron\data\scanimage\' mouseName];
sessionList = dir(fullfile(mousePath, [mouseName '*']));
session = find(contains({sessionList.name},'AK001_170625'));% can loop over session here

%% Initialize Acq2P object
sessionPath=fullfile(mousePath,sessionList(session).name);
output_dir = fullfile(mousePath, '\sessions\',sessionList(session).name);
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
cd(sessionPath);
acq = Acquisition2P([],@ak_init); 


%% Produce motion corrected movies
tic
acq.motionRefChannel = 2; % red channel 
acq.motionCorrect; % Produce motion corrected 
toc
%load([sessionList(session).name '_001']);
%eval(['acq = ' sessionList(session).name '_001;'])
%acq.metaDataSI = metaDataSI;
channelNum = 1; % green?

sizeBit = acq.correctedMovies.slice(1).channel(1).size;
sizeArray = repmat(sizeBit,[size(acq.Movies,2) 1]);
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
