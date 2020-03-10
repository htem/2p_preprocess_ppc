%% 180826 AK working through 2P analysis code
% Mouse AK001, data from Loki Rig

% From Harvey Lab & Laura Driscoll

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
session = find(contains({sessionList.name},'AK001_170626'));% can loop over session here


%% Initialize Acq2P object
sessionPath=fullfile(mousePath,sessionList(session).name);
output_dir = fullfile(mousePath, '\sessions\',sessionList(session).name);
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
cd(sessionPath);
acq = Acquisition2P([],@ak_init); 

%% Prohttps://github.com/epnev/constrained-foopsiduce motion corrected movies

acq.motionRefChannel = 2; % red channel 
acq.motionCorrect; % Produce motion corrected 
%acq.linkCorrectedMovies;
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
%load('Z:\Lee_Lab\Aaron\data\scanimage\AK001\sessions\AK001_170626\AK001_170626.mat')
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
sliceNum = 1; %Choose a slice to analyze
channelNum = 1; %Choose the GCaMP channel
actImg = acq.roiInfo.slice(sliceNum).covFile.activityImg;
img = img/2 + actImg/2;

% Note that the reference image is only used for display purposes, and has no impact
% on the segmentation algorithm itself.

% Now start the ROI selection GUI. This tool is complex enough to have its
% own tutorial, located in the same folder as this file. Again, all
% arguments are optional, provided here just for clarity.
smoothWindow = 15; % Gaussian window with std = smoothWin/5, for displaying traces
excludeFrames = []; %List of frames that need to be excluded, e.g. if they contain artifacts


channeNum = 1; % green

acq.selectROIs(img,sliceNum,channelNum,smoothWindow,excludeFrames);

% Take the time to read through 'Using the ROI selection tool', and then
% select your cells. Once you've selected ROIs, be sure to save the acquisition again
acq.save;
% AK171117, name of Acq2P object is AK001 after save, not acq, leads to
% some bookkeeping issues.
%load('Z:\Lee_Lab\Aaron\data\scanimage\AK001\AK001_170626\AK001.mat')

%% get dF
% AK171117 replacing script get_dF.m
sliceNum = 1;
channelNum = 1;
[dF,traces,rawF,roiList] = extractROIsBin(acq,[],sliceNum,channelNum);
z_dF = zscore(dF,0,2);

% Plot dF
figure; hold on;
frames = 1:length(dF);
time = frames * acq.metaDataSI.SI4.fastZPeriod;
for i = 1:size(dF,1)
    plot(time,dF(i,:)+2*i);
end
xlabel('time (sec)');
ylabel('dF/F');

%% Deconvolution 
%[c,b,c1,g,sn,sp] = constrained_foopsi(dF);

%% Laura's get_dF code, needs some work to get it to run
%{

cd('Z:\Lee_Lab\Aaron\data\scanimage\AK001\sessions')
mother_dir = pwd;
sessionList = dir('AK*');

cd('Z:\Laura\code_workspace\AK001\sessions')
mother_dir = pwd;et dF, etc

sessionList = dir('AK*');
for session = 12%1:size(sessionList,1)%3%15%37%%%%LD186_10%%%LD183_36
    sessionDir = fullfile(mother_dir,sessionList(session).name);
    cd(sessionDir)
    
    clear acq
    acqName = dir('selected_rois*');
    if ~isempty(acqName)
        load(acqName.name)
        for sliceNum = 2:4
            filename = sprintf('%s_slice%02.0f_deconv', 'dF_struct', sliceNum);
            dF_struct.filename = filename;
            
            [~,covName,~] = fileparts(acq.roiInfo.slice(sliceNum).covFile.fileName);
            [~,indName,~] = fileparts(acq.indexedMovie.slice(sliceNum).channel.fileName);
            acq.roiInfo.slice(sliceNum).covFile.fileName = ...
                fullfile(sessionDir,[covName '.bin']);
            acq.indexedMovie.slice(sliceNum).channel.fileName = ...
                fullfile(sessionDir,[indName '.bin']);
            
            [dF_struct.dF, dF_struct.traces, dF_struct.rawF, ...
                dF_struct.roi, dF_struct.traceNeuropil, dF_struct.rawF_nuc] = ...
                extractROIsBin(acq,[],sliceNum,1);
            
            dF_struct.dF_adj = baseline_adj(dF_struct.dF);
            
            [c,b,c1,g,sn,sp] = run_constrained_foopsi(dF_struct.dF_adj);
            
            for cel = find(~isnan(sum(sp,2)))'
                maxPulse = max(impulseAR(g{cel,1}));
                unitaryDF = 0.1;
                deconvInSpikes(cel,:) = sp(cel,:)*maxPulse/unitaryDF;
            end
            
            save(fullfile('Z:\Laura\code_workspace',acq.acqName(1:5),'sessions',...
                acq.acqName(1:12),filename),'dF_struct','c','b','c1','g','sn','sp');
        end
    end
end
%}
%% Line up pclamp, virmen and scanimage data
% for now, just one session passed
% need to update to accept list of sessions
syncVirmenData = lineUpAll(sessionList_all,session);


%% Quick parsing of virmen output
% %1 time
%         outputDATA(2,frameNumber) = mean(combinedVirmenFiles(2,itsUSE)); %2 x position
%         outputDATA(3,frameNumber) = mean(combinedVirmenFiles(3,itsUSE)); %3 y position
%         outputDATA(4,frameNumber) = mean(combinedVirmenFiles(4,itsUSE)); %4 heading angle
%         outputDATA(5,frameNumber) = mean(combinedVirmenFiles(5,itsUSE)); %5 x velocity
%         outputDATA(6,frameNumber) = mean(combinedVirmenFiles(6,itsUSE)); %6 y velocity
%         outputDATA(7,frameNumber) = combinedVirmenFiles(7,itsUSE(end)); %7 vr.cuePos 2 or 3
%         outputDATA(8,frameNumber) = max(combinedVirmenFiles(8,itsUSE)); %8 vr.isReward
%         outputDATA(9,frameNumber) = max(combinedVirmenFiles(9,itsUSE)); %9 vr.inITI
%         outputDATA(10,frameNumber) = mode(combinedVirmenFiles(10,itsUSE));%10 vr.greyFac (always .5)
%         outputDATA(11,frameNumber) = mode(combinedVirmenFiles(10,itsUSE)); %11 vr.breakFlag (always 0)

virmen_time = syncVirmenData(1,:);
x_pos = syncVirmenData(2,:);
y_pos = syncVirmenData(3,:);
y_vel = syncVirmenData(6,:);
mazeLength = max(y_pos);
cuePos = syncVirmenData(7,:);
isReward = syncVirmenData(8,:);
inITI = syncVirmenData(9,:);
frames = 1:length(time);
inTrial = 1-inITI;

%% Delineate Trial Landmarks
[pks,ITIstarts] = findpeaks(inITI);
[pks,trialStarts] = findpeaks(inTrial); trialStarts = [1 trialStarts];
[pks,rewards] = findpeaks(isReward);
numRewards = length(rewards);
trialEnd = ITIstarts;
numTrials = length(ITIstarts);

for i = 1:numTrials
    if i < numTrials
        success(i) = max(isReward(trialStarts(i):trialStarts(i+1)));
        delayStart(i) = find(y_pos(trialStarts(i):trialStarts(i+1)) > 0.5*mazeLength,1)+trialStarts(i)-1;
        runOnset(i) = find(y_vel(trialStarts(i):trialStarts(i+1)) > 1 ,1) + trialStarts(i) - 1; % what is the right y-vel cutoff?
    else
        success(i) = max(isReward(trialStarts(i):end));
        delayStart(i) = find(y_pos(trialStarts(i):end) > 0.5*mazeLength,1);
        runOnset(i) = find(y_vel(trialStarts(i):end) > 1,1) + trialStarts(i) - 1;
    end
    cueType(i) = cuePos(ITIstarts(i));
   
end




%% fill in synched trial data per roi, using zscored dF for now
trial_dF = nan(size(dF,1),numTrials,76); % trial_dF(cell, trial, frames)
for i = 1:numTrials
    for roiIdx = 1:size(dF,1)
        if i == 1
            trial_dF(roiIdx,i,1:12) = zeros(1,1,12);
            trial_dF(roiIdx,i,13)= dF(roiIdx,trialStarts(i));
            trial_dF(roiIdx,i,14:26) = dF(roiIdx,runOnset(i):runOnset(i)+12);
            trial_dF(roiIdx,i,27:51) = dF(roiIdx,delayStart(i)-12:delayStart(i)+12);
            trial_dF(roiIdx,i,52:76) = dF(roiIdx,trialEnd(i)-12:trialEnd(i)+12);
        else
            trial_dF(roiIdx,i,1:13) = dF(roiIdx,trialStarts(i)-12:trialStarts(i));
            trial_dF(roiIdx,i,14:26) = dF(roiIdx,runOnset(i):runOnset(i)+12);
            trial_dF(roiIdx,i,27:51) = dF(roiIdx,delayStart(i)-12:delayStart(i)+12);
            trial_dF(roiIdx,i,52:76) = dF(roiIdx,trialEnd(i)-12:trialEnd(i)+12);
        end
    end
end

trialAvg_dF = squeeze(mean(trial_dF,2));
figure; imagesc(trialAvg_dF); colormap gray;

% Plot trial avg dF
figure; hold on;
frames = 1:size(trialAvg_dF,2);
time = frames;% * acq.metaDataSI.SI4.fastZPeriod;
for i = 1:size(dF,1)
    plot(time,trialAvg_dF(i,:));
end
xlabel('time (frames)');
ylabel('dF/F');
    %}
    
    %
% 
% frame 13 is trial onset with 12 frames before  1-12,13: 1-12, no cue, no image? Fr 13 cue goes on?
% frame 14 is running onset with 12 frames after 14,15-26
% frame 39 is cue offset (delay period onset) with 12 frames before and after 27-38,39,40-51: Fr 39 is when the mouse gets far enough to turn off the cue?
% frame 54 is trial end with 12 frames before and after Do you mean 64? 52-63,64,65-76: Fr 54 is end of trial. Is this when the mouse turns? Is the reward given here?


 %% Sequence Plot (for intro)


% For now, average over all sessions
trialMean = trialAvg_dF;
zTrialMean = zscore(trialMean,0,2);
%figure; imagesc(zTrialMean(:,:,1));

% find peaks
[peaks, tPeaks] = max(trialMean(:,:,:),[],2); tPeaks = squeeze(tPeaks);
[temp,cellOrder] = sort(tPeaks,1);

figure; colormap jet; % Just a figure to look at sequences averaged over all sessions
imagesc(zTrialMean(cellOrder,:)); title('Sorted');

line([13 0],[-1 -1],'Linewidth',2);
h(1) = text(12.5,-2,'Trial Begin');
line([14 26],[-1 -1],'Linewidth',2);
h(2) = text(14.5,-2,'Running Onset');
line([27 51],[-1 -1],'Linewidth',2);
h(3) = text(39,-2,'Cue Offset');
line([52 76],[-1 -1],'Linewidth',2);
h(4) = text(64,-2,'Trial End');
set(h,'Rotation',90);
ylim([-20 size(zTrialMean,1)]);
% subplot(2,2,2);imagesc(trialMean(cellOrderR,12:64,1));title('L sorted by R');
% subplot(2,2,3);imagesc(trialMean(cellOrderL,12:64,2));title('R sorted by L');
% subplot(2,2,4);imagesc(trialMean(cellOrderL,12:64,1));title('L sorted by L');


%% Sequence Plot separated R, L

correctR = find(cueType == 2 & success);
correctL = find(cueType == 3 & success);

RtrialAvg_dF = squeeze(mean(trial_dF(:,correctR,:),2));
LtrialAvg_dF = squeeze(mean(trial_dF(:,correctL,:),2));

zTrialMean_R = zscore(RtrialAvg_dF,0,2);
zTrialMean_L = zscore(LtrialAvg_dF,0,2);
%figure; imagesc(zTrialMean(:,:,1));


% find peaks
[peaksR, tPeaksR] = max(zTrialMean_R(:,:,:),[],2); tPeaksR = squeeze(tPeaksR);
[peaksL, tPeaksL] = max(zTrialMean_L(:,:,:),[],2); tPeaksL = squeeze(tPeaksL);
[tempR,cellOrderR] = sort(tPeaksR,1);
[tempL,cellOrderL] = sort(tPeaksL,1);

figure; colormap jet; % Just a figure to look at sequences averaged over all sessions
%imagesc(zTrialMean(cellOrder,:)); title('Sorted');
subplot(2,2,1);imagesc(zTrialMean_R(cellOrderR,:));title('R sorted by R');
subplot(2,2,2);imagesc(zTrialMean_L(cellOrderR,:));title('L sorted by R');
subplot(2,2,3);imagesc(zTrialMean_R(cellOrderL,:));title('R sorted by L');
subplot(2,2,4);imagesc(zTrialMean_L(cellOrderL,:));title('L sorted by L');


%%
% Plot trial avg dF
figure; hold on;
frames = 1:size(trialAvg_dF,2);
time = frames(12:64);% * acq.metaDataSI.SI4.fastZPeriod;
for i = 1:size(dF,1)
    plot(time,zTrialMean(cellOrder(i),12:64)+0.3*i);
end
xlabel('time (frames)');
ylabel('dF/F, zscored');
%%
 actCellsR = find(selectIndex>0);
 actCellsL = find(selectIndex<0);
 
 cellOrderRSelect = cellOrderR(ismember(cellOrderR,actCellsR));
 cellOrderLSelect = cellOrderL(ismember(cellOrderL,actCellsL));

 figure; colormap jet; % Just a figure to look at sequences averaged over all sessions
subplot(2,1,1);imagesc(trialMean(cellOrderRSelect,12:64,2));title('R sorted by R');
subplot(2,1,2);imagesc(trialMean(cellOrderLSelect,12:64,1));title('L sorted by L');

%[temp, cellOrderRSelect] = sort(tPeaks(selectIndex>.2,2),1);
%[temp, cellOrderLSelect] = sort(tPeaks(selectIndex<-.2,1),1);
%colormap gray
