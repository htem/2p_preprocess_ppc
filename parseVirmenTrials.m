function trialAlignedData = parseVirmenTrials(syncVirmenData,syncCaData)

%% Quick parsing of virmen output
% 1 time
% outputDATA(2,frameNumber) = mean(combinedVirmenFiles(2,itsUSE)); %2 x position
% outputDATA(3,frameNumber) = mean(combinedVirmenFiles(3,itsUSE)); %3 y position
% outputDATA(4,frameNumber) = mean(combinedVirmenFiles(4,itsUSE)); %4 heading angle
% outputDATA(5,frameNumber) = mean(combinedVirmenFiles(5,itsUSE)); %5 x velocity
% outputDATA(6,frameNumber) = mean(combinedVirmenFiles(6,itsUSE)); %6 y velocity
% outputDATA(7,frameNumber) = combinedVirmenFiles(7,itsUSE(end)); %7 vr.cuePos 
% outputDATA(8,frameNumber) = max(combinedVirmenFiles(8,itsUSE)); %8 vr.isReward
% outputDATA(9,frameNumber) = max(combinedVirmenFiles(9,itsUSE)); %9 vr.inITI
% outputDATA(10,frameNumber) = mode(combinedVirmenFiles(10,itsUSE));%10 vr.greyFac (always .5)
% outputDATA(11,frameNumber) = mode(combinedVirmenFiles(10,itsUSE)); %11 vr.breakFlag (always 0)
dF = syncCaData;

virmen_time = syncVirmenData(1,:);
x_pos = syncVirmenData(2,:);
y_pos = syncVirmenData(3,:);
y_vel = syncVirmenData(6,:);
mazeLength = max(y_pos);
cuePos = syncVirmenData(7,:);
isReward = syncVirmenData(8,:);
inITI = syncVirmenData(9,:);
frames = 1:length(virmen_time);
inTrial = 1-inITI;


% Delineate Trial Landmarks
[pks,ITIstarts] = findpeaks(inITI);
[pks,trialStarts] = findpeaks(inTrial); trialStarts = [1 trialStarts];
[pks,rewards] = findpeaks(isReward);
numRewards = length(rewards);
trialEnds = ITIstarts;

numTrials = length(trialStarts);
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
    cueType(i) = cuePos(trialEnds(1,i));
    % filter for erroneous trial starts where mouse is not at beginning of maze
    if y_pos(trialStarts(i)) > 50
        cueType(i) = nan;
        disp(['Removing erroneous trial ' num2str(i)])
    end
end

% fill in synched trial data per roi, using zscored dF for now
% frame 13 is trial onset with 12 frames before  1-12,13: 1-12, no cue, no image? Fr 13 cue goes on?
% frame 14 is running onset with 12 frames after 14,15-26
% frame 39 is cue offset (delay period onset) with 12 frames before and after 27-38,39,40-51: Fr 39 is when the mouse gets far enough to turn off the cue?
% frame 54 is trial end with 12 frames before and after Do you mean 64? 52-63,64,65-76: Fr 54 is end of trial. Is this when the mouse turns? Is the reward given here?


trial_dF = nan(size(dF,1),numTrials,76); % trial_dF(cell, trial, frames)
trial_type = nan(1,numTrials);
numVirVars = size(syncVirmenData,1);
trial_virmen = nan(numVirVars,numTrials,76);
for i = 1:numTrials
    trial_type(1,i) = cueType(i);
    if i == 1
        trial_virmen(:,i,1:12) = zeros(numVirVars,i,12);
        trial_virmen(:,i,13)= syncVirmenData(:,trialStarts(i));
    else
        trial_virmen(:,i,14:26) = syncVirmenData(:,runOnset(i):runOnset(i)+12);
        trial_virmen(:,i,27:51) = syncVirmenData(:,delayStart(i)-12:delayStart(i)+12);
        trial_virmen(:,i,52:76) = syncVirmenData(:,trialEnds(i)-12:trialEnds(i)+12);      
    end
    for roiIdx = 1:size(dF,1)
        if i == 1
            trial_dF(roiIdx,i,1:12) = zeros(1,1,12);
            trial_dF(roiIdx,i,13)= dF(roiIdx,trialStarts(i));
        else
            trial_dF(roiIdx,i,1:13) = dF(roiIdx,trialStarts(i)-12:trialStarts(i));
            trial_dF(roiIdx,i,1:13) = dF(roiIdx,trialStarts(i)-12:trialStarts(i));
        end   
        trial_dF(roiIdx,i,14:26) = dF(roiIdx,runOnset(i):runOnset(i)+12);
        trial_dF(roiIdx,i,27:51) = dF(roiIdx,delayStart(i)-12:delayStart(i)+12);
        trial_dF(roiIdx,i,52:76) = dF(roiIdx,trialEnds(i)-12:trialEnds(i)+12);
    end
end

trialAvg_dF = squeeze(mean(trial_dF,2));


% Export data struct
trialAlignedData = struct;
trialAlignedData.CaData = trial_dF;
trialAlignedData.trialType = trial_type;
trialAlignedData.virmenData = trial_virmen;
 