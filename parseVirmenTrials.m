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

numTrials = length(trialEnds);
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
    % filter for erroneous trials
    if delayStart(i) < 12
        cueType(i) = nan;
        disp(['Removing erroneous trial ' num2str(i)])
    elseif y_pos(trialStarts(i)) > 50
        cueType(i) = nan;
        disp(['Removing erroneous trial ' num2str(i)])
    elseif max(abs(x_pos(trialStarts(i):trialStarts(i)+12))) > 0.2 
        cueType(i) = nan;
        disp(['Removing erroneous trial ' num2str(i)])
    elseif max(abs(x_pos(delayStart(i)-12:delayStart(i)+12))) > 0.2 
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
    trial_virmen(:,i,:) = nan(numVirVars,76);
    if trialStarts(i)>12
        trial_virmen(:,i,1:13) = syncVirmenData(:,trialStarts(i)-12:trialStarts(i));
    end
    if runOnset(i)>12
        trial_virmen(:,i,14:26) = syncVirmenData(:,runOnset(i):runOnset(i)+12);
    end
    if delayStart(i)>12 
        trial_virmen(:,i,27:51) = syncVirmenData(:,delayStart(i)-12:delayStart(i)+12);
    end
    if trialEnds(i)>12
        trial_virmen(:,i,52:76) = syncVirmenData(:,trialEnds(i)-12:trialEnds(i)+12);
    end
    for roiIdx = 1:size(dF,1)
        trial_dF(roiIdx,i,:) = nan(76,1);
        if trialStarts(i)>12
            trial_dF(roiIdx,i,1:13) = dF(roiIdx,trialStarts(i)-12:trialStarts(i));
        end
        if runOnset(i)>12
            trial_dF(roiIdx,i,14:26) = dF(roiIdx,runOnset(i):runOnset(i)+12);
        end
        if delayStart(i)>12
            trial_dF(roiIdx,i,27:51) = dF(roiIdx,delayStart(i)-12:delayStart(i)+12);
        end
        if trialEnds(i)>12
            trial_dF(roiIdx,i,52:76) = dF(roiIdx,trialEnds(i)-12:trialEnds(i)+12);
        end
    end
end

trialAvg_dF = squeeze(mean(trial_dF,2));


% Export data struct
trialAlignedData = struct;
trialAlignedData.CaData = trial_dF;
trialAlignedData.trialType = trial_type;
trialAlignedData.virmenData = trial_virmen;
trialAlignedData.numCells = size(syncCaData,1);

%% Calculate derived quantities
trialTypes = {'newL_trials','bR_trials','wL_trials','newR_trials'};
for i = 1:length(trialTypes)
    % Ca activity and virmen activity for each trial type individually
    trialAlignedData.(trialTypes{i}) = struct;
    trialAlignedData.(trialTypes{i}).Ca = trialAlignedData.CaData(:,trialAlignedData.trialType==i,:);
    trialAlignedData.(trialTypes{i}).virmen = trialAlignedData.virmenData(:,trialAlignedData.trialType==i,:);

    % Calculate trial average activities
    trialAlignedData.(trialTypes{i}).Ca_trialMean = squeeze(mean(trialAlignedData.(trialTypes{i}).Ca,2));
end

%% Sanity check with trialwise behavioral data
beh_idx = 2;

figure; hold on;
pltC = {'m-','b-','r-','c-'};
for j = 1:length(trialTypes)
    if size(trialAlignedData.(trialTypes{j}).virmen,2) > 0
        for i = 1:size(trialAlignedData.(trialTypes{j}).virmen,2)
            plot(squeeze(trialAlignedData.(trialTypes{j}).virmen(beh_idx,i,:)),pltC{j})
        end
    end
end
%% Calculate R/L selectivity indices
thisTrial = 13:76;
% exclude ITI-before (which is unrelated to this trial), 13:76
r = nanmean(trialAlignedData.bR_trials.Ca_trialMean(:,thisTrial),2);
l = nanmean(trialAlignedData.wL_trials.Ca_trialMean(:,thisTrial),2);
trialAlignedData.RL_selectIdx = (r - l)./(r + l);
figure; histogram(trialAlignedData.RL_selectIdx)
% definition fo cellSelectIdx from Harvey 2012 

% sanity-check activity for very selective cells
figure; hold on; 
subplot(1,2,1);
plot(thisTrial,trialAlignedData.bR_trials.Ca_trialMean(trialAlignedData.RL_selectIdx>0.4,thisTrial)')
ylim([0,25]);
subplot(1,2,2); 
plot(thisTrial,trialAlignedData.wL_trials.Ca_trialMean(trialAlignedData.RL_selectIdx>0.4,thisTrial)')
ylim([0,25]);



%% Calculate pairwise correlation matrices

trialAlignedData.corr_all = corrcoef(syncCaData');
for i = 1:length(trialTypes)
    if size(trialAlignedData.(trialTypes{i}).virmen,2) > 0
        disp(trialTypes{i})
        trialAlignedData.(trialTypes{i}).corrcoef = corrcoef(trialAlignedData.(trialTypes{i}).Ca_trialMean','rows','complete');
    end
end
% sanity check by comparing
figure;
plot(trialAlignedData.corr_all(:),trialAlignedData.bR_trials.corrcoef(:),'.');
xlabel('total correlation'); ylabel('bR trial correlation');

figure;
plot(trialAlignedData.bR_trials.corrcoef(:),trialAlignedData.wL_trials.corrcoef(:),'.');
xlabel('bR trials correlations'); ylabel('wL trials correlation');

%% Calculate timing metrics


for i = 1:length(trialTypes)
    if size(trialAlignedData.(trialTypes{i}).virmen,2) > 0
        Ca = trialAlignedData.(trialTypes{i}).Ca_trialMean;
        % time center-of-mass 
        trialAlignedData.(trialTypes{i}).tCOM = (Ca(:,thisTrial)*thisTrial')./sum(Ca(:,thisTrial),2);
        % time of max
        tMax = nan(trialAlignedData.numCells);
        for cid = 1:trialAlignedData.numCells
            tMax(cid) = find(Ca(cid,:) == max(Ca(cid,thisTrial)),1,'last');
        end
        trialAlignedData.(trialTypes{i}).tMax = tMax;
    end
end
% sanity check by plotting high COM 
figure;
plot(trialAlignedData.bR_trials.Ca_trialMean(trialAlignedData.bR_trials.tCOM>48,:)');

% sanity check by plotting high tmax
figure;
plot(trialAlignedData.bR_trials.Ca_trialMean(trialAlignedData.bR_trials.tMax>70,:)');

%% Calculate activity timing indices
% ATK sanity check doesn't seem right, hold off for now
%{
% Define time blocks
cueBlockEarly = 14:26; % 14 is running onset + 12 frames after
cueBlockLate = 27:38; % 12 frames before cue offset (frame 39)
delayBlockEarly = 39:51; % 39 is cue offset + 12 frames after
delayTurnBlock = 52:64; % 12 frames before end of trial (turn a certain amt)
turnITIblock = 65:76; % Trial end (reward given) an 12 frames after (dark, ITI)

for i = 1:length(trialTypes)
    if size(trialAlignedData.(trialTypes{i}).virmen,2) > 0
        disp(['type ' trialTypes{i}])
        meanCueEarly = squeeze(mean(trialAlignedData.(trialTypes{i}).Ca_trialMean(:,cueBlockEarly),2));
        meanCueLate = squeeze(mean(trialAlignedData.(trialTypes{i}).Ca_trialMean(:,cueBlockLate),2)); 
        meanDelayEarly = squeeze(mean(trialAlignedData.(trialTypes{i}).Ca_trialMean(:,delayBlockEarly),2));
        meanDelayTurn = squeeze(mean(trialAlignedData.(trialTypes{i}).Ca_trialMean(:,delayTurnBlock),2));
        meanTurnITI = squeeze(mean(trialAlignedData.(trialTypes{i}).Ca_trialMean(:,turnITIblock),2));

        normSpRate = meanCueEarly + meanCueLate + meanDelayEarly + meanDelayTurn + meanTurnITI;
        
        trialAlignedData.(trialTypes{i}).t_selectIdx = [meanCueEarly meanCueLate meanDelayEarly meanDelayTurn meanTurnITI]./normSpRate;
    end
end

% sanity check by plotting very selective traces
figure;
hold on;
colors = jet(5);
for tIdx = 5:5
    disp(sum(trialAlignedData.wL_trials.t_selectIdx(:,tIdx)>.5))
    %plot(trialAlignedData.wL_trials.Ca(trialAlignedData.wL_trials.t_selectIdx(:,tIdx)>.6,thisTrial)','color',colors(tIdx,:));
    %plot(thisTrial,trialAlignedData.wL_trials.Ca(trialAlignedData.wL_trials.t_selectIdx(:,tIdx)>.5,thisTrial));
    imagesc(trialAlignedData.wL_trials.Ca(trialAlignedData.wL_trials.t_selectIdx(:,tIdx)>.5,thisTrial));
end
%}


