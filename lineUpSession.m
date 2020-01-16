function lineUpSession(session)
%lines up imaging frames with virmen. requires pclamp and virmen files to run. output is binned virmen data with length the same as the number of frames.
% 200115 ATK updated for use with suite2p

data_path = '/home/atk13/ppc/2P_data/';
mouse = session(1:5);
frameClockThresh = 4.5; % ATK in Volts

 
VirmenData1 = [];
VirmenData2 = [];
VirmenData3 = [];
VirmenData4 = [];
allRewardsStart = [];
allRewardsEnd = [];
VirmenDataMazeID = [];
maxRowsData = 11;
file = session(1:12);
acqSize = 2000;

pclamp_directory = fullfile(data_path, 'pclamp', mouse, file);
virmenList = dir(fullfile(data_path, 'virmen', mouse, [file '*.mat']));

%check for multiple virmen files
dataSize = 0;
virmenFilesUSE = [];
for fileNumber = 1:size(virmenList,1)
    if isempty(strfind(virmenList(fileNumber,1).name,'Cell'))
        virmenFilesUSE = [virmenFilesUSE fileNumber]; %#ok<AGROW>
        load(fullfile(data_path, 'virmen', mouse, ...
        virmenList(fileNumber,1).name),'data')
        dataSize = cat(2,dataSize,size(data,2));
    end
end

if size(dataSize,2)>2
    error('check number of virmen files')
end

% combine multiple virmen files as default but will stop with error
% warning first
numVirmenIts = sum(dataSize);
combinedVirmenFiles = zeros(maxRowsData,numVirmenIts);
mazeName = cell(1,numVirmenIts);
dataItInd = 1;
for fileNumber = virmenFilesUSE
    load(fullfile(data_path, 'virmen', mouse, ...
        virmenList(fileNumber,1).name),'data')
    combinedVirmenFiles(1:size(data,1),...
        dataItInd:dataItInd+size(data,2)-1) = data;
    dataItInd = dataItInd+size(data,2);
end
numVirmenIts = dataItInd-1;

%set up variable for frame associated with each variable
itCount = zeros(1,numVirmenIts);

d = dir(fullfile(pclamp_directory, [session(7:end) '*.abf']));
useFile = ones(size(d));
for f = 1:size(d,1)
    %remove orientation tuning files for V1 data
    useFile(f) = isempty(strfind(d(f).name,'vis'));
end
d = d(useFile==1);

clear last_It
virmenFileShiftedItInd = 1;
It_ind =1;
figure;
for trial = 1:size(d,1) % for each pclamp file

    %load pclamp files
    pclamp_file_string = fullfile(pclamp_directory, d(trial,:).name);
    [raw_traces,~,~] = abfload(pclamp_file_string);

    %visualize 
    subplot(3,ceil(size(d,1)/3),trial)
    hold on
    plot(raw_traces(:,2))   %frame clock
    plot(raw_traces(:,1))   %-1 on virmen its, -5 otherwise

    real_length = min(5900000,size(raw_traces,1));
    [lastInd] = find(raw_traces(1:real_length,2)>frameClockThresh,1,'last'); % ATK changed thresh to 4.5 V - wasn't reaching 5 
    raw_traces = raw_traces(1:lastInd+200,:);

    %crop file to end of imaging
    diff_raw_traces = diff(raw_traces(:,1));
    indPulse= find(raw_traces(:,1)>0);                          %find indices of virmen iteration pulses in pclamp
    maxIt = round(10*max(raw_traces(indPulse,1)))*1e4;          %exact frame number for relative counting 

    if maxIt < 1000000
        first_pulse_time = min(find(diff(raw_traces(:,1))<-2));
        total_virmen_its = size(find(diff(raw_traces(:,1))<-2),1);
        It_ind = maxIt + dataSize(virmenFileShiftedItInd); %set iteration index of largest virmen pulse
        %label all following pulses with iteration number
        for pClamp_ind = max(indPulse):size(raw_traces,1)-1
            if diff_raw_traces(pClamp_ind) < -2
                itCount(It_ind) = pClamp_ind+1;
                It_ind = It_ind + 1;
            end
        end
        last_It = It_ind;
        first_It = last_It - total_virmen_its;
        It_ind = first_It;
        %now starting from the front, label all pulses with iteration number
        for pClamp_ind = first_pulse_time:max(indPulse)
            if diff_raw_traces(pClamp_ind) < -2
                if diff_raw_traces(pClamp_ind) > -10
                    itCount(It_ind) = pClamp_ind+1;
                    It_ind = It_ind + 1;
                end
            end
        end
    else
        %if virmen starts during pclamp file
        It_ind =1+ dataSize(virmenFileShiftedItInd);
        itCount(It_ind:end) = 0;
        for pClamp_ind = min(indPulse):find(diff(raw_traces(:,1))<-2,1,'last')
            if diff_raw_traces(pClamp_ind) < -2
                itCount(It_ind) = pClamp_ind+1;
                It_ind = It_ind + 1;
            end
        end
    end
    pulses = raw_traces(:,1)>-1.9;
    pulseStart = find(diff(pulses)==1);
    pulseEnd = find(diff(pulses)==-1);
    if min(pulseEnd)<min(pulseStart)
        pulseEnd = pulseEnd(2:end);
    end
    if max(pulseEnd)<max(pulseStart)
        pulseStart = pulseStart(1:end-1);
    end
    if size(pulseEnd) ~= size(pulseStart)
        error('mismatch number of pulses')
    end
    pulseLength = pulseEnd - pulseStart;
    rewardIndStart = pulseStart(pulseLength>900);
    rewardIndEnd = pulseEnd(pulseLength>900);
    allRewardsStart = [allRewardsStart; rewardIndStart];%#ok<AGROW>
    allRewardsEnd = [allRewardsEnd; rewardIndEnd];%#ok<AGROW>

end
legend('imaging frame clock (scanimage)', 'behavior iteration clock (virmen/pclamp)');


figure;plot(itCount,'k.')
totalFrames = (acqSize)*size(d,1);
frame_times = zeros(totalFrames,2,size(d,1));

for aquisitionNumber = 1:size(d,1)
    pclamp_file_string = fullfile(pclamp_directory, d(aquisitionNumber,:).name);
    if size(d,1)>1
    [raw_traces,samplinginterval,headerinfo] = abfload(pclamp_file_string);

    %crop file to end of imaging
    real_length = min(5900000,size(raw_traces,1));
    [lastInd] = find(raw_traces(1:real_length,2)>frameClockThresh,1,'last');
    raw_traces = raw_traces(1:lastInd+2%% Quick parsing of virmen output
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

syncVirmenData = vData.VirmenCombined;

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
% frame 13 is trial onset with 12 frames before  1-12,13: 1-12, no cue, no image? Fr 13 cue goes on?
% frame 14 is running onset with 12 frames after 14,15-26
% frame 39 is cue offset (delay period onset) with 12 frames before and after 27-38,39,40-51: Fr 39 is when the mouse gets far enough to turn off the cue?
% frame 54 is trial end with 12 frames before and after Do you mean 64? 52-63,64,65-76: Fr 54 is end of trial. Is this when the mouse turns? Is the reward given here?

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
figure; imagesc(trialAvg_dF); colormap gray;0,:);

    end
    y = diff(round(raw_traces(:,2)/5));
    [~,frameStart] = findpeaks(y,'Threshold',0.9);
    [~,frameStop] = findpeaks(-y,'Threshold',0.9);
    frame_times(1:size(frameStart,1),1,aquisitionNumber) = frameStart;
    frame_times(1:size(frameStop,1),2,aquisitionNumber) = frameStop;
    if aquisitionNumber==1
        numFrames = size(frameStart,1);
    end
end
peaksItCount = [itCount(1:end-1) 0];
[pclamp_end,it_end] = findpeaks(peaksItCount,'minpeakheight',1000000);

if size(pclamp_end) ~= size(d,1)
    error('number of itCount peaks and pClamp files arent the same')
end

virmenFrameNumber = zeros(1,size(combinedVirmenFiles,2));
for aquisitionNumber = 1:size(d,1)
    if aquisitionNumber>1
        it_framestart=it_end(aquisitionNumber-1)+1;
    else
        it_framestart=1;
    end
    for frameNumber = 1:numFrames
        iterationNumber = find(itCount(it_framestart:it_end(aquisitionNumber)) > frame_times(frameNumber,1,aquisitionNumber)...
            & itCount(it_framestart:it_end(aquisitionNumber))  < frame_times(frameNumber,2,aquisitionNumber));
        virmenFrameNumber(iterationNumber+it_framestart-1) = (aquisitionNumber-1)*numFrames+ frameNumber;
    end
end

%figure;  plot (virmenFrameNumber,'linewidth',3)
figure;  plot (virmenFrameNumber,'k.')
xlabel('virmen iteration')
ylabel('scanimage frame number')

[ outputSlice1, outputSlice2, outputSlice3, outputSlice4, outputMazeName ] =...
    saveVirmenData(combinedVirmenFiles, mazeName,virmenFrameNumber);
VirmenDataID(1,size(VirmenData1,2)+1:size(VirmenData1,2)+size(outputSlice1,2)) = {file};
VirmenData1 = cat(2,VirmenData1,outputSlice1(1:10,:));
VirmenData2 = cat(2,VirmenData2,outputSlice2(1:10,:));
VirmenData3 = cat(2,VirmenData3,outputSlice3(1:10,:));
VirmenData4 = cat(2,VirmenData4,outputSlice4(1:10,:));
VirmenDataMazeID = cat(2,VirmenDataMazeID,outputMazeName(1:size(outputSlice1,2)));

combineVirmen
output_dir = fullfile(data_path, 'code_workspace',mouse,'virmen');   
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
save(fullfile(output_dir,file),'VirmenCombined');

