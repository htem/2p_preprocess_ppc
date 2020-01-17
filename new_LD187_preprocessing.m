% ATK 200115 new_LD187_preprocessing for use with suite2p
% Needs constrained_foopsi https://github.com/epnev/constrained-foopsi
% Needs HarveyLab helper functions (abfload) https://github.com/HarveyLab/helperFunctions
% Needs matlab signal processing toolbox

%% make sessionList_all variable for each mouse and save 
masterPath = '/home/atk13/ppc/2P_data/';
mouse = 'LD187';

sessionList = dir('/home/atk13/NeuroShare/Lee_Lab/Aaron/data/scanimage/LD187/LD187*');
sessionList_all = cell(size(sessionList,1),1);
for s = 1:size(sessionList,1)
    sessionList_all{s,1} = sessionList(s,1).name;
end

%% Select Session

%%%%%%%%%
session = sessionList_all{36,1};
disp(session)
trialAlignedData = struct;
%%%%%%%%%
%% line up 
lineUpSession(session) % for now run 141213
% Exports VirminCombined to ppc/2P_data/code_workspace/LD187/virmen/FILE.mat

%% Load suite2p and aligned virmen data

% Load suite2p output
suite2pPath = fullfile(masterPath,'scanimage',mouse,session,'suite2p');
suite2pOutput = fullfile(suite2pPath,'combined/Fall.mat');
s2p = load(suite2pOutput);
numCells = size(s2p.F,1);

% Load aligned virmen data
linedUpPath = fullfile(masterPath,'code_workspace',mouse,'virmen',[session '.mat']);
vData = load(linedUpPath);

%% Deconvolution
neucoeff = 0.7;
s2p.Fsub = s2p.F - neucoeff * s2p.Fneu;

% for now use suite2p's deconvolution
spks = s2p.spks;
% dF = zscore(s2p.Fsub,0,2);
%[c, s, options] = deconvolveCa(s2p.Fsub(1,:));
%[c,b,c1,g,sn,sp] = constrained_foopsi(s2p.Fsub(1:2,:));

%% Plot example traces
figure; hold on;
frames = 1:length(spks);
dt = 1/5.3; % 5.3 Hz sampling rate
time = frames * dt;
for i = 1:30%size(mySignal,1)
    plot(time,spks(i,:)+100*i);
end
xlabel('time (sec)');
ylabel('spike rate');


%% Parse virmen trials

trialAlignedData = parseVirmenTrials(vData.VirmenCombined, spks);
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
beh_idx = 4;

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
r = mean(trialAlignedData.bR_trials.Ca_trialMean(:,thisTrial),2);
l = mean(trialAlignedData.wL_trials.Ca_trialMean(:,thisTrial),2);
trialAlignedData.RL_selectIdx = (r - l)./(r + l);
figure; histogram(trialAlignedData.RL_selectIdx)
% definition fo cellSelectIdx from Harvey 2012 

% sanity-check activity for very selective cells
figure; hold on; 
subplot(1,2,1);
plot(thisTrial,trialAlignedData.bR_trials.Ca_trialMean(trialAlignedData.RL_selectIdx>0.45,13:76)')
ylim([0,25]);
subplot(1,2,2); 
plot(thisTrial,trialAlignedData.wL_trials.Ca_trialMean(trialAlignedData.RL_selectIdx>0.45,13:76)')
ylim([0,25]);



%% Calculate pairwise correlation matrices

trialAlignedData.corr_all = corrcoef(spks');
for i = 1:length(trialTypes)
    if size(trialAlignedData.(trialTypes{i}).virmen,2) > 0
        trialAlignedData.(trialTypes{i}).corrcoef = corrcoef(trialAlignedData.(trialTypes{i}).Ca_trialMean');
        trialAlignedData.(trialTypes{i}).corrcoef = corrcoef(trialAlignedData.(trialTypes{i}).Ca_trialMean');
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
        tMax = nan(numCells);
        for cid = 1:numCells
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


