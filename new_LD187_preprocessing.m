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
session = sessionList_all{36,1};
%% line up 
lineUpSession(session) % for now run 141213
% Exports VirminCombined to ppc/2P_data/code_workspace/LD187/virmen/FILE.mat

%% Load suite2p and aligned virmen data

% Load suite2p output
suite2pPath = fullfile(masterPath,'scanimage',mouse,session,'suite2p');
suite2pOutput = fullfile(suite2pPath,'combined/Fall.mat');
s2p = load(suite2pOutput);

% Load aligned virmen data
linedUpPath = fullfile(masterPath,'code_workspace',mouse,'virmen',[session '.mat']);
vData = load(linedUpPath);

%% Deconvolution
neucoeff = 0.7;
s2p.Fsub = s2p.F - neucoeff * s2p.Fneu;

%[c, s, options] = deconvolveCa(s2p.Fsub(1,:));
%[c,b,c1,g,sn,sp] = constrained_foopsi(s2p.Fsub(1:2,:));

% for now skip deconvolution
dF = zscore(s2p.Fsub,0,2);

%% Plot example traces
figure; hold on;
frames = 1:length(dF);
dt = 1/5.3; % 5.3 Hz sampling rate
time = frames * dt;
for i = 1:30%size(mySignal,1)
    plot(time,dF(i,:)+3*i);
end
xlabel('time (sec)');
ylabel('dF/F');

%% Parse virmen trials

trialAlignedData = parseVirmenTrials(vData.VirmenCombined, dF);

% break into individual trial types
trialAlignedData.bR_trials_Ca = trialAlignedData.CaData(:,trialAlignedData.trialType==2,:);
trialAlignedData.wL_trials_Ca = trialAlignedData.CaData(:,trialAlignedData.trialType==3,:);
trialAlignedData.new_trials_Ca = trialAlignedData.CaData(:,trialAlignedData.trialType==4,:);
trialAlignedData.bR_trials_virmen = trialAlignedData.virmenData(:,trialAlignedData.trialType==2,:);
trialAlignedData.wL_trials_virmen = trialAlignedData.virmenData(:,trialAlignedData.trialType==3,:);
trialAlignedData.new_trials_virmen = trialAlignedData.virmenData(:,trialAlignedData.trialType==4,:);

%% Sanity check with trialwise behavioral data
beh_idx = 2;
figure; hold on;
for i = 1:size(trialAlignedData.bR_trials_virmen,2)
    plot(squeeze(trialAlignedData.bR_trials_virmen(beh_idx,i,:)),'b-')
end
for i = 1:size(trialAlignedData.wL_trials_virmen,2)
    plot(squeeze(trialAlignedData.wL_trials_virmen(beh_idx,i,:)),'r-')
end    
for i = 1:size(trialAlignedData.new_trials_virmen,2)
    plot(squeeze(trialAlignedData.new_trials_virmen(beh_idx,i,:)),'k-')
end    
%% Plot trial avg dF
figure; hold on;
frames = 1:size(trialAvg_dF,2);
time = frames;% * acq.metaDataSI.SI4.fastZPeriod;
for i = 1:size(dF,1)
    plot(time,trialAvg_dF(i,:));
end
xlabel('time (frames)');
ylabel('dF/F');

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

%% Plot trial avg dF
figure; hold on;
frames = 1:size(trialAvg_dF,2);
time = frames(12:64);% * acq.metaDataSI.SI4.fastZPeriod;
for i = 1:size(dF,1)
    plot(time,zTrialMean(cellOrder(i),12:64)+0.3*i);
end
xlabel('time (frames)');
ylabel('dF/F, zscored');

