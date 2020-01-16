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
trialAlignedData = struct;

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
    plot(time,spks(i,:)+3*i);
end
xlabel('time (sec)');
ylabel('dF/F');


%% Parse virmen trials

trialAlignedData = parseVirmenTrials(vData.VirmenCombined, spks);

% break into individual trial types
trialAlignedData.bR_trials_Ca = trialAlignedData.CaData(:,trialAlignedData.trialType==2,:);
trialAlignedData.wL_trials_Ca = trialAlignedData.CaData(:,trialAlignedData.trialType==3,:);
trialAlignedData.new_trials_Ca = trialAlignedData.CaData(:,trialAlignedData.trialType==4,:);
trialAlignedData.bR_trials_virmen = trialAlignedData.virmenData(:,trialAlignedData.trialType==2,:);
trialAlignedData.wL_trials_virmen = trialAlignedData.virmenData(:,trialAlignedData.trialType==3,:);
trialAlignedData.new_trials_virmen = trialAlignedData.virmenData(:,trialAlignedData.trialType==4,:);

%% Sanity check with trialwise behavioral data
beh_idx = 4;
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


%% Calculate trial average activities

trialAlignedData.bR_Ca_trialMean = squeeze(mean(trialAlignedData.bR_trials_Ca,2));
trialAlignedData.wL_Ca_trialMean = squeeze(mean(trialAlignedData.wL_trials_Ca,2));
trialAlignedData.new_Ca_trialMean = squeeze(mean(trialAlignedData.new_trials_Ca,2));

%% Calculate R/L selectivities
% exclude ITI-before (which is unrelated to this trial), 13:76
thisTrial = 13:76;
r = mean(trialAlignedData.bR_Ca_trialMean(:,thisTrial),2);
l = mean(trialAlignedData.wL_Ca_trialMean(:,thisTrial),2);
trialAlignedData.cellSelectIdx = (r - l)./(r + l);
figure; histogram(trialAlignedData.cellSelectIdx)
% definition fo cellSelectIdx from Harvey 2012 

% sanity-check activity for very selective cells
figure; hold on; 
subplot(1,2,1);
plot(thisTrial,trialAlignedData.bR_Ca_trialMean(trialAlignedData.cellSelectIdx>0.45,13:76)')
ylim([0,25]);
subplot(1,2,2); 
plot(thisTrial,trialAlignedData.wL_Ca_trialMean(trialAlignedData.cellSelectIdx>0.45,13:76)')
ylim([0,25]);

%% Calculate pairwise correlation matrices

trialAlignedData.corr_all = corrcoef(spks');
trialAlignedData.corr_bR_trials = corrcoef(trialAlignedData.bR_Ca_trialMean');
trialAlignedData.corr_wL_trials = corrcoef(trialAlignedData.wL_Ca_trialMean');

% sanity check by comparing
figure;
plot(trialAlignedData.corr_all(:),trialAlignedData.corr_bR_trials(:),'.');
xlabel('total correlation'); ylabel('bR trial correlation');

figure;
plot(trialAlignedData.corr_bR_trials(:),trialAlignedData.corr_wL_trials(:),'.');
xlabel('bR trials correlations'); ylabel('wL trials correlation');

%% Calculate time center-of-mass 

rmean = trialAlignedData.bR_Ca_trialMean(:,thisTrial);
lmean = trialAlignedData.wL_Ca_trialMean(:,thisTrial);
trialAlignedData.bR_tCOM = (rmean*thisTrial')./sum(rmean,2);
trialAlignedData.wL_tCOM = (lmean*thisTrial')./sum(lmean,2);

% sanity check by plotting high COM 
figure;
plot(trialAlignedData.bR_Ca_trialMean(trialAlignedData.bR_tCOM>48,:)');

% sanity check by RL correlation
figure; plot(trialAlignedData.bR_tCOM,trialAlignedData.wL_tCOM,'o');
xlabel('bR_tCOM');ylabel('wL_tCOM');

%% Calculate time maximum
for i = 1:size(rmean,1)
    bR_tMax(i) = find(rmean(i,:) == max(rmean(i,:)),1,'last');
    wL_tMax(i) = find(lmean(i,:) == max(lmean(i,:)),1,'last');
end
%wL_tMax = find(lmean == max(lmean,2));
%[m, wL_tMax] = max(lmean,2);

% sanity check by plotting high COM 
figure;
plot(trialAlignedData.bR_Ca_trialMean(bR_tMax>55,:)');

% sanity check by plotting
figure; plot(bR_tMax,wL_tMax,'.');
xlabel('bR_tCOM');ylabel('wL_tCOM');

%% Plot example activity as test

s2p_cid = 1351;
matlab_cid = s2p_cid+1;
figure;

subplot(1,2,1); hold on;
plot(cueBlockEarly,trialAlignedData.bR_Ca_trialMean(matlab_cid,cueBlockEarly),'color',cueEarlyColor);
plot(cueBlockLate,trialAlignedData.bR_Ca_trialMean(matlab_cid,cueBlockLate),'color',cueLateColor);
plot(delayBlockEarly,trialAlignedData.bR_Ca_trialMean(matlab_cid,delayBlockEarly),'color',delayEarlyColor);
plot(delayTurnBlock,trialAlignedData.bR_Ca_trialMean(matlab_cid,delayTurnBlock),'color',delayTurnColor);
plot(turnITIblock,trialAlignedData.bR_Ca_trialMean(matlab_cid,turnITIblock),'color',turnITIcolor);
title('right');

subplot(1,2,2); hold on;
plot(cueBlockEarly,trialAlignedData.wL_Ca_trialMean(matlab_cid,cueBlockEarly),'color',cueEarlyColor);
plot(cueBlockLate,trialAlignedData.wL_Ca_trialMean(matlab_cid,cueBlockLate),'color',cueLateColor);
plot(delayBlockEarly,trialAlignedData.wL_Ca_trialMean(matlab_cid,delayBlockEarly),'color',delayEarlyColor);
plot(delayTurnBlock,trialAlignedData.wL_Ca_trialMean(matlab_cid,delayTurnBlock),'color',delayTurnColor);
plot(turnITIblock,trialAlignedData.wL_Ca_trialMean(matlab_cid,turnITIblock),'color',turnITIcolor);
title('left');
%%
figure; hold on;
subplot(1,2,1);
plot(trialAlignedData.bR_Ca_trialMean(matlab_cid,thisTrial));
title('right');
subplot(1,2,2);
plot(trialAlignedData.wL_Ca_trialMean(matlab_cid,thisTrial));
title('left');

%%
plot(cueBlockEarly+toffset,spData_all(cIdx,cueBlockEarly,sIdx,k)+j*yoffset,'color',cueEarlyColor);
plot(cueBlockLate+2*toffset,spData_all(cIdx,cueBlockLate,sIdx,k)+j*yoffset,'color',cueLateColor);
plot(delayBlockEarly+2*toffset,spData_all(cIdx,delayBlockEarly,sIdx,k)+j*yoffset,'color',delayEarlyColor);
plot(delayTurnBlock+3*toffset,spData_all(cIdx,delayTurnBlock,sIdx,k)+j*yoffset,'color',delayTurnColor);
plot(turnITIblock+3*toffset,spData_all(cIdx,turnITIblock,sIdx,k)+j*yoffset,'color',turnITIcolor);
            
for k = 1:2
    subplot(1,2,k); hold on; axis off; % title(lr(k));
    %text(-12,i*yoffset, num2str(cIdx));
    line([runOnset+toffset runOnset+toffset],[0 yoffset]+yoffset*(earliestSession-1),'color','k');
    text(runOnset+toffset, yoffset * (earliestSession-1), 'Run','Rotation',90,'HorizontalAlignment','right','color','g');
    line([cueOff+2*toffset cueOff+2*toffset],[0 yoffset]+yoffset*(earliestSession-1),'color','k');
    text(cueOff+2*toffset, yoffset * (earliestSession-1), 'Cue Off','Rotation',90,'HorizontalAlignment','right','color','r');
    line([trialEnd+3*toffset trialEnd+3*toffset],[0 yoffset]+yoffset*(earliestSession-1),'color','k');
    text(trialEnd+3*toffset, yoffset * (earliestSession-1), 'End','Rotation',90,'HorizontalAlignment','right','color','b');
    for j = 1:length(mySessions)
        sIdx = mySessions(j);
        if ~isnan(spData_tAvg(cIdx,mySessions(j),k))
            if k == 2
                text(-10,j*yoffset, ['Session ' num2str(mySessions(j))],'HorizontalAlignment','center');
            end
            %timingColor = [cueCoef(cIdx,j,k) delayCoef(cIdx,j,k) turnCoef(cIdx,j,k)].^cpower;
            timingColor = 'b';
            %plot(ITIbeforeBlock,spData_all(cIdx,ITIbeforeBlock,sIdx,k)+j*yoffset,'color','k');

            % Display selectivity indices (colored text)
            s = cellSelectIdx(cIdx,j);
            lr_color = lr_cmap(round(s*100)+101,:);
            subplot(1,2,2); text(t(end)+10*toffset,j*yoffset,num2str(s,2),...
                'color',lr_color,'BackgroundColor',[.5 .5 .5]);
            subplot(1,2,k); hold on; axis off;
        end
    end
    axis([0 max(t)+3*toffset (earliestSession-2)*yoffset (latestSession+3)*yoffset]);
end


%% Quantify signal-to-noise

% From raw (non-trial aligned) traces.
snr_ratio_raw = prctile(s2p.Fsub,90,2)./abs(prctile(s2p.Fsub,10,2));
figure; histogram(snr_ratio_raw);

%% Sanity check by plotting high and low snr 
figure;
plot(s2p.Fsub(snr_ratio_raw>1.3,:)');
%%
figure; hold on;
plot(trialAlignedData.bR_Ca_trialMean(1:30,:)')
%% Plot trial avg dF
figure; hold on;
frames = 1:size(trialAvg_dF,2);
time = frames;% * acq.metaDataSI.SI4.fastZPeriod;
for i = 1:size(spks,1)
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
for i = 1:size(spks,1)
    plot(time,zTrialMean(cellOrder(i),12:64)+0.3*i);
end
xlabel('time (frames)');
ylabel('dF/F, zscored');

