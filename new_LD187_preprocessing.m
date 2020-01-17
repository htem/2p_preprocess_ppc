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

%% Select Session and lineup

%%%%%%%%%
session = sessionList_all{33,1};
disp(session)
%%%%%%%%%

%% line up 
lineUpSession(session) 
% Exports VirminCombined to ppc/2P_data/code_workspace/LD187/virmen/FILE.mat

%% Load suite2p and aligned virmen data

% Load suite2p output
suite2pPath = fullfile(masterPath,'scanimage',mouse,session,'suite2p');
suite2pOutput = fullfile(suite2pPath,'combined/Fall.mat');
s2p = load(suite2pOutput);
numCells = size(s2p.F,1);
disp(['Suite2p video data has ' num2str(size(s2p.F,2)) ' frames']);

% Load aligned virmen data
linedUpPath = fullfile(masterPath,'code_workspace',mouse,'virmen',[session '.mat']);
vData = load(linedUpPath);
disp(['Virmen Data has ' num2str(size(vData.VirmenCombined,2)) ' frames']);
%% Deconvolution
neucoeff = 0.7;
s2p.Fsub = s2p.F - neucoeff * s2p.Fneu;

% for now use suite2p's deconvolution
spks = s2p.spks;
% dF = zscore(s2p.Fsub,0,2);
%[c, s, options] = deconvolveCa(s2p.Fsub(1,:));
%[c,b,c1,g,sn,sp] = constrained_foopsi(s2p.Fsub(1:2,:));


%% Parse virmen trials

trialAlignedData = parseVirmenTrials(vData.VirmenCombined, spks);
output_dir = fullfile(masterPath, 'code_workspace',mouse,'syncedData');   
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
save(fullfile(output_dir,session),'trialAlignedData');

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


%% Plot example activity as test

s2p_cid = 1351;
matlab_cid = s2p_cid+1;

% Define colors
cueEarlyColor = [0 161 75]/255; 
cueLateColor = [255 222 23]/255;
delayEarlyColor = [237 28 36]/255;
delayTurnColor = [127 63 152]/255;
turnITIcolor = [33 64 154]/255;

figure;

subplot(1,2,1); hold on;
plot(cueBlockEarly,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,cueBlockEarly),'color',cueEarlyColor);
plot(cueBlockLate,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,cueBlockLate),'color',cueLateColor);
plot(delayBlockEarly,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,delayBlockEarly),'color',delayEarlyColor);
plot(delayTurnBlock,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,delayTurnBlock),'color',delayTurnColor);
plot(turnITIblock,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,turnITIblock),'color',turnITIcolor);
title('right');

subplot(1,2,2); hold on;
plot(cueBlockEarly,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,cueBlockEarly),'color',cueEarlyColor);
plot(cueBlockLate,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,cueBlockLate),'color',cueLateColor);
plot(delayBlockEarly,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,delayBlockEarly),'color',delayEarlyColor);
plot(delayTurnBlock,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,delayTurnBlock),'color',delayTurnColor);
plot(turnITIblock,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,turnITIblock),'color',turnITIcolor);
title('left');


