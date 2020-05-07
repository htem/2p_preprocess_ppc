% ATK 200115 new_LD187_preprocessing for use with suite2p
% Needs constrained_foopsi https://github.com/epnev/constrained-foopsi
% Needs HarveyLab helper functions (abfload) https://github.com/HarveyLab/helperFunctions
% Needs matlab signal processing toolbox

%% make sessionList_all variable for each mouse and save 
masterPath = '/home/atk13/ppc/2P_data/';
mouse = 'LD187';

sessionList = dir('/n/groups/htem/temcagt/datasets/ppc/2P_data/scanimage/LD187/LD187*');
sessionList_all = cell(size(sessionList,1),1);
for s = 1:size(sessionList,1)
    sessionList_all{s,1} = sessionList(s,1).name;
end

%% Select Session and lineup

%%%%%%%%%
session = sessionList_all{38,1};
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
neucoeff = 0.7;
s2p.Fsub = s2p.F - neucoeff * s2p.Fneu;

%%
tic
[c,b,c1,g,sn,sp] = run_constrained_foopsi(double(s2p.Fsub));
toc

output_dir_deconv = fullfile(masterPath, 'code_workspace',mouse,'deconvData');   
if ~exist(output_dir_deconv,'dir')
    mkdir(output_dir_deconv)
end
save(fullfile(output_dir_deconv,session),'c','b','c1','g','sn','sp');

%% Load deconv data

output_dir_deconv = fullfile(masterPath, 'code_workspace',mouse,'deconvData');   
load(fullfile(output_dir_deconv,session),'c','b','c1','g','sn','sp');

%% Calculate SNR metrics
snr = nan([],length(sn));
for i = 1:length(sn)
    snr(i) = std(c(i,:))^2/sn(i)^2;
end

%% Plot SNR 
figure; plot(snr); xlabel('mask ID'); ylabel('SNR');
figure; histogram(snr, 100); xlabel('SNR'); ylabel('# of masks');  xlim([0,150]);
%%
cid =674;
dF = s2p.Fsub;
sorted_dF = sort(dF(cid,:));
baseline_adj = mode(round(100*sorted_dF(.05*size(sorted_dF,2):...
        .95*size(sorted_dF,2)))/100);
dF_zeroed = dF(cid,:) - baseline_adj;

figure; hold on; 
plot(dF_zeroed);
%plot(c(cid,:));
title(['SNR: '  num2str(snr(cid))])

%% Parse virmen trials
% Hack for 141210 %% trialAlignedData = parseVirmenTrials(vData.VirmenCombined(:,6001:16000), sp);
trialAlignedData = parseVirmenTrials(vData.VirmenCombined, sp);
output_dir = fullfile(masterPath, 'code_workspace',mouse,'syncedData');   
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
trialAlignedData.SNR_raw = snr;
save(fullfile(output_dir,session),'trialAlignedData');

%% Plot example traces
figure; hold on;
frames = 1:length(sp);
dt = 1/5.3; % 5.3 Hz sampling rate
time = frames * dt;
for i = 1:30%size(mySignal,1)
    plot(time,sp(i,:)+100*i);
    
end
xlabel('time (sec)');
ylabel('spike rate');

%% Plot example activity as test

s2p_cid = 8;
matlab_cid = s2p_cid+1;

% Define colors
cueEarlyColor = [0 161 75]/255; 
cueLateColor = [255 222 23]/255;
delayEarlyColor = [237 28 36]/255;
delayTurnColor = [127 63 152]/255;
turnITIcolor = [33 64 154]/255;
% Define blocks
cueBlockEarly = 14:26; % 14 is running onset + 12 frames after
cueBlockLate = 27:38; % 12 frames before cue offset (frame 39)
delayBlockEarly = 39:51; % 39 is cue offset + 12 frames after
delayTurnBlock = 52:64; % 12 frames before end of trial (turn a certain amt)
turnITIblock = 65:76; % Trial end (reward given) an 12 frames after (dark, ITI)

figure;
tiledlayout(1,2)

ax1 = nexttile;
%subplot(1,2,1); 
hold on;
plot(cueBlockEarly,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,cueBlockEarly),'color',cueEarlyColor);
plot(cueBlockLate,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,cueBlockLate),'color',cueLateColor);
plot(delayBlockEarly,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,delayBlockEarly),'color',delayEarlyColor);
plot(delayTurnBlock,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,delayTurnBlock),'color',delayTurnColor);
plot(turnITIblock,trialAlignedData.bR_trials.Ca_trialMean(matlab_cid,turnITIblock),'color',turnITIcolor);
title('right');

ax2 = nexttile;
%subplot(1,2,2); 
hold on;
plot(cueBlockLate,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,cueBlockLate),'color',cueLateColor);
plot(delayBlockEarly,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,delayBlockEarly),'color',delayEarlyColor);
plot(delayTurnBlock,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,delayTurnBlock),'color',delayTurnColor);
plot(turnITIblock,trialAlignedData.wL_trials.Ca_trialMean(matlab_cid,turnITIblock),'color',turnITIcolor);
title('left');
linkaxes([ax1 ax2],'xy')

%% compare suite2p and constrained_foopsi deconv
cell_idx = find(s2p.iscell(:,1)==1);

figure;
num_ex = 3;
start = 100;
num_pts = 400;
tsec = (1:num_pts)/5.6;
for i = start:start+num_ex-1
    disp(cell_idx(i))
    subplot(num_ex,3,3*(i-start)+1);
    plot(tsec,s2p.Fsub(cell_idx(i),1:num_pts));
    title('raw');
    subplot(num_ex,3,3*(i-start) + 2);
    plot(tsec,s2p.spks(cell_idx(i),1:num_pts));
    title('suite2p');
    subplot(num_ex,3,3*(i-start)+3);    
    plot(tsec,sp(cell_idx(i),1:num_pts));
    title('const foopsi');
end