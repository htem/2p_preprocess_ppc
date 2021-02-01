%% ATK 201030 for loading from HD
% Adds ROC AUC metrics (from Najafi et al. 2019)
% updated 210201 to also include significance testing from
% older RL selectivity measure from Harvey et al. 2012
masterPath = "/Volumes/Aaron's 5TB HD/ppc/2P_data/";
mouse = 'LD187';


 %mySessions = {'LD187_141210'};
   
mySessions = {'LD187_141216','LD187_141215','LD187_141214','LD187_141213','LD187_141212',...
             'LD187_141211','LD187_141209','LD187_141208','LD187_141207','LD187_141206'};        
         
for i = 1:length(mySessions)
    session = cell2mat(mySessions(i));

    %% Select Session and lineup

    %%%%%%%%%
    %session = sessionList_all{36};
    %session = 'LD187_141216';
    disp(session)
    %%%%%%%%%
    % Load aligned virmen data
    linedUpPath = fullfile(masterPath,'code_workspace',mouse,'virmen',[session '.mat']);
    vData = load(linedUpPath);
    disp(['Virmen Data has ' num2str(size(vData.VirmenCombined,2)) ' frames']);

    % Load deconv data

    disp(['Loading deconv data for ' session]);
    %disp(session)
    output_dir_deconv = fullfile(masterPath, 'code_workspace',mouse,'deconvData');   
    load(fullfile(output_dir_deconv,session),'c','b','c1','g','sn','sp');

    % Calculate SNR metrics
    disp(['Calculating SNR for ' session]);
    snr = nan([],length(sn));
    for i = 1:length(sn)
        snr(i) = std(c(i,:))^2/sn(i)^2;
    end

    % Parse virmen trials
    disp(['Saving trial aligned data for ' session]);
    % Hack for 141210 below 
    % trialAlignedData = parseVirmenTrials(vData.VirmenCombined(:,6001:16000), sp);
    trialAlignedData = parseVirmenTrials(vData.VirmenCombined, sp);
    output_dir = fullfile(masterPath, 'code_workspace',mouse,'syncedData');   
    if ~exist(output_dir,'dir')
        mkdir(output_dir)
    end
    trialAlignedData.SNR_raw = snr;

    %% testing new R/L metrics
    lr_trials = ismember(trialAlignedData.trialType,[2,3]);
    trial_types = trialAlignedData.trialType(lr_trials);
    num_trials = length(trial_types);
    fullTrial_Act = nanmean(trialAlignedData.CaData(:,lr_trials,14:76),3);
    %
    % define blocks
    num_blocks = 5;

    blocks_Act = NaN([size(fullTrial_Act),5]);
    blocks_Act(:,:,1) =  nanmean(trialAlignedData.CaData(:,lr_trials,14:26),3); % 14 is running onset + 12 frames after
    blocks_Act(:,:,2) =  nanmean(trialAlignedData.CaData(:,lr_trials,27:38),3); % 12 frames before cue offset (frame 39)
    blocks_Act(:,:,3)=  nanmean(trialAlignedData.CaData(:,lr_trials,39:51),3); % 39 is cue offset + 12 frames after
    blocks_Act(:,:,4)= nanmean(trialAlignedData.CaData(:,lr_trials,52:64),3); % 12 frames before end of trial (turn a certain amt)
    blocks_Act(:,:,5) = nanmean(trialAlignedData.CaData(:,lr_trials,65:76),3); % Trial end (reward given) an 12 frames after (dark, ITI)
    blocks_labels = {'cueEarly_Act','cueLate_Act','delayEarly_Act','delayTurn_Act','turnITI_Act'};

    % create shuffle distribution of trial types
    num_shuf = 20;
    trial_types_shuf = NaN(num_shuf, num_trials);
    for shuf_idx = 1:num_shuf
        trial_types_shuf(shuf_idx,:) = trial_types(randperm(length(trial_types)));
    end
    
    % Calculate ROC AUC metrics
    tic
    num_masks = size(fullTrial_Act,1);

    AUC_fullTrial = NaN(num_masks,1);
    AUC_fullTrial_shuf = NaN(num_masks, num_blocks, num_shuf);
    selectivity_AUC = NaN(num_masks,1);
    AUC_blocks = NaN(num_masks,num_blocks);
    AUC_blocks_shuf = NaN(num_masks,num_blocks,num_shuf);
    selectivity_AUC_blocks = NaN(num_blocks,num_blocks);

    for mask_idx = 1:size(fullTrial_Act,1)
        [X,Y,T,AUC_fullTrial(mask_idx)] = perfcurve(trial_types,fullTrial_Act(mask_idx,:),2);
        for blk_idx = 1:5
            [X,Y,T,AUC_blocks(mask_idx,blk_idx)] = perfcurve(trial_types,blocks_Act(mask_idx,:,blk_idx),2);
        end

        for shuf_idx = 1:num_shuf
            [X,Y,T,AUC_fullTrial_shuf(mask_idx,shuf_idx)] = perfcurve(trial_types_shuf(shuf_idx,:),fullTrial_Act(mask_idx,:),2);

            for blk_idx = 1:num_blocks
                [X,Y,T,AUC_blocks_shuf(mask_idx,blk_idx,shuf_idx)] = perfcurve(trial_types_shuf(shuf_idx,:),blocks_Act(mask_idx,:,blk_idx),2);
            end
        end
        disp(['AUC for mask ' num2str(mask_idx) ' is ' num2str(AUC_fullTrial(mask_idx)) ' vs shuffle: ' num2str(nanmean(AUC_fullTrial_shuf(mask_idx,:)))]);
        if prctile(AUC_fullTrial_shuf(mask_idx,:),2.5)>AUC_fullTrial(mask_idx)
            selectivity_AUC(mask_idx) = 2;
        elseif prctile(AUC_fullTrial_shuf(mask_idx,:),97.5)<AUC_fullTrial(mask_idx)    
            selectivity_AUC(mask_idx) = 3;
        else
            selectivity_AUC(mask_idx) = 0;
        end
        for blk_idx = 1:num_blocks
            if prctile(AUC_blocks_shuf(mask_idx,blk_idx,:),2.5)>AUC_blocks(mask_idx,blk_idx)
                selectivity_AUC_blocks(mask_idx,blk_idx) = 2;
            elseif prctile(AUC_blocks_shuf(mask_idx,blk_idx,:),97.5)<AUC_blocks(mask_idx,blk_idx)
                selectivity_AUC_blocks(mask_idx,blk_idx) = 3;
            else
                selectivity_AUC_blocks(mask_idx,blk_idx) = 0;
            end
        end         
    end
    toc
    
    % Calculate old RL metrics
    tic
    num_masks = size(fullTrial_Act,1);

    RL_fullTrial = NaN(num_masks,1);
    RL_fullTrial_shuf = NaN(num_masks, num_blocks, num_shuf);
    selectivity_RL = NaN(num_masks,1);
    RL_blocks = NaN(num_masks,num_blocks);
    RL_blocks_shuf = NaN(num_masks,num_blocks, num_shuf);
    selectivity_RL_blocks = NaN(num_blocks,num_blocks);

    for mask_idx = 1:size(fullTrial_Act,1)
        R_tr_mean = nanmean(fullTrial_Act(mask_idx,trial_types==2));
        L_tr_mean = nanmean(fullTrial_Act(mask_idx,trial_types==3));
        RL_fullTrial(mask_idx) = (R_tr_mean - L_tr_mean)/(R_tr_mean + L_tr_mean);
        for blk_idx = 1:5
            R_tr_mean_blocks = nanmean(blocks_Act(mask_idx,trial_types==2,blk_idx));
            L_tr_mean_blocks = nanmean(blocks_Act(mask_idx,trial_types==3,blk_idx));
            RL_blocks(mask_idx,blk_idx) = (R_tr_mean_blocks - L_tr_mean_blocks)/(R_tr_mean_blocks + L_tr_mean_blocks);
        end

        for shuf_idx = 1:num_shuf
            R_tr_mean = nanmean(fullTrial_Act(mask_idx,trial_types_shuf(shuf_idx,:)==2));
            L_tr_mean = nanmean(fullTrial_Act(mask_idx,trial_types_shuf(shuf_idx,:)==3));
            RL_fullTrial_shuf(mask_idx,shuf_idx) = (R_tr_mean - L_tr_mean)/(R_tr_mean + L_tr_mean);
            for blk_idx = 1:num_blocks
                R_tr_mean_blocks = nanmean(blocks_Act(trial_types_shuf(shuf_idx,:)==2,blk_idx));
                L_tr_mean_blocks = nanmean(blocks_Act(trial_types_shuf(shuf_idx,:)==3,blk_idx));
                RL_blocks_shuf(mask_idx,blk_idx) = (R_tr_mean_blocks - L_tr_mean_blocks)/(R_tr_mean_blocks + L_tr_mean_blocks);
            end
        end
        disp(['RL for mask ' num2str(mask_idx) ' is ' num2str(RL_fullTrial(mask_idx)) ' vs shuffle: ' num2str(nanmean(RL_fullTrial_shuf(mask_idx,:)))]);
        if prctile(RL_fullTrial_shuf(mask_idx,:),2.5)>RL_fullTrial(mask_idx)
            selectivity_RL(mask_idx) = 2;
        elseif prctile(AUC_fullTrial_shuf(mask_idx,:),97.5)<AUC_fullTrial(mask_idx)    
            selectivity_RL(mask_idx) = 3;
        else
            selectivity_RL(mask_idx) = 0;
        end
        for blk_idx = 1:num_blocks
            if prctile(RL_blocks_shuf(mask_idx,blk_idx,:),2.5)>RL_blocks(mask_idx,blk_idx)
                selectivity_RL_blocks(mask_idx,blk_idx) = 2;
            elseif prctile(RL_blocks_shuf(mask_idx,blk_idx,:),97.5)<RL_blocks(mask_idx,blk_idx)
                selectivity_RL_blocks(mask_idx,blk_idx) = 3;
            else
                selectivity_RL_blocks(mask_idx,blk_idx) = 0;
            end
        end         
    end
    toc
    
    output_dir = '/Users/akuan/repos/ppc_project_analysis/new_2p_analysis/210201_update/';
    save(fullfile(output_dir,session),'AUC_fullTrial', 'AUC_blocks', 'RL_fullTrial', 'RL_blocks',...
        'selectivity_AUC','selectivity_AUC_blocks','selectivity_RL_blocks','selectivity_RL_blocks','-v7.3');

    figure; plot(sum(selectivity_AUC_blocks == 2),'o-'); hold on;
    plot(sum(selectivity_AUC_blocks == 3),'o-');
    plot(sum(selectivity_AUC_blocks == 0),'o-');
    legend({'left','right','non selective'});
end
%%
figure;
plot(trialAlignedData.RL_selectIdx, AUC_fullTrial,'.');
xlabel('RL select Idx');
ylabel('AUC full Trial');

figure;
plot(trialAlignedData.RL_selectIdx, RL_fullTrial,'.');
xlabel('RL select Idx');
ylabel('RL full Trial');
%%
