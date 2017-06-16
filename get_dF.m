cd('Z:\Laura\code_workspace\LD187\sessions')
mother_dir = pwd;
sessionList = dir('LD*');
for session = 32%1:size(sessionList,1)%3%15%37%%%%LD186_10%%%LD183_36
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