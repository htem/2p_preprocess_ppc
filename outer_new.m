% needs https://github.com/HarveyLab/helperFunctions (abfload)



%cd('Z:\HarveyLab\Aaron\scanimage\AK001')
cd('/home/atk13/NeuroShare/Lee_Lab/Aaron/data/scanimage/LD187')
mouseDir = pwd;
sessionList = dir('LD187*');
%
%{
for session = 1:size(sessionList,1)
cd(fullfile(mouseDir,sessionList(session).name))
acq = Acquisition2P([],@ak_init);
acq.motionCorrect;
output_dir = fullfile('Z:\HarveyLab\Aaron\scanimage\',sessionList(session).name(1:5),...
    'sessions',sessionList(session).name);
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
load([sessionList(session).name '_001']);
eval(['acq = ' sessionList(session).name '_001;'])
channelNum = 1;
acq.metaDataSI = metaDataSI;
sizeArray = repmat(sizeBit,[size(acq.Movies,2) 1]);
for sliceNum = 1:4
acq.correctedMovies.slice(sliceNum).channel(1).size = sizeArray;
acq.indexMovie(sliceNum,channelNum,output_dir);
acq.calcPxCov([],[],[],sliceNum,channelNum,output_dir);
save(fullfile(output_dir,sessionList(session).name),'acq','-v7.3')
end
clear acq
end
%}
%% make sessionList_all variable for each mouse and save 

sessionList_all = cell(size(sessionList,1),1);
for s = 1:size(sessionList,1)
    sessionList_all{s,1} = sessionList(s,1).name;
end