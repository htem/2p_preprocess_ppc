% needs https://github.com/HarveyLab/helperFunctions (abfload)



%cd('Z:\HarveyLab\Aaron\scanimage\AK001')
cd('/home/atk13/NeuroShare/Lee_Lab/Aaron/data/scanimage/LD187')
mouseDir = pwd;
sessionList = dir('LD187*');
%

%% make sessionList_all variable for each mouse and save 
sessionList_all = cell(size(sessionList,1),1);
for s = 1:size(sessionList,1)
    sessionList_all{s,1} = sessionList(s,1).name;
end


% line up 
lineUpAll_new(sessionList_all)