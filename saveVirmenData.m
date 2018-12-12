function [syncVirmenData, outputDATA, outputSlice1, outputSlice2, ...
    outputSlice3, outputSlice4, outputSlice5,outputSlice6,outputMazeName ] = ...
    saveVirmenData(combinedVirmenFiles, mazeName,virmenFrameNumber)

numSlices = 6;
outputDATA = zeros(size(combinedVirmenFiles,1), max(virmenFrameNumber));
outputMazeName = cell(1,max(virmenFrameNumber));
    for frameNumber = (min(virmenFrameNumber(virmenFrameNumber~=0)))+1:max(virmenFrameNumber)
        if mod(frameNumber-1,5000) ~= 0 && ~isempty(find(virmenFrameNumber==frameNumber, 1))
            itsUSE = find(virmenFrameNumber==frameNumber);
            outputDATA(1,frameNumber) = mean(combinedVirmenFiles(1,itsUSE));
            outputDATA(2,frameNumber) = mean(combinedVirmenFiles(2,itsUSE));
            outputDATA(3,frameNumber) = mean(combinedVirmenFiles(3,itsUSE));
            outputDATA(4,frameNumber) = mean(combinedVirmenFiles(4,itsUSE));
            outputDATA(5,frameNumber) = mean(combinedVirmenFiles(5,itsUSE));
            outputDATA(6,frameNumber) = mean(combinedVirmenFiles(6,itsUSE));
            outputDATA(7,frameNumber) = combinedVirmenFiles(7,itsUSE(end));
            outputDATA(8,frameNumber) = max(combinedVirmenFiles(8,itsUSE));
            outputDATA(9,frameNumber) = max(combinedVirmenFiles(9,itsUSE));
            outputDATA(10,frameNumber) = mode(combinedVirmenFiles(10,itsUSE));
            outputMazeName(frameNumber) = {mazeName(itsUSE(1))};
        end
    end

empties = find(outputDATA(1,:)==0);
for x = (min(virmenFrameNumber(virmenFrameNumber~=0)))+1:size(empties,2)
    outputDATA(:,empties(x))= outputDATA(:,empties(x)-1);
end
length = size(outputDATA,2);
outputSlice1 = zeros(size(outputDATA,1),size(1:numSlices:length - mod(length, numSlices),2));
outputSlice2 = zeros(size(outputDATA,1),size(2:numSlices:length - mod(length, numSlices),2));
outputSlice3 = zeros(size(outputDATA,1),size(3:numSlices:length - mod(length, numSlices),2));
outputSlice4 = zeros(size(outputDATA,1),size(4:numSlices:length - mod(length, numSlices),2));
outputSlice5 = zeros(size(outputDATA,1),size(5:numSlices:length - mod(length, numSlices),2));
outputSlice6 = zeros(size(outputDATA,1),size(6:numSlices:length - mod(length, numSlices),2));
syncVirmenData = zeros(size(outputDATA,1),size(4:numSlices:length - mod(length, numSlices),2));

outputSlice1(:,1:floor(length/numSlices)) = outputDATA(:,1:numSlices:length - mod(length, numSlices));
outputSlice2(:,1:floor(length/numSlices)) = outputDATA(:,2:numSlices:length - mod(length, numSlices));
outputSlice3(:,1:floor(length/numSlices)) = outputDATA(:,3:numSlices:length - mod(length, numSlices));
outputSlice4(:,1:floor(length/numSlices)) = outputDATA(:,4:numSlices:length - mod(length, numSlices));
outputSlice5(:,1:floor(length/numSlices)) = outputDATA(:,5:numSlices:length - mod(length, numSlices));
outputSlice6(:,1:floor(length/numSlices)) = outputDATA(:,6:numSlices:length - mod(length, numSlices));
outputMazeName(1:floor(length/numSlices)) = outputMazeName(1:numSlices:length - mod(length, numSlices));

% syncVirmenData is volumewise average (approx 2 fps on Loki)
syncVirmenData = mean(cat(3, outputSlice1, outputSlice2, outputSlice3, outputSlice4),3);

end