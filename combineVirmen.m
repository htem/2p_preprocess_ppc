% VirmenCombined = nan(10,size(VirmenData1,2));
VirmenCombined = nan(10,size(VirmenData1,2));
for frameNumber = 1:size(VirmenData1,2)
VirmenCombined(1,frameNumber) = mean([VirmenData1(1,frameNumber) VirmenData2(1,frameNumber) VirmenData3(1,frameNumber) VirmenData4(1,frameNumber)]);
VirmenCombined(2,frameNumber) = mean([VirmenData1(2,frameNumber) VirmenData2(2,frameNumber) VirmenData3(2,frameNumber) VirmenData4(2,frameNumber)]);
VirmenCombined(3,frameNumber) = mean([VirmenData1(3,frameNumber) VirmenData2(3,frameNumber) VirmenData3(3,frameNumber) VirmenData4(3,frameNumber)]);
VirmenCombined(4,frameNumber) = mean([VirmenData1(4,frameNumber) VirmenData2(4,frameNumber) VirmenData3(4,frameNumber) VirmenData4(4,frameNumber)]);
VirmenCombined(5,frameNumber) = mean([VirmenData1(5,frameNumber) VirmenData2(5,frameNumber) VirmenData3(5,frameNumber) VirmenData4(5,frameNumber)]);
VirmenCombined(6,frameNumber) = mean([VirmenData1(6,frameNumber) VirmenData2(6,frameNumber) VirmenData3(6,frameNumber) VirmenData4(6,frameNumber)]);
VirmenCombined(7,frameNumber) = mode([VirmenData1(7,frameNumber) VirmenData2(7,frameNumber) VirmenData3(7,frameNumber) VirmenData4(7,frameNumber)]);
VirmenCombined(8,frameNumber) = max([VirmenData1(8,frameNumber) VirmenData2(8,frameNumber) VirmenData3(8,frameNumber) VirmenData4(8,frameNumber)]);
VirmenCombined(9,frameNumber) = max([VirmenData1(9,frameNumber) VirmenData2(9,frameNumber) VirmenData3(9,frameNumber) VirmenData4(9,frameNumber)]);
VirmenCombined(10,frameNumber) = mode([VirmenData1(10,frameNumber) VirmenData2(10,frameNumber) VirmenData3(10,frameNumber) VirmenData4(10,frameNumber)]);
end
