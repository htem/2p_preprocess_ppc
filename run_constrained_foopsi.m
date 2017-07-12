function [c,b,c1,g,sn,sp] = run_constrained_foopsi(dF)

numCels = size(dF,1);
c = nan(size(dF));
sp = nan(size(dF));
b = nan(numCels,1);
c1 = nan(numCels,1);
sn = nan(numCels,1);
g = cell(numCels,1);

for cel = find(~isnan(sum(dF,2))')
    sorted_dF = sort(dF(cel,:));
    baseline_adj = mode(round(100*sorted_dF(.05*size(sorted_dF,2):...
        .95*size(sorted_dF,2)))/100);
    dF_zeroed = dF(cel,:) - baseline_adj;
    [c(cel,:),b(cel),c1(cel),g{cel},sn(cel),sp(cel,:)] = constrained_foopsi(dF_zeroed);
end
end

% for cel = find(~isnan(sum(dF,2))')
%     figure('position',[50 50 1600 300]);
%     plot(dF(cel,:))
%     hold on
%     plot(sp(cel,:),'lineWidth',2)
%     plot(c(cel,:),'lineWidth',2)
%     xlim([2000 3000])
%     ylim([-.4 1])
%     close all
% end
