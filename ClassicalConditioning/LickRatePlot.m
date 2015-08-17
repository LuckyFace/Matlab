clc; clear all; close all;

load('D:\Cheetah_data\classical_conditioning\PV1\2015-08-11_14-00-54\Events.mat');

fHandle = figure('PaperUnits', 'centimeters', 'PaperPosition', [2 2 8.9/2 6.88/2]);
lineColor = {[0.8 0.4 0],[0.4 0 0.8],[0.8 0 0.4],[0 0.4 0.8]};

nCue = 4;
yMax = ceil(max(lickSemConv(:)));
axes('Position',[0.1 0.1 0.8 0.4]);
for iCue = 1:nCue
    hold on;
    plot([0.5 0.5], [0 100], 'LineStyle', ':', 'LineWidth', 0.2, 'Color', [0.8 0.8 0.8]);
    plot([1.5 1.5], [0 100], 'LineStyle', ':', 'LineWidth', 0.2, 'Color', [0.8 0.8 0.8]);
    plot([4 4], [0 100], 'LineStyle', ':', 'LineWidth', 0.2, 'Color', [0.8 0.8 0.8]);
    rectangle('Position', [0.5 yMax*0.95 1 yMax*0.05], 'LineStyle', 'none', 'FaceColor', [0 0 0]);
    rectangle('Position', [4 yMax*0.95 0.1 yMax*0.05], 'LineStyle', 'none', 'FaceColor', [0 0 0]);
    fill(lickSemBin, lickSemConv(iCue,:), lineColor{iCue}, 'LineStyle', 'none');
    plot(lickBin, lickMeanConv(iCue,:), ...
        'Color', lineColor{iCue}, ...
        'LineWidth', 1);
    set(gca, 'box', 'off', 'TickDir', 'out', 'LineWidth', 0.2, 'FontSize', 5, ...
        'XLim', [0 8], 'XTick', [0 0.5 1.5 4 8], ...
        'YLim', [0 yMax], 'YTick', [0 yMax]);
end
print(fHandle,'-depsc','D:\Cloud\project\classical_conditioning\fig\LickRatePlot.eps');