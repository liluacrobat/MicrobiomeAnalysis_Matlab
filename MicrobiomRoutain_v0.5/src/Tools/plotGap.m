function [h, x, y] = plotGap(eva_Gap)
% Plot Gap statistic result

Gap = eva_Gap.ExpectedLogW-eva_Gap.LogW;
Gap_SE = Gap-eva_Gap.SE;
delta_Gap = Gap(1:end-1)-Gap_SE(2:end);

% Data to plot
x = eva_Gap.InspectedK(1:end-1);
y = delta_Gap;

h = figure; hold on;
bar(x, y);
xlabel('Number of Clusters')
ylabel('Gap(k)-[Gap(k+1)-SE]')


