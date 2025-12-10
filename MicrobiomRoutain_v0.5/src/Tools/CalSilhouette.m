function [AveSilh, Silh, numSilhouette, meanSilhouette,fig]=CalSilhouette(X, cidx,geoDis)
%% ------------------------------------------------------------------------
% apply K-means clustering and confirm the robustness using Silhouette
%% ------------------------------------------------------------------------
if nargin<3
    Silh = silhouette(X',cidx,'cityblock');
else
    reGeo = @(x,Y,w)w(x,Y);
    Silh = silhouette(X',cidx,reGeo,geoDis);
end

% % order the clustering label
fig = figure;
hold on
n = 0;
temp = Silh;
cluster_legend = [];
c_num = 0;
numClusters = max(cidx);
facecolor = defaultColor(numClusters);
tf = facecolor;
for k = numClusters:-1:1
    temp_c = sort(temp(cidx==k),'ascend');
    numSilhouette(k) = length(temp_c);
    meanSilhouette(k) = mean(temp_c);
    barh(n+(1:length(temp_c)),temp_c,'FaceColor',tf(k,:),'EdgeColor',...
        tf(k,:));
    n = n+length(temp_c);
    c_num = c_num+1;
    cluster_legend{c_num} = ['Cluster ' num2str(k)];
end
% plot([mean(temp) mean(temp)],[1 n],'--r','linewidth',2);
xlabel('Silhouette Width');
ylabel('');
legend(cluster_legend);
length(find(temp<0))
AveSilh = mean(Silh);
% xticks(-0.1:0.1:0.5);
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ycolor',[1 1 1])
ax = axis;
ax(3) = -1;
ax(4)= length(cidx)+1;
axis(ax);
% axis([-0.101 0.5 -3 275]);
pbaspect([1 3 1])
end
