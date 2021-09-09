function [row_sel,col_sel] = plotHeatmap(X,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heatmap analysis
% X: m x n data, m features, n samples
% para: parameter for analysis
%   norm - normalization of data
%          '01': normalized to 0 to 1(default)|'std': normalized to 0 mean 1 std
%   rlinkage - linkage for hierarchcal clustering of columns (default: 'average')
%   clinkage - linkage for hierarchcal clustering of columns (default: 'average')
%   drow - distance used for hierarchcal clustering of rows (default: 'spearman')
%   dcol - distance used for hierarchcal clustering of columns (default: 'euclidean')
%   plotGroup - whether plot colorbar of groups (default: 0)
%   Group - labels of samples (only used when plotLabel=1)
%   Group_legend - group legend (only used when plotLabel=1)
%   colormap - color map of the heatmap
%   cluster - dimension for cluster 'all'(default)|'column'|'row'
%   plotrow - plot row names (default: 1)
%   rowname - row names
%   colname - rowcolumn names
%   plotcol - plot column names (default: 1)
%   rpath - path to store results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
if ~isfield(para,'plotcol')
    para.plotcol = 0;
end
if ~isfield(para,'group_col')
    para.group_col = 1;
end
if ~isfield(para,'plotrow')
    para.plotrow = 0;
end
if ~isfield(para,'rlinkage')
    para.rlinkage = 'average';
end
if ~isfield(para,'clinkage')
    para.clinkage = 'average';
end
if ~isfield(para,'drow')
    para.drow = 'spearman';
end
if ~isfield(para,'dcol')
    para.dcol = 'euclidean';
end
if ~isfield(para,'norm')
    para.norm = 'std';
end
if ~isfield(para,'plotGroup')
    para.plotGroup = 0';
end

X_vis = dataNorm(X,para.norm);

if ~isfield(para,'colormap')
    cmap = [165,0,38
        215,48,39
        244,109,67
        253,174,97
        254,224,144
        255,255,191
        224,243,248
        171,217,233
        116,173,209
        69,117,180
        49,54,149]/255;
    cmap = flip(cmap);
end

if ~isfield(para,'cluster')
    para.cluster = 'all';
end
[~,n] = size(X);

pd_c = pdist(X',para.dcol);
Z_c = linkage(pd_c,para.clinkage);
[~,~,outperm_c] = dendrogram(Z_c,n);
plotPDF(gcf,strcat(para.rpath,'_dendrogram'));

flag_sym = strcmp(para.norm,'std')==1;
if para.group_col==1
    CGobj = clustergram(X,'Cluster',para.cluster,'Colormap',cmap,'Symmetric',false,...
        'OptimalLeafOrder',false,'RowPDist',para.drow,'ColumnPDist',para.dcol,...
        'Linkage',{para.rlinkage,para.clinkage});
    % keyboard
    plot(CGobj);
    plotPDF(gcf,strcat(para.rpath,'_clustergram_w_den'));
    
    [~,col_sel] = heatmapLabels(CGobj);
else
    CGobj = [];
    col_sel = 1:size(X,2);
end

CGobj_vis = clustergram(X_vis(:,col_sel),'Cluster','column','Colormap',cmap,'Symmetric',flag_sym,...
    'OptimalLeafOrder',false,'RowPDist',para.drow,'ColumnPDist',para.dcol,...
    'Linkage',{para.rlinkage,para.clinkage},'DisplayRatio',[0.001 0.15]);

[row_sel,~] = heatmapLabels(CGobj_vis);
a=para.rowname(row_sel);
% if strcmp(para.norm,'std')
%     fig =figure;
%     colormap(fig,cmap);
%     image(X_vis(row_sel,col_sel),'CDataMapping','scaled');
% end

if para.plotrow==1
    set(CGobj_vis,'RowLabels',para.rowname);
else
    set(CGobj_vis,'RowLabels',[]);
end
if para.plotcol==1
    set(CGobj_vis,'ColumnLabels',para.colname(col_sel));
else
    set(CGobj_vis,'ColumnLabels',[]);
end
plot(CGobj_vis);
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 12)
keyboard
plotPDF(gcf,strcat(para.rpath,'_clustergram_normalized'));

save(strcat(para.rpath,'_clustergram_normalizedresults'),'CGobj','CGobj_vis',...
    'col_sel','row_sel','Z_c','para','outperm_c');


if para.plotGroup==1
    ng = size(para.Group,2);
    num_gc = zeros(1,ng);
    for i=1:ng
        num_gc(i) = length(para.Group_legend{i});
    end
    facecolorGroup = defaultGroupColor(num_gc);
    genGroupBar(para.Group(col_sel,:),facecolorGroup,para.Group_legend,para.rpath);
    
end
end
function genGroupBar(X,ColorG,Name,rpath)
[m,n] = size(X);

gap = 0.2;
step = round(1/gap);
img = zeros(n*step,m,3);
for k=1:n
    
    facecolor = ColorG{k};
    
    for i=1:3
        for j=1:m
            img(((k-1)*step+1):(k*step-1),j,i)= facecolor(X(j,k),i);
        end
    end
    for i=1:3
        for j=1:m
            img(k*step,j,i)= 255;
        end
    end
    
end
figure,
image(img);
axis off
plotPDF(gcf,strcat(rpath,'_group_bar'));
for k=1:length(Name)
    Name_t = Name{k};
    n = length(Name_t);
    img_t = zeros(n*step,1,3);
    facecolor = ColorG{k};
    
    figure,hold on
    
    for p=1:n
        plot(1,1,'s','MarkerSize',16,'MarkerFaceColor',facecolor(p,:),'MarkerEdgeColor',facecolor(p,:));
        legend(Name_t,'FontSize',16);
    end
    plotPDF(gcf,strcat(rpath,'_group_bar_legend_',num2str(k)));
end
end

function colorLine(h,linecolor)
for i=1:length(h)
    h(i).Color = linecolor;
end
end
function Y = dataNorm(X,method)
%% normalize the data so that each feature is comparable
switch method
    case '01'
        % Normalize to [0, 1]
        [MIN,~] = min(X,[],2);
        [MAX,~] = max(X,[],2);
        STA = MAX-MIN;
        Y = X-MIN*ones(1,size(X,2));
        Y(STA~=0,:) = Y(STA~=0,:)./(STA(STA~=0,1)*ones(1,size(Y,2)));
        Y(STA==0,:)=0;
    case 'std'
        % Normalize to Z-score
        M = mean(X,2);
        STD=std(X,[],2);
        Y = (X-M*ones(1,size(X,2)))./(STD*ones(1,size(X,2)));
        Y(STD==0,:) = 0;
end
end
function [row_sel_num,col_sel_num] = heatmapLabels(CGobj)
row_sel = CGobj.RowLabels;
row_sel = flipud(row_sel);
col_sel = CGobj.ColumnLabels;
for i=1:length(row_sel)
    row_sel_num(i)=str2num(row_sel{i});
end
for i=1:length(col_sel)
    col_sel_num(i)=str2num(col_sel{i});
end
end