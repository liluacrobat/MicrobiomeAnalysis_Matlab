function [h,composition,fig]= plotComposition(subtypeIndex, clusterIndex, type,FaceColor)
% Plot subytype composition in each cluster
%
%
%% =======================================================

composition = crosstab(clusterIndex, subtypeIndex);
composition=composition./repmat(sum(composition,2),1,size(composition,2))*100;
[numCluster, numSubtype] = size(composition);
if nargin<4
    FaceColor=ColorSel(1);
    switch length(unique(subtypeIndex))
        case 6
            temp=FaceColor;
            FaceColor(4,:)=temp(6,:);
            FaceColor(5:6,:)=temp(4:5,:);
        case 3
            temp=FaceColor;
            FaceColor(3,:)=temp(4,:);
            FaceColor(5:6,:)=temp(4:5,:);
        case 2
            temp=FaceColor;
            FaceColor(2,:)=temp(4,:);
            FaceColor(5:6,:)=temp(4:5,:);
        case 4
            tf=cbrewer('qual', 'Paired',12);
            tf2=cbrewer('qual', 'Set3',12);
            FaceColor(1,:)=tf(4,:);
            FaceColor(2,:)=tf2(3,:);
            FaceColor(3,:)=tf2(5,:);
            FaceColor(4,:)=tf2(4,:);
            
            %         FaceColor(1,:)=tf(4,:);
            %         FaceColor(2,:)=tf(10,:);
            %         FaceColor(3,:)=tf(1,:);
            %         FaceColor(4,:)=tf(6,:);
        otherwise
            FaceColor =  cbrewer('qual', 'Set2',length(unique(subtypeIndex)));
    end
end
subtypeColor=FaceColor;
switch type
    
    case 'barh'
        % Calculate cluster composition
        fig = figure;
        h = barh(flipud(composition), 'stacked');
        %         h = barh(flipud(composition./repmat(sum(composition,2),1,length(subtypeLabel))), 'stacked');
        %         xlabel('# of Samples')
        xlabel('Percentage of Samples (%)');
        ylabel('Cluster ID');
        set(gca, 'YTickLabel', cellstr(num2str((numCluster:-1:1)')));
        % Set color
        for n = 1:numSubtype; set(h(n), 'FaceColor', subtypeColor(n,:));end
        %         legend(subtypeLabel, 'Location', 'southeast');
        %         boldify2;
        % Percentage plot
        %         figure;
        %
        %         set(gca, 'YTickLabel', cellstr(num2str((numClusters:-1:1)')));
        %         for n = 1:length(subtypeLabel)
        %             set(h(n), 'FaceColor', subtypeColor(n,:));
        %             baseline_handle = get(h(n),'BaseLine');
        %             set(baseline_handle,'LineStyle','none');
        %         end
        %         boldify2;
        %         legend(subtypeLabel, 'Location', 'BestOutside');
        
    case 'bar'
        % Calculate cluster composition
        fig = figure;
        h = bar(composition, 'stacked');
        ylabel('Percentage of Samples (%)');
        xlabel('Cluster ID');
        %         set(gca, 'XTickLabel', cellstr(num2str((numCluster:-1:1)')));
        % Set color
        for n = 1:numSubtype; set(h(n), 'FaceColor', subtypeColor(n,:));end
        %         boldify2;
        
        
    case 'pie'
        numCluster = length(unique(clusterIndex));
        
        for n = 1:numCluster
            h = figure;
            %             subplot(1, numCluster, n)
            hp=pie(composition(n,:));
            color_temp=FaceColor(composition(n,:)>0,:);
            for k=1:length(hp)/2
                set(hp(k*2-1), 'FaceColor', color_temp(k,:));
            end
            boldify_pca
        end
        
end