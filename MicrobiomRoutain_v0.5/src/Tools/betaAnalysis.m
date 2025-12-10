function [mapped_data,power] = betaAnalysis(fileName,measure,measure_title,groups_y,header,d,para)
groups_y = groups_y(:);
if ~isfield(para,'Ellipse')
    para.Ellipse=1;
end
if ~isfield(para,'confidence')
    para.confidence=0.95;
end
if ~isfield(para,'facecolor')
    para.facecolor = defaultColor(length(header));
end

facecolor = para.facecolor;
d = min(d,3);
if para.pca==0
    dis = measure;
    pcoa = f_pcoa(dis,1);
    mapped_data = pcoa.scores(:,1:d)';
    power = pcoa.expl(:,1);
    plotFIGURE(mapped_data,groups_y,facecolor,para);
    set(gca,'FontSize',14);
    xlabel(['PCo1 (' num2str(round(power(1)*10)/10) '%)']);
    ylabel(['PCo2 (' num2str(round(power(2)*10)/10) '%)']);
    if d==3
        zlabel(['PCo3 (' num2str(round(power(2)*10)/10) '%)']);
    end
else
    [mapped_data,power] = plotPCA(measure,groups_y,facecolor,d,para);
end
legend(header);
set(gca,'FontSize',14);
pbaspect([1 1 1]);
aa=axis;
aa_x = mean(aa(1:2));
aa_y = mean(aa(3:4));
aa_wx = aa(2)-aa(1);
aa_wy = aa(4)-aa(3);

aa = [aa_x-aa_wx*0.55 aa_x+aa_wx*0.55 aa_y-aa_wy*0.55 aa_y+aa_wy*0.55];
axis(aa);
box on
% hold on
% uy=unique(groups_y);
% for i=1:length(unique(groups_y))
%     id = find(groups_y==uy(i))
% plot(mapped_data(1,id)',mapped_data(2,id)','.-k');
% end
keyboard
if para.pca==1
    plotPDF(gcf,strcat(fileName,'_PCA'));
else
    plotPDF(gcf,strcat(fileName,'_PCoA'));
end
% PERMANOVA test
if para.pca~=0
    dis = f_dis(measure','euc');
end
result_PERMANOVA = f_npManova(dis,groups_y,10000);
result_PERMANOVA_pw = f_npManovaPW(dis,groups_y,10000);

fid = fopen(strcat(fileName,'_PERMANOVA.txt'),'w');
fprintf(fid,'Distance: %s\n',measure_title);
fprintf(fid,'PERMANOVA F: %f\n',result_PERMANOVA(1).F);
fprintf(fid,'PERMANOVA p-value: %f\n',result_PERMANOVA(1).p);

fprintf(fid,'\n');

fprintf(fid,'Pair-wise comparisons\n');
fprintf(fid,'Groups\tGroup_ID\t# of samples\n');
for i=1:length(header)
    fprintf(fid,'%d\t%s\t%d\n',i,header{i},sum(groups_y==i));
end
fprintf(fid,'\n');

for i=1:length(header)
    fprintf(fid,'\t%s',header{i});
end
fprintf(fid,'\n');
fprintf(fid,'Pairwise PERMANOVA with Bonferroni adjusted p-value:\n');
p_mtx = zeros(length(header));
ncmp = length(result_PERMANOVA_pw.p);
for i=1:ncmp
    p = result_PERMANOVA_pw.pairList(i,1);
    q = result_PERMANOVA_pw.pairList(i,2);
    p_mtx(p,q) = result_PERMANOVA_pw.p_bon(i);
end

for i=1:length(header)
    fprintf(fid,'%s',header{i});

    for j=1:length(header)
        if i==j
            fprintf(fid,'\t%f',1);
        else
            if i>j
                fprintf(fid,'\t%f',p_mtx(j,i));
            else
                fprintf(fid,'\t%f',p_mtx(i,j));
            end
        end
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'Pairwise PERMANOVA with Dunn-Sidak adjusted p-value:\n');
p_mtx = zeros(length(header));
ncmp = length(result_PERMANOVA_pw.p);
for i=1:ncmp
    p = result_PERMANOVA_pw.pairList(i,1);
    q = result_PERMANOVA_pw.pairList(i,2);
    p_mtx(p,q) = result_PERMANOVA_pw.p_ds(i);
end

for i=1:length(header)
    fprintf(fid,'%s',header{i});

    for j=1:length(header)
        if i==j
            fprintf(fid,'\t%f',1);
        else
            if i>j
                fprintf(fid,'\t%f',p_mtx(j,i));
            else
                fprintf(fid,'\t%f',p_mtx(i,j));
            end
        end
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'Pairwise PERMANOVA without correction:\n');
p_mtx = zeros(length(header));
ncmp = length(result_PERMANOVA_pw.p);
for i=1:ncmp
    p = result_PERMANOVA_pw.pairList(i,1);
    q = result_PERMANOVA_pw.pairList(i,2);
    p_mtx(p,q) = result_PERMANOVA_pw.p(i);
end

for i=1:length(header)
    fprintf(fid,'%s',header{i});

    for j=1:length(header)
        if i==j
            fprintf(fid,'\t%f',1);
        else
            if i>j
                fprintf(fid,'\t%f',p_mtx(j,i));
            else
                fprintf(fid,'\t%f',p_mtx(i,j));
            end
        end
    end
    fprintf(fid,'\n');
end
end
function h = plotFIGURE(X,Y,FaceColor,para)
% X: principal coordinate m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
if nargin<2||isempty(Y)
    Y=ones(size(X,2),1);
end


U=sort(unique(Y));
str=cell(1,length(U));
if nargin<3
    FaceColor=cbrewer('qual', 'Set1',9);
end
D=size(X,1);
if D>3
    X = X(1:3,:);
    D=3;
end
mapped_data=X;
switch D
    case 2
        if nargin<3
            h = figure;
        end
        for i=1:length(U)
            plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
            hold on
            str{i}=num2str(U(i));
        end
        for i=1:length(U)
            if sum(Y==U(i))>=3 && para.Ellipse==1
                plotEllipse(mapped_data(1:2,Y==U(i)),FaceColor(i,:),para.confidence)
            end
        end
        %         grid
    case 3
        if nargin<3
            h = figure;
        end
        for i=1:length(U)
            plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
            hold on
            str{i}=num2str(U(i));
        end

        grid
end
% boldify
set(gca,'FontSize',14);
% set(gca,'FontSize',14,'FontWeight','Bold');
pbaspect([1 1 1]);
end
function [mapped_data,power,fig,FaceColor]=plotPCA(X,Y,FaceColor,d,para)
% X: m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
% flag: plot ellipses
if nargin<2
    Y = ones(size(X,2),1);
end
if nargin<4
    d = 2;
end
if nargin<5
    flag = 1;
end
U=sort(unique(Y));
str=cell(1,length(U));

if nargin<3
    FaceColor = defaultColor(length(unique(Y)));

end

[mapped_data,~,power]=compute_mapping(X','PCA',d);
mapped_data=mapped_data';

fig = figure;
hold on
if d==3
    for i=1:length(U)
        plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
    end
    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
    zlabel(['PC3 (' num2str(round(power(3)*10000)/100) '%)']);
else
    for i=1:length(U)
        plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
    end
    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
    if para.Ellipse==1
        for i=1:length(U)
            if sum(Y==U(i))>=3
                plotEllipse(mapped_data(1:2,Y==U(i)),FaceColor(i,:),para.confidence)
            end
        end
    end
end
set(gca,'FontSize',14);
% set(gca,'FontSize',14,'FontWeight','Bold');
pbaspect([1 1 1]);
box on
end
function plotEllipse(X,FaceColor,confidence_level)
if nargin<3
    confidence_level = 0.95;
end
data = X';
% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);

% Get the 67% confidence interval error ellipse

chisquare_val = sqrt(chi2inv(confidence_level,2));
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-','color',FaceColor,'linewidth',2);
hold on;
plot(X0,Y0,'+','color',FaceColor,'linewidth',2);
% 
% % Plot the original data
% plot(data(:,1), data(:,2), '.');
% mindata = min(min(data));
% maxdata = max(max(data));
% Xlim([mindata-3, maxdata+3]);
% Ylim([mindata-3, maxdata+3]);
% hold on;

% % Plot the eigenvectors
% quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
% quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);

end