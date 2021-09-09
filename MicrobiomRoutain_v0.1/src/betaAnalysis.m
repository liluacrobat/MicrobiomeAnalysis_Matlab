function betaAnalysis(fileName,measure,measure_title,groups_y,header,d,flag_pca)
facecolor = defaultColor(length(header));
d = min(d,3);
if flag_pca==0
    dis = measure;
    pcoa = f_pcoa(dis,1);
    mapped_data = pcoa.scores(:,1:d)';
    power = pcoa.expl(:,1);
    plotFIGURE(mapped_data,groups_y,facecolor);
    set(gca,'FontSize',14);
    xlabel(['PCo1 (' num2str(round(power(1)*10)/10) '%)']);
    ylabel(['PCo2 (' num2str(round(power(2)*10)/10) '%)']);
    if d==3
        ylabel(['PCo3 (' num2str(round(power(2)*10)/10) '%)']);
    end
else
    plotPCA(measure,groups_y,facecolor,d);
end
legend(header);
set(gca,'FontSize',14);
pbaspect([1 1 1]);
box on
keyboard
if flag_pca==1
    plotPDF(gcf,strcat(fileName,'_PCA'));
else
    plotPDF(gcf,strcat(fileName,'_PCoA'));
end
% PERMANOVA test
if flag_pca~=0
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
end