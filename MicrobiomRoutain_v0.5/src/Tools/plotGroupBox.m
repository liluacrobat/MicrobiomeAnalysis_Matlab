function plotGroupBox(visit_number,measure,groups_y,header,medianValue)
nVisit = max(visit_number);


key_tmp = header(groups_y);
meta.pattern = categorical(key_tmp',header);
UPattern = unique(meta.pattern);

for i=1:length(UPattern)
    for j=1:nVisit
        medianValue(i,j) = median(measure(meta.pattern==UPattern(i) & visit_number==j),'omitmissing');
    end
end

nGroup= length(unique(meta.pattern));
stepp = 1/nGroup;
for i=1:nVisit
    vv{i} = ['Visit ' num2str(i)];
end
meta.visit_number = vv(visit_number);
for i=1:nVisit
    vv{i*2-1} = ['Visit ' num2str(i)];
    vv{i*2} = ['Gap ' num2str(i)];
end
vv = vv(1:end-1);
meta.visit_number = categorical(meta.visit_number,vv);

figure,hold on
FaceColor=defaultColor(5);
colormap(FaceColor);
bb = boxchart(meta.visit_number,measure,'GroupByColor',key_tmp,'BoxWidth',0.5);



for i=1:nGroup
    bbcolor(i,:) = bb(i).BoxEdgeColor;
end
xt = xticks;
xt([2 4]) = categorical(" "); % Invisible character (ATL+255), not a space!
xticklabels(xt)
al = legend;
pbaspect([2 1 1]);

baseX = 1:2:nVisit*2;
for i=1:nGroup
    xm = baseX-stepp+stepp*(i-1);
    ym= medianValue(i,:);
    plot(xm,ym,'-o','Color',bbcolor(i,:),'MarkerFaceColor',bbcolor(i,:),'LineWidth',1);
end
set(gca,'FontSize',18);
legend(UPattern);
end