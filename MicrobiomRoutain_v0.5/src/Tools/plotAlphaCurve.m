function plotAlphaCurve(x,depth,y,groups)
figure,
hold on
u = unique(y);
n = length(depth);
facecolor = defaultColor;
for i=1:length(u)
    xs{i} = x(y==u(i),:);
    xs_m{i} = mean(xs{i},1,'omitnan');
    xs_s{i} = std(xs{i},1,'omitnan');
    plot(depth,xs_m{i},'linewidth',1.5,'color',facecolor(mod(i-1,7)+1,:));
end

% for i=1:length(u)
%     errorbar(depth,xs_m{i},xs_s{i},'linewidth',1.5,'color',facecolor(mod(i-1,7)+1,:))
% end
set(gca,'FontSize',14);
xlabel('Sequencing Depth');
legend(groups,'location','southeast');
pbaspect([1 1 1])
end
