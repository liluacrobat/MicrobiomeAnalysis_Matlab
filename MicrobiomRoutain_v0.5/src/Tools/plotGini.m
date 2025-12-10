function plotGini(gini,name,m)
figure,
hold on
sel = gini>0.01;
gini = gini(sel);
name = name(sel);
% m = [zeros(d,size(m,2));m];
[~,gc] = max(m,[],2);
% gc(1:d) = max(gc)+1;

gc = gc(sel,:);
[gini,idx] = sort(gini,'descend');
name = name(idx);
% if length(idx)>30
%     idx = idx(1:30);
%     gini = gini(1:30);
%     name = name(idx);
% end
% name{2} = 'Peptostreptococcaceae saphenum';
% name{8} = 'Anaerolineae bacterium HMT439';
% name{16} = 'Peptostreptococcaceae bacterium HMT369';
% name{11} = 'Peptostreptococcaceae nodatum';
% name{22} = 'Peptoniphilaceae bacterium HMT113';
% name{28} = 'Peptostreptococcaceae brachy';

% name{4} = 'Peptostreptococcaceae nodatum';
% name{7} = 'Smoking years';
% name{14} = 'Mchol';
% 
% name{21} = 'Anaerolineae bacterium HMT 439';
% name{26} = 'Peptostreptococcaceae saphenum';


% name{4} = 'Anaerolineae bacterium HMT 439';
% name{9} = 'Peptostreptococcaceae nodatum';
% name{11} = 'Peptostreptococcaceae saphenum';
% name{13} = 'Lifetime smoked packs per year';
% name{16} = 'Saccharibacteria (TM7) bacterium HMT 347';
% name{17} = 'Veillonellaceae bacterium HMT 132';
% name{27} = 'Lifetime total ounces of ethanol for usu'; 
gc = gc(idx,:);

u = unique(gc);
gini(gc==u(1)) = gini(gc==u(1))*(-1);

[gini,idx2] = sort(gini,'ascend');
gc = gc(idx2,:);
name = name(idx2);
x = 1:length(idx);
facecolor = defaultColor;
yticks(x);

% xtickangle(45);

axis([min(gini)*1.1 max(gini)*1.1 0.5 max(x)+0.5 ]);
for i=length(u):-1:1
    ax = x(gc==u(i));
    ay = gini(gc==u(i));
    for d=1
        barh(ax(d),ay(d),'FaceColor',facecolor(mod(i-1,7)+1,:),'BarWidth',0.5);
    end
end
for i=length(u):-1:1
    ax = x(gc==u(i));
    ay = gini(gc==u(i));
    for d=1:length(ax)
        barh(ax(d),ay(d),'FaceColor',facecolor(mod(i-1,7)+1,:),'BarWidth',0.5);
    end
end

yticklabels(name);
% yticks([-50 0 50 100]);
% yticklabels({'50','0','50','100'});
% pbaspect([2 1 1]);
pbaspect([1 1.2 1 ]);
xticks(-2:2:10);
a =axis;
a(1)=-3;
axis(a)
set(gca,'FontSize',14);
xlabel('Importance Score','FontSize',16);
end