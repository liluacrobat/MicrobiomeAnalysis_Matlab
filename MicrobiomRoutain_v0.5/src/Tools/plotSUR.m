function plotSUR(t,c,idx,thr)
figure,
hold on
t=t/365;
% thr=thr*365;
c(t>thr)=1;
t(t>thr)=thr;
str=cell(1,length(unique(idx)));
d=0;
rng default
U=unique(idx);
FaceColor =  distinguishable_colors(length(U));
for i=1:length(U)
    if ~isempty(t(idx==U(i)))
        d=d+1;
    str{d}=num2str(U(i));
    iT=t(idx==U(i));
    i_cen=c(idx==U(i));
    [f,x]=ecdf(iT','censoring',i_cen,'function','survivor');
    stairs(x,f,'Color',FaceColor(i,:),'linewidth',2);
    end
end
str=str(1:d);
legend(str);
grid;
end