function [Y,header,FaceColor]= fun_catColor(mm)
% X: principal coordinate m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
mm = mm(:)';
Y=mm-min(mm)+1;
U=sort(unique(Y));
ColorG = defaultColor(length(U));
ColorG(ColorG>1)=1;
ColorG(ColorG<0)=0;
header=cell(1,length(U));
for i=1:length(header)
    header{i}=['Level' num2str(i)];
end
FaceColor=ColorG;

end