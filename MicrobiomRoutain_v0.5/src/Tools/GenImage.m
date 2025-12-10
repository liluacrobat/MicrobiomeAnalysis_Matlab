function fig =GenImage(X,Color)
Y=zeros(1,length(X),3);
for i=1:3
    for j=1:length(X)
        Y(1,j,i)=Color(X(j),i);
    end
end
fig = figure;
image(Y);
end