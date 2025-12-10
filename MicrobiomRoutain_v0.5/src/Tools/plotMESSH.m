function plotMESSH(X,Z)
y=X(2,:);
x=X(1,:);
Z=Z(:)';
[xq,yq] = meshgrid(min(x):0.1:max(x),min(y):0.1:max(y));
vq = griddata(x,y,Z,xq,yq);
figure,
mesh(xq,yq,vq);
colormap(jet);
figure,
contourf(xq,yq,vq);
end