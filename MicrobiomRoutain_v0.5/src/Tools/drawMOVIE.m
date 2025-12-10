function drawMOVIE(k)
A(1:length(1:k:360))=struct('cdata',[],'colormap',[]);
for az=1:k:360
    view(az,30);
   drawnow
    f=getframe;
    f=frame2im(f);
    [X,map]=rgb2ind(f,256);
    if az==1
        
        imwrite(X,map,'ex_imwrite.gif');
    else
        imwrite(X,map,'ex_imwrite.gif','WriteMode','Append');
    end
end
end