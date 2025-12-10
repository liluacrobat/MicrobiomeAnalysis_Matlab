function plotGIF(filename,fig)
n=50;
if nargin<1
    filename = 'testnew51.gif';
end
if nargin<2
    fig=1;
end
az = linspace(0,360,n);
for i = 1:n
      view(az(i),30);
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end
end