function plotDistribution(measure)
alpha = 0.3;
x_values = linspace(min(measure(:)),max(measure(:)),100);
pd = pdfKernel(measure);
y_values = pdf(pd,x_values);
area(x_values,y_values,'FaceColor','g','FaceAlpha',alpha);
end
function pd = pdfKernel(x)
% pd = fitdist(x,'Kernel','Kernel','epanechnikov');
pd = fitdist(x,'Kernel');
end
