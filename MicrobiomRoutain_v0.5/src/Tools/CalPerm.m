function [sen,spe] = CalPerm(x,y)
pos = max(y);
neg = min(y);
sen = sum(x(y==pos)==pos)/sum(y==pos);
spe = sum(x(y==neg)==neg)/sum(y==neg);
end