function cmap = genColormap(climit,L)
climit_low = climit(2,:);
climit_high = climit(1,:);
lstick  = linspace(0,1,11)';
l1 = lstick;
l2 = flip(lstick);
w = [1 1 1];
cmap = (1-l1)*climit_low + l1*climit_high + w-0.5*(climit_low+climit_high);
% mat1 = repmat(lstick,1,3);
% mat2 = repmat(flip(lstick),1,3);
% 
% cmap = [
%     l2.*(w-climit_high)+climit_high;
%     1 1 1;
%     l1.*(w-climit_low)+climit_low;
% ];
end