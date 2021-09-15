function y = selBeta(x,sel)
if isempty(x)
    y = x;
else
    
    y = x;
    y.bc = x.bc(sel,sel);
    y.jc = x.jc(sel,sel);
    y.thetayc = x.thetayc(sel,sel);
end
end