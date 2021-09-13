function y = selBeta(x,sel)
y = x;
y.bc = x.bc(sel,sel);
y.jc = x.jc(sel,sel);
y.thetayc = x.thetayc(sel,sel);
end