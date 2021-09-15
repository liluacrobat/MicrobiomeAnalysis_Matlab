function y = selBeta(x,sel)

y = x;
if ~isempty(x)
    if isempty(x.bc)
        y.bc = [];
    else
        y.bc = x.bc(sel,sel);
    end
    if isempty(x.jc)
        y.jc = [];
    else
        y.jc = x.jc(sel,sel);
    end
    if isempty(x.thetayc)
        y.thetayc = [];
    else
        y.thetayc = x.thetayc(sel,sel);
    end
end
end