function R=isclose(a,b)
atol = 10^(-8);
rtol = 10^(-5);
R=zeros(length(b),1);
for i=1:length(b)
    R(i)=abs(a- b(i)) <= (atol + rtol * abs(b(i)));
end
end