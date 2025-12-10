function dis = f_dis_yc(x)
[~,n] = size(x);
x = CalRel(x);
dis = zeros(n);
for i=1:n
    for j=i+1:n
        x_a = x(:,i);
        x_b = x(:,j);
        sab = x_a'*x_b;
        bot = sum((x_a-x_b).^2)+sum(x_a.*x_b);
        dis(i,j) = 1-sab/bot;
        if sab==0
            dis(i,j) = 1;
        end
        dis(j,i) = dis(i,j);
    end
end
end
