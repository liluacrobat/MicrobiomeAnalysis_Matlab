function plot_PCA(Data,Label,dim)
L_t=Label;
if dim==3
    [mapped_data,mapping]=compute_mapping(Data','PCA',3);
    mapped_data=mapped_data';
    figure,
%     i=0;
%     index=find(L_t==i);
%     plot3(mapped_data(1,index),mapped_data(2,index),mapped_data(3,index),'xr');
%     hold on
    for i=1:max(Label)
    index=find(L_t==i);
    plot3(mapped_data(1,index),mapped_data(2,index),mapped_data(3,index),'o');
    hold on
    end
    colormap(colorcube);
    grid;
else
    [mapped_data,mapping]=compute_mapping(Data','PCA',2);
    mapped_data=mapped_data';
    figure,
    for i=1:max(Label)
    index=find(L_t==i);
    plot(mapped_data(1,index),mapped_data(2,index),'x');
    hold on
    end
    colormap(colorcube);
    grid;
end