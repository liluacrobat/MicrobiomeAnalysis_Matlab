function [Mean, E_vector,Projection, E_value] = PrincipalComAna(DATA,row,col,dim,plotfig,percent)
 
%====================================================================================
%This fuction is to reduce the data dimesionality using Principal Components Analysis
%Useage [Mean, E_vector,Projection]=PCA_general(DATA,row,col,dim)
%Input:
%        DATA     : the normalized real data in the form of data=[u1 u2,...] 
%        Row/Col  : an optional to vasulize the eigenvector
%        dim      :   The desirable or finally demensionality
%        plotfig    : 1,plot the image; 0, not plot the image 
%        percent    : if percent~=0; dimensionality is determined by energy
%Output:
%        Mean       :  The mean of the pool data
%        E_vector   : The resulting eigenvector
%        Projection : the projected data
%=====================================================================================
%Copyright(c) 2002, Yijun Sun 
%SAL, The University of Florida
Mean=mean(DATA,2);           % a column vector
DATA=DATA-Mean*ones(1,size(DATA,2));
if size(DATA,1)>size(DATA,2)
   covmatrix=(DATA')*DATA;    % U'U
   [V,dia]=eig(covmatrix);
   E_vector=DATA*V;           %eigenvector matrix
   for n = 1:size(E_vector,2)
      E_vector(:,n) =  E_vector(:,n)/norm(E_vector(:,n));
   end
end
 
if size(DATA,2)>=size(DATA,1)
   covmatrix=(DATA)*(DATA');% UU'
   [V,dia]=eig(covmatrix);
   E_vector=V;
end
 
for n=1:size(dia,1)
   E_value(n)=abs(real(dia(n,n)));
end
 
[E_value,index]=sort(E_value);
E_value=fliplr(E_value);
index=fliplr(index);
E_vector=E_vector(:,index);
 
% figure(1)
% set(gca,'fontsize',13)
% %E_value1=fliplr(E_value);
% plot(E_value,'r-o');grid on;
%  
% if percent~=0
%    totalE=sum(E_value);
%    count=0;
%    Per=0;
%    while Per<percent
%       count=count+1;
%       Per=sum(E_value(1:count))/totalE;
%    end
%    dim=count;
% end
% hold on
% stem(dim,E_value(dim));hold off;
%-----visionlize the eigenvector-----
 
if plotfig==1
   Nu=10;
   %t=sin(linspace(0,1,31).');
   %[dum,tau,xi]=ambifunb(t);
   for n=1:Nu
      eigImage=reshape(E_vector(:,n),row,col);%keyboard
      figure(20)
      subplot(ceil(Nu/3),3,n)
      imagesc(eigImage);shading interp;
      axis equal; axis tight;
      axis off
      title([ 'eigenvector ' num2str(n)])
      %contour(2*tau,xi,eigImage,16); 
      %grid
      %xlabel('Delay'); ylabel('Doppler'); shading interp      
      %mesh(eigImage)
      %dispimage(abs(eigImage),-100,0.03* [0 size(eigImage,2)-1],0.03* [0 size(eigImage,1)-1],jet)
   end
end
 
Projection=(((E_vector(:,1:dim))')*(DATA));
% if fs==1 %feature selection
%     Projection=DATA'*E_vector;%(((E_vector)')*(DATA));
%     NumberF=dim;
%     FurtherStep=10;
%     optional=[1,8];
%     disp('>>> Fisher criteria + SFFS') 
%     SFFSdata={Projection(1:No(1),:),Projection(1:No(2),:)};
%     feature=SFFS(NumberF,FurtherStep,SFFSdata , optional,'LDA');
%     feature
%     E_vector=E_vector(:,feature);
%     Projection=(E_vector')*(DATA);
% end