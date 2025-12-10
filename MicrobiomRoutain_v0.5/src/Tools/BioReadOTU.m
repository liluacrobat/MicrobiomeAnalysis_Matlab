function OTU_table=BioReadOTU(filename)
fid=fopen(filename);
n=0;
line=fgetl(fid);
X=[];
while ~feof(fid)
    n=n+1;
    line=fgetl(fid);
    S=strsplit(line,'\t');
    if n==1
    ID=S(2:end);
    else
        OTU{n-1,1}=S{1};
        temp=zeros(1,length(S)-1);
        for i=1:length(S)-1
            temp(i)=str2num(S{i+1});
        end
        X=[X;temp];
    end
end
OTU_table.X=X;
OTU_table.ID=ID;
OTU_table.OTU=OTU;
end