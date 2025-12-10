function [ SAMPLE_ID,OTU_ID,OTU_TABLE,Type,MAP_ID]=readMap(map_name, SAMPLE_ID,OTU_ID,OTU_TABLE)
clc;
m_id=fopen(map_name,'r');
n=0;
line=fgetl(m_id);
S=strsplit(line,'\t');
while ~feof(m_id)
    n=n+1;
    line=fgetl(m_id);
    S=strsplit(line,'\t');
    MAP_ID{n}=S{1};
    Type(n,1)=str2num(S{4});
end
[IDX1,IDX2]=ID_alig(MAP_ID,SAMPLE_ID);
SAMPLE_ID=SAMPLE_ID(IDX1(IDX1~=0));
OTU_TABLE=OTU_TABLE(:,IDX1(IDX1~=0));
MAP_ID=MAP_ID(IDX1~=0);
Type=Type(IDX1~=0);
IDcheck(MAP_ID,SAMPLE_ID);
end
function IDcheck(X,Y)
flag=0;
if length(X)~=length(Y)
    flag=1;
    display('mismatch');
end
for i=1:length(X)
    if strcmpi(X{i},Y{i})~=1
        flag=1;
        display('mismatch');
    end
end
if flag==0
    display('Pass');
end
end






