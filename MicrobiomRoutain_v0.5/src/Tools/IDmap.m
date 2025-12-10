function nd=IDmap(id,or,new)
nd=id;
nd(id~=0)=nan;
for i=1:length(or)
    nd(id==or(i))=new(i);
end
end