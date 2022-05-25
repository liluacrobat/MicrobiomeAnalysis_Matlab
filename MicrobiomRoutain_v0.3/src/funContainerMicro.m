classdef funContainerMicro
    methods
        function R=CalRel(obj,X,p)
            if nargin<3
                R=X./repmat(sum(X,1),size(X,1),1);
            else
                Q=zeros(1,size(X,2));
                for i=1:length(Q)
                    temp=X(:,i);
                    S=cumsum(sort(temp,'ascend'));
                    Q(i)=quantile(S,p);
                    if  Q(i)==0
                        Q(i)=eps;
                    end
                end
                R=X./repmat(Q,size(X,1),1);
            end
        end
        function R1=CalCLR(obj,X)
            R1=X;
            if min(X(:))==0
                X=X+1;
            end
            for i=1:size(X,2)
                temp=X(:,i);
                temp2=temp(temp~=0);
                %     temp(temp==0)=0.1;
                m = geomean(temp2);
                R1(:,i)=log(temp)-log(m);
            end
        end
        function f_ls = addPrefix(obj,id_ls,prefix,num_flag)
            if nargin<4
                num_flag=0;
            end
            
            n = length(id_ls);
            f_ls = cell(n,1);
            if num_flag==1
                
                for i=1:n
                    f_ls{i} = strcat(prefix,num2str(id_ls(i)));
                end
            else
                for i=1:n
                    f_ls{i} = strcat(prefix,id_ls(i));
                end
                
            end
        end
        function OTU_tbl = loadOTUtbl(obj,tbl_name,tax_flag)
            % read OTU table and extract sample id, taxonomy, and counts
            if nargin<3
                tax_flag = 0;
            else
                tax_flag = 1;
            end
            tbl = readtable(tbl_name,'delimiter','\t');
            OTU_tbl.tax = table2array(tbl(:,1));
            OTU_tbl.sample_id = tbl.Properties.VariableNames(2:end-tax_flag);
            OTU_tbl.counts = table2array(tbl(:,2:end-tax_flag));
            if tax_flag==1
                OTU_tbl.otu = OTU_tbl.tax;
                OTU_tbl.tax = table2array(tbl(:,end));
            end
        end
        function f_tbl = filterBysample(obj,OTU_tbl,idx)
            if length(unique(idx))~=2
                idx = idx(idx~=0);
            end
            f_tbl.sample_id = OTU_tbl.sample_id(idx);
            counts = OTU_tbl.counts(:,idx);
            otu_sel = sum(counts,2)>0;
            f_tbl.counts = counts(otu_sel,:);
            f_tbl.tax = OTU_tbl.tax(otu_sel);
            if isfield(OTU_tbl,'otu')
                f_tbl.otu = OTU_tbl.otu(otu_sel);
            end
            if isfield(OTU_tbl,'clr')
                f_tbl.clr = OTU_tbl.clr(otu_sel,idx);
            end
            if isfield(OTU_tbl,'rel')
                f_tbl.rel = OTU_tbl.rel(otu_sel,idx);
            end
        end
        function plotPDF(obj,fig,filename)
            savefig(fig, strcat(filename, '.fig'))
            set(fig,'Units','Inches');
            pos = get(fig,'Position');
            set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(fig,strcat(filename, '.pdf'),'-dpdf', '-r600');
            % print(fig,strcat(filename, '.eps'),'-depsc', '-r600')
        end
        function writeResultTable(obj,map_name,OTU_ID,Sample_ID,Data,tax)
            fid=fopen([map_name '.txt'],'w');
            fprintf(fid,'OTU ID');
            for i=1:length(Sample_ID)
                fprintf(fid,'\t%s',Sample_ID{i});
            end
            if nargin==6
                fprintf(fid,'\ttaxonomy');
            end
            fprintf(fid,'\n');
            for i=1:length(OTU_ID)
                fprintf(fid,'%s',OTU_ID{i});
                for j=1:length(Sample_ID)
                    fprintf(fid,'\t%f',Data(i,j));
                end
                if nargin==6
                    fprintf(fid,'\t%s',tax{i});
                end
                fprintf(fid,'\n');
            end
        end
        
        
    end
end
