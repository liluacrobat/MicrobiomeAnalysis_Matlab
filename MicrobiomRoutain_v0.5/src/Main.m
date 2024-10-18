function Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routain analysis of microbiome study
% Matlab 2023b
% Analysis
%  Basic summary of samples
%  Top 20 OTUs
%  Alpha diversity
%  Beta diversity
%  Group comparison using LEfSe
%  Heqatmap of abundance of taxa selected by LEfSe
%
% Note:
% lines having comments starting with '% ***' are ones need to change based on
% specific study
%
% Author: Lu Li
% Version:0.5
% Last updated 10/18/2024
% Required R packages: R.matlab, fossil, vegan
% Setting:
% Create an environment of Mothur using anconda: conda create mothur mothur,biom-format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;close all hidden
% Setting environment for the R, may need to change accordingly in Windows
% For example: '/' in paths to '\'
% system
PATH = getenv('PATH');
if ~contains(PATH,':/Library/Frameworks/R.framework/Resources/bin','IgnoreCase',true)
    setenv('PATH', [PATH ':/Library/Frameworks/R.framework/Resources/bin']);
end

%% Preprocessing
%% Specify project name
proj = 'MgPreventive';  % *** project name

%% Parameters for data preprocessing (1:Yes,0:No)
% Need to specify the parameters and run only once. Then, you can set
% flag.pre = 0, the program will use the calculated files

flag.pre = 0;      % *** Perform pre-processsing
flag.seq_plot = 0; % *** Summarize the # of reads in each sample
flag.alpha = 0;    % *** Calculate alpha diversity
flag.beta = 0;     % *** Calculate beta diversity

%% Set the working directory
% !!! Need to change to where the MicrobiomRoutain_v0.4 folder locates
WD = '/Users/luli/Documents/MATLAB/OrganizePipeline/Pit_Mg/MicrobiomRoutain_v0.5'; % *** Wprking directory

Rdir = [WD '/result'];
Ddir = [WD '/data'];
Cdir = [WD '/src'];

addpath(genpath(WD));
mkdir(Rdir);
mkdir(Ddir );
cd(Cdir);

%% Initial file name of OTU table and clinical meta data
OTU_tbl_name = 'plate_16S.OTU.reformated_L7.txt'; % *** Change to your OTU table name. Make sure the table is put in "MicrobiomRoutain_v0.4/data" folder
meta_tbl_name = 'MgPreventive_meta.txt'; % *** Change to your meta table name. Make sure the table is put in "MicrobiomRoutain_v0.4/data" folder

meta_key = 1; % Specify the column in meta file used as the sample ID to match the OTU table. No need to change unless the first colum of meta table is not used to match the sample ID in OTU table
dc = [0.00,0.45,0.74]; % default color (blue)

nthr = 3000; % *** Sequencing depth threshold. Change as necessary. Low depth: 1000, moderate: 3000, high depth: 10000

mkdir(strcat(Rdir,'/',proj)); % Create the directory to store results

%% Pre-processing
if flag.pre==1

    % Load OTU table and mapping file
    OTU_tbl = loadOTUtbl(strcat(Ddir,'/',OTU_tbl_name));
    meta_tbl = readtable(strcat(Ddir,'/',meta_tbl_name),'delimiter','\t');

    % Match meta data and OTU table
    key = table2array(meta_tbl(:,meta_key));

    %     key = addPrefix(key,'sample',1); % Add prefix if the ID
    %     used in mapping file is different

    [match_idx,~] = AlignMem(key,OTU_tbl.sample_id);
    meta_tbl = meta_tbl(match_idx~=0,:);
    x_tbl = filterBysample(OTU_tbl,match_idx);

    % List the samples not included in OTU table
    disp(['# of samples not included in OTU table: ' num2str(sum(match_idx==0))]);
    disp(key(match_idx==0));

    total_counts = sum(x_tbl.counts,1);

    if flag.seq_plot==1

        % Summarize the sequencing depth
        % Plot with all data. The sample name was displayed along x axis
        figure,semilogy(total_counts,'ko','MarkerFaceColor',dc,'MarkeredgeColor',dc,'MarkerSize',8,'linewidth',1);
        ax = axis;
        ax(1) = 0;
        ax(2) = length(total_counts)+1;
        ax(3) = max(ax(3)/10,1);
        hold on
        plot([-10 length(total_counts)+10],[nthr nthr]);
        axis(ax);
        pbaspect([3,1,1]);
        set(gca,'FontSize',9);
        ylabel('Total reads','Fontweight','Normal')
        xticks(1:length(total_counts));
        xticklabels(strrep(x_tbl.sample_id,'_',' '));
        keyboard % pause so that you could modify the plot before it is saved.
        plotPDF(gcf,strcat(Rdir,'/',proj,'/',proj,'_reads_per_sample_dot'));

        % Plot highlight the samples below depth threshold
        figure,semilogy(total_counts,'ko','MarkerFaceColor',dc,'MarkeredgeColor',dc,'MarkerSize',8,'linewidth',1);
        ax = axis;
        ax(1) = 0;
        ax(2) = length(total_counts)+1;
        ax(3) = max(ax(3)/10,1);
        hold on
        plot([-10 length(total_counts)+10],[nthr nthr]);
        axis(ax);
        pbaspect([2,1,1]);
        set(gca,'FontSize',12);
        xlabel('Samples','Fontweight','Normal');
        ylabel('Total reads','Fontweight','Normal')
        xticks([]);
        s = rmStr(x_tbl.sample_id,'_');
        for i=1:length(total_counts)
            if total_counts(i)<nthr
                text(i,total_counts(i),s{i},'FontSize',14);
            end
        end
        keyboard
        plotPDF(gcf,strcat(Rdir,'/',proj,'/',proj,'_reads_per_sample_dot_w_sample'));
    end

    % Filter samples by sequencing depth
    depth_sel = total_counts>=nthr;
    meta_ready = meta_tbl(depth_sel,:);
    x_ready = filterBysample(x_tbl,depth_sel);
    x_ready.clr = CalCLR(x_ready.counts);
    x_ready.rel = CalRel(x_ready.counts);

    disp(strcat('Minimum total reads: ',num2str(min(sum(x_ready.counts,1)))));
    writeResultTable(strcat(Rdir,'/',proj,'/',proj,'_table_ready'),x_ready.tax,x_ready.sample_id,x_ready.counts);

    if flag.alpha==1

        % Calculate alpha diversity
        para4alpha.rare = 2; % Use rarefaction in alpha diversity caluclation
        para4alpha.num = min(sum(x_ready.counts,1)); % Rarefaction to the lowest depth of samples passed filtering
        para4alpha.rep = 100;  % Repeat 100 times, report the average value

        [observed_otu, shannon_index, chao1_index,simposon,Result] = CalAlphaDiversity(...
            x_ready.counts,para4alpha);
        mkdir(strcat(Rdir,'/',proj));
        save(strcat(Rdir,'/',proj,'/',proj,'_alpha'),'observed_otu', 'shannon_index', ...
            'chao1_index', 'simposon','Result')
    else
        observed_otu = [];
        shannon_index = [];
        chao1_index = [];
        simposon = [];
        if isfile(strcat(Rdir,'/',proj,'/',proj,'_alpha.mat'))
            load(strcat(Rdir,'/',proj,'/',proj,'_alpha.mat'),'observed_otu', 'shannon_index', ...
                'chao1_index', 'simposon','Result');
            disp('Load alpha diversity... Succeed.')
        else
            disp('Error: Alpha diversity has not been calculated yet.');
            keyboard;
        end
    end

    if flag.beta==1

        % Calculate beta diversity use Mothur

        cd(strcat(Rdir,'/',proj));
        genBetaSH_Gen(strcat(Rdir,'/',proj,'/',proj,'_table_ready')) % Generate shell script
        system('chmod 777 tmp_cal_beta.sh');
        keyboard
        % Open terminal and change working directory to the folder of 'tmp_cal_beta.sh' using 'cd' command
        % Run the code in script  'tmp_cal_beta.sh' with Mothur. Please modify according to your personal settings
        cd(Cdir);

        % Load the beta diversity
        if isfile(strcat(Rdir,'/',proj,'/',proj,'_table_ready.braycurtis.userLabel.lt.ave.dist'))
            beta_dist.bc = loadBeta(strcat(Rdir,'/',proj,'/',proj,'_table_ready.braycurtis.userLabel.lt.ave.dist'),x_ready.sample_id);
        end
        if isfile(strcat(Rdir,'/',proj,'/',proj,'_table_ready.jclass.userLabel.lt.ave.dist'))
            beta_dist.jc = loadBeta(strcat(Rdir,'/',proj,'/',proj,'_table_ready.jclass.userLabel.lt.ave.dist'),x_ready.sample_id);
        end
        if isfile(strcat(Rdir,'/',proj,'/',proj,'_table_ready.thetayc.userLabel.lt.ave.dist'))
            beta_dist.thetayc = loadBeta(strcat(Rdir,'/',proj,'/',proj,'_table_ready.thetayc.userLabel.lt.ave.dist'),x_ready.sample_id);
        end
    else
        beta_dist = [];
    end

    save(strcat(Rdir,'/',proj,'/',proj,'_data_ready'),'observed_otu', 'shannon_index', ...
        'chao1_index', 'simposon','beta_dist','meta_ready','x_ready');
else
    load(strcat(Rdir,'/',proj,'/',proj,'_data_ready'),'observed_otu', 'shannon_index', ...
        'chao1_index', 'simposon','beta_dist','meta_ready','x_ready');
end

%% Code of Analysis
para.Rdir = strcat(Rdir,'/',proj);
para.proj = proj;
para.test = 'Mg'; % Specify the test name
para.prim = meta_key; % Sample ID column in mapping file

%% Specify the samples for analysis
% If sel = [], it will use all the samples for comparison
key_treatment = catLabel(meta_ready,3); % use colum 3 of meta_ready for grouping
% Select a subgroup of samples, use & (and), | (or) operators to creat a selection
sel =  key_treatment.y==2 | key_treatment.y==3 | key_treatment.y==4 ;

%% Specify the column for comparison
para.cmp = 3; % Column in mapping used to define comparison groups

%% Specify the tested included in analysis (1:Yes,0:No)
para.test_top = 0;         % Calculate the Top20 OTUs
para.test_alpha = 0;       % Alpha diversity comparison
para.test_beta = 0;        % Beta diversity comparison
para.test_lefse = 0;       % LEfSe analysis
para.test_heatmap = 1;        % Hierarchial clustering

%% Load calculated apha and beta diversity
alpha.ob_otu = observed_otu;
alpha.shannon = shannon_index;
alpha.simposon = simposon;
alpha.chao1 = chao1_index;
beta = beta_dist;

routainAnalysis(x_ready,alpha,beta,meta_ready,sel,para);
end
function routainAnalysis(x_input,alpha,beta,meta_input,sel,para)
cfolder = strcat(para.proj,'_',para.test);
mkdir(para.Rdir,cfolder);
cfolder = strcat(para.Rdir,'/',cfolder);
if ~isfield(para,'prim')
    para.prim = 1;
end
if isempty(sel)
    sel = ones(size(x_input.counts,2),1)==1;
end

x = filterBysample(x_input,sel);
meta = meta_input(sel,:);

alpha = selAlpha(alpha,sel);
beta = selBeta(beta,sel);

n_key = length(para.cmp); % number of variables for comparisons
n_s = size(x.counts,2);   % number of samples

key_ls = cell(1,n_key);
label = zeros(n_s,n_key);
legend_ls = cell(1,n_key);
for i=1:n_key
    key_tmp = catLabel(meta,para.cmp(i));
    key_ls{i} = key_tmp;
    label(:,i) = key_tmp.y;
    legend_ls{i} = key_tmp.legend;
    if ~iscell(key_tmp.legend)
        tmp = cell(size(key_tmp.legend));
        for k=1:length(tmp)
            tmp{k} = num2str(key_tmp.legend(k));
        end
        legend_ls{i} = tmp;
    end
end
if isfield(para,'legend_ls')
    legend_ls = para.legend_ls;
end
[groups_y, header] = genGroups(label,legend_ls);

%% Change the order of headers as necessary
% [header,groups_y] = reOrder([2 1 3 4],header,groups_y);

if para.test_top ~= 0
    %% Top20 OTUs with highest average relative abunacne of all samples

    m_rel = mean(x.rel,2);
    [~,rel_rank] = sort(m_rel,'descend');
    Top20.rel = x.rel(rel_rank(1:20),:);
    Top20.tax = x.tax(rel_rank(1:20));

    % Relative abundance of Top 20 OTUs per samples
    mkdir(cfolder,'Top20');
    fileName = strcat(cfolder,'/Top20/',para.proj,'_',para.test,'_top20rel_per_sample');
    write2Table(fileName,[Top20.tax(:); 'Others'],x.sample_id,[Top20.rel;1-sum(Top20.rel,1)],'num');

    % Average relative abundance of Top 20 OTUs per group

    MRel = zeros(20,length(header));
    for i=1:length(header)
        MRel(:,i) = mean(Top20.rel(:,groups_y==i),2);
    end
    fileName = strcat(cfolder,'/Top20/',para.proj,'_',para.test,'_top20rel_per_group');
    write2Table(fileName,[Top20.tax(:); 'Others'], header,[MRel;1-sum(MRel,1)],'num');
end
if para.test_alpha ~= 0

    %% Alpha diversity
    mkdir(cfolder,'Alpha_div');
    measure_alpha = alpha.ob_otu;
    measure_alpha_title = 'Obsearved OTUs';
    fileName = strcat(cfolder,'/Alpha_div/',para.proj,'_',para.test,'_alpha_observed_otu');
    alphaAnalysis_w_dot(fileName,measure_alpha,measure_alpha_title,groups_y,header)

    measure_alpha = alpha.shannon;
    measure_alpha_title = 'Shannon Index';
    fileName = strcat(cfolder,'/Alpha_div/',para.proj,'_',para.test,'_alpha_shannon');
    alphaAnalysis_w_dot(fileName,measure_alpha,measure_alpha_title,groups_y,header)


    measure_alpha = alpha.chao1;
    measure_alpha_title = 'Chao1 Index';
    fileName = strcat(cfolder,'/Alpha_div/','/',para.proj,'_',para.test,'_alpha_chao1');
    alphaAnalysis_w_dot(fileName,measure_alpha,measure_alpha_title,groups_y,header)


    measure_alpha = alpha.simposon;
    measure_alpha_title = 'Simpson Index';
    fileName = strcat(cfolder,'/Alpha_div/','/',para.proj,'_',para.test,'_alpha_simpson');
    alphaAnalysis_w_dot(fileName,measure_alpha,measure_alpha_title,groups_y,header)
end

if para.test_beta ~= 0
    %% Beta diversity
    mkdir(cfolder,'Beta_div');
    measure_beta = x.clr;
    measure_beta_title = 'euclidean';
    para_beta.pca=1;
    para_beta.Ellipse=1; % Draw ellipse
    fileName = strcat(cfolder,'/Beta_div/',para.proj,'_',para.test,'_',measure_beta_title);
    betaAnalysis(fileName,measure_beta,measure_beta_title,groups_y,header,2,para_beta)

    measure_beta = beta.bc;
    measure_beta_title = 'bray-curtis';
    para_beta.pca=0;
    para_beta.Ellipse=1;
    fileName = strcat(cfolder,'/Beta_div/',para.proj,'_',para.test,'_',measure_beta_title);
    betaAnalysis(fileName,measure_beta,measure_beta_title,groups_y,header,2,para_beta)

    measure_beta = beta.jc;
    measure_beta_title = 'jaccard';
    para_beta.pca=0;
    para_beta.Ellipse=1;
    fileName = strcat(cfolder,'/Beta_div/',para.proj,'_',para.test,'_',measure_beta_title);
    betaAnalysis(fileName,measure_beta,measure_beta_title,groups_y,header,2,para_beta)

    measure_beta = beta.thetayc;
    measure_beta_title = 'ThetaYC';
    para_beta.pca=0;
    para_beta.Ellipse=1;
    fileName = strcat(cfolder,'/Beta_div/',para.proj,'_',para.test,'_',measure_beta_title);
    betaAnalysis(fileName,measure_beta,measure_beta_title,groups_y,header,2,para_beta)
end

if para.test_lefse~=0
    %% LEfSe analysis
    % Pairwise comparison
    mkdir(cfolder,'LEfSe');

    % Specify the selection criteria
    para4sel.thr_fre = 0.15; % Prevelence >15% in at leat one group
    para4sel.thr_rel = 2*10^-4; % Average relative abundance >0.2% in at leat one group
    com = nchoosek(1:length(header),2); % Determine pairwise comparisons
    for k=1:size(com,1)
        ki = com(k,1);
        kj = com(k,2);
        sel_com = groups_y==ki | groups_y==kj;
        % Select samples of compariosn group pairs
        x_sub = filterBysample(x,sel_com);
        para4sel.y = groups_y(sel_com);

        % Filter OTU by the defined criteria of prevelence and abundance
        sel4lefse = filterOTU(x_sub.rel,para4sel);
        x4lefse = x_sub.rel(sel4lefse,:)*10^6;
        y4lefse = para4sel.y;



        %% Following can be used to simplify taxa name
        para4lefse.tax = simpleName(x_sub.tax(sel4lefse));
        %% Otherwise use the full name
        % para4lefse.tax = x_sub.tax(sel4lefse);

        para4lefse.y_legend = header;
        para4lefse.rpath = strcat(cfolder,'/LEfSe/',para.proj,'_',para.test,'_LEfSe_V1_',header{ki},'_V2',header{kj});
        para4lefse.plot = 1;
        [order,score,Enrich_lefse,tax_lefse,rpath_lefse] = CalLEfSe(x4lefse,y4lefse,para4lefse);
    end
end

if isfield(para,'test_heatmap')
    if para.test_heatmap~=0
        %% Heatmap analysis
        % Filtering by relative abundances
        mkdir(cfolder,'heatmap');

        %% Specify the group of taxa used for heatmap

        % Example of picking the top 20 varate taxa
        taxa_std = std(x.rel,[],2);
        [~,tax_sorted_idx] = sort(taxa_std,'descend');

        picked = zeros(length(taxa_std),1);
        picked(tax_sorted_idx(1:20)) = 1;

        sel4heat = picked==1;

        disp(['# of OTUs selected for heatmap analysis: ' num2str(sum(sel4heat))]);
        tax = simpleName(x.tax);

        % Heatmap, each row normalized to z-score for visualization
        % !!! Normalized data is not used for clustering !!!
        para4heat.rpath = strcat(cfolder,'/heatmap/',para.proj,'_',para.test,'_heatmap_std');
        para4heat.plotrow = 1;
        para4heat.rowname = tax(sel4heat);
        para4heat.plotcol = 1;
        para4heat.colname = x.sample_id;
        para4heat.plotGroup =1;
        para4heat.Group = label;
        para4heat.Group_legend = legend_ls;
        para4heat.norm = 'std';
        [row_sel,col_sel] = plotHeatmap(x.clr(sel4heat,:),para4heat);

        % % Heatmap, each row normalized to 0 to 1 for visualization
        % para4heat.rpath = strcat(cfolder,'/heatmap/',para.proj,'_',para.test,'_heatmap_01');
        % para4heat.plotrow = 1;
        % para4heat.rowname = tax(sel4heat);
        % para4heat.plotcol = 1;
        % para4heat.colname = x.sample_id;
        % para4heat.plotGroup =1;
        % para4heat.Group = label;
        % para4heat.Group_legend = legend_ls;
        % para4heat.norm = '01';
        % [row_sel,col_sel] = plotHeatmap(x.clr(sel4heat,:),para4heat);
    end
end
end













