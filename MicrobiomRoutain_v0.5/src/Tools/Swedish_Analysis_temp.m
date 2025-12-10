function Swedish_Analysis_temp
%% ========================================================================
% Exploring the dynamics of Crohn's disease based on behavior
% Lu Li
% update history: 2/27/2018
%% ========================================================================
clc;clear;close all;
%% loading

addpath('/Users/Lu/Documents/MATLAB/MicroDynamic/NatureMicrobiology/SwedishData/data');
addpath('/Users/Lu/Documents/MATLAB/MicroDynamic/NatureMicrobiology/SwedishData/KEGG_map');
load('SwedishData');
Behavior_legend = {'HC','Inflammatory (B1)','Stricturing (B2)','Penetrating (B3)'};
DIA_legend = {'HC','CCD','ICD'};
DIA_SUB_legend = {'HC','CCD','ICD-nr','ICD-r'};
%% ------------------------------------------------------------------------
% SAMPLE_ID: sample's ID
% PATIENT: patient ID
% Behavior: CD behavior
%            1:not applicable (HC)
%            2:Non-stricturing, non-penetrating (B1)
%            3:Stricturing (B2)
%            4:Penetrating (B3)
% Location: disease location
%            0:not applicable (HC)
%            1:Ileal (L1)
%            2:Colonic (L2)
%            3:Ileocolonic (L3)
%            4:Ileal and Upper-GI (L1+L4)
%            5:Ileocolonic and Upper-GI (L3+L4)
% DIA: CD diagnosis
%       1:HC
%       2:colonic Crohn's disease (CCD)
%       3:ileal Crohn's disease (ICD)
% DIA_SUB: CD diagnosis with surgery information
%           1:HC
%           2:CCD
%           3:ICD without ileocaecal resection (ICD-nr)
%           4:ICD with ileocaecal resection (ICD-r)
% bmi: Body mass index (BMI)
% sex: gender
%       1:male
%       2:female
% dur: duration of CD
% perianal: perianal disease
%            0:not applicable (HC)
%            1:no
%            2:yes
% Time: sampling time
% calprotectin: fecal calprotecin concentration
% counts_per_sample: number of counts for each sample
% TABLE_L6: OTU table in genus level, D*N matrix (D:number of genus;
%                                                 N:number of smaples )
% Tax_L6: taxonomy in genus level
% TABLE_L7: OTU table in species level, D*N matrix (D:number of species;
%                                                   N:number of smaples )
% Tax_L7: taxonomy in species level
% TABLE_OTU: OTU table, D*N matrix (D:number of OTUs; N:number of smaples )
% Tax_OTU: taxonomy of OTUs
% FUN: funation table of KEGG orthology (KO), D*N matrix (D:number of KO;
%                                                         N:number of smaples )
% FUN_name: KEGG orthology (KO)
% FUN_desp: description of KO
% PATHWAY: table of KEGG pathways, D*N matrix (D:number of pathways;
%                                              N:number of smaples )
% PATHWAY_name: KEGG pathways
%
% FUN_desp: description of KO
% PATHWAY: table of KEGG pathways, D*N matrix (D:number of pathways;
%                                              N:number of smaples )
% EC_Table: EC number (KO)
% EC_ID: ID of EC
%
% KO_Table: KEGG abundance
% KO_ID: KO ID
%
% ECP_Table: pathway abundance (KO)
% ECP_ID: ID of pathways
%% ------------------------------------------------------------------------

%% pre-processing
option.count = 0;
if option.count == 1
    %% plot the sample counts
    color_map = cbrewer('qual', 'Paired',12);
    figure,
    semilogy(counts_per_sample,'o','MarkerFaceColor',color_map(1,:),'MarkerSize',...
        8,'MarkerEdgeColor',color_map(2,:));
    hold on
    semilogy([0 length(counts_per_sample)],[10000 10000],'--r','linewidth',3);
    xlabel('Samples');
    ylabel('OTU Counts')
    boldify_line
    axis([0 length(counts_per_sample) 10^1 10^7])
end

% exclude the samples with counts<10000
count_sel = counts_per_sample>=10000;
SAMPLE_ID = SAMPLE_ID(count_sel);
PATIENT = PATIENT(count_sel);
Time = Time(count_sel);
Behavior = Behavior(count_sel);
Location = Location(count_sel);
DIA = DIA(count_sel);
DIA_SUB = DIA_SUB(count_sel);
bmi = bmi(count_sel);
sex = sex(count_sel);
dur = dur(count_sel);
perianal = perianal(count_sel);
calprotectin = calprotectin(count_sel);
counts_per_sample = counts_per_sample(count_sel);
TABLE_L6 = TABLE_L6(:,count_sel);
TABLE_L7 = TABLE_L7(:,count_sel);
TABLE_OTU = TABLE_OTU(:,count_sel);
FUN = FUN(:,count_sel);
PATHWAY = PATHWAY(:,count_sel);
EC_Table = EC_TABLE(:,count_sel);
ECP_Table = ECP_TABLE(:,count_sel);
KO_Table = KO_TABLE(:,count_sel);
% calculate relative abundance
Rel_OTU = CalRel(TABLE_OTU);

% remove the OTUs with counts <=20 to facilitate the feature selection
[~,idx20] = RmSingle(TABLE_OTU,20);
Rel_filtered = Rel_OTU(idx20>0,:);

% apply log10 transformation
LogRel = log10(Rel_filtered+10^(-6));
training = LogRel;
Label = Behavior;
Label(Behavior==4)=3;

option.cross = 0;
if option.cross == 1
    %% apply 10-fold crossvalidation to select the parameter for LOGO
    % set the range of parameter
    lam =logspace((-5),2,15);
    % apply 10-fold crossvalidation
    ACC = cell(1,length(lam));
    Ave_ACC = zeros(1,length(lam));
    for pi = 1:length(lam)
        Para4logo.distance = 'block';
        Para4logo.plotfigure = 0;
        Para4logo.sigma = 11;
        Para4logo.lambda = lam(pi);
        Para4logo.soft = 0;
        Para4logo.kernel = 'parabolic';
        Para4logo.seed = 98;
        Para4logo.maxit = 10;
        Para4logo.folds = 10;
        
        %         [ACC{pi}, ~, ~] = CrossVail_LOGO(training, Label,...
        %             Para4logo);
        %         Ave_ACC(pi) = mean(ACC{pi});
        %         save('SwedishCrossValidationResult','ACC','Ave_ACC','lam');
    end
end

option.cross_plot = 0;
if option.cross_plot == 1
    %% plot the 10-fold crossvalidation results
    load('SW_CV_Result_fold_10_or_sep98','Ave_ACC','lam');
    Ave_ACC = Ave_ACC(:,98);
    Testing_Error=100-Ave_ACC;
    figure,
    semilogx(lam,Testing_Error,'-o','Color','b','MarkerFaceColor','b','MarkerSize',8,'MarkerEdgeColor','b','linewidth',2);
    hold on
    [~, idx_best] = max(Ave_ACC);
    plot([lam(idx_best) lam(idx_best)],[0 80],'--r','linewidth',3);
    xlabel('Regulation Parameter \lambda');
    ylabel('Cross-validation Error (%)');
    axis([10^(-5) 100 0 80]);
    boldify_line;
end

option.feature = 0;
if option.feature == 1
    %% feature selection;
    load('SwedishCrossValidationResult','Ave_ACC','lam');
    [~, idx_best] = max(Ave_ACC);
    
    Para4logo.distance = 'block';
    Para4logo.plotfigure = 0;
    Para4logo.sigma = 11;
    Para4logo.lambda = lam(idx_best);
    Para4logo.soft = 0;
    Para4logo.kernel = 'parabolic';
    Para4logo.maxit = 10;
    
    [Weight,~,~] = Logo_kernel(training,Label,Para4logo);
    save('SwedishFeatureSelectionResult', 'Weight');
end

load('SwedishFeatureSelectionResult', 'Weight');
% normalize the feature weight
Weight = Weight/sum(Weight);
threshold = 0.001;
ids = Weight>threshold;

option.feature_plot = 0;
if option.feature_plot == 1
    %% plot normalized feature weight
    [Weight_sorted,~]=sort(Weight,'descend');
    figure,semilogx(Weight_sorted,'-b');
    hold on
    boldify_pca
    plot([sum(ids) sum(ids)],[0 max(Weight)],'--r','linewidth',3);
    ylabel('Normalized Feature Weight');
    xlabel('OTUs');
    axis([1 10^5 0 max(Weight)])
end

% data after feature selection
X = LogRel(ids,:);
option.pca_plot = 0;
if option.pca_plot == 1
    %% plot PCA with tunnel
    % PCA analysis
    Y = Behavior;
    [mapped_data,~,power_pca]=compute_mapping(X','PCA',3);
    mapped_data=mapped_data';
    X_PCA = mapped_data(1:3,:);
    
    % estimate the backbone using PrincipalTree method
    para4PTree = [];
    para4PTree.maxIter = 100;
    %     %% estimate the parameter using elbow method
    %     lambda_range = logspace(-2,log10(1000),20);
    %     ln = length(lambda_range);
    %     Error4PTree = zeros(1,ln);
    %     CurveLength4PTree = zeros(1,ln);
    %     for i = 1:ln
    %         para4PTree.lambda = lambda_range(i);
    %         [~,history4PTree] = principal_tree(X_PCA, para4PTree);
    %         Error4PTree(i)=history4PTree.mse(end);
    %         CurveLength4PTree(i)= history4PTree.length(end);
    %     end
    %     figure,plot(CurveLength4PTree,Error4PTree,'o');
    %     lambda4PTree_opt = ElbowPosition(Error4PTree, CurveLength4PTree,...
    %         lambda_range,1,0,2);
    %             save('SwedishPTreeElbow', 'lambda4PTree_opt', 'Error4PTree', ...
    %                 'CurveLength4PTree', 'lambda_range')
    
    load('SwedishPTreeElbow')
    lambda4PTree_opt = ElbowPosition(Error4PTree, CurveLength4PTree,...
        lambda_range,1,0,2);
    para4PTree.lambda = lambda4PTree_opt;%5;
    [model4PTree,~] = principal_tree(X_PCA, para4PTree);
    % postprocessing
    subtype=Label2Type(Y,Behavior_legend);
    [PCurve4PTree,stree4PTree,ProjDis4PTree] = PrincipalTreePost(model4PTree.MU',...
        X_PCA', subtype);
    model4PTree.MU = PCurve4PTree;
    model4PTree.stree = stree4PTree;
    % exclude the 15% samples farest to the backbone
    [~,idx]=sort(ProjDis4PTree);
    ratio = 0.90;
    num = round(length(DIA)*ratio);
    
    X_PCA_filtered = X_PCA(:,idx(1:num));
    Y_filtered = Y(idx(1:num));
    ProjDis_filtered = ProjDis4PTree(idx(1:num));
    % plot the tunnel
    backbone = 1; % 0:don't plot backbone; 1:plot backbone
    save('SwedishTunnelResult', 'X_PCA_filtered', 'model4PTree',...
        'Y_filtered', 'power_pca', 'PCurve4PTree', 'stree4PTree',...
        'ProjDis_filtered');
    load('SwedishTunnelResult');
    plotPTree(X_PCA_filtered, model4PTree, Y_filtered, power_pca,[],backbone);
    drewTunnel_Surf(X_PCA_filtered,PCurve4PTree,stree4PTree,model4PTree.sigma,...
        ProjDis_filtered)
    legend(Behavior_legend);
    view(12,22)
end

option.ddr_elbow = 0;
if option.ddr_elbow == 1
    %% estimating parameters of DDR using elbow method
    Para4ddr.gamma = 2; % Regulation parameter for square error
    Para4ddr.dim = 3;
    Para4ddr.maxIter = 100;
    Para4ddr.eps = 1e-9;
    Para4ddr.knn = 1;
    
    % fix lambda=1, estimate sigma
    f_sig = 1;
    offset_sig = 10;
    Para4ddr.lambda =1;  % tree strength
    %     rang_sig=[0.0001 0.001 0.005 0.01 0.05 0.1 0.15 0.5 1.5 2.5 4 5 10 15 20];
    %     [CurveLength_sig,Error_sig,rang_sig,sigma_opt] = ElbowDDR(X, Para4ddr,...
    %         f_sig, rang_sig, offset_sig);
    %             save('Swedish_DDR_ELBOW_sig','Error_sig', 'CurveLength_sig', 'rang_sig');
    load('Swedish_DDR_ELBOW_sig');
    sigma_opt = ElbowPosition(Error_sig, CurveLength_sig,rang_sig,1, f_sig, offset_sig);
    
    % fix sigma=0.5, estimate lambda
    f_sig = 0;
    offset_lam = 10;
    %     Para4ddr.sigma = 1;
    %     rang_lam=sort([2 3 5 7 10 20 50 150 500 1000 2000 3000 5000  7000 10000]);
    %     [CurveLength_lam,Error_lam,rang_lam,lambda_opt] = ElbowDDR(X, Para4ddr,...
    %         f_sig, rang_lam, offset_lam);
    %         save('Swedish_DDR_ELBOW_lam','Error_lam', 'CurveLength_lam', 'rang_lam');
    load('Swedish_DDR_ELBOW_lam');
    lambda_opt = ElbowPosition(Error_lam, CurveLength_lam,rang_lam,1, f_sig, offset_lam);
end
option.cluster = 0;
if option.cluster == 1
    %% clustering analysis using K-means
    
    % determine the number of clusters using gap statistics
    %     numClusters = gapKmeans(X,['Gap/Swedish_Gap_Clustering_Raw'],1);
    load('Swedish_Gap_Clustering.mat','eva_Gap','numClusters');
    [h, gapX, gapY] = plotGap(eva_Gap);
    boldify_line
    grid
    xticks(1:9);
    ylabh = get(gca,'ylabel');
    get(ylabh,'position')
    set(ylabh, 'Position', get(ylabh,'position')+ ...
        [-0.15, 0, 0]);
    
    % confirm the robustness of clustering using Silhouette method
    [cidx, AveSilh, Silh] = SilhouetteKmeans(X, numClusters);
    axis([-0.15 0.5 0 275]);
    plotPCA3D_color(X,cidx);
    save('SwedishClusteringResult','cidx');
end
load('SwedishClusteringResult','cidx');
cluster_legend = {'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', ...
    'Cluster 5'};
option.ddr = 0;
if option.ddr == 1
    %% structure learning using DDR method
    Para4ddr.gamma = 2;
    Para4ddr.dim = 3;
    Para4ddr.maxIter = 100;
    Para4ddr.eps = 1e-9;
    
    Para4ddr.lambda = 150;
    Para4ddr.sigma = 0.5;
    %         [W, Z,stree, PTree, ~] = DDRTree(X, Para4ddr);
    %         subtype = Label2Type(Behavior,Behavior_legend);
    %         [extracted_path, extracted_pathDist, PTree_etd, extracted_curve,...
    %             stree_etd] = StruPostProcess(PTree, Z, subtype);
    
    %         save('SwedishDDR_Reslt','W', 'Z', 'stree', 'PTree', 'extracted_path',...
    %             'extracted_pathDist', 'PTree_etd', 'extracted_curve', 'stree_etd');
    
    load('SwedishDDR_Reslt');
    flag_PTree = 1;
    plotDDRtree(Z, PTree_etd, stree_etd, Behavior,Behavior_legend, flag_PTree)
    boldify_line
    grid
    view(-100,10);
    flag_PTree = 1;
    plotDDRtree(Z, PTree_etd, stree_etd, cidx,cluster_legend, flag_PTree)
    boldify_line
    grid
    view(-100,10);
    % organize the order of progression paths
    seq = [3 1 4 2];
    for i=1:length(seq)
        Path{i} = extracted_path{seq(i)};
        PathDis{i} = extracted_pathDist{seq(i)};
        Y_path = ones(size(DIA))*1;
        Y_path(Path{i})=2;
        plotDDRtree(Z, PTree_etd, stree_etd,Y_path, {'Other',['Path' num2str(i)]});
        grid
        boldify_line
    end
    %         save('SwedishDDR_Path','Path','PathDis');
    plotDDRtree(Z, PTree_etd, stree_etd, cidx,cluster_legend, flag_PTree)
end
load('SwedishDDR_Path','Path', 'PathDis');
% save SW_Merge4Analy_data
M_OTU = mean(Rel_OTU,2);
Fre_OTU = CalFre(Rel_OTU);

Rel_L7 = CalRel(TABLE_L7);
M_L7 = mean(Rel_L7,2);
Fre_L7 = CalFre(Rel_L7);

Rel_EC = CalRel(EC_Table);
M_EC = mean(EC_Table,2);
Fre_EC = CalFre(Rel_EC);

Rel_ECP = CalRel(ECP_Table);
M_ECP = mean(ECP_Table,2);
Fre_ECP = CalFre(Rel_ECP);

Rel_KO = CalRel(KO_Table);
M_KO = mean(KO_Table,2);
Fre_KO = CalFre(Rel_KO);

option.plot_cluster = 0;
if option.plot_cluster == 1
    %% Clustering Composition
    FaceColor=ColorList('ibd_pca_be2');
    % FaceColor1 = FaceColor*255;
    [h,comp] = plotComposition(Behavior, cidx, 'bar',FaceColor);
    C_legend = {'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5'};
    set(gca,'xticklabel',C_legend);
    xtickangle(-30)
    xlabel('');
    boldify_pca
    legend(Behavior_legend);
    
    FaceColor=ColorList('ibd_pca_dia');
    [h,comp]= plotComposition(DIA_SUB, cidx, 'bar');
    set(gca,'xticklabel',C_legend);
    xtickangle(-30)
    xlabel('');
    axis([0.2 5.8 0 100])
    boldify_pca
    legend(DIA_SUB_legend);
end
load('SW_KOP')
[idx12,~] = AlignID(SAMPLE_ID,KOP_Sample);
KOP = KOP(:,idx12);
Rel_KOP = CalRel(KOP);
M_KOP = mean(KOP,2);
Fre_KOP = CalFre(Rel_KOP);
option.correlation = 1;
if option.correlation == 1
    %% correlation analysis
%     for i=1:length(Path)
%         [corr_rho_KOP{i},corr_p_KOP{i},corr_q_KOP{i},corr_M_KOP{i},...
%             corr_Fre_KOP{i}] = PathTest( Path{i}, PathDis{i}, Rel_KOP,'dsfdr',0.05);
%     end
%     save('IBD_spearman_KOP_dsfdr1','corr_rho_KOP','corr_p_KOP','corr_q_KOP','corr_M_KOP',...
%         'corr_Fre_KOP','M_KOP','Rel_KOP','Path','PathDis');
%     % species level
%     for i=1:length(Path)
%         [corr_rho_L7{i},corr_p_L7{i},corr_q_L7{i},corr_M_L7{i},corr_Fre_L7{i}]=PathTest( Path{i}, PathDis{i},Rel_L7,'dsfdr',0.05);
%     end
%     save('IBD_spearman_L7_dsfdr2','corr_rho_L7','corr_p_L7','corr_q_L7','corr_M_L7','corr_Fre_L7','M_L7','Rel_L7','Path','PathDis');
    

%
    %OTUs
    
        for i=1:length(Path)
            [corr_rho_OTU{i},corr_p_OTU{i},corr_q_OTU{i},corr_M_OTU{i},corr_Fre_OTU{i}]=PathTest( Path{i}, PathDis{i},Rel_OTU,'dsfdr',0.05);
        end
        save('IBD_spearman_OTU_dsfdr0','corr_rho_OTU','corr_p_OTU','corr_q_OTU','corr_M_OTU','corr_Fre_OTU','M_OTU','Rel_OTU','Path','PathDis');
    
    
    %     % functional
    %             for i=1:length(Path)
    %                 [corr_rho_fun{i},corr_p_fun{i},corr_q_fun{i},corr_M_fun{i},corr_Fre_fun{i}]=PathTest( Path{i}, PathDis{i},Rel_fun,'dsfdr',0.05);
    %             end
    %             save('IBD_spearman_fun_dsfdr','corr_rho_fun','corr_p_fun','corr_q_fun','corr_M_fun','corr_Fre_fun','M_fun','Rel_fun','Path','PathDis');
    
    %% pathway activities
    %         for i=1:length(Path)
    %             [corr_rho_pathway{i},corr_p_pathway{i},corr_q_pathway{i},corr_M_pathway{i},corr_Fre_pathway{i}]=PathTest( Path{i}, PathDis{i}, Rel_pathway,'dsfdr',0.05);
    %         end
    %         save('IBD_spearman_pathway_dsfdr','corr_rho_pathway','corr_p_pathway','corr_q_pathway','corr_M_pathway','corr_Fre_pathway','M_pathway','Rel_pathway','Path','PathDis');
    
    
    %     %% ECP
    %     for i=1:length(Path)
    %         [corr_rho_ECP{i},corr_p_ECP{i},corr_q_ECP{i},corr_M_ECP{i},...
    %             corr_Fre_ECP{i}] = PathTest( Path{i}, PathDis{i}, Rel_ECP,'dsfdr',0.05);
    %     end
    %     save('IBD_spearman_ECP_dsfdr2','corr_rho_ECP','corr_p_ECP','corr_q_ECP','corr_M_ECP','corr_Fre_ECP','M_ECP','Rel_ECP','Path','PathDis');
    %% EC
%     for i=1:length(Path)
%         [corr_rho_EC{i},corr_p_EC{i},corr_q_EC{i},corr_M_EC{i},...
%             corr_Fre_EC{i}] = PathTest( Path{i}, PathDis{i}, Rel_EC,'dsfdr',0.05);
%     end
%     save('IBD_spearman_EC_dsfdr2','corr_rho_EC','corr_p_EC','corr_q_EC','corr_M_EC',...
%         'corr_Fre_EC','M_EC','Rel_EC','Path','PathDis');
    %% KO
%     for i=1:length(Path)
%         [corr_rho_KO{i},corr_p_KO{i},corr_q_KO{i},corr_M_KO{i},...
%             corr_Fre_KO{i}] = PathTest( Path{i}, PathDis{i}, Rel_KO,'dsfdr',0.05);
%     end
%     save('IBD_spearman_KO_dsfdr2','corr_rho_KO','corr_p_KO','corr_q_KO','corr_M_KO',...
%         'corr_Fre_KO','M_KO','Rel_KO','Path','PathDis');
end
%
% load('IBD_spearman_L7_dsfdr');
% load('IBD_spearman_OTU_dsfdr');
% load('IBD_spearman_ECP_dsfdr');
% load('IBD_spearman_EC_dsfdr');
% load('IBD_spearman_KO_dsfdr');
% %% plot correlated OTUs
% TaxSimplified = SimplifyNameRefine(Tax_OTU);
%
% for i=1:length(Path)
% %     sel_idx_OTU{i} = corr_q_OTU{i}<0.05 &  M_OTU>0.0001 ;%& Fre_OTU>0.3;
%     sel_idx_OTU{i} = corr_q_OTU{i}<0.01 &  corr_M_OTU{i}>0.001;% & corr_Fre_OTU{i}>0.15;
%     display(['# of selected OTUs along Path ' num2str(i) ': ' num2str(sum(sel_idx_OTU{i}))]);
%     [Ptax_OTU{i},Ntax_OTU{i},PSCoe_OTU{i},NSCoe_OTU{i}]=plotOTU(TaxSimplified,...
%         corr_rho_OTU{i}, sel_idx_OTU{i},0);
% %     a = Ptax_OTU{i};
%     Tax_corr_OTU{i} = union(Ptax_OTU{i},Ntax_OTU{i});
% %     idx = find(sel_idx_OTU{i}==1);
% %     temp = corr_rho_OTU{i};
% %     for j=1:length(idx)
% %         if  temp(idx(j))>0
% %             figure,plot(PathDis{i},Rel_OTU(idx(j),Path{i}),'o','markersize',8);
% %             title(TaxSimplified{idx(j)});
% %         end
% %     end
% end
% P1 = setdiff(Ptax_OTU{1},Ptax_OTU{2});
% P2 = setdiff(Ptax_OTU{2},Ptax_OTU{1});
% I=intersect(Ptax_OTU{1},Ptax_OTU{2});
%
% P3 = setdiff(Ptax_OTU{3},Ptax_OTU{4});
% P4 = setdiff(Ptax_OTU{4},Ptax_OTU{3});
% I3=intersect(Ptax_OTU{3},Ptax_OTU{4});
%
% PathDis = PathMatch(Path,PathDis);
% TotalDis = zeros(1,length(Label))*nan;
% TotalDis(Path{3}) = PathDis{3};
% TotalDis(Path{2}) = PathDis{2};
% TotalDis(Path{1}) = PathDis{1};

% %% progression diversity
% option.progression_diversity = 0;
% if option.progression_diversity == 1
%     num = max(Time);
%     diff_pat = cell(1,num);
%     U_pat = unique(PATIENT);
%     dis_pat_pool = [];
%         pat_pool = [];
%         dia_pool = [];
%         pat_l_pool = [];
%     for i=2:num
%         for j=1:length(U_pat)
%             temp_t = Time(PATIENT == U_pat(j));
%             if ~isempty(find(temp_t==1, 1)) && ~isempty(find(temp_t==i, 1))
%                 temp_dis = TotalDis(PATIENT == U_pat(j));
%                 pat_dia = DIA(PATIENT == U_pat(j));
%                 if pat_dia(1)~=1
%                     pidx = find(temp_t==i);
%                     for k=1:length(pidx)
%                         dis_pat_pool = [dis_pat_pool temp_dis(pidx(k))- temp_dis(temp_t==1)];
%                         pat_pool = [pat_pool U_pat(j)];
%                         dia_pool = [dia_pool pat_dia(1)];
%                         pat_l_pool = [pat_l_pool i];
%                     end
%                 end
%             end
%         end
%
%     end
%     num = max(pat_l_pool);
%     for i = 2:num
%          m_diff(i-1) = mean(dis_pat_pool(pat_l_pool==i));
%          s_diff(i-1) = std(dis_pat_pool(pat_l_pool==i));
%          m_diff_ccd(i-1) = mean(dis_pat_pool(pat_l_pool==i&dia_pool==2));
%          m_diff_icd(i-1) = mean(dis_pat_pool(pat_l_pool==i&dia_pool==3));
%          s_diff_ccd(i-1) = std(dis_pat_pool(pat_l_pool==i&dia_pool==2));
%          s_diff_icd(i-1) = std(dis_pat_pool(pat_l_pool==i&dia_pool==3));
%     end
%
%     figure,boxplot(dis_pat_pool(:),pat_l_pool(:));
%     figure,errorbar(2:max(pat_l_pool),m_diff,s_diff/2,'color','k');
%     hold on
%     errorbar(2:max(pat_l_pool),m_diff_ccd,s_diff_ccd/2,'color','b');
%     errorbar(2:max(pat_l_pool),m_diff_icd,s_diff_icd/2,'color','r');
%     figure,boxplot(TotalDis(DIA_SUB~=1),Time(DIA_SUB~=1));
%     figure,plot(Time(DIA_SUB~=1),TotalDis(DIA_SUB~=1),'o','markersize',8);
%     figure,plot(pat_l_pool(:)',dis_pat_pool(:)','o','markersize',8);
%
%     dis_pat_pool1 = dis_pat_pool(pat_l_pool<=5);
%     pat_l_pool1=pat_l_pool(pat_l_pool<=5);
%     pp=polyfit(pat_l_pool1',dis_pat_pool1',1);
%     xin = [unique(pat_l_pool1) 6:10];
%     yin=polyval(pp,xin);
%     figure,plot(xin,yin);hold on;plot(2:5,m_diff(1:4),'o','markersize',8)
% end
% corr_idx = sel_idx_OTU{1};
% for i=2:length(Path)
%     corr_idx = corr_idx|sel_idx_OTU{i};
% end
%
% s = idx20;
% s(idx20==1)=ids;
% CalHeatmap(Rel_OTU(s==1,:),TotalDis,cidx,TaxSimplified(s==1),Behavior,DIA_SUB)
% CalHeatmap(Rel_OTU(corr_idx==1,:),TotalDis,cidx,TaxSimplified(corr_idx==1),Behavior,DIA_SUB)
%%
% %% plot correlated speciese level
% TaxSimplified = SimplifyNameRefine(Tax_L7);
% for i=1:length(Path)
%     sel_idx_L7{i} = corr_q_L7{i}<0.01 &  corr_M_L7{i}>0.0001 & Fre_L7>0.30;
%     display(['# of selected OTUs along Path ' num2str(i) ': ' num2str(sum(sel_idx_L7{i}))]);
%     [Ptax_L7{i},Ntax_L7{i},PSCoe_L7{i},NSCoe_L7{i}]=plotOTU(TaxSimplified,...
%         corr_rho_L7{i}, sel_idx_L7{i});
%     Tax_corr_L7{i} = union(Ptax_L7{i},Ntax_L7{i});
% end
%
% corr_idx = sel_idx_L7{1};
% for i=2:length(Path)
%     corr_idx = corr_idx|sel_idx_L7{i};
% end
% CalHeatmap(Rel_L7(corr_idx==1,:),TotalDis,cidx,TaxSimplified(corr_idx==1),Behavior,DIA_SUB)

%%
a=mean( Rel_KOP(:,cidx>3),2);
b=mean( Rel_ECP(:,cidx>3),2);
option.cluster_comparision = 1;
if option.cluster_comparision == 1
     %% KOP
    [~,pair_PVAL_KOP,pair_result_KOP] = dsfdrN(...
        Rel_KOP(:,cidx>3),cidx(cidx>3),'ranksum',0.05,'dsfdr');
    save('IBD_cluster4vs5_KOP','pair_PVAL_KOP','pair_result_KOP');
    
    %% pairwise cluster analysis
    
    %% species level
    [~,pair_PVAL_L7,pair_result_L7] = dsfdrN(Rel_L7(:,cidx>3),cidx(cidx>3),...
        'ranksum',0.05,'dsfdr');
    save('IBD_cluster4vs5_L7','pair_PVAL_L7','pair_result_L7');
    
    %% OTUs
    [~,pair_PVAL_OTU,pair_result_OTU] = dsfdrN(Rel_OTU(:,cidx>3),cidx(cidx>3),...
        'ranksum',0.05,'dsfdr');
    save('IBD_cluster4vs5_otu','pair_PVAL_OTU','pair_result_OTU');
    
    %% KO
    [~,pair_PVAL_KO,pair_result_KO] = dsfdrN(...
        Rel_KO(:,cidx>3),cidx(cidx>3),'ranksum',0.05,'dsfdr');
    save('IBD_cluster4vs5_KO','pair_PVAL_KO','pair_result_KO');
    %% EC
    [~,pair_PVAL_EC,pair_result_EC] = dsfdrN(...
        Rel_EC(:,cidx>3),cidx(cidx>3),'ranksum',0.05,'dsfdr');
    save('IBD_cluster4vs5_EC','pair_PVAL_EC','pair_result_EC');
    %% ECP
    [~,pair_PVAL_ECP,pair_result_ECP] = dsfdrN(...
        Rel_ECP(:,cidx>3),cidx(cidx>3),'ranksum',0.05,'dsfdr');
    save('IBD_cluster4vs5_ECP','pair_PVAL_ECP','pair_result_ECP');
end
load('IBD_cluster4vs5_otu');
LogFC_OTU = CalFC(Rel_OTU(:,cidx>3),cidx(cidx>3));
LogFC_OTU(LogFC_OTU>10)=10;
TaxSimplified = SimplifyNameRefine(Tax_OTU);
sel_idx_OTU = pair_result_OTU.Q<0.05& mean(Rel_OTU(:,cidx>3),2)>0.0005 & abs(LogFC_OTU)>0.5850;
sel_idx_OTU = pair_result_OTU.Q<0.05& M_OTU>0.0001 & abs(LogFC_OTU)>0.5850;
display(['# of selected OTUs: ' num2str(sum(sel_idx_OTU))]);
[Ptax_com_OTU,Ntax_com_OTU,PSCoe_com_OTU,NSCoe_com_OTU]=plotOTU(TaxSimplified,...
    LogFC_OTU, sel_idx_OTU,1);
Tax_com_OTU = union(Ptax_com_OTU,Ntax_com_OTU);
Com_P = setdiff(Ptax_com_OTU,Ntax_com_OTU);
Com_N = setdiff(Ntax_com_OTU,Ptax_com_OTU);
option.alpha = 0;
if option.alpha == 1
    %% alpha diversity analysis
    % calculate the alpha diversity without rarefaction
    para4alpha.rare=0;
    %     [observed_otu, shannon_index, Chao1_index] = CalAlphaDiversity(TABLE_OTU,para4alpha);
    %     save('SwedishAlphaDiversity_direct','observed_otu', ...
    %         'shannon_index', 'Chao1_index', 'para4alpha');
    
    % calculate the alpha diversity with proportionality rarefaction
    para4alpha.rare=1;
    para4alpha.num=10000;
    %     [observed_otu, shannon_index, Chao1_index] = CalAlphaDiversity(TABLE_OTU,para4alpha);
    %     save('SwedishAlphaDiversity_prop_rare','observed_otu', ...
    %         'shannon_index', 'Chao1_index', 'para4alpha');
    
    % calculate the alpha diversity with random rarefaction
    para4alpha.rare=2;
    para4alpha.num=10000;
    %     [observed_otu, shannon_index, Chao1_index] = CalAlphaDiversity(TABLE_OTU,para4alpha);
    %     save('SwedishAlphaDiversity_random_rare','observed_otu', ...
    %         'shannon_index', 'Chao1_index', 'para4alpha');
end


option.calprotecin = 0;
if option.calprotecin == 1
    calP_no_nan = calprotectin(~isnan(calprotectin));
    cidx_no_nan = cidx(~isnan(calprotectin));
    calP_C4 = calP_no_nan(cidx_no_nan==4);
    calP_C5 = calP_no_nan(cidx_no_nan==5);
    calP_C5_quan = quantile(calP_C4,[.05 .50 .95]);
    calP_C5_quan = quantile(calP_C5,[.05 .50 .95]);
    calP_p_ranksum = ranksum(calP_C4, calP_C5);
    
    calP_select = calP_no_nan(cidx_no_nan>3);
    cidx_select = cidx_no_nan(cidx_no_nan>3);
    
    [table_calP,chi2_calP,p_calP] = crosstab(calP_select > 150, cidx_select);
    
end

option.bmi = 0;
if option.bmi == 1
    bmi_no_nan = bmi(~isnan(bmi));
    cidx_no_nan = cidx(~isnan(bmi));
    
    bmi_C4 = bmi_no_nan(cidx_no_nan==4);
    bmi_C5 = bmi_no_nan(cidx_no_nan==5);
    
    bmi_C4_quan = quantile(bmi_C4,[.05 .50 .95]);
    bmi_C5_quan = quantile(bmi_C5,[.05 .50 .95]);
    [~,p_bmi] = ttest2(bmi_C4,bmi_C5);
end

option.duration = 0;
if option.duration == 1
    duration_no_nan = dur(~isnan(dur));
    cidx_no_nan = cidx(~isnan(dur));
    
    duration_C4 = duration_no_nan(cidx_no_nan==4);
    duration_C5 = duration_no_nan(cidx_no_nan==5);
    
    duration_C4_quan = quantile(duration_C4,[.05 .50 .95]);
    duration_C5_quan = quantile(duration_C5,[.05 .50 .95]);
    
    [~,p_duration] = ttest2(duration_C4,duration_C5);
end

option.sex = 0;
if option.sex == 1
    sex_select = sex(cidx>3);
    cidx_select = cidx(cidx>3);
    [table_sex, chi2_sex, p_sex] = crosstab(sex_select,cidx_select);
    
    % check patient
    
    sex_C4 = sex(cidx==4);
    Pat_C4 = PATIENT(cidx==4);
    [~,ia]=unique(Pat_C4);
    sex_C4 = sex_C4(ia);
    
    sex_C5 = sex(cidx==5);
    Pat_C5 = PATIENT(cidx==5);
    [~,ia]=unique(Pat_C5);
    sex_C5 = sex_C5(ia);
    
    [table_sex_pat, chi2_sex_pat, p_sex_pat] = crosstab([sex_C4;sex_C5],...
        [ones(size(sex_C4))*4;ones(size(sex_C5))*5]);
end

option.perianal = 0;
if option.perianal == 1
    % 1:no; 2:yes
    perianal_select = perianal(cidx>3);
    cidx_select = cidx(cidx>3);
    [table_Per,chi2_Per,p_Per] = crosstab(perianal_select,cidx_select);
    
    % check patient
    
    perianal_C4 = perianal(cidx==4);
    Pat_C4 = PATIENT(cidx==4);
    [~,ia]=unique(Pat_C4);
    perianal_C4 = perianal_C4(ia);
    
    perianal_C5 = perianal(cidx==5);
    Pat_C5 = PATIENT(cidx==5);
    [~,ia]=unique(Pat_C5);
    perianal_C5 = perianal_C5(ia);
    
    [table_Per_pat,chi2_Per_pat,p_Per_pat] = crosstab([perianal_C4;perianal_C5],...
        [ones(size(perianal_C4))*4;ones(size(perianal_C5))*5]);
end

end

function [model,history] = principal_tree(X,params)

[D, N] = size(X);

% initialize MU, lmabda,
if ~isfield(params,'MU')
    K = N;
    MU = X;
else
    MU = params.MU;
    K = size(MU,2);
end

lambda = params.lambda * N; % scale the parameter by N

if ~isfield(params,'bandwidth')
    distsqX = sqdist(X,X);
    sigma = 0.1 * sum(sum(distsqX))/ (N*N);
else
    sigma = params.bandwidth;
end

for iter=1:params.maxIter
    
    % Kruskal method to find a spanning tree
    distsqMU = sqdist(MU,MU);
    stree = graphminspantree(sparse(tril(distsqMU)),'Method','Kruskal');
    stree = stree + stree';
    e = stree ~= 0;
    
    history.mu{iter} = MU;
    history.stree{iter} = stree;
    
    % update data assignment matrix
    distMUX = sqdist(MU,X);
    min_dist = repmat(min(distMUX,[],1),K,1);
    tmp_distMUX = distMUX - min_dist;
    tmp_R = exp(-tmp_distMUX' ./ sigma);
    
    %     for i=1:size(tmp_R,1)
    %         temp = tmp_R(i,:);
    %         [~,idx] = sort(temp,'descend');
    %         tmp_R(i,idx(51:end))=0;
    %     end
    
    R = tmp_R ./ repmat( sum(tmp_R,2), 1, K );
    
    % compute objective function
    obj1 = - sigma * sum( log( sum( exp(- tmp_distMUX./sigma) ,1) ) ...
        - min_dist(1,:) ./ sigma );
    reg = sum(sum(full(stree)));
    obj = (obj1 + 0.5 * lambda * reg)/N;
    history.objs(iter) = obj;
    
    % projected mean square error
    %     mse = sum(sum( R .* distMUX'));
    %     mse = obj1 + sigma *N* log(K);
    
    projd = min(distMUX,[],1);
    mse = mean(projd);
    
    history.mse(iter) = mse;
    
    % length of the structure
    history.length(iter)=reg;
    
    % terminate condition
    fprintf('iter=%d, obj=%f, mse=%f, len=%f\n',iter,obj,mse,reg);
    if iter > 1
        if abs((obj - old_obj)/old_obj) < 1e-3
            break;
        end
    end
    
    % compute the means
    L = diag(sum(e,1)) - e;
    MU = X * R / ( lambda * L + diag(sum(R,1)) );
    
    old_obj = obj;
end

model.MU = MU;
model.stree = stree;
model.sigma = sigma;
model.lambda = params.lambda;
end
function [R,P,Q,M,C]=PathTest(extracted_path,extracted_pathDist,D,Type,ph)
% Select OTUs which are significantly correlated with the progression path
% in terms of the Spearman correlation test.
% Input
% extracted_path: extracted path
% extracted_pathDist: progression distance of extracted path
% D: relative abudance
% Type: multiple testing correction
%        S: estimate the FDR using Storey's method
%        BH1: Benjamini and Hochberg implemented in matlab
%        BH2: Benjamini and Hochberg implemented by others
% ph: FDR level
% Output
% R: Correlation coefficient
% P: Correlation p value
% Q: Adjusted p value
% M: Average abundance
if nargin <5
    ph=0.05;
end
Q=ones(size(D,1),1);
if strcmp(Type,'dsfdr')~=1
    for kk=1
        for jj=1:size(D,1)
            ty=D(jj,extracted_path); % Relative abundance
            tx=extracted_pathDist;   % Progression distance
            %         [RHO,PVAL] = corr(tx',ty','type','Pearson');
            [RHO,PVAL] = corr(tx',ty','type','Spearman');
            
            M(jj,kk)=mean(ty);
            C(jj,kk)=length(find(ty>0))/length(ty);
            if C(jj,kk)>0
                R(jj,kk)=RHO;
                P(jj,kk)=PVAL;
            else
                R(jj,kk)=0;
                P(jj,kk)=1;
            end
        end
        switch Type
            case 'S'
                adjust_p=ones(size(C(:,kk)));
                ids=C(:,kk)>0;
                adjust_p(ids==1)=mafdr(P(ids==1,kk));
            case 'BH1'
                adjust_p=ones(size(C(:,kk)));
                ids=C(:,kk)>0;
                adjust_p(ids==1)=mafdr(P(ids==1,kk),'BHFDR',true);
            case 'BH2'
                adjust_p=ones(size(C(:,kk)));
                ids=C(:,kk)>0;
                [~, ~, ~, adjust_p(ids==1)]=fdr_bh(P(ids==1,kk),ph,'pdep','yes');
        end
        Q(:,kk)=adjust_p;
    end
else
    for kk=1
        for jj=1:size(D,1)
            ty=D(jj,extracted_path); % Relative abundance
            C(jj,kk)=length(find(ty>0))/length(ty);
        end
        ty=D(C>0,extracted_path); % Relative abundance
        tx=extracted_pathDist;   % Progression distance
        %         [RHO,PVAL] = corr(tx',ty','type','Pearson');
        %             [RHO,PVAL] = corr(tx',ty','type','Spearman');
        [~,PVAL,result]=dsfdrN(ty,tx,'spearman',0.05,'dsfdr');
        R=zeros(size(C));
        P=ones(size(C));
        M=zeros(size(C));
        R(C>0,kk)=result.rho;
        P(C>0,kk)=result.pvals;
        M(C>0,kk)=mean(ty,2);
        Q(C>0,kk)=result.Q;
    end
end
end
function [reject,pvals,result]=dsfdrN(data,labels,method,alpha,fdr_method)
%%
% calculate the Discrete FDR for the data
%
%     input:
%     data : N x S numpy array
%         each column is a sample (S total), each row an OTU (N total)
%     labels : a 1d numpy array (length S)
%         the labels of each sample (same order as data) with the group
%         (0/1 if binary, 0-G-1 if G groups, or numeric values for correlation)
%
%
%     transform_type : str or None
%         transformation to apply to the data before caluculating
%         the test statistic
%         'rankdata' : rank transfrom each OTU reads
%         'log2data' : calculate log2 for each OTU using minimal cutoff of 2
%         'normdata' : normalize the data to constant sum per samples
%         'binarydata' : convert to binary absence/presence
%         'clrdata' : clr transformation of data (after replacing 0 with 1)
%          None : no transformation to perform
%
%     method : str or function
%         the method to use for calculating test statistics:
%         'meandiff' : mean(A)-mean(B) (binary)
%         'mannwhitney' : mann-whitney u-test (binary)
%         'kruwallis' : kruskal-wallis test (multiple groups)
%         'stdmeandiff' : (mean(A)-mean(B))/(std(A)+std(B)) (binary)
%         'spearman' : spearman correlation (numeric)
%         'pearson' : pearson correlation (numeric)
%         'nonzerospearman' : spearman correlation only non-zero entries
%                             (numeric)
%         'nonzeropearson' : pearson correlation only non-zero entries (numeric)
%         function : use this function to calculate the test statistic
%         (input is data,labels, output is array of float)
%
%     alpha : float
%         the desired FDR control level
%     numperm : int
%         number of permutations to perform
%
%     fdr_method : str
%         the FDR procedure to determine significant bacteria
%         'dsfdr' : discrete FDR method
%         'bhfdr' : Benjamini-Hochberg FDR method
%         'byfdr' : Benjamini-Yekutielli FDR method
%         'filterBH' : Benjamini-Hochberg FDR method with filtering
%         'gilbertBH' : Benjamini-Hochberg FDR method with Gilbert (2005) pre-filtering
%
%     output:
%     reject : np array of bool (length N)
%         True for OTUs where the null hypothesis is rejected
%     tstat : np array of float (length N)
%         the test statistic value for each OTU (for effect size)
%     pvals : np array of float (length N)
%         the p-value for each OTU
%%
rng(100)
numperm=1000;
idx_s=zeros(size(data,1),1);
for i=1:size(data,1)
    if length(unique(data(i,:)))==1
        idx_s(i)=1;
    end
end
data=data(idx_s==0,:);
orig_numbact=size(data,1);
filtered_order =1:orig_numbact;
if strcmp(fdr_method,'filterBH')==1
    index=[];
    U=unique(labels);
    temp=labels;
    labels(temp==U(1))=0;
    labels(temp==U(2))=1;
    n0=sum(labels==0);
    n1=sum(labels==1);
    for i=1:orig_numbact
        nonzeros = length(find(data(i, :)~=0));
        if nonzeros<min(n0,n1)
            pval_min =(length(combnk(n0, nonzeros))+length(combnk(n1, nonzeros)))...
                /length(combnk(n0+n1, nonzeros));
            if pval_min<=alpha
                index=[index i];
            end
        else
            index=[index i];
        end
    end
    data=data(index,:);
    filtered_order=filtered_order(index);
else
    if strcmp(fdr_method,'gilbertBH')==1
        % caluclate the Gilbert alpha* per feature (minimal ibtainable p-value)
        alpha_star = [];
        U=unique(labels);
        temp=labels;
        labels(temp==U(1))=0;
        labels(temp==U(2))=1;
        n0=sum(labels==0);
        n1=sum(labels==1);
        for i=1:orig_numbact
            % test if all values are identical, max p-val is 1 (need to filter)
            if length(unique(data(i,:)))==1
                alpha_star=[alpha_star 1];
                continue
            end
            cdat = sort(data(i,:)); % sort in acending order
            rdat = sort(data(i,:),'descend'); % sort in decending order
            
            [p1,tbl1,stats1] = kruskalwallis(cdat(:)',[zeros(1,n0) ones(1,n1)],'off');
            [p2,tbl2,stats2] = kruskalwallis(rdat(:)',[zeros(1,n0) ones(1,n1)],'off');
            s1=tbl1{2,5};
            s2=tbl2{2,5};
            alpha_star=[alpha_star min(p1,p2)];
            %find the smallest K which is big enough for Bonferoni (that's how it's done in Gilbert)
            for ck = 1:orig_numbact
                num_ok = sum(alpha_star < alpha / ck);
                if num_ok <= ck
                    break
                end
                
                
            end
            %and keep only the features which match it
            index = (alpha_star < alpha / ck);
            data = data(index, :);
            filtered_order = filtered_order(index);
        end
        
    end
end

% transform the data

if isempty(alpha)
    alpha=0.05;
end

if isempty(fdr_method)
    fdr_method='dsfdr';
end

if isempty(method)
    method='meandiff';
end
[numbact,S]=size(data);

labels=labels(:);
switch method
    case 'kruwallis'
        tstat=StaTest(data,labels,method);
        t=abs(tstat);
        u=zeros(numbact,numperm);
        for i=1:numperm
            display(['Perm: ' num2str(i)]);
            rlabels=labels(randperm(S));
            rt=StaTest(data,rlabels,method);
            u(:,i)=rt;
        end
    case 'ranksum'
        tstat=StaTest(data,labels,method);
        t=abs(tstat);
        u=zeros(numbact,numperm);
        for i=1:numperm
            display(['Perm: ' num2str(i)]);
            rlabels=labels(randperm(S));
            rt=StaTest(data,rlabels,method);
            u(:,i)=rt;
        end
    case 'ttest'
        tstat=StaTest(data,labels,method);
        t=abs(tstat);
        u=zeros(numbact,numperm);
        for i=1:numperm
            display(['Perm: ' num2str(i)]);
            rlabels=labels(randperm(S));
            rt=StaTest(data,rlabels,method);
            u(:,i)=rt;
        end
        
    case 'spearman'
        [tstat,rho]=StaTest(data,labels,method);
        %         data=para_spearman.x;
        %         labels=para_spearman.y;
        t = abs(tstat);
        permlabels = zeros([length(labels), numperm]);
        for cperm =1: numperm
            display(['Perm: ' num2str(cperm)]);
            rnd_num = randperm(length(labels));
            rlabels = labels(rnd_num);
            permlabels(:, cperm) = rlabels(:);
        end
        %         u = abs(data, permlabels);
        u=StaTest(data,permlabels,method);
        result.rho=zeros(size(idx_s))*nan;
        result.rho(idx_s==0)=rho;
end

pvals=tstat;
pvals_u=u;



% calculate FDR
% pvals=round(pvals/10^(-8))*10^(-8);
pvals_unique = unique(pvals);
sortp = sort(pvals_unique,'descend');

% find a data-dependent threshold for the p-value
foundit = 0;
allfdr = [];
allt = [];
realcp=2;
for i=1:length(sortp)
    cp=sortp(i);
    realnum = sum(pvals <= cp);
    fdr = (realnum + length(find(pvals_u <= cp))) / ...
        (realnum * (numperm + 1));
    allfdr=[allfdr fdr];
    allt=[allt cp];
    if fdr <alpha&&foundit ==0
        realcp = cp;
        foundit = 1;
    end
end

if foundit ==0
    reject=zeros(numbact,1);
else
    reject=pvals <= realcp;
end
adjp=pvals;
for i=1:length(allt)
    adjp(pvals==allt(i))=allfdr(i);
end
result.tstat=zeros(size(idx_s))*nan;
result.tstat(idx_s==0)=tstat;
result.s=idx_s;

result.pvals=zeros(size(idx_s))*nan;
result.pvals(idx_s==0)=pvals;

result.pvals_u=zeros(size(idx_s,1),numperm)*nan;
result.pvals_u(idx_s==0,:)=pvals_u;

result.realcp=realcp;
result.allfdr=allfdr;
result.allt=allt;
result.Q=zeros(size(idx_s))*nan;
result.Q(idx_s==0)=adjp;
end

function [tstat,para]=StaTest(X,Y,method)
tstat=zeros(size(X,1),size(Y,2));
para=[];
U=unique(Y);
switch method
    case 'kruwallis'
        for i=1:size(X,1)
            [p,tbl,stats] = kruskalwallis(X(i,:),Y,'off');
            %             p=ranksum(X(i,Y==U(1)),X(i,Y==U(2)));
            tstat(i)=p;%tbl{2,5};
        end
    case 'ranksum'
        for i=1:size(X,1)
            p=ranksum(X(i,Y==U(1)),X(i,Y==U(2)));
            tstat(i)=p;%tbl{2,5};
        end
    case 'ttest'
        U=unique(Y);
        for i=1:size(X,1)
            [~,p] = ttest2(X(i,Y==U(1)),X(i,Y==U(2)));
            tstat(i)=p;%tbl{2,5};
        end
    case 'spearman'
        %         Y=transform_rankdata(Y);
        %         Y=Y-mean(Y);
        %         X=transform_rankdata(X);
        %         meanval=mean(X,2);
        %         X=X-repmat(meanval,1,size(X,2));
        %         tstat=X*Y(:);
        %         para.x=X;
        %         para.x=Y;
        [para,tstat]=corr(X',Y,'type','Spearman');
        
    case 'nonzerospearman'
        for i=1:size(X,1)
            idx0=find(X(i,:)~=0);
            Xt=X(i,idx0);
            Yt=Yt(idx0);
            
            Xt=transform_rankdata(Xt);
            Yt=transform_rankdata(Yt);
            Yt=Yt-mean(Yt);
            
            meanval=mean(Xt);
            Xt=Xt-meanval;
            tstat(i)=Xt*Yt(:);
        end
end
end