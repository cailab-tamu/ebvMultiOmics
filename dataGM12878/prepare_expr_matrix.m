%[X1,genelist1,barcodelist]=sc_readmtxfile('matrixRNA.mtx','genesRNA.txt','barcodesRNA.txt',1);
%[X2,genelist2]=sc_readmtxfile('ebvGM12878_matrix.mtx','ebvGM12878_geneList.txt',[],1);
X=[X1; X2];
genelist=[genelist1;genelist2];

%[Xnorm]=sc_norm(X1,'type','deseq');            % normalize using DeSeq method
Xnorm=sc_transform(X1);
[s_tsne]=sc_tsne(Xnorm,3,false,false);
%[T,Xhvg]=sc_hvg(Xnorm,genelist1,true);         % identify highly variable genes (HVGs)
%[s_tsne]=sc_tsne(Xhvg(1:1000,:),3,false,false);   % using expression of top 2000 HVGs for tSNE
sce=SingleCellExperiment(X,genelist,s_tsne);      % make SCE class
% sce=sce.estimatepotency(2);                   % estimate differentiation potency (1-human; 2-mouse)
% sce=sce.estimatecellcycle;                  % estimate cell cycle phase using R/Seurat
%id=sc_cluster_s(s_tsne,10);                   % clustering on tSNE coordinates using k-means
%sce.c_cluster_id=id;                          % assigning cluster Ids to SCE class
sc_scatter(sce)         