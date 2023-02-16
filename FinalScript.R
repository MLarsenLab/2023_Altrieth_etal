#Version 2
#Uses Seurat packages for analysis: https://satijalab.org/seurat/index.html
#This version is using data that is deconvoluted using CellRanger
#Includes function for finding number/percentage of cells expressing a specific gene
#Major parameter and statistical decisions made: 
###reprocesses individual samples to have % mt gene removed and have all samples processed with the same settings
###integrates mock sample from 2/15/22 with ligated sample from 11/30/21
###integrates mock samples from 2/15/22 and 6/29/21 with ligated samples from 6/30/21 and 11/30/21
##subclusters endothelial cells for mock 2/15/22 and ligated 11/30/21
##subclusters endothelial cells for mock 2/15/22 and 6/29/21 and ligated 11/30/21 and 6/30/21

###########################################SESSION INFO##############################
#R version 3.6.3 (2020-02-29)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Debian GNU/Linux 10 (buster)

#Matrix products: default
#BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so

#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
#[4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
# [1] metap_1.3           pheatmap_1.0.12     BiocManager_1.30.10 stringr_1.4.0      
#[5] ggplot2_3.3.0       patchwork_1.0.0     dplyr_0.8.5         cowplot_1.0.0      
#[9] Seurat_3.1.5       

#loaded via a namespace (and not attached):
# [1] nlme_3.1-144        tsne_0.1-3          bitops_1.0-6        RcppAnnoy_0.0.16   
#[5] RColorBrewer_1.1-2  httr_1.4.1          numDeriv_2016.8-1.1 sctransform_0.2.1  
#[9] tools_3.6.3         R6_2.4.1            irlba_2.3.3         KernSmooth_2.23-16 
#[13] BiocGenerics_0.32.0 sn_1.6-1            uwot_0.1.8          lazyeval_0.2.2     
#[17] colorspace_1.4-1    npsurv_0.4-0        withr_2.2.0         mnormt_1.5-6       
#[21] tidyselect_1.0.0    gridExtra_2.3       compiler_3.6.3      Biobase_2.46.0     
#[25] TFisher_0.2.0       sandwich_2.5-1      plotly_4.9.2.1      labeling_0.3       
#[29] caTools_1.18.0      scales_1.1.0        mvtnorm_1.1-0       lmtest_0.9-37      
#[33] ggridges_0.5.2      pbapply_1.4-2       rappdirs_0.3.1      digest_0.6.25      
#[37] pkgconfig_2.0.3     htmltools_0.4.0     bibtex_0.4.2.2      plotrix_3.7-8      
#[41] htmlwidgets_1.5.1   rlang_0.4.5         rstudioapi_0.11     farver_2.0.3       
#[45] zoo_1.8-7           jsonlite_1.6.1      ica_1.0-2           gtools_3.8.2       
#[49] magrittr_1.5        Matrix_1.2-18       Rcpp_1.0.4.6        munsell_0.5.0      
#[53] ape_5.3             reticulate_1.15     lifecycle_0.2.0     multcomp_1.4-13    
#[57] stringi_1.4.6       gbRd_0.4-11         MASS_7.3-51.5       gplots_3.0.3       
#[61] Rtsne_0.15          plyr_1.8.6          grid_3.6.3          parallel_3.6.3     
#[65] gdata_2.18.0        listenv_0.8.0       ggrepel_0.8.2       crayon_1.3.4       
#[69] lattice_0.20-38     splines_3.6.3       multtest_2.42.0     pillar_1.4.3       
#[73] igraph_1.2.5        stats4_3.6.3        future.apply_1.5.0  reshape2_1.4.4     
#[77] codetools_0.2-16    leiden_0.3.3        mutoss_0.1-12       glue_1.4.0         
#[81] lsei_1.2-0          data.table_1.12.8   png_0.1-7           vctrs_0.2.4        
#[85] Rdpack_0.11-1       gtable_0.3.0        RANN_2.6.1          purrr_0.3.4        
#[89] tidyr_1.0.2         future_1.17.0       assertthat_0.2.1    rsvd_1.0.3         
#[93] RSpectra_0.16-0     survival_3.1-8      viridisLite_0.3.0   tibble_3.0.1       
#[97] cluster_2.1.0       globals_0.12.5      TH.data_1.0-10      fitdistrplus_1.0-14
#[101] ellipsis_0.3.0      ROCR_1.0-7   

###############################Functions#####################################
#Fucntion for quickly making UMAP and violin plots for all genes from a gene list
#Does not split by protocol! Loop to create UMAPS, Violin plots, and Expression Levels per cluster for GenesLists, can pass any gene-list to generate plots

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(stringr)
library(magrittr)

setwd("Insert_Path")

##Initial QC Analysis:
##Reads the 10x Generated data
agg.mock215.data <- Read10X(data.dir = "/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_021522_Mock_Repeat/ML_2_Mock_021522/outs/filtered_feature_bc_matrix") 
agg.lig1130.data <- Read10X(data.dir = "/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_216_Ligated_Repeat_02/ML_1_ligated_113021/outs/filtered_feature_bc_matrix")

##Creates Seurat object of the 10x data, looks for a minimum of 3 cells with at least 200 features(AKA genes)
agg.lig1130 <- CreateSeuratObject(agg.lig1130.data, project = "agg.lig1130", min.cells = 3, min.features = 200) 
agg.mock215 <- CreateSeuratObject(agg.mock215.data, project = "agg.mock215", min.cells = 3, min.features = 200) 

##Creates new parameter selecting for mitochondrial genes
agg.lig1130[["percent.mt"]] <- PercentageFeatureSet(agg.lig1130, pattern = "^mt-")
agg.mock215[["percent.mt"]] <- PercentageFeatureSet(agg.mock215, pattern = "^mt-") 

#generates 3 different violin plots based on the Features, Count RNA, and percent mitochrondrial genes, each dot is a cell
VlnPlot(agg.lig1130, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(agg.mock215, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

##Filtering of Data and Scaling:
#filters out based on at least 200 features but less than 9000 and less than 5% mitochondrialDNA
agg.lig1130 <- subset(agg.lig1130, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)
agg.mock215 <- subset(agg.mock215, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)

#normalizes the data, allows for accurate comparison of clusters and to be able to perform differential gene expression between clusters
agg.lig1130<- NormalizeData(agg.lig1130) 
agg.mock215<- NormalizeData(agg.mock215) 

#identifies the most highly variable genes and returns the top 2000 by default between the different clusters     
agg.lig1130 <- FindVariableFeatures(agg.lig1130, selection.method = "vst", nfeatures = 2000) 
agg.mock215 <- FindVariableFeatures(agg.mock215, selection.method = "vst", nfeatures = 2000) 

#apply a linear transformation that is a standard pre-processing step prior to dimensional reduction techniques like PCA, standard preprocessing step for UMAP and tSNA plots
agg.lig1130 <- ScaleData(agg.lig1130, feature = all.geneslig1130) 
agg.mock215 <- ScaleData(agg.mock215, feature = all.genesmock215) 

##Clustering:
#Construct a graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity), and takes as input the previously defined dimensionality of the dataset (first 40 PCs), want to choose based on your elbow plot!    
agg.lig1130 <- FindNeighbors(agg.lig1130, dims = 1:21) 
agg.mock215 <- FindNeighbors(agg.mock215, dims = 1:21) 

#Generate the clusters, resolution increases the specificity of each cluster (assigns individual cells to clusters based on sim/diff)
agg.lig1130<- FindClusters(agg.lig1130, resolution = 1)
agg.mock215<- FindClusters(agg.mock215, resolution = 1)

#Generate a UMAP dataset, use the same dimensionality you used for FindNeighborsGenerate a UMAP dataset, use the same dimensionality you used for FindNeighbors
agg.lig1130 <- RunUMAP(agg.lig1130, dims = 1:21) 
agg.mock215 <- RunUMAP(agg.mock215, dims = 1:21) 

#Saves Formated dataset after all major computational processing
saveRDS(agg.lig1130, file = "Lig1130_(SEURAT_v3)_02.rds")
saveRDS(agg.mock215, file = "mock215_(SEURAT_v3)_02.rds")

#SessionInfo
#R version 4.1.2 (2021-11-01)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19043)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] magrittr_2.0.3     stringr_1.4.1      clustree_0.5.0     ggraph_2.0.6       dplyr_1.0.9       
#[6] patchwork_1.1.2    ggplot2_3.3.6      sp_1.5-0           SeuratObject_4.1.0 Seurat_4.1.1      

#loaded via a namespace (and not attached):
#  [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6          ellipsis_0.3.2       
#[5] ggridges_0.5.3        spatstat.data_2.2-0   leiden_0.4.2          listenv_0.8.0        
#[9] farver_2.1.0          graphlayouts_0.8.1    ggrepel_0.9.1         fansi_1.0.3          
#[13] codetools_0.2-18      splines_4.1.2         polyclip_1.10-0       jsonlite_1.8.0       
#[17] ica_1.0-3             cluster_2.1.2         png_0.1-7             rgeos_0.5-9          
#[21] uwot_0.1.11           ggforce_0.3.4         shiny_1.7.2           sctransform_0.3.3    
#[25] spatstat.sparse_2.1-1 compiler_4.1.2        httr_1.4.3            assertthat_0.2.1     
#[29] Matrix_1.3-4          fastmap_1.1.0         lazyeval_0.2.2        cli_3.3.0            
#[33] later_1.3.0           tweenr_2.0.1          htmltools_0.5.2       tools_4.1.2          
#[37] igraph_1.3.4          gtable_0.3.0          glue_1.6.2            RANN_2.6.1           
#[41] reshape2_1.4.4        Rcpp_1.0.8.3          scattermore_0.8       vctrs_0.4.1          
#[45] nlme_3.1-153          progressr_0.10.1      lmtest_0.9-40         spatstat.random_2.2-0
#[49] globals_0.15.1        mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.1      
#[53] irlba_2.3.5           goftest_1.2-3         future_1.26.1         MASS_7.3-54          
#[57] zoo_1.8-10            scales_1.2.0          tidygraph_1.2.2       spatstat.core_2.4-4  
#[61] promises_1.2.0.1      spatstat.utils_2.3-1  parallel_4.1.2        RColorBrewer_1.1-3   
#[65] reticulate_1.25       pbapply_1.5-0         gridExtra_2.3         rpart_4.1-15         
#[69] stringi_1.7.6         rlang_1.0.2           pkgconfig_2.0.3       matrixStats_0.62.0   
#[73] lattice_0.20-45       ROCR_1.0-11           purrr_0.3.4           tensor_1.5           
#[77] htmlwidgets_1.5.4     cowplot_1.1.1         tidyselect_1.1.2      parallelly_1.32.1    
#[81] RcppAnnoy_0.0.19      plyr_1.8.7            R6_2.5.1              generics_0.1.3       
#[85] DBI_1.1.3             pillar_1.7.0          withr_2.5.0           mgcv_1.8-38          
#[89] fitdistrplus_1.1-8    survival_3.2-13       abind_1.4-5           tibble_3.1.7         
#[93] future.apply_1.9.0    crayon_1.5.1          KernSmooth_2.23-20    utf8_1.2.2           
#[97] spatstat.geom_2.4-0   plotly_4.10.0         viridis_0.6.2         grid_4.1.2           
#[101] data.table_1.14.2     digest_0.6.29         xtable_1.8-4          tidyr_1.2.0          
#[105] httpuv_1.6.5          munsell_0.5.0         viridisLite_0.4.0

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(stringr)
library(magrittr)

genes = c('Abi3bp','Aebp1','Cilp','Creld2','Crispld1','Efemp2','Egflam','Eln','Emilin1','Fbn2',
          'Gas6','Igfbp7','Igfbp6','Ltbp2','Mfap2','Mfap4','Mfap5','Mgp','Postn',
          'Slit2','Sparc','Spp1','Svep1','Thbs2','Thbs4','Col11a1','Col12a1',
          'Col16a1','Col28a1','Col4a3','Col4a4','Col4a5','Col8a1','Col8a2','Bgn',
          'Fmod','Omd','Dpp4')

#

Matrisome = c('5430419D17Rik','Abi3bp','Adipoq','Aebp1','Agrn',
              'Ambn','Amelx','AW551984','Bglap2','Bglap3','Bmper','Bsph1','Bsph2',
              'Cdcp2','Cilp','Cilp2','Coch','Colq','Comp','Creld1','Creld2','Crim1',
              'Crispld1','Crispld2','Ctgf','Cthrc1','Cyr61','Ddx26b','Dmbt1','Dmp1','Dpt','Dspp',
              'Ecm1','Ecm2','Edil3','Efemp1','Efemp2','Egfem1','Egflam','Eln','Emid1','Emilin1','Emilin2','Emilin3',
              'Fbln1','Fbln2','Fbln5','Fbln7','Fbn1','Fbn2','Fga','Fgb','Fgg','Fgl1',
              'Fgl2','Fn1','Fndc1','Fndc7','Fndc8','Fras1','Gas6','Gldn','Hmcn1','Hmcn2',
              'Ibsp','Igfals','Igfbp1','Igfbp2','Igfbp3','Igfbp4','Igfbp5','Igfbp6','Igfbp7',
              'Igfbpl1','Igsf10','Kcp','Lama1','Lama2','Lama3','Lama4','Lama5','Lamb1',
              'Lamb2','Lamb3','Lamc1','Lamc2','Lamc3','Lgi1','Lgi2','Lgi3','Lgi4','Lrg1',
              'Ltbp1','Ltbp2','Ltbp3','Ltbp4','Matn1','Matn2','Matn3','Matn4','Mepe','Mfap1a',
              'Mfap1b','Mfap2','Mfap3','Mfap4','Mfap5','Mfge8','Mgp','Mmrn1','Mmrn2','Ndnf',
              'Nell1','Nell2','Nid1','Nid2','Nov','Npnt','Ntn1','Ntn3','Ntn4','Ntn5','Ntng1',
              'Ntng2','Oit3','Otog','Otogl','Otol1','Papln','Pcolce','Pcolce2','Postn','Pxdn',
              'Reln','Rspo1','Rspo2','Rspo3','Rspo4','Sbspon','Slamf6','Slit1','Slit2',
              'Slit3','Smoc1','Smoc2','Sned1','Sparc','Sparcl1','Spon1','Spon2','Spp1',
              'Srpx','Srpx2','Sspo','Svep1','Tecta','Tectb','Tgfbi','Thbs1','Thbs2',
              'Thbs3','Thbs4','Thsd4','Tinag','Tinagl1','Tnc','Tnfaip6','Tnn','Tnr',
              'Tnxb','Tsku','Tspear','Vit','Vtn','Vwa1','Vwa2','Vwa3a','Vwa3b','Vwa5a',
              'Vwa5b1','Vwa5b2','Vwa7','Vwa9','Vwce','Vwde','Vwf','Wisp1','Wisp2','Wisp3',
              'Zp1','Zp2','Zp3','Zp3r','Zpld1','Col10a1','Col11a1','Col11a2','Col12a1','Col13a1','Col14a1','Col15a1','Col16a1','Col17a1','Col18a1','Col19a1','Col1a1','Col1a2','Col20a1',
              'Col22a1','Col23a1','Col24a1','Col25a1','Col26a1','Col27a1','Col28a1','Col2a1','Col3a1',
              'Col4a1','Col4a2','Col4a3','Col4a4','Col4a5','Col4a6','Col5a1','Col5a2','Col5a3','Col6a1',
              'Col6a2','Col6a3','Col6a4','Col6a5','Col6a6','Col7a1','Col8a1','Col8a2','Col9a1','Col9a2',
              'Col9a3','Acan','Aspn','Bcan','Bgn','Chad','Chadl','Dcn','Epyc','Esm1','Fmod','Hapln1','Hapln2','Hapln3',
              'Hapln4','Hspg2','Impg1','Impg2','Kera','Lum','Ncan','Nepn',
              'Nyx','Ogn','Omd','Optc','Podn','Podnl1','Prelp','Prg2','Prg3','Spock2','Spock3','Srgn','Vcan')

#Read Files
Mock = readRDS("Set_Path/mock215_(SEURAT_v3)_02.rds")
Ligated = readRDS("Set_Path/Lig1130_(SEURAT_v3)_02.rds")

#Labeling Meta data based on condition
Mock@meta.data[,"protocol"] = "Mock"
Ligated@meta.data[,"protocol"] = "Ligated"

#merging datasets
Lig_Mock_Merge = merge(Mock, y = Ligated, add.cell.ids = c("Mock","Ligated"), project = "protocol")

#Subsetting based of Pdgfr and Gli1 expression
pdgfr_Gli_subset = subset(Lig_Mock_Merge, subset = Pdgfrb > .5 | Pdgfra > .5 | Gli1 > 0)

#Subset to compare Gli1 expressing cells by condition
Gli_subset_norm = subset(Lig_Mock_Merge, subset = Gli1 > 0)

#Data Scaling and Normalization
pdgfr_Gli_subset <- ScaleData(pdgfr_Gli_subset, verbose = TRUE)
pdgfr_Gli_subset = NormalizeData(pdgfr_Gli_subset)

#
Lig_Mock_Merge = ScaleData(Lig_Mock_Merge)
Lig_Mock_Merge = NormalizeData(Lig_Mock_Merge)

#
Gli_subset_norm = ScaleData(Gli_subset_norm)
Gli_subset_norm = NormalizeData(Gli_subset_norm)

#PCA

pdgfr_Gli_subset = FindVariableFeatures(pdgfr_Gli_subset, selection.method = "vst", nfeatures = 500)
pdgfr_Gli_subset <- RunPCA(pdgfr_Gli_subset, npcs = 20, verbose = TRUE)
#pdgfr_Gli_subset = JackStraw(pdgfr_Gli_subset, dims = 20, num.replicate = 100)
#pdgfr_Gli_subset = ScoreJackStraw(pdgfr_Gli_subset, dims = 1:20)
#JackStrawPlot(pdgfr_Gli_subset, dims = 1:20, ymax = .5)
#DimHeatmap(pdgfr_Gli_subset, dims = 1:20)

#

Lig_Mock_Merge = FindVariableFeatures(Lig_Mock_Merge, selection.method = "vst", nfeatures = 500)
Lig_Mock_Merge = RunPCA(Lig_Mock_Merge, npcs = 35)
#Lig_Mock_Merge = JackStraw(Lig_Mock_Merge, dims = 35, num.replicate = 100)
#Lig_Mock_Merge = ScoreJackStraw(Lig_Mock_Merge, dims = 1:35)
#JackStrawPlot(Lig_Mock_Merge, dims = 1:35, ymax = 1)
#DimHeatmap(Lig_Mock_Merge, dims = 1:15)
#DimHeatmap(Lig_Mock_Merge, dims = 16:35)

#

Gli_subset_norm = FindVariableFeatures(Gli_subset_norm, selection.method = "vst", nfeatures = 500)
Gli_subset_norm = RunPCA(Gli_subset_norm, npcs = 35)
#Gli_subset_norm = JackStraw(Gli_subset_norm, dims = 35, num.replicate = 100)
#Gli_subset_norm = ScoreJackStraw(Gli_subset_norm, dims = 1:35)
#JackStrawPlot(Gli_subset_norm, dims = 1:35, ymax = 1)
#DimHeatmap(Gli_subset_norm, dims = 1:15)
#DimHeatmap(Gli_subset_norm, dims = 16:35)

#

pdgfr_Gli_subset = FindNeighbors(pdgfr_Gli_subset, dims = 1:14)
#pdgfr_Gli_subset = FindClusters(pdgfr_Gli_subset, resolution = (seq(0, .03, by= .005)))
#clustree(pdgfr_Gli_subset, prefix = "RNA_snn_res.")
pdgfr_Gli_subset = FindClusters(pdgfr_Gli_subset, resolution = .025)

#

Lig_Mock_Merge = FindNeighbors(Lig_Mock_Merge, dims = 1:14)
#Lig_Mock_Merge = FindClusters(Lig_Mock_Merge, resolution = (seq(0, 1, by = .05)))
#clustree(Lig_Mock_Merge, prefix = "RNA_snn_res.")
Lig_Mock_Merge = FindClusters(Lig_Mock_Merge, resolution = .3)

#

Gli_subset_norm = FindNeighbors(Gli_subset_norm, dims = 1:14)
#Gli_subset_norm = FindClusters(Gli_subset_norm, resolution = (seq(0, .3, by = .01)))
#clustree(Gli_subset_norm, prefix = "RNA_snn_res.")
Gli_subset_norm = FindClusters(Gli_subset_norm, resolution = .1)

#Spatial visualization of genes

Lig_Mock_Merge <- RunUMAP(Lig_Mock_Merge, reduction = "pca", dims = 1:14)
FeaturePlot(Lig_Mock_Merge, features = c("Ptprcap","Pdgfra","Pdgfrb","Pecam1","Adgre1","Col1a1","Epcam","Cd74"))
VlnPlot(Lig_Mock_Merge, features = c("Ptprcap","Pdgfra","Pdgfrb","Pecam1","Adgre1","Col1a1","Epcam","Cd74","Gli1"))
VlnPlot(Lig_Mock_Merge, features = "Gli1")
#VlnPlot(Lig_Mock_Merge, features = "Adgre1", split.by = "protocol")
#VlnPlot(Lig_Mock_Merge, features = "Pdgfra", split.by = "protocol")
#VlnPlot(Lig_Mock_Merge, features = "Pdgfrb", split.by = "protocol")
#VlnPlot(Lig_Mock_Merge, features = "Pecam1", split.by = "protocol")
DimPlot(Lig_Mock_Merge,reduction = "umap",label = TRUE,label.size = 4,pt.size = .5,repel = TRUE)
DimPlot(Lig_Mock_Merge, reduction = "umap", split.by = "protocol")
p1 <- DimPlot(Lig_Mock_Merge, reduction = "umap", group.by = "protocol",)
p2 <- DimPlot(Lig_Mock_Merge, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
p1

#Generates Table for total gene expression, used for identifying cluster identities in whole dataset
Pos_genes_all = FindAllMarkers(Lig_Mock_Merge, only.pos = TRUE, return.thresh =  .000002498043)

#

pdgfr_Gli_subset <- RunUMAP(pdgfr_Gli_subset, reduction = "pca", dims = 1:14)
DimPlot(pdgfr_Gli_subset,reduction = "umap",label = TRUE,label.size = 10,pt.size = 1.5,repel = TRUE)+theme(text = element_text(size = 25),axis.text = element_text(size = 20))
VlnPlot(pdgfr_Gli_subset, features = c("Pdgfra","Pdgfrb","Col1a1","Gli1"))
VlnPlot(pdgfr_Gli_subset, features = c("Pdgfra","Pdgfrb","Col1a1","Gli1","Postn","Dpp4"),split.by = "protocol")
p3 <- DimPlot(pdgfr_Gli_subset, reduction = "umap", group.by = "protocol",)
p3
DimPlot(pdgfr_Gli_subset, reduction = "umap", split.by = "protocol")
FeaturePlot(pdgfr_Gli_subset, features = c("Pdgfra","Pdgfrb","Col1a1","Gli1","Adgre1","Pecam1"))

#

Gli_subset_norm <- RunUMAP(Gli_subset_norm, reduction = "pca", dims = 1:14)
DimPlot(Gli_subset_norm,reduction = "umap",label = TRUE,label.size = 4,pt.size = 2,repel = TRUE)
VlnPlot(Gli_subset_norm, features = c("Pdgfra","Pdgfrb","Col1a1","Gli1","Adgre1","Pecam1"))
VlnPlot(Gli_subset_norm, features = c("Pdgfra","Pdgfrb","Col1a1","Gli1","Adgre1","Pecam1"),split.by = "protocol")
p3 <- DimPlot(Gli_subset_norm, reduction = "umap", group.by = "protocol",pt.size = 2)
p3
DimPlot(Gli_subset_norm, reduction = "umap", split.by = "protocol")
FeaturePlot(Gli_subset_norm, features = c("Pdgfra","Pdgfrb","Col1a1","Gli1","Adgre1","Pecam1"))

#Changing identities of clusters for mock vs ligated analysis
Mod_merge = Lig_Mock_Merge
Mod_merge$cluster.protocol = paste(Idents(Mod_merge),Mod_merge$protocol,sep =  "_")
Mod_merge$cluster = Idents(Mod_merge)
Idents(Mod_merge) = 'cluster.protocol'

#

Mod_subset = pdgfr_Gli_subset
Mod_subset$cluster.protocol = paste(Idents(Mod_subset),Mod_subset$protocol,sep =  "_")
Mod_subset$cluster = Idents(Mod_subset)
Idents(Mod_subset) = 'cluster.protocol'

#

Mod_merge_gli1 = Gli_subset_norm
Mod_merge_gli1$cluster.protocol = paste(Idents(Mod_merge_gli1),Mod_merge_gli1$protocol,sep = "_")
Mod_merge_gli1$cluster = Idents(Mod_merge_gli1)
Idents(Mod_merge_gli1) = 'cluster.protocol'

#Markers
Mod_Gli1_Markers = FindAllMarkers(Mod_merge_gli1, return.thresh =  .000002498043,)
Mod_subset_Markers = FindAllMarkers(Mod_subset,return.thresh =  .000002498043)
pdgfr_Gli_subset_markers = FindAllMarkers(pdgfr_Gli_subset, only.pos = TRUE, return.thresh =  .000002498043)
Mod_Markers = FindAllMarkers(Mod_merge, return.thresh =  .000002498043,)

#cell count per cluster

table(Idents(Lig_Mock_Merge))
table(Idents(Mod_merge))
table(Idents(pdgfr_Gli_subset))
table(Idents(Gli_subset_norm))
table(Idents(Mod_subset))
table(Idents(Mod_merge_gli1))

#Dotplots
#DotPlot(Lig_Mock_Merge, features = c("Mgp","Dcn",'Bgn','Col1a2','Col3a1','Spp1','Col1a1','Cilp','Lum','Thbs4','Mfap5','Col6a1','Col6a3'),split.by = 'protocol')
#DotPlot(pdgfr_Gli_subset, features = c("Mgp","Dcn",'Bgn','Col1a2','Col3a1','Spp1','Col1a1','Cilp','Lum','Thbs4','Mfap5','Col6a1','Col6a3'),split.by = 'protocol')
#DotPlot(Mod_subset, features = c("Mgp","Dcn",'Bgn','Col1a2','Col3a1','Spp1','Col1a1','Cilp','Lum','Thbs4','Mfap5','Col6a1','Col6a3'))

#cluster_0 = subset(pdgfr_Gli_subset, idents = 0)
DotPlot(pdgfr_Gli_subset, features = genes,cols = "RdBu", split.by = 'protocol')+RotatedAxis()+ 
  theme(text = element_text(size = 15),axis.text = element_text(size = 15))#+theme(panel.background = 
# element_rect(fill = 'gray',colour = 'gray'))


levels(Mod_merge) = c("0_Mock","0_Ligated","1_Mock","1_Ligated","2_Mock","2_Ligated",
                      "3_Mock","3_Ligated", "4_Mock","4_Ligated",
                      "5_Mock","5_Ligated","6_Mock","6_Ligated","7_Mock",
                      "7_Ligated","8_Mock","8_Ligated","9_Mock","9_Ligated",
                      "10_Mock","10_Ligated","11_Mock","11_Ligated")
DoHeatmap(Mod_merge, features = Matrisome, size = 1.75)+ theme(axis.text.y = element_text(size = 5))

levels(Mod_subset) = c("0_Mock","0_Ligated","1_Mock","1_Ligated","2_Mock","2_Ligated",
                       "3_Mock","3_Ligated", "4_Mock","4_Ligated")
DoHeatmap(Mod_subset, features = Matrisome, size = 1.75)+ theme(axis.text.y = element_text(size = 5))




