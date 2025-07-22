#### Load environment ####
library("cytoreason.ccm.pipeline.qc")
library("cytoreason.individual.variation")
library(cytoreason.validator.apps.client)
library(cytoreason.cc.client)
library(Hmisc)
library(testthat)
library(scales)
library(igraph)
library(visNetwork)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(patchwork)
library(ggplot2)
library(cytoreason.assets, lib.loc = "/opt/R/4.4.2/lib/R/library")
source("~/UC_NW_2025/Vis_functions.R")
devtools::load_all("~/analysis-p03-poc1/")
.get_esets_ccm <- function(ccmfit, esets_names=NULL) {
  if (is.null(esets_names)) {
    esets_names <- names(ccmfit$datasets)
  }
  
  esets <- lapply(esets_names, function(nam) assayDataExpression(ccmfit, nam))
  names(esets) <- esets_names
  esets
}

#### Load Inputs ####
##### Load input example #####
# moriah_inputs <- sapply(cytoreason.cc.client::get_task_inputs("wf-2ee7856ff0", "0"), readRDS)
# esets_example <- moriah_inputs$esets
# group_table_example <- moriah_inputs$group_table

##### Load UC Input - 6 datasets #####
UC_model <- as_ccm_fit('wf-2309373e83')
UC_config = modelMetadata(UC_model)
hc_datasets <- subset(UC_config,comparison %in% c("UC_vs_HC", "hc_vs_HC", "UC_non_inf_vs_HC"))$dataset_id
hc_datasets <- unique(hc_datasets)
# check
hc_datasets %in% names(UC_model$datasets)
# subset to wanted datasets
hc_esets <- .get_esets_ccm(UC_model, esets_names = c(hc_datasets))
hc_esets$GSE83687__rnaseq <- NULL # tumor adjacent
hc_esets$GSE53306__GPL14951 <- NULL # an outlier in the inflamed NW

##### Generate group table #####
group_table <- data.frame(matrix(ncol = 3, nrow = length(hc_esets)))
colnames(group_table) <- c('experiment_id', 'group_variable', 'reference_category')

for (dataset_num in 1:length(hc_esets)){
  group_table$experiment_id[dataset_num] <- names(hc_esets)[dataset_num]
  group_table$group_variable[dataset_num] <- 'all_dz'
  group_table$reference_category[dataset_num] <- 'DZ'
  metadata <- pData(hc_esets[[dataset_num]])
  metadata$all_dz <- 'DZ'
  pData(hc_esets[[dataset_num]]) <- metadata
  
}


# Filter the samples to only hc samples that passed QC:
esets_nw <- list()
for (dataset_num in 1:length(hc_esets)){
  sub_e <- hc_esets[[dataset_num]][, pData(hc_esets[[dataset_num]])$condition %in% c("healthy","control")]
  sub_e
  exclude_samp <- unique(UC_config$exclude_samples)[!is.na(unique(UC_config$exclude_samples))]
  if(any(pData(sub_e)$sample_id %in% exclude_samp)) {
    sub_e <- sub_e[,!pData(sub_e)$sample_id %in% exclude_samp]
  }
  esets_nw[[dataset_num]] <- sub_e
}
names(esets_nw) <- names(hc_esets)

## keep 1 per patient
index_to_keep <- seq_len(length(esets_nw$GSE206171__GPL13667$subject_id)) %%2 != 1 
esets_nw$GSE206171__GPL13667 <- esets_nw$GSE206171__GPL13667[,sampleNames(esets_nw$GSE206171__GPL13667 ) %in% esets_nw$GSE206171__GPL13667$sample_id[index_to_keep]]


print(sapply(esets_nw, dim))


##### Generate full corr matrix
net_id <- "hc_nw_UC"
IMAGE <- "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest"
default_memory_request <- "200Gi"

##### Step 1 - generate full correlation matrix
full_corr_matrix_res <- get_full_corr_matrix(esets = esets_nw,
                                             group_table = group_table,
                                             IMAGE = IMAGE,
                                             net_id = net_id,
                                             gpu = TRUE,
                                             mem_request = default_memory_request)
full_corr_matrix_res <- read_asset('wf-3283822bf6')
##### Step 2 - check outliers, if any is above 0.2 - remove and repeat step 1

outliers <- check_outliers(reference_edges = "wf-3283822bf6", image_url = "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest", mem_request="20Gi", cyto_cc_force_execution=FALSE)
corrplot::corrplot(as.matrix(outliers), order = "hclust", addCoef.col = "black",
                   addCoefasPercent=FALSE, method="color",
                   col=rev(RColorBrewer::brewer.pal(n=10, name="RdBu")))

##### Step 3 -  Calculate meta gene-gene correlation matrix ##### 
edge_min_count <- 0
# if chunk_size is NULL, a default value (1000) will be used
chunk_size <- NULL
edge_meta_res <- meta_pvalue_edge(nets ='wf-3283822bf6', # nets = full_corr_matrix_res$wf_id
                                  edge_min_count = edge_min_count,
                                  IMAGE = IMAGE,
                                  mem_req = default_memory_request,
                                  chunk_size = chunk_size)
# https://cyto-cc.cytoreason.com/workflow/wf-f8cc81df81
get_workflow(edge_meta_res[["edge_meta"]][["wf_id"]], wait = TRUE)

###### plot edges q_pvalue_fdr ###### 
all_edges <- cytoreason.io::read_data(cytoreason.cc.client::get_task_outputs(res=edge_meta_res[["edge_meta"]][["wf_id"]], task_id = NULL, task_name = edge_meta_res[["edge_meta"]][["task_name"]]))
#all_edges <- cytoreason.io::read_data(cytoreason.cc.client::get_task_outputs(res= 'wf-f8cc81df81', task_id = NULL, task_name = 'meta-pvalue-post'))

library(ggplot2)
p <- ggplot(all_edges, aes(x = q_pvalue_fdr)) +
  geom_histogram(binwidth = 0.002, fill = "skyblue", color = "black") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))
p + ylim(0,(5*10^6)) +geom_vline(xintercept = 0.02)


##### Step 4 - Perform edge pruning analysis on meta gene-gene correlation matrix ##### 
edge_pruning_corr <- seq(0, 1, 0.05)
fdr_thresholds <- list(c(fdr_th=0.05, q_fdr_th=0.02))
important_genes <- NULL
pruning_info_wfid <- 
  run_function_dist(function(full_corr_matrix_res, ...) {
    pruning_info_res <- cytoreason.individual.variation::edge_pruning_th_comb(...)
    return(pruning_info_res)
  },
  nets = 'wf-f8cc81df81', 
  nets_task_name = 'meta-pvalue-post',
  cor_th = edge_pruning_corr, 
  fdr_ths = fdr_thresholds, 
  full_corr_matrix_res = 'wf-3283822bf6',
  IMAGE = IMAGE,
  mem_req = default_memory_request,
  cytocc = TRUE,
  plot = FALSE, 
  important_genes = important_genes,
  image = IMAGE,
  memory_request = default_memory_request,
  force_execution = FALSE, replace_image_tags = TRUE, 
  tags = list(list("name" = "NW_notebook", "value" =  "edge_pruning_analysis")))
# https://cyto-cc.cytoreason.com/workflow/wf-6db2ba42b9


###### print and check pruning info ###### 
get_workflow(pruning_info_wfid, wait = TRUE)
pruning_info <- get_outputs_dist('wf-6db2ba42b9')$output.rds$edge_pruning_info
pruning_info_pos <- pruning_info[pruning_info$co]
p_paretofront(pruning_info)
p_nodes(pruning_info)
p_edges(pruning_info)
p_components(pruning_info)
unique(pruning_info$num_nodes[pruning_info$min_corr_threshold == '0.45']);
unique(pruning_info$num_edges[pruning_info$min_corr_threshold == '0.45']);
length(unique(pruning_info$node_measures.components[pruning_info$min_corr_threshold == '0.45']));

###### count connected component ###### 
conn_comp <- as.data.frame(table(pruning_info$node_measures.components[pruning_info$min_corr_threshold == '0.45']))
conn_comp
conn_comp$Freq[1]/sum(conn_comp$Freq)

##### Step 5 - postprocess edge meta output #####
nets = consume_workflow_task_output(wf = 'wf-f8cc81df81', # wf = edge_meta_res$wf_id
                                    task_id = NULL,
                                    task_name = 'meta-pvalue-post')

# nets = consume_workflow_task_output(wf = edge_meta_res[["edge_meta"]][["wf_id"]], 
#                                     task_id = NULL,
#                                     task_name = edge_meta_res[["edge_meta"]][["task_name"]])

# See orchestrated_edge_pruning doc for details regarding parameters 
prune <- c("min_abs_corr" = 0.45, "max_fdr" = 0.05, "min_q_fdr" = 0.02)
corr_direction <- "positive_only"
subset_gene_list <- NULL
min_connected_component <- 4
filter_components <- NULL
post_processes_nw <- 
  run_function_dist(cytoreason.individual.variation::orchestrated_edge_pruning,
                    nets = nets, 
                    prune = prune,
                    corr_direction=corr_direction,
                    subset_gene_list = subset_gene_list,
                    min_connected_component = min_connected_component,
                    filter_components = filter_components,
                    image = IMAGE,
                    memory_request = default_memory_request,
                    force_execution = FALSE, replace_image_tags = TRUE, 
                    tags = list(list("name" = "uc_hc_network", "value" = "orchestrated_edge_pruning")))
# wf-96cbcdc4ae - Arnon fixed


##### Step 6 - Compute NW Measurements per node ##### 
nw_measures <- 
  run_function_dist(function(nets, ...) {
    graph <- cytoreason.individual.variation::get_as_graph(nets, nodes_columns = c("Var1","Var2"))
    cytoreason.individual.variation::compute_measures(graph, ...)
  },
  nets = consume_workflow_task_output(wf = post_processes_nw), 
  cytocc = TRUE, 
  mem_req = default_memory_request,
  IMAGE = IMAGE,
  image = IMAGE,
  memory_request = default_memory_request, 
  force_execution = FALSE, replace_image_tags = TRUE, 
  tags = list(list("name" = "NW_notebook", "value" = "compute_measures")))
# https://cyto-cc.cytoreason.com/workflow/wf-b887714056 - Arnon fixed

get_workflow(nw_measures, wait = TRUE)
#nw_measures_res <- readRDS(get_task_outputs(nw_measures, task_id = 0, files_names_grepl_pattern = "output.rds"))
nw_measures_res <- readRDS(get_task_outputs('wf-b887714056', task_id = 0, files_names_grepl_pattern = "output.rds"))
#write.csv(nw_measures_res, '~/UC_NW_2025/outputs/hc_nw_measures_res_0.45.csv')

##### Internal Visualization - Ariel & Gali ##### 
###### Centrality histograms ###### 
# degree hist
library(ggplot2)
ggplot(nw_measures_res, aes(x = degree)) +
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "black") +
  labs(title = paste0("Degree distribution - UC healthy Network"),
       x = "Node degree",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))

# eigen_centrality_hist
ggplot(nw_measures_res, aes(x = eigen_centrality)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = paste0("Eigen centrality distribution - UC healthy Network"),
       x = "Node eigen centrality",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))

# betweenness_hist
ggplot(nw_measures_res, aes(x = betweenness)) +
  geom_histogram(binwidth = 0.0005, fill = "skyblue", color = "black") +
  labs(title = paste0("Betweenness distribution - UC healthy Network"),
       x = "Node betweenness",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))

# Page rank
ggplot(nw_measures_res, aes(x = page_rank)) +
  geom_histogram(binwidth = 0.00005, fill = "skyblue", color = "black") +
  labs(title = paste0("Page rank distribution - UC healthy Network"),
       x = "Node page rank",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))

# closeness_hist 
ggplot(nw_measures_res, aes(x = closeness)) +
  geom_histogram(binwidth = 0.0005, fill = "skyblue", color = "black") +
  labs(title = paste0("Closeness distribution - UC healthy Network"),
       x = "Node closeness",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))



# ### Check max distance (for unweighted network)
# igraph object
nw_bulk_pruned_graph <- igraph::graph_from_data_frame(d=nw_measures_res, directed = FALSE)
# pairwise distances
gene_gene_distances <- igraph::distances(nw_bulk_pruned_graph, mode = "all", algorithm = "unweighted")
# maximum distance
max_distance <- max(gene_gene_distances[gene_gene_distances != Inf])  # Exclude infinite values (if there are disconnected components)
print(max_distance)  # results show 4


###### Bio QC for sub-graphs  ###### 
customSignatures <- read_asset("wf-13096e7be3")
combo_genesets <- customSignatures$dupiMono
mono_genesets <- customSignatures$monotherapies
mono_genesets
post_processes_nw_tmp <- readRDS(get_task_outputs(post_processes_nw$workflow_id, task_id = 0, files_names_grepl_pattern = "output.rds"))
#post_processes_nw_tmp <- readRDS(get_task_outputs('wf-96cbcdc4ae', task_id = 0, files_names_grepl_pattern = "output.rds"))
colnames(post_processes_nw_tmp)[colnames(post_processes_nw_tmp) == 'Var1'] <- 'From'
colnames(post_processes_nw_tmp)[colnames(post_processes_nw_tmp) == 'Var2'] <- 'To'
post_processes_nw_as_graph <- igraph::graph_from_data_frame(d=post_processes_nw_tmp[,2:10])
# positive - IL23, TNF, TSLP, IL4, IL6
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[19], nw_measures_res, centrality_measure = "page_rank") #IL6
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[22], nw_measures_res, centrality_measure = "page_rank") #TNFa
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[12], nw_measures_res, centrality_measure = "page_rank") #IL23
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[18], nw_measures_res, centrality_measure = "page_rank") #IL4
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[23], nw_measures_res, centrality_measure = "page_rank") #TSLP

# negative - BMP7, KLK5, KLK7, GPR15L
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[33], nw_measures_res, centrality_measure = "page_rank") #BMP7
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[26], nw_measures_res, centrality_measure = "page_rank") #GPR15L
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[28], nw_measures_res, centrality_measure = "page_rank") #KLK5
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[29], nw_measures_res, centrality_measure = "page_rank") #KLK7


###### Compare bio-exp to random- Run Shiran's code ###### 
devtools::load_all("~/analysis-p03-poc1/")
customSignatures <- read_asset("wf-13096e7be3") # top 50 monotherapies + Dupi with random
combined_customSignatures <- list(monotherapies=customSignatures$monotherapies)
fullList = unlist(combined_customSignatures, recursive = FALSE)
names(fullList) = stringr::str_replace(names(fullList),"\\.","__")

# background
#nw_measures_res <- readRDS(get_task_outputs('wf-b887714056', task_id = 0, files_names_grepl_pattern = "output.rds"))
EntrezInput = unique(as.character(nw_measures_res$feature_id))
nIter <- 100 #updating number of iterations

bulk_centrality_wf <- "wf-b887714056" # nw_measures_res

bulk_nw_cent <-  cytoreason.cc.client::run_method_dist(
  method = "NWParam_VersusRandom_modified",  
  ns = "analysis.p03.POC1",
  NW_WF = bulk_centrality_wf, 
  Signaturelist = fullList,
  BackgroundGeneInput = EntrezInput,
  memory_request = "40Gi",                                  
  Image= "eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest", 
  image="eu.gcr.io/cytoreason/ci-analysis.p03.poc1-package:shiran_API_latest",
  Niter = nIter)

# wf-f034d1b873

#process bulk geneset centrality
library(tidyr)
library(dplyr)
library(scales)
library(ggpubr, lib.loc = "/opt/R/4.4.2/lib/R/library")
BulkAD_NW_cent <- read_asset("wf-f034d1b873")
BulkAD_NW_cent_pval <- BulkAD_NW_cent$Significance
BulkAD_NW_cent_pval$Criteria <- paste(BulkAD_NW_cent_pval$Criteria, BulkAD_NW_cent_pval$Measure, sep="_")
BulkAD_NW_cent_pval$Criteria[BulkAD_NW_cent_pval$Criteria %in% "OverlapProbability_NA"] <- "OverlapPercent"

BulkAD_NW_cent_values <- BulkAD_NW_cent$FullParam
BulkAD_NW_cent_values <- BulkAD_NW_cent_values[BulkAD_NW_cent_values$ListType %in% "gene_list",]

#add values in addition to pval
BulkAD_NW_cent_values_long <- reshape2::melt(BulkAD_NW_cent_values , id.vars = c("ListType", "ListName"), variable.name = "Criteria")
BulkAD_NW_cent <- merge(BulkAD_NW_cent_pval, BulkAD_NW_cent_values_long, by=c("ListName", "Criteria"), all.x=T)


centralitySignificance_ALL = cbind(BulkAD_NW_cent,list(Type = "bulk"))
centralitySignificance_ALL <-   centralitySignificance_ALL %>% dplyr::mutate(pval = ifelse(pval == 0, 1/(nIter+1), pval))

centralitySignificance_ALL <- centralitySignificance_ALL %>%
  separate(ListName, into = c("collection", "pathway"), sep = "__", extra = "merge", remove = FALSE)
centralitySignificance_ALL$Identifier <- centralitySignificance_ALL$pathway


#convert to final table structure
#-----raw measure*significance---#

nwCentrality <-  centralitySignificance_ALL %>%
  dplyr::mutate(Identifier = ListName) %>%
  dplyr::mutate(Subset = "L",
                DataType = "subnetwork_centrality") %>%
  dplyr::mutate(Disease = "UC") %>%
  dplyr::mutate(NL.normalize = case_when(Subset == "L" ~ "no",
                                         Subset == "L_vs_NL" ~ "yes")) %>%
  dplyr::mutate(adjusted = case_when(Type == "bulk" ~ "no",
                                     Type == "adjusted" ~ "yes")) %>%
  dplyr::group_by(Criteria,Type, collection) %>%
  dplyr::mutate(
    # Check if all 'value' or 'pval' are identical
    all_values_identical = n_distinct(value) == 1,
    all_pvals_identical = n_distinct(pval) == 1,
    # Rescale but if all identical values, set the scaled value to 1
    scaled_value = if_else(all_values_identical, 1, rescale(x = value, to = c(0.1, 1))),
    scaled_pval = if_else(all_pvals_identical, 1, rescale(x = (-log10(pval)), to = c(0.1, 1))),
    # Compute the score as product of scaled value and scaled pval
    Score = scaled_value * scaled_pval) %>%
  #  dplyr::mutate(scaled_value = rescale(x = value, to=c(0.1,1)), 
  #                scaled_pval = rescale(x = (-log10(pval)), to=c(0.1,1)), 
  #                Score = scaled_value*scaled_pval) %>%
  dplyr::mutate(SigScaled = Score) %>%
  # dplyr::mutate(SigScaled = rescale(x = Score)) %>%
  dplyr::select(Identifier,Criteria,Type,Subset,Score, SigScaled,DataType,value, Disease, pval, NL.normalize, adjusted, scaled_value,scaled_pval, collection)

nwCentrality$Type <- factor(nwCentrality$Type, levels=c("bulk", "adjusted"))

targetColors = read_asset("wf-c588727d1c") #bottom BTLA & SLAMF6

#targetColors['PD1'] <- "#FFE4C4"
names(targetColors)[names(targetColors) %in% "IL13.Dupi"] <- "IL13.minFDR"
names(targetColors)[names(targetColors) %in% "IL4.Dupi"] <- "IL4.minFDR"

nwCentrality <- nwCentrality[!grepl("3sd", nwCentrality$Identifier),]
# nwCentrality <- nwCentrality[nwCentrality$BroadCategory %in% c("Target", "negativeControls"),]
nwCentrality <- nwCentrality[!grepl("avg", nwCentrality$Identifier),]
nwCentrality$SubCategory <- nwCentrality$Identifier
nwCentrality$SubCategory <- gsub("monotherapies__","",nwCentrality$SubCategory)
#nwCentrality$SubCategory[nwCentrality$Identifier == 'monotherapies__ad-permutations'] <- 'ad-permutations'
nwCentrality$log10pval <- -log10(nwCentrality$pval)

nwCentrality_page_rank <- nwCentrality[nwCentrality$Criteria %in% "page_rank_median",]
nwCentrality_page_rank[nwCentrality_page_rank$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"

#ZeroCenteredProbability = -(2*(0.05 - .5))
#NegLogPVAL_thr = -log10(1-abs(ZeroCenteredProbability))

nwCentrality_page_rank$SubCategory <- factor(nwCentrality_page_rank$SubCategory, levels=rev(unique(filtered_data$SubCategory)))
nwCentrality_page_rank <- nwCentrality_page_rank[!is.na(nwCentrality_page_rank$SubCategory),]

tail_probability__page_rank_g <- ggplot(nwCentrality_page_rank, aes(x= reorder(SubCategory, value), y=-log10(pval), fill=SubCategory)) + 
  facet_wrap(~Type, ncol=1)+
  geom_bar(position=position_dodge(), stat="identity",color = "black") + 
  scale_fill_manual(values = targetColors)+
  coord_flip()+
  border() +
  ylab("Permutation Tail \n Probability") + xlab("")+
  geom_hline(yintercept=1.0, linetype="dashed")+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")

tail_probability__page_rank_g

###########################3
# AD_NW_cent_df_sub_combined <- read_asset("wf-c8c1f9bf0f")
BulkAD_NW_cent <- read_asset("wf-f034d1b873") # check what is missing
BulkAD_NW_cent <- BulkAD_NW_cent$FullParam
BulkAD_NW_cent$Type <- 'bulk'
BulkAD_NW_cent$BroadCategory <- 'Target'
BulkAD_NW_cent$Identifier <-  gsub("monotherapies__","",BulkAD_NW_cent$ListName)
BulkAD_NW_cent$Identifier[BulkAD_NW_cent$Identifier == 'ad-permutations'] <-  'UC-permutations'

BulkAD_NW_cent$CleanName <- BulkAD_NW_cent$Identifier
BulkAD_NW_cent$SubCategory <-  BulkAD_NW_cent$Identifier
BulkAD_NW_cent$BroadCategory[BulkAD_NW_cent$Identifier %in% c('KLK5', 'KLK7', 'GPR15L', 'BMP7','random-fc' ,'UC-permutations')] <- 'negativeControls'
# colnames(AD_NW_cent_df_sub_combined)[!colnames(AD_NW_cent_df_sub_combined) %in% colnames(BulkAD_NW_cent)]


#remove largest component criterion

##### page rank gene_list vs random#####
AD_NW_cent_df_sub_combined <- BulkAD_NW_cent
AD_NW_cent_df_sub_combined_Target <- AD_NW_cent_df_sub_combined[!AD_NW_cent_df_sub_combined$Identifier %in% c("steadyArm__IL13.avgFDR","steadyArm__IL4.avgFDR"),]

AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[!grepl("3sd", AD_NW_cent_df_sub_combined$ListName),]
AD_NW_cent_df_sub_combined$Type <- factor(AD_NW_cent_df_sub_combined$Type)
AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[!grepl("avg", AD_NW_cent_df_sub_combined$SubCategory),]

#fix the factor (order) of the subcategories
filtered_data <- AD_NW_cent_df_sub_combined %>%
  dplyr::filter(ListType == "gene_list" & BroadCategory %in% c("Target", "negativeControls")) %>%
  dplyr::mutate(Type = factor(Type, levels = c("bulk", "adjusted"))) %>%
  # Ensure BroadCategory has "negativeControls" last
  dplyr::mutate(BroadCategory = factor(BroadCategory, levels = c("Target", "negativeControls"))) %>%
  # Order SubCategory by page_rank_median within each Type
  dplyr::group_by(Type, BroadCategory) %>%
  dplyr::arrange(Type, BroadCategory, desc(page_rank_median)) %>%
  dplyr::mutate(SubCategory = factor(SubCategory, levels = rev(unique(SubCategory)))) %>%
  dplyr::ungroup()

filtered_data$SubCategory <- as.character(filtered_data$SubCategory)
filtered_data[filtered_data$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"
filtered_data$SubCategory <- factor(filtered_data$SubCategory, levels=rev(unique(filtered_data$SubCategory)))

library(patchwork)

AD_NW_cent_df_sub_combined <- AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$Type %in% "bulk",]

max_values <- AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$ListType %in% "gene_list",] %>%
  dplyr::filter(BroadCategory %in% "negativeControls") %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(max_page_rank_median = max(page_rank_median, na.rm = TRUE)) 

AD_NW_cent_df_sub_combined[AD_NW_cent_df_sub_combined$SubCategory %in% "PD1",'SubCategory'] <- "PDL1"

page_rank_median_g = AD_NW_cent_df_sub_combined %>%
  dplyr::filter(BroadCategory %in% c("Target", "negativeControls")) %>%
  dplyr::mutate(Type = factor(Type,levels = c("bulk","adjusted"))) %>%
  dplyr::select(c(ListName,page_rank_median,Type,SubCategory,CleanName,ListType)) %>%
  dplyr::mutate(SubCategory = factor(SubCategory,levels = rev(as.character(unique(filtered_data$SubCategory))))) %>%
  unique() %>%
  ggplot(aes(x = SubCategory,y = page_rank_median)) +
  facet_wrap(~ Type ,scales = "free",ncol = 1)+
  geom_point(data = function(x) subset(x,ListType != "gene_list"),aes(color = SubCategory), alpha=0.5 ,color = "grey") +
  geom_point(data = function(x) subset(x,ListType == "gene_list"),aes(color = SubCategory), size=3) +
  geom_hline(data = max_values, aes(yintercept = max_page_rank_median), linetype = "dashed", color = "black") +  # Add hline here
  coord_flip()+
  theme_minimal() + 
  theme(legend.position = "none")+
  #  theme(aspect.ratio = 1.5) + 
  border() +ylab("Median page rank - HC UC Targets") + xlab("") +
  #  scale_x_log10() +
  scale_color_manual(values = targetColors)
page_rank_median_g
