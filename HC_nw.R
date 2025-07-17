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
min_connected_component <- NULL
filter_components <- as.integer(c(3))
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
# wf-e9b7254e63

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
# https://cyto-cc.cytoreason.com/workflow/wf-58b9314909
get_workflow(nw_measures, wait = TRUE)
#nw_measures_res <- readRDS(get_task_outputs(nw_measures, task_id = 0, files_names_grepl_pattern = "output.rds"))
nw_measures_res <- readRDS(get_task_outputs('wf-58b9314909', task_id = 0, files_names_grepl_pattern = "output.rds"))
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
#post_processes_nw_tmp <- readRDS(get_task_outputs('wf-e9b7254e63', task_id = 0, files_names_grepl_pattern = "output.rds"))
colnames(post_processes_nw_tmp)[colnames(post_processes_nw_tmp) == 'Var1'] <- 'From'
colnames(post_processes_nw_tmp)[colnames(post_processes_nw_tmp) == 'Var2'] <- 'To'
post_processes_nw_as_graph <- igraph::graph_from_data_frame(d=post_processes_nw_tmp[,2:10])
# IL23, TNF, TSLP, IL4, IL6
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[19], nw_measures_res, centrality_measure = "page_rank") #IL6
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[33], nw_measures_res, centrality_measure = "page_rank") #TNFa
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[12], nw_measures_res, centrality_measure = "page_rank") #IL23
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[18], nw_measures_res, centrality_measure = "page_rank") #IL4
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[23], nw_measures_res, centrality_measure = "page_rank") #TSLP



###### Compare bio-exp to random- Run Shiran's code ###### 
devtools::load_all("~/analysis-p03-poc1/")
customSignatures <- read_asset("wf-13096e7be3") # top 50 monotherapies + Dupi with random
combined_customSignatures <- list(monotherapies=customSignatures$monotherapies)
fullList = unlist(combined_customSignatures, recursive = FALSE)
names(fullList) = stringr::str_replace(names(fullList),"\\.","__")

# background
#nw_measures_res <- readRDS(get_task_outputs('wf-58b9314909', task_id = 0, files_names_grepl_pattern = "output.rds"))
EntrezInput = unique(as.character(nw_measures_res$feature_id))
nIter <- 100 #updating number of iterations

bulk_centrality_wf <- "wf-58b9314909" # nw_measures_res

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

# wf-65a5e3bfa0

