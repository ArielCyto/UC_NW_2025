
#### Load environment ####
library("cytoreason.ccm.pipeline.qc")
library(cytoreason.assets, lib.loc = "/opt/R/4.4.2/lib/R/library")
library("cytoreason.individual.variation")
library(cytoreason.validator.apps.client)
library(cytoreason.cc.client)
library(Hmisc)
library(stringr)
library(dplyr)
library(checkmate)
library(testthat)
library(scales)
library(igraph)
library(visNetwork)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(patchwork)
library(ggplot2)
source("~/UC_NW_2025/Vis_functions.R")

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

##### Load UC Input - all datasets #####
UC_model <- as_ccm_fit('wf-2309373e83')
UC_config = modelMetadata(UC_model)
inf_datasets <- subset(UC_config,comparison %in% c("inf_vs_HC", "Inflamed_vs_non_inflamed", "network_baseline"))$dataset_id
inf_datasets <- unique(inf_datasets)

# Subset to mucosal datasets only
mucosal_inf_datasets <- c('E-MTAB-184','E-MTAB-2967','E-MTAB-7604__rnaseq','E-MTAB-7845__rnaseq','GSE16879__GPL570',
                        'GSE38713__GPL570','GSE47908__GPL570','GSE75214__GPL6244','GSE92415__GPL13158')
                        #,'GSE9452__GPL570' - removed due to low correlation,
                        # ,'GSE73661__GPL6244' - removed due to duplications with GSE75214__GPL6244)
inf_esets <- .get_esets_ccm(UC_model, esets_names = c(mucosal_inf_datasets))

##### Generate group table #####
group_table <- data.frame(matrix(ncol = 3, nrow = length(inf_esets)))
colnames(group_table) <- c('experiment_id', 'group_variable', 'reference_category')

for (dataset_num in 1:length(inf_esets)){
  group_table$experiment_id[dataset_num] <- names(inf_esets)[dataset_num]
  group_table$group_variable[dataset_num] <- 'all_dz'
  group_table$reference_category[dataset_num] <- 'DZ'
  metadata <- pData(inf_esets[[dataset_num]])
  metadata$all_dz <- 'DZ'
  pData(inf_esets[[dataset_num]]) <- metadata
  
}


# Filter the samples to only inf samples that passed QC:
esets_nw <- list()
for (dataset_num in 1:length(inf_esets)){
  sub_e <- inf_esets[[dataset_num]][,pData(inf_esets[[dataset_num]])$sample_classification %in% "Inflamed" & 
                                      pData(inf_esets[[dataset_num]])$condition %in% "ulcerative colitis"]
  sub_e
  exclude_samp <- unique(UC_config$exclude_samples)[!is.na(unique(UC_config$exclude_samples))]
  if(any(pData(sub_e)$sample_id %in% exclude_samp)) {
    sub_e <- sub_e[,!pData(sub_e)$sample_id %in% exclude_samp]
  }
  esets_nw[[dataset_num]] <- sub_e
}
names(esets_nw) <- names(inf_esets)


# manual changes

# esets_nw$GSE107593__rnaseq <- esets_nw$GSE107593__rnaseq[,sampleNames(esets_nw$GSE107593__rnaseq) %in% 
#                                                            esets_nw$GSE107593__rnaseq$sample_id[esets_nw$GSE107593__rnaseq$is_pediatric == 'ADULT']]

esets_nw$GSE16879__GPL570 <- esets_nw$GSE16879__GPL570[,sampleNames(esets_nw$GSE16879__GPL570) %in% 
                                                         esets_nw$GSE16879__GPL570$sample_id[esets_nw$GSE16879__GPL570$tissue == 'colon' & esets_nw$GSE16879__GPL570$time == 'W0']]

# esets_nw$GSE206285__GPL13158 <- inf_esets$GSE206285__GPL13158[,sampleNames(inf_esets$GSE206285__GPL13158) %in% 
#                                                                 inf_esets$GSE206285__GPL13158$sample_id[inf_esets$GSE206285__GPL13158$sample_classification == 'No specific classification found']]

esets_nw$`E-MTAB-7604__rnaseq` <- inf_esets$`E-MTAB-7604__rnaseq` [,sampleNames(inf_esets$`E-MTAB-7604__rnaseq` ) %in% 
                                                                     inf_esets$`E-MTAB-7604__rnaseq` $sample_id[inf_esets$`E-MTAB-7604__rnaseq`$condition == 'inflammatory bowel disease']]

esets_nw$GSE38713__GPL570 <- esets_nw$GSE38713__GPL570[,sampleNames(esets_nw$GSE38713__GPL570 ) %in% 
                                                         esets_nw$GSE38713__GPL570 $sample_id[esets_nw$GSE38713__GPL570$Condition_comment == 'active disease']]

# esets_nw$GSE193677__GPL16791 <- esets_nw$GSE193677__GPL16791[,sampleNames(esets_nw$GSE193677__GPL16791 ) %in% 
#                                                                esets_nw$GSE193677__GPL16791 $sample_id[esets_nw$GSE193677__GPL16791$Tissu_type == 'Colon']]

# esets_nw$GSE206171__GPL13667 <- esets_nw$GSE206171__GPL13667[,sampleNames(esets_nw$GSE206171__GPL13667 ) %in% 
#                                                                esets_nw$GSE206171__GPL13667 $sample_id[esets_nw$GSE206171__GPL13667$Condition_comment == 'active disease']]

# esets_nw$GSE73661__GPL6244 <- esets_nw$GSE73661__GPL6244[,sampleNames(esets_nw$GSE73661__GPL6244 ) %in% 
#                                                            esets_nw$GSE73661__GPL6244 $sample_id[esets_nw$GSE73661__GPL6244$time == 'W0']]



print(sapply(esets_nw, dim))
# upload to cyto-cc
# library('cytoreason.ccm.pipeline')
# tags <- list(list(name="owner", value="Ariel Simon"),
#              list(name="data", value="NW-UC-INF"),
#              list(name="details", value="esets of infalmed uc samples"),
#              list(name="n_datasets", value='22'))
# 
# wf <- run_function_dist(function(obj){return(obj)},
#                         obj=esets_nw, tags = tags)
# wf
# # https://cyto-cc.cytoreason.com/workflow/wf-0f7f0d29ed

##### Generate full corr matrix
net_id <- "UC_mucosal_inflamed_NW_9_datasets"
IMAGE <- "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest"
default_memory_request <- "200Gi"

##### Step 1 - generate full correlation matrix
full_corr_matrix_res <- get_full_corr_matrix(esets = esets_nw,
                                             group_table = group_table,
                                             IMAGE = IMAGE,
                                             net_id = net_id,
                                             gpu = TRUE,
                                             mem_request = default_memory_request)

#################################################################################

##### Step 3 -  Calculate meta gene-gene correlation matrix
edge_min_count <- 0
# if chunk_size is NULL, a default value (1000) will be used
chunk_size <- NULL
edge_meta_res <- meta_pvalue_edge(nets = full_corr_matrix_res[["reference_edges_wid"]],
                                  edge_min_count = edge_min_count,
                                  IMAGE = IMAGE,
                                  mem_req = default_memory_request,
                                  chunk_size = chunk_size)
get_workflow(edge_meta_res[["edge_meta"]][["wf_id"]], wait = TRUE)

all_edges <- cytoreason.io::read_data(cytoreason.cc.client::get_task_outputs(res=edge_meta_res[["edge_meta"]][["wf_id"]], task_id = NULL, task_name = edge_meta_res[["edge_meta"]][["task_name"]]))
q_pvalue_plot(all_edges)



##### Step 4 - Perform edge pruning analysis on meta gene-gene correlation matrix
edge_pruning_corr <- seq(0.2, 0.8, 0.05)
fdr_thresholds <- list(c(fdr_th=0.05, q_fdr_th=0.05))
important_genes <- NULL
pruning_info_wfid <- 
  run_function_dist(function(full_corr_matrix_res, ...) {
    pruning_info_res <- cytoreason.individual.variation::edge_pruning_th_comb(...)
    return(pruning_info_res)
  },
  nets = edge_meta_res[["edge_meta"]][["wf_id"]], 
  nets_task_name = edge_meta_res[["edge_meta"]][["task_name"]],
  cor_th = edge_pruning_corr, 
  fdr_ths = fdr_thresholds, 
  full_corr_matrix_res = full_corr_matrix_res,
  IMAGE = IMAGE,
  mem_req = default_memory_request,
  cytocc = TRUE,
  plot = FALSE,
  important_genes = important_genes,
  image = IMAGE,
  memory_request = default_memory_request,
  force_execution = FALSE, replace_image_tags = TRUE, 
  tags = list(list("name" = "NW_notebook", "value" =  "edge_pruning_analysis")))
# wf-08d37a2968

#get_workflow(pruning_info_wfid, wait = TRUE)
pruning_info <- get_outputs_dist('wf-08d37a2968')$output.rds$edge_pruning_info
p_nodes(pruning_info)
p_edges(pruning_info)
p_components(pruning_info)
p_paretofront(pruning_info)
if (!is.null(important_genes)){
  p_important_genes(pruning_info)
}
# print and check
unique(pruning_info$num_nodes[pruning_info$min_corr_threshold == '0.6'])
unique(pruning_info$num_nodes[pruning_info$min_corr_threshold == '0.55'])
unique(pruning_info$num_nodes[pruning_info$min_corr_threshold == '0.5'])

unique(pruning_info$num_edges[pruning_info$min_corr_threshold == '0.6'])
unique(pruning_info$num_edges[pruning_info$min_corr_threshold == '0.55'])
unique(pruning_info$num_edges[pruning_info$min_corr_threshold == '0.5'])

length(unique(pruning_info$node_measures.components[pruning_info$min_corr_threshold == '0.6']));
length(unique(pruning_info$node_measures.components[pruning_info$min_corr_threshold == '0.55']));
length(unique(pruning_info$node_measures.components[pruning_info$min_corr_threshold == '0.5']));

##### Step 5
# use case 1: postprocess edge meta output
nets = consume_workflow_task_output(wf = edge_meta_res[["edge_meta"]][["wf_id"]], 
                                    task_id = NULL,
                                    task_name = edge_meta_res[["edge_meta"]][["task_name"]])

# use case 2: postprocess on previous orchestrated_edge_pruning results

# See orchestrated_edge_pruning doc for details regarding parameters 
prune <- c("min_abs_corr" = 0.6, "max_fdr" = 0.05, "min_q_fdr" = 0.05)
corr_direction <- "none"
subset_gene_list <- NULL
min_connected_component <- NULL
filter_components <- NULL
post_processes_nw <- 
  run_function_dist(cytoreason.individual.variation::orchestrated_edge_pruning,
                    nets = nets, 
                    prune = prune, 
                    corr_direction = corr_direction, 
                    subset_gene_list = subset_gene_list,
                    min_connected_component = min_connected_component,
                    filter_components = filter_components,
                    image = IMAGE,
                    memory_request = default_memory_request,
                    force_execution = FALSE, replace_image_tags = TRUE, 
                    tags = list(list("name" = "NW_notebook", "value" = "orchestrated_edge_pruning")))

# https://cyto-cc.cytoreason.com/workflow/wf-c3b19c53be
# wf-ae634ffd00
post_processes_nw <- readRDS(get_task_outputs('wf-c3b19c53be', task_id = 0, files_names_grepl_pattern = "output.rds"))

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
# https://cyto-cc.cytoreason.com/workflow/wf-f92e284fb6

#get_workflow(nw_measures, wait = TRUE)
nw_measures_res <- readRDS(get_task_outputs('wf-f92e284fb6', task_id = 0, files_names_grepl_pattern = "output.rds"))
# nw_measures_res <- readRDS(get_task_outputs(nw_measures, task_id = 0, files_names_grepl_pattern = "output.rds"))
#head(nw_measures_res)
write.csv(nw_measures_res, '~/UC_NW_2025/outputs/inf_mucosal_nw_measures_res_0.6.csv')


######## Vis
RA_nw_bulk_pruned_0.5_degree_hist <- ggplot(nw_measures_res, aes(x = degree)) +
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "black") +
  labs(title = paste0("Degree distribution - UC Inflamed Mucosal Network"),
       x = "Node degree",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))
RA_nw_bulk_pruned_0.5_degree_hist

RA_nw_bulk_pruned_0.5_eigen_centrality_hist <- ggplot(nw_measures_res, aes(x = eigen_centrality)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = paste0("Eigen centrality distribution - UC Inflamed Mucosal Network"),
       x = "Node eigen centrality",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))
RA_nw_bulk_pruned_0.5_eigen_centrality_hist

RA_nw_bulk_pruned_0.5_betweenness_hist <- ggplot(nw_measures_res, aes(x = betweenness)) +
  geom_histogram(binwidth = 0.0005, fill = "skyblue", color = "black") +
  labs(title = paste0("Betweenness distribution - UC Inflamed Mucosal Network"),
       x = "Node betweenness",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))
RA_nw_bulk_pruned_0.5_betweenness_hist


RA_nw_bulk_pruned_0.5_page_rank_hist <- ggplot(nw_measures_res, aes(x = page_rank)) +
  geom_histogram(binwidth = 0.00005, fill = "skyblue", color = "black") +
  labs(title = paste0("Page rank distribution - UC Inflamed Mucosal Network"),
       x = "Node page rank",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))
RA_nw_bulk_pruned_0.5_page_rank_hist

RA_nw_bulk_pruned_0.5_closeness_hist <- ggplot(nw_measures_res, aes(x = closeness)) +
  geom_histogram(binwidth = 0.0005, fill = "skyblue", color = "black") +
  labs(title = paste0("Closeness distribution - UC Inflamed Mucosal Network"),
       x = "Node closeness",
       y = "No. nodes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5,face='bold'))
RA_nw_bulk_pruned_0.5_closeness_hist


### Check max distance (for unweighted network)
# igraph object
nw_bulk_pruned_graph <- igraph::graph_from_data_frame(d=nw_measures_res, directed = FALSE)

# pairwise distances
gene_gene_distances <- igraph::distances(nw_bulk_pruned_graph, mode = "all", algorithm = "unweighted")

# maximum distance
max_distance <- max(gene_gene_distances[gene_gene_distances != Inf])  # Exclude infinite values (if there are disconnected components)
print(max_distance)  # results show 6


### Bio QC for Gali
# post_processes_nw_short <- (post_processes_nw[,c(2,3,9)])
# colnames(post_processes_nw_short) <- c('from', 'to', 'cor')
post_processes_nw_chec <- post_processes_nw
colnames(post_processes_nw_chec)[colnames(post_processes_nw_chec) == 'Var1'] <- 'From'
colnames(post_processes_nw_chec)[colnames(post_processes_nw_chec) == 'Var2'] <- 'To'
post_processes_nw_as_graph <- igraph::graph_from_data_frame(d=post_processes_nw_chec[,2:10])
visualize_gene_set_subgraphs(post_processes_nw_as_graph , mono_genesets[5], nw_measures_res, centrality_measure = "page_rank")

