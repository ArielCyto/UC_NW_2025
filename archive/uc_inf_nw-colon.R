#### Load environment ####
library("cytoreason.ccm.pipeline.qc")
library("cytoreason.individual.variation")
library(cytoreason.validator.apps.client)
library(cytoreason.cc.client)
library(Hmisc)
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

# Subset to colon datasets only
colon_inf_datasets <- c('GSE107593__rnaseq','GSE174159__rnaseq','GSE206285__GPL13158','GSE193677__GPL16791',
                        'GSE206171__GPL13667','GSE179285__GPL6480','GSE72819__rnaseq','GSE83687__rnaseq',
                          'GSE23597__GPL570','GSE66407__GPL13667') #,'GSE95437__rnaseq' - removed due to low corr and low sample size (11))
inf_esets <- .get_esets_ccm(UC_model, esets_names = c(colon_inf_datasets))

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

esets_nw$GSE107593__rnaseq <- esets_nw$GSE107593__rnaseq[,sampleNames(esets_nw$GSE107593__rnaseq) %in%
                                                           esets_nw$GSE107593__rnaseq$sample_id[esets_nw$GSE107593__rnaseq$is_pediatric == 'ADULT']]

# esets_nw$GSE16879__GPL570 <- esets_nw$GSE16879__GPL570[,sampleNames(esets_nw$GSE16879__GPL570) %in% 
#                                                          esets_nw$GSE16879__GPL570$sample_id[esets_nw$GSE16879__GPL570$tissue == 'colon' & esets_nw$GSE16879__GPL570$time == 'W0']]

esets_nw$GSE206285__GPL13158 <- inf_esets$GSE206285__GPL13158[,sampleNames(inf_esets$GSE206285__GPL13158) %in%
                                                                inf_esets$GSE206285__GPL13158$sample_id[inf_esets$GSE206285__GPL13158$sample_classification == 'No specific classification found']]

# esets_nw$`E-MTAB-7604__rnaseq` <- inf_esets$`E-MTAB-7604__rnaseq` [,sampleNames(inf_esets$`E-MTAB-7604__rnaseq` ) %in% 
#                                                                      inf_esets$`E-MTAB-7604__rnaseq` $sample_id[inf_esets$`E-MTAB-7604__rnaseq`$condition == 'inflammatory bowel disease']]

# esets_nw$GSE38713__GPL570 <- esets_nw$GSE38713__GPL570[,sampleNames(esets_nw$GSE38713__GPL570 ) %in% 
#                                                          esets_nw$GSE38713__GPL570 $sample_id[esets_nw$GSE38713__GPL570$Condition_comment == 'active disease']]

# esets_nw$GSE193677__GPL16791 <- esets_nw$GSE193677__GPL16791[,sampleNames(esets_nw$GSE193677__GPL16791 ) %in%
#                                                                esets_nw$GSE193677__GPL16791 $sample_id[esets_nw$GSE193677__GPL16791$Tissu_type == 'Colon']]

esets_nw$GSE206171__GPL13667 <- esets_nw$GSE206171__GPL13667[,sampleNames(esets_nw$GSE206171__GPL13667 ) %in%
                                                               esets_nw$GSE206171__GPL13667 $sample_id[esets_nw$GSE206171__GPL13667$Condition_comment == 'active disease']]

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
net_id <- "UC_colon_inflamed_NW_10_datasets"
IMAGE <- "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest"
default_memory_request <- "200Gi"

##### Step 1 - generate full correlation matrix
full_corr_matrix_res <- get_full_corr_matrix(esets = esets_nw,
                                             group_table = group_table,
                                             IMAGE = IMAGE,
                                             net_id = net_id,
                                             gpu = TRUE,
                                             mem_request = default_memory_request)



##### Step 2 - check outliers, if any is above 0.2 - remove and repeat step 1

# outliers <- check_outliers(reference_edges = "wf-ab1a4d302d", image_url = "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest", mem_request="50Gi", cyto_cc_force_execution=FALSE)
outliers <- check_outliers_python(reference_edges = full_corr_matrix_res[["reference_edges_wid"]])
corrplot::corrplot(as.matrix(outliers), order = "hclust", addCoef.col = "black",
                   addCoefasPercent=FALSE, method="color",
                   col=rev(RColorBrewer::brewer.pal(n=10, name="RdBu")))
#######################################################################################################

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


get_workflow(pruning_info_wfid, wait = TRUE)
pruning_info <- get_outputs_dist(pruning_info_wfid)$output.rds$edge_pruning_info
p_nodes(pruning_info)
p_edges(pruning_info)
p_components(pruning_info)
p_paretofront(pruning_info)
if (!is.null(important_genes)){
  p_important_genes(pruning_info)
}


##### Step 5
# use case 1: postprocess edge meta output
nets = consume_workflow_task_output(wf = edge_meta_res[["edge_meta"]][["wf_id"]], 
                                    task_id = NULL,
                                    task_name = edge_meta_res[["edge_meta"]][["task_name"]])

# use case 2: postprocess on previous orchestrated_edge_pruning results

# See orchestrated_edge_pruning doc for details regarding parameters 
prune <- c("min_abs_corr" = 0.45, "max_fdr" = 0.05, "min_q_fdr" = 0.05)
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

# https://cyto-cc.cytoreason.com/workflow/wf-5ee339dc90

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

get_workflow(nw_measures, wait = TRUE)
nw_measures_res <- readRDS(get_task_outputs(nw_measures, task_id = 0, files_names_grepl_pattern = "output.rds"))
head(nw_measures_res)
write.csv(nw_measures_res, '~/UC_NW_2025/outputs/inf_colon_nw_measures_res_0.45.csv')



edges <- get_outputs_dist(post_processes_nw)$output.rds
edges$Var1 <- as.numeric(edges$Var1)
edges$Var2 <- as.numeric(edges$Var2)
shortest_path_matrix_wfid <- calc_distance_matrix_gpu(edges = edges, force_abs = TRUE, weight_col = "cor", is_correlation = TRUE)
shortest_path_matrix <- get_distance_matrix_gpu(wfid = shortest_path_matrix_wfid)

write.csv(shortest_path_matrix, '~/UC_NW_2025/outputs/inf_colon_shortest_path_matrix_0.45.csv')

