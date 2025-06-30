
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
# check
inf_datasets %in% names(UC_model$datasets)
# subset to wanted datasets
inf_esets <- .get_esets_ccm(UC_model, esets_names = c(inf_datasets))
inf_esets$GSE53306__GPL14951 <- NULL # Allprep RNA/DNA FFPE kit protocol, outlier

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

esets_nw$GSE16879__GPL570 <- esets_nw$GSE16879__GPL570[,sampleNames(esets_nw$GSE16879__GPL570) %in% 
                                                         esets_nw$GSE16879__GPL570$sample_id[esets_nw$GSE16879__GPL570$tissue == 'colon' & esets_nw$GSE16879__GPL570$time == 'W0']]

esets_nw$GSE206285__GPL13158 <- inf_esets$GSE206285__GPL13158[,sampleNames(inf_esets$GSE206285__GPL13158) %in% 
                                                                inf_esets$GSE206285__GPL13158$sample_id[inf_esets$GSE206285__GPL13158$sample_classification == 'No specific classification found']]

esets_nw$`E-MTAB-7604__rnaseq` <- inf_esets$`E-MTAB-7604__rnaseq` [,sampleNames(inf_esets$`E-MTAB-7604__rnaseq` ) %in% 
                                                                     inf_esets$`E-MTAB-7604__rnaseq` $sample_id[inf_esets$`E-MTAB-7604__rnaseq`$condition == 'inflammatory bowel disease']]

esets_nw$GSE38713__GPL570 <- esets_nw$GSE38713__GPL570[,sampleNames(esets_nw$GSE38713__GPL570 ) %in% 
                                                         esets_nw$GSE38713__GPL570 $sample_id[esets_nw$GSE38713__GPL570$Condition_comment == 'active disease']]

esets_nw$GSE193677__GPL16791 <- esets_nw$GSE193677__GPL16791[,sampleNames(esets_nw$GSE193677__GPL16791 ) %in% 
                                                               esets_nw$GSE193677__GPL16791 $sample_id[esets_nw$GSE193677__GPL16791$Tissu_type == 'Colon']]

esets_nw$GSE193677__GPL16791 <- esets_nw$GSE193677__GPL16791[,sampleNames(esets_nw$GSE193677__GPL16791 ) %in% 
                                                               esets_nw$GSE193677__GPL16791 $sample_id[esets_nw$GSE193677__GPL16791$severity_condition %in% c('Moderate','Severe')]]


esets_nw$GSE206171__GPL13667 <- esets_nw$GSE206171__GPL13667[,sampleNames(esets_nw$GSE206171__GPL13667 ) %in% 
                                                               esets_nw$GSE206171__GPL13667 $sample_id[esets_nw$GSE206171__GPL13667$Condition_comment == 'active disease']]

esets_nw$GSE73661__GPL6244 <- esets_nw$GSE73661__GPL6244[,sampleNames(esets_nw$GSE73661__GPL6244 ) %in% 
                                                           esets_nw$GSE73661__GPL6244 $sample_id[esets_nw$GSE73661__GPL6244$time == 'W0']]

esets_nw$GSE95437__rnaseq <- esets_nw$GSE95437__rnaseq[,sampleNames(esets_nw$GSE95437__rnaseq ) %in% 
                                                               esets_nw$GSE95437__rnaseq $sample_id[esets_nw$GSE95437__rnaseq$Severity %in% c('moderate','severe')]]



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
net_id <- "inflamed_no_mild_UC"
IMAGE <- "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest"
default_memory_request <- "200Gi"

##### Step 1 - generate full correlation matrix
full_corr_matrix_res <- get_full_corr_matrix(esets = esets_nw,
                                             group_table = group_table,
                                             IMAGE = IMAGE,
                                             net_id = net_id,
                                             gpu = TRUE,
                                             mem_request = default_memory_request)

# https://cyto-cc.cytoreason.com/workflow/wf-e57066d676

##### Step 2 - check outliers, if any is above 0.2 - remove and repeat step 1

outliers <- check_outliers(reference_edges = "wf-e57066d676", image_url = "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest", mem_request="20Gi", cyto_cc_force_execution=FALSE)
corrplot::corrplot(as.matrix(outliers), order = "hclust", addCoef.col = "black",
                   addCoefasPercent=TRUE, method="color",
                   col=rev(RColorBrewer::brewer.pal(n=10, name="RdBu")))

##### Rerun Step 1
# we found 5 outliers ao we will rerun step 1 without them:
outlier_datasets <- c('GSE95437', 'GSE9452', 'GSE38713', 'GSE107593', 'GSE179285')
esets_nw_final <- esets_nw
for (dataset_num in 1:length(names(esets_nw_final))){
  if (experimentID(esets_nw_final[[dataset_num]]) %in% outlier_datasets){
    esets_nw_final[[dataset_num]] <- NULL
  }
    
}
esets_nw_final$GSE95437__rnaseq <- NULL
esets_nw_final$GSE73661__GPL6244 <- NULL # dup

group_table_final <- group_table
group_table_final <- group_table_final[-c(21,22,5,8,13,17),]
group_table_final$experiment_id %in% names(esets_nw_final)

net_id <- "UC_16_datasets_no_outliers"
full_corr_matrix_res <- get_full_corr_matrix(esets = esets_nw_final,
                                             group_table = group_table_final,
                                             IMAGE = IMAGE,
                                             net_id = net_id,
                                             gpu = TRUE,
                                             mem_request = default_memory_request)
# https://cyto-cc.cytoreason.com/workflow/wf-b9e52e3ba9

##### Step 3 -  Calculate meta gene-gene correlation matrix
edge_min_count <- 0
# if chunk_size is NULL, a default value (1000) will be used
chunk_size <- NULL
edge_meta_res <- meta_pvalue_edge(nets = full_corr_matrix_res[["reference_edges_wid"]],
                                  edge_min_count = edge_min_count,
                                  IMAGE = IMAGE,
                                  mem_req = default_memory_request,
                                  chunk_size = chunk_size)
# https://cyto-cc.cytoreason.com/workflow/wf-9654ac7614
get_workflow(edge_meta_res[["edge_meta"]][["wf_id"]], wait = TRUE)

all_edges <- cytoreason.io::read_data(cytoreason.cc.client::get_task_outputs(res=edge_meta_res[["edge_meta"]][["wf_id"]], task_id = NULL, task_name = edge_meta_res[["edge_meta"]][["task_name"]]))
q_pvalue_plot <- function (res) 
{
  ecdf_func <- ecdf(res[["q_pvalue"]])
  p <- plot(function(x) (1 - ecdf_func(x)) * nrow(res), xlim = c(0, 
                                                                 1), xlab = "q_pvalue", ylab = "Number of edges above", 
            main = "Number of edges above q_pvalue")
}

FDR_plot <- function (res) 
{
  ecdf_func <- ecdf(res[["q_pvalue_fdr"]])
  p <- plot(function(x) (1 - ecdf_func(x)) * nrow(res), xlim = c(0, 
                                                                 1), xlab = "q_pvalue_fdr", ylab = "Number of edges above", 
            main = "Number of edges above q_pvalue_fdr")
}


q_pvalue_plot(all_edges)
FDR_plot(all_edges)



##### Step 4 - Perform edge pruning analysis on meta gene-gene correlation matrix
edge_pruning_corr <- seq(0, 1, 0.05)
fdr_thresholds <- list(c(fdr_th=0.05, q_fdr_th=0.02))
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

# https://cyto-cc.cytoreason.com/workflow/wf-dd6a6ddf1e

get_workflow(pruning_info_wfid, wait = TRUE)
pruning_info <- get_outputs_dist(pruning_info_wfid)$output.rds$edge_pruning_info
p_nodes(pruning_info)
p_edges(pruning_info)
p_components(pruning_info)
p_paretofront(pruning_info)
if (!is.null(important_genes)){
  p_important_genes(pruning_info)
}
