
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

esets_nw$GSE206171__GPL13667 <- esets_nw$GSE206171__GPL13667[,sampleNames(esets_nw$GSE206171__GPL13667 ) %in% 
                                                               esets_nw$GSE206171__GPL13667 $sample_id[esets_nw$GSE206171__GPL13667$Condition_comment == 'active disease']]

esets_nw$GSE73661__GPL6244 <- esets_nw$GSE73661__GPL6244[,sampleNames(esets_nw$GSE73661__GPL6244 ) %in% 
                                                           esets_nw$GSE73661__GPL6244 $sample_id[esets_nw$GSE73661__GPL6244$time == 'W0']]



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
net_id <- "first_test_UC"
IMAGE <- "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest"
default_memory_request <- "200Gi"

##### Step 1 - generate full correlation matrix
full_corr_matrix_res <- get_full_corr_matrix(esets = esets_nw,
                                             group_table = group_table,
                                             IMAGE = IMAGE,
                                             net_id = net_id,
                                             gpu = TRUE,
                                             mem_request = default_memory_request)

# https://cyto-cc.cytoreason.com/workflow/wf-8f453e4460

##### Step 2 - check outliers, if any is above 0.2 - remove and repeat step 1

outliers <- check_outliers(reference_edges = "wf-8f453e4460", image_url = "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest", mem_request="20Gi", cyto_cc_force_execution=FALSE)
corrplot::corrplot(as.matrix(outliers), order = "hclust", addCoef.col = "black",
                   addCoefasPercent=FALSE, method="color",
                   col=rev(RColorBrewer::brewer.pal(n=10, name="RdBu")))

##### Rerun Step 1
# we found 4 outliers ao we will rerun step 1 without them:
outlier_datasets <- c('GSE95437', 'GSE9452', 'GSE38713', 'GSE107593')
esets_nw_final <- esets_nw
for (dataset_num in 1:length(names(esets_nw_final))){
  if (experimentID(esets_nw_final[[dataset_num]]) %in% outlier_datasets){
    esets_nw_final[[dataset_num]] <- NULL
  }
    
}

group_table_final <- group_table
group_table_final <- group_table_final[-c(21,22,5,13),]
group_table_final$experiment_id %in% names(esets_nw_final)

net_id <- "UC_18_datasets_no_outliers"
full_corr_matrix_res <- get_full_corr_matrix(esets = esets_nw_final,
                                             group_table = group_table_final,
                                             IMAGE = IMAGE,
                                             net_id = net_id,
                                             gpu = TRUE,
                                             mem_request = default_memory_request)

