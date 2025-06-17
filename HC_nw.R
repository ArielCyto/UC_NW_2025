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
