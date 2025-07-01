edge_pruning_corr <- seq(0, 1, 0.05)
fdr_thresholds <- list(c(fdr_th=0.05, q_fdr_th=0.07))
important_genes <- NULL
pruning_info_wfid <- 
  run_function_dist(function(full_corr_matrix_res, ...) {
    pruning_info_res <- cytoreason.individual.variation::edge_pruning_th_comb(...)
    return(pruning_info_res)
  },
  nets = 'wf-d896d61ba5' , # 'wf-9654ac7614', 
  nets_task_name = 'meta-pvalue-post',
  cor_th = edge_pruning_corr, 
  fdr_ths = fdr_thresholds, 
  full_corr_matrix_res = 'wf-b9e52e3ba9',
  IMAGE = IMAGE,
  mem_req = default_memory_request,
  cytocc = TRUE,
  plot = FALSE, 
  important_genes = important_genes,
  image = IMAGE,
  memory_request = default_memory_request,
  force_execution = FALSE, replace_image_tags = TRUE, 
  tags = list(list("name" = "test_010725", "value" =  "edge_pruning_analysis")))
#  wf-6c50137329 - test with yasmin
#  wf-fe4962f065 - positive cor only