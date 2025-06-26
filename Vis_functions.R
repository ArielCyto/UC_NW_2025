### Visualization functions fron Gil's code

get_node_positions <- function(network) {
  layout <- layout_with_fr(network)  # Fruchterman-Reingold layout for better visualization
  positions <- data.frame(id = names(V(network)), x = layout[,1], y = layout[,2])
  return(positions)
}

visualize_gene_set_subgraphs <- function(whole_network, gene_sets, centrality_df, centrality_measure = "betweenness") {
  for (i in seq_along(gene_sets)) {
    gene_set <- gene_sets[[i]]
    gene_set <- gene_set[gene_set %in% names(V(whole_network))]
    
    node_positions <- get_node_positions(whole_network)
    
    vids_df <- which(V(whole_network)$name %in% gene_set)
    # Create subgraph for the gene set
    subgraph <- induced_subgraph(whole_network, vids = vids_df)
    
    # Extract node names from the subgraph
    nodes <- V(subgraph)$name
    
    # Merge subgraph nodes with centrality measures
    subgraph_centrality <- centrality_df %>%
      dplyr::filter(feature_id %in% nodes) %>%
      dplyr::arrange(feature_id)
    
    # Prepare data for visNetwork
    nodes_df <- data.frame(
      id = subgraph_centrality$feature_id,
      label = subgraph_centrality$feature_id,
      #color = scales::col_numeric("viridis", domain = NULL)(subgraph_centrality[[centrality_measure]])
      color = scales::col_numeric(palette = c("#f0e3ab", "#f49d77", "#ce5f5f", "#b23030"), domain = NULL)(subgraph_centrality[[centrality_measure]]),
      font.size = 32  # Increase font size
    )
    
    # Add positions to nodes_df
    nodes_df <- merge(nodes_df, node_positions, by = "id", all.x = TRUE)
    
    ann_df_gene_name <- mapIds(keys=as.character(nodes_df$id), x = org.Hs.eg.db, keytype = c("ENTREZID"), column = c("SYMBOL")) %>%
      tibble::enframe(name = "ENTREZID", value = "SYMBOL")
    ann_df_gene_name <- data.frame(ann_df_gene_name, stringsAsFactors = F)
    
    nodes_df$label <- ann_df_gene_name[match(nodes_df$id, ann_df_gene_name$ENTREZID),'SYMBOL']
    
    edges_df <- igraph::as_data_frame(subgraph, what = "edges")
    
    #edges_df$from <- ann_df_gene_name[match(edges_df$from, ann_df_gene_name$ENTREZID),'SYMBOL']
    #edges_df$to <- ann_df_gene_name[match(edges_df$to, ann_df_gene_name$ENTREZID),'SYMBOL']
    
    vis_net <- visNetwork(nodes_df, edges_df) %>%
      visNodes(font = list(size = 20)) %>%  # Increase font size
      visEdges(color = list(color = "gray")) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visLayout(randomSeed = 123) %>%
      visPhysics(stabilization = FALSE, 
                 solver = "barnesHut",
                 barnesHut = list(gravitationalConstant = -8000,
                                  centralGravity = 0.3,
                                  springLength = 100,
                                  springConstant = 0.04,
                                  damping = 0.09))  # Disable physics movement
    
    clean_name <- gsub("[^A-Za-z0-9_]", "",  names(gene_sets)[i])
    
    # Save visualization as HTML file
    file_name <- paste0("Gene_Set_", clean_name, "_", centrality_measure, ".html")
    visSave(vis_net, file = file_name)
    cat("Saved visualization to", file_name, "\n")
  }
}


customSignatures <- read_asset("wf-13096e7be3")
combo_genesets <- customSignatures$dupiMono
mono_genesets <- customSignatures$monotherapies