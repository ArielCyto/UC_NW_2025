
#install.packages("Polychrome")
library(Polychrome)
library(CHBUtils)
library(cytoreason.validator.apps.client)
library(cytoreason.ccm.pipeline)


inf_esets <- .get_esets_ccm(UC_model, esets_names = c(inf_datasets))

P36 = createPalette(36,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P36)

dataset <- inf_esets[[9]]
dataset <- dataset[,dataset$sample_classification == 'Inflamed']
dataset$plot_data <- paste(dataset$tissue_class,dataset$sample_classification)
PCAplot.eset(eset=dataset, categories="tissue_class",title=paste("PCAplot - inflamed tissues", inf_datasets[9]),colorpalette=P36[1:length(unique(dataset$plot_data))], alpha=0.8, numcomponents=2)        
PCAplot.eset(eset=dataset, categories="tissue_severity_condition",title=paste("PCAplot - inflamed tissues - severity", inf_datasets[9]),colorpalette=P36[1:length(unique(dataset$tissue_severity_condition))], alpha=0.8, numcomponents=2)             

dataset_ann <- cytoreason.curator.annotations::get_annotations('GSE193677')
dataset$Severity <- NA
for (sample_id in unique(dataset$geo_accession)){
  dataset$Severity[dataset$geo_accession == sample_id] <- dataset_ann$severity[dataset_ann$sample_id == sample_id]
}
PCAplot.eset(eset=dataset, categories="Severity",title=paste("PCAplot - inflamed tissues - severity", inf_datasets[9]),colorpalette=P36[1:length(unique(dataset$Severity))], alpha=0.8, numcomponents=2)             
dataset$plot_data <- paste(dataset$tissue_class,dataset$Severity)
PCAplot.eset(eset=dataset, categories="plot_data",title=paste("PCAplot - inflamed tissues - severity", inf_datasets[9]),colorpalette=P36[1:length(unique(dataset$plot_data))], alpha=0.8, numcomponents=2)             
dataset_pruned <- dataset
dataset_pruned <- dataset_pruned[,dataset_pruned$Severity %in% c('moderate', 'severe')]
PCAplot.eset(eset=dataset_pruned, categories="plot_data",title=paste("PCAplot - inflamed - moderate and severe only", inf_datasets[9]),colorpalette=P36[1:length(unique(dataset_pruned$plot_data))], alpha=0.8, numcomponents=2)             




dataset <- inf_esets[[16]]
dataset <- dataset[,dataset$sample_classification == 'Inflamed']
dataset$plot_data <- paste(dataset$`tissue:ch1`,dataset$sample_classification)
PCAplot.eset(eset=dataset, categories="plot_data",title=paste("PCAplot - inflamed tissues", inf_datasets[6]),colorpalette=P36[1:length(unique(dataset$plot_data))], alpha=0.8, numcomponents=2)        
#PCAplot.eset(eset=dataset, categories="severity_condition",title=paste("PCAplot - inflamed tissues - severity", inf_datasets[6]),colorpalette=P36[1:length(unique(dataset$severity_condition))], alpha=0.8, numcomponents=2)             

dataset <- inf_esets[[23]]
dataset <- dataset[,dataset$sample_classification == 'Inflamed']
PCAplot.eset(eset=dataset, categories="Severity",title=paste("PCAplot - inflamed tissues - severity", inf_datasets[23]),colorpalette=P36[1:length(unique(dataset$Severity))], alpha=0.8, numcomponents=2)             
PCAplot.eset(eset=dataset, categories="Severity",title=paste("PCAplot - inflamed tissues - severity", inf_datasets[23]),colorpalette=P36[1:length(unique(dataset$Severity))], alpha=0.8, numcomponents=2)