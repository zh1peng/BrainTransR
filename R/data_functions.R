# load data functions
library(dplyr)
#================================================================================================== 
#Functions to load relevant data
#==================================================================================================

load_GeneExp <- function (data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneExp',
                          atlas=c('desikan',
                                  'schaefer100',
                                  'schaefer200',
                                  'schaefer300'),
                          rdonor=c('r0.2',
                                    'r0.4',
                                    'r0.6'),
                          hem=c('L','R','B')) {
  atlas=match.arg(atlas)
  rdonor=match.arg(rdonor)
  hem=match.arg(hem)

  if (atlas=='desikan'){
  search_pattern=sprintf('^%s_',hem)
  }else if (atlas %in% c('schaefer100','schaefer200','schaefer300')){
    search_pattern=sprintf('%sH', hem)
  }

# get csv file of the gene expression data
  GeneExpCSV=sprintf('%s/%s_%s.csv',data_path, atlas,rdonor)
# check if csvfile exist otherwise stop
  if (!file.exists(GeneExpCSV)) {
    stop(sprintf('GeneExp file %s does not exist',GeneExpCSV))
  }
gene.df=read.csv(GeneExpCSV, stringsAsFactors = F)
 if (hem=='L'|hem=='R'){
    gene.df=gene.df %>% filter(grepl(search_pattern,Region))
  }                
  gene.df=gene.df %>% filter(complete.cases(.)) %>% tibble::column_to_rownames('Region')
  gene.mx=as.matrix(gene.df)
  return(gene.mx)
}


load_BrainDat <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainDat',
                          atlas=c('desikan',
                                   'schaefer100',
                                   'schaefer200'),
                           type=c('sample_data',
                                  'sim_nospatial',
                                  'sim_spatial0.03',
                                  'sim_spatial0.02',
                                  'sim_spatial0.01',
                                  'real_brain'),
                          col_idx=1){
  atlas=match.arg(atlas)
  type=match.arg(type)
  BrainDatCSV=sprintf('%s/%s_%s.csv',data_path,atlas,type)
  
  #check if BrainDatFile exist
  if (!file.exists(BrainDatCSV)) {
    stop(sprintf('BrainDat file %s does not exist',BrainDatCSV))
  }

  df=read.csv(BrainDatCSV,row.names = 1) 
  df.mx=as.matrix(df)
  # if col_idx is NA, use all columns
  if (col_idx=='all'){col_idx=c(1:dim(df.mx)[2])} 

  # select columns
  df.mx=df.mx[,col_idx,drop=F] 

  # if only one column, rename it as 'phenotype'
  if (length(col_idx)==1){
  colnames(df.mx)='phenotype'
  }
  return(df.mx)
}


load_GeneSets <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneSets',
                          atlas=c('desikan',
                                  'schaefer100',
                                  'schaefer200'),
                            rdonor=c('r0.2',
                                      'r0.4',
                                      'r0.6'),
                            gs_type=c('MF',
                                      'Sim',
                                      'SynGO')){
  atlas=match.arg(atlas)
  rdonor=match.arg(rdonor)
  gs_type=match.arg(gs_type)

  GeneSetsRDS=sprintf('%s/%s_%s_%s.rds',data_path, atlas,rdonor,gs_type)

  #check if GeneSetsRDS exist
  if (!file.exists(GeneSetsRDS)) {
    stop(sprintf('GeneSetsRDS file %s does not exist',GeneSetsRDS))
  }
geneSetList=readRDS(GeneSetsRDS)
return(geneSetList)
}




load_sampled_GeneSets <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/sampled_GeneSets',
                          atlas=c('desikan',
                                  'schaefer100',
                                  'schaefer200'),
                            rdonor=c('r0.2',
                                      'r0.4',
                                      'r0.6'),
                            gs_type=c('MF',
                                      'Sim',
                                      'SynGO'),
                            constrain=c('match_coexp','gene_subset')){
  atlas=match.arg(atlas)
  rdonor=match.arg(rdonor)
  gs_type=match.arg(gs_type)

  GeneSetsRDS=sprintf('%s/%s_%s_%s_%s.rds',data_path, atlas,rdonor,gs_type,constrain)

  #check if GeneSetsRDS exist
  if (!file.exists(GeneSetsRDS)) {
    stop(sprintf('GeneSetsRDS file %s does not exist',GeneSetsRDS))
  }
geneSetList=readRDS(GeneSetsRDS)
return(geneSetList)
}