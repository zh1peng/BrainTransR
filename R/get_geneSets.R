



#' Retrieve Gene Sets Based on GO Terms
#'
#' This function fetches gene sets from the Gene Ontology (GO) data, 
#' filters them based on correlation with synthetic brain data, and applies size constraints.
#'
#' @param gene_data A matrix or data frame of gene expression data.
#' @param ont2use A character vector specifying the ontology categories to use: 'MF' (Molecular Function),
#'        'BP' (Biological Process), or 'CC' (Cellular Component).
#' @param minGSSize Minimum size of gene sets to consider.
#' @param maxGSSize Maximum size of gene sets to consider.
#'
#' @return A list of gene sets filtered based on the specified size and correlation criteria.
#' @importFrom clusterProfiler get_GO_data
#' @importFrom DOSE geneSet_filter getGeneSet
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom stats cor
#' @export
#'
#' @examples
#' gene_data <- matrix(rnorm(1000), ncol=10)
#' gene_sets <- get_GO_geneSetList(gene_data)
get_GO_geneSetList <- function(gene_data,
                               ont2use = c('MF', 'BP', 'CC'),
                               minGSSize = 20,
                               maxGSSize = 200) {
  # Retrieve GO data using predefined ontology and symbol
  GO_env <- get_GO_data('org.Hs.eg.db', ont = ont2use, 'SYMBOL')
  USER_DATA <- getGeneSet(GO_env)

  # Create synthetic brain data and compute correlation with gene data
  brain_data <- rnorm(nrow(gene_data)) 
  geneList_tmp <- cor(gene_data, brain_data)

  # Filter gene sets based on correlation and size constraints
  geneSetList <- geneSet_filter(USER_DATA, geneList_tmp[, 1], minGSSize = minGSSize, maxGSSize = maxGSSize)

  return(geneSetList)
}



build_Anno=getFromNamespace("build_Anno", "DOSE")

get_synGO_geneSetList <- function(gene_data,
                                  synGO_file='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneSets/syngo_ontologies.xlsx',
                                  minGSSize = 20,
                                  maxGSSize = 200) {
#check if synGO_file exist
if (!file.exists(synGO_file)) {
  stop(sprintf('synGO_file %s does not exist',synGO_file))
}

GS_df=readxl::read_xlsx(synGO_file)
TERM2GENE=GS_df%>%dplyr::select(`GO term ID`, `genes - hgnc_symbol`) %>% # TERM2GENE
                  dplyr::rename(cLabel=`GO term ID`,
                                geneID=`genes - hgnc_symbol`) %>% 
                  mutate(geneID = strsplit(geneID, ';')) %>%
                  unnest(cols = c(geneID))  
TERM2NAME=GS_df%>%dplyr::select(`GO term ID`, `GO term name`) %>% 
                  dplyr::rename(cLabel=`GO term ID`,
                                description=`GO term name`)
USER_DATA = build_Anno(TERM2GENE,TERM2NAME)
USER_GS=getGeneSet(USER_DATA)
brain_data=rnorm(nrow(gene_data)) # create a fake brain data
geneList.tmp=cor(gene_data, brain_data)
geneSetList=geneSet_filter(USER_GS, geneList.tmp[,1], minGSSize=minGSSize,maxGSSize=maxGSSize)
return(geneSetList)
}


load_geneSetList <- function(gene_data,
                            type=c('Cell_Martins2021_lake',
                                    'Cell_Martins2021_pooled',
                                    'Cell_Seidlitz2020',
                                    'GO_MF',
                                    'GO_BP',
                                    'GO_CC',
                                   'GO_GSEA',
                                   'AUD_GWAS_GSEA',
                                   'GO_GSEA_MF'),
                            minGSSize=10,
                            maxGSSize=200){
type=match.arg(type)
  
  if (type=='Cell_Martins2021_lake') {
    # glist from Lake
    df.glist = clusterProfiler::read.gmt('F:/Google Drive/post-doc/GSVA/code/cell_types/geneset_LAKE.gmt')
    df.glist = df.glist %>% filter(!gene == '')
    geneSetList = split(df.glist$gene, df.glist$term)
  } else if (type == 'Cell_Martins2021_pooled') {
    # glist from pooled
    df.glist = clusterProfiler::read.gmt('F:/Google Drive/post-doc/GSVA/code/cell_types/geneset_Pooled.gmt')
    df.glist = df.glist %>% filter(!gene == '')
    geneSetList = split(df.glist$gene, df.glist$term)
  } else if (type == 'Cell_Seidlitz2020') {
    df.glist = read.csv('F:/Google Drive/post-doc/GSVA/code/cell_types/2020NC.txt')
    df.glist = df.glist %>% filter(!gene == '') %>% rename(term = class)
    geneSetList = split(df.glist$gene, df.glist$term)
  } else if (type == 'GO_BP') {
    # ont2use='BP'
    # go_env=get_GO_data('org.Hs.eg.db',ont2use,'SYMBOL')
    # geneSetList=getGeneSet(go_env)
    # saveRDS(geneSetList,'F:/Google Drive/post-doc/GSVA/code/GO/GO_BP.Rds')
    geneSetList = readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_BP.Rds')
  } else if (type == 'GO_CC') {
    # ont2use='CC'
    # go_env=get_GO_data('org.Hs.eg.db',ont2use,'SYMBOL')
    # geneSetList=getGeneSet(go_env)
    # saveRDS(geneSetList,'F:/Google Drive/post-doc/GSVA/code/GO/GO_CC.Rds')
    geneSetList = readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_CC.Rds')
  } else if (type == 'GO_MF') {
    # ont2use='MF'
    # go_env=get_GO_data('org.Hs.eg.db',ont2use,'SYMBOL')
    # geneSetList=getGeneSet(go_env)
    # saveRDS(geneSetList,'F:/Google Drive/post-doc/GSVA/code/GO/GO_MF.Rds')
    geneSetList = readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_MF.Rds')
  } else if (type == 'GO_GSEA') {
    # geneSetList_CC=readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_CC.Rds')
    # geneSetList_BP=readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_BP.Rds')
    # geneSetList_MF=readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_MF.Rds')
    # geneSetList_all=c(geneSetList_CC,geneSetList_BP,geneSetList_MF)
    # geneSetList=geneSetList_all[df.res$ID]
    # saveRDS(geneSetList,'F:/Google Drive/post-doc/GSVA/code/GO/GO_GSEA.Rds')
    geneSetList = readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_GSEA.Rds')
  } else if (type == 'GO_GSEA_MF'){
    # geneSetList_MF=readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_MF.Rds')
    # df.res=read.csv('F:/Google Drive/post-doc/Structural_subtype_new/all_res/GSEA/MF.csv')
    # geneSetList=geneSetList_MF[df.res$ID]
    # saveRDS(geneSetList,'F:/Google Drive/post-doc/GSVA/code/GO/GO_GSEA_MF.Rds')
    geneSetList=readRDS('F:/Google Drive/post-doc/GSVA/code/GO/GO_GSEA_MF.Rds')
  }
  else if (type == 'AUD_GWAS_GSEA') {
    df.glist = read.csv('F:/Google Drive/post-doc/GSVA/code/GWAS/magma_gsa_IDconverted.csv')
    df.glist = df.glist %>% filter(complete.cases(SYMBOL)) %>% rename(term = class)
    geneSetList = split(df.glist$SYMBOL, df.glist$term)
  }
  # df.glist=read.csv('F:/Google Drive/post-doc/GSVA/code/GWAS/GWAS.csv') %>% filter(journal=='NC')
  # geneSetList=split(df.glist$geneID,df.glist$term)
  
  # make a fake geneList to do the filter
  bg_genes=c(1:ncol(gene_data))
  names(bg_genes)=colnames(gene_data)
  geneSetList_filtered=geneSet_filter(geneSetList, bg_genes, minGSSize, maxGSSize)
  return(geneSetList_filtered)
}

get_GO_data=getFromNamespace('get_GO_data','clusterProfiler')
geneSet_filter=getFromNamespace('geneSet_filter','DOSE')
build_Anno=getFromNamespace("build_Anno", "DOSE")
getGeneSet=getFromNamespace('getGeneSet','DOSE')
get_DGN_data=getFromNamespace('get_DGN_data','DOSE')
mapIds=getFromNamespace('mapIds','AnnotationDbi')
EXTID2NAME=getFromNamespace('EXTID2NAME','DOSE')

EXTID2NAME(OrgDb='org.Hs.eg.db',geneID=c("1468", "4210","99999"),keytype="ENTREZID")


get_geneSetList <- function(type=c('GO_BP',
                                   'GO_MF',
                                   'GO_CC',
                                   'DO',
                                   'KEGG',
                                   'WikiPathways',
                                   'Reactome',
                                   'MeSH',
                                   'SynGO',
                                   'CellType')){
  type=match.arg(type)
  if (type=='GO_BP') {
    annoData <- get_GO_data('org.Hs.eg.db', ont = 'BP', 'SYMBOL')
    geneSetList <- getGeneSet(annoData)
  } else if (type=='GO_MF') {
    annoData <- get_GO_data('org.Hs.eg.db', ont = 'MF', 'SYMBOL')
    geneSetList <- getGeneSet(annoData)
  } else if (type=='GO_CC') {
    annoData <- get_GO_data('org.Hs.eg.db', ont = 'CC', 'SYMBOL')
    geneSetList <- getGeneSet(annoData)
  } else if (type=='DisGeNet') {
  annoData=get_DGN_data()
  geneSetList <- getGeneSet(annoData)
  geneSetNames <- get("PATHID2NAME", envir = annoData)
  # convert entrezid to symbol
  geneSetList <- lapply(geneSetList, function(x) entrezid2symbol(x))
  


  } else if (type=='KEGG') {


  }
  return(geneSetList)
  }



fitler_geneSetList <- function(bg_genes, geneSetList, minGSSize, maxGSSize) {
  # Create synthetic brain data and compute correlation with gene data
bg_genes=c(1:ncol(gene_data))
  names(bg_genes)=colnames(gene_data)
  geneSetList_filtered=geneSet_filter(geneSetList, bg_genes, minGSSize, maxGSSize)
  return(geneSetList_filtered)
}




entrezid2symbol <- function(entrezid) {
  # Convert the input to character to ensure compatibility with mapIds
  entrezid <- as.character(entrezid)

  # Corrected function call
  mappedSymbol <- suppressMessages(
    mapIds(x = org.Hs.eg.db,
           keys = entrezid,
           keytype = "ENTREZID",
           column = "SYMBOL",
           multiVals = "first")
  )
  mappedSymbol <- na.omit(mappedSymbol)
  # Remove all attributes from the result
    attributes(mappedSymbol) <- NULL

  return(unname(mappedSymbol))
}



prepare_WP_data
prepare_KG_data






































