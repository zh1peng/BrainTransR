

get_GO_data=getFromNamespace('get_GO_data','clusterProfiler')
geneSet_filter=getFromNamespace('geneSet_filter','DOSE')
build_Anno=getFromNamespace("build_Anno", "DOSE")
getGeneSet=getFromNamespace('getGeneSet','DOSE')
get_DGN_data=getFromNamespace('get_DGN_data','DOSE')
mapIds=getFromNamespace('mapIds','AnnotationDbi')
EXTID2NAME=getFromNamespace('EXTID2NAME','DOSE')
prepare_KEGG=getFromNamespace('prepare_KEGG','clusterProfiler')
prepare_WP_data=getFromNamespace('prepare_WP_data','clusterProfiler')
get_Reactome_DATA=getFromNamespace('get_Reactome_DATA','ReactomePA')
get_MeSH_data = getFromNamespace('get_MeSH_data', 'meshes')

get_geneSetList <- function(type=c('GO',
                                   'DO',
                                   'KEGG',
                                   'WikiPathways',
                                   'Reactome',
                                   'MeSH',
                                   'SynGO',
                                   'CellType'),
                              parameter){
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
    # convert entrezid to symbol
    geneSetList <- lapply(geneSetList, function(x) entrezid2symbol(x))
  } else if (type=='KEGG') {
    annoData <- prepare_KEGG('hsa', 'MKEGG', 'kegg')
    geneSetList <- getGeneSet(annoData)
    geneSetList <- lapply(geneSetList, function(x) entrezid2symbol(x))
  } else if (type=='WikiPathways'){
    wpdata <- prepare_WP_data('Homo sapiens')
    TERM2GENE = wpdata$WPID2GENE
    TERM2NAME = wpdata$WPID2NAME
    annoData = build_Anno(TERM2GENE, TERM2NAME)
    geneSetList <- getGeneSet(annoData)
    geneSetList <- lapply(geneSetList, function(x) entrezid2symbol(x))
  } else if(type=='Reactome'){
    annoData=get_Reactome_DATA('human')
    geneSetList <- getGeneSet(annoData)
  }
  return(geneSetList)
}


# From BioC 3.14 (Nov. 2021, with R-4.2.0)
library(AnnotationHub)
library(MeSHDbi)
ah <- AnnotationHub(localHub=TRUE)
hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
file_hsa <- hsa[[1]]
db <- MeSHDbi::MeSHDb(file_hsa)



#' Filter Gene Set List
#'
#' This function filters a list of gene sets based on the background genes and specified size constraints.
#'
#' @param bg_genes A vector of background gene symbols to be used for filtering.
#' @param geneSetList A list of gene sets to be filtered.
#' @param minGSSize Minimum gene set size for filtering.
#' @param maxGSSize Maximum gene set size for filtering.
#' @return A filtered list of gene sets that meet the size constraints and background genes criteria.
#' @importFrom DOSE geneSet_filter
#' @export
filter_geneSetList <- function(bg_genes, geneSetList, minGSSize, maxGSSize) {
  # Create temporary values to name the background genes
  tmp.val <- seq_along(bg_genes)
  names(tmp.val) <- bg_genes
  
  # Filter the gene set list
  geneSetList_filtered <- geneSet_filter(geneSetList, tmp.val, minGSSize, maxGSSize)
  
  return(geneSetList_filtered)
}




#' Convert Entrez IDs to Gene Symbols
#'
#' This function converts a vector of Entrez IDs to gene symbols using the org.Hs.eg.db annotation package.
#'
#' @param entrezid A vector of Entrez IDs to be converted to gene symbols.
#' @return A vector of gene symbols corresponding to the input Entrez IDs.
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @examples
#' \dontrun{
#' entrez_ids <- c("2451", "3142", "66666")
#' gene_symbols <- entrezid2symbol(entrez_ids)
#' print(gene_symbols)
#' }
#'
#' @export
entrezid2symbol <- function(entrezid) {
  # Ensure input is character vector
  entrezid <- as.character(entrezid)
  
  # Map Entrez IDs to gene symbols
  mappedSymbol <- suppressMessages(
    mapIds(
      x = org.Hs.eg.db,
      keys = entrezid,
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )
  )
  
  # Remove NA values and attributes
  mappedSymbol <- na.omit(mappedSymbol)
  attributes(mappedSymbol) <- NULL
  
  return(unname(mappedSymbol))
}



#' Download and process SynGO data
#'
#' This function downloads the SynGO data from the specified URL, extracts the `syngo_ontologies.xlsx` file,
#' processes the data to create TERM2GENE and TERM2NAME data frames, builds the annotation using `build_Anno`,
#' and removes the downloaded and temporary files afterward.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Downloads the SynGO data zip file from the provided URL.
#'   \item Extracts the `syngo_ontologies.xlsx` file from the zip archive.
#'   \item Reads the Excel file and processes the data to create TERM2GENE and TERM2NAME data frames.
#'   \item Uses the `build_Anno` function to build the final annotation data.
#'   \item Cleans up by removing the temporary zip file and extracted Excel file.
#' }
#'
#' @param url The URL to download the SynGO data zip file from. Default is the latest release URL.
#' @return A data frame containing the processed SynGO data.
#' @import readxl
#' @importFrom dplyr select rename mutate %>%
#' @importFrom tidyr unnest
#' @importFrom utils download.file unzip
#' @importFrom DOSE build_Anno
#'
#' @examples
#' \dontrun{
#' # Example usage of get_SynGO_data function
#' data <- get_SynGO_data()
#' print(head(data))  # Prints the first few rows of the processed data
#' }
#'
#' @export
get_SynGO_data <- function(url = "https://www.syngoportal.org/data/SynGO_bulk_download_release_20231201.zip") {
  # Create temporary file and directory
  zip_path <- tempfile(fileext = ".zip")
  temp_dir <- tempdir()
  extracted_file <- file.path(temp_dir, "syngo_ontologies.xlsx")

  # Clean-up function to ensure temp files are removed
  on.exit({
    if (file.exists(zip_path)) unlink(zip_path)
    if (file.exists(extracted_file)) unlink(extracted_file)
  }, add = TRUE)

  tryCatch({
    # Download the file
    message("Downloading SynGO data...")
    download.file(url, zip_path, mode = "wb")

    # Unzip the required file
    message("Unzipping SynGO data...")
    unzip(zip_path, files = "syngo_ontologies.xlsx", exdir = temp_dir)

    # Read the data from the extracted file
    if (!file.exists(extracted_file)) {
      stop("Failed to extract the required file.")
    }

    data <- read_xlsx(extracted_file)

    # Process TERM2GENE
    TERM2GENE <- data %>%
      select(id, hgnc_symbol) %>%
      rename(cLabel = id, geneID = hgnc_symbol) %>%
      mutate(geneID = strsplit(geneID, ";")) %>%
      unnest(cols = c(geneID))

    # Process TERM2NAME
    TERM2NAME <- data %>%
      select(id, name) %>%
      rename(cLabel = id, description = name)

    # Use build_Anno function from DOSE package
    USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)

    message("SynGO data has been processed and temporary files removed.")
    return(USER_DATA)
  }, error = function(e) {
    stop("An error occurred while processing SynGO data: ", e$message)
  })
}

#' Download and process cell type gene set data
#'
#' This function downloads cell type gene set data from specified resources, processes the data,
#' and returns a list of gene sets based on the specified type.
#'
#' @details
#' The function performs the following steps based on the specified type:
#' \itemize{
#'   \item Downloads the cell type gene set data from the provided URL.
#'   \item Reads the data from the file and processes it to create a gene set list.
#'   \item Supports the following types:
#'     \itemize{
#'       \item \strong{Seidlitz2020}: \url{https://github.com/jms290/PolySyn_MSNs/blob/master/Data/AHBA/celltypes_PSP.csv}
#'       \itemize{
#'         \item Seidlitz, J., Nadig, A., Liu, S., Bethlehem, R. A., Vértes, P. E., Morgan, S. E., ... & Raznahan, A. (2020). Transcriptomic and cellular decoding of regional brain vulnerability to neurogenetic disorders. Nature communications, 11(1), 3358.
#'       }
#'       \item \strong{Lake2018}: \url{https://github.com/molecular-neuroimaging/Imaging_Transcriptomics/blob/main/imaging_transcriptomics/data/geneset_LAKE.gmt}
#'       \itemize{
#'         \item Lake, B. B., Chen, S., Sos, B. C., Fan, J., Kaeser, G. E., Yung, Y. C., ... & Zhang, K. (2018). Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain. Nature biotechnology, 36(1), 70-80.
#'       }
#'       \item \strong{Martins2021}: \url{https://github.com/molecular-neuroimaging/Imaging_Transcriptomics/blob/main/imaging_transcriptomics/data/geneset_Pooled.gmt}
#'       \itemize{
#'         \item Imaging transcriptomics: Convergent cellular, transcriptomic, and molecular neuroimaging signatures in the healthy adult human brain. Daniel Martins, Alessio Giacomel, Steven CR Williams, Federico Turkheimer, Ottavia Dipasquale, Mattia Veronese, PET templates working group. Cell Reports.
#'       }
#'     }
#' }
#'
#' @param type The type of cell type gene set data to download and process. Must be one of 'Seidlitz2020', 'Lake2018', 'Martins2021'.
#' @return A list of gene sets where each element is a vector of genes associated with a specific term.
#' @importFrom dplyr filter rename
#' @importFrom utils download.file
#' @importFrom clusterProfiler read.gmt
#' @importFrom tibble as_tibble
#' @examples
#' \dontrun{
#' # Example usage of get_celltype_data function
#' geneSetList <- get_celltype_data('Seidlitz2020')
#' print(geneSetList)
#' }
#'
#' @export
get_celltype_data <- function(type = c('Seidlitz2020', 'Lake2018', 'Martins2021')) {
  type <- match.arg(type)
  
  # Define URLs for each type
  urls <- list(
    Seidlitz2020 = "https://github.com/jms290/PolySyn_MSNs/blob/master/Data/AHBA/celltypes_PSP.csv?raw=true",
    Lake2018 = "https://github.com/molecular-neuroimaging/Imaging_Transcriptomics/raw/main/imaging_transcriptomics/data/geneset_LAKE.gmt",
    Martins2021 = "https://github.com/molecular-neuroimaging/Imaging_Transcriptomics/raw/main/imaging_transcriptomics/data/geneset_Pooled.gmt"
  )

  # Create temporary file
  temp_file <- tempfile()
  
  # Clean-up function to ensure temp files are removed
  on.exit({
    if (file.exists(temp_file)) unlink(temp_file)
  }, add = TRUE)
  
  # Download the file
  message("Downloading cell type data...")
  download.file(urls[[type]], temp_file, mode = "wb")

  tryCatch({
    if (type == 'Seidlitz2020') {
      # Read CSV file
      df.glist <- read.csv(temp_file)
      df.glist <- df.glist %>% filter(!is.na(gene)) %>% rename(term = class)
    } else {
      # Read GMT file
      df.glist <- read.gmt(temp_file)
      df.glist <- df.glist %>% filter(!is.na(gene))
    }
    
    # Create gene set list
    geneSetList <- split(df.glist$gene, df.glist$term)

    message("Cell type data has been processed and temporary files removed.")
    return(geneSetList)
  }, error = function(e) {
    stop("An error occurred while processing cell type data: ", e$message)
  })
}




































