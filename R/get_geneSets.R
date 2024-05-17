#' Get Gene Set List based on type and parameters
#'
#' This function retrieves a gene set list based on the specified type and additional parameters.
#'
#' @details
#' The function supports various types of gene sets, including 'GO', 'KEGG', 'WikiPathways', 'Reactome', 'SynGO', and 'CellType'.
#' For 'GO', the additional parameter can be 'BP', 'MF', or 'CC'.
#' For 'CellType', the additional parameter can be 'Seidlitz2020', 'Lake2018', or 'Martins2021'.
#'
#' @param type The type of gene set data to retrieve. Must be one of 'GO', 'KEGG', 'WikiPathways', 'Reactome', 'SynGO', 'CellType'.
#' @param parameter Additional parameter for specific types like 'GO' and 'CellType'. For 'GO', use 'BP', 'MF', or 'CC'. For 'CellType', use 'Seidlitz2020', 'Lake2018', or 'Martins2021'.
#' @return A list of gene sets where each element is a vector of genes associated with a specific term.
#' @examples
#' \dontrun{
#' # Example usage of get_geneSetList function
#' geneSetList <- get_geneSetList('GO', 'BP')
#' print(geneSetList)
#' }
#'
#' @export
get_geneSetList <- function(type = c('GO', 'KEGG', 'WikiPathways', 'Reactome', 'SynGO', 'CellType'),
                            parameter = NULL) {
  type <- match.arg(type)
  
  # Main processing based on type
  annoData <- NULL
  convert_to_symbol <- FALSE
  
  if (type == 'GO') {
    if (is.null(parameter) || !parameter %in% c('BP', 'MF', 'CC')) {
      stop("For type 'GO', parameter must be one of 'BP', 'MF', or 'CC'.")
    }
    get_GO_data=getFromNamespace('get_GO_data','clusterProfiler')
    annoData <- get_GO_data('org.Hs.eg.db', ont = parameter, 'SYMBOL')
  } else if (type == 'DisGeNet') {
    get_DGN_data=getFromNamespace('get_DGN_data','DOSE')
    annoData <- get_DGN_data()
    convert_to_symbol <- TRUE
  } else if (type == 'KEGG') {
    prepare_KEGG=getFromNamespace('prepare_KEGG','clusterProfiler')
    annoData <- prepare_KEGG('hsa', 'MKEGG', 'kegg')
    convert_to_symbol <- TRUE
  } else if (type == 'WikiPathways') {
    prepare_WP_data=getFromNamespace('prepare_WP_data','clusterProfiler')
    wpdata <- prepare_WP_data('Homo sapiens')
    TERM2GENE <- wpdata$WPID2GENE
    TERM2NAME <- wpdata$WPID2NAME
    build_Anno=getFromNamespace("build_Anno", "DOSE")
    annoData <- build_Anno(TERM2GENE, TERM2NAME)
    convert_to_symbol <- TRUE
  } else if (type == 'Reactome') {
    get_Reactome_DATA=getFromNamespace('get_Reactome_DATA','ReactomePA')
    annoData <- get_Reactome_DATA('human')
    convert_to_symbol <- TRUE
  } else if (type == 'CellType') {
    if (is.null(parameter) || !parameter %in% c('Seidlitz2020', 'Lake2018', 'Martins2021')) {
      stop("For type 'CellType', parameter must be one of 'Seidlitz2020', 'Lake2018', or 'Martins2021'.")
    }
    annoData <- get_celltype_data(parameter)
  } else if (type=='SynGO'){
    annoData <- get_SynGO_data()
  }

  getGeneSet=getFromNamespace('getGeneSet','DOSE')
  geneSetList <- getGeneSet(annoData)
  
  if (length(geneSetList) > 10000) {
    warning("The geneSetList is quite large (>10k), this may take some time to process.")
  }
  
  if (convert_to_symbol) {
    geneSetList <- lapply(geneSetList, entrezid2symbol)
  }
  
  return(geneSetList)
}


# # From BioC 3.14 (Nov. 2021, with R-4.2.0)
# library(AnnotationHub)
# library(MeSHDbi)
# ah <- AnnotationHub(localHub=TRUE)
# hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
# file_hsa <- hsa[[1]]
# db <- MeSHDbi::MeSHDb(file_hsa)




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
#' @importFrom readxl read_xlsx
#' @importFrom dplyr select rename mutate %>%
#' @importFrom tidyr unnest
#' @importFrom utils download.file unzip
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
      rename(term = id,gene=hgnc_symbol) %>%
      mutate(gene = strsplit(gene, ";")) %>%
      unnest(cols = c(gene))

    # Process TERM2NAME
    TERM2NAME <- data %>%
      select(id, name) %>%
      rename(term = id, description = name)

    build_Anno=getFromNamespace("build_Anno", "DOSE")
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
#' @importFrom dplyr filter rename select %>%
#' @importFrom utils download.file
#' @importFrom clusterProfiler read.gmt
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
      TERM2GENE <- read.csv(temp_file) %>% 
              mutate(term=class) %>% 
              filter(!gene == '') %>% 
              select(term, gene)
      TERM2NAME <- TERM2GENE %>% 
                    mutate(description=term) %>% 
                    select(term, description)
  
    } else {
      # Read GMT file
      TERM2GENE <- suppressWarnings({read.gmt(temp_file) %>% 
                    filter(!gene == '') %>%
                    select(term, gene)})
      TERM2NAME <- TERM2GENE %>% mutate(description=term) %>% select(term, description)
    }

    build_Anno=getFromNamespace("build_Anno", "DOSE")
    # Create gene set list
    USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)


    message("Cell type data has been processed and temporary files removed.")
    return(USER_DATA)
  }, error = function(e) {
    stop("An error occurred while processing cell type data: ", e$message)
  })
}




#' Filter Gene Set List
#'
#' This function filters a list of gene sets based on the background genes and specified size constraints.
#'
#' @param bg_genes A vector of background gene symbols to be used for filtering.
#' @param geneSetList A list of gene sets to be filtered.
#' @param minGSSize Minimum gene set size for filtering.
#' @param maxGSSize Maximum gene set size for filtering.
#' @return A filtered list of gene sets that meet the size constraints and background genes criteria.
#' @export
filter_geneSetList <- function(bg_genes, geneSetList, minGSSize, maxGSSize) {
  # Create temporary values to name the background genes
  tmp.val <- seq_along(bg_genes)
  names(tmp.val) <- bg_genes
  # Filter the gene set list
  geneSet_filter=getFromNamespace('geneSet_filter','DOSE')
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


































