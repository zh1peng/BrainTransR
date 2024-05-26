# Fetch annotation data from BrainEnrich package and save it
annoData_synGO <- BrainEnrich::get_annoData(type='SynGO')
usethis::use_data(annoData_synGO, overwrite = TRUE, compress = "xz")

# Process PC1 data, filter for regions starting with 'L_' and set row names
brain_data_PC1 <- read.csv('data-raw/desikan_PC1_data.csv') %>% 
            filter(stringr::str_detect(Region, '^L_')) %>%
            tibble::column_to_rownames('Region')
usethis::use_data(brain_data_PC1, overwrite = TRUE, compress = "xz")

# Process Cushing's data, rename columns, filter, and reorder based on PC1 data
Cushing_data_tmp <- read.csv('data-raw/desikan_cushing_data.csv') %>% 
            rename(Region = 'Dependent.Variable',
                   cohend = 'Cohen.s.d') %>%
            filter(stringr::str_detect(Region, '^L_')) %>%
            select(Region, cohend)
# Reorder rows based on variable order in PC1 data
var_order <- read.csv('data-raw/desikan_PC1_data.csv') %>% 
            filter(stringr::str_detect(Region, '^L_')) %>%
            select(Region)
brain_data_Cushing <- left_join(var_order, Cushing_data_tmp, by = 'Region') %>%
            tibble::column_to_rownames('Region')
usethis::use_data(brain_data_Cushing, overwrite = TRUE, compress = "xz")

# Load gene expression data and save it
gene_data_sample <- BrainEnrich::get_geneExp(atlas = 'desikan', rdonor = 'r0.6', hem = 'L')
usethis::use_data(gene_data_sample, overwrite = TRUE, compress = "xz")

# Load coordinate data and save it
coord_dk_lh <- read.csv('data-raw/desikan_centroid.csv', stringsAsFactors = FALSE) %>%
tibble::column_to_rownames('Row')
usethis::use_data(coord_dk_lh, overwrite = TRUE, compress = "xz")


