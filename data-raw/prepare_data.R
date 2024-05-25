library(dplyr)
# Fetch annotation data from BrainEnrich package and save it
syngo_annoData=BrainEnrich::get_annoData(type='SynGO')
usethis::use_data(syngo_annoData, overwrite = TRUE,compress = "xz")

# Read permutation ID data from raw data directory and save it
perm_id=readRDS('data-raw/desikan_spin_brain_perm_id.rds')
usethis::use_data(perm_id, overwrite = TRUE,compress = "xz")

# Process PC1 data, filter for regions starting with 'L_' and set row names
PC1_data=read.csv('data-raw/desikan_PC1_data.csv') %>% 
            filter(stringr::str_detect(Region, '^L_'))%>%
            tibble::column_to_rownames('Region')
usethis::use_data(PC1_data, overwrite = TRUE,compress = "xz")


# Process Cushing's data, rename columns, filter, and reorder based on PC1 data
Cushing_data.tmp=read.csv('data-raw/desikan_cushing_data.csv') %>% 
            rename(Region = 'Dependent.Variable',
                    cohend='Cohen.s.d') %>%
            filter(stringr::str_detect(Region, '^L_'))%>%
            select(Region, cohend)
# Reorder rows based on variable order in PC1 data
var_order=read.csv('data-raw/desikan_PC1_data.csv') %>% 
            filter(stringr::str_detect(Region, '^L_'))%>%
            select(Region)
Cushing_data=left_join(var_order,Cushing_data.tmp,by='Region')%>%
            tibble::column_to_rownames('Region')
usethis::use_data(Cushing_data, overwrite = TRUE,compress = "xz")

# Load gene expression data and save it
gene_data <- BrainEnrich::get_geneExp(atlas = 'desikan', rdonor = 'r0.6', hem = 'L')
usethis::use_data(gene_data, overwrite = TRUE,compress = "xz")

# load coord data and save it
coord=read.csv('data-raw/desikan_centroid.csv',stringsAsFactors = FALSE) %>% 
  tibble::column_to_rownames('Row')
usethis::use_data(coord, overwrite = TRUE,compress = "xz")

