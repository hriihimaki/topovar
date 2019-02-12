# Set environment for the lsp_wrapper.R
library(RSAGA)
library(parallel)

# Set SAGA GIS location:
SAGA_path <- "C:/Program Files/QGIS 3.4/apps/saga-ltr" 

# Set cores 
n_cores <- detectCores()#-1

# Set my environment
myenv <- rsaga.env(workspace=getwd(), path=SAGA_path, cores = n_cores)