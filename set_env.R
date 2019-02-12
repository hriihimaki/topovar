# Set environment for the lsp_wrapper.R
library(RSAGA)
library(parallel)

# Set workspace
my_workspace <-  getwd() #or e.g. D:/topovar 

# Set SAGA GIS location:
SAGA_path <- "C:/Program Files/QGIS 3.4/apps/saga-ltr" 

# Set my environment
myenv <- rsaga.env(workspace=my_workspace, path=SAGA_path, parallel = T)
