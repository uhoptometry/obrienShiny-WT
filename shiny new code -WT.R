reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

devtools::install_github("SGDDNB/ShinyCell")
install.packages("devtools")
library(Seurat)
library(ShinyCell)

#getExampleData()                       # Download example dataset (~200 MB)
seu = readRDS("WT_3_new2.rds")
scConf = createConfig(seu)
makeShinyApp(seu, scConf, gene.mapping = TRUE,
             shiny.title = "Zebrafish Retina WT dataset") 

# second dataset


#code for shiny app
citation = list(
  author  = "Liu X., Ouyang J.F., Rossello F.J. et al.",
  title   = "",
  journal = "Nature",
  volume  = "586",
  page    = "101-107",
  year    = "2020", 
  doi     = "10.1038/s41586-020-2734-6",
  link    = "https://www.nature.com/articles/s41586-020-2734-6")
makeShinyCodesMulti(
  shiny.title = "Multi-dataset Tutorial", shiny.footnotes = citation,
  shiny.prefix = c("sc1", "sc2"),
  shiny.headers = c("Integrated", "P23H"), 
  shiny.dir = "shinyAppMulti/") 