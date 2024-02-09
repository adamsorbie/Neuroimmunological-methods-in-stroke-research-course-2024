### LIBRARIES 
if(!require("pacman")){
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
packages <- c("BiocManager",
               "tidyverse",
               "EnhancedVolcano",
               "DESeq2",
               "clusterProfiler",
               "AnnotationDbi",
               "org.Mm.eg.db",
               "ggpubr",
               "ggsci",
               "rstatix",
               "phyloseq",
               "vegan",
               "ANCOMBC",
               "zoo")

pacman::p_load(char = packages)

uninstalled <- packages[!(packages %in% installed.packages()[,"Package"])]

if (length(uninstalled) > 0){
  print(paste0("installation not completed successfully, please install: ", uninstalled))
} else {
  print("Installed packages successfully :)")
}
