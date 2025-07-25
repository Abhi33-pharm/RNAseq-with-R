# Install packages 
BiocManager::install("biomaRt",force = TRUE)

#Load Packages 
library(biomaRt)
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(readr)
library(stringr)

#Import Raw data 

data<-read_csv("node_score.csv",
               col_types=cols(.default = col_double(),
                              node_name=col_character(),
                              ClusteringCoefficient =col_character()))
data<-data |> 
  mutate(
    ClusteringCoefficient=as.numeric(str_remove(ClusteringCoefficient,",$"))
  )


#Provide gene symbol

data<-data |> 
  mutate(protein_id=sub("^\\d+\\.(ENSP\\d+)$", "\\1", node_name))
  
mart<-useDataset("hsapiens_gene_ensembl",useMart("ensembl"))

symbols<-getBM(attributes =c("ensembl_peptide_id","hgnc_symbol"),
               filters = "ensembl_peptide_id",
              values = unique(data$protein_id),
              mart = mart) |> 
  distinct(ensembl_peptide_id,.keep_all = TRUE)

data_final<-left_join(data,symbols,
                      by=c("protein_id"="ensembl_peptide_id"))


write.csv(data_final, "node_score_with_symbols.csv", row.names = FALSE)

glimpse(data)

