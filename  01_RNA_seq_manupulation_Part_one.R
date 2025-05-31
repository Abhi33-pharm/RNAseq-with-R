# Install required R packages 

install.packages(c("tidyverse", "ggpubr", "openxlsx"))

# Install bioconductor packages 

BiocManager::install(c("GEOquery", "TCGAbiolinks","DESeq2", "airway"),force = TRUE)

#Load required packages 
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(GEOquery)
library(readxl)
library(naniar)
#Data Manipulation

# 1. Import data
data<- read.csv("GSE183947_fpkm_long.csv")

#Data Exploration
#Dimension of the data 
dim(data)
nrow(data)
ncol(data)


#Examine first few rows
head(data)
head(data,20)
head(data,n=30)

#Examine the last few rows

tail(data)
tail(data, 20)

#sampling 

sample(data)
sample_n(data,100)
sample_frac(data,0.25)

#Is there any missing value
is.na(data)
sum(is.na(data))
miss_var_summary(data)
gg_miss_var(data)
miss_var_which(data)

#missing value removing from data frame
drop_na(data)

#Select Function (for sub-setting Column); 2 methods
# 1. Select single column by col number

select(data, 1)

# 1.1 select multiple column by col number
select(data,c(1,3,5))

# 1.2 select range of column by using ':' operator

select(data,1:4)

# 1.3 Select single column by col name
select(data, samples)

# 1.4 Select multiple column by col name 
select(data,c(gene, samples,tissue))

#2Filter data 

names(data)

#2.1 Filter data using (==)
filter(data, metastasis=="yes")
names(data)

# 2.2 Filter data using (<)
filter(data, fpkm<10)
#2.3 Filter data using (>)

filter(data, fpkm>10)

# 2.4 Filter data using (>=)
filter(data, fpkm>=10)

#2.5 Filter data using (<=)

filter(data, fpkm<=10)

#2.6 Filter data using ("==", "&")
head(data)

filter(data, metastasis=="yes" & fpkm>10)

#Filter data using ("==",|)

filter(data, metastasis=="yes" |fpkm<10)

# Select and Filter 
head(data)
data_new<- select(data, gene, fpkm, metastasis)
data_new<-filter(data_new, metastasis=="yes" & fpkm<10)

# Select and filter using chaining method (pipe operator= Ctrl+shift+M)

data |> 
  select(gene, fpkm, metastasis) |> 
  filter(metastasis=="yes" & fpkm>10) |> 
  head()

#Multiple filtering criteria

data |> 
  filter(gene %in% c("BRCA1", "BRCA2", "FUCA2", "TP53", "MYCN")) |> 
  head()


# 3. Mutate (Creating column )

data |> mutate(log_fpkm=log(fpkm)) |> 
  head()

#4. Grouping and summarizing 


head(data)

data |> 
  filter(gene=="BRCA1" | gene=="BRCA2") |> 
  group_by(tissue) |> 
  summarise(mean(fpkm))

# 4.1 Which gene

data |> 
  filter(gene=="BRCA1" | gene=="BRCA2") |> 
  group_by(tissue, gene) |> 
  summarise(mean(fpkm))



# 4.2 Giving summary column name

data |> 
  filter(gene=="BRCA1" | gene=="BRCA2") |> 
  group_by(tissue, gene) |> 
  summarise(mean_fpkm= mean(fpkm), median_fpkm=median(fpkm))

#5. Arrange 

data |> 
  filter(gene=="BRCA1" | gene=="BRCA2") |> 
  group_by(tissue, gene) |> 
  summarise(mean_fpkm= mean(fpkm), median_fpkm=median(fpkm)) |> 
  arrange(mean_fpkm)

# Arrange in descending order

data |> 
  filter(gene=="BRCA1" | gene=="BRCA2") |> 
  group_by(tissue, gene) |> 
  summarise(mean_fpkm= mean(fpkm), median_fpkm=median(fpkm)) |> 
  arrange(desc(mean_fpkm))











