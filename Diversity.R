#Install programs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
install.packages('knitr')

#Open required programs
library("phyloseq")
library(microbiome)
library(knitr)

#OTU have to be raw, not normalized, import the data
Reads_aa <- read.csv(file = 'OTU_1.csv', row.names=1)
Metadata_aa <- read.csv(file = 'Metadata_1.csv', row.names=1)
Taxonomy_aa <- read.csv(file = 'Tax.csv', row.names=1)
Taxonomy_aa1 <- as.matrix(Taxonomy_aa)
Reads_aa1 <- as.matrix(Reads_aa)


#Create tables
OTU <- otu_table(Reads_aa1, taxa_are_rows = TRUE)
TAX <- tax_table(Taxonomy_aa1)
physeq = phyloseq(OTU, TAX)

#Assign Richness to the variable
Richness <- (estimate_richness(physeq, split = TRUE, measures = NULL))

#Write csv file
write.csv(Richness, "Richness.csv")
