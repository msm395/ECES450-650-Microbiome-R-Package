# import packages along with some example data sets

options(max.print = 1000) # can adjust this value to alter the amount of info that is printed onto the console
library(microbiome)
library(knitr)
library(vegan)
data("atlas1006") 
data("dietswap")
data("peerj32")

# ------------------------------------------
#assign variables to the data sets in phyloseq format
pseq <- atlas1006
pseq2 <- dietswap
pseq3 <- peerj32$phyloseq

#------------------------------------------
# Useful Functions to use to better understand datasets (Optional)

OTUData <- otu_table(pseq2) #provides taxa and its abundances
Sampleinfo <- sample_data(pseq2) #provides metrics of each sample
Taxoinfo <- tax_table(pseq2) #provides list of all taxa

print("OTU Table")
print(OTUData)
print("Sample Metrics")
print(Sampleinfo)
print("Taxa Info")
print(Taxoinfo)

# ----------------------------------------------
#Specific functions available to demonstrate certain variations in alpha diversity

alpdiv <- alpha(pseq3, index = "shannon", zeroes = TRUE) #calculates shannon diversity
dom <- dominant(pseq3,level = NULL) # gives a list of the most abundant taxa starting from 1
evexample <- evenness(pseq2,index = "all", zeroes = TRUE, detection = 0) #calculates all methods of measuring evenness
richexample <- estimate_richness(pseq, split = TRUE, measures = NULL) #gives estimator functions such as chao1 and ace
dom2 <- dominance(pseq,index = "all") #calculates dominance of the inputted data set

#-------------------------------------------
#print variable calculated above
print("Shannon Diversity")
print(alpdiv)
print("Dominance Calculation")
print(dom)
print("Evenness Calculation")
print(evexample)
print("Richness Calculation")
print(richexample)
print("Second Dominance Calculation")
print(dom2)

#--------------------------------------------
#perform alpha diversity and assign diversity factors into table

tab1 <- microbiome::alpha(pseq2, index = "all") # calculates all diversity metrics and puts into table
kable(head(tab1))