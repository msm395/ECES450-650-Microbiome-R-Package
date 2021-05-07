# import packages along with some example data sets
options(max.print = 1000) # can adjust this value to alter the amount of info that is printed onto the console
library(microbiome)
library(vegan)
library(dplyr)
library(tidyverse)
library(MicrobeDS)

data("atlas1006")
data("dietswap")
data("peerj32")
data("MovingPictures")

# ------------------------------------------
#assign variables to the data sets in phyloseq format
pseq <- atlas1006
pseq2 <- dietswap
pseq3 <- peerj32$phyloseq
pseq4 <- MovingPictures

#-------------------------------------------
#Simple Divergence Function
#compares a subset of samples within the atlas1006 study. primarily adults that are lean are compared to the 
#whole sample set for dissimilarity based on different species abundances (0 = 100% similar, 1 = 0% similar). 

pseqsub <- subset_samples(atlas1006, bmi_group == "lean")
ref <- apply(abundances(pseqsub), 1, median)
divSUB <- divergence(pseq, ref, method = "bray")
print("lean group divergence")
print(divSUB)

#-------------------------------------------
# simple divergence function part 2 using different example data set (dietswap)
#compares a subset of samples within the dietswap study. primarily adults that are males are compared to the whole sample set
# for dissimilarity based on the species compositions (0 = share the exact same species, 1 = no species shared  )

pseq2sub <- subset_samples(dietswap, sex == "male")
ref2 <- apply(abundances(pseq2sub),1, median)
divSUB2 <- divergence(pseq2, ref2, method = "jaccard")
print("male group divergence")
print(divSUB2)

#---------------------------------------------
#distance function could also be used for jaccard and is probably preferable since jaccard does not factor
#abundances like bray. Used the third data set for this one.

distcalc <- distance(pseq3, "jaccard", binary = TRUE)
print("jaccard distance function")
print(distcalc)

#-------------------------------------------
#more complex use cases for divergence. 

betas <- list()
groups <- as.character(unique(meta(pseq3)$group))
for (g in groups) {
  
  df <- subset(meta(pseq3), group == g)
  beta <- c()
  
  for (subj in df$subject) {
    # Pick the samples for this subject
    dfs <- subset(df, subject == subj)
    # Check that the subject has two time points
    if (nrow(dfs) == 2) {
      s <- as.character(dfs$sample)
      # Here with just two samples we can calculate the
      # beta diversity directly
      beta[[subj]] <- divergence(abundances(pseq3)[, s[[1]]],abundances(pseq3)[, s[[2]]],method = "bray")
    }
  }
  betas[[g]] <- beta
}
# boxplot
df <- as.data.frame(unlist(betas))
s<- rownames(df)
si<- as.data.frame(s)
si<- separate(si, s, into = c('names','s'))
df1<- bind_cols(df, si)
rownames(df1)<- df1$s ; df1$s<- NULL

p<- ggplot(df1, aes(x = names, y = `unlist(betas)`))+ geom_boxplot() + ylab('') + xlab('')

plot(p) 

##---------------------------------------------------------------------------
## Divergence experienced to a single subject may chance over time (increase). The blue best fit
# curve demonstrates the change in divergence as the days progress. Some days have high similarity = 49
# others not so much like day ~180 

s <- "F4" # Selected subject
b <- "UBERON:feces" # Selected body site

# Let us pick a subset
pseq4 <- subset_samples(MovingPictures, host_subject_id == s & body_site == b) 

# Rename variables
sample_data(pseq4)$subject <- sample_data(pseq4)$host_subject_id
sample_data(pseq4)$sample <- sample_data(pseq4)$X.SampleID

# Tidy up the time point information (convert from dates to days)
sample_data(pseq4)$time <- as.numeric(as.Date(gsub(" 0:00", "", as.character(sample_data(pseq4)$collection_timestamp)), "%m/%d/%Y") - as.Date("10/21/08", "%m/%d/%Y"))

# Order the entries by time
df <- meta(pseq4) %>% arrange(time)

# Calculate the beta diversity between each time point and
# the baseline (first) time point
beta <- c() # Baseline similarity
s0 <- subset(df, time == 0)$sample
# Let us transform to relative abundance for Bray-Curtis calculations
a <-microbiome::abundances(microbiome::transform(pseq4, "compositional")) 
for (tp in df$time[-1]) {
  # Pick the samples for this subject
  # If the same time point has more than one sample,
  # pick one at random
  st <- sample(subset(df, time == tp)$sample, 1)
  # Beta diversity between the current time point and baseline
  b <- vegdist(rbind(a[, s0], a[, st]), method = "bray")
  # Add to the list
  beta <- rbind(beta, c(tp, b))
}
colnames(beta) <- c("time", "beta")
beta <- as.data.frame(beta)

theme_set(theme_bw(20))
library(ggplot2)
p <- ggplot(beta, aes(x = time, y = beta)) +
  geom_point() +
  geom_line() +
  geom_smooth() +
  labs(x = "Time (Days)", y = "Beta diversity (Bray-Curtis)")
print(p)


