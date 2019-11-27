# Hena R. Ramay
# This code is modified from the HMP2 bitbucket repository 
#https://bitbucket.org/biobakery/hmp2_analysis/src/default/common/

# It attempts calculates the dysbiosis score from 
# "Multi-omics of the gut microbial ecosystem 
# in inflammatory bowel diseases",(Lloyd-Price et al.)
# https://doi.org/10.1038/s41586-019-1237-9
# and generate fig 2C
# for patients present in the metagenomics dataset.



library(tidyverse)
library(ggplot2)

# Taken from
# https://bitbucket.org/biobakery/hmp2_analysis/src/default/common/pcl_utils.r

pcl.normalize <- function(dat, s = rowSums(dat, na.rm = T)) {
  # Normalize samples to sum to 1
  s[s == 0] <- 1
  #if (dat$nf >= 1)
  dat <- dat / s
  return (dat)
}


# Read in the metadata file
# File taken from 
# https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv

metadata <- read.csv("data/hmp2_metadata.csv", stringsAsFactors = F)

# Read in the taxonoimc profiles from the metagenomics file

min_reads = 1000000

# Reads in the metagenomics taxonomic_profiles file taken from
# (https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products/taxonomic_profiles.tsv.gz )
# Joins metadata file to this data and finds data_type == metagenomics based on the External.IDs
# It also keeps only samples with > 1 M reads.

metag <-
  as.data.frame(t(
    read.csv(
      "data/taxonomic_profiles.tsv",
      sep = "\t",
      row.names = 1,
      stringsAsFactors = F
    )
  ))  %>%
  tibble::rownames_to_column(var = "External.ID") %>%
  left_join(by = "External.ID", x = ., y = metadata) %>%
  filter(., data_type == "metagenomics") %>%  filter(reads_filtered >= min_reads) %>%
  select(-External.ID)


# Reference nonIBD set with samples only from  > 20th week
ref_set <- (metag$diagnosis == "nonIBD") & (metag$week_num >= 20)

# Taken and modified from pcl.only code in
# https://bitbucket.org/biobakery/hmp2_analysis/src/default/common/pcl_utils.r

temp <-
  metag[, grepl(sprintf("(^|\\|)%s__[^\\|]+$", 's'), colnames(metag))] %>% pcl.normalize(.)


# Modified from
# https://bitbucket.org/biobakery/hmp2_analysis/src/default/common/disease_activity.r

# calculate bray-curtis dissimilarity
dist<-as.matrix(vegdist(x = temp, method = "bray")) 

# Calculate the activity index
metag$dysbiosis_score <- sapply(seq_along(ref_set), function(i) {
  print(i)
  median(dist[i, ref_set &
                (metag$Participant.ID != metag$Participant.ID[i])])
})


disease_activity_threshold <-
  quantile(metag$dysbiosis_score[metag$diagnosis == "nonIBD"], 0.9)
eubiosis_lower_threshold <-
  quantile(metag$dysbiosis_score[metag$diagnosis == "nonIBD"], 0.1)

metag$dysbiotic <- metag$dysbiosis_score >= disease_activity_threshold



fig2c <- ggplot(metag,
                aes(
                  x = dysbiosis_score,
                  group = diagnosis,
                  color = diagnosis,
                  fill = diagnosis
                )) +
  geom_density() + scale_fill_manual(values = alpha(c("red", "blue", "yellow"), .3)) +
  scale_color_manual(values = alpha(c("red", "blue", "yellow"), .3)) +
  geom_vline(xintercept = disease_activity_threshold) +
  xlab("Dysbiosis Score")

ggsave(
  "results/Fig2c_hmp2_ibd.png",
  plot = fig2c,
  width = 9,
  height = 5
)


table(metag$active)

write_csv(path =  "results/dysbiosis_Score.csv",
          x = metag[, c("Participant.ID", "dysbiosis_score", "dysbiotic")],
          col_names = T)
