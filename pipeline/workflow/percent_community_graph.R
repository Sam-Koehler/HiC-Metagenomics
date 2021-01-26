# Percent Community Composition Figure 

# arguments 
args = commandArgs(trailingOnly=TRUE)

# gtdbtk files
# 1 : gtdbtk input file1
# 2 : gtdbtk input file2
# 3 : gtdbtk input file3
# 4 : gtdbtk input file4

# checkm files
# 5 : gtdbtk input file5
# 6 : gtdbtk input file6
# 7 : gtdbtk input file7
# 8 : gtdbtk input file8

# 9 : graph classification
# 10 : number of bins
# 11 : annotated only boolean

# libraries
library(ggplot2)


# read in data
### HEALTHY PONDS
## algae52 ###### ###### ###### ###### ###### ######
algae_52_gtbtk_path <- args[1]
algae_52_checkm_path <- args[5]

# data
algae_52_gtbtk_data <- read.table(algae_52_gtbtk_path, sep = "\t")
algae_52_checkm_data <- read.table(algae_52_checkm_path, header = TRUE, sep = "\t")

# merge dataframes
algae52_combined <- merge(algae_52_gtbtk_data, algae_52_checkm_data, by.x = "V1", by.y = "Bin.Id")
colnames(algae52_combined) <- c("user_genome","domain", "Phylum", "class", "order", "family", "genus", "species","bin_size","mapped_reads","pecent_mapped","percent_binned","percent_community")

# sort dataframe
algae52_sorted <- algae52_combined[order(algae52_combined$percent_community, decreasing = TRUE),]

## algae53 ###### ###### ###### ###### ###### ######
algae_53_gtbtk_path <- args[2]
algae_53_checkm_path <- args[6]

# data
algae_53_gtbtk_data <- read.table(algae_53_gtbtk_path, sep = "\t")
algae_53_checkm_data <- read.table(algae_53_checkm_path, header = TRUE, sep = "\t")

# merge dataframes
algae53_combined <- merge(algae_53_gtbtk_data, algae_53_checkm_data, by.x = "V1", by.y = "Bin.Id")
colnames(algae53_combined) <- c("user_genome","domain", "Phylum", "class", "order", "family", "genus", "species","bin_size","mapped_reads","pecent_mapped","percent_binned","percent_community")

# sort dataframe
algae53_sorted <- algae53_combined[order(algae53_combined$percent_community, decreasing = TRUE),]



### SICK PONDS
## algae114 ###### ###### ###### ###### ###### ######
algae_114_gtbtk_path <- args[3]
algae_114_checkm_path <- args[7]

# data
algae_114_gtbtk_data <- read.table(algae_114_gtbtk_path, sep = "\t")
algae_114_checkm_data <- read.table(algae_114_checkm_path, header = TRUE, sep = "\t")

# merge dataframes
algae114_combined <- merge(algae_114_gtbtk_data, algae_114_checkm_data, by.x = "V1", by.y = "Bin.Id")
colnames(algae114_combined) <- c("user_genome","domain", "Phylum", "class", "order", "family", "genus", "species","bin_size","mapped_reads","pecent_mapped","percent_binned","percent_community")

# sort dataframe
algae114_sorted <- algae114_combined[order(algae114_combined$percent_community, decreasing = TRUE),]

## algae115 ###### ###### ###### ###### ###### ######
algae_115_gtbtk_path <- args[4]
algae_115_checkm_path <- args[8]

# data
algae_115_gtbtk_data <- read.table(algae_115_gtbtk_path, sep = "\t")
algae_115_checkm_data <- read.table(algae_115_checkm_path, header = TRUE, sep = "\t")

# merge dataframes
algae115_combined <- merge(algae_115_gtbtk_data, algae_115_checkm_data, by.x = "V1", by.y = "Bin.Id")
colnames(algae115_combined) <- c("user_genome","domain", "Phylum", "class", "order", "family", "genus", "species","bin_size","mapped_reads","pecent_mapped","percent_binned","percent_community")

# sort dataframe
algae115_sorted <- algae115_combined[order(algae115_combined$percent_community, decreasing = TRUE),]

# plotting function
graph_community_comp <- function(num_bins, classification, annotated_only){
  
  all_ponds <- data.frame()
  
  if (annotated_only == FALSE){
    # algae 52 top two species and percentages
    algae52_tops <- algae52_sorted[1:num_bins,c(classification, "percent_community")]
    
    # sum what percent of community was selected
    percent_selected <- 0
    for (i in seq(num_bins)){
      percent_selected <- percent_selected + as.numeric(algae52_tops[i,"percent_community"])
    }
    # calculate and assign remaining percent community
    remaining <- 100 - percent_selected
    algae52_tops[num_bins + 1, c(classification, "percent_community")] <- c("remaining", remaining) 
    algae52_tops$pond <- rep(c("Healthy1"), num_bins + 1)
    
    # algae 53 top two species and percentages
    algae53_tops <- algae53_sorted[1:num_bins,c(classification, "percent_community")]
    
    # sum what percent of community was selected
    percent_selected <- 0
    for (i in seq(num_bins)){
      percent_selected <- percent_selected + as.numeric(algae53_tops[i,"percent_community"])
    }
    # calculate and assign remaining percent community
    remaining <- 100 - percent_selected
    algae53_tops[num_bins + 1, c(classification, "percent_community")] <- c("remaining", remaining) 
    algae53_tops$pond <- rep(c("Healthy2"), num_bins + 1)
    
    # algae 114 top two species and percentages
    algae114_tops <- algae114_sorted[1:num_bins,c(classification, "percent_community")]
    
    # sum what percent of community was selected
    percent_selected <- 0
    for (i in seq(num_bins)){
      percent_selected <- percent_selected + as.numeric(algae114_tops[i,"percent_community"])
    }
    # calculate and assign remaining percent community
    remaining <- 100 - percent_selected
    algae114_tops[num_bins + 1, c(classification, "percent_community")] <- c("remaining", remaining) 
    algae114_tops$pond <- rep(c("Sick1"), num_bins + 1)
    
    # algae 115 top two species and percentages
    algae115_tops <- algae115_sorted[1:num_bins,c(classification, "percent_community")]
    
    # sum what percent of community was selected
    percent_selected <- 0
    for (i in seq(num_bins)){
      percent_selected <- percent_selected + as.numeric(algae114_tops[i,"percent_community"])
    }
    # calculate and assign remaining percent community
    remaining <- 100 - percent_selected
    algae115_tops[num_bins + 1, c(classification, "percent_community")] <- c("remaining", remaining) 
    algae115_tops$pond <- rep(c("Sick2"), num_bins + 1)
    
    # combine sick and healthy ponds
    # add algae52 
    start <- 1
    end <- num_bins + 1
    all_ponds[start:end, 1:3] <- algae52_tops
    
    # add algae53
    start <- end + 1
    end <- start + num_bins 
    all_ponds[start:end, 1:3] <- algae53_tops
    
    # add algae114
    start <- end + 1
    end <- start + num_bins 
    all_ponds[start:end,1:3] <- algae114_tops
    
    # algae 115
    start <- end + 1
    end <- start + num_bins 
    all_ponds[start:end,1:3] <- algae115_tops
    
    # make scale 0-1
    all_ponds$percent_community <- as.numeric(all_ponds$percent_community)/100
    
  }
  
  else{
    # algae 52
    pond_total <- sum(algae52_sorted[1:num_bins,"percent_community"])
    algae52_tops <- algae52_sorted[1:num_bins,c(classification, "percent_community")]
    algae52_tops$percent_community <- algae52_tops$percent_community/pond_total
    algae52_tops$pond <- rep(c("Healthy1"), num_bins)
    
    # algae 53
    pond_total <- sum(algae53_sorted[1:num_bins,"percent_community"])
    algae53_tops <- algae53_sorted[1:num_bins,c(classification, "percent_community")]
    algae53_tops$percent_community <- algae53_tops$percent_community/pond_total
    algae53_tops$pond <- rep(c("Healthy2"), num_bins)
    
    # algae 114
    pond_total <- sum(algae114_sorted[1:num_bins,"percent_community"])
    algae114_tops <- algae114_sorted[1:num_bins,c(classification, "percent_community")]
    algae114_tops$percent_community <- algae114_tops$percent_community/pond_total
    algae114_tops$pond <- rep(c("Sick1"), num_bins)
    
    # algae 115
    pond_total <- sum(algae115_sorted[1:num_bins,"percent_community"])
    algae115_tops <- algae115_sorted[1:num_bins,c(classification, "percent_community")]
    algae115_tops$percent_community <- algae115_tops$percent_community/pond_total
    algae115_tops$pond <- rep(c("Sick2"), num_bins)
    
    # combine sick and healthy ponds
    # add algae52 
    start <- 1
    end <- num_bins 
    all_ponds[start:end, 1:3] <- algae52_tops
    
    # add algae53
    start <- end + 1
    end <- start + num_bins - 1
    all_ponds[start:end, 1:3] <- algae53_tops
    
    # add algae114
    start <- end + 1
    end <- start + num_bins - 1
    all_ponds[start:end,1:3] <- algae114_tops
    
    # algae 115
    start <- end + 1
    end <- start + num_bins - 1
    all_ponds[start:end,1:3] <- algae115_tops
  }
  # Stacked + percent
  ggplot(all_ponds, aes_string(fill=classification, y="percent_community", x="pond")) + 
    geom_bar(position="fill", stat="identity") + 
    labs( y = "Percent Community", x = "Pond Sample") +
    theme(title = element_text(size = 25), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15))  
    #scale_fill_manual(values=c("#ABD08B", "#8BCDD0","#E2C990", "#D08E8B", "#B08BD0", "#90A9E2"))
  
}


# function call and arguments

classy <- args[9]
bins <- args[10]
annotated_only <- args[11]
graph_community_comp(bins, classy, annotated_only)





