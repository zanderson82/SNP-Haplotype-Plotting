#!/usr/bin/env Rscript
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(readr)
#library(optparse)

args <- commandArgs(trailingOnly = TRUE)
demographic_file <- args[1]
output_plot <- args[2]
file_list_path <- args [3]
title <- args [4]


#demographic_file <- snakemake@params[["demographic_data"]]
#output_plot <- snakemake@output[["plots"]]
#file_list_path <- snakemake@input[["filtered_files"]]
demographic_df <- read_delim(demographic_file, delim = "\t")
file_list <- read_delim(file_list_path, delim = "\t", col_names = FALSE)$X1
print(file_list)
# Initialize an empty list to store all unique SNPs
all_snps <- list()
# First loop to get all unique SNPs

print("starting first loop")
for (bed_file in file_list) {
  lines <- readLines(bed_file)
  
  # Check if the file contains any data beyond the header
  if (length(lines) > 1) {
    # Extract the column names from the first line
    column_names <- strsplit(lines[1], "\t")[[1]]
    
    # Read the rest of the lines
    bed <- do.call(rbind, strsplit(lines[-1], "\t"))
    bed <- as.data.frame(bed, stringsAsFactors = FALSE)
    
    # Set the column names
    colnames(bed) <- column_names
    
    # Remove any rows that are actually the column names repeated
    bed <- bed[bed$CHROM != "CHROM", ]
    
    # Inspect the data
    print(head(bed))
    
    # Ensure the GENE and EXON columns are correctly interpreted
    bed$GENE <- as.character(bed$GENE)
    bed$EXON <- as.numeric(bed$EXON)
    
    # Check for problematic rows
    problematic_rows <- bed[is.na(bed$POS) | is.na(bed$QUAL) | is.na(bed$EXON), ]
    print(problematic_rows)
    
    # Clean the data if necessary
    bed <- bed[!is.na(bed$POS) & !is.na(bed$QUAL) & !is.na(bed$EXON), ]
    
    # Convert columns to numeric
    bed$POS <- as.numeric(bed$POS)
    bed$QUAL <- as.numeric(bed$QUAL)
    
    ref <- bed$REF
    alt <- bed$ALT
    
    # Check if the ref and alt are both 1 character
    for (i in seq_len(length(ref))) {
      if (nchar(ref[i]) == 1 && nchar(alt[i]) == 1) {
        snp_tag <- paste0("SNP", bed$POS[i], "_", ref[i], "_", alt[i])
        # Generate SNP tag if it hasn't been made yet
        if (!snp_tag %in% names(all_snps)) {
          all_snps[[snp_tag]] <- list("GENE" = bed$GENE[i], "EXON" = bed$EXON[i])
        }
      }
    }
  } else {
    message(paste("The file", bed_file, "contains no SNPs."))
  }
}

snp_info <- lapply(all_snps, function(x) list())

print(snp_info)
str(snp_info)
print("starting second loop")
# Second loop to get genotype, PS tag, GENE, and EXON for each SNP in each sample
for (bed_file in file_list) {
  lines <- readLines(bed_file)
  print("extracting column names")
  
  # Check if the file contains any data beyond the header
  if (length(lines) > 1) {
    # Extract the column names from the first line
    column_names <- strsplit(lines[1], "\t")[[1]]
    
    # Read the rest of the lines
    bed <- do.call(rbind, strsplit(lines[-1], "\t"))
    bed <- as.data.frame(bed, stringsAsFactors = FALSE)
    print("applying the column names to each bed")
    
    # Set the column names
    colnames(bed) <- column_names
    
    # Remove any rows that are actually the column names repeated
    bed <- bed[bed$CHROM != "CHROM", ]
    
    # Inspect the data
    print(head(bed))
    
    # Ensure the GENE and EXON columns are correctly interpreted
    bed$GENE <- as.character(bed$GENE)
    bed$EXON <- as.numeric(bed$EXON)
    
    # Check for problematic rows
    problematic_rows <- bed[is.na(bed$POS) | is.na(bed$QUAL) | is.na(bed$EXON), ]
    print(problematic_rows)
    
    # Clean the data if necessary
    bed <- bed[!is.na(bed$POS) & !is.na(bed$QUAL) & !is.na(bed$EXON), ]
    
    # Convert columns to numeric
    bed$POS <- as.numeric(bed$POS)
    bed$QUAL <- as.numeric(bed$QUAL)
    
    ref <- bed$REF
    alt <- bed$ALT

    geno_data <- substr(bed$SAMPLE, start = 1, stop = 3)

    for (i in seq_len(nrow(bed))) {
      # Check if "PS" is present in the FORMAT column in the bed dataframe
      if ("PS" %in% strsplit(bed$FORMAT, ":")[[1]]) {
        # Retrieve the PS tag from the SAMPLE column as anything after the last colon
        ps_tag <- sapply(strsplit(as.character(bed$SAMPLE), ":"), function(x) tail(x, n = 1))
        # If the PS tag is not a number, set it to NA
        ps_tag[!grepl("^\\d+$", ps_tag)] <- NA
      } else {
        ps_tag <- NA
      }
    }
    
    # Extract the sample name from the BED file name
    sample_name <- basename(bed_file)
    # Extract the sample name by removing anything after the first hyphen
    sample_name <- sub("-.*$", "", sample_name)
    
    # Initialize a vector to keep track of which SNPs are present in this sample
    snps_present <- rep(FALSE, length(all_snps))
    
    # Generating snp_tags
    for (i in seq_len(length(ref))) {
      if (nchar(ref[i]) == 1 && nchar(alt[i]) == 1) {
        snp_tag <- paste0("SNP", bed$POS[i], "_", ref[i], "_", alt[i])
        if (snp_tag %in% names(all_snps)) {
          # Extract the genotype, PS tag, GENE, and EXON from the SAMPLE column
          genotype <- geno_data[i]
          if (genotype %in% c("./.", "0/0")) {
            ps_tag <- NA
          }
          if (is.null(snp_info[[snp_tag]][[sample_name]]$PS)) {
            snp_info[[snp_tag]][[sample_name]]$PS <- ps_tag
          }
          gene <- bed$GENE[bed$POS == bed$POS[i]]
          exon <- bed$EXON[bed$POS == bed$POS[i]]
          snp_info[[snp_tag]][[sample_name]] <- list("genotype" = genotype, "PS" = ps_tag, "ALT" = alt[i], "GENE" = gene, "EXON" = exon)
          snps_present[names(all_snps) == snp_tag] <- TRUE
        }
      }
    }
    
    # For any SNPs not present in this sample, add the default genotype
    for (snp_tag in names(all_snps)[!snps_present]) {
      gene <- all_snps[[snp_tag]]$GENE
      exon <- all_snps[[snp_tag]]$EXON
      snp_info[[snp_tag]][[sample_name]] <- list("genotype" = "0/0", "PS" = NA, "ALT" = ".", "GENE" = gene, "EXON" = exon)
    }
  } else {
    # Handle the case where the file contains no SNPs
    sample_name <- basename(bed_file)
    sample_name <- sub("-.*$", "", sample_name)
    for (snp_tag in names(all_snps)) {
      gene <- all_snps[[snp_tag]]$GENE
      exon <- all_snps[[snp_tag]]$EXON
      snp_info[[snp_tag]][[sample_name]] <- list("genotype" = "0/0", "PS" = NA, "ALT" = ".", "GENE" = gene, "EXON" = exon)
    }
  }
}

# Convert the snp_info list to a data frame
snp_df <- do.call(rbind, lapply(names(snp_info), function(snp_tag) {
  sample_data <- lapply(names(snp_info[[snp_tag]]), function(sample_name) {
    genotype <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$genotype), "0/0", snp_info[[snp_tag]][[sample_name]]$genotype)
    ps <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$PS), NA, snp_info[[snp_tag]][[sample_name]]$PS)
    alt <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$ALT), ".", snp_info[[snp_tag]][[sample_name]]$ALT)
    gene <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$GENE), NA, snp_info[[snp_tag]][[sample_name]]$GENE)
    exon <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$EXON), NA, snp_info[[snp_tag]][[sample_name]]$EXON)
    alleles <- strsplit(genotype, "/")[[1]]
    data.frame(
      SNP = rep(snp_tag, 2),
      Sample = rep(sample_name, 2),
      Haplotype = 1:2,
      Genotype = alleles,
      PS = rep(ps, 2),
      ALT = rep(alt, 2),
      GENE = rep(gene, 2),
      EXON = rep(exon, 2),
      stringsAsFactors = FALSE
    )
  })
  if (length(sample_data) > 0) {
    do.call(rbind, sample_data)
  }
}))

head(snp_df)
final_df <- left_join(snp_df, demographic_df, by = c('Sample' = 'Sample name'))
final_df$Sex <- as.factor(final_df$Sex)
# function to print the alt allele under under a column named "Allele". Also specifies which haplotype, which corresponds to the genotype. 
# Define a function to map genotypes to alleles
get_allele <- function(genotype, alt, ps, haplotype, Sex = NULL) {
  if (!is.null(Sex) && Sex == "male") {
    if (genotype == "1|0" && haplotype == 1) {
      return(alt)
    } else if (genotype == "0|1" && haplotype == 1 | haplotype == 2) {
      return(alt)
    } else if (genotype == "1" || genotype == "1|1") {
      return(alt)
    }
  } else {
    if (genotype == "1|0" && haplotype == 1) {
      return(alt)
    } else if (genotype == "0|1" && haplotype == 2) {
      return(alt)
    } else if (genotype == "1" || genotype == "1|1") {
      return(alt)
    }
  }
  return(NA)
}

final_df$Allele <- mapply(get_allele, final_df$Genotype, final_df$ALT, final_df$PS, final_df$Haplotype, final_df$Sex)

#function to return return a color that corresponds to each alternate allele and white for na which means the sample is the same as the reference at that SNP position. 
get_color <- function(allele) {
  if (is.na(allele)) {
    return("white")
  }  else if (allele == "A") {
    return("red")
  } else if (allele == "T") {
    return("lightblue")
  } else if (allele == "C") {
    return("green")
  } else if (allele == "G") {
    return("yellow")
  } else {
    return("white")  # return white if allele is not A, T, C, or G
  }
}

# Run this code chunk if you want to plot samples from 1000g with demographics
final_df$Color <- sapply(final_df$Allele, get_color)
final_df <- final_df %>% arrange(Sample, desc(Haplotype))
final_df$Sample <- factor(final_df$Sample, levels = unique(final_df$Sample))
final_df$Haplotype <- factor(final_df$Haplotype, levels = c(1, 2))
# make a gene_exon name for each pair. this allows you to facet by gene and exon in one step
final_df$GENE_EXON <- paste(final_df$GENE, final_df$EXON, sep = "_")

# Run this chunk when plotting 1000g snps. This will allow you to include information such as superpopulation and Sex.
# Join demographic information with each sample (Sex, superpop, etc)
##final_df <- left_join(snp_df, demographic_df, by = c('Sample' = 'Sample name'))
##final_df$Sex <- as.factor(final_df$Sex)
# Add XY for samples labeled as 'male; and XX for samples labeled 'female'
final_df <- final_df %>% mutate(Sample = ifelse(Sex == "male", paste0(Sample, " (XY)"), paste0(Sample, " (XX)")))
superpop_list <- split(final_df, final_df$'Superpop. code')

manual_color_mapping <- c("A" = "#8ff035", "T" = "#f70a0a", "C" = "#75a2f0", "G" = "#ffc445", "NA" = "white")
# You can change which gene_exon you want to visualize or just gene if intersect the vcf is filtered with a bed file that does not contain exon information. 
include_gene_exon <- c("OPN1LW_3", "OPN1LW_5", "OPN1MW_3", "OPN1MW_5")

#if (opt$xchrom) {
#run this if you are plotting SNPs on X chromosome and if joining snps and samples with demographic data (output is called final_df) 
df <- final_df %>% mutate(Haplotype = ifelse(Sex == "male" & Haplotype == 2, NA, Haplotype))
#Run this if just making one plot
filtered_df <- df %>% dplyr::filter(GENE_EXON %in% include_gene_exon)
p <- ggplot(filtered_df, aes(x = SNP, y = Haplotype, fill = Allele)) +
    geom_tile(color = "black") +
    geom_text(aes(label = Allele), size = 2, na.rm = TRUE) +
    scale_fill_manual(values = manual_color_mapping, na.value = "white") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 8)) +
    theme(strip.text.y = element_text(angle = 0, size = 6)) +
    theme(legend.text = element_text(size = 8)) +
    theme(legend.position = "right") +
    ylab("Haplotype (1/2)") +
    facet_grid(Sample ~ GENE_EXON, scales = "free", space = "free") +
    scale_y_discrete(breaks = c(1, 2)) +
    ggtitle(title)
ggsave(file = output_plot, plot = p, width = 16, height = 12)


#run this to plot make multiple plots
#create_plot <- function(filtered_df){
  #filtered_df <- df %>% dplyr::filter(GENE_EXON %in% include_gene_exon)
  #m <- ggplot(filtered_df, aes(x = SNP, y = Haplotype, fill = Allele)) +
      #geom_tile(color = "black") +
      #geom_text(aes(label = Allele), size = 2, na.rm = TRUE) +
      #scale_fill_manual(values = manual_color_mapping, na.value = "white") +
      #theme_classic() +
     # theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 8)) +
      #theme(strip.text.y = element_text(angle = 0, size = 7)) +
      #theme(legend.text = element_text(size = 8)) +
      #theme(legend.position = "right") +
      #ylab("Haplotype (1/2)") +
      #facet_grid(Sample ~ GENE_EXON, scales = "free", space = "free") +
      #ggtitle(unique(df$'Superpop. code')) +
      #scale_y_discrete(breaks = c(1, 2)) 


  #m <- ggplotly (m)
  #return(m)
#}

#create a plot for each 'Superpop. code'
#plots <- lapply(superpop_list, create_plot)
#for (i in seq_along(plots)) {
  #superpop_code <- unique(superpop_list[[i]]$'Superpop. code')
  #svg
#}