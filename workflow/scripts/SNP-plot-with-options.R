#!/usr/bin/env Rscript
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(readr)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-d", "--demographic"), type="character", default=NULL,
              help="Path to demographic file [default= %default]"),
  make_option(c("-o", "--output"), type="character", default="plot.pdf",
              help="Output plot file [default= %default]"),
  make_option(c("-f", "--file_list"), type="character",
              help="Path to file list"),
  make_option(c("-t", "--title"), type="character", default="SNP Plot",
              help="Plot title [default= %default]"),
  make_option(c("-x", "--xchrom"), action="store_true", default=FALSE,
              help="Process as X chromosome data [default= %default]"),
  make_option(c("-e", "--exon"), action="store_true", default=FALSE,
              help="Include exon information in visualization [default= %default]"),
  make_option(c("-w", "--width"), type="numeric", default=12,
              help="Plot width [default= %default]"),
  make_option(c("-h", "--height"), type="numeric", default=9,
              help="Plot height [default= %default]"),
  make_option(c("-g", "--genes"), type="character", default=NULL,
              help="Comma-separated list of genes to include"),
  make_option(c("-c", "--gene_exons"), type="character", default=NULL,
              help="Comma-separated list of gene_exons to include (e.g., 'OPN1LW_3,OPN1LW_5')")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read input files
file_list <- read_delim(opt$file_list, delim = "\t", col_names = FALSE)$X1
print(file_list)

# Read demographic data if provided
if (!is.null(opt$demographic)) {
  demographic_df <- read_delim(opt$demographic, delim = "\t")
}

# Initialize an empty list to store all unique SNPs
all_snps <- list()

print("starting first loop")
for (bed_file in file_list) {
  lines <- readLines(bed_file)
  
  if (length(lines) > 1) {
    column_names <- strsplit(lines[1], "\t")[[1]]
    bed <- do.call(rbind, strsplit(lines[-1], "\t"))
    bed <- as.data.frame(bed, stringsAsFactors = FALSE)
    colnames(bed) <- column_names
    bed <- bed[bed$CHROM != "CHROM", ]
    
    print(head(bed))
    
    bed$GENE <- as.character(bed$GENE)
    if (opt$exon) {
      bed$EXON <- as.numeric(bed$EXON)
      bed <- bed[!is.na(bed$POS) & !is.na(bed$QUAL) & !is.na(bed$EXON), ]
    } else {
      bed <- bed[!is.na(bed$POS) & !is.na(bed$QUAL), ]
    }
    
    bed$POS <- as.numeric(bed$POS)
    bed$QUAL <- as.numeric(bed$QUAL)
    
    ref <- bed$REF
    alt <- bed$ALT
    
    for (i in seq_len(length(ref))) {
      if (nchar(ref[i]) == 1 && nchar(alt[i]) == 1) {
        snp_tag <- paste0("SNP", bed$POS[i], "_", ref[i], "_", alt[i])
        if (!snp_tag %in% names(all_snps)) {
          if (opt$exon) {
            all_snps[[snp_tag]] <- list("GENE" = bed$GENE[i], "EXON" = bed$EXON[i])
          } else {
            all_snps[[snp_tag]] <- list("GENE" = bed$GENE[i])
          }
        }
      }
    }
  } else {
    message(paste("The file", bed_file, "contains no SNPs."))
  }
}

snp_info <- lapply(all_snps, function(x) list())

print("starting second loop")
for (bed_file in file_list) {
  lines <- readLines(bed_file)
  
  if (length(lines) > 1) {
    column_names <- strsplit(lines[1], "\t")[[1]]
    bed <- do.call(rbind, strsplit(lines[-1], "\t"))
    bed <- as.data.frame(bed, stringsAsFactors = FALSE)
    colnames(bed) <- column_names
    bed <- bed[bed$CHROM != "CHROM", ]
    
    bed$GENE <- as.character(bed$GENE)
    if (opt$exon) {
      bed$EXON <- as.numeric(bed$EXON)
      bed <- bed[!is.na(bed$POS) & !is.na(bed$QUAL) & !is.na(bed$EXON), ]
    } else {
      bed <- bed[!is.na(bed$POS) & !is.na(bed$QUAL), ]
    }
    
    bed$POS <- as.numeric(bed$POS)
    bed$QUAL <- as.numeric(bed$QUAL)
    
    ref <- bed$REF
    alt <- bed$ALT
    geno_data <- substr(bed$SAMPLE, start = 1, stop = 3)

    for (i in seq_len(nrow(bed))) {
      if ("PS" %in% strsplit(bed$FORMAT, ":")[[1]]) {
        ps_tag <- sapply(strsplit(as.character(bed$SAMPLE), ":"), function(x) tail(x, n = 1))
        ps_tag[!grepl("^\\d+$", ps_tag)] <- NA
      } else {
        ps_tag <- NA
      }
    }
    
    sample_name <- basename(bed_file)
    sample_name <- sub("-.*$", "", sample_name)
    
    snps_present <- rep(FALSE, length(all_snps))
    
    for (i in seq_len(length(ref))) {
      if (nchar(ref[i]) == 1 && nchar(alt[i]) == 1) {
        snp_tag <- paste0("SNP", bed$POS[i], "_", ref[i], "_", alt[i])
        if (snp_tag %in% names(all_snps)) {
          genotype <- geno_data[i]
          if (genotype %in% c("./.", "0/0")) {
            ps_tag <- NA
          }
          if (is.null(snp_info[[snp_tag]][[sample_name]]$PS)) {
            snp_info[[snp_tag]][[sample_name]]$PS <- ps_tag
          }
          gene <- bed$GENE[bed$POS == bed$POS[i]]
          if (opt$exon) {
            exon <- bed$EXON[bed$POS == bed$POS[i]]
            snp_info[[snp_tag]][[sample_name]] <- list("genotype" = genotype, "PS" = ps_tag, "ALT" = alt[i], "GENE" = gene, "EXON" = exon)
          } else {
            snp_info[[snp_tag]][[sample_name]] <- list("genotype" = genotype, "PS" = ps_tag, "ALT" = alt[i], "GENE" = gene)
          }
          snps_present[names(all_snps) == snp_tag] <- TRUE
        }
      }
    }
    
    for (snp_tag in names(all_snps)[!snps_present]) {
      gene <- all_snps[[snp_tag]]$GENE
      if (opt$exon) {
        exon <- all_snps[[snp_tag]]$EXON
        snp_info[[snp_tag]][[sample_name]] <- list("genotype" = "0/0", "PS" = NA, "ALT" = ".", "GENE" = gene, "EXON" = exon)
      } else {
        snp_info[[snp_tag]][[sample_name]] <- list("genotype" = "0/0", "PS" = NA, "ALT" = ".", "GENE" = gene)
      }
    }
  }
}

# Create the data frame with appropriate columns based on options
snp_df <- do.call(rbind, lapply(names(snp_info), function(snp_tag) {
  sample_data <- lapply(names(snp_info[[snp_tag]]), function(sample_name) {
    genotype <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$genotype), "0/0", snp_info[[snp_tag]][[sample_name]]$genotype)
    ps <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$PS), NA, snp_info[[snp_tag]][[sample_name]]$PS)
    alt <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$ALT), ".", snp_info[[snp_tag]][[sample_name]]$ALT)
    gene <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$GENE), NA, snp_info[[snp_tag]][[sample_name]]$GENE)
    alleles <- strsplit(genotype, "/")[[1]]
    
    df <- data.frame(
      SNP = rep(snp_tag, 2),
      Sample = rep(sample_name, 2),
      Haplotype = 1:2,
      Genotype = alleles,
      PS = rep(ps, 2),
      ALT = rep(alt, 2),
      GENE = rep(gene, 2),
      stringsAsFactors = FALSE
    )
    
    if (opt$exon) {
      exon <- ifelse(is.null(snp_info[[snp_tag]][[sample_name]]$EXON), NA, snp_info[[snp_tag]][[sample_name]]$EXON)
      df$EXON <- rep(exon, 2)
      df$GENE_EXON <- paste(gene, exon, sep = "_")
    }
    
    return(df)
  })
  if (length(sample_data) > 0) {
    do.call(rbind, sample_data)
  }
}))

# Process with demographic data if available
if (!is.null(opt$demographic)) {
  final_df <- left_join(snp_df, demographic_df, by = c('Sample' = 'Sample name'))
  final_df$Sex <- as.factor(final_df$Sex)
  if (opt$xchrom) {
    final_df <- final_df %>% 
      mutate(Sample = ifelse(Sex == "male", paste0(Sample, " (XY)"), paste0(Sample, " (XX)"))) %>%
      mutate(Haplotype = ifelse(Sex == "male" & Haplotype == 2, NA, Haplotype))
  }
} else {
  final_df <- snp_df
}

# Define allele mapping function
get_allele <- function(genotype, alt, ps, haplotype, sex = NULL) {
  if (!is.null(sex) && sex == "male" && opt$xchrom) {
    if (genotype == "1|0" && haplotype == 1) {
      return(alt)
    } else if (genotype == "0|1" && (haplotype == 1 | haplotype == 2)) {
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

# Apply allele mapping
if (!is.null(opt$demographic) && opt$xchrom) {
  final_df$Allele <- mapply(get_allele, final_df$Genotype, final_df$ALT, 
                           final_df$PS, final_df$Haplotype, final_df$Sex)
} else {
  final_df$Allele <- mapply(get_allele, final_df$Genotype, final_df$ALT, 
                           final_df$PS, final_df$Haplotype)
}

# Color mapping
manual_color_mapping <- c("A" = "#8ff035", "T" = "#f70a0a", "C" = "#75a2f0", 
                         "G" = "#ffc445", "NA" = "white")

# Filter data based on genes or gene_exons if specified
if (!is.null(opt$genes)) {
  genes_to_include <- unlist(strsplit(opt$genes, ","))
  final_df <- final_df %>% filter(GENE %in% genes_to_include)
} else if (!is.null(opt$gene_exons) && opt$exon) {
  gene_exons_to_include <- unlist(strsplit(opt$gene_exons, ","))
  final_df <- final_df %>% filter(GENE_EXON %in% gene_exons_to_include)
}

# Prepare final data for plotting
final_df <- final_df %>% 
  arrange(Sample, desc(Haplotype)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)),
         Haplotype = factor(Haplotype, levels = c(1, 2)))

# Create the plot
p <- ggplot(final_df, aes(x = SNP, y = Haplotype, fill = Allele)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Allele), size = 2, na.rm = TRUE) +
  scale_fill_manual(values = manual_color_mapping, na.value = "white") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 8),
        strip.text.y = element_text(angle = 0, size = 6),
        legend.text = element_text(size = 8),
        legend.position = "right") +
  ylab("Haplotype (1/2)") +
  {if (opt$exon) 
    facet_grid(Sample ~ GENE_EXON, scales = "free", space = "free")
   else 
    facet_grid(Sample ~ GENE, scales = "free", space = "free")
  } +
  scale_y_discrete(breaks = c(1, 2)) +
  ggtitle(opt$title)

ggsave(file = opt$output, plot = p, width = opt$width, height = opt$height)