#!/usr/bin/env Rscript

# VCF Metrics Line Chart Generator
# This script generates a line chart for VCF quality metrics (DP, FS, QD, MQ)
# across chromosomal positions

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(tools)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input VCF file", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output image file (PNG/JPEG). If not specified, uses the input filename with .png extension")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if input file is provided
if (is.null(opt$input)) {
  stop("Input VCF file must be provided. Use --help for more information.")
}

# Set output filename if not provided
if (is.null(opt$output)) {
  opt$output <- paste0(file_path_sans_ext(opt$input), ".png")
}

# Check file extension
output_ext <- tolower(file_ext(opt$output))
if (!(output_ext %in% c("png", "jpg", "jpeg"))) {
  warning("Output format not recognized. Using PNG format.")
  opt$output <- paste0(file_path_sans_ext(opt$output), ".png")
  output_ext <- "png"
}

cat("Processing file:", opt$input, "\n")
cat("Output will be saved as:", opt$output, "\n")

# Function to parse VCF files and extract relevant information
parse_vcf <- function(vcf_file) {
  # Read VCF file, skipping header lines
  lines <- readLines(vcf_file)
  header_end <- which(startsWith(lines, "#CHROM"))
  
  if (length(header_end) == 0) {
    stop("Invalid VCF format: CHROM header line not found")
  }
  
  # Read data lines
  data_lines <- lines[(header_end+1):length(lines)]
  
  # Create a temporary file with data lines
  temp_file <- tempfile()
  writeLines(data_lines, temp_file)
  
  # Read as data.table
  vcf_data <- fread(temp_file, sep="\t", header=FALSE, 
                    col.names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"),
                    fill=TRUE, na.strings=c(".", ""), stringsAsFactors=FALSE)
  
  # Clean up
  unlink(temp_file)
  
  return(vcf_data)
}

# Function to extract metrics from INFO field
extract_metrics <- function(vcf_data) {
  # Extract DP, FS, QD, MQ from INFO field
  metrics_data <- vcf_data %>%
    select(CHROM, POS, INFO) %>%
    mutate(
      DP = as.numeric(gsub(".*DP=([0-9]+).*", "\\1", INFO)),
      FS = as.numeric(gsub(".*FS=([0-9.]+).*", "\\1", INFO)),
      QD = as.numeric(gsub(".*QD=([0-9.]+).*", "\\1", INFO)),
      MQ = as.numeric(gsub(".*MQ=([0-9.]+).*", "\\1", INFO))
    ) %>%
    select(CHROM, POS, DP, FS, QD, MQ)
  
  return(metrics_data)
}

# Function to create genomic coordinates (for better visualization across chromosomes)
create_genomic_coordinates <- function(metrics_data) {
  # Get chromosome lengths (approximation based on data)
  chrom_lengths <- metrics_data %>%
    group_by(CHROM) %>%
    summarize(max_pos = max(POS, na.rm = TRUE)) %>%
    mutate(offset = lag(cumsum(max_pos), default = 0))
  
  # Join with original data to create genome coordinates
  metrics_data <- metrics_data %>%
    left_join(chrom_lengths, by = "CHROM") %>%
    mutate(genome_pos = POS + offset)
  
  # Create chromosome boundaries for plotting
  chrom_boundaries <- chrom_lengths %>%
    mutate(
      boundary = offset + max_pos/2,
      chrom_label = CHROM
    )
  
  return(list(metrics_data = metrics_data, chrom_boundaries = chrom_boundaries))
}

# Main function to create the visualization
create_vcf_metrics_plot <- function(vcf_file, output_file) {
  # Parse VCF file
  vcf_data <- parse_vcf(vcf_file)
  
  # Extract metrics
  metrics_data <- extract_metrics(vcf_data)
  
  # Create genomic coordinates
  coord_data <- create_genomic_coordinates(metrics_data)
  metrics_data <- coord_data$metrics_data
  chrom_boundaries <- coord_data$chrom_boundaries
  
  # Convert to long format for plotting
  plot_data <- metrics_data %>%
    select(CHROM, POS, genome_pos, DP, FS, QD, MQ) %>%
    pivot_longer(cols = c(DP, FS, QD, MQ), names_to = "Metric", values_to = "Value")
  
  # Handle potential outliers for better visualization
  plot_data <- plot_data %>%
    group_by(Metric) %>%
    mutate(
      q1 = quantile(Value, 0.01, na.rm = TRUE),
      q99 = quantile(Value, 0.99, na.rm = TRUE),
      Value = pmin(pmax(Value, q1), q99)
    ) %>%
    select(-q1, -q99)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = genome_pos, y = Value, color = Metric)) +
    geom_line(alpha = 0.7, size = 0.5) +
    geom_vline(data = chrom_boundaries, aes(xintercept = offset), 
               linetype = "dashed", color = "darkgray", alpha = 0.7) +
    scale_color_manual(values = c("DP" = "blue", "FS" = "red", "QD" = "green", "MQ" = "purple")) +
    labs(
      title = "VCF Quality Metrics Across Chromosomal Positions",
      x = "Chromosomal Position",
      y = "Value",
      color = "Metric"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    scale_x_continuous(
      breaks = chrom_boundaries$boundary,
      labels = chrom_boundaries$chrom_label
    ) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 1)
  
  # Save the plot
  ggsave(output_file, plot = p, width = 12, height = 10, dpi = 300)
  cat("Plot saved to:", output_file, "\n")
}

# Execute the main function
create_vcf_metrics_plot(opt$input, opt$output)