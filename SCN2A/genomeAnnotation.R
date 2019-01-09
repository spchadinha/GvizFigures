#!/usr/bin/env Rscript

# Developer: Spencer Chadinha

# Loads the Gviz package and all dependencies
if (!suppressMessages(require(Gviz))) {
	# suppressMessages(
	# 	if (!requireNamespace("BiocManager", quietly = TRUE))
 #    		install.packages("BiocManager")
	# 	BiocManager::install("Gviz", version = "3.8")
	# )
}
# Load RColorBrewer for color assignments
if (!suppressMessages(require(RColorBrewer))) {
	suppressMessages(install.packages("RColorBrewer"))
}

# Read in gene model data and additional data track information
# Gene model
# Note: data was pulled from ENSEMBL BioMart web interface and has 8
# columns described by model_cnames. This could be done with the 
# biomaRt package, but it was giving me a hard time when I tried.
print("Reading gene model data")
model <- read.csv("./SCN2A/SCN2A_biomart.tsv", sep='\t')
# Gviz requires specific column names for GeneRegionTrack, which
# I use to make the gene model
model_cnames <- c("chromosome", "start", "end", "strand", 
				"gene", "exon", "transcript", "symbol")
colnames(model) <- model_cnames

# Create a color column from color blind friendly color palletes for each transcript
# Work in progress
q <- brewer.pal.info[brewer.pal.info$category == 'qual' & brewer.pal.info$colorblind,]
col_vector = unlist(mapply(brewer.pal, q$maxcolors, rownames(q)))
cols <- col_vector[1:length(unique(model[,'transcript']))]
transcript_ids <- as.character(model[,'transcript'])
col_list <- character(length(transcript_ids))
at = 1
for (i in transcript_ids) {
	col_list[at] <- cols[i]
	at = at + 1
}
model$color=col_list

# Canonical gene model data
# Note: This model was derived from the above model using a python
# script to conglomerate all unique exons into a "super" transcript.
print("Reading union model data")
cmodel <- read.csv("./SCN2A/SCN2A_union.tsv", sep='\t')
colnames(cmodel) <- model_cnames

# H3K4Me1 Data
# Note: this data is pulled from UCSC, but the original file contained
# comment lines that needed to be remove for Gviz to use the data
print("Reading methylation data")
mdata <- read.csv("./SCN2A/SCN2A_H3K4Me1.txt", sep='\t')

# Declare genomic ranges and gene names
genome <- "hg38"
chr <- "chr2"

# # Create a chromosome ideogram track
# itrack <- IdeogramTrack(
# 	genome = genome, 
# 	chromosome = chr)

# Create an axis track to show location on the chromosome
gtrack <- GenomeAxisTrack()

# Define the SNP location and track
SNP_name = "rs17183814"
SNP_start = 165295879
SNP_end = 165295879

snptrack <- AnnotationTrack(
	start = SNP_start, 
	end = SNP_end, 
	chromosome = chr,  
	genome = genome, 
	transcriptAnnotation = SNP_name, 
	name = SNP_name)

# Create gene model track from the model Data.Frame
genetrack <- GeneRegionTrack(model, 
	genome = genome, 
	chromosome = chr,
	transcriptAnnotation = "transcript", 
	name = "SCN2A")

# Create canonical gene model track and highlight junction of interest 
jxn_start = 165239641
jxn_end = 165295879

union <- GeneRegionTrack(cmodel, 
	genome = genome, 
	chromosome = chr,
	transcriptAnnotation = "gene", 
	name = "union")

jxn_union <- HighlightTrack(
	trackList = c(union), 
	start = c(jxn_start), 
	end = c(jxn_end), 
	chromosome = chr)

# Create H3K4Me1 histogram track
methyltrack <- DataTrack(mdata, 
	genome=genome, 
	chromosome=chr, 
	type="histogram", 
	name="H3K4Me1")

# Create a list of tracks and plot them with Gviz
# tracks <- list(itrack, gtrack, snptrack, genetrack, jxn_canonical, methyltrack)
tracks <- list(gtrack, snptrack, genetrack, jxn_union, methyltrack)

pdf("./SCN2A/SCN2A_figure_script.pdf")
	plotTracks(tracks)
dev.off()