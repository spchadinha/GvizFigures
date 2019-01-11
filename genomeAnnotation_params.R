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
# Read depths data
print("Reading sequencing data")
rddata_dlpfc <- read.csv("./SCN2A/SCN2A_DLPFC_U.txt", sep='\t')
rddata_amyg <- read.csv("./SCN2A/SCN2A_amyg_U.txt", sep='\t')
rddata_sacc <- read.csv("./SCN2A/SCN2A_sacc_U.txt", sep='\t')
# Gviz requires specific column names for GeneRegionTrack, which
# I use to make the gene model
model_cnames <- c("chromosome", "start", "end", "strand", 
				"gene", "exon", "transcript", "symbol", "transcript_type")
colnames(model) <- model_cnames

#######################################################################################
# Create a color column from color blind friendly color palletes for each transcript
# Work in progress

# column to seperate transcripts by color
seperator <- 'transcript_type'
# Generate a data frame of color panels
q <- brewer.pal.info[brewer.pal.info$category == 'qual' & brewer.pal.info$colorblind,]
# Create a list of color blind friendly colors
col_vector = unlist(mapply(brewer.pal, q$maxcolors, rownames(q)))
# Subset the list to get a number of colors equal to the number of unique transcripts
transcript_cols <- col_vector[1:length(unique(model[,seperator]))]

########################################################################################

# Canonical gene model data
# Note: This model was derived from the above model using a python
# script to conglomerate all unique exons into a "super" transcript.
print("Reading union model data")
cmodel <- read.csv("./SCN2A/SCN2A_union.tsv", sep='\t')
colnames(cmodel) <- model_cnames[1:8]

# Declare genomic ranges and gene names
genome <- "hg38"
chr <- "chr2"

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
	group = model$transcript,
	feature = model$transcript,
	genome = genome, 
	chromosome = chr,
	transcriptAnnotation = "transcript", 
	name = "SCN2A",
	ENST00000636985 = transcript_cols[2],
	ENST00000637266 = transcript_cols[2],
	ENST00000375437 = transcript_cols[2],
	ENST00000486878 = transcript_cols[2],
	ENST00000631182 = transcript_cols[2],
	ENST00000424833 = transcript_cols[2],
	ENST00000640791 = transcript_cols[1],
	ENST00000283256 = transcript_cols[2],
	ENST00000480032 = transcript_cols[3],
	ENST00000636662 = transcript_cols[4],
	ENST00000636135 = transcript_cols[4],
	ENST00000636769 = transcript_cols[4],
	ENST00000638151 = transcript_cols[1],
	ENST00000375427 = transcript_cols[2],
	ENST00000635945 = transcript_cols[1],
	ENST00000636384 = transcript_cols[4],
	ENST00000636071 = transcript_cols[2],
	ENST00000637367 = transcript_cols[4],
	col = NULL)

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

# Create read depths histogram track
# Each brain region has a separate track - stacking info makes a messy graph
rdtrack_d <- DataTrack(rddata_dlpfc, 
                       genome=genome, 
                       chromosome=chr, 
                       type="histogram", 
                       name="DLPFC")

rdtrack_a <- DataTrack(rddata_amyg, 
                       genome=genome, 
                       chromosome=chr, 
                       type="histogram", 
                       name="Amygdala")

rdtrack_s <- DataTrack(rddata_sacc, 
                       genome=genome, 
                       chromosome=chr, 
                       type="histogram", 
                       name="sACC")

# Create a list of tracks and plot them with Gviz
# tracks <- list(itrack, gtrack, snptrack, genetrack, jxn_canonical, methyltrack)
tracks <- list(gtrack, snptrack, genetrack, jxn_union, rdtrack_d)

pdf("./SCN2A/SCN2A_figure_script_tst.pdf")
	plotTracks(tracks)
dev.off()