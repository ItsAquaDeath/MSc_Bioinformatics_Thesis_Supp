pacman::p_load(Rsubread)

# Set wd and load names of sorted BAM files (contained in a txt file)
work_dir = paste0(getwd(), "/")
BAM_files = read.table("bamfiles.txt")
BAM_files = paste0(work_dir, "bam_files/" , BAM_files[,1])

# Execute featureCounts to create counts matrix using the genomic.gff
fc_RNAseq <- featureCounts(BAM_files, 
                           annot.ext = paste0(work_dir, "gff_files/genomic.gff"),
                           isGTFAnnotationFile = TRUE,
                           chrAliases = paste0(work_dir, "chr.csv"),
                           GTF.featureType = c("CDS", "exon"),
                           GTF.attrType = "gene",
                           GTF.attrType.extra = "pseudogene",
                           countMultiMappingReads = FALSE,
                           isPairedEnd = TRUE,
                           useMetaFeatures = TRUE,
                           nthreads = 12
)

# Save de SummarizedExperiment object
save(fc_RNAseq, file="C13RNAseq.rda")
