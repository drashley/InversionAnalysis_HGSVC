# inversion detection for HGSVC paper

Collaborators: Ashley Sanders, and David Porubsky

The main purpose of this R-package is to locate inversions from a collection of Strand-seq bam files.
1. Generates a directional composite file from single cell data
2. Uses breakpointR functionalities to locate strand state changes
3. locates and genotypes putative inversions

To generate composite files:
breakpointer("./bamFiles", dataDirectory ="./CompositeFile", createCompositeFile=T, pairedEndReads=T, scaleWindowSize=T, WC.cutoff=0.9, chromosomes=paste0('chr', c(1:22, "X", "Y")), windowsize=250, minReads=50)

To locate inversions in composite files:

