# inversion detection for HGSVC

Collaborators: Ashley D. Sanders and David Porubsky

The main purpose of this R-package is to locate inversions from a collection of Strand-seq bam files.
1. Generates a directional composite file from single cell Strand-seq libraries 
2. Uses breakpointR functionalities to locate strand state changes in the composite file
3. Predicts and genotypes putative inversions based on type of strand state change localized

To generate composite files:
breakpointer("./bamFiles", dataDirectory ="./CompositeFile", createCompositeFile=T, pairedEndReads=T, scaleWindowSize=T, WC.cutoff=0.9, chromosomes=paste0('chr', c(1:22, "X", "Y")), windowsize=250, minReads=50)

