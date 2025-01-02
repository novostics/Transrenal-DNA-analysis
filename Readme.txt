#read me
This pipeline is designed for fragmentomic analysis of cell-free DNA (cfDNA). The fragmentomic features analyzed include cfDNA end motifs, fragment sizes, observed/expected (O/E) ratios, and end density. The basic workflow for this analysis is as follows:
Step 1: Use 01.bam2bed_hg19.py to extract genomic coordinates, fragment sizes, and end motifs of each sequenced cfDNA fragment based on the human reference genome. Input: .bam file; Output: .bed file.
Step 2: Use 02.bed2motif.py to calculate end motif frequencies, defined as the frequencies of DNA fragments carrying specific ending motifs (combinations of the four nucleotides, e.g., CCCA) among all cfDNA fragments in a sample. Input: .bed file from Step 1.
Step 3: Use 03.bed2size.pl to calculate size frequencies, defined as the frequencies of DNA fragments with specific lengths (e.g., 100 bp) among all cfDNA fragments in a sample. Input: .bed file from Step 1.
Step 4: Use 04.OEratio.demo.py to calculate the O/E ratio, defined as the ratio of the observed percentage of cfDNA fragments from a genomic region (O) to the expected percentage deduced from simulated reads mapped to that region (E). Input: .bed file from Step 1.
Step 5: Use 05.bed2enddensity.pl to calculate end density, defined as the occurrence of fragment ends at genomic loci relative to the center of regions of interest (ROI), normalized by the median of these values across loci spanning 1 kb upstream or downstream of the ROI center. Input: .bed file from Step 1.

