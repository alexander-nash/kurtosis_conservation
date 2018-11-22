# A Novel Measure of Non-coding Genome Conservation Identifies Genomic Regulatory Blocks Within Primates
### Alexander J. Nash<sup>1, 2</sup> and Boris Lenhard<sup>1, 2, *</sup>
<sup>1</sup>Computational Regulatory Genomics Group, MRC London Institute of Medical Sciences, Du Cane Road, London W12 0NN, UK, <sup>2</sup>Institute of Clinical Sciences, Faculty of Medicine, Imperial College London, Hammersmith Campus, Du Cane Road, London W12 0NN, UK.
<sup>*</sup>To whom correspondence should be addressed. 
## Abstract
__Motivation:__ Clusters of extremely conserved non-coding elements (CNEs) mark genomic regions devoted to cis-regulation of key developmental genes in Metazoa. We have recently show that their span coincides with that of topologically associating domains (TADs), making them useful for estimating conserved TAD boundaries in the absence of Hi-C data. The standard approach - detecting CNEs in genome alignments and then establishing the boundaries of their clusters - requires tuning of several parameters and breaks down when comparing closely related genomes.
__Results:__ We present a novel, kurtosis-based measure of pairwise non-coding conservation that requires no pre-set thresholds for conservation level and length of CNEs. We show that it performs robustly across a large span of evolutionary distances, including across the closely related genomes of primates for which standard approaches fail. The method is straightforward to implement and enables detection and comparison of clusters of CNEs and estimation of underlying TADs across a vastly increased range of Metazoan genomes.




## Scripts and data used for publication

General kurtosis scripts are located in the base directory, while all Data is found in /Data/

## How to use the provided scripts

At this point, paths are hardcoded in the scripts. In time these will be updated, but for now the paths will need to be manually changed by the user. 

### get_identical_seq_locations.R

This is the script used to identify all the runs of perfect sequence identity between the two species of interest. It takes __4 arguments__, __the first__ is genome assembly abbreviation of the reference species of interest (e.g. hg38, mm10, dm6), while __the second__ is the genome assembly of the query species which will be compared to the reference species. __The third__ argument is a chromosome argument. This is used to parallelise the script such that each chromosome is run as a separate job. This isn't necessary for smaller genomes, but significantly improves speed for large genomes. The chromosome argument should be provided in UCSC style (e.g. chr1, chr20). __The fourth__ argument is a logical (TRUE/FALSE) indicating whether the job should run on an individual chromosome or not. If this argument is FALSE, then the chromosome argument is ignored. 

For this script to run, __you must have two .axt alignments__, one for reference species to query species, and one for query species to reference species. These can be generated using LASTZ. __You must also have .2bit files__ for both genome assemblies. At the moment the two species arguments serve as a way to locate the .axt and .2bit files in hard coded directories. For now, these directories will have to be manually changed in the script. 

__The script produces a .tsv file containing the coordinates (in the reference species genome) of all runs of sequence identity between the two species.__ 

### binned_quantile_kurtosis_genome_ave.R

This is the script that calculates the kurtosis-based conservation score in bins across the reference genome. It takes the output of the previous script and calculates the kurtosis of the distribution of lengths of runs of prefect sequence identity in bins across the genome. 

The script takes __4 arguments__. As with the previous script, __the first__ is the genome assembly abbreviation of the reference species of interest (e.g. hg38, mm10, dm6), while __the second__ is the genome assembly of the query species which will be compared to the reference species. __The third__ argument is the bin size in which to calculate kurtosis-based conservation. This argument must be provided in bp. __The fourth__ argument is the chromosome argument. Once again, in large genomes it is much faster to parallelise this script and run it on each chromosome separately. The chromosome argument should be provided in UCSC style (e.g. chr1, chr20). If you are using a small genome and want to run this script on the entire genome in one job, the fourth argument should be "all". 

This script requires a __.bed file containing regions to be filtered out before kurtosis calculation.__ Since this is a method of measuring non-coding conservation, this filter file should contain the coordinates of exons, and all known repeats. A file containing the chromosome sizes of the reference species is also required. 

__The script produces a .wig file containing the binned kurtosis-based conservation values across the reference genome.__

### binned_quantile_kurtosis_genome_ave_grbs.R

This script takes the genome-wide kurtosis-based conservation scores generated in the previous script and locates the boundaries of GRBs. 

It takes __5 arguments__. __The first__ argument is  the genome assembly abbreviation of the reference species of interest (e.g. hg38, mm10, dm6), while __the second__ is the genome assembly of the query species which will be compared to the reference species. __The third__ argument is the bin size (in bp) in which the kurtosis-based conservation was calculated in the previous script. __The fourth__ argument is the ARL0 parameter used in the change point detection - see [the CPM R library for more details](https://rdrr.io/cran/cpm/man/processStream.html). In all previous cases we have used a value of 370 for this argument, as this provides the most sensitive change point detection. __The fifth__ argument is the merging percentile for the final GRB identification step. It controls the stringency of the merging of adjacent ranges after boundary identification. A good starting point for this parameter is usually 0.7. 

__This script produces a .bed file containing the coordinates of all identified GRBs in the reference species genome.__
