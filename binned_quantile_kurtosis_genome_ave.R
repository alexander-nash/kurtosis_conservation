library(GenomicRanges)
library(rtracklayer)
library(data.table)

.libPaths("/mnt/biggles/csc_home/an1413/Rlibs/")

library(npsm)
library(rphast)


##################################################################################################

args<-commandArgs(TRUE)

window_size<-as.numeric(args[3]) 	##in bp

species1<-args[1]			## UCSC stle abbreviation (hg19, canFam3, etc)
species2<-args[2]			
chrom<-args[4]				## either "all" or UCSC style notation (chr1, chr2, chrX, etc)

#no.points=5
#one_d_window_size<-as.numeric(args[5])

outDir<-paste0("/mnt/biggles/csc_projects/an1413/kurtosis/", species1, "/binned_quantile_kurtosis_genome_ave/")

if(!dir.exists(file.path(outDir))){
	dir.create(file.path(outDir))
}

## read in the identical sequence locations - output of get_identical_seq_locations.R

identical_seqs_fn<-paste0("/mnt/biggles/csc_projects/an1413/kurtosis/", species1, "/identical_seqs/")

if(chrom == "all"){
	identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", species2, ".*all*"), full.names=T)
} else {
	identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", species2, ".*chr*"), full.names=T)
}

identical_seqs.dt<-rbindlist(lapply(identical_seqs_files, fread))

###get the range of the middle 50% of the distribution of all runs of conserved seqs across the genome - used later in kurtosis calculation

denominator_quants<-quantile(identical_seqs.dt$first.width, c(0.75, 0.25))
denominator<-denominator_quants[1]-denominator_quants[2]
	
identical_seqs.gr<-GRanges(seqnames=identical_seqs.dt$first.seqnames, 
		   ranges=IRanges(start=identical_seqs.dt$first.start,
				  end=identical_seqs.dt$first.end),
		   strand="*")

if(chrom != "all"){
	identical_seqs.gr<-subset(identical_seqs.gr, seqnames(identical_seqs.gr)==chrom)
}

### read in the filter regions - usually a list of exon and repeat coordinates

filter_regions <- import.bed(dir(paste0("/export/data/CNEs/", species1, "/filters/"), pattern="filter_region.*.bed", full.names=T)) ##Filter regions for reference genome - repeats and exons
filtered_identical_seqs.gr<-setdiff(identical_seqs.gr, filter_regions)

####tile the genome

genome_info<-read.table(paste0("/export/data/goldenpath/", species1, "/assembly.sizes"))

genome_seqlengths<-genome_info$V2
names(genome_seqlengths)<-genome_info$V1

tiled_genome<-tileGenome(genome_seqlengths, tilewidth=window_size, cut.last.tile.in.chrom=T)

if(chrom=="all"){
	tiled_chr<-tiled_genome
} else {
	tiled_chr<-tiled_genome[seqnames(tiled_genome)==chrom]
}

####get distribution in each window and calculate kurtosis

tileHits<-findOverlaps(tiled_chr, filtered_identical_seqs.gr)
tileSplit<-split(subjectHits(tileHits), queryHits(tileHits))

kurt_by_tile<-lapply(as.list(1:length(tiled_chr)) , function(x){
			     obj<-tileSplit[[as.character(x)]]
			     if(is.null(obj)){
				     return(0)
			     } else {
				     runs<-filtered_identical_seqs.gr[obj]
				     dist<-width(runs)
				     score<-(quantile(dist, 0.99)-quantile(dist, 0.01))/denominator
				     if(score=="NaN") score = 0
				     return(score)
			     }
			   })
				

tiled_chr$score<-unlist(kurt_by_tile)

#drop windows that arent the correct width - the ends of chromosomes usually

tiled_chr<-tiled_chr[width(tiled_chr)==window_size]

export.wig(tiled_chr, paste0(outDir, species1, "_", species2, "_", chrom, "_", window_size/1000, "kb_min_window.wig"))

print("kurtosis calculated")

