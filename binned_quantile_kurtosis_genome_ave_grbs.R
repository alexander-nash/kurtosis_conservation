library(GenomicRanges)
library(rtracklayer)
library(data.table)

.libPaths("/mnt/biggles/csc_home/an1413/Rlibs/")

library(cpm)
library(randtests)
library(npsm)
library(rphast)


##########Define variables

args<-commandArgs(TRUE)

species1<-args[1]			## UCSC stle abbreviation (hg19, canFam3, etc)
species2<-args[2]
window_size<-as.numeric(args[3])	## in bp
ARL0=as.numeric(args[4])		## one of  370, 500, 600, 700, ..., 1000, 2000, 3000, ..., 10000, 20000, ..., 50000. see ?processStream
percentile=as.numeric(args[5])		## between 0-1, usually start with 0.7


#####################GRBS##########################

inDir<-paste0("/mnt/biggles/csc_projects/an1413/kurtosis/", species1, "/binned_quantile_kurtosis_genome_ave/")

score_files<-dir(inDir, pattern=paste0(species2, ".*", window_size/1000, "kb_min_window"), full.names=T)

tiled_chr<-lapply(score_files, import.wig, genome=species1)

tiled_chr_all<-Reduce(c, tiled_chr)

edges<-lapply(tiled_chr, function(x){
		      fwrd.edges<-processStream(x$score, cpmType="Mann-Whitney", ARL0=ARL0)
		      rev.edges<-processStream(rev(x$score), cpmType="Mann-Whitney", ARL0=ARL0)
		      return(list("fwrd.edges"=fwrd.edges,
				  "rev.edges"=rev.edges))
})



boundaries<-mapply(x = tiled_chr,
		   y = edges, 
		   FUN = function(x,y){
			   x$fwrd.boundary=0
			   x$rev.boundary=0
			   x$fwrd.boundary[y$fwrd.edges$changePoints]=1
			   revx<-rev(x)
			   revx$rev.boundary[y$rev.edges$changePoints]=1
			   x$rev.boundary=rev(revx$rev.boundary)
			   x$score<-NULL
			   fwrd<-x
			   mcols(fwrd)<-NULL
			   fwrd$score<-x$fwrd.boundary
			   fwrd=subset(fwrd, score==1)
			   rev<-x
			   mcols(rev)<-NULL
			   rev$score<-x$rev.boundary
			   rev=subset(rev, score==1)
			   out<-list("fwrd"=fwrd, "rev"=rev)
			   return(out)
		   }, SIMPLIFY=FALSE)

boundaries<-list("fwrd"=lapply(boundaries, function(x) x$fwrd),
		 "rev"=lapply(boundaries, function(x) x$rev))

if(length(boundaries[[1]]) == 1){
	boundaries<-lapply(boundaries, function(x) split(x[[1]], seqnames(x[[1]])))
}

ranges<-lapply(boundaries, function(x){
		       out<-Reduce(c, lapply(x, function(y){
						     if(length(y) < 2) return(GRanges())
						     ranges<-GRanges()
						     for(i in 1:length(y)){
							     if(i == length(y)) break
							     range<-GRanges(seqnames=seqnames(y)[i],
									    ranges=IRanges(start=start(y)[i],
											   end=end(y)[(i+1)]),
									    strand="*")
							     ranges<-c(ranges, range)
						     }
						     return(ranges)
				  }))
		       return(out)
		 })


scores_by_ranges<-lapply(ranges, function(x){
				 out<-lapply(x, function(y){
						     scores<-subsetByOverlaps(tiled_chr_all, y)$score
						     return(scores)
})
				 return(out)
		 })

mean_scores_by_ranges<-lapply(scores_by_ranges, function(x){
				      out<-do.call(c, lapply(x, mean))
		 })

grbs_raw<-list()

grbs_raw[["fwrd"]]<-reduce(ranges$fwrd[which(mean_scores_by_ranges$fwrd > quantile(mean_scores_by_ranges$fwrd, percentile))])

grbs_raw[["rev"]]<-reduce(ranges$rev[which(mean_scores_by_ranges$rev > quantile(mean_scores_by_ranges$rev, percentile))])


combined_raw_grbs<-intersect(Reduce(c, grbs_raw$fwrd), Reduce(c, grbs_raw$rev))

#if(species1 != "OikDioicaNorway"){
#	combined_raw_grbs<-keepStandardChromosomes(combined_raw_grbs)
#}


outDir<-paste0("/mnt/biggles/csc_projects/an1413/kurtosis/", species1, "/binned_quantile_kurtosis_genome_ave_grbs/")

dir.create(outDir, showWarnings=F)

export.bed(combined_raw_grbs, paste0(outDir, species1, "_", species2, "_binned_kurtosis_arl_", ARL0, "_", percentile, "_pc_", window_size/1000, "_kb_windows_grbs.bed"))

message("Script Complete")
