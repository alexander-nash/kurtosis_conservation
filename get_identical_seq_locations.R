library(rtracklayer)
library(Hmisc)
library(GenomicRanges)
library(data.table)
library(stringr)
library(Biostrings)
library(CNEr)

###################

getLengthsOfIdenticalSeqs<-function(x){
	#x must be an Axt object imported by CNEr
	qr<-DNAStringSet(apply(as.data.frame(querySeqs(x)), 2, function(z){
				       tmp<-paste0("M", z)
				       out<-paste0(tmp, "N")
				}))
	tr<-DNAStringSet(apply(as.data.frame(targetSeqs(x)), 2, function(z){
				       tmp<-paste0("K", z)
				       out<-paste0(tmp, "V")
				}))

	identical_bases <- mapply(x=qr, y=tr, function(x,y){
					 out<-which(as.raw(x)==as.raw(y))
				})

	relative_ranges <- lapply(identical_bases, function(x){
					  reduce(IRanges(x-2, x-2))
				})

	absolute_ranges <- mapply(z=relative_ranges, y=as.list(x), function(z,y){
					  if(length(z) == 0){
					    out<-list(data.frame())
					  } else {
					    out<-list(data.frame("first.seqnames"=seqnames(first(y)),
					                    "first.start"=(start(first(y)) + start(z)),
					                    "first.end"=(start(first(y)) + end(z)),
					                    "first.width"=(end(z) - start(z) + 1),
					                    "first.strand"="*",
					                    "second.seqnames"=seqnames(second(y)),
					                    "second.start"=(start(second(y)) + start(z)),
					                    "second.end"=(start(second(y)) + end(z)),
					                    "second.width"=(end(z) - start(z) + 1),
					                    "second.strand"="*"))
					  }
				})

	return(absolute_ranges)
}


####################################Start of Script #####################



args<-commandArgs(T)

species1<-args[1]
species2<-args[2]
chrom<-args[3]


speciesDir<-paste0("/mnt/biggles/csc_projects/an1413/kurtosis/", species1)
outDir<-paste0("/mnt/biggles/csc_projects/an1413/kurtosis/", species1, "/identical_seqs/")

if(!dir.exists(file.path(speciesDir))){
	dir.create(file.path(speciesDir))
}


if(!dir.exists(file.path(outDir))){
	dir.create(file.path(outDir))
}

###
if(species2 != "lepOcu1"){
	axts<-lapply(species2, function(x) {
		     lel<-dir(paste0("/mnt/biggles/data/UCSC/axtNet/", species1, "/"), pattern=paste0(x, ".*.axt"), full.names=T)
		     lel<-lel[!grepl("Exon", lel)]
		     lel<-lel[!grepl("broken", lel)]
		     tfn<-paste0("/export/data/goldenpath/", species1, "/assembly.2bit")
		     if(!file.exists(tfn)) tfn<-paste0("/mnt/biggles/data/UCSC/goldenpath/", species1, "/bigZips/", species1, ".2bit")
		     qfn<-paste0("/export/data/goldenpath/", species2, "/assembly.2bit")
		     if(!file.exists(qfn)) qfn<-paste0("/mnt/biggles/data/UCSC/goldenpath/", species2, "/bigZips/", species2, ".2bit")
		     out<-readAxt(lel, tAssemblyFn=tfn, qAssemblyFn=qfn)
		     })
} else {
	axts<-lapply("LepOcu1", function(x) {
		     lel<-dir(paste0("/mnt/biggles/data/UCSC/axtNet/", species1, "/"), pattern=paste0(x, ".*.axt"), full.names=T)
		     lel<-lel[!grepl("Exon", lel)]
		     lel<-lel[!grepl("broken", lel)]
		     tfn<-paste0("/export/data/goldenpath/", species1, "/assembly.2bit")
		     qfn<-paste0("/export/data/goldenpath/", species2, "/assembly.2bit")
		     out<-readAxt(lel, tAssemblyFn=tfn, qAssemblyFn=qfn)
		     })
}

names(axts)<-species2

print("Axts Read")

if(args[4] == "TRUE"){
	chr_axts<-axts[[1]][which(seqnames(first(axts[[1]])) == chrom)]
} else {
	chr_axts<-axts[[1]]
}
rm(axts)

axtSubIndex<-cut2(1:length(chr_axts), m=5000)
print(paste0("number of cuts: ", length(levels(axtSubIndex))))

chr_axts_split<-split(chr_axts, axtSubIndex)

rm(chr_axts)


for(i in 1:length(chr_axts_split)){
  seqs<-getLengthsOfIdenticalSeqs(chr_axts_split[[i]])
  for(j in 1:length(seqs)){
    if(j==1){
	    write.table(seqs[[j]], paste0(outDir, species1, "_", species2, "_", chrom, "_", i, "_identical_seqs.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
    } else {
	    write.table(seqs[[j]], paste0(outDir, species1, "_", species2, "_", chrom, "_", i, "_identical_seqs.tsv"), sep="\t", col.names=F, row.names=F, append=T, quote=F)
    }
  }
  rm(seqs)
  gc()
  print(i)
}

print("Script Complete")
