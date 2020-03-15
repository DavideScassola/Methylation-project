library(BSgenome.Mmusculus.UCSC.mm10)
require(Biostrings)
require(parallel)

get_sites <- function(pattern = "CG")
{
  Find <- function(Genome, Cores){
    if (class(Genome) != "BSgenome") stop("Genome must be a BSgenome!")
    
    XpY <- mclapply(seqlevels(Genome), function(x) start(matchPattern(pattern, Genome[[x]])), mc.cores = Cores)
    return(
      suppressWarnings(
        do.call(c, mclapply(1:length(seqlevels(Genome)), function(x) GRanges(names(Genome)[x], 
                                                                             IRanges(XpY[[x]], width = 2)
        ), mc.cores=Cores))
      )
    )
  }
  
  mm10.XpY <- Find(Genome = BSgenome.Mmusculus.UCSC.mm10, Cores = 2)
  
  total_length <- sum(mm10.XpY@seqnames@lengths[1:21])
  chr_name <- rep_len("",total_length)
  j <- 1
  
  for(i in 1:21) {
    l   <- mm10.XpY@seqnames@lengths[i]
    chr <- mm10.XpY@seqinfo@seqnames[i]
    chr_name[j:(j+l-1)] = chr
    j = j + l
  }
  
  position <- mm10.XpY@ranges@start[1:total_length]
  
  sites <- data.frame(chr_name, position)
  return(sites)
}


pattern_positions <- function(chromosome, pattern, Genome = BSgenome.Mmusculus.UCSC.mm10)
{
  ranges <- matchPattern(pattern,(Genome[[chromosome]]))
  return(IRanges(ranges)@start)
}






Find_CpG <- function(Genome, Cores){
  if (class(Genome) != "BSgenome") stop("Genome must be a BSgenome!")
  
  CpG <- mclapply(seqlevels(Genome), function(x) start(matchPattern("CG", Genome[[x]])), mc.cores = Cores)
  return(
    suppressWarnings(
      do.call(c, mclapply(1:length(seqlevels(Genome)), function(x) GRanges(names(Genome)[x], 
                                                                           IRanges(CpG[[x]], width = 2)
      ), mc.cores=Cores))
    )
  )
}

mm10.CpG <- Find_CpG(Genome = BSgenome.Mmusculus.UCSC.mm10, Cores = 2)

setwd("~/Scrivania/Tesi")
load("mm10.CpG.RData")

total_length <- sum(mm10.CpG@seqnames@lengths[1:21])
chr_name <- rep_len("",total_length)
j <- 1

for(i in 1:21) {
  l   <- mm10.CpG@seqnames@lengths[i]
  chr <- mm10.CpG@seqinfo@seqnames[i]
  chr_name[j:(j+l-1)] = chr
  j = j + l
}

position <- mm10.CpG@ranges@start[1:total_length]

sites <- data.frame(chr_name, position)
