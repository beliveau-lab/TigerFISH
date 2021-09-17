
#!/usr/bin/Rscript

install.packages("chromoMap",repos="http://cran.us.r-project.org")

library(chromoMap) 
library(htmlwidgets)
library(optparse)
library(tidyverse)
library(splitstackshape)

option_list = list(
  make_option(c("-c", "--chrom_sizes"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-r", "--repeat_annotation"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$repeat_annotation)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


#read in probe file
probe_repeat <- read.table(opt$repeat_annotation, 
                           sep = "\t", header=FALSE)

#name columns
colnames(probe_repeat)<-c("probe_coords","repeat_coords","probe")

probe_repeat$probe_coords<-NULL
probe_repeat$probe<-NULL

out<-cSplit(probe_repeat, "repeat_coords", sep="-")

out<-cSplit(out, "repeat_coords_1", sep=":")

renamed_regions <- out %>% 
  rename(
    chrom = repeat_coords_1_1,
    start = repeat_coords_1_2,
    stop = repeat_coords_2
    
  )

renamed_regions$annot<-"annot"

#drop and reorder cols
renamed_regions<-renamed_regions[,c("annot","chrom","start","stop")]
renamed_regions <- as.data.frame(renamed_regions)


#read in chrom length file
chrom_sizes <- read.table(opt$chrom_sizes, 
                           sep = "\t",header=FALSE)

#name columns
colnames(chrom_sizes)<-c("chrom","stop")

#add start of chrom
chrom_sizes$start<-1

chrom_sizes<-renamed_regions[,c("chrom","start","stop")]

#generate chromomap
map<-chromoMap(list(chrom_sizes),list(renamed_regions))

#isave as html
saveWidget(map,opt$out)
