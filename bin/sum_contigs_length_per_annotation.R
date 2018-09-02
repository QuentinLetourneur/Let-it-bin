library(plyr)

args = commandArgs(TRUE)

blast_best_hit = args[1]
out = args[2]

blast_bh = read.table(blast_best_hit, header=F,stringsAsFactors = F)

sum_contig_length_per_annotation = ddply(blast_bh, c("V5"), function(x) sum(x$V2))

sum_contig_gt1000_per_annotation = ddply(blast_bh[which(blast_bh$V2 >= 1000),], c("V5"), function(x) sum(x$V2))

res = merge(sum_contig_length_per_annotation, sum_contig_gt1000_per_annotation, by = 1)

names(res) = c("species","sum contigs length","sum contigs length >= 1kb")

write.table(res, file=out, quote=F, row.names=F, sep="\t")