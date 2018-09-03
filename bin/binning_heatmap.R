for(package in c('ggplot2','plyr','reshape2','cowplot')) {
    suppressPackageStartupMessages(do.call(library, list(package)))
}

args = commandArgs(TRUE)

binning_folder = args[1]
cov_stats_folder = args[2]
ani_species = args[3]
conta_thresh = as.numeric(args[4])
# min_bin_size = as.numeric(args[5])

conta_thresh_checkm = conta_thresh * 100

ani_max_per_species = read.table(ani_species,header=T,stringsAsFactors = F,sep=",")

prec_rec_files = system(sprintf("ls %s/*/prec_rec_dominant.tsv",binning_folder),intern=T)

cov_files = system(sprintf("ls %s/Ech*stat.txt",cov_stats_folder),intern=T)

prec_rec_comb = data.frame(bin_name=character(),Annot_maj_sp=character(),
                           Contamination=double(),Completeness=double(),soft=character())

for(file in prec_rec_files) {
    prec_rec = read.table(file,header=T, stringsAsFactors = F)
    prec_rec = prec_rec[which(prec_rec$Contamination <= conta_thresh ),] #& prec_rec$length_bin_annotated_sequence >= min_bin_size
    prec_rec = prec_rec[order(prec_rec$Annot_maj_sp, -prec_rec$Completeness),]
    prec_rec_uniq_annot = prec_rec[!duplicated(prec_rec$Annot_maj_sp),]
    binning_soft = unlist(strsplit(file,"/"))[length(unlist(strsplit(file,"/"))) - 1]
    soft = rep(binning_soft,nrow(prec_rec_uniq_annot))
    prec_rec_uniq_annot = cbind(prec_rec_uniq_annot,soft)
    prec_rec_comb = rbind(prec_rec_comb, prec_rec_uniq_annot)
}

prec_rec_comb = prec_rec_comb[which(prec_rec_comb$soft != "SCIMM"),]

stats_comb = data.frame(Genome=character(), Size=integer(), Coverage=character(), Reads=integer())

for(file in cov_files) {
    stats = read.table(file, header=T, stringsAsFactors = F)
    stats_comb = rbind(stats_comb, stats)
}

stats_comb$Coverage = as.numeric(gsub("X","",stats_comb$Coverage))

sum_cov_per_genome = ddply(stats_comb, c("Genome"), function(x) sum(x$Coverage))

names(sum_cov_per_genome)[2] = "cov_sum"

prec_rec_cov = merge(prec_rec_comb, sum_cov_per_genome, by.x="Annot_maj_sp", by.y="Genome")
prec_rec_cov_ani = merge(prec_rec_cov, ani_max_per_species[,c("genome1", "ANI", "grp_abun")],
                         by.x="Annot_maj_sp", by.y="genome1")

# if(length(unique(sort(prec_rec_cov$soft))) == 8) {
#     prec_rec_cov$soft = factor(prec_rec_cov$soft, levels = unique(sort(prec_rec_cov$soft))[c(2,8,1,3,4,6,5)])
# } else {
#     prec_rec_cov$soft = factor(prec_rec_cov$soft, levels = unique(sort(prec_rec_cov$soft))[c(1,2,3,4,6,5)])
# }

names(prec_rec_cov_ani)[6] = "Softwares"

# png(filename=paste0(getwd(), "/plot_completeness_fct_cov_filt_conta_10.png"), width=750, height=700)
# 
# ggplot(prec_rec_cov, aes(x=Annot_maj_sp, y=Completeness, group=Softwares, color=Softwares)) +
#     geom_point(aes(shape=Softwares), size=2) +
#     geom_line(aes(linetype=Softwares)) +
#     scale_x_discrete(limits=sum_cov_per_genome[order(sum_cov_per_genome$cov_sum),]$Genome) +
#     scale_shape_manual(values = c(1, 2, 0, 4, 5, 7, 8, 6)) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle=75, hjust=1), axis.title.x = element_blank())
#     #geom_line(aes(x=Annot_maj_sp, y=cov_sum))
# dev.off()

# ggplot(prec_rec_cov, aes(x=Annot_maj_sp, y=cov_sum)) +
#     geom_point() + 
#     geom_line() +
#     scale_x_discrete(limits=sum_cov_per_genome[order(sum_cov_per_genome$cov_sum),]$Genome) +
#     theme(axis.text.x = element_text(angle=75, hjust=1))

# heatmap view

prec_rec_cov_ani[which(prec_rec_cov_ani$ANI > 0), "Annot_maj_sp"] = 
    paste0(prec_rec_cov_ani[which(prec_rec_cov_ani$ANI > 0), "Annot_maj_sp"], " (", prec_rec_cov_ani[which(prec_rec_cov_ani$ANI > 0), "ANI"], ")")
prec_rec_cov_ani[which(prec_rec_cov_ani$ANI == 0), "Annot_maj_sp"] = 
    paste0(prec_rec_cov_ani[which(prec_rec_cov_ani$ANI == 0), "Annot_maj_sp"], " (", prec_rec_cov_ani[which(prec_rec_cov_ani$ANI == 0), "cov_sum"], ")")

png(filename = paste0(getwd(), "/heatmap_completeness_fct_cov_filt_conta_10.png"), width=580, height=650)

#If simulated data have groups of abundance
if(regexpr('random', cov_stats_folder)[1] == -1) {
    levT = as.character(unique(prec_rec_cov_ani[which(prec_rec_cov_ani$grp_abun == "T"), "Annot_maj_sp"]))
    levANI = prec_rec_cov_ani[which(prec_rec_cov_ani$ANI > 0 & prec_rec_cov_ani$grp_abun != "T"),]
    levANI = as.character(unique(levANI[order(levANI$ANI, decreasing = T), "Annot_maj_sp"]))
    levCov = prec_rec_cov_ani[which(prec_rec_cov_ani$ANI == 0 & prec_rec_cov_ani$grp_abun != "T"),]
    levCov = as.character(unique(levCov[order(levCov$cov_sum, decreasing = T), "Annot_maj_sp"]))
    
    prec_rec_cov_ani$Annot_maj_sp = factor(prec_rec_cov_ani$Annot_maj_sp, levels = rev(c(levT, levANI, levCov)))
    
     plot1 = ggplot(prec_rec_cov_ani, aes(x=Softwares, y=Annot_maj_sp)) +
        geom_tile(aes(fill=Completeness), show.legend = F) +
        theme_minimal() +
        scale_x_discrete(limits=c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2")) +
        scale_fill_gradient2(low="#FEFFD5", mid="#FFFF00", high="#FF0000FF", midpoint = 0.5) +
        theme(panel.grid=element_blank(), axis.text.x= element_text(angle=45, hjust=1)) +
        annotate("segment", x=0, xend=8.5, y=37.5, yend=37.5) +
        annotate("segment", x=0, xend=8.5, y=13.5, yend=13.5) +
        ggtitle("BLAST") +
        ylab("") +
        theme(plot.title = element_text(hjust= 0.5))
} else {
    levANI = prec_rec_cov_ani[which(prec_rec_cov_ani$ANI > 0),]
    levANI = as.character(unique(levANI[order(levANI$ANI, decreasing = T), "Annot_maj_sp"]))
    levCov = prec_rec_cov_ani[which(prec_rec_cov_ani$ANI == 0),]
    levCov = as.character(unique(levCov[order(levCov$cov_sum, decreasing = T), "Annot_maj_sp"]))
    
    prec_rec_cov_ani$Annot_maj_sp = factor(prec_rec_cov_ani$Annot_maj_sp, levels = rev(c(levANI, levCov)))
    
    plot1 = ggplot(prec_rec_cov_ani, aes(x=Softwares, y=Annot_maj_sp)) +
        geom_tile(aes(fill=Completeness)) +
        theme_minimal() +
        scale_x_discrete(limits=c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2")) +
        scale_fill_gradient2(low="#FEFFD5", mid="#FFFF00", high="#FF0000FF", midpoint = 0.5) +
        theme(panel.grid=element_blank(), axis.text.x= element_text(angle=45, hjust=1)) +
        annotate("segment", x=0, xend=8.5, y=15.5, yend=15.5) +
        ylab("")
}

dev.off()

files_checkm = system(sprintf("ls %s/*/chkm_res/tree_qa+qa.tsv", binning_folder), intern=T)

tab_rec_uniq_checkm_blast_annot = c()

for(file in files_checkm) {
    soft = strsplit(file, "/")[[1]][length(strsplit(file, "/")[[1]])-2]
    prec_rec_checkm = read.table(file, sep="\t", as.is=T, header=T, comment.char = "")
    prec_rec_blast = read.table(paste0(binning_folder, "/", soft, "/prec_rec_dominant.tsv"), sep="\t", header=T, stringsAsFactors = F)
    prec_rec_checkm_clean = prec_rec_checkm[which(prec_rec_checkm$Contamination <= conta_thresh_checkm),]
    prec_rec_checkm_clean_blast_annotation = merge(prec_rec_checkm_clean, prec_rec_blast[,c("bin_name", "Annot_maj_sp")], by.x = "Bin.id", by.y = "bin_name")
    if(nrow(prec_rec_checkm_clean_blast_annotation) > 0) {
        # annot_sp = sapply(prec_rec_clean$Taxonomy, annot_sp_fct)
        # if(sum(is.na(annot_sp)) < length(annot_sp)) {
        # prec_rec_clean = cbind(prec_rec_clean, annot_sp)
        prec_rec_checkm_clean_ordered = prec_rec_checkm_clean_blast_annotation[order(
            prec_rec_checkm_clean_blast_annotation$Annot_maj_sp, -prec_rec_checkm_clean_blast_annotation$Completeness),]
        prec_rec_checkm_clean_ordered_uniq = prec_rec_checkm_clean_ordered[!duplicated(prec_rec_checkm_clean_ordered$Annot_maj_sp, incomparables = NA),]
        software = rep(soft, nrow(prec_rec_checkm_clean_ordered_uniq))
        tab_rec_uniq_checkm_blast_annot = rbind(tab_rec_uniq_checkm_blast_annot, cbind(prec_rec_checkm_clean_ordered_uniq, software))
        #} else {
        #     annot_sp = rep(NA, nrow(prec_rec_clean))
        #     software = rep(soft, nrow(prec_rec_clean))
        #     tab_rec_uniq_checkm = rbind(tab_rec_uniq_checkm, cbind(prec_rec_clean, annot_sp, software))
        # }
    }
}

tab_rec_uniq_checkm_blast_annot = tab_rec_uniq_checkm_blast_annot[which(tab_rec_uniq_checkm_blast_annot$software != "SCIMM"),]

tab_rec_uniq_checkm_blast_annot_cov = merge(tab_rec_uniq_checkm_blast_annot, sum_cov_per_genome, by.x="Annot_maj_sp", by.y="Genome")
tab_rec_uniq_checkm_blast_annot_cov_ani = merge(tab_rec_uniq_checkm_blast_annot_cov, ani_max_per_species[,c("genome1", "ANI", "grp_abun")],
                         by.x="Annot_maj_sp", by.y="genome1")

names(tab_rec_uniq_checkm_blast_annot_cov_ani)[17] = "Softwares"

tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI > 0), "Annot_maj_sp"] = 
    paste0(tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI > 0), "Annot_maj_sp"],
           " (", tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI > 0), "ANI"], ")")
tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI == 0), "Annot_maj_sp"] = 
    paste0(tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI == 0), "Annot_maj_sp"], 
           " (", tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI == 0), "cov_sum"], ")")

png(filename = paste0(getwd(), "/heatmap_completeness_CheckM_fct_cov_filt_conta_10.png"), width=580, height=650)

#If simulated data have groups of abundance
if(regexpr('random', cov_stats_folder)[1] == -1) {
    levT = as.character(unique(tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$grp_abun == "T"), "Annot_maj_sp"]))
    levANI = tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI > 0 & tab_rec_uniq_checkm_blast_annot_cov_ani$grp_abun != "T"),]
    levANI = as.character(unique(levANI[order(levANI$ANI, decreasing = T), "Annot_maj_sp"]))
    levCov = tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI == 0 & tab_rec_uniq_checkm_blast_annot_cov_ani$grp_abun != "T"),]
    levCov = as.character(unique(levCov[order(levCov$cov_sum, decreasing = T), "Annot_maj_sp"]))
    
    tab_rec_uniq_checkm_blast_annot_cov_ani$Annot_maj_sp = factor(tab_rec_uniq_checkm_blast_annot_cov_ani$Annot_maj_sp, levels = rev(c(levT, levANI, levCov)))
    
    plot2 = ggplot(tab_rec_uniq_checkm_blast_annot_cov_ani, aes(x=Softwares, y=Annot_maj_sp)) +
        geom_tile(aes(fill=Completeness)) +
        theme_minimal() +
        scale_x_discrete(limits=c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2")) +
        scale_fill_gradient2(low="#FEFFD5", mid="#FFFF00", high="#FF0000FF", midpoint = 0.5) +
        theme(panel.grid=element_blank(), axis.text.x= element_text(angle=45, hjust=1)) +
        annotate("segment", x=0, xend=8.5, y=37.5, yend=37.5) +
        annotate("segment", x=0, xend=8.5, y=13.5, yend=13.5) +
        ylab("") +
        ggtitle("CheckM") +
        theme(plot.title = element_text(hjust= 0.5), axis.text.y = element_blank())
} else {
    levANI = tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI > 0),]
    levANI = as.character(unique(levANI[order(levANI$ANI, decreasing = T), "Annot_maj_sp"]))
    levCov = tab_rec_uniq_checkm_blast_annot_cov_ani[which(tab_rec_uniq_checkm_blast_annot_cov_ani$ANI == 0),]
    levCov = as.character(unique(levCov[order(levCov$cov_sum, decreasing = T), "Annot_maj_sp"]))
    
    tab_rec_uniq_checkm_blast_annot_cov_ani$Annot_maj_sp = factor(tab_rec_uniq_checkm_blast_annot_cov_ani$Annot_maj_sp, levels = rev(c(levANI, levCov)))
    
    plot2 = ggplot(tab_rec_uniq_checkm_blast_annot_cov_ani, aes(x=Softwares, y=Annot_maj_sp)) +
        geom_tile(aes(fill=Completeness)) +
        theme_minimal() +
        scale_x_discrete(limits=c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2")) +
        scale_fill_gradient2(low="#FEFFD5", mid="#FFFF00", high="#FF0000FF", midpoint = 0.5) +
        theme(panel.grid=element_blank(), axis.text.x= element_text(angle=45, hjust=1)) +
        annotate("segment", x=0, xend=8.5, y=15.5, yend=15.5) +
        ylab("")
}

dev.off()

# grid.arrange(plot1, plot2, ncol=2)
plot_grid(plot1, plot2, rel_widths = c(1.4,1))
