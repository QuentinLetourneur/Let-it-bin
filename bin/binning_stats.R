for(package in c('ggplot2','plyr','reshape2')) {
    suppressPackageStartupMessages(do.call(library, list(package)))
}

args = commandArgs(TRUE)

# type = 'MetaBAT'
bin_file = args[1] #"/mnt/atlas/my_home/binning_wf/Bacteria_V4/bins/bin"
# bin_file = "/mnt/atlas/my_home/binning_wf/Bacteria_V3/bins/bin"
#file = "/pasteur/homes/qletourn/binning_wf/Bacteria_V3/bins/bin"

out = args[2]
annotDir = args[3]
refsinfo = args[4]
min_bin_size_for_plot = as.numeric(args[5])

#bin files names
files <- system(sprintf("ls %s.* | egrep '\\.[0-9]+\\.fa$'", bin_file), intern=T)
if(length(files) == 0)
    stop(sprintf("Cannot find bins: %s.*", file))

refs_info = read.table(refsinfo, as.is=T, header=T, 
                       sep="\t", strip.white = T, stringsAsFactors = F)

# annot strain
annot_maj_fun = function(refs_info, annot) {
    subset_genome = annot[annot[,3] %in% trimws(refs_info[3]),]
    subset_plasmids = annot[annot[,3] %in% strsplit(refs_info[4],";")[[1]],]
    res=c(refs_info[2], sum(subset_genome[,2]), sum(subset_plasmids[,2]), 
          sum(subset_genome[,2]) + sum(subset_plasmids[,2]))
    return(res)
}

# class_fun = function(annot, res_stat, refs) {
#     if(annot %in% refs) {
#         if(nrow(res_stat[which(res_stat[,2] == annot),]) > 1) {
#             class = c("1", rep("2", nrow(res_stat[which(res_stat == annot),]) - 1))
#         }
#         else if(nrow(res_stat[which(res_stat == annot),]) == 1) {
#             class = "1"
#         }
#     }
#     else {
#         class = "2"
#     }
#     return(class)
# }

# res_stat_strain = data.frame(bin_name=character(), Annot_maj=character(),Prec=double(),
                             # Rec=double(),TP_genome=integer(), TP_plasmids=integer(), 
                             # TF=integer(), stringsAsFactors=F)

# for (file in files) {
#     bin_name = system(sprintf("echo %s | egrep -o 'bin\\.[0-9]+'",file), intern=T)
#     # annot = read.table(paste("/mnt/atlas/my_home/binning_wf/Bacteria_V3/Blast_on_refs/Annot/",
#     #                          bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")
#     annot = read.table(paste("/mnt/atlas/my_home/binning_wf/Bacteria_V4/Annot/",
#                              bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")
#     #annot = read.table(paste("/mnt/atlas/my_home/binning_wf/Bacteria_V3/Annot/",
#     #                         bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")
#     
#     annot_maj_strain = as.data.frame(t(apply(refs_info, 1, function(x,y) annot_maj_fun(x,y),
#                                              y=annot)), stringsAsFactors=F)
#     
#     annot_other = annot[!(annot[,3] %in% strsplit(paste(refs_info[,3], 
#                                                         refs_info[,4], collapse=";",
#                                                         sep=";"), ";")[[1]]),c(2,3)]
#     
#     #annot_maj_other = ddply(annot_other, .(V3), function(x) sum(x$V2))
#     #colnames(annot_maj_other) = c("V1","V2")
#     #annot_maj_other = data.frame(annot_maj_other, V3=rep(0, nrow(annot_maj_other)),
#     #                             V4=rep(0, nrow(annot_maj_other)))
#     
#     annot_maj = annot_maj_strain
#     #annot_maj = rbind(annot_maj_strain, annot_maj_other)
#     annot_maj[,2] = as.numeric(annot_maj[,2])
#     annot_maj[,3] = as.numeric(annot_maj[,3])
#     annot_maj[,4] = as.numeric(annot_maj[,4])
#     annot_maj = annot_maj[order(annot_maj[,4], decreasing = T),]
#     
#     #print(bin_name)
#     #print(annot_maj)
#     
#     TP_genome = annot_maj[1,2]
#     TP_plasmids = annot_maj[1,3]
#     TP = TP_genome + TP_plasmids
#     TF = sum(annot_maj[2:nrow(annot_maj),2]) + sum(annot_maj[2:nrow(annot_maj),3])
#     Prec = TP / (sum(annot_maj[,2]) + sum(annot_maj[,3]))
#     if(annot_maj[1,1] %in% refs_info[,2]) {
#         Rec = TP / as.numeric(refs_info[which(annot_maj[1,1] == refs_info[,2]),7])
#     }
#     else {
#         Rec = 0
#     }
#     
#     res_stat_strain = rbind(res_stat_strain, data.frame("bin_name"=bin_name, 
#                                                         "Annot_maj"=annot_maj[1,1],
#                                                         "Prec"=Prec, "Rec"=Rec, 
#                                                         "TP_genome"=TP_genome, 
#                                                         "TP_plasmids"=TP_plasmids, 
#                                                         "TF"=TF,
#                                                         stringsAsFactors=F))
# }
# 
# res_stat_strain = data.frame(res_stat_strain, 
#                         fac=unlist(apply(as.matrix(unique(res_stat_strain[,2])), 1, 
#                                      function(x,y,z) class_fun(x,y,z),
#                                      y=res_stat_strain, z=refs_info[,2])))

# res_stat_strain$TF[is.na(res_stat_strain$TF)] = 0

# annot species
# if ls fasta_.../strain* | wc -l == 1 :
# calc Rec & Prec
res_stat_sp = data.frame(bin_name=character(), Annot_maj_sp=character(), Prec_sp=double(),
                         Rec_sp=double(), Prec_strain1=double(), Prec_plasm1=double(),
                         Prec_strain2=double(), Prec_plasm2=double(), Rec_strain1=double(),
                         Rec_plasm1=double(), Rec_strain2=double(), Rec_plasm2=double(),
                         TF_sp=double(), bin_size=integer(), mean_genome_size=integer(),
                         stringsAsFactors=F)

for (file in files) {
    sp_annot_maj = "NA"
    col = 11
    
    bin_name = system(sprintf("echo %s | egrep -o 'bin\\.[0-9]+'",file), intern=T)
    # annot = read.table(paste("/mnt/atlas/my_home/binning_wf/Bacteria_V3/Blast_on_refs/Annot/",
    #                          bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")
    
    annot = read.table(paste(annotDir, "/", bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")
    while (sp_annot_maj == "NA" && col >= 4) {
        
        annot_maj_sp = ddply(annot, c(paste("V", as.character(col), sep="")), 
                             function(x) sum(x$V2))
        annot_maj_sp = annot_maj_sp[order(annot_maj_sp[,2], decreasing=T),]
        sp_annot_maj = annot_maj_sp[1,1]
        
        col = col - 1
    }
    
    sub_refs_info = refs_info[which(refs_info[,1] == gsub(" ", "_", sp_annot_maj)),]
    
    count_by_strain = 
        as.data.frame(t(apply(sub_refs_info, 1, function(x,y) annot_maj_fun(x,y),
                              y=annot[which(annot$V12 >= 95),])), stringsAsFactors=F)
                              # to concider only match that are annotated at species level
    count_by_strain = count_by_strain[order(as.numeric(count_by_strain$V4), decreasing = T),]
    
    bin_size = sum(annot_maj_sp[,2])
    mean_genome_size = mean(sub_refs_info$total_length)
    
    count_by_strain_cp = count_by_strain
    
    count_by_strain$V2 = (as.double(count_by_strain$V2) / bin_size) * 100
    count_by_strain$V3 = (as.double(count_by_strain$V3) / bin_size) * 100
    
    count_by_strain_cp$V2 = (as.double(count_by_strain_cp$V2) / mean_genome_size) * 100
    count_by_strain_cp$V3 = (as.double(count_by_strain_cp$V3) / mean_genome_size) * 100
    
    #print(bin_name)
    #print(count_by_strain)
    TP_sp = annot_maj_sp[1,2]
    TF_sp = (sum(annot_maj_sp[2:nrow(annot_maj_sp),2]) / bin_size) * 100
    Prec_sp = TP_sp / bin_size
    Rec_sp = TP_sp / mean_genome_size
    #as.numeric(system(sprintf("grep -v '^>' /mnt/atlas/my_home/fasta_genome_bact/%s*complete_genome.fa | wc -c",gsub(" ", "_",sp_annot_maj)), intern=T))
    if(nrow(count_by_strain) == 1) {
        res_stat_sp = rbind(res_stat_sp, 
                            data.frame("bin_name"=bin_name, "Annot_maj_sp"=sp_annot_maj,
                                       "Prec_sp"=Prec_sp, "Rec_sp"=Rec_sp, "Prec_strain1"=as.numeric(count_by_strain[1,2]),
                                       "Prec_plasm1"=count_by_strain[1,3], "Prec_strain2"=0, "Prec_plasm2"=0,
                                       "Rec_strain1"=count_by_strain_cp[1,2], "Rec_plasm1"=count_by_strain[1,3],
                                       "Rec_strain2"=0, "Rec_plasm2"=0, "TF_sp"=TF_sp, "bin_size"=bin_size,
                                       "mean_genome_size"=mean_genome_size, stringsAsFactors=F))
    } else {
        res_stat_sp = rbind(res_stat_sp,
                            data.frame("bin_name"=bin_name, "Annot_maj_sp"=sp_annot_maj,
                                       "Prec_sp"=Prec_sp, "Rec_sp"=Rec_sp, "Prec_strain1"=count_by_strain[1,2],
                                       "Prec_plasm1"=count_by_strain[1,3], "Prec_strain2"=count_by_strain[2,2],
                                       "Prec_plasm2"=count_by_strain[2,3], "Rec_strain1"=count_by_strain_cp[1,2],
                                       "Rec_plasm1"=count_by_strain_cp[1,3], "Rec_strain2"=count_by_strain_cp[2,2],
                                       "Rec_plasm2"=count_by_strain_cp[2,3], "TF_sp"=TF_sp, "bin_size"=bin_size,
                                       "mean_genome_size"=mean_genome_size, stringsAsFactors=F))
    }
    
}

# res_stat_sp = cbind(res_stat_sp, 
#                     unlist(apply(as.matrix(unique(res_stat_sp[,2])), 1, 
#                                  function(x,y,z) class_fun(x,y,z),
#                                  y=res_stat_sp, z=refs_info$V1)))

res_stat_sp = res_stat_sp[order(res_stat_sp$Annot_maj_sp, res_stat_sp$bin_size,
                                decreasing=T),]
res_stat_sp$TF_sp[is.na(res_stat_sp$TF_sp)] = 0

write.table(res_stat_sp[,c(1,2,3,4)], paste0(out, "/prec_rec.tsv"), sep = "\t", row.names = F)

res_sp_trunc = res_stat_sp[which(res_stat_sp$bin_size >= min_bin_size_for_plot),]

res_sp_trunc$bin_name = paste(res_sp_trunc$bin_name, res_sp_trunc$Annot_maj_sp, sep = " ")

# barplot_prec = melt(res_stat_strain,id.vars="bin_name", measure.vars = c("TP_genome","TP_plasmids","TF"),
#                     variable.name = "fac", value.name = "nb_base")
# barplot_prec = barplot_prec[order(barplot_prec[,2], barplot_prec[,1]),]
# sum_TP_TF = ddply(barplot_prec, "bin_name", function(x) sum(x[,3]))
# sum_TP_TF = sum_TP_TF[order(sum_TP_TF$bin_name),]

# res_stat_strain = data.frame(res_stat_strain, sum_TPTF=sum_TP_TF[,2])
# res_stat_strain = res_stat_strain[order(res_stat_strain$TP_genome + res_stat_strain$TP_plasmids,
                                        # decreasing=T),]

#barplot with annot
# ggplot(data=barplot_prec, aes(x=bin_name, y=nb_base, fill=fac)) + 
#     geom_col(position = position_stack(reverse = TRUE)) + 
#     geom_text(aes(label=c(res_stat_strain[order(res_stat_strain$sum_TPTF, decreasing = T),2],rep("",nrow(res_stat_strain)*2)),
#                  y=max(res_stat_strain$sum_TPTF)-((max(res_stat_strain$sum_TPTF) * 20)/100))) +
    # geom_text(aes(label=c(as.character(round(res_stat_strain[order(res_stat_strain$sum_TPTF, decreasing = T),]$Prec,2)),
    #                       rep("",nrow(res_stat_strain)*2)),
    #              y=c(res_stat_strain[order(res_stat_strain$sum_TPTF, decreasing = T),]$TP/2,rep(0,nrow(res_stat_strain)*2)))) +
    # coord_flip() +
    # scale_fill_manual(values=c("#30AEA5","#30AEA5","#F94D54")) +
    # guides(fill=guide_legend(title=NULL)) +
    # theme(axis.text.x= element_text(size="11"), axis.text.y= element_text(size="12")) +
    # scale_x_discrete(limits=rev(as.factor(res_stat_strain[order(res_stat_strain$sum_TPTF,
                                                                # decreasing = T),]$bin_name)))

# pos_text_rec=0.03
# pos_text_annot=0.75

# barplot_rec = melt(res_stat_strain[order(res_stat_strain$Rec, decreasing = T),],
#                    id.vars=c("bin_name","Annot_maj"), measure.vars = c("TP_genome", "TP_plasmids"),
#                    variable.name = "fac", value.name = "TP")

# barplot_rec = data.frame(barplot_rec, order=seq(nrow(barplot_rec)))

# barplot_rec = merge(barplot_rec, refs_info[,c(2,7)], by.x=2, by.y=1)
# barplot_rec = barplot_rec[order(barplot_rec$fac),]
# 
# Completeness = (barplot_rec$TP / barplot_rec$total_length) * 100
# 
# barplot_rec=data.frame(barplot_rec,Completeness=Completeness)

# ggplot(data=barplot_rec[,2:6], aes(x=bin_name,y=Completeness, fill=fac)) +
#     geom_col(position = position_stack(reverse = TRUE)) +
    # geom_text(aes(label=as.matrix(round(res_stat_strain$Rec,2)), y=pos_text_rec)) +
    # geom_text(aes(label=res_stat_strain$Annot_maj),y=pos_text_annot) +
    # coord_flip() +
    # theme(axis.text.x= element_text(size="11"), axis.text.y= element_text(size="12")) +
    # scale_x_discrete(limits=rev(as.factor(res_stat_strain[order(res_stat_strain$Rec,
    #                                                             decreasing = T),]$bin_name)))

barplot_prec_sp = melt(res_sp_trunc,id.vars="bin_name",
                       measure.vars = c("Prec_strain1", "Prec_plasm1", "Prec_strain2",
                                        "Prec_plasm2","TF_sp"),
                       variable.name = "fac", value.name = "nb_base")

# barplot_prec_sp = barplot_prec_sp[order(barplot_prec_sp[,2], 
#                                         barplot_prec_sp[,6] + barplot_prec_sp[,8],
#                                         decreasing = T),]

sum_TP_TF_sp = ddply(barplot_prec_sp, "bin_name", function(x) sum(x[,3]))[,2]
#barplot with annot

png(filename=paste0(getwd(), "/contamination.png"), width=750, height=600)

ggplot(data=barplot_prec_sp, aes(x=bin_name, y=nb_base, fill=fac)) + 
    geom_col(position = position_stack(reverse = TRUE)) + 
    # geom_text(aes(label=c(res_stat_sp[order(res_stat_sp[,1]),2],rep("",nrow(res_stat_sp))),
    #               y=max(sum_TP_TF_sp)-((max(sum_TP_TF_sp) * 20)/100))) + 
    # geom_text(aes(label=c(as.character(round(res_stat_sp[order(res_stat_sp$bin_name),]$Prec,2)),
    #                       rep("",nrow(res_stat_sp))),
    #               y=c(res_stat_sp[order(res_stat_sp$bin_name),]$TP/2,rep(0,nrow(res_stat_sp))))) +
    coord_flip() +
    scale_fill_manual(values=c("#2A4BE9","#1E90FF","#11A121","#5DDD6B","#F94D54"),
                      labels = c("Strain 1","Plasmids 1","Strain 2","Plasmids 2","Contamination")) +
    guides(fill=guide_legend(title=NULL)) +
    theme(axis.text.x= element_text(size="15"), axis.text.y= element_text(size="14"),
          legend.position = "bottom", legend.text=element_text(size=11)) + 
    scale_x_discrete(limits=as.factor(res_sp_trunc[order(res_sp_trunc$Rec_sp),]$bin_name))

dev.off()

barplot_rec_sp = melt(res_sp_trunc, id.vars=c("bin_name", "mean_genome_size"), 
                      measure.vars = c("Rec_strain1","Rec_plasm1","Rec_strain2","Rec_plasm2"),
                      variable.name="fac", value.name = "TP")

barplot_rec_sp$TP = as.numeric(barplot_rec_sp$TP)

Completeness_sp = (barplot_rec_sp$TP / barplot_rec_sp$mean_genome_size) * 100

barplot_rec_sp = data.frame(barplot_rec_sp, Completeness_sp=Completeness_sp)

# fill_col_name = as.factor(apply(barplot_rec_sp,1,function(x) paste(x[1],x[3],sep="")))
# 
# fill_col=c("gold2","gold2","gold2","darkgoldenrod4","bisque3","bisque3","chartreuse4",
#            "cyan3","darkorange3","darkgoldenrod4","blue","gold2","darkorchid4",
#            "darkred","dodgerblue","yellow","yellow","yellow","darkgoldenrod1",
#            "antiquewhite","antiquewhite","chartreuse","cadetblue","darkorange",
#            "darkgoldenrod1","aquamarine","yellow","darkorchid","red","deepskyblue",
#            "dodgerblue2","dodgerblue2","dodgerblue2","bisque3","darkgoldenrod4",
#            "darkgoldenrod4","NA","NA","NA","bisque3","NA","dodgerblue2","NA","NA",
#            "gold2","deepskyblue","deepskyblue","deepskyblue","antiquewhite","darkgoldenrod1",
#            "darkgoldenrod1","NA","NA","NA","antiquewhite","NA","deepskyblue","NA",
#            "NA","yellow")

# names(fill_col) = fill_col_name

png(filename=paste0(getwd(), "/completeness.png"), width=750, height=600)

ggplot(data=barplot_rec_sp, aes(x=bin_name,y=TP, fill=fac)) + #, colour=fill_cols 
    geom_col(position = position_stack(reverse = TRUE)) + #, size=2, width=0.7
    #geom_text(aes(label=as.matrix(round(res_stat_sp$Rec_sp,2)), y=pos_text_rec)) +
    #geom_text(aes(label=res_stat_sp$Annot_maj_sp),y=pos_text_annot) +
    coord_flip() +
    guides(fill=guide_legend(title=NULL)) +
    theme(axis.text.x= element_text(size="14"), axis.text.y= element_text(size="15"),
          legend.position="bottom", legend.text=element_text(size=11)) +
    labs(y="Completeness") +
    scale_x_discrete(limits=as.factor(res_sp_trunc[order(res_sp_trunc$Rec_sp),]$bin_name)) +
    scale_fill_manual(values=c("#2A4BE9","#1E90FF","#11A121","#5DDD6B"),
                      labels = c("Strain 1","Plasmids 1","Strain 2","Plasmids 2")) #+
    # scale_colour_manual(values=fill_col)

dev.off()

# c("blue","aquamarine","NA","NA","charteuse4","chartreuse",
#   "NA","NA","darkorchid4","darkorchid1","darkorange3",
#   "darkorange","NA","NA","darkred","red","NA","NA",
#   "cyan3","cadetblue1","NA","NA","darkgoldenroad4",
#   "darkgoldenroad1","bisque3","antiquewhite","gold2",
#   "yellow","dodgerblue2","deepskyblue","dodgerblue2",
#   "deepskyblue","gold2","yellow","gold2","yellow",
#   "dodgerblue2","deepskyblue","gold2","yellow",
#   "dodgerblue2","deepskyblue","gold2","yellow",
#   "dodgerblue2","deepskyblue","darkgoldenroad4",
#   "darkgoldenroad1","bisque3","antiquewhite","bisque3",
#   "antiquewhite","darkgoldenroad4","darkgoldenroad1")

# compo plasmids
# bin_name = "bin.15"
# annot = read.table(paste("/mnt/atlas/my_home/binning_wf/Bacteria_V3/Blast_on_refs/Annot/",
#                          bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")
#annot = read.table(paste("/mnt/atlas/my_home/binning_wf/Bacteria_V3/Annot/",
#                         bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")

# subset = annot[annot[,3] %in% strsplit(refs_info[7,3],";")[[1]]
#                [2:length(strsplit(refs_info[7,3],";")[[1]])],]
# sum(subset[,2])
# 
# for (file in files) {
#     bin_name = system(sprintf("echo %s | egrep -o 'bin\\.[0-9]+'",file), intern=T)
#     annot = read.table(paste("/mnt/atlas/my_home/binning_wf/Bacteria_V3/Blast_on_refs/Annot/",
#                              bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")
    #annot = read.table(paste("/mnt/atlas/my_home/binning_wf/Bacteria_V3/Annot/",
    #                         bin_name, "_annotation.txt", sep=""), as.is=T, sep="\t")
    
#     annot_maj_strain = as.data.frame(t(apply(refs_info, 1, function(x,y) annot_maj_fun(x,y),
#                                              y=annot)),stringsAsFactors=F)
#     annot_other = annot[!(annot[,3] %in% strsplit(paste(refs_info[,3],collapse=";"),";")[[1]]),c(2,3)]
#     annot_maj_other = ddply(annot_other, .(V3), function(x) sum(x$V2))
#     colnames(annot_maj_other) = c("V1","V2")
#     annot_maj = rbind(annot_maj_strain, annot_maj_other)
#     annot_maj[,2] = as.numeric(annot_maj[,2])
#     annot_maj = annot_maj[order(annot_maj$V2, decreasing=T),]
#     print(bin_name)
#     print(sum(annot_maj[,2]))
# }


# plot(c(2,4,6,9), type ="l", col="red")
