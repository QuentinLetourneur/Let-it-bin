for(package in c('ggplot2','plyr','reshape2')) {
    suppressPackageStartupMessages(do.call(library, list(package)))
}

args = commandArgs(TRUE)

# Binnings folder
input = args[1]

# contamination threshold
conta_thresh_blast = as.numeric(args[2])
comp_thresh_blast = as.numeric(args[3])

conta_thresh_checkm = conta_thresh_blast * 100
comp_thresh_checkm = comp_thresh_blast * 100

# specify if the summary of results is done based on BLAST ('b') or CheckM ('c')
# if no argument given, do them both

if(args[4] == "T") {
    sim_data = TRUE
} else {
    sim_data = FALSE
}

nb_ref = as.numeric(args[5])

choice = args[6]

if(is.na(args[7])) {
    output = getwd()
} else {
    output = args[7]
}

message("output : ", output)

# print(args)
# i=1
# for(arg in args) {
#     message("arg ",i," : ", arg)
#     i=i+1
# }

if ( !is.na(choice) & choice != "b" & choice != "c") {
    stop(message("If you want to only do the results summary based on BLAST give 'b' as the third argument or 'c' if you want to take CheckM results.\nIF no arguments are given summary will be done based on BLAST and CheckM"))
}

# count the number of bin per completeness interval
nb_bin_per_cat_score_fun = function (thresh, res, weight) {
    counts = c()
    score = c()
    for(i in seq(length(thresh) - 1)) {
        nb = length(which(res >= thresh[i] & res < thresh[i + 1]))
        counts = c(counts, nb)
        sc = nb * weight[i]
        score = c(score, sc)
    }
    output = data.frame(counts,score)
    return(output)
}

nb_bin_per_cat_fun = function (thresh, res) {
    counts = c()
    for(i in seq(length(thresh) - 1)) {
        nb = length(which(res >= thresh[i] & res < thresh[i + 1]))
        counts = c(counts, nb)
    }
    return(counts)
}

annot_sp_fct = function(x) {
    bst_annot = strsplit(x, ";")[[1]][length(strsplit(x, ";")[[1]])]
    if(! is.na(bst_annot) && strsplit(bst_annot, "_")[[1]][1] == "s") {
        return(bst_annot)
    } else {
        return(NA)
    }
}

nb_bin_per_soft = c()

weight_rec = c(1,2,2.5,3,3.5)
weight_prec = c(2.5,2,1.5,1,0.5,-0.5)

#BLAST
if (is.na(choice) | choice == '' | choice == "" | choice == "b") {
    
    eval_blast = system(sprintf("ls %s/*/prec_rec*.tsv", input), intern=T)
    
    thresh_rec_blast = c(0.6,0.7,0.8,0.9,0.95,1.00001)
    thresh_rec_combined_blast = c(0.6,0.7,0.8,0.9,0.95,1.00001,1.100001)
    thresh_prec_blast = c(0,0.05,0.10,0.15,0.2,0.3,1)
    
    tab_rec_uniq_blast = c()
    tab_prec_uniq_blast = c()
    tab_rec_combined_blast = c()
    tab_rec_score = c()
    tab_prec_score = c()
    
    # message("BLAST")
    
    # count the number of bin per contamination threshold per binning softwares
    for(file in eval_blast) {
    
        prec_rec_tab = read.table(file, sep="\t", as.is=T, header=T)
        soft = strsplit(file, "/")[[1]][length(strsplit(file, "/")[[1]])-1]
        nb_bin_per_soft = rbind(nb_bin_per_soft, c(soft, nrow(prec_rec_tab)))
        rec_clean = prec_rec_tab[which(prec_rec_tab$Contamination <= conta_thresh_blast),]
        prec_clean = prec_rec_tab[which(prec_rec_tab$Completeness >= comp_thresh_blast),]
        # message(soft)
        if(nrow(rec_clean) > 0) {
            rec_clean_ordered = rec_clean[order(rec_clean$Annot_maj_sp, -rec_clean$Completeness),]
            rec_clean_ordered_uniq = rec_clean_ordered[!duplicated(rec_clean_ordered$Annot_maj_sp),]
            rec_clean_combined = aggregate(Completeness ~ Annot_maj_sp, rec_clean_ordered, sum)
            
            count_score_rec_uniq = nb_bin_per_cat_score_fun(thresh_rec_blast, rec_clean_ordered_uniq[,"Completeness"], weight_rec)
            count_rec_uniq = count_score_rec_uniq$counts
            count_rec_combined = nb_bin_per_cat_fun(thresh_rec_combined_blast, rec_clean_combined[,"Completeness"])
            
            tab_rec_uniq_blast = rbind(tab_rec_uniq_blast,c(soft,count_rec_uniq))
            tab_rec_combined_blast = rbind(tab_rec_combined_blast,c(soft,count_rec_combined))
            tab_rec_score = rbind(tab_rec_score, c(soft, count_score_rec_uniq$score))
        } else {
            tab_rec_uniq_blast = rbind(tab_rec_uniq_blast,c(soft,rep(0, length(thresh_rec_blast) - 1)))
            tab_rec_combined_blast = rbind(tab_rec_combined_blast,c(soft,rep(0, length(thresh_rec_combined_blast) - 1)))
        }
        if(nrow(prec_clean) > 0) {
            prec_clean_ordered = prec_clean[order(prec_clean$Annot_maj_sp, -prec_clean$Completeness),]
            prec_clean_ordered_uniq = prec_clean_ordered[!duplicated(prec_clean_ordered$Annot_maj_sp),]
            
            count_score_prec_uniq = nb_bin_per_cat_score_fun(thresh_prec_blast, prec_clean_ordered_uniq[,"Contamination"], weight_prec)
            count_prec_uniq = count_score_prec_uniq$counts
            
            tab_prec_uniq_blast = rbind(tab_prec_uniq_blast, c(soft, count_prec_uniq))
            tab_prec_score = rbind(tab_prec_score, c(soft, count_score_prec_uniq$score))
        } else {
            tab_prec_uniq_blast = rbind(tab_prec_uniq_blast, c(soft,rep(0, length(thresh_prec_blast) - 1)))
        }
    }
    
    tab_rec_uniq_blast = as.data.frame(tab_rec_uniq_blast, stringsAsFactors = F)
    tab_rec_score = as.data.frame(tab_rec_score, stringsAsFactors = F)
    tab_rec_combined_blast = as.data.frame(tab_rec_combined_blast, stringsAsFactors = F)
    tab_prec_uniq_blast = as.data.frame(tab_prec_uniq_blast, stringsAsFactors = F)
    tab_prec_score = as.data.frame(tab_prec_score, stringsAsFactors = F)
    
    colnames(tab_rec_uniq_blast) = c('soft',as.character(thresh_rec_blast[1:length(thresh_rec_blast) - 1]))
    colnames(tab_rec_score) = c('soft',as.character(thresh_rec_blast[1:length(thresh_rec_blast) - 1]))
    colnames(tab_rec_combined_blast) = c('soft',as.character(thresh_rec_combined_blast[1:length(thresh_rec_combined_blast) - 1]))
    colnames(tab_prec_uniq_blast) = c('soft', as.character(thresh_prec_blast[1:length(thresh_prec_blast) - 1]))
    colnames(tab_prec_score) = c('soft', as.character(thresh_prec_blast[1:length(thresh_prec_blast) - 1]))
    
    tab_rec_uniq_blast = tab_rec_uniq_blast[which(tab_rec_uniq_blast$soft != "SCIMM"),]
    tab_rec_score = tab_rec_score[which(tab_rec_score$soft != "SCIMM"),]
    tab_rec_combined_blast = tab_rec_combined_blast[which(tab_rec_combined_blast$soft != "SCIMM"),]
    tab_prec_uniq_blast = tab_prec_uniq_blast[which(tab_prec_uniq_blast$soft != "SCIMM"),]
    tab_prec_score = tab_prec_score[which(tab_prec_score$soft != "SCIMM"),]
    
    if( 'MetaGen' %in% tab_rec_uniq_blast$soft) {
        tab_rec_score = tab_rec_score[which(tab_rec_score$soft != 'MetaGen'),]
        tab_prec_score = tab_prec_score[which(tab_prec_score$soft != 'MetaGen'),]
    }
    # else {
    #     tab_rec_score = tab_rec_uniq_blast
    #     tab_score_prec = tab_prec_uniq_blast
    # }
    
    # to obtain the 'score' for the binnings on a dataset
    score_rec = (mean(apply(as.matrix.data.frame(tab_rec_score[,-1]), 1, function(x) sum(as.numeric(x)))) / nb_ref) / max(weight_rec)
    score_prec = (mean(apply(as.matrix.data.frame(tab_prec_score[,-1]), 1, function(x) sum(as.numeric(x)))) / nb_ref) / max(weight_prec)
    
    barplot_rec_uniq_blast = melt(tab_rec_uniq_blast, id.vars='soft',
                             measures.vars=as.character(thres_rec_blast[1:length(thresh_rec_blast) - 1]), 
                             variable.name = 'Completeness', value.name = "number_of_bin")
    barplot_rec_combined_blast = melt(tab_rec_combined_blast, id.vars='soft', 
                                 measures.vars=as.character(thresh_rec_combined_blast[1:length(thresh_rec_combined_blast) - 1]), 
                                 variable.name = 'Completeness', value.name = "number_of_bin")
    barplot_prec_uniq_blast = melt(tab_prec_uniq_blast, id.vars='soft',
                              measures.name=as.character(thresh_prec_blast[1:length(thresh_prec_blast - 1)]),
                              variable.name = 'Contamination', value.name = "number_of_bins")
    
    barplot_rec_uniq_blast$number_of_bin = as.numeric(barplot_rec_uniq_blast$number_of_bin)
    barplot_rec_combined_blast$number_of_bin = as.numeric(barplot_rec_combined_blast$number_of_bin)
    barplot_prec_uniq_blast$number_of_bin = as.numeric(barplot_prec_uniq_blast$number_of_bin)
    
    png(filename=paste0(output, "/barplot_rec_uniq_BLAST.png"), width=468, height=571)
    p = ggplot(data=barplot_rec_uniq_blast, aes(x=soft,y=number_of_bin,fill=Completeness)) +
        geom_col(position = position_stack()) +
        scale_fill_brewer(palette = "Blues",name="Completeness (%)", 
                          labels=c("[60-70[","[70-80[","[80-90[","[90-95[","[95-100]")) +
        scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2", "Dastool")) +
        labs(y="Number of bins", x="Softwares") + #+ ylim(0,40)
        theme_minimal() + 
        theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
              legend.position = "bottom", legend.text=element_text(size=12), 
              axis.title = element_text(size=14), legend.title = element_text(size=14)) +
        annotate("text", label=paste0("Score : ", round(score_rec, 3)), x=4.5, y=nb_ref + 1, size = 6)
    print(p)
    dev.off()
    
    png(filename=paste0(output, "/barplot_rec_combined_BLAST.png"), width=468, height=571)
    p = ggplot(data=barplot_rec_combined_blast, aes(x=soft,y=number_of_bin,fill=Completeness)) +
        geom_col(position = position_stack()) +
        scale_fill_brewer(palette = "Blues",name="Completeness (%)", 
                          labels=c("[60-70[","[70-80[","[80-90[","[90-95[","[95-100]","]100-110]")) +
        scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2")) +
        labs(y="Number of bins", x="Softwares") + ylim(0, nb_ref) + 
        theme_minimal() +
        theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
              legend.position = "bottom", legend.text=element_text(size=12), 
              axis.title = element_text(size=14), legend.title = element_text(size=14))
    print(p)
    dev.off()
    
    png(filename=paste0(output, "/barplot_prec_uniq_BLAST.png"), width=468, height=571)
    p = ggplot(data=barplot_prec_uniq_blast, aes(x=soft,y=number_of_bin,fill=Contamination)) +
        geom_col(position = position_stack(reverse = T)) +
        scale_fill_brewer(palette = "Blues",name="Contamination (%)", 
                          labels=c("[0-5]","]5-10]","]10-15]","]15-20]","]20-30]","> 30"),
                          direction = -1) +
        scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2", "Dastool")) +
        labs(y="Number of bins", x="Softwares") + 
        theme_minimal() +
        theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
              legend.position = "bottom", legend.text=element_text(size=12), 
              axis.title = element_text(size=14), legend.title = element_text(size=14)) +
        annotate("text", label=paste0("Score : ", round(score_prec, 3)), x=4.5, y=nb_ref + 1, size = 6)
    print(p)
    dev.off()
}

#CheckM
if (is.na(choice) | choice == '' | choice == "" | choice == "c") {
    
    files_checkm = system(sprintf("ls %s/*/chkm_res/tree_qa+qa.tsv", input), intern=T)
    
    thresh_rec_chkm = c(60,70,80,90,95,100.00001)
    thresh_rec_combined_chkm = c(60,70,80,90,95,100.00001,110.00001)
    thresh_prec_chkm = c(0,5,10,15,20,30,500)
    
    tab_rec_uniq_checkm = c()
    tab_rec_combined_checkm = c()
    tab_prec_uniq_checkm = c()
    
    # message("Checkm")
    
    for(file in files_checkm) {
        soft = strsplit(file, "/")[[1]][length(strsplit(file, "/")[[1]])-2]
        prec_rec = read.table(file, sep="\t", as.is=T, header=T, comment.char = "")
        if(choice == 'c') {
            nb_bin_per_soft = rbind(nb_bin_per_soft, c(soft, nrow(prec_rec)))
        }
        rec_clean = prec_rec[which(prec_rec$Contamination <= conta_thresh_checkm),]
        prec_clean = prec_rec[which(prec_rec$Completeness >= comp_thresh_checkm),]
        if(nrow(rec_clean) > 0) {
            annot_sp = sapply(rec_clean$Taxonomy, annot_sp_fct)
            if(sum(is.na(annot_sp)) < length(annot_sp)) {
                rec_clean = cbind(rec_clean, annot_sp)
                rec_clean_ordered = rec_clean[order(rec_clean$annot_sp, -rec_clean$Completeness),]
                rec_clean_ordered_uniq = rec_clean_ordered[!duplicated(rec_clean_ordered$annot_sp, incomparables = NA),]
                rec_clean_combined = aggregate(Completeness ~ annot_sp, rec_clean_ordered, sum)
                rec_clean_combined = rbind(rec_clean_combined, rec_clean_ordered[which(is.na(rec_clean_ordered$annot_sp) == TRUE), c("annot_sp","Completeness")])
                
            } else {
                rec_clean_ordered_uniq = rec_clean
                rec_clean_combined = rec_clean
            }
            count_rec_uniq = nb_bin_per_cat_fun(thresh_rec_chkm, rec_clean_ordered_uniq[,c("Completeness")])
            count_rec_combined = nb_bin_per_cat_fun(thresh_rec_combined_chkm, rec_clean_combined[,c("Completeness")])
            tab_rec_uniq_checkm = rbind(tab_rec_uniq_checkm,c(soft,count_rec_uniq))
            tab_rec_combined_checkm = rbind(tab_rec_combined_checkm,c(soft,count_rec_combined))
        } else {
            tab_rec_uniq_checkm = rbind(tab_rec_uniq_checkm,c(soft,rep(0,6)))
            tab_rec_combined_checkm = rbind(tab_rec_combined_checkm,c(soft,rep(0,6)))
        }
        if(nrow(prec_clean) > 0) {
            annot_sp = sapply(prec_clean$Taxonomy, annot_sp_fct)
            if(sum(is.na(annot_sp)) < length(annot_sp)) {
                prec_clean = cbind(prec_clean, annot_sp)
                prec_clean_ordered = prec_clean[order(prec_clean$annot_sp, -prec_clean$Completeness),]
                prec_clean_ordered_uniq = prec_clean_ordered[!duplicated(prec_clean_ordered$annot_sp, incomparables = NA),]
            } else {
                prec_clean_ordered_uniq = prec_clean
            }
            
            count_score_prec_uniq = nb_bin_per_cat_score_fun(thresh_prec_chkm, prec_clean_ordered_uniq[,"Contamination"], weight_prec)
            count_prec_uniq = count_score_prec_uniq$counts
            
            tab_prec_uniq_checkm = rbind(tab_prec_uniq_checkm, c(soft, count_prec_uniq))
        } else {
            tab_prec_uniq_checkm = rbind(tab_prec_uniq_checkm, c(soft,rep(0, length(thresh_prec_checkm) - 1)))
        }
    }
    
    tab_rec_uniq_checkm = as.data.frame(tab_rec_uniq_checkm, stringsAsFactors = F)
    tab_rec_combined_checkm = as.data.frame(tab_rec_combined_checkm, stringsAsFactors = F)
    tab_prec_uniq_checkm = as.data.frame(tab_prec_uniq_checkm, stringsAsFactors = F)
    
    colnames(tab_rec_uniq_checkm) = c('soft',as.character(thresh_rec_chkm[1:length(thresh_rec_chkm) - 1]))
    colnames(tab_rec_combined_checkm) = c('soft',as.character(thresh_rec_combined_chkm[1:length(thresh_rec_combined_chkm) - 1]))
    colnames(tab_prec_uniq_checkm) = c('soft',as.character(thresh_prec_chkm[1:length(thresh_prec_chkm) - 1]))
    
    tab_rec_uniq_checkm = tab_rec_uniq_checkm[which(tab_rec_uniq_checkm$soft != "SCIMM"),]
    tab_rec_combined_checkm = tab_rec_combined_checkm[which(tab_rec_combined_checkm$soft != "SCIMM"),]
    tab_prec_uniq_checkm = tab_prec_uniq_checkm[which(tab_prec_uniq_checkm$soft != "SCIMM"),]
    
    if( ! 'MetaGen' %in% tab_rec_uniq_checkm$soft) {
        tab_rec_uniq_checkm = rbind(tab_rec_uniq_checkm, c('MetaGen', rep(0,6)))
        tab_rec_combined_checkm = rbind(tab_rec_combined_checkm, c('MetaGen', rep(0,6)))
        tab_prec_uniq_checkm = rbind(tab_prec_uniq_checkm, c('MetaGen', rep(0,6)))
    }
    
    barplot_rec_uniq_checkm = melt(tab_rec_uniq_checkm, id.vars='soft', 
                              measures.vars=as.character(thresh_rec_chkm[1:length(thresh_rec_chkm) - 1]), 
                              variable.name = 'Completeness', value.name = "number_of_bin")
    barplot_rec_combined_checkm = melt(tab_rec_combined_checkm, id.vars='soft', 
                                  measures.vars=as.character(thresh_rec_combined_chkm[1:length(thresh_rec_combined_chkm) - 1]), 
                                  variable.name = 'Completeness', value.name = "number_of_bin")
    barplot_prec_uniq_checkm = melt(tab_prec_uniq_checkm, id.vars='soft',
                               measures.vars=as.character(thresh_prec_chkm[1:length(thresh_prec_chkm) - 1]),
                               variable.name = 'Contamination', value.name = "number_of_bin")
    
    barplot_rec_uniq_checkm$number_of_bin = as.numeric(barplot_rec_uniq_checkm$number_of_bin)
    barplot_rec_combined_checkm$number_of_bin = as.numeric(barplot_rec_combined_checkm$number_of_bin)
    barplot_prec_uniq_checkm$number_of_bin = as.numeric(barplot_prec_uniq_checkm$number_of_bin)
    
    if(sim_data) {
        png(filename=paste0(output, "/barplot_rec_uniq_CheckM.png"), width=468, height=571)
        p = ggplot(data=barplot_rec_uniq_checkm, aes(x=soft,y=number_of_bin,fill=Completeness)) +
            geom_col(position = position_stack()) +
            scale_fill_brewer(palette = "Blues",name="Completeness (%)", 
                              labels=c("[60-70[","[70-80[","[80-90[","[90-95[","[95-100]")) +
            scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2", "Dastool")) +
            labs(y="Number of bins", x="Softwares") + ylim(0, nb_ref) +
            theme_minimal() +
            theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
                  legend.position = "bottom", legend.text=element_text(size=12), 
                  axis.title = element_text(size=14), legend.title = element_text(size=14))
        print(p)
        dev.off()
        
        png(filename=paste0(output, "/barplot_rec_combined_CheckM.png"), width=468, height=571)
        p = ggplot(data=barplot_rec_combined_checkm, aes(x=soft,y=number_of_bin,fill=Completeness)) +
            geom_col(position = position_stack()) +
            scale_fill_brewer(palette = "Blues",name="Completeness (%)", 
                              labels=c("[60-70[","[70-80[","[80-90[","[90-95[","[95-100]","]100-110]")) +
            scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2", "Dastool")) +
            labs(y="Number of bins", x="Softwares") + ylim(0, nb_ref) +
            theme_minimal() +
            theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
                  legend.position = "bottom", legend.text=element_text(size=12), 
                  axis.title = element_text(size=14), legend.title = element_text(size=14))
        print(p)
        dev.off()
        
        png(filename=paste0(output, "/barplot_prec_uniq_CheckM.png"), width=468, height=571)
        p = ggplot(data=barplot_prec_uniq_checkm, aes(x=soft,y=number_of_bin,fill=Contamination)) +
            geom_col(position = position_stack(reverse=T)) +
            scale_fill_brewer(palette = "Blues",name="Contamination (%)", 
                              labels=c("[0-5]","]5-10]","]10-15]","]15-20]","]20-30]","> 30"),
                              direction = -1) +
            scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2", "Dastool")) +
            labs(y="Number of bins", x="Softwares") + ylim(0, nb_ref) +
            theme_minimal() +
            theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
                  legend.position = "bottom", legend.text=element_text(size=12), 
                  axis.title = element_text(size=14), legend.title = element_text(size=14))
        print(p)
        dev.off()
    } else {
        png(filename=paste0(output, "/barplot_rec_uniq_CheckM.png"), width=468, height=571)
        p = ggplot(data=barplot_rec_uniq_checkm, aes(x=soft,y=number_of_bin,fill=Completeness)) +
            geom_col(position = position_stack()) +
            scale_fill_brewer(palette = "Blues",name="Completeness (%)", 
                              labels=c("[60-70[","[70-80[","[80-90[","[90-95[","[95-100]")) +
            scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2", "Dastool")) +
            labs(y="Number of bins", x="Softwares") +
            theme_minimal() +
            theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
                  legend.position = "bottom", legend.text=element_text(size=12), 
                  axis.title = element_text(size=14), legend.title = element_text(size=14))
        print(p)
        dev.off()
        
        png(filename=paste0(output, "/barplot_rec_combined_CheckM.png"), width=468, height=571)
        p = ggplot(data=barplot_rec_combined_checkm, aes(x=soft,y=number_of_bin,fill=Completeness)) +
            geom_col(position = position_stack()) +
            scale_fill_brewer(palette = "Blues",name="Completeness (%)", 
                              labels=c("[60-70[","[70-80[","[80-90[","[90-95[","[95-100]","]100-110]")) +
            scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2", "Dastool")) +
            labs(y="Number of bins", x="Softwares") +
            theme_minimal() +
            theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
                  legend.position = "bottom", legend.text=element_text(size=12), 
                  axis.title = element_text(size=14), legend.title = element_text(size=14))
        print(p)
        dev.off()
        
        png(filename=paste0(output, "/barplot_prec_uniq_CheckM.png"), width=468, height=571)
        p = ggplot(data=barplot_prec_uniq_checkm, aes(x=soft,y=number_of_bin,fill=Contamination)) +
            geom_col(position = position_stack(reverse=T)) +
            scale_fill_brewer(palette = "Blues",name="Contamination (%)", 
                              labels=c("[0-5]","]5-10]","]10-15]","]15-20]","]20-30]","> 30"),
                              direction = -1) +
            scale_x_discrete(limits = c("Canopy", "MetaGen", "BinSanity-wf", "Cocacola", "Concoct", "MaxBin", "Metabat", "Metabat2", "Dastool")) +
            labs(y="Number of bins", x="Softwares") +
            theme_minimal() +
            theme(axis.text.x= element_text(size="12", angle=45, hjust=1), axis.text.y= element_text(size="13"),
                  legend.position = "bottom", legend.text=element_text(size=12), 
                  axis.title = element_text(size=14), legend.title = element_text(size=14))
        print(p)
        dev.off()
    }
}

nb_bin_per_soft = as.data.frame(nb_bin_per_soft)
names(nb_bin_per_soft) = c('Softwares', 'number of bins')
write.csv(nb_bin_per_soft, file = paste(output, "nb_bin_per_soft.csv", sep="/"), quote = F, row.names = F)
# barplot_prec = melt(tab_prec, id.vars='soft', 
#                    measures.vars=as.character(thres_prec[1:length(thres_prec) - 1]), 
#                    variable.name = 'Contamination', value.name = "number_of_bin")
# 
# barplot_prec$number_of_bin = as.numeric(barplot_prec$number_of_bin)
# 
# ggplot(data=barplot_prec, aes(x=soft,y=number_of_bin,fill=Contamination)) +
#     geom_col(position = position_stack()) +
#     scale_fill_brewer(palette = "Blues")
    
# nb_per_thresh_prec_fun = function (thres, res) {
#     output = c()
#     for(i in seq(length(thres) - 1)) {
#         nb = length(which(res >= thres[i] & res < thres[i + 1]))
#         output=c(output, nb)
#     }
#     return(output)
# }

# tab_prec = c()