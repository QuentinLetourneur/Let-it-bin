#!/usr/bin/env nextflow

/*workflow.onComplete {
    def subject = 'Benchmark Binning'
    def recipient = 'quentin.letourneur@pasteur.fr'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}*/

workflow.onComplete {
    //any worlflow property can be used here
    if ( workflow.success ) {
        println "Pipeline Complete"
    }
    println "Command line: $workflow.commandLine"
}

workflow.onError {
    println "Oops .. something went wrong"
}

params.help=false

def usage() {
    println("\nLet-it-bin has been developped to analyse metagenomic short paired-end reads data from multiple samples\nIt can be used on simulated data to benchmark binners performance or real data where binning results will only be evaluated by a reference free approach.\nThe pipeline take raw paired-end reads from multiple samples as primary input (the pairs names MUST finish with _{1,2}.fq/fastq). It comprise 4 major steps, reads preprocessing, assembly, binning and evaluation.\nThe pipeline can be started from the second or third step by adding arguments to the command line, provided that you have the needed inputs.")
    println("You have to select the binning softwares that will run in the following list :\n  binsanity, canopy, concoct, cocacola, maxbin, metabat, metabat2 and metagen\nYou just have to prefix the name of the wanted programms with '--' (Ex : --concoct).\nIf you want to use them all just use --all T\n\nWhen path are needed please give full path. For directories don't place a trailing '/' \n")
    println("    GENERAL ARGUMENTS :\n")
    println("  --reads [PATH] Directory containing unzipped paired reads files.\nOnly needed if you start from raw reads or reads from which contaminant have been removed")
    println("  --nb_samples [INT] Number of samples you are using")
    println("  --out [PATH] Directory were will be stored the results of the pipeline")
    println("  --sim_data [CHAR] Can be either T (Default) or F. Will change the execution of the pipeline depending on the analysed data (simulated or not).")
    println("  --cpus [INT] Number of cpus used for task with multiprocessing (Default 4)")
    println("  --min_contigs_length [INT] Minimum contigs length in base to be passed to binning programms (Default 1000)")
    println("  --nb_ref [INT] If you use simulated data specify the total number of different genomes present in the samples")
    println("  --dastool [CHAR] Can be either T (Default) or F. If you use multiple binning softwares you can use dastool to combine these results and try to extract best bins corresponding to the same microorganism.")
    println("  --tmp_checkm [PATH] Directory were will be stored CheckM temporary files. The path length souhldn't exeed 65 chars")
    println("  --help Print the help message")
    println("\n    READS PREPROCESSING :\n")
    println("  --contaminant [PATH] Path and prefix of the bowtie2 index files of contaminant sequences (HAVE TO be computed before lauching the pipeline)")
    println("  --minlength [INT] Minimum length for trimed contigs (Default 45)")
    println("  --alienseq [PATH] Fasta file containing adaptaters sequence for AlienTrimmer")
    println("  --cleaned_readsDir [PATH] Folder were will be stored reads that have been filtered to eliminate contaminant and trimmed (Default [out]/cleaned_reads)\nIf there are already fastq files in the folder it will take it as input for khmer and skip the cleaning step")
    println("  --filt_readsDir [PATH] Directory containing Khmer results.\nIF SPECIFIED the workflow will start at the assembly by taking as input the filtered fastq files in the directory.")
    println("\n    ASSEMBLY :\n")
    println("  --cpuassembly [INT] Number of cpus used for reads assembly (Default 10)")
    println("  --memassembly [INT] Quantity of RAM in Mb used for reads assembly. Default 160000 Mb with 2 retries. If the wanted memory is <= 160000, the allocated memory will grow according to the following formula : number of retry * [memmapping]\nElse no retry will be done")
    println("  --qos_assembly [STRING] If you run the pipeline on a cluster with SLURM you can specify the queue of submision (qos) to use : fast, normal (Default) or long (on TARS you can only use 5 cpus in long)")
    println("  --mode [STRING] Name of the assembler to be used. Can be either spades (Default), clc, megahit or ray")
    println("  --multi_assembly [CHAR] By default a co-assembly of the samples is done. If you want to change that behaviour and do an assembly by sample set this parameter to T (Default F). The generated contigs will then be pulled in one file and filtered to lessen the redundancy but eliminating it is hard so there migth still be redundancy in the filtered contigs")
    println("  --contigs [PATH] If the assembly has already been done you can specify the path to the fasta file containing contigs. If provided the assembly steps will be skipped")
    println("  --refs_info [PATH] If you use a simulated dataset and specify the --contigs option. Give the path to the sum_contigs_length_per_annotation.tsv file contained in [out]/assembly/")
    println("\n    CONTIGS ANNOTATION :\n")
    println("  --blast_db [PATH] If you use simulated data. Path to the BLAST database containing reference sequences (HAVE TO be computed before running the pipeline)")
    println("  --coverage [INT] Coverage threshold used to filter alignments of contigs on reference genomes or a public database (Default 90)")
    println("  --identity [INT] Identity threshold used to filter alignments of contigs on reference genomes or a public database (Default 95)")
    println("  --mismatch [INT] Number of mismatch allowed in the seed aligment of BLAST (Default 1)")
    println("  --evalue [INT] E-value used for BLAST (Default 10)")
    println("  --hit [INT] Maximum number of hits for each querry in the BLAST output (Default 10)")
    println("  --link_ref_id_species [PATH] For simulated dataset tab-delimited file containing contigs IDs of reference sequence and the species to which they belong. Used to identify the target in the BLAST of contigs against the references")
    println("  --contigs_annotation [PATH] If you use simulated data and have specified the --contigs option specify the path to the (bh_blast_contigs_on_refs.tsv)(to be updated) file in [out]/assembly/")
    println("\n    MAPPING :\n")
    println("  --cpumapping [INT] Number of cpus used for mapping reads on the assembly (Default 4)\n  One mapping per sample")
    println("  --memmapping [INT] Quantity of RAM in Mb for the mapping of each sample (Default 5000 )")
    println("  --bowtie2_indexDir [PATH] Directory were will be stored bowtie2 index of assembled sequences (Default [out]/bowtie2_index)")
    println("  --index_prefix [STRING] Prefix for the index files generated by bowtie2-build")
    println("  --bamDir [PATH] Directory were will be stored sorted BAM files generated by the mapping and if you use MetaGen indexed BAM files (Default [out]/mapping/bam). If there are sorted BAM files in the folder the mapping step will be skipped and they will be taken as input for the binning programs")
    println("  --count_matrix [PATH] If you use Canopy and have specified the --bamDir option. Provide the path to the count_matrix file contained by default in the [out]/mapping/comptage folder")
    println("\n    BINNING :\n")
    println("  --cpubinning [INT] Number of cpus for binning programs (Default same value as --cpus option)")
    println("  --nb_cluster [INT] Needed if you use Concoct. Rougth estimation of the expected number of species in the samples. This value will be a starting point for the program that will then refine it")
    println("  --bic_step [INT] MetaGen parameter corresponding to the step for the search of the optimal number of clusters (Default 5)")
    println("  --auto_method [INT] MetaGen parameter can be either 1 or 2. Recommended to be 2 for large datasets and is the default here")
    println("\n    BINNING EVALUATION :\n")
    println("  --min_bin_size [INT] Minimum size of bins in base (sum of the length of the sequences it contains) to be places in plots (Default 500000)")
    println("  --conta_threshold [FLOAT] Maximum contamination percent for good quality bins [0-1] (Default 0.1)")
    println("  --comp_threshold [FLOAT] Minimum completeness percent for good quality bins [0-1] (Default 0.6)")
    println("\n    EXAMPLES :\n")
    println("  For real data starting from raw reads or reads filtered from contaminant and co-assembly with Spades")
    println("  ~/nextflow -c ~/let-it-bin/nextflow_slurm_singularity_common.config run -w [directory to store temporary files] let-it-bin.nf --reads ~/data/reads --out ~/results --sim_data F --cpus 4\n   --tmp_checkm ~/tmp --metabat --canopy --maxbin \n   --index_prefix spades_contigs")
    println("\n  For simulated data starting from raw reads and co-assembly with megahit")
    println("  ~/nextflow -c ~/let-it-bin/nextflow_slurm_singularity_common.config run -w [directory to store temporary files] let-it-bin.nf --reads ~/data/reads --out ~/results --sim_data T --cpus 4\n  --nb_ref 50 --tmp_checkm ~/tmp\n   --metabat2 --cocacola --metagen --mode megahit\n   --blast_db ~/blast_db/refs_seq --link_ref_id_species ~/link_id_species.tsv\n   --index_prefix spades_contigs")
}

if( params.help ) {
    usage()
    exit(1)
}

params.reads = " "
params.nb_samples = ""
params.cpus = 4
params.memtrimming = '5g'
params.cpumapping = 4
params.memmapping = '5000'
params.cpuassembly = 10
params.memassembly = 160000
params.cpubinning = params.cpus
params.out = " "
params.binDir = "${params.out}/Binnings"
params.coverage = 90
params.identity = 95
params.mismatch = 1
params.alienseq = "${baseDir}/alienTrimmerPF8contaminants.fasta"
params.minlength = 45
params.cleaned_readsDir = "${params.out}/cleaned_reads"
params.mode = "spades"
params.contaminant = "" // path to the prefix of the index of bowtie2 for contaminants
params.bowtie2_indexDir = "${params.out}/bowtie2_index"
params.index_prefix = "" // prefix for the index of bowtie2 for analysed contigs
params.bamDir = "${params.out}/mapping/bam"
params.evalue = 10
params.gitaxidnucl = "/path/to/taxodb/gi_taxid_nucl.dmp"
params.names = "/path/to/taxodb/names.dmp"
params.nodes = "/path/to/taxodb/nodes.dmp"
params.hit = 10
params.nb_cluster = ""
params.metabat = " "
params.concoct = " "
params.cocacola = " "
params.maxbin = " "
params.canopy = " "
params.likelybin = " "
params.metagen = " "
params.binsanity = " "
params.metabat2 = " "
params.all = "F"
params.metabatDir = "${params.binDir}/Metabat"
params.concoctDir = "${params.binDir}/Concoct"
params.cocacolaDir = "${params.binDir}/Cocacola"
params.maxbinDir = "${params.binDir}/MaxBin"
params.canopyDir = "${params.binDir}/Canopy"
params.metagenDir = "${params.binDir}/MetaGen"
params.binsanityDir = "${params.binDir}/BinSanity-wf"
params.metabat2Dir = "${params.binDir}/Metabat2"
params.multi_assembly = "F"
params.contigs = ""
params.count_matrix = ""
// params.interleaved = ""
params.refs_info = ""
params.link_ref_id_species = "" // take the absolute path to the file containing the link between reference genome IDs and species name
params.blast_db = "/path/to/nt"
params.tmp_checkm = "${params.binDir}/tmp"
params.min_bin_size = 500000
params.conta_threshold = 0.1 // at max 10% contamination
params.comp_threshold = 0.6 // min 60% completeness
params.qos_assembly = "normal" //should be either fast, normal (default) or long
params.filt_readsDir = ""
params.plot_width = "750"
params.plot_heigth = "1000"
params.min_contig_length = 1000
params.sim_data = "F"
params.blast_ncbi = "F"
params.contigs_annotation = " "
params.bic_step = 5
params.auto_method = 2
params.nb_ref = ""
params.dastool = "T"
params.dastoolDir = "${params.binDir}/Dastool"

if( params.sim_data == "T" && params.nb_ref == "" ) {
    exit 1, "If you are using a simulated dataset you have to give the number of genomes in the dataset. It will be used to scale plot axis"
}

if( params.nb_samples == "" ) {
    exit 1, "You have to specify the number of samples you are using with the --nb_samples options"
}

myDir = file(params.out)
myDir.mkdirs()

cleanDir = file("${params.cleaned_readsDir}")
cleanDir.mkdirs()

if( params.contaminant != "" && params.reads != " ") {

    readChannel = Channel.fromFilePairs("${params.reads}/*{1,2}.{fastq,fq}")
                 .ifEmpty { exit 1, "Cannot find any reads file in ${params.reads}" }
                 //.subscribe { println it }

    // Map reads on the genome(s) of organisms that could have been sequenced along the wanted DNA.
    // Only reads that don't map on it will be kept
    process filtering {
        cpus params.cpus

        input:
        set pair_id, file reads from readChannel

        output:
        set pair_id, file("unmapped/*.1"), file("unmapped/*.2") into unmappedChannel

        shell:
        """
        export LC_ALL=C
        mkdir unmapped
        bowtie2 -q -N !{params.mismatch} -1 !{reads[0]} -2 !{reads[1]} -x ${params.contaminant} --un-conc unmapped/ -S /dev/null -p !{params.cpus} --very-sensitive-local
        """
    }

    // reads cleaning
    process trimming_rr {

        input:
        set pair_id, file(forward), file(reverse) from unmappedChannel

        output:
        set pair_id, file("*_paired_1.fastq"), file("*_paired_2.fastq") into trimChannel
        set pair_id, file("*_paired_1.fastq"), file("*_paired_2.fastq") into mappingChannel
        file("*paired_?.fastq") into outChannel mode flatten

        shell:
        """
        AlienTrimmer -if !{forward} -ir !{reverse} -of !{pair_id}_alien_1.fastq \
        -or !{pair_id}_alien_2.fastq -os un-conc-mate_sgl.fastq -c !{params.alienseq} \
        -l !{params.minlength} -p 80

        awk 'NR % 4 == 1 {gsub("/1\$","",\$0);gsub(" 1:.*","",\$0); print \$0}' !{pair_id}_alien_1.fastq > r1_headers.txt
        awk 'NR % 4 == 1 {gsub("/2\$","",\$0);gsub(" 2:.*","",\$0); print \$0}' !{pair_id}_alien_2.fastq > r2_headers.txt

        if [ `diff -q r1_headers r2_headers | wc -l` != "0" ];then
            #filter out single reads
            /home/bbmap/repair.sh in=!{pair_id}_alien_1.fastq in2=!{pair_id}_alien_2.fastq out=!{pair_id}_paired_1.fastq out2=!{pair_id}_paired_2.fastq outs=!{pair_id}_single.fastq -Xmx!{params.memtrimming}
            rm !{pair_id}_alien_1.fastq !{pair_id}_alien_2.fastq
        else
            mv !{pair_id}_alien_1.fastq !{pair_id}_paired_1.fastq
            mv !{pair_id}_alien_2.fastq !{pair_id}_paired_2.fastq
        fi
        """
    }
    outChannel.subscribe { it.copyTo(cleanDir) }

}
else if( file("${params.cleaned_readsDir}/*.f*q").size < params.nb_samples && params.reads != " " ) {

    readChannel = Channel.fromFilePairs("${params.reads}/*{1,2}.{fastq,fq}")
                 .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }

    // reads cleaning
    process trimming_fr {

        input:
        set pair_id, file(reads) from readChannel

        output:
        set pair_id, file("*_paired_1.fastq"), file("*_paired_2.fastq") into trimChannel
        set pair_id, file("*_paired_1.fastq"), file("*_paired_2.fastq") into mappingChannel
        file("*paired_?.fastq") into outChannel mode flatten

        shell:
        """
        AlienTrimmer -if !{reads[0]} -ir !{reads[1]} -of !{pair_id}_alien_1.fastq \
        -or !{pair_id}_alien_2.fastq -os un-conc-mate_sgl.fastq -c !{params.alienseq} \
        -l !{params.minlength} -p 80

        awk 'NR % 4 == 1 {gsub("/1\$","",\$0);gsub(" 1:.*","",\$0); print \$0}' !{pair_id}_alien_1.fastq > r1_headers.txt
        awk 'NR % 4 == 1 {gsub("/2\$","",\$0);gsub(" 2:.*","",\$0); print \$0}' !{pair_id}_alien_2.fastq > r2_headers.txt

        if [ `diff -q r1_headers r2_headers | wc -l` != "0" ];then
            #filter out single reads
            /home/bbmap/repair.sh in=!{pair_id}_alien_1.fastq in2=!{pair_id}_alien_2.fastq out=!{pair_id}_paired_1.fastq out2=!{pair_id}_paired_2.fastq outs=!{pair_id}_single.fastq -Xmx!{params.memtrimming}
            rm !{pair_id}_alien_1.fastq !{pair_id}_alien_2.fastq
        else
            mv !{pair_id}_alien_1.fastq !{pair_id}_paired_1.fastq
            mv !{pair_id}_alien_2.fastq !{pair_id}_paired_2.fastq
        fi
        """
    }
    outChannel.subscribe { it.copyTo(cleanDir) }
}
else {
    trimChannel = Channel.fromFilePairs("${params.cleaned_readsDir}/*{1,2}.{fastq,fq}").ifEmpty { exit 1, "No cleaned reads were found in ${params.cleaned_readsDir}"}.map{it.flatten()}
    mappingChannel = Channel.fromFilePairs("${params.cleaned_readsDir}/*{1,2}.{fastq,fq}").ifEmpty { exit 1, "No cleaned reads were found in ${params.cleaned_readsDir}"}.map{it.flatten()}
}

// Perform a co-assembly
if( params.contigs == "" && params.multi_assembly == "F" ) {

    if( params.filt_readsDir == "" ) {
        // filter reads redoudancy and even their aboudance
        process khmer_ca {
            cpus params.cpus

            input:
            set pair_id, file(fw), file(rv) from trimChannel

            output:
            file("*_filt_1.fastq") into R1_Channel
            file("*_filt_2.fastq") into R2_Channel

            script:
            """
            interleave-reads.py ${fw} ${rv} --output interleaved.pe
            normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct  interleaved.pe --output output.pe.keep
            rm interleaved.pe
            filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
            rm output.pe.keep
            extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
            rm output.pe.filter
            split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
            rm output.dn.pe
            mkdir -p ${params.out}/khmer_res
            cp *filt*fastq ${params.out}/khmer_res
            """
        }
    }
    else {
        R1_Channel = Channel.fromPath("${params.filt_readsDir}/*filt_1.{fastq,fq}").ifEmpty { exit 1, "No forward filt reads found in ${params.filt_readsDir}"}
        R2_Channel = Channel.fromPath("${params.filt_readsDir}/*filt_2.{fastq,fq}").ifEmpty { exit 1, "No reverse filt reads found in ${params.filt_readsDir}"}
    }

    // Pull all the forward/reverse reads together to make a co-assembly
    process merge {

        input:
        file R1 from R1_Channel.toList()
        file R2 from R2_Channel.toList()

        output:
        set val("mergesamples"), file("mergesamples_1.fastq"), file("mergesamples_2.fastq") into mergeChannel

        shell:
        """
        cat !{R1} > mergesamples_1.fastq
        cat !{R2} > mergesamples_2.fastq
        """
    }

    process assembly_co {

        if(params.mode != "ray") {
            cpus params.cpuassembly
            if(params.memassembly <= 160000) {
                errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
                maxRetries 2
                memory { params.memassembly * task.attempt }
            }
            else {
                memory params.memassembly
            }

            if(params.mode == "clc") {
                queue 'clcbio'
                clusterOptions="--qos=clcbio"
            }
            else {
                clusterOptions="--qos=${params.qos_assembly}"
            }
        }
        else {
            cpus 1
            clusterOptions="--qos=${params.qos_assembly}"
            exitReadTimeout = '604800sec'
        }

        input:
        set pair_id, file(forward), file(reverse) from mergeChannel

        output:
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_2
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_3
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_4
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_5
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_6
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_7
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_8
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_9
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_10
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_11
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_12
        file("assembly/*_{spades,clc,megahit,ray}.fasta") into assemblyChannel_13

        shell:
        """
        #!/bin/bash
        if [ !{params.mode} == "spades" ];then
            mkdir assembly
            spades.py --meta --only-assembler -1 !{forward} -2 !{reverse} -t !{params.cpuassembly} -o assembly/ -m 500
            #mv assembly/scaffolds.fasta assembly/!{pair_id}_spades.fasta
            mv assembly/contigs.fasta assembly/contigs_spades.fasta
        elif [ !{params.mode} ==  "clc" ];then
            mkdir assembly
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpuassembly}
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpuassembly} -n
            mv assembly/contigs.fasta assembly/!{pair_id}_clc.fasta
        elif [ !{params.mode} == "megahit" ];then
            megahit --mem-flag 0 -t !{params.cpuassembly} -1 !{forward} -2 !{reverse} -o assembly --out-prefix !{pair_id}_megahit
            mv assembly/!{pair_id}_megahit.contigs.fa assembly/!{pair_id}_megahit.fasta
            sed -i "s/ /_/g" assembly/!{pair_id}_megahit.fasta
        elif [ !{params.mode} == "ray" ];then
            sbatch -n !{params.cpuassembly} --mem-per-cpu 30000 --qos=!{params.qos_assembly} -J ray --wait -o ray.log --wrap="mpirun -n !{params.cpuassembly} Ray -k 31 -p !{forward} !{reverse} -o assembly"
            #mv assembly/Scaffolds.fasta assembly/!{pair_id}_ray.fasta
            mv assembly/Contigs.fasta assembly/!{pair_id}_ray.fasta
            sed -i "s/ /_/g" assembly/!{pair_id}_ray.fasta
        else
            echo "ERROR : the mode should be either spades, minia, megahit, ray or clc if you have it on your system"
            exit 1
        fi
        cp -r assembly !{params.out}
        """
    }

    if( params.sim_data == "T" && params.blast_db != "" && params.link_ref_id_species != "" ) {

        process contigs_annotation {
            cpus params.cpus

            input:
            file assembly from assemblyChannel_12

            output:
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_2
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_3
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_4
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_5
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_6
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_7
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_8
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_9
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_10

            shell:
            """
            infile=`echo !{assembly} | cut -f 1 -d "."`
            blastn -query !{assembly} -out "blast_"\$infile"_on_refs.txt" -outfmt \
                   "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
                   -db !{params.blast_db} \
                   -evalue !{params.evalue} -num_threads !{params.cpus} \
                   -max_target_seqs !{params.hit}

            ExtractAnnotation.py -f "blast_"\$infile"_on_refs.txt" -a !{params.link_ref_id_species} -o "annotation_"\$infile".txt" -nb 1 -r . -fc !{params.coverage} -fi !{params.identity}

            Rscript /usr/local/bin/sum_contigs_length_per_annotation.R "annotation_"\$infile".txt" sum_contigs_length_per_annotation.tsv

            cp "blast_"\$infile"_on_refs.txt" "annotation_"\$infile".txt" sum_contigs_length_per_annotation.tsv !{params.out}/assembly
            """
        }
    }
    else if ( params.sim_data == "F" && params.blast_db != "" ) {
        // annotate contigs on a public database (in progress)
        process contigs_annotation_rd {
            cpus params.cpus

            input:
            file assembly from assemblyChannel_12

            output:
            file("contigs_annotation.txt") into annotationChannel

            shell:
            """
            blastn -query !{assembly} -out blast_contigs_on_refs.txt -outfmt \
                       "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
                       -db !{params.blast_db} \
                       -evalue !{params.evalue} -num_threads !{params.cpus} \
                       -max_target_seqs !{params.hit}

            # Annot ncbi
            fasta_bn=`echo \$infile | cut -f 1,2 -d "." | cut -f 1 -d "_"`
            python get_taxonomy.py -i !{assembly} \
                    -o \$fasta_bn"_taxonomy.txt" -t !{params.gitaxidnucl} \
                    -n !{params.names} -d !{params.nodes}
            python ExtractNCBIDB.py \
                    -f !{assembly} -g \$fasta_bn"_taxonomy.txt" -fc !{params.coverage} \
                    -o \$fasta_bn"_annotation.txt" -nb 1 -fi

            # Get not annotated sequences
            if [ -f \$fasta_bn"_catalogue_annotation.txt" ]; then
                cat  \$fasta_bn"_annotation.txt" \$fasta_bn"_catalogue_annotation.txt"\
                 > annotated
                python extract_fasta.py -q annotated \
                    -t \$fasta_bn.fa -n -o \$fasta_bn"_not_annotated.fasta"
            else
                python extract_fasta.py \
                    -q \$fasta_bn"_annotation.txt" -t \$fasta_bn.fa -n \
                    -o \$fasta_bn"_not_annotated.fasta"
            fi
            """
        }
    }
    else if( params.sim_data == "F" && params.blast_db == "" ) {
        process decoy_2 {

            input:
            file assemby from assemblyChannel_12

            output:
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_2
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_3
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_4
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_5
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_6
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_7
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_8
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_9
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_10

            shell:
            """
            touch sum_contigs_length_per_annotation.tsv decoy_annotation.txt
            """
        }
    }
    else {
       exit 1, "If you are running the pipeline on simulated data, make sure you have specified --blast_db path to the BLAST database files and --lcs a file containing the association between sequence ids, the corresponding taxon annotation and a uniq number (see file given as exemple)"
    }
}

// perform a multi-assembly (to finish)
else if( params.contigs == "" && params.multi_assembly == "T" ) {

    if( params.filt_readsDir == "" ) {
        // filter reads singleton and reduce abundance of highly represented reads
        process khmer_ma {
            cpus params.cpus

            input:
            set pair_id, file(fw), file(rv) from trimChannel

            output:
            set pair_id, file("*_filt_1.fastq"), file("*_filt_2.fastq") into khmerChannel
            // file("interleaved.pe") into conv_fastqChannel

            script:
            """
            interleave-reads.py ${fw} ${rv} --output interleaved.pe
            normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct  interleaved.pe --output
            rm interleaved.pe
            filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
            rm output.pe.keep
            extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
            rm output.pe.filter
            split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
            rm output.dn.pe
            mkdir -p ${params.out}/khmer_res
            cp *filt*fastq ${params.out}/khmer_res
            """
        }
    }
    else {
        R1_Channel = Channel.fromPath("${params.filt_readsDir}/*filt_1.{fastq,fq}").ifEmpty { exit 1, "No forward filt reads found in ${params.filt_readsDir}"}
        R2_Channel = Channel.fromPath("${params.filt_readsDir}/*filt_2.{fastq,fq}").ifEmpty { exit 1, "No reverse filt reads found in ${params.filt_readsDir}"}
    }

    // make an assembly per sample
    process assembly_mult {
        publishDir "$myDir", mode: 'copy'

        if(params.mode != "ray") {
            cpus params.cpuassembly
            memory "70G"
        }
        else {
            cpus 1
        }

        if(params.mode == "clc") {
            queue 'clcbio'
            clusterOptions='--qos=${params.qos_assembly}'
        }
        else {
            clusterOptions="--qos=${params.qos_assembly}"
        }
        //~ else if(params.mode == "minia" || params.mode == "megahit" || params.mode == "ray"){
            //~ clusterOptions='--qos=long'
        //~ }

        input:
        //~ set pair_id, file(forward), file(reverse) from khmerChannel
        set pair_id, file(reads) from skiptrimChannel

        output:
        file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into contigsChannel
        //file("assembly/*_{spades,clc,minia}.fasta") into contigsChannel_2

        shell:
        """
        mkdir assembly
        if [ !{params.mode} == "spades" ];then
            spades.py --meta -1 !{reads[0]} -2 !{reads[1]} -t !{params.cpuassembly} -o assembly/
            mv assembly/scaffolds.fasta assembly/!{pair_id}_spades.fasta
        elif [ !{params.mode} ==  "clc" ];then
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{reads[0]} !{reads[1]} --cpus !{params.cpuassembly}
            mv assembly/contigs.fasta assembly/!{pair_id}_clc.fasta
        elif [ !{params.mode} == "minia" ];then
            /pasteur/projets/policy01/Matrix/metagenomics/entvirus/bin/gatb-minia-pipeline/gatb -1 !{forward} -2 !{reverse} -o assembly/!{pair_id}_minia.fasta
        elif [ !{params.mode} == "megahit" ];then
            megahit -m 0.5 --mem-flag 0 -t !{params.cpuassembly} -1 !{forward} -2 !{reverse} -o assembly/ --out-prefix !{pair_id}_megahit
            mv assembly/!{pair_id}_megahit.contigs.fa assembly/!{pair_id}_megahit.fasta
            sed -i "s/ /_/g" assembly/!{pair_id}_megahit.fasta
        elif [ !{params.mode} == "ray" ];then
            sbatch -n !{params.cpuassembly} --mem-per-cpu 10000 --qos=!{params.qos_assembly} -J ray --wait -o ray.log --wrap="mpirun -n !{params.cpuassembly} Ray -k 31 -p !{forward} !{reverse} -o assembly"
            mv assembly/Scaffolds.fasta assembly/!{pair_id}_ray.fasta
        else
            echo "ERROR : the mode should be either spades, clc, minia, megahit or ray"
            exit 1
        fi
        """
        //~ #spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpuassembly} -o assembly/
        //~ #clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpuassembly}
        //~ #interleave-reads.py !{forward} !{reverse} --output assembly/!{pair_id}.pe
            //~ #minia -in assembly/!{pair_id}.pe -out assembly/!{pair_id} -nb-cores !{params.cpus}
            //~ #!{baseDir}/gatb-minia-pipeline/gatb -1 !{forward} -2 !{reverse} -o assembly/!{pair_id}_gatb
            //~ #mv  assembly/!{pair_id}.contigs.fa assembly/!{pair_id}_minia.fasta
    }

    // Pull all the generated contigs together and filter their redundancy
    process cdhit {
        publishDir "$myDir/assembly", mode: 'copy'
        cpus params.cpus
        clusterOptions='--qos=normal -p common'

        input:
        file contigs from contigsChannel.toList()

        output:
        file("cata_contig_nr.fasta") into assemblyChannel
        file("cata_contig_nr.fasta") into assemblyChannel_2
        file("cata_contig_nr.fasta") into assemblyChannel_3
        file("cata_contig_nr.fasta") into assemblyChannel_4
        file("cata_contig_nr.fasta") into assemblyChannel_5
        file("cata_contig_nr.fasta") into assemblyChannel_6
        file("cata_contig_nr.fasta") into assemblyChannel_7
        file("cata_contig_nr.fasta") into assemblyChannel_8
        file("cata_contig_nr.fasta") into assemblyChannel_9
        file("cata_contig_nr.fasta") into assemblyChannel_10
        file("cata_contig_nr.fasta") into assemblyChannel_11
        file("cata_contig_nr.fasta") into assemblyChannel_12
        file("cata_contig_nr.fasta") into assemblyChannel_13

        shell:
        """
        cat !{contigs} > cata_contigs.fasta
        cd-hit-est -i cata_contigs.fasta -o cata_contig_nr.fasta -c 0.95 -T !{params.cpus} -aS 0.9 -d 0 -g 1 -M 0
        """
    }
}
else {
    if( file("${params.contigs}").isFile() ) {
        assemblyChannel = Channel.fromPath("${params.contigs}")
                         .ifEmpty { exit 1, "Cannot find the contigs file : ${params.contigs}" }
        assemblyChannel_2 = Channel.fromPath("${params.contigs}")
        assemblyChannel_3 = Channel.fromPath("${params.contigs}")
        assemblyChannel_4 = Channel.fromPath("${params.contigs}")
        assemblyChannel_5 = Channel.fromPath("${params.contigs}")
        assemblyChannel_6 = Channel.fromPath("${params.contigs}")
        assemblyChannel_7 = Channel.fromPath("${params.contigs}")
        assemblyChannel_8 = Channel.fromPath("${params.contigs}")
        assemblyChannel_9 = Channel.fromPath("${params.contigs}")
        assemblyChannel_10 = Channel.fromPath("${params.contigs}")
        assemblyChannel_11 = Channel.fromPath("${params.contigs}")
        assemblyChannel_12 = Channel.fromPath("${params.contigs}")
        assemblyChannel_13 = Channel.fromPath("${params.contigs}")
    }
    else {
        exit 1, "${params.contigs} is not a file"
    }

    if( params.sim_data == "T" ) {
        if( file("${params.refs_info}").isFile() && file("${params.contigs_annotation}").isFile() ) {
            annotationChannel = Channel.fromPath("${params.contigs_annotation}").ifEmpty { exit 1, "Cannot find the contigs annotation file : ${params.contigs_annotation}"}
            annot_by_ref_infosChannel = Channel.fromPath("${params.refs_info}").ifEmpty { exit 1, "Cannot find the references informations file : ${params.refs_info}"}
            annotationAndInfosByRefChannel = annot_by_ref_infosChannel.concat(annotationChannel).collect()

            annotationChannel_2 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_2 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_2 = annot_by_ref_infosChannel_2.concat(annotationChannel_2).collect()

            annotationChannel_3 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_3 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_3 = annot_by_ref_infosChannel_3.concat(annotationChannel_3).collect()

            annotationChannel_4 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_4 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_4 = annot_by_ref_infosChannel_4.concat(annotationChannel_4).collect()

            annotationChannel_5 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_5 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_5 = annot_by_ref_infosChannel_5.concat(annotationChannel_5).collect()

            annotationChannel_6 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_6 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_6 = annot_by_ref_infosChannel_6.concat(annotationChannel_6).collect()

            annotationChannel_7 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_7 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_7 = annot_by_ref_infosChannel_7.concat(annotationChannel_7).collect()

            annotationChannel_8 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_8 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_8 = annot_by_ref_infosChannel_8.concat(annotationChannel_8).collect()

            annotationChannel_9 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_9 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_9 = annot_by_ref_infosChannel_9.concat(annotationChannel_9).collect()

            annotationChannel_10 = Channel.fromPath("${params.contigs_annotation}")
            annot_by_ref_infosChannel_10 = Channel.fromPath("${params.refs_info}")
            annotationAndInfosByRefChannel_10 = annot_by_ref_infosChannel_10.concat(annotationChannel_10).collect()
        }
    }
    else {
        process decoy_1 {

            input:
            file assemby from assemblyChannel_12

            output:
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_2
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_3
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_4
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_5
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_6
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_7
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_8
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_9
            set file("sum_contigs_length_per_annotation.tsv"), file("*_annotation.txt") into annotationAndInfosByRefChannel_10

            shell:
            """
            touch sum_contigs_length_per_annotation.tsv decoy_annotation.txt
            """
        }
    }
}

bowt_refDir = file(params.bowtie2_indexDir)
bowt_refDir.mkdirs()

if( file("${params.bamDir}/*.bam").size == 0 &&  params.index_prefix != "" && (! file("${params.bowtie2_indexDir}/${params.index_prefix}.1.bt2").exists() && ! file("${params.bowtie2_indexDir}/${params.index_prefix}.rev.2.bt2").exists() )) {

    // generate bowtie2 index of the contigs
    process bowtie2_index {
        publishDir "$bowt_refDir", mode: 'copy'
        cpus params.cpus

        input:
        file assembly from assemblyChannel

        output:
        file("*.bt2") into indexChannel

        shell:
        """
        if [ `grep -v "^>" !{assembly} | wc -c` -lt \$((2**32)) ];then
            bowtie2-build !{assembly} !{params.index_prefix} --threads !{params.cpus}
        else
            bowtie2-build !{assembly} !{params.index_prefix} --threads !{params.cpus} --large-index
        fi
        """
    }
}
else {
    indexChannel = Channel.fromPath("${params.bowtie2_indexDir}/${params.index_prefix}*.bt2").ifEmpty { exit 1, "Can't find the index files : ${params.bowtie2_indexDir}/${params.index_prefix}*.bt2"}
}

if( file("${params.bamDir}/*.bam").size < params.nb_samples && params.index_prefix != "" ) {

    // map cleaned reads on the contigs
    process mapping_count {
        cpus params.cpumapping
        memory "${params.memmapping} MB"

        input:
        set pair_id, file(reads_1), file(reads_2) from mappingChannel
        file idx from indexChannel.first()

        output:
        file("bam/sorted*.bam") into sortedChannel mode flatten
        file("bam/sorted*.bam") into sortedChannel_2 mode flatten
        file("bam/sorted*.bam") into sortedChannel_3
        file("bam/sorted*.bam") into sortedChannel_4
        file("bam/sorted*.bam") into sortedChannel_5
        file("bam/sorted*.bam") into sortedChannel_6
        file("comptage/*.txt") into countsChannel

        shell:
        """
        #!/bin/bash

        export LC_ALL=C

        mkdir -p !{params.out}/mapping/{bam,comptage,out,unmapped}

        function timer()
        {
            if [[ \$# -eq 0 ]]; then
                echo \$(date '+%s')
            else
                local  stime=\$1
                etime=\$(date '+%s')
                if [[ -z '\$stime' ]]; then stime=\$etime; fi
                dt=\$((etime - stime))
                ds=\$((dt % 60))
                dm=\$(((dt / 60) % 60))
                dh=\$((dt / 3600))
                printf '%d:%02d:%02d' \$dh \$dm \$ds
            fi
        }

        mkdir unmapped out comptage bam sam

        echo "Mapping files: !{reads_1} & !{reads_2}"
        echo "Mapping with bowtie2 started at \$(date +"%T")"
        start_time=\$(timer)

        bowtie2 -p !{params.cpumapping} -x !{params.bowtie2_indexDir}/!{params.index_prefix} -q --local --sensitive-local -1 !{reads_1} -2 !{reads_2} --un-conc unmapped/!{pair_id} -S sam/!{pair_id}.sam > out/!{pair_id}.txt 2>&1

        echo "Elapsed time : \$(timer \$start_time)"

        start_time_count=\$(timer)
        echo "Count table creation started at \$(date +"%T")"

        counting.py sam/!{pair_id}.sam comptage/!{pair_id}.txt --best -m !{params.memmapping} -t !{params.cpumapping}

        echo "Elapsed time : \$(timer \$start_time_count)"
        echo "Total duration :  \$(timer \$start_time)"

        rm bam/filtered*
        find ./bam -maxdepth 1 -type f ! -iname sorted*.bam -delete
        cp bam/* !{params.out}/mapping/bam/
        cp out/* !{params.out}/mapping/out/
        cp comptage/* !{params.out}/mapping/comptage/
        cp unmapped/* !{params.out}/mapping/unmapped/
        """
    }
}
else if( file("${params.bamDir}/sorted*.bam").size == params.nb_samples.toInteger() && (params.concoct != " " || params.cocacola != " " || params.maxbin!= " " || params.metabat != " " || params.metagen != " " || params.binsanity != " " || params.metabat2 != " " || params.all == "T") ) {

    sortedChannel = Channel.fromPath("${params.bamDir}/sorted*.bam").ifEmpty { exit 1, "Cannot find any sorted BAM files in : ${params.bamDir}" }
    sortedChannel_2 = Channel.fromPath("${params.bamDir}/sorted*.bam")
    sortedChannel_3 = Channel.fromPath("${params.bamDir}/sorted*.bam")
    sortedChannel_4 = Channel.fromPath("${params.bamDir}/sorted*.bam")
    sortedChannel_5 = Channel.fromPath("${params.bamDir}/sorted*.bam")
    sortedChannel_6 = Channel.fromPath("${params.bamDir}/sorted*.bam")
    countsChannel = Channel.fromPath("${params.out}/mapping/comptage/*.txt").ifEmpty {exit 1, "Cannot find any count data in : ${params.out}/mapping/comptage" }
}
else {
    exit 1, "If you want to start at the binning steps you have to ensure that in the ${params.bamDir} folder are present all the sorted BAM files. Also make sure to have specified the wanted binning software to run or specified --all T"
}

if( (params.canopy != " " || params.all == "T") && params.count_matrix == "" ) {

    process count_matrix {

        input:
        file counts from countsChannel.toList()

        output:
        file("count_matrix.txt") into count_matrixChannel

        shell:
        """
        if [ -e !{params.out}/mapping/comptage/count_matrix.txt ];then
            rm !{params.out}/mapping/comptage/count_matrix.txt
        fi

        count_matrix.py -d !{params.out}/mapping/comptage -o count_matrix.txt
        summarise_mapping_PE.sh !{params.out}/mapping !{params.out}/mapping/stats_mapping.tsv

        cp count_matrix.txt !{params.out}/mapping/comptage
        """
    }
}
else if( (params.canopy != " " || params.all == "T") && params.count_matrix != "" ) {
    if( file("${params.count_matrix}").isFile() ) {
        count_matrixChannel = Channel.fromPath("${params.count_matrix}").ifEmpty { exit 1, "Can't find the count matrix file : ${params.count_matrix}" }
    }
    else {
        exit 1, "${params.count_matrix} is not a file"
    }
}
// else {
//     exit 1, "You have to specify the --count_matrix option if you want to start after the mapping (+ BAM sorting) step and use Canopy."
// }

if( params.metagen != " " || params.binsanity != " " || params.all == "T" ) {

    if( file("${params.bamDir}/*.bai").size == 0 || file("${params.bamDir}/*.bai").size < params.nb_samples ) {
        process index_bam {

            input:
            file sort_bam from sortedChannel

            output:
            file("*.bai") into indexedChannel
            file("*.bai") into indexedChannel_2
            file("*.bai") into indexedChannel_3

            shell:
            """
            samtools index !{sort_bam} !{sort_bam}.bai
            cp !{sort_bam}.bai !{params.out}/mapping/bam/
            """
        }
    }
    else if( file("${params.bamDir}/*.bai").size > params.nb_samples ) {
        exit 1, "  ERROR : There are more bam.bai files than than the given number of samples"
    }
    else {
        indexedChannel = Channel.fromPath("${params.bamDir}/*.bam.bai").ifEmpty { exit 1, "No indexed bam files found in ${params.bamDir}" }
        indexedChannel_2 = Channel.fromPath("${params.bamDir}/*.bam.bai").ifEmpty { exit 1, "No indexed bam files found in ${params.bamDir}" }
        indexedChannel_3 = Channel.fromPath("${params.bamDir}/*.bam.bai").ifEmpty { exit 1, "No indexed bam files found in ${params.bamDir}" }
    }
}

if( params.metagen != " " || params.all == "T" ) {

    process idxstat {

        input:
        file sort_bam from sortedChannel_6
        file idx_bam from indexedChannel.toList()

        output:
        file("*.stat") into statChannel

        shell:
        """
        samtools idxstats !{sort_bam} > !{sort_bam}.stat
        """
    }
}

if( params.concoct != " " || params.cocacola != " " || params.maxbin!= " " || params.all == "T" ) {

    // generate converage files from the alignment files
    process gen_cov_bed {

        input:
        file bam_sorted from sortedChannel_2

        output:
        file("*.gcbout") into gen_cov_bedChannel

        shell:
        """
        bn=`echo !{bam_sorted} | cut -f1 -d "."`
        genomeCoverageBed -ibam !{bam_sorted} > \$bn.gcbout
        """
    }

    // make abundance tables based on the coverage files
    process abun_and_link_profile {

        input:
        file bed_file from gen_cov_bedChannel.toList()
        file assembly from assemblyChannel_2

        output:
        file("covTable.tsv") into abundanceProfileChannel
        file("covTable.tsv") into abundanceProfileChannel_2
        file("covTable.tsv") into abundanceProfileChannel_3

        shell:
        """
        gen_input_table.py !{assembly} !{bed_file} --isbedfiles > covTable.tsv
        """
    }
}

if( params.metabat != " " || params.metabat2 != " " || params.all == "T" ) {

    // make abundance file for metabat
    process jgi_summa_depths {

        input:
        file bams from sortedChannel_3.toList()

        output:
        file("depth.txt") into abunChannel
        file("depth.txt") into abunChannel_2

        shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth depth.txt !{bams}
        """
    }
}

if( params.canopy != " " || params.all == "T" ) {
    file(params.canopyDir).mkdirs()

    process Canopy {
        cpus params.cpubinning

        input:
        file count_mat from count_matrixChannel
        file assembly from assemblyChannel_8
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_2

        output:
        set val("Canopy"), file("${refs_info}"), file("${contigs_annot}"), file("bin*") into canopyChannel
        // val("a") into canopyChannel_2
        val("Canopy") into canopyChannel_3
        file("*clustering.tsv") into canopyChannel_4

        shell:
        '''
        sed -r s/'(\t[0-9]+)\\.[0-9]+'/'\\1'/g !{count_mat} > tmp.txt
        awk -v min_size=!{params.min_contig_length} '\$2 >= min_size {print \$0}' tmp.txt > tmp2.txt
        tail -n +2 tmp2.txt | cut -f 1,3- -d $'\t' > count_matrix_gt1000.tsv
        #tail -n +2 tmp.txt | cut -f 1,3- -d $'\t' > count_matrix.tsv

        cc.bin -n !{params.cpubinning} -i count_matrix_gt1000.tsv -o cluster.out

        sed -i "s/>//" cluster.out
        nb_clust=`cut -f 1 cluster.out | uniq | wc -l`
        for clust in `seq $nb_clust`;do
            grep -P "CAG0*$clust\t" cluster.out | cut -f 2 > list.txt
            extract_fasta_from_list.py !{assembly} list.txt bin.$clust.fa
        done

        awk -F \$'\t' '{print \$2,\$1}' OFS=\$'\t' cluster.out > canopy_clustering.tsv
        sed -i 's/CAG/Canopy_bin./' canopy_clustering.tsv

        cp bin* !{params.canopyDir}
        '''
    }
}
else {
    canopyChannel = Channel.empty()
    // canopyChannel_2 = Channel.from("a")
    canopyChannel_3 = Channel.empty()
    canopyChannel_4 = Channel.from("")
}

if( params.maxbin != " " || params.all == "T" ) {
    file(params.maxbinDir).mkdirs()

    process maxbin {
        cpus params.cpubinning

        input:
        file assembly from assemblyChannel_7
        file abun from abundanceProfileChannel_3
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_3

        output:
        set val("MaxBin"), file("${refs_info}"), file("${contigs_annot}"), file("bin*") into maxbinChannel
        // val("a") into maxbinChannel_2
        val("MaxBin") into maxbinChannel_3
        file("*clustering.tsv") into maxbinChannel_4

        shell:
        """
        nb_col=`awk '{print NF; exit}' !{abun}`
        split_abun_file.sh !{abun} \$nb_col
        ls cov_mean* > list_abun_files.txt

        run_MaxBin.pl -contig !{assembly} -out out -abund_list list_abun_files.txt -thread !{params.cpubinning} -min_contig_length !{params.min_contig_length}

        for file in `ls out.*.fasta`;do
            name=`echo \$file | cut -f 1,2 -d "." | sed 's/out/bin/'`
            mv \$file \$name.fa
        done

        for file in `ls *.fa`;do
            bn=`basename \$file .fa`
            grep "^>" \$file | sed 's/>//' > tmp_\$bn
            sed -i "s/\$/\tMaxbin_\$bn/" tmp_\$bn
        done
        cat tmp* > maxbin_clustering.tsv
        rm tmp*

        cp bin* !{params.maxbinDir}
        """
    }
}
else {
    maxbinChannel = Channel.empty()
    // maxbinChannel_2 = Channel.from("a")
    maxbinChannel_3 = Channel.empty()
    maxbinChannel_4 = Channel.from("")
}

if( params.metabat != " " || params.all == "T" ) {
    file(params.metabatDir).mkdirs()

    process Metabat {
        cpus params.cpubinning

        input:
        file assembly from assemblyChannel_3
        file depth from abunChannel
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_4

        output:
        set val("Metabat"), file("${refs_info}"), file("${contigs_annot}"), file("bin*") into metabatChannel
        // val("a") into metabatChannel_2
        val("Metabat") into metabatChannel_3
        file("*clustering.tsv") into metabatChannel_4

        shell:
        """
        if [ "!{params.min_contig_length}" -le "1500" ];then
            metabat1 -i !{assembly} -a depth.txt -o bin -t !{params.cpubinning} --seed 10 -m 1500
        else
            metabat1 -i !{assembly} -a depth.txt -o bin -t !{params.cpubinning} -m !{params.min_contig_length} --seed 10
        fi

        for file in `ls *.fa`;do
            bn=`basename \$file .fa`
            grep "^>" \$file | sed 's/>//' > tmp_\$bn
            sed -i "s/\$/\tMetabat_\$bn/" tmp_\$bn
        done
        cat tmp* > metabat_clustering.tsv
        rm tmp*

        cp bin* !{params.metabatDir}
        """
    }
}
else {
    metabatChannel = Channel.empty()
    // metabatChannel_2 = Channel.from("a")
    metabatChannel_3 = Channel.empty()
    metabatChannel_4 = Channel.from("")
}

if( params.metabat2 != " " || params.all == "T" ) {
    file(params.metabat2Dir).mkdirs()

    process Metabat2 {
        cpus params.cpubinning

        input:
        file assembly from assemblyChannel_13
        file depth from abunChannel_2
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_5

        output:
        set val("Metabat2"), file("${refs_info}"), file("${contigs_annot}"), file("bin*") into metabat2Channel
        // val("a") into metabat2Channel_2
        val("Metabat2") into metabat2Channel_3
        file("*clustering.tsv") into metabat2Channel_4

        shell:
        """
        if [ "!{params.min_contig_length}" -le "1500" ];then
            metabat2 -i !{assembly} -a depth.txt -o bin -t !{params.cpubinning} --seed 10
        else
            metabat2 -i !{assembly} -a depth.txt -o bin -t !{params.cpubinning} -m !{params.min_contig_length} --seed 10
        fi

        for file in `ls *.fa`;do
            bn=`basename \$file .fa`
            grep "^>" \$file | sed 's/>//' > tmp_\$bn
            sed -i "s/\$/\tMetabat2_\$bn/" tmp_\$bn
        done
        cat tmp* > metabat2_clustering.tsv
        rm tmp*

        cp bin* !{params.metabat2Dir}
        """
    }
}
else {
    metabat2Channel = Channel.empty()
    // metabat2Channel_2 = Channel.from("a")
    metabat2Channel_3 = Channel.empty()
    metabat2Channel_4 = Channel.from("")
}

if( params.concoct != " " || params.all == "T" ) {
    file(params.concoctDir).mkdirs()

    process CONCOCT {

        input:
        file assembly from assemblyChannel_4
        file abun_prof from abundanceProfileChannel
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_6

        output:
        set val("Concoct"), file("${refs_info}"), file("${contigs_annot}"), file("bin*") into concoctChannel
        // val("a") into concoctChannel_2
        val("Concoct") into concoctChannel_3
        file("*clustering.tsv") into concoctChannel_4

        shell:
        """
        cut -f1,3- !{abun_prof} > covTableR.tsv
        concoct -c !{params.nb_cluster} --coverage_file covTableR.tsv --composition_file !{assembly} -b Concoct/ -l !{params.min_contig_length}

        nb_clust=`sort -t "," -rnk2,2 Concoct/clustering_gt1000.csv | head -n 1 | cut -f 2 -d ","`
        for clust in `seq 0 \$nb_clust`;do
            grep ",\$clust\$" Concoct/clustering_gt1000.csv | cut -f 1 -d "," > list.txt
            python /usr/local/bin/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        done

        sed 's/,/\tConcoct_bin./' Concoct/clustering_gt1000.csv > concoct_clustering.tsv

        cp bin* !{params.concoctDir}
        """
    }
}
else {
    concoctChannel = Channel.empty()
    // concoctChannel_2 = Channel.from("a")
    concoctChannel_3 = Channel.empty()
    concoctChannel_4 = Channel.from("")
}

if( params.cocacola != " " || params.all == "T" ) {
    file(params.cocacolaDir).mkdirs()

    process COCACOLA {
        cpus params.cpubinning

        input:
        file abun_table from abundanceProfileChannel_2
        file assembly from assemblyChannel_5
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_7

        output:
        set val("Cocacola"), file("${refs_info}"), file("${contigs_annot}"), file("bin*") into cocacolaChannel
        // val("a") into cocacolaChannel_2
        val("Cocacola") into cocacolaChannel_3
        file("*clustering.tsv") into cocacolaChannel_4

        shell:
        """
        echo "Generating coverage table ..."
        head -n1 !{abun_table} > covtable_gt"!{params.min_contig_length}".tsv
        tail -n +2 !{abun_table} | awk -v min_size=!{params.min_contig_length} '\$2 >= min_size {print \$0}' !{abun_table} >> covtable_gt"!{params.min_contig_length}".tsv
        tail -n +2 covtable_gt"!{params.min_contig_length}".tsv | awk '{print \$1}' > headers.txt
        echo "Done"

        echo "Generating kmer table ..."
        nb_contigs=`grep "^>" !{assembly} | wc -l`
        fasta_to_features.py !{assembly} \$nb_contigs 4 kmer_4.csv
        head -n1 kmer_4.csv > kmer_4_gt"!{params.min_contig_length}".csv
        grep -Fwf headers.txt kmer_4.csv >> kmer_4_gt"!{params.min_contig_length}".csv
        echo "Done"

        echo "Executing COCAOLA"
        cocacola.py --contig_file !{assembly} --abundance_profiles covtable_gt"!{params.min_contig_length}".tsv --composition_profiles kmer_4_gt"!{params.min_contig_length}".csv --out result_gt"!{params.min_contig_length}".csv --aux_dir /home/COCACOLA-python --threads !{params.cpubinning}
        echo "Done"

        num_clust=`sort -unk2,2 -t "," result_gt"!{params.min_contig_length}".csv | cut -f 2 -d ","`
        for clust in \$num_clust; do
            grep ",\$clust\$" result_gt"!{params.min_contig_length}".csv | cut -f 1 -d "," > list.txt
            python /usr/local/bin/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        done

        sed 's/,/\tCocacola_bin./' result_gt"!{params.min_contig_length}".csv > cocacola_clustering.tsv

        cp bin* !{params.cocacolaDir}
        """
    }
}
else {
    cocacolaChannel = Channel.empty()
    // cocacolaChannel_2 = Channel.from("a")
    cocacolaChannel_3 = Channel.empty()
    cocacolaChannel_4 = Channel.from("")
}

if( params.metagen != " " || params.all == "T" ) {

    file(params.metagenDir).mkdirs()

    process MetaGen {
        cpus params.cpubinning

        input:
        file assembly from assemblyChannel_10
        file count from statChannel.toList()
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_8

        output:
        set val("MetaGen"), file("${refs_info}"), file("${contigs_annot}"), file("bin*") into metagenChannel
        // val("a") into metagenChannel_2
        val("MetaGen") into metagenChannel_3
        file("*clustering.tsv") into metagenChannel_4

        shell:
        """
        mkdir map output contigs reads
        ln !{assembly} contigs/Contigs.fasta

        for f in !{count};do
            bn=`basename \$f .bam.stat`
            mkdir map/\$bn
            ln \$f map/\$bn/count.dat
        done

        /home/MetaGen-master/script/combine-counts.sh -p .
        /home/MetaGen-master/script/sum-reads.sh -p !{params.cleaned_readsDir} .
        Rscript /home/MetaGen-master/R/metagen.R -w . -n !{params.cpubinning} -l !{params.min_contig_length} -o !{params.auto_method} -s !{params.bic_step} -m /home/MetaGen-master

        sed -i -e 's/"//g' -e 's/>//' -e 's/ /,/' output/segs.txt
        num_clust=`sort -unk2,2 -t "," output/segs.txt | cut -f 2 -d ","`
        for clust in \$num_clust;do
            grep ",\$clust\$" output/segs.txt | cut -f 1 -d "," > list.txt
            extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        done

        sed 's/,/\tMetagen_bin./' output/segs.txt > metagen_clustering.tsv

        cp bin* !{params.metagenDir}
        """
    }
}
else {
    metagenChannel = Channel.empty()
    // metagenChannel_2 = Channel.from("a")
    metagenChannel_3 = Channel.empty()
    metagenChannel_4 = Channel.from("")
}

if( params.binsanity != " " || params.all == "T" ) {

    file(params.binsanityDir).mkdirs()

    process Binsanity_wf {
        cpus params.cpubinning

        input:
        file assembly from assemblyChannel_11
        file sorted_bam from sortedChannel_5.toList()
        file bam_index from indexedChannel_3.toList()
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_9

        output:
        set val("BinSanity-wf"), file("${refs_info}"), file("${contigs_annot}"), file("*.fa") into binsanityChannel
        // val("a") into binsanityChannel_2
        val("BinSanity-wf") into binsanityChannel_3
        file("*clustering.tsv") into binsanityChannel_4

        shell:
        """
        #source /local/gensoft2/exe/conda/3.19.0/conda/bin/activate envbinsanity
        get-ids -f "" -l !{assembly} -o ids.txt -x !{params.min_contig_length}
        Binsanity-profile -i !{assembly} -s . --ids ids.txt -c matrix -T !{params.cpubinning}
        Binsanity-wf -f "" -l !{assembly} -o . --threads !{params.cpubinning} --binPrefix "" -c matrix.cov

        for infile in `ls BinSanity-Final-bins/*.fna`; do
            nn=`echo \$infile | cut -f 2 -d "/" | sed -e 's/^_//' -e 's/\\.fna/\\.fa/'`
            mv \$infile ./\$nn
        done

        for infile in `ls *.fa`;do
            bn=`basename \$file .fa`
            grep "^>" \$file | sed 's/>//' > tmp_\$bn
            sed -i "s/\$/\tBinsanity-wf_\$bn/" tmp_\$bn
        done
        cat tmp* > binsanity-wf_clustering.tsv
        rm tmp*

        cp *.fa !{params.binsanityDir}
        """
    }
}
else {
    binsanityChannel = Channel.empty()
    // binsanityChannel_2 = Channel.from("a")
    binsanityChannel_3 = Channel.empty()
    binsanityChannel_4 = Channel.from("")
}

if( params.dastool == "T" ) {

    file(params.dastoolDir).mkdirs()

    process dastool {
        cpus params.cpus

        input:
        file inp1 from canopyChannel_4
        file inp2 from maxbinChannel_4
        file inp3 from metabatChannel_4
        file inp4 from metabat2Channel_4
        file inp5 from concoctChannel_4
        file inp6 from cocacolaChannel_4
        file inp7 from metagenChannel_4
        file inp8 from binsanityChannel_4
        file assembly from assemblyChannel_6
        set file(refs_info), file(contigs_annot) from annotationAndInfosByRefChannel_10

        output:
        set val("Dastool"), file("${refs_info}"), file("${contigs_annot}"), file("out_DASTool_bins/*.fa") into dastoolChannel
        val("Dastool") into dastoolChannel_2

        shell:
        """
        input=""
        for file in !{inp1} !{inp2} !{inp3} !{inp4} !{inp5} !{inp6} !{inp7} !{inp8};do
            if [[ "\$file" =~ "clustering" ]];then
                if [[ \$input == "" ]];then
                    input=\$file
                else
                    input="\$input,\$file"
                fi
            fi
        done

        DAS_Tool -i \$input -c !{assembly} -o ./out --score_threshold 0 --write_bins 1 -t !{params.cpus} --search_engine blast
        """
    }
}
else {
    dastoolChannel = Channel.empty()
    dastoolChannel_2 = Channel.empty()
}

checkmInputChannel = canopyChannel_3.mix(maxbinChannel_3, metabatChannel_3, metagenChannel_3, metabat2Channel_3, binsanityChannel_3, concoctChannel_3, cocacolaChannel_3, dastoolChannel_2)

// evaluate the quality of the bins generated by the chosen binning softwares
process checkm {
    cpus params.cpus

    input:
    val(soft) from checkmInputChannel

    output:
    file("chkm_res/tree_qa+qa.tsv") into checkmChannel
    val("plup") into resumeChannel_2

    shell:
    '''
    if ! mkdir chkm_res 2>/dev/null ; then
        rm -r chkm_res
        mkdir chkm_res
    fi

    mkdir -p !{params.tmp_checkm}/!{soft}

    checkm tree -t !{params.cpus} -x fa --tmpdir !{params.tmp_checkm}/!{soft} !{params.binDir}/!{soft} chkm_res
    checkm tree_qa -f chkm_res/tree_qa.tsv --tab_table --tmpdir !{params.tmp_checkm}/!{soft} chkm_res
    checkm lineage_set --tmpdir !{params.tmp_checkm}/!{soft} chkm_res chkm_res/lineage.ms
    checkm analyze -t !{params.cpus} --tmpdir !{params.tmp_checkm}/!{soft} -x fa chkm_res/lineage.ms !{params.binDir}/!{soft} chkm_res
    checkm qa -t !{params.cpus} -f chkm_res/qa_res.tsv --tab_table --tmpdir !{params.tmp_checkm}/!{soft} chkm_res/lineage.ms chkm_res

    echo -e 'Bin id\tTaxonomy\tMarker lineage\t# genomes\t# markers\tmarker sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain heterogeneity' > chkm_res/tree_qa+qa.tsv
    join -t $'\t' -1 1 -2 1 -o 1.1,1.4,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 <(tail -n +2 chkm_res/tree_qa.tsv | sort -k1,1) <(tail -n +2 chkm_res/qa_res.tsv | sort -k1,1) >> chkm_res/tree_qa+qa.tsv

    cp -r chkm_res !{params.binDir}/!{soft}
    '''
}

if( params.sim_data == "T" ) {

    AnnotInputChannel = canopyChannel.mix(maxbinChannel, metabatChannel, metagenChannel, metabat2Channel, binsanityChannel , concoctChannel, cocacolaChannel, dastoolChannel)

    process annotation_by_bin {

        input:
        set soft, file(refs_info), file(contigs_annot), file(bins) from AnnotInputChannel

        output:
        set val(soft), file("${refs_info}"), file("*_annotation.txt") into annot_binChannel

        script:
        """
        #!/bin/bash
        mkdir -p ${params.binDir}/${soft}/Annotation

        for infile in ${bins};do
            bin_name=`basename \$infile .fa`
            tax_count=`wc -l < \$infile`
            if [ "\$tax_count" -gt "0" ]; then
                grep "^>" \$infile | sed "s/>//" > headers.txt
                grep -wFf headers.txt ${contigs_annot} > \$bin_name"_annotation.txt"
            else
                touch \$bin_name"_annotation.txt"
            fi
        done

        cp *_annotation.txt ${params.binDir}/${soft}/Annotation
        """
    }
}
else if( params.sim_data == "T" && params.link_ref_id_species == "" ) {
    exit 1, "With simulated data, if you want to extract BLAST best hits and replace contigs name by the name of the species to which they belong specify the --lcs option with the path to the file annotation_num.tsv"
}

if( params.sim_data == "T" ) {
    process eval_complet_conta {

        input:
        set val(soft), file(refs_info), file(annotations) from annot_binChannel

        output:
        file("contamination.png") into evalChannel
        val("done") into resumeChannel

        shell:
        """
        Rscript /usr/local/bin/binning_stats.R !{params.binDir}/!{soft}/ !{params.binDir}/!{soft} !{params.binDir}/!{soft}/Annotation !{refs_info} !{params.min_bin_size} !{params.plot_width} !{params.plot_heigth}
        cp contamination.png completeness.png !{params.binDir}/!{soft}
        """
    }
}
else {
    resumeChannel = Channel.from("done")
}

process resume_res {

    input:
    file prec_rec_blast from resumeChannel.toList()
    file prec_rec_checkm from resumeChannel_2.toList()

    output:
    file("barplot*.png") into finishChannel

    shell:
    """
    if [ "!{params.sim_data}" == "T" ];then
        Rscript /usr/local/bin/nb_bin_per_threshold.R !{params.binDir} !{params.conta_threshold} !{params.comp_threshold} !{params.sim_data} !{params.nb_ref}
    else
        Rscript /usr/local/bin/nb_bin_per_threshold.R !{params.binDir} !{params.conta_threshold} !{params.comp_threshold} !{params.sim_data} "" c
    fi
    cp -f barplot*.png !{params.binDir}
    cp nb_bin_per_soft.csv !{params.binDir}
    """
}
