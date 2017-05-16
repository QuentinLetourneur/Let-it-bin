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
// to be completed
def usage() {
    println("benchmark_binning.nf --in <contigs_dir> --out <output_dir> --cpus <nb_cpus>")
}


if(params.help){
    usage()
    exit(1)
}

//baseDir = folder where is located the nextflow script
params.in="$baseDir/../../Simulated_data/Bacteria/"
readChannel = Channel.fromFilePairs("${params.in}/*_R{1,2}.fastq")
                     .ifEmpty { exit 1, "Cannot find any reads matching: ${params.in}" }
                     //.subscribe { println it }

/*contigsChannel = Channel
                .fromPath("${params.in}/*.fasta")
                .map {file -> tuple(file.baseName.replaceAll(".fasta",""),file)}*/

params.cpus = 4
params.out = "$baseDir/../../binning_wf"
params.coverage = 80
params.mismatch = 1
params.alienseq = "$baseDir/../../genomes_databases/alienTrimmerPF8contaminants.fasta"
params.minlength = 45
params.cleaned_reads = "${params.out}/cleaned_reads"
params.mode = "spades"
params.contaminant = "" // prefix for the index of bowtie2 for contaminants
params.bowt_index = "$baseDir/../../bowtie_ref"
params.index_prefix = "" // prefix for the index of bowtie2 for analysed contigs
params.bamDir = ""
params.mappingDir = "${params.out}/mapping" // were mapping results will be stored
params.metabatDir = "${params.out}/Metabat"
params.skiptrim = "F"
params.nt = "/pasteur/projets/policy01/BioIT/amine/catalogue/nt"
params.evalue = 10
params.rvdb = "/pasteur/projets/policy01/Biomics/metagenomics/catalogue/rVDBv10.2.fasta"
params.gitaxidnucl = "/local/databases/release/taxodb/gi_taxid_nucl.dmp"
params.annotDir = "${params.out}/Annot"
params.names = "/local/databases/release/taxodb/names.dmp"
params.nodes = "/local/databases/release/taxodb/nodes.dmp"
params.hit = 10
params.catalogue = "/pasteur/projets/policy01/Biomics/metagenomics/catalogue/metabat_bin.fa"
params.nb_cluster = ""
params.metabat = " "
params.concoct = " "
params.cocacola = " "
params.mycc = " "
params.groopm = " "
params.maxbin = " "
params.metacluster5 = " "
params.mbbc = " "
params.canopy = " "
params.abundancebin = " "
params.scimm = " "
params.all = "F"
params.concoctDir = "${params.out}/Concoct"
params.cocacolaDir = "${params.out}/Cocacola"
params.myccDir = "${params.out}/MyCC"
params.groopmDir = "${params.out}/GroopM"
params.maxbinDir = "${params.out}/MaxBin"
params.metacluster5Dir "${params.out}/Metacluster5"
params.mbbcDir = "${params.out}/MBBC"
params.canopyDir = "${params.out}/Canopy"
params.abundancebinDir = "${params.out}/AbundanceBin"
params.scimmDir = "${params.out}/SCIMM"
params.multi_assembly = "false"
params.contigs = ""
params.count_matrix = ""
params.scripts = "/pasteur/homes/qletourn/scripts"
//~ params.chkmDir = "${params.out}/chkm_res"
//~ params.vp1 = "$baseDir/databases/vp1_seq.fasta"
//~ params.ncbi = "$baseDir/databases/ncbi_viruses.fna"
//~ params.rvdb = "$baseDir/databases/rVDBv10.2.fasta"
//~ params.uniprot = "$baseDir/databases/uniprot_taxonomy.fasta"
//~ params.uniref = "$baseDir/databases/uniref_uniprot.fasta"
//~ params.viral = "$baseDir/databases/viral_catalogue_poltson.fna"


myDir = file(params.out)
myDir.mkdirs()

cleanDir = file("${params.cleaned_reads}")
cleanDir.mkdirs()

if( params.contaminant != "" ) {
    process filtering {
        //publishDir "$myDir", mode: 'copy'
        cpus params.cpus
        
        input:
        set pair_id, file reads from readChannel
        
        output:
        set pair_id, file("unmapped/*.1"), file("unmapped/*.2") into unmappedChannel
        
        shell:
        """
        mkdir unmapped
        bowtie2 -q -N !{params.mismatch} -1 !{reads[0]} -2 !{reads[1]} -x ${params.contaminant} --un-conc unmapped/ -S /dev/null -p !{params.cpus}
        """
    }
    
    
    process trimming {
        //publishDir "$cleanDir", mode: 'copy'
        cpus params.cpus
        
        input:
        set pair_id, file(forward), file(reverse) from unmappedChannel
        
        output:
        set pair_id, file("*_1.fastq"), file("*_2.fastq") into trimChannel
        file("*.fastq") into mappingChannel mode flatten
        file("*.fastq") into metaclustinputChannel mode flatten
        
        script:
        """
        AlienTrimmer -if ${forward} -ir ${reverse} -of ${pair_id}_1.fastq \
        -or ${pair_id}_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} \
        -l ${params.minlength}
        """
    }
    mappingChannel.subscribe { it.copyTo(cleanDir) }
    
    //~ process khmer {
        //~ cpus params.cpus
        //~ memory "50000"
        
        //~ input:
        //~ set pair_id, file(forward), file(reverse) from trimChannel
        
        //~ output:
        //~ set pair_id, file("*_1.fastq"), file("*_2.fastq") into khmerChannel
        
        //~ script:
        //~ """
        //~ interleave-reads.py ${forward} ${reverse} --output interleaved.pe
        //~ normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct interleaved.pe --output output.pe.keep
        //~ filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
        //~ extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
        //~ split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
        //~ """
    //~ }
}
else if( params.skiptrim == "F" ) {
    process trimming {
        //publishDir "$myDir", mode: 'copy'
        
        input:
        set pair_id, file(reads) from readChannel
        
        output:
        set pair_id, file("*_1.fastq"), file("*_2.fastq") into trimChannel
        file("*.fastq") into mappingChannel mode flatten
        file("*.fastq") into metaclustinputChannel mode flatten
        
        script:
        """
        AlienTrimmer -if ${reads[0]} -ir ${reads[1]} -of ${pair_id}_1.fastq \
        -or ${pair_id}_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} \
        -l ${params.minlength}
        """
    }
    mappingChannel.subscribe { it.copyTo(cleanDir) }
    // modif process khmer pr chaque condition
    //~ process khmer {
        //~ cpus params.cpus
        //~ memory "100000"
        
        //~ input:
        //~ set pair_id, file(forward), file(reverse) from trimChannel
        
        //~ output:
        //~ set pair_id, file("*_1.fastq"), file("*_2.fastq") into khmerChannel
        
        //~ script:
        //~ """
        //~ interleave-reads.py ${forward} ${reverse} --output interleaved.pe
        //~ normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct interleaved.pe --output output.pe.keep
        //~ filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
        //~ extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
        //~ split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
        //~ """
    //~ }
}
else if( params.contigs == "" && params.multi_assembly == "false" ) {
    mergeChannel = Channel.fromFilePairs("${params.cleaned_reads}/*_{1,2}.fastq").ifEmpty { exit 1, "No clean reads were found"}
    metaclustinputChannel = Channel.fromFilePairs("${params.cleaned_reads}/*_{1,2}.fastq").ifEmpty { exit 1, "No clean reads were found"}
    // modify so that there is on script to merge apply on the list of R1 files then on the R2
    // To render it usable independently of the number of samples
    // Add this box before each khmer process !
    process merge {
        
        input:
        set fA, fB, fC, fD, fE, fF, fG from mergeChannel.toList()
        
        output:
        set val("mergefstq"), file("*_1.fastq"), file("*_2.fastq") into trimChannel
        
        shell:
        """
        bash !{params.scripts}/merge_fastq.sh !{fA[1][0]} !{fA[1][1]} !{fB[1][0]} !{fB[1][1]} !{fC[1][0]} !{fC[1][1]} !{fD[1][0]} !{fD[1][1]} !{fE[1][0]} !{fE[1][1]} !{fF[1][0]} !{fF[1][1]} !{fG[1][0]} !{fG[1][1]}
        """
    }
}

if( params.contigs == "" ) {
    process khmer {
        cpus params.cpus
        memory "100000"
        clusterOptions='--qos=normal'
        
        input:
        set pair_id, file(fw), file(rv) from trimChannel
        
        output:
        set pair_id, file("*t_1.fastq"), file("*t_2.fastq") into khmerChannel
        file("output.dn.pe") into conv_fastqChannel
        
        script:
        """
        interleave-reads.py ${fw} ${rv} --output interleaved.pe
        normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct  interleaved.pe --output output.pe.keep
        filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
        extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
        split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
        """
    }
}

if( params.contigs == "" && params.multi_assembly == "false" ) {
    process assembly {
        // permet de copier dans myDir les fichiers de l'output
        publishDir "$myDir", mode: 'copy'
        
        cpus params.cpus
        if(params.mode == "clc"){
            clusterOptions='--qos=fast -C clcbio' 
        }
        input:
        set pair_id, file(forward), file(reverse) from khmerChannel
        
        output:
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_2
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_3
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_4
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_5
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_6
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_7
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_8
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_9
        
        shell:
        """
        #!/bin/bash
        mkdir assembly
        if [ !{params.mode} == "spades" ];then
            spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpus} -o assembly/
            mv assembly/scaffolds.fasta assembly/!{pair_id}_spades.fasta
        elif [ !{params.mode} ==  "clc" ];then
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpus}
            mv assembly/contigs.fasta assembly/!{pair_id}_clc.fasta
        else
            #interleave-reads.py !{forward} !{reverse} --output assembly/!{pair_id}.pe
            #minia -in assembly/!{pair_id}.pe -out assembly/!{pair_id} -nb-cores !{params.cpus}
            !{baseDir}/gatb-minia-pipeline/gatb -1 !{forward} -2 !{reverse} -o assembly/!{pair_id}_gatb
            #mv  assembly/!{pair_id}.contigs.fa assembly/!{pair_id}_minia.fasta
        fi
        """
    }
}
else if( params.contigs == "" && params.multi_assembly == "true" ) {
    
    process assembly {
        // permet de copier dans myDir les fichiers de l'output
        publishDir "$myDir", mode: 'copy'
        
        cpus params.cpus
        if(params.mode == "clc"){
            clusterOptions='--qos=fast -C clcbio' 
        }
        input:
        set pair_id, file(forward), file(reverse) from khmerChannel
        
        output:
        file("assembly/*_{spades,clc,minia}.fasta") into contigsChannel
        //file("assembly/*_{spades,clc,minia}.fasta") into contigsChannel_2
        
        shell:
        """
        mkdir assembly
        if [ !{params.mode} == "spades" ]
        then
            spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpus} -o assembly/
            mv assembly/scaffolds.fasta assembly/!{pair_id}_spades.fasta
        elif [ !{params.mode} ==  "clc" ]
        then
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpus}
            mv assembly/contigs.fasta assembly/!{pair_id}_clc.fasta
        else
            #interleave-reads.py !{forward} !{reverse} --output assembly/!{pair_id}.pe
            #minia -in assembly/!{pair_id}.pe -out assembly/!{pair_id} -nb-cores !{params.cpus}
            !{baseDir}/gatb-minia-pipeline/gatb -1 !{forward} -2 !{reverse} -o assembly/!{pair_id}_gatb
            #mv  assembly/!{pair_id}.contigs.fa assembly/!{pair_id}_minia.fasta
        fi
        """
    }
    
    
    process cdhit {
        publishDir "${myDir}/assembly", mode: 'copy'
        cpus params.cpus
        clusterOptions='--qos=normal'
        
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
        
        shell:
        """
        bash !{params.scripts}/merge_fasta.sh !{contigs}
        cd-hit-est -i cata_contigs.fasta -o cata_contig_nr.fasta -c 0.95 -T !{params.cpus} -aS 0.9 -d 0 -g 1 -M 0
        """
    }
}
else {
    assemblyChannel = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_2 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_3 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_4 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_5 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_6 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_7 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_8 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_9 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
}


bowt_refDir = file(params.bowt_index)
bowt_refDir.mkdirs()

//file(params.mappingDir).mkdirs()

if( params.bamDir == "" && params.index_prefix != "" && ! file("${params.bowt_index}/${params.index_prefix}.1.bt2").exists() ) {
    
    process index {
        publishDir "$bowt_refDir", mode: 'copy'
        cpus params.cpus
        
        input:
        file assembly from assemblyChannel
        
        output:
        file("*.bt2") into indexChannel
        
        shell:
        """
        bowtie2-build !{assembly} !{params.index_prefix} --threads !{params.cpus}
        """
    }
    
    process mapping_count {

        //cpus params.cpus
        
        input:
        file idx from indexChannel.first()
        
        output:
        file("mapping/bam/*.bam") into bamChannel mode flatten
        file("mapping/comptage/count_matrix.txt") into countChannel
        
        shell:
        """
        mbma.py mapping -i !{cleanDir} -o mapping -db !{params.bowt_index}/!{params.index_prefix} -t 6 -q fast --bowtie2 --shared -e quentin.letourneur@pasteur.fr
        bash !{params.scripts}/summarise_mapping_PE.sh mapping mapping/stats_mapping.tsv
        rm -r mapping/sam
        cp -r mapping/ !{params.out}/
        """
    }
    
    
    process sort_bam {
    
        cpus params.cpus
        
        input:
        file bam from bamChannel
        
        output:
        file("*.sorted") into sortedChannel
        file("*.sorted") into sortedChannel_2
        file("*.sorted") into sortedChannel_3
        file("*.sorted") into sortedChannel_4
        
        shell:
        """
        samtools sort -o !{bam}.sorted -@ !{params.cpus} -O bam !{bam}
        """
    }
}
else if (params.bamDir == "" && params.index_prefix != "") {
    
    process mapping_count {

        //cpus params.cpus
        
        input:
        file assembly from assemblyChannel
        
        output:
        file("mapping/bam/*.bam") into bamChannel mode flatten
        file("mapping/comptage/count_matrix.txt") into countChannel
        
        shell:
        """
        mbma.py mapping -i !{cleanDir} -o mapping -db !{params.bowt_index}/!{params.index_prefix} -t 6 -q fast --bowtie2 --shared -e quentin.letourneur@pasteur.fr
        bash !{params.scripts}/summarise_mapping_PE.sh mapping mapping/stats_mapping.tsv
        rm -r mapping/sam
        cp -r mapping/ !{params.out}/
        """
    }
    
    
    process sort_bam {
    
        cpus params.cpus
        
        input:
        file bam from bamChannel
        
        output:
        file("*.sorted") into sortedChannel
        file("*.sorted") into sortedChannel_2
        file("*.sorted") into sortedChannel_3
        file("*.sorted") into sortedChannel_4
        
        shell:
        """
        samtools sort -o !{bam}.sorted -@ !{params.cpus} -O bam !{bam}
        """
    }
}
else if ( params.bamDir != "" && params.count_matrix != "" ) {
    
    bamChannel = Channel.fromPath("${params.bamDir}/*.bam").ifEmpty { exit 1, "Cannot find any bams in the directory : ${params.in}" }
    
    process sort_bam {
        
        cpus params.cpus
        
        input:
        file bam from bamChannel
        
        output:
        file("*.sorted") into sortedChannel
        file("*.sorted") into sortedChannel_2
        file("*.sorted") into sortedChannel_3
        file("*.sorted") into sortedChannel_4
        
        shell:
        """
        samtools sort -o !{bam}.sorted -@ !{params.cpus} -O bam !{bam}
        """
    }
}
else {
    exit 1, "If no prefix for bowtie index is given you have to give the path to the bams directory and the count_matrix.txt file"
}

if ( params.concoct != " " || params.cocacola != " " || params.groopm != " " ||  params.all == "T" ) {
    
    process index_bam {
        
        //~ cpus params.cpus
        
        input:
        file sort_bam from sortedChannel_2
        
        output:
        file("*.bai") into indexedChannel
        file("*.bai") into indexedChannel_2
        
        shell:
        """
        samtools index !{sort_bam} !{sort_bam}.bai
        """
    }
    
    process abun_and_link_profile {
    
        //~ cpus params.cpus
        
        input:
        file bam_sorted from sortedChannel_3.toList()
        file bam_idx from indexedChannel.toList()
        file assembly from assemblyChannel_2
        
        output:
        file("covTable.tsv") into abundanceProfileChannel
        file("covTable.tsv") into abundanceProfileChannel_2
        file("covTable.tsv") into abundanceProfileChannel_3
        file("linkTable.tsv") into linkageChannel
        file("linkTable.tsv") into linkageChannel_2
        
        shell:
        """
        python /pasteur/homes/qletourn/tools/CONCOCT-0.4.0/scripts/gen_input_table.py !{assembly} !{bam_sorted} > covTable.tsv
        python /pasteur/homes/qletourn/tools/CONCOCT-0.4.0/scripts/bam_to_linkage.py -m !{params.cpus} --regionlength 500 --fullsearch !{assembly} !{bam_sorted} > linkTable.tsv
        """
    }
}

//~ process fastq_to_fasta {
    
    //~ input:
    //~ file fastq 
    
//~ }

if( params.metabat != " " || params.mycc != " " || params.all == "T" ) {
    
    process jgi_summa_depths {
        
        input:
        file bams from sortedChannel.toList()
        
        output:
        file("depth.txt") into abunChannel
        file("depth.txt") into abunChannel_2
        
        shell:
        """
        /pasteur/homes/qletourn/tools/metabat/jgi_summarize_bam_contig_depths --outputDepth depth.txt !{bams}
        """
    }
}

//~ if( params.scimm != " " || params.all == "T" ) {
    //~ file(params.scimmDir).mkdirs()
    
        //~ process SCIMM {
        
        //~ input:
        
        //~ output:
        
        //~ shell:
        //~ """
        
        //~ """
    //~ }
//~ }

//~ if( params.abundancebin != " " || params.all == "T" ) {
    //~ file(params.abundancebinDir).mkdirs()
    
        //~ process AbundanceBin {
        
        //~ input:
        
        //~ output:
        
        //~ shell:
        //~ """
        
        //~ """
    //~ }
//~ }

if( params.canopy != " " || params.all == "T" ) {
    file(params.canopyDir).mkdirs()
    
        process Canopy {
        cpus params.cpus
        memory "100G"
        
        input:
        file count_mat from countChannel
        
        output:
        output file("*") into canopyChannel
        
        shell:
        """
        sed -r "s/(\t[0-9]+)\.[0-9]+/\1/g" !{count_mat} > tmp.txt
        tail -n +2 tmp.txt | cut -f 1,3- -d $'\t' > count_matrix.tsv
        cc.bin -n !{params.cpus} -i count_matrix.tsv -o bin -c //...
        //...
        """
    }
}

//~ if( params.mbbc != " " || params.all == "T" ) {
    //~ file(params.mbbcDir).mkdirs()
    
        //~ process MBBC {
        //~ memory "60000"
        
        //~ input:
        
        
        //~ output:
        
        
        //~ shell:
        //~ """
        //~ java -jar -Xmx50g -i 
        //~ """
    //~ }
//~ }

if( params.metacluster5 != " " || params.all == "T" ) {
    file(params.metacluster5Dir).mkdirs()
    
    process Metacluster5 {
        
        input:
        file assembly from assemblyChannel_9
        file reads from metaclustinputChannel.toList()
        
        output:
        file("*") into metacluster5Channel
        
        shell:
        """
        Metacluster5_1 !{assembly} !{reads} 
        Metacluster5_2 !{assembly} !{reads}
        """
    }
}

if( params.maxbin != " " || params.all == "T" ) {
    file(params.maxbinDir).mkdirs()
    
    process maxbin {
        cpus params.cpus
        
        input:
        file assembly from assemblyChannel_8
        file abun from abundanceProfileChannel_3
        
        output:
        file("bin*") into maxbinChannel
        
        shell:
        """
        nb_col=`awk '{print NF; exit}' !{abun}`
        bash !{params.scripts}/split_abun_file.sh !{abun} \$nb_col
        ls cov_mean* > list_abun_files.txt
        perl run_MaxBin.pl -contigs !{assembly} -out bin -abund_list list_abun_files.txt -thread !{params.cpus}
        for file in `ls bin.0*`;do
            basename=`cut -f 1,2 -d "." $file`
            mv $file $basename.fa
        done
        """
    }
}

if( params.groopm != " " || params.all == "T" ) {
    file(params.groopmDir).mkdirs()
    
    process GroopM {
        cpus params.cpus
        
        input:
        file assembly from assemblyChannel_7
        file bams from sortedChannel_4.toList()
        file bams_idx from indexedChannel_2
        
        output:
        file("*") into groopmChannel
        
        shell:
        """
        groopm parse -t !{params.cpus} -c 1000 database.gm !{assembly} !{bams}
        groopm core database.gm
        groopm recruit -c 1000
        groopm extract -t !{params.cpus} database.gm !{assembly} -o GroopM
        """
    }
}

if( params.metabat != " " || params.all == "T" ) {
    file(params.metabatDir).mkdirs()
    
    process Metabat {
        publishDir "${params.binDir}", mode: 'copy'
        cpus params.cpus
        
        input:
        file assembly from assemblyChannel_3
        file depth from abunChannel
        
        output:
        file("bin*") into metabatChannel
        file("bin*") into metabatChannel_2 mode flatten
        
        shell:
        """
        /pasteur/homes/qletourn/tools/metabat/metabat -i !{assembly} -a depth.txt -o bin -t !{params.cpus} --minSamples 5
        """
    }
}

if( params.concoct != " " || params.all == "T" ) {
    file(params.concoctDir).mkdirs()
    
    process CONCOCT {
        
        //~ cpus params.cpus
        
        input:
        file assembly from assemblyChannel_4
        file abun_prof from abundanceProfileChannel
        file link from linkageChannel
        
        output:
        file("Concoct/*") into concoctChannel
        
        shell:
        """
        cut -f1,3- !{abun_prof} > covTableR.tsv
        concoct -c !{params.nb_cluster} --coverage_file covTableR.tsv --composition_file !{assembly} -b Concoct/
        
        nb_clust=`sort -t "," -rnk2,2 Concoct/clustering_gt1000.csv | head -n 1 | cut -f 2 -d ","`
        for clust in `seq \$nb_clust`;do
            grep ",\$clust" Concoct/clustering_gt1000.csv | cut -f 1 -d "," > list.txt
            python !{params.scripts}/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        done
        """
    }
}

if( params.cocacola != " " || params.all == "T" ) {
    file(params.cocacolaDir).mkdirs()
    
    process COCACOLA {
        
        //~ cpus params.cpus
        
        input:
        file abun_table from abundanceProfileChannel_2
        file link from linkageChannel_2
        file assembly from assemblyChannel_5
        
        output:
        file("*") into cocacolaChannel
        
        shell:
        """
        nb_contigs=`grep "^>" !{assembly} | wc -l`
        python /pasteur/homes/qletourn/tools/CONCOCT-0.4.0/scripts/fasta_to_features.py !{assembly} \$nb_contigs 4 kmer_4.csv
        python /pasteur/homes/qletourn/tools/COCACOLA-python/cocacola.py --contig_file !{assembly} --abundance_profiles !{abun_table} --composition_profiles kmer_4.csv --out result.csv
        """
    }
}

if( params.mycc != " " || params.all == "T" ) {
    file(params.myccDir).mkdirs()
    
    process MyCC {
        input:
        file assembly from assemblyChannel_6
        file coverage from abunChannel_2
        
        output:
        
        shell:
        """
        cut -f1,4- !{coverage} > abun.txt
        python /pasteur/homes/qletourn/tools/MyCC/MyCC.py !{assembly} -a abun.txt -meta
        """
    }
}


//~ process checkm {
    //~ cpus params.cpus
    //~ memory "100G"
    
    //~ input:
    //~ file(bins) from binChannel.first()
    
    //~ output:
    //~ file("chkm_res/tree_qa+qa.tsv") into annotChannel
    
    //~ shell:
    //~ '''
    //~ if ! mkdir chkm_res 2>/dev/null ; then
        //~ rm -r chkm_res
        //~ mkdir chkm_res
    //~ fi
    
    //~ checkm tree -t !{params.cpus} -x fa --tmpdir /pasteur/homes/qletourn/tmp_chkm !{params.binDir} chkm_res
    //~ checkm tree_qa -f chkm_res/tree_qa.tsv --tab_table --tmpdir /pasteur/homes/qletourn/tmp_chkm chkm_res
    //~ checkm lineage_set --tmpdir /pasteur/homes/qletourn/tmp_chkm chkm_res chkm_res/lineage.ms
    //~ checkm analyze -t !{params.cpus} --tmpdir /pasteur/homes/qletourn/tmp_chkm -x fa chkm_res/lineage.ms !{params.binDir} chkm_res
    //~ checkm qa -t !{params.cpus} -f chkm_res/qa_res.tsv --tab_table --tmpdir /pasteur/homes/qletourn/tmp_chkm chkm_res/lineage.ms chkm_res
    
    //~ echo -e 'Bin id\tTaxonomy\tMarker lineage\t# genomes\t# markers\tmarker sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain heterogeneity' > chkm_res/tree_qa+qa.tsv
    //~ join -t $'\t' -1 1 -2 1 -o 1.1,1.4,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 <(tail -n +2 chkm_res/tree_qa.tsv | sort -k1,1) <(tail -n +2 chkm_res/qa_res.tsv | sort -k1,1) >> chkm_res/tree_qa+qa.tsv
    
    //~ cp -r chkm_res !{params.out}
    //~ '''
//~ }


//~ process dotplot {
    //~ publishDir "${params.out}", mode: 'copy'
    
    //~ input:
    //~ file annot from annotChannel
    
    //~ output:
    //~ file ("*.png") into dpChannel
    
    //~ shell:
    //~ """
    //~ bash !{params.scripts}/dotplot.sh !{annot} !{params.binDir}
    //~ """
//~ }

//~ file(params.annotDir).mkdirs()

//~ process blast {
    //~ publishDir "$params.annotDir", mode: 'copy'
    //~ cpus params.cpus
    //~ memory "20G"
    
    //~ input:
    //~ file(fasta) from binChannel_2
    
    //~ output:
    //~ set file(fasta), file("*_nt.txt") into blastChannel //, file("*_catalogue.txt")
    
    //~ script:
    //~ fasta_bn = fasta.baseName
    //~ """
    //~ #!/bin/bash
    //~ # local DB
    //~ blastn -query ${fasta} -out ${fasta_bn}_nt.txt -outfmt \
           //~ "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
           //~ -db /pasteur/homes/qletourn/Blastdb/Bact_Anita \
           //~ -evalue ${params.evalue} -num_threads ${params.cpus} \
           //~ -max_target_seqs ${params.hit}
    //~ # NCBI
    //~ #blastn -query ${fasta} -out ${fasta_bn}_nt.txt -outfmt \
           //~ "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
           //~ -db ${params.nt} \
           //~ -evalue ${params.evalue} -num_threads ${params.cpus} \
           //~ -max_target_seqs ${params.hit}
    //~ # rVDB
    //~ #blastn -query ${fasta} -out ${fasta_bn}_rvdb.txt -outfmt \
           //~ "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
           //~ -db ${params.rvdb} \
           //~ -evalue ${params.evalue} -num_threads ${params.cpus}
    //~ #cat ${fasta_bn}_rvdb.txt >> ${fasta_bn}_nt.txt
    //~ # Microbial catalogue
    //~ #blastn -query ${fasta} -out ${fasta_bn}_catalogue.txt -outfmt \
           //~ "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
           //~ -db ${params.catalogue}  \
           //~ -evalue ${params.evalue} -num_threads ${params.cpus}
    //~ """
//~ }

//~ blastChannel = Channel.fromPath("${baseDir}/../../binning_wf/Bacteria_V3/Annot/bin.*_nt.txt")
//~ binChannel_2 = Channel.fromPath("${baseDir}/../../binning_wf/Bacteria_V3/bins/bin.*.fa")

//~ process taxonomy {
    //~ publishDir "$params.annotDir", mode: 'copy'
    //~ memory "20G"
    
    //~ //echo true
    
    //~ input:
    //~ set file(fasta), file(nt) from blastChannel //, file(catalogue)   
    //~ //~~ file(nt) from blastChannel //, file(catalogue)   
    //~ //~~ file(fasta) from binChannel_2
    
    //~ output:
    //~ file("*_krona.txt") into taxChannel
    //~ file("*_not_annotated.fasta") into notAnnotatedChannel
    //~ //file("log.txt") into logChannel
    //~ file("*_{taxonomy,annotation}.txt") into resChannel mode flatten
    
    //~ script:
    //~ fasta_bn = fasta.baseName
    //~ """
    //~ #!/bin/bash
    //~ tax_count=`wc -l ${nt} | cut -f 1 -d " "`
    //~ if [ "\$tax_count" -gt "0" ]; then
        //~ # Annot ncbi
        //~ python !{params.scripts}/get_taxonomy.py -i ${nt} \
                //~ -o ${fasta_bn}_taxonomy.txt -t ${params.gitaxidnucl} \
                //~ -n ${params.names} -d ${params.nodes}
        //~ python !{params.scripts}/ExtractNCBIDB.py \
                //~ -f ${nt} -g ${fasta_bn}_taxonomy.txt -fc ${params.coverage} \
                //~ -o ${fasta_bn}_annotation.txt -nb 1
        //~ # Interest column for krona
        //~ cut -s -f 3-10 ${fasta_bn}_annotation.txt > ${fasta_bn}_annotation_interest.txt
        
        //~ # Get sequence not annotated
        //~ if [ -f ${fasta_bn}_catalogue_annotation.txt ]; then
            //~ cat  ${fasta_bn}_annotation.txt ${fasta_bn}_catalogue_annotation.txt\
             //~ > annotated
            //~ python !{params.scripts}/extract_fasta.py -q annotated \
                //~ -t ${fasta} -n -o ${fasta_bn}_not_annotated.fasta
        //~ else
            //~ python !{params.scripts}/extract_fasta.py \
                //~ -q ${fasta_bn}_annotation.txt -t ${fasta} -n \
                //~ -o ${fasta_bn}_not_annotated.fasta
        //~ fi
        
        
        //~ # Create Krona annotation
        //~ while read line; do
            //~ echo -e "1\t\${line}"
        //~ done < ${fasta_bn}_annotation_interest.txt > ${fasta_bn}_krona.txt
        //~ annot=`wc -l ${fasta_bn}_krona.txt | cut -f 1 -d ' '`
        //~ #echo "\$count_reads" > log.txt
        //~ #echo "\$annot" >>log.txt
        //~ # Count not annoted elements
        //~ if [ "\$count_reads" -gt "\$annot" ]; then
            //~ val=\$(( count_reads - annot))
        
            //~ echo -e "\$val\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> ${fasta_bn}_krona.txt
        //~ fi
    //~ else 
        //~ cat ${fasta} > ${fasta_bn}_not_annotated.fasta
        //~ touch ${fasta_bn}_krona.txt ${fasta_bn}_annotation.txt \
            //~ ${fasta_bn}_annotation.txt
    //~ fi
    //~ """
//~ }

