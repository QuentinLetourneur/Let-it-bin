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
readChannel = Channel.fromFilePairs("${params.in}/*_{1,2}.fastq")
                     .ifEmpty { exit 1, "Cannot find any reads matching: ${params.in}" }
                     //.subscribe { println it }

/*contigsChannel = Channel
                .fromPath("${params.in}/*.fasta")
                .map {file -> tuple(file.baseName.replaceAll(".fasta",""),file)}*/

params.cpus = 4
params.cpusassembly = 12
params.out = "$baseDir/../../binning_wf"
params.binDir = "${params.out}/Binnings"
params.coverage = 80
params.mismatch = 1
params.alienseq = "/pasteur/projets/policy01/Biomics/metagenomics/alienTrimmerPF8contaminants.fasta"
params.minlength = 45
params.cleaned_reads = "${params.out}/cleaned_reads"
params.mode = "spades"
params.contaminant = "" // prefix for the index of bowtie2 for contaminants
params.bowt_index = "$baseDir/../../bowtie_ref"
params.index_prefix = "" // prefix for the index of bowtie2 for analysed contigs
params.bamDir = ""
params.mappingDir = "${params.out}/mapping" // were mapping results will be stored
params.skiptrim = "F"
params.nt = "/pasteur/projets/policy01/Biomics/metagenomics/catalogue/nt"
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
params.likelybin = " "
params.all = "F"
params.metabatDir = "${params.binDir}/Metabat"
params.concoctDir = "${params.binDir}/Concoct"
params.cocacolaDir = "${params.binDir}/Cocacola"
params.myccDir = "${params.binDir}/MyCC"
params.groopmDir = "${params.binDir}/GroopM"
params.maxbinDir = "${params.binDir}/MaxBin"
params.metacluster5Dir = "${params.binDir}/Metacluster5"
params.mbbcDir = "${params.binDir}/MBBC"
params.canopyDir = "${params.binDir}/Canopy"
params.abundancebinDir = "${params.binDir}/AbundanceBin"
params.scimmDir = "${params.binDir}/SCIMM"
params.likelybinDir = "${params.binDir}/LikelyBin"
params.multi_assembly = "F"
params.contigs = ""
params.count_matrix = ""
params.scripts = "/pasteur/projets/policy01/BioIT/quentin/scripts"
params.interleaved = ""
params.refs_info = "/pasteur/projets/policy01/BioIT/quentin/refs_info.tsv"
params.genomes = "/pasteur/projets/policy01/BioIT/quentin/Blastdb/Bact_Anita"
params.tools = "/pasteur/projets/policy01/BioIT/quentin/tools"
params.email = "quentin.letourneur@pasteur.fr"
params.tmp_checkm = "/pasteur/scratch/amine/tmp_chkm"
params.min_bin_size = 1000000
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
        bowtie2 -q -N !{params.mismatch} -1 !{reads[0]} -2 !{reads[1]} -x ${params.contaminant} --un-conc unmapped/ -S /dev/null -p !{params.cpus} --very-sensitive-local
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

        script:
        """
        AlienTrimmer -if ${forward} -ir ${reverse} -of ${pair_id}_1.fastq \
        -or ${pair_id}_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} \
        -l ${params.minlength}
        """
    }
    mappingChannel.subscribe { it.copyTo(cleanDir) }

}
else if( params.skiptrim == "F" ) {
    process trimming {
        //publishDir "$myDir", mode: 'copy'

        input:
        set pair_id, file(reads) from readChannel

        output:
        set pair_id, file("*_alien_1.fastq"), file("*_alien_2.fastq") into trimChannel
        //~ file("*_1.fastq") into R1Channel
        //~ file("*_2.fastq") into R2Channel
        file("*.fastq") into mappingChannel mode flatten

        script:
        """
        AlienTrimmer -if ${reads[0]} -ir ${reads[1]} -of ${pair_id}_alien_1.fastq \
        -or ${pair_id}_alien_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} \
        -l ${params.minlength}
        """
    }
    mappingChannel.subscribe { it.copyTo(cleanDir) }
}
else {
    skiptrimChannel = Channel.fromFilePairs("${params.cleaned_reads}/*_{1,2}.fastq").ifEmpty { exit 1, "No clean reads were found"}
}

if( params.contigs == "" && params.multi_assembly == "F" ) {

    if( params.skiptrim == "F" ) {
        process khmer {
            cpus params.cpus
            //memory "100000"
            clusterOptions='--qos=normal -p common'
    
            input:
            set pair_id, file(fw), file(rv) from trimChannel
    
            output:
            //~ set pair_id, file("*t_1.fastq"), file("*t_2.fastq") into khmerChannel
            file("*_filt_1.fastq") into R1Channel
            file("*_filt_2.fastq") into R2Channel
            file("interleaved.pe") into conv_fastqChannel
    
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
    else {
        process khmer {
            cpus params.cpus
            //memory "100000"
            clusterOptions='--qos=normal -p common'
    
            input:
            set pair_id, file(clean_reads) from skiptrimChannel
    
            output:
            //~ set pair_id, file("*t_1.fastq"), file("*t_2.fastq") into khmerChannel
            file("*_filt_1.fastq") into R1Channel
            file("*_filt_2.fastq") into R2Channel
            file("interleaved.pe") into conv_fastqChannel
    
            script:
            """
            interleave-reads.py ${clean_reads[0]} ${clean_reads[1]} --output interleaved.pe
            normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct  interleaved.pe --output output.pe.keep
            filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
            extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
            split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
            """
        }
    }

    process merge {

        input:
        file R1 from R1Channel.toList()
        file R2 from R2Channel.toList()

        output:
        set val("mergefstq"), file("*_1.fastq"), file("*_2.fastq") into mergeChannel

        shell:
        """
        bash !{params.scripts}/merge_fastq.sh !{R1} !{R2}
        """
    }

    process assembly {
        // permet de copier dans myDir les fichiers de l'output
        publishDir "$myDir", mode: 'copy'

        cpus params.cpusassembly
        if(params.mode == "clc"){
            clusterOptions='--qos=fast -C clcbio'
        }
        input:
        set pair_id, file(forward), file(reverse) from mergeChannel

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
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_10
        file("assembly/*_{spades,clc,minia}.fasta") into assemblyChannel_11

        shell:
        """
        #!/bin/bash
        mkdir assembly
        if [ !{params.mode} == "spades" ];then
            spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpusassembly} -o assembly/
            mv assembly/scaffolds.fasta assembly/!{pair_id}_spades.fasta
        elif [ !{params.mode} ==  "clc" ];then
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpusassembly}
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
// ERROR when the script is run with --multi_assembly true but don't understand where it comes from
else if( params.contigs == "" && params.multi_assembly == "T" ) {

    if( params.skiptrim == "F" ) {
        process khmer {
            cpus params.cpus
            memory "60000"
            clusterOptions='--qos=normal'
            
            input:
            set pair_id, file(fw), file(rv) from trimChannel
            
            output:
            set pair_id, file("*_filt_1.fastq"), file("*_filt_2.fastq") into khmerChannel
            file("interleaved.pe") into conv_fastqChannel
            
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
    else {
        process khmer {
            cpus params.cpus
            memory "60000"
            clusterOptions='--qos=normal'
    
            input:
            set pair_id, file(clean_reads) from skiptrimChannel
    
            output:
            set pair_id, file("*_filt_1.fastq"), file("*_filt_2.fastq") into khmerChannel
            file("interleaved.pe") into conv_fastqChannel
            
            script:
            """
            interleave-reads.py ${clean_reads[0]} ${clean_reads[1]} --output interleaved.pe
            normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct  interleaved.pe --output output.pe.keep
            filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
            extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
            split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
            """
        }
    }
    
    process assembly {
        // permet de copier dans myDir les fichiers de l'output
        publishDir "$myDir", mode: 'copy'

        cpus params.cpusassembly
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
            spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpusassembly} -o assembly/
            mv assembly/scaffolds.fasta assembly/!{pair_id}_spades.fasta
        elif [ !{params.mode} ==  "clc" ]
        then
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpusassembly}
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
        file("cata_contig_nr.fasta") into assemblyChannel_11

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
    assemblyChannel_10 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_11 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    if( params.metacluster5 !=  " " ) {
        conv_fastqChannel = Channel.fromPath("${params.interleaved}")
                         .ifEmpty { exit 1, "You have to give the path to the interleaved.pe file if you want to run metacluster5 without doing the assembly" }
    }
}


bowt_refDir = file(params.bowt_index)
bowt_refDir.mkdirs()

//file(params.mappingDir).mkdirs()

if( params.bamDir == "" && params.index_prefix != "" && ! file("${params.bowt_index}/${params.index_prefix}.1.bt2").exists() ) {

    process index {
        publishDir "$bowt_refDir", mode: 'copy'
        // cpus params.cpus

        input:
        file assembly from assemblyChannel

        output:
        file("*.bt2") into indexChannel

        shell:
        """
        bowtie2-build !{assembly} !{params.index_prefix}
        """
    }

    process mapping_count {

        //cpus params.cpus

        modules = 'Python/2.7.8:samtools/1.3:mbma/tars'

        input:
        file idx from indexChannel.first()

        output:
        file("mapping/bam/*.bam") into bamChannel mode flatten
        file("mapping/comptage/count_matrix.txt") into countChannel

        shell:
        """
        #!/bin/bash
        mbma.py mapping -i !{cleanDir} -o mapping -db !{params.bowt_index}/!{params.index_prefix} -t 4 -q normal --bowtie2 --best -e !{params.email}
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
        mbma.py mapping -i !{cleanDir} -o mapping -db !{params.bowt_index}/!{params.index_prefix} -t 6 -q normal --bowtie2 --best -e !{params.email}
        bash !{params.scripts}/summarise_mapping_PE.sh mapping mapping/stats_mapping.tsv
        rm -r mapping/sam
        cp -r mapping/ !{params.out}/
        """
    }

    if( params.concoct != " " || params.cocacola != " " || params.groopm != " " || params.maxbin!= " " || params.metabat != " " || params.mycc != " " || params.all == "T" ) {
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
}
else if ( params.bamDir != "" && params.count_matrix != "" && (params.concoct != " " || params.cocacola != " " || params.groopm != " " || params.maxbin!= " " || params.metabat != " " || params.mycc != " " || params.all == "T") ) {

    bamChannel = Channel.fromPath("${params.bamDir}/*.bam").ifEmpty { exit 1, "Cannot find any bams in the directory : ${params.in}" }
    countChannel = Channel.fromPath("${params.count_matrix}").ifEmpty { exit 1, "Cannot find the count_matrix in the directory : ${params.count_matrix}" }

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
else if( params.likelybin == " " && params.scimm ==  " " && params.metacluster5 == " " && params.mbbc == " " && params.canopy ==  " " && params.abundancebin == " " ) {
    exit 1, "If no prefix for bowtie index is given you have to give the path to the bams directory and the count_matrix.txt file"
}

if ( params.concoct != " " || params.cocacola != " " || params.groopm != " " || params.maxbin!= " " || params.all == "T" ) {

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
}

if( params.concoct != " " || params.cocacola != " " || params.maxbin!= " " || params.all == "T" ) {
    process abun_and_link_profile {

        cpus params.cpus

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
        python !{params.tools}/CONCOCT-0.4.0/scripts/gen_input_table.py !{assembly} !{bam_sorted} > covTable.tsv
        python !{params.tools}/CONCOCT-0.4.0/scripts/bam_to_linkage.py -m !{params.cpus} --regionlength 500 --fullsearch !{assembly} !{bam_sorted} > linkTable.tsv
        """
    }
}

if( params.metacluster5 != " " ) {
    process fastq_to_fasta {

        input:
        file fastq from conv_fastqChannel.toList()

        output:
        file("merged.fa") into metaclustinputChannel

        shell:
        """
        args=!{fastq}
        if [ `awk '{print NF; exit}' \$args` -gt 1 ];then
            for file in \$args;do
                cat \$file >> merged.fastq
            done
            egrep "^@[a-z]+" -A 1 merged.fastq | sed -e 's/^@/>/' -e '/^--/d' > merged.fa
        else
            egrep "^@[a-z]+" -A 1 !{fastq} | sed -e 's/^@/>/' -e '/^--/d' > merged.fa
        fi
        """
    }
}

if( params.metabat != " " || params.mycc != " " || params.all == "T" ) {

    process jgi_summa_depths {

        input:
        file bams from sortedChannel.toList()

        output:
        file("depth.txt") into abunChannel
        file("depth.txt") into abunChannel_2

        shell:
        """
        !{params.tools}/metabat/jgi_summarize_bam_contig_depths --outputDepth depth.txt !{bams}
        """
    }
}

if( params.likelybin != " " ) {

    process LikelyBin {

        input:
        file assembly from assemblyChannel_10

        output:
        file("*") into likelybinChannel

        shell:
        """
        perl !{params.tools}/likelybin-master/mcmc.pl !{assembly} -num_source=!{params.nb_cluster}
        """
    }
}
else { likelybinChannel = Channel.from("a") }

if( params.scimm != " " || params.all == "T" ) {
    file(params.scimmDir).mkdirs()

    process SCIMM {
        cpus params.cpus

        input:
        file assembly from assemblyChannel_11

        output:
        file("bin*") into scimmChannel
        set val("SCIMM"), file("bin*") into scimmChannel_2
        val("a") into scimmChannel_3

        shell:
        """
        python !{params.tools}/scimm/bin/scimm.py -s !{assembly} -k !{params.nb_cluster} -p 2 --ln 1000 --cs 0 --ct 2
        for file in `ls cluster*.fa`;do
            name=`echo \$file | sed -e 's/cluster/bin/' -e 's/-/./'`
            mv \$file \$name
        done
        cp bin*.fa !{params.scimmDir}
        """
    }
}
else {
    scimmChannel = Channel.from("a")
    scimmChannel_2 = Channel.empty()
    scimmChannel_3 = Channel.from("a")
}

//~ if( params.abundancebin != " " || params.all == "T" ) {
    //~ file(params.abundancebinDir).mkdirs()

    //~ process AbundanceBin {

        //~ input:

        //~ output:
        //~ file ("*") into abundancebinChannel

        //~ shell:
        //~ """

        //~ """
    //~ }
//~ }
//~ else { abundancebinChannel = Channel.from("a") }

if( params.canopy != " " || params.all == "T" ) {
    file(params.canopyDir).mkdirs()

    process Canopy {
        cpus params.cpus
        //memory "100G"

        input:
        file count_mat from countChannel
        file assembly from assemblyChannel_9

        output:
        file("bin*") into canopyChannel
        set val("Canopy"), file("bin*") into canopyChannel_2
        val("a") into canopyChannel_3

        shell:
        '''
        sed -r s/'(\t[0-9]+)\\.[0-9]+'/'\\1'/g !{count_mat} > tmp.txt
        tail -n +2 tmp.txt | cut -f 1,3- -d $'\t' > count_matrix.tsv
        cc.bin -n !{params.cpus} -i count_matrix.tsv -o cluster.out
        nb_clust=`sort -rnk1,1 cluster.out | head -n 1 | egrep -o "CAG[0-9]+" | egrep -o "[0-9]+"`
        for clust in `seq $nb_clust`;do
            egrep "CAG0*\$clust" cluster.out | cut -f 2 > list.txt
            python !{params.scripts}/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        done
        cp bin* !{params.canopyDir}
        '''
    }
}
else {
    canopyChannel = Channel.from("a")
    canopyChannel_2 = Channel.empty()
    canopyChannel_3 = Channel.from("a")
}

//~ if( params.mbbc != " " || params.all == "T" ) {
    //~ file(params.mbbcDir).mkdirs()

    //~ process MBBC {
        //~ memory "60000"

        //~ input:


        //~ output:
        //~ file("*") into mbbcChannel

        //~ shell:
        //~ """
        //~ java -jar -Xmx50g -i
        //~ """
    //~ }
//~ }
//~ else { mbbcChannel = Channel.from("a") }

//~ if( params.metacluster5 != " " || params.all == "T" ) {
    //~ file(params.metacluster5Dir).mkdirs()

    //~ process Metacluster5 {

        //~ input:
        //~ file reads from metaclustinputChannel.toList()

        //~ output:
        //~ file("*") into metacluster5Channel

        //~ shell:
        //~ """
        //~ MetaCluster5_1 !{reads}
        //~ MetaCluster5_2 !{reads}
        //~ """
    //~ }
//~ }
//~ else { metacluster5Channel = Channel.from("a") }

if( params.maxbin != " " || params.all == "T" ) {
    file(params.maxbinDir).mkdirs()

    process maxbin {
        cpus params.cpus
        module = 'gcc/4.9.0:bowtie2/2.2.3:hmmer/3.1b1:FragGeneScan/1.30:idba/1.1.1:MaxBin/2.2.3'

        input:
        file assembly from assemblyChannel_8
        file abun from abundanceProfileChannel_3

        output:
        file("bin*") into maxbinChannel
        set val("MaxBin"), file("bin*") into maxbinChannel_2
        val("a") into maxbinChannel_3

        shell:
        """
        nb_col=`awk '{print NF; exit}' !{abun}`
        bash !{params.scripts}/split_abun_file.sh !{abun} \$nb_col
        ls cov_mean* > list_abun_files.txt
        run_MaxBin.pl -contig !{assembly} -out out -abund_list list_abun_files.txt -thread !{params.cpus}
        for file in `ls out.0*`;do
            name=`echo \$file | cut -f 1,2 -d "." | sed 's/out/bin/'`
            mv \$file \$name.fa
        done
        cp bin* !{params.maxbinDir}
        """
    }
}
else {
    maxbinChannel = Channel.from("a")
    maxbinChannel_2 = Channel.empty()
    maxbinChannel_3 = Channel.from("a")
}

//~ if( params.groopm != " " || params.all == "T" ) {
    //~ file(params.groopmDir).mkdirs()

    //~ process GroopM {
        //~ cpus params.cpus

        //~ input:
        //~ file assembly from assemblyChannel_7
        //~ file bams from sortedChannel_4.toList()
        //~ file bams_idx from indexedChannel_2.toList()

        //~ output:
        //~ file("*") into groopmChannel

        //~ shell:
        //~ """
        //~ groopm2 parse -t !{params.cpus} -c 1000 database.gm !{assembly} !{bams}
        //~ groopm2 core database.gm
        //~ groopm2 extract -t !{params.cpus} database.gm !{assembly} -o GroopM
        //~ """
        //~~ groopm parse -t !{params.cpus} -c 1000 database.gm !{assembly} !{bams}
        //~~ groopm core database.gm
        //~~ groopm recruit -c 1000
        //~~ groopm extract -t !{params.cpus} database.gm !{assembly} -o GroopM
    //~ }
//~ }
//~ else { groopmChannel = Channel.from("a") }

if( params.metabat != " " || params.all == "T" ) {
    file(params.metabatDir).mkdirs()

    process Metabat {
        publishDir "${params.metabatDir}", mode: 'copy'
        cpus params.cpus

        input:
        file assembly from assemblyChannel_3
        file depth from abunChannel

        output:
        file("bin*") into metabatChannel
        set val("Metabat"), file("bin*") into metabatChannel_2
        val("a") into metabatChannel_3
        
        shell:
        """
        !{params.tools}/metabat/metabat -i !{assembly} -a depth.txt -o bin -t !{params.cpus} --minSamples 5
        """
    }
}
else {
    metabatChannel = Channel.from("a")
    metabatChannel_2 = Channel.empty()
    metabatChannel_3 = Channel.from("a")
}

if( params.concoct != " " || params.all == "T" ) {
    //~ file(params.concoctDir).mkdirs()

    //~ process CONCOCT {
        //~ cpus params.cpus

        //~ input:
        //~ file assembly from assemblyChannel_4
        //~ file abun_prof from abundanceProfileChannel
        //~ file link from linkageChannel

        //~ output:
        //~ file("bin*") into concoctChannel
        //~ set val("Concoct"), file("bin*") into concoctChannel_2
        //~ val("a") into concoctChannel_3

        //~ shell:
        //~ """
        //~ cut -f1,3- !{abun_prof} > covTableR.tsv
        //~ concoct -c !{params.nb_cluster} --coverage_file covTableR.tsv --composition_file !{assembly} -b Concoct/

        //~ nb_clust=`sort -t "," -rnk2,2 Concoct/clustering_gt1000.csv | head -n 1 | cut -f 2 -d ","`
        //~ for clust in `seq \$nb_clust`;do
            //~ grep ",\$clust" Concoct/clustering_gt1000.csv | cut -f 1 -d "," > list.txt
            //~ python !{params.scripts}/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        //~ done
        //~ cp bin* !{params.concoctDir}
        //~ """
    //~ }
}
else {
    concoctChannel = Channel.from("a")
    concoctChannel_2 = Channel.empty()
    concoctChannel_3 = Channel.from("a")
}

if( params.cocacola != " " || params.all == "T" ) {
    //~ file(params.cocacolaDir).mkdirs()

    //~ process COCACOLA {
        //~~ cpus params.cpus

        //~ input:
        //~ file abun_table from abundanceProfileChannel_2
        //~ file link from linkageChannel_2
        //~ file assembly from assemblyChannel_5

        //~ output:
        //~ file("bin*") into cocacolaChannel
        //~ set val("Cocacola"), file("bin*") into cocacolaChannel_2
        //~ val("a") into cocacolaChannel_3

        //~ shell:
        //~ """
        //~ nb_contigs=`grep "^>" !{assembly} | wc -l`
        //~ python !{params.tools}/CONCOCT-0.4.0/scripts/fasta_to_features.py !{assembly} \$nb_contigs 4 kmer_4.csv
        //~ python !{params.tools}/COCACOLA-python/cocacola.py --contig_file !{assembly} --abundance_profiles !{abun_table} --composition_profiles kmer_4.csv --out result.csv --aux_dir !{params.tools}/COCACOLA-python
        //~ nb_clust=`sort -t "," -rnk2,2 Concoct/clustering_gt1000.csv | head -n 1 | cut -f 2 -d ","`
        //~ for clust in `seq \$nb_clust`;do
            //~ grep ",\$clust" Concoct/clustering_gt1000.csv | cut -f 1 -d "," > list.txt
            //~ python !{params.scripts}/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        //~ done
        //~ cp bin* !{params.cocacolaDir}
        //~ """
    //~ }
}
//~ else {
    //~ cocacolaChannel = Channel.from("a")
    //~ cocacolaChannel_2 = Channel.empty()
    //~ cocacolaChannel_3 = Channel.from("a")
//~ }

if( params.mycc != " " || params.all == "T" ) {
    //~ file(params.myccDir).mkdirs()

    //~ process MyCC {

        //~ input:
        //~ file assembly from assemblyChannel_6
        //~ file coverage from abunChannel_2

        //~ output:
        //~ file ("bin*") into myccChannel
        //~ set val("MyCC"), file("bin*") into myccChannel_2
        //~ val("a") into myccChannel_3

        //~ shell:
        //~ """
        //~ cut -f1,4- !{coverage} > abun.txt
        //~ python !{params.tools}/MyCC/MyCC.py !{assembly} -a abun.txt -meta
        //~ """
    //~ }
}
else {
    myccChannel = Channel.from("a")
    myccChannel_2 = Channel.empty()
    myccChannel_3 = Channel.from("a")
}

process checkm {
    cpus params.cpus
    memory "100G"
    clusterOptions='--qos=normal -p common'


    input:
    file bins from scimmChannel.first()
    //~~ file bins_1 from abundancebinChannel.first()
    file bins_2 from canopyChannel.first()
    //~~ file bins_3 from mbbcChannel.first()
    //~~ file bins_4 from metacluster5Channel.first()
    file bins_5 from maxbinChannel.first()
    //~~ file bins_6 from groopmChannel.first()
    file bins_7 from metabatChannel.first()
    //~ file bins_8 from concoctChannel.first()
    //~ file bins_9 from cocacolaChannel.first()
    //~ file bins_10 from myccChannel.first()
    //~~ file bins_11 from likelybinChannel.first()

    //~ output:
    //~ file("chkm_res/tree_qa+qa.tsv") into annotChannel

    shell:
    '''
    for tool in `ls !{params.binDir}/`;do
        if ! mkdir chkm_res 2>/dev/null ; then
            rm -r chkm_res
            mkdir chkm_res
        fi


        checkm tree -t !{params.cpus} -x fa --tmpdir !{params.tmp_checkm} !{params.binDir}/\$tool chkm_res
        checkm tree_qa -f chkm_res/tree_qa.tsv --tab_table --tmpdir !{params.tmp_checkm} chkm_res
        checkm lineage_set --tmpdir !{params.tmp_checkm} chkm_res chkm_res/lineage.ms
        checkm analyze -t !{params.cpus} --tmpdir !{params.tmp_checkm} -x fa chkm_res/lineage.ms !{params.binDir}/\$tool chkm_res
        checkm qa -t !{params.cpus} -f chkm_res/qa_res.tsv --tab_table --tmpdir !{params.tmp_checkm} chkm_res/lineage.ms chkm_res

        echo -e 'Bin id\tTaxonomy\tMarker lineage\t# genomes\t# markers\tmarker sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain heterogeneity' > chkm_res/tree_qa+qa.tsv
        join -t $'\t' -1 1 -2 1 -o 1.1,1.4,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 <(tail -n +2 chkm_res/tree_qa.tsv | sort -k1,1) <(tail -n +2 chkm_res/qa_res.tsv | sort -k1,1) >> chkm_res/tree_qa+qa.tsv

        cp -r chkm_res !{params.binDir}/\$tool
    done
    '''
}


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

process binnings_done {

    input:
    val inp1 from scimmChannel_3
    val inp2 from canopyChannel_3
    val inp3 from maxbinChannel_3
    //~ val inp4 from concoctChannel_3
    val inp5 from metabatChannel_3
    //~ val inp6 from cocacolaChannel_3
    //~ val inp7 from myccChannel_3

    output:
    val("") into concatChannel

    shell:
    """
    echo "plop"
    """
}

blastinputChannel = concatChannel.concat( scimmChannel_2, canopyChannel_2, maxbinChannel_2, metabatChannel_2 ) //, concoctChannel_2, cocacolaChannel_2, myccChannel_2

process blast {
    //publishDir "$params.annotDir", mode: 'copy'
    publishDir "$params.annotDir", mode: 'copy'
    cpus params.cpus
    memory "20G"

    input:
    set val(soft), file(fasta) from blastinputChannel

    output:
    set val(soft), file(fasta), file("*_nt.txt") into blastChannel //, file("*_catalogue.txt")

    when:
    soft != ""

    script:
    // fasta_bn = fasta.baseName
    """
    #!/bin/bash
    mkdir -p ${params.binDir}/${soft}/Annotation
    for infile in ${fasta};do
        fasta_bn=`echo \$infile | cut -f 1,2 -d "."`
        # local DB
        blastn -query \$infile -out \$fasta_bn"_nt.txt" -outfmt \
               "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
               -db ${params.genomes} \
               -evalue ${params.evalue} -num_threads ${params.cpus} \
               -max_target_seqs ${params.hit}
        # NCBI
        #blastn -query \$fasta_bn.fa -out \$fasta_bn"_nt.txt" -outfmt \
               "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
               -db ${params.nt} \
               -evalue ${params.evalue} -num_threads ${params.cpus} \
               -max_target_seqs ${params.hit}
        # rVDB
        #blastn -query \$fasta_bn.fa -out \$fasta_bn"_rvdb.txt" -outfmt \
               "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
               -db ${params.rvdb} \
               -evalue ${params.evalue} -num_threads ${params.cpus}
        #cat \$fasta_bn"_rvdb.txt" >> \$fasta_bn"_nt.txt"
        # Microbial catalogue
        #blastn -query \$fasta_bn.fa -out \$fasta_bn"_catalogue.txt" -outfmt \
               "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
               -db ${params.catalogue}  \
               -evalue ${params.evalue} -num_threads ${params.cpus}
    done
    cp *_nt.txt ${params.binDir}/${soft}/Annotation
    """
}

//~~ blastChannel = Channel.fromPath("${baseDir}/../../binning_wf/Bacteria_V3/Annot/bin.*_nt.txt")
//~~ binChannel_2 = Channel.fromPath("${baseDir}/../../binning_wf/Bacteria_V3/bins/bin.*.fa")

// process taxonomy {
//     //publishDir "$params.annotDir", mode: 'copy'
//     memory "20G"

//     //~ clusterOptions='--qos=normal'

//     input:
//     set soft, file(fasta), file(nt) from blastChannel //, file(catalogue)
//     //~~ file(nt) from blastChannel //, file(catalogue)
//     //~~ file(fasta) from binChannel_2

//     output:
//     //file("*_krona.txt") into taxChannel
//     //file("*_not_annotated.fasta") into notAnnotatedChannel
//     //file("log.txt") into logChannel
//     file("*_{taxonomy,annotation}.txt") into resChannel mode flatten

//     script:
    // fasta_bn = fasta.baseName
    // """
    // #!/bin/bash
    // for infile in ${nt};do
    //     tax_count=`wc -l \$infile | cut -f 1 -d " "`
    //     fasta_bn=`echo \$infile | cut -f 1,2 -d "." | cut -f 1 -d "_"`
    //     if [ "\$tax_count" -gt "0" ]; then
    //         # Annot ncbi
    //         python ${params.scripts}/get_taxonomy.py -i \$infile \
    //                 -o \$fasta_bn"_taxonomy.txt" -t ${params.gitaxidnucl} \
    //                 -n ${params.names} -d ${params.nodes}
    //         python ${params.scripts}/ExtractNCBIDB.py \
    //                 -f \$infile -g \$fasta_bn"_taxonomy.txt" -fc ${params.coverage} \
    //                 -o \$fasta_bn"_annotation.txt" -nb 1
    //         # Interest column for krona
    //         cut -s -f 3-10 \$fasta_bn"_annotation.txt" > \$fasta_bn"_annotation_interest.txt"
    //         count_reads=\$(grep "^>" -c \$fasta_bn.fa)

    //         # Get sequence not annotated
    //         if [ -f \$fasta_bn"_catalogue_annotation.txt" ]; then
    //             cat  \$fasta_bn"_annotation.txt" \$fasta_bn"_catalogue_annotation.txt"\
    //              > annotated
    //             python ${params.scripts}/extract_fasta.py -q annotated \
    //                 -t \$fasta_bn.fa -n -o \$fasta_bn"_not_annotated.fasta"
    //         else
    //             python ${params.scripts}/extract_fasta.py \
    //                 -q \$fasta_bn"_annotation.txt" -t \$fasta_bn.fa -n \
    //                 -o \$fasta_bn"_not_annotated.fasta"
    //         fi


    //         # Create Krona annotation
    //         while read line; do
    //             echo -e "1\t\${line}"
    //         done < \$fasta_bn"_annotation_interest.txt" > \$fasta_bn"_krona.txt"
    //         annot=`wc -l \$fasta_bn"_krona.txt" | cut -f 1 -d ' '`
    //         #echo "\$count_reads" > log.txt
    //         #echo "\$annot" >>log.txt
    //         # Count not annoted elements
    //         if [ "\$count_reads" -gt "\$annot" ]; then
    //             val=\$(( count_reads - annot))

    //             echo -e "\$val\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> \$fasta_bn"_krona.txt"
    //         fi
    //     else
    //         cat \$fasta_bn.fa > \$fasta_bn"_not_annotated.fasta"
    //         touch \$fasta_bn"_krona.txt" \$fasta_bn"_annotation.txt"
    //     fi
    // done
    // cp *_krona.txt ${params.binDir}/${soft}/Annotation
    // cp *_not_annotated.fasta ${params.binDir}/${soft}/Annotation
    // cp *_taxonomy.txt ${params.binDir}/${soft}/Annotation
    // cp *_annotation.txt ${params.binDir}/${soft}/Annotation
    // """
// }

// process evaluation {

//     input:
//     file annotations from resChannel.toList()

//     output:
//     file("contamination.png") into finishChannel

//     shell:
//     """
//     for tool in `ls !{params.binDir}/`;do
//         Rscript !{params.scripts}/binning_stats.R !{params.binDir}/\$tool/bin !{params.binDir}/\$tool !{params.binDir}/\$tool/Annotation !{params.refs_info} !{params.min_bin_size}
//         cp contamination.png completeness.png !{params.binDir}/\$tool
//     done
//     """
// }

