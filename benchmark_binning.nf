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
// a completer
def usage() {
    println("benchmark_binning.nf --in <contigs_dir> --out <output_dir> --cpus <nb_cpus>")
}


if(params.help){
    usage()
    exit(1)
}

//baseDir = dossier ou se trouve le script nextflow
params.in="$baseDir/../../simulated_data/Bacteria/"
readChannel = Channel.fromFilePairs("${params.in}/*_R{1,2}.fastq")
                     .ifEmpty { exit 1, "Cannot find any reads matching: ${params.in}" }
                     //.subscribe { println it }

/*contigsChannel = Channel
                .fromPath("${params.in}/*.fasta")
                .map {file -> tuple(file.baseName.replaceAll(".fasta",""),file)}*/

params.cpus = 6
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
params.binDir = "${params.out}/bins"
params.skiptrim = "F"
//~ params.chkmDir = "${params.out}/chkm_res"
//~ params.vp1 = "$baseDir/databases/vp1_seq.fasta"
//~ params.ncbi = "$baseDir/databases/ncbi_viruses.fna"
//~ params.rvdb = "$baseDir/databases/rVDBv10.2.fasta"
//~ params.uniprot = "$baseDir/databases/uniprot_taxonomy.fasta"
//~ params.uniref = "$baseDir/databases/uniref_uniprot.fasta"
//~ params.viral = "$baseDir/databases/viral_catalogue_poltson.fna"
//~ params.nt = "/pasteur/projets/policy01/BioIT/amine/catalogue/nt"
//~ params.gitaxidnucl = "/local/databases/release/taxodb/gi_taxid_nucl.dmp"
//~ params.names = "/local/databases/release/taxodb/names.dmp"
//~ params.nodes = "/local/databases/release/taxodb/nodes.dmp"


myDir = file(params.out)
myDir.mkdirs()

cleanDir = file("${params.cleaned_reads}")
cleanDir.mkdirs()

if(params.contaminant != "") {
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
        
        script:
        """
        AlienTrimmer -if ${forward} -ir ${reverse} -of ${pair_id}_1.fastq \
        -or ${pair_id}_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} \
        -l ${params.minlength}
        """
    }
    mappingChannel.subscribe { it.copyTo(cleanDir) }
    
    process khmer {
        cpus params.cpus
        memory "50000"
        
        input:
        set pair_id, file(forward), file(reverse) from trimChannel
        
        output:
        set pair_id, file("*_1.fastq"), file("*_2.fastq") into assemblyChannel
        
        script:
        """
        interleave-reads.py ${forward} ${reverse} --output interleaved.pe
        normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct  interleaved.pe --output output.pe.keep
        filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
        extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
        split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
        """
    }
}
else if (params.skiptrim == "F") {
    process trimming {
        //publishDir "$myDir", mode: 'copy'
        
        input:
        set pair_id, file(reads) from readChannel
        
        output:
        set pair_id, file("*_1.fastq"), file("*_2.fastq") into trimChannel
        file("*.fastq") into mappingChannel mode flatten
        
        script:
        """
        AlienTrimmer -if ${reads[0]} -ir ${reads[1]} -of ${pair_id}_1.fastq \
        -or ${pair_id}_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} \
        -l ${params.minlength}
        """
    }
    mappingChannel.subscribe { it.copyTo(cleanDir) }
    
    process khmer {
        cpus params.cpus
        memory "50000"
        
        input:
        set pair_id, file(forward), file(reverse) from trimChannel
        
        output:
        set pair_id, file("*_1.fastq"), file("*_2.fastq") into assemblyChannel
        
        script:
        """
        interleave-reads.py ${forward} ${reverse} --output interleaved.pe
        normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct  interleaved.pe --output output.pe.keep
        filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
        extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
        split-paired-reads.py output.dn.pe -1 ${pair_id}_filt_1.fastq -2 ${pair_id}_filt_2.fastq
        """
    }
}
else {
    trimChannel = Channel.fromFilePairs("${params.cleaned_reads}/*_{1,2}.fastq").ifEmpty { exit 1, "No clean reads were found"}
    
    process khmer {
        cpus params.cpus
        memory "50000"
        
        input:
        set pair_id, file(read) from trimChannel
        
        output:
        set pair_id, file("*_1.fastq"), file("*_2.fastq") into assemblyChannel
        
        script:
        """
        interleave-reads.py ${read[0]} ${read[1]} --output interleaved.pe
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
    
    cpus params.cpus
    if(params.mode == "clc"){
        clusterOptions='--qos=fast -C clcbio' 
    }
    input:
    set pair_id, file(forward), file(reverse) from assemblyChannel
    
    output:
    file("assembly/*_{spades,clc,minia}.fasta") into contigsChannel
    //file("assembly/*_{spades,clc,minia}.fasta") into contigsChannel_2
    
    shell:
    """
    #!/bin/bash
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
    file("cata_contig_nr.fasta") into cdhitChannel
    file("cata_contig_nr.fasta") into cdhitChannel_2
    
    shell:
    """
    bash /pasteur/homes/qletourn/scripts/merge_fasta.sh !{contigs}
    cd-hit-est -i cata_contigs.fasta -o cata_contig_nr.fasta -c 0.95 -T !{params.cpus} -aS 0.9 -d 0 -g 1 -M 0
    """
}

bowt_refDir = file(params.bowt_index)
bowt_refDir.mkdirs()

file(params.mappingDir).mkdirs()

if (params.bamDir == "" && params.index_prefix != "" && ! file("${params.bowt_index}/${params.index_prefix}.1.bt2").exists() ) {
    
    process index {
        publishDir "$bowt_refDir", mode: 'copy'
        cpus params.cpus
        
        input:
        file assembly from cdhitChannel
        
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
        file("res_mapping/bam/*.bam") into bamChannel mode flatten
        file("res_mapping/comptage/count_matrix.txt") into countChannel
        
        shell:
        """
        mbma.py mapping -i !{cleanDir} -o res_mapping -db !{params.bowt_index}/!{params.index_prefix} -t 6 -q fast --bowtie2 --shared -e quentin.letourneur@pasteur.fr
        cp -r res_mapping/* !{params.mappingDir}/
        rm -r !{params.mappingDir}/sam
        bash /pasteur/homes/qletourn/scripts/summarise_mapping_PE.sh !{params.mappingDir} !{params.mappingDir}/res_mapping.tsv
        """
    }
    
    
    process sort_bam {
    
        cpus params.cpus
        
        input:
        file bam from bamChannel
        
        output:
        file("*.sorted") into sortedChannel
        
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
        file contig from cdhitChannel
        
        output:
        file("res_mapping/bam/*.bam") into bamChannel mode flatten
        file("res_mapping/comptage/count_matrix.txt") into countChannel
        
        shell:
        """
        mbma.py mapping -i !{cleanDir} -o res_mapping -db !{params.bowt_index}/!{params.index_prefix} -t 6 -q fast --bowtie2 --shared -e quentin.letourneur@pasteur.fr
        cp -r res_mapping/* !{params.mappingDir}/
        rm -r !{params.mappingDir}/sam
        bash /pasteur/homes/qletourn/scripts/summarise_mapping_PE.sh !{params.mappingDir} !{params.mappingDir}/res_mapping.tsv
        """
    }
    
    
    process sort_bam {
    
        cpus params.cpus
        
        input:
        file bam from bamChannel
        
        output:
        file("*.sorted") into sortedChannel
        
        shell:
        """
        samtools sort -o !{bam}.sorted -@ !{params.cpus} -O bam !{bam}
        """
    }
}
else if ( params.bamDir != "" ) {
    
    bamChannel = Channel.fromPath("${params.bamDir}/*.bam")
    
    process sort_bam {
        
        cpus params.cpus
        
        input:
        file bam from bamChannel
        
        output:
        file("*.sorted") into sortedChannel
        
        shell:
        """
        samtools sort -o !{bam}.sorted -@ !{params.cpus} -O bam !{bam}
        """
    }
}
else {
    exit 1, "If no bam file path is specified you have to give a prefix name for the bowtie2 index files"
}
// maybe give choice to use jgi_summurized... or the count_matrix of mbma for the read count 

file(params.binDir).mkdirs()

process binning {
    publishDir "${params.binDir}", mode: 'copy'
    cpus params.cpus
    
    input:
    file assembly from cdhitChannel_2
    file bams from sortedChannel.toList()
    
    output:
    file("bin*") into binChannel
    
    shell:
    """
    /pasteur/homes/qletourn/tools/metabat/jgi_summarize_bam_contig_depths --outputDepth depth.txt !{bams}
    /pasteur/homes/qletourn/tools/metabat/metabat -i !{assembly} -a depth.txt -o bin -t !{params.cpus} --minSamples 5
    """
}


process annotaion {
    cpus params.cpus
    memory "100G"
    
    input:
    file(bins) from binChannel.first()
    
    output:
    file("chkm_res/tree_qa+qa.tsv") into annotChannel
    
    shell:
    """
    if ! mkdir chkm_res 2>/dev/null ; then
        rm -r chkm_res
        mkdir chkm_res
    fi
    
    checkm tree -t !{params.cpus} -x fa --tmpdir /pasteur/homes/qletourn/tmp_chkm !{params.binDir} chkm_res
    checkm tree_qa -f chkm_res/tree_qa.tsv --tab_table --tmpdir /pasteur/homes/qletourn/tmp_chkm chkm_res
    checkm lineage_set --tmpdir /pasteur/homes/qletourn/tmp_chkm chkm_res chkm_res/lineage.ms
    checkm analyze -t !{params.cpus} --tmpdir /pasteur/homes/qletourn/tmp_chkm -x fa chkm_res/lineage.ms !{params.binDir} chkm_res
    checkm qa -t !{params.cpus} -f chkm_res/qa_res.tsv --tab_table --tmpdir /pasteur/homes/qletourn/tmp_chkm lineage.ms chkm_res
     
    join -t $'\t' -1 1 -2 1 -o 1.1,1.4,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 <(tail -n +2 chkm_res/tree_qa.tsv | sort -k1,1) <(tail -n +2 chkm_res/qa_res.tsv | sort -k1,1) | sed '1s/^/Bin id\tTaxonomy\tMarker lineage\t# genomes\t# markers\tmarker sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain heterogeneity\n/' > chkm_res/tree_qa+qa.tsv
    
    cp -r chkm_res !{params.out}
    """
}

