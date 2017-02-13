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

workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}


workflow.onError = {
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

params.cpus = 8
params.out = "$baseDir/../../binning_wf/"
params.coverage = 80
params.mismatch = 1
params.alienseq = "$baseDir/../../genomes_databases/alienTrimmerPF8contaminants.fasta"
params.minlength = 45
params.cleaned_reads = "${params.out}/cleaned_reads"
params.mode = "spades"
params.contaminant = ""
//~ params.vp1 = "$baseDir/databases/vp1_seq.fasta"
//~ params.ncbi = "$baseDir/databases/ncbi_viruses.fna"
//~ params.rvdb = "$baseDir/databases/rVDBv10.2.fasta"
//~ params.uniprot = "$baseDir/databases/uniprot_taxonomy.fasta"
//~ params.uniref = "$baseDir/databases/uniref_uniprot.fasta"
//~ params.viral = "$baseDir/databases/viral_catalogue_poltson.fna"
//params.nt = "/local/databases/fasta/nt"
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
	    set pair_id, file(reads) from readChannel
	
	    output:
	    set pair_id, file("unmapped/*.1"), file("unmapped/*.2") into unmappedChannel
	    
	    shell:
	    """
	    mkdir unmapped
	    bowtie2 -q -N !{params.mismatch} -1 !{reads[0]} -2 !{reads[1]} -x ${params.contaminant} --un-conc unmapped/ -S /dev/null -p !{params.cpus}
	    """
	}


	process trimming {
	    //publishDir "$myDir", mode: 'copy'
	    cpus params.cpus
	    
	    input:
	    set pair_id, file(forward), file(reverse) from unmappedChannel
	    
	    
	    output:
	    set pair_id, file("*_1.fastq"), file("*_2.fastq") into trimChannel
	    
	    script:
		"""
		AlienTrimmer -if ${forward} -ir ${reverse} -of un-conc-mate_1.fastq 
		-or un-conc-mate_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} 
		-l ${params.minlength}
		"""
	}
}
else {
	process trimming {
	    //publishDir "$myDir", mode: 'copy'
	    cpus params.cpus
	    
	    input:
	    set set pair_id, file(reads) from readChannel
	    
	    
	    output:
	    set pair_id, file("*_1.fastq"), file("*_2.fastq") into trimChannel
	    
	    script:
		"""
		AlienTrimmer -if ${forward} -ir ${reverse} -of un-conc-mate_1.fastq 
		-or un-conc-mate_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} 
		-l ${params.minlength}
	    """
	}
}

process khmer {
    cpus params.cpus
    memory "50000"
    input:
    set pair_id, file(forward), file(reverse) from trimChannel
    
    output:
    file("*.fastq") into mappingChannel mode flatten
    set pair_id, file("*_1.fastq"), file("*_2.fastq") into assemblyChannel
    
    script:
    """
    interleave-reads.py ${forward} ${reverse} --output interleaved.pe
    normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct  interleaved.pe --output output.pe.keep
    filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter -T ${params.cpus}
    extract-paired-reads.py output.pe.filter --output-paired output.dn.pe  --output-single output.dn.se
    split-paired-reads.py output.dn.pe -1 ${pair_id}_1.fastq -2 ${pair_id}_2.fastq
    """
}

/*mappingChannel = Channel.create()
assemblyChannel = Channel.create()
khmerChannel.into { mappingChannel ; assemblyChannel }*/

/*Channel
  .fromFilePairs('*-{1,2}.txt') 
    .flatMap { key, files -> [[key], files].combinations() }
      .collectFile()*/


mappingChannel.subscribe { it.copyTo(cleanDir) }

process assembly {
    publishDir "$myDir", mode: 'copy'
    cpus params.cpus
    if(params.mode == "clc"){
        clusterOptions='--qos=fast -C clcbio' 
    }
    input:
    set pair_id, file(forward), file(reverse) from assemblyChannel

    output:
    set pair_id, file("assembly/*_{spades,clc,minia}.fasta") into contigsChannel mode flatten //test mode flatten

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
//voir si stockage des assemblages
//contigsChannel.subscribe { it.copyTo() }

