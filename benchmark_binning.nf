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
params.cpusmapping = 4
params.cpusassembly = 12
params.cpusscimm = params.cpus
params.out = "$baseDir/../../binning_wf"
params.binDir = "${params.out}/Binnings"
params.coverage = 90
params.identity = 95
params.mismatch = 1
params.alienseq = "/pasteur/projets/policy01/Biomics/metagenomics/alienTrimmerPF8contaminants.fasta"
params.minlength = 45
params.cleaned_reads = "${params.out}/cleaned_reads"
params.mode = "spades"
params.contaminant = "" // prefix for the index of bowtie2 for contaminants
params.bowt_index = "/pasteur/projets/policy01/BioIT/quentin/bowtie_ref"
params.index_prefix = "" // prefix for the index of bowtie2 for analysed contigs
params.bamDir = ""
params.mappingDir = "${params.out}/mapping" // were mapping results will be stored
params.skiptrim = "F"
params.nt = "/pasteur/projets/policy01/Biomics/metagenomics/catalogue/nt"
params.evalue = 10
params.rvdb = "/pasteur/projets/policy01/Biomics/metagenomics/catalogue/rVDBv10.2.fasta"
params.gitaxidnucl = "/local/databases/release/taxodb/gi_taxid_nucl.dmp"
//~ params.annotDir = "${params.out}/Annot"
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
params.metagen = " "
params.binsanity = " "
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
params.metagenDir = "${params.binDir}/MetaGen"
params.binsanityDir = "${params.binDir}/BinSanity"
params.multi_assembly = "F"
params.contigs = ""
params.count_matrix = ""
params.scripts = "/pasteur/projets/policy01/BioIT/quentin/scripts"
params.interleaved = ""
params.refs_info = "/pasteur/projets/policy01/BioIT/quentin/refs_info.tsv"
params.lcs = "" // take the absolute path to the file containing the link between reference contigs/genome ID and species name
params.blast_db = ""
params.tools = "/pasteur/projets/policy01/BioIT/quentin/tools"
params.email = "quentin.letourneur@pasteur.fr"
params.tmp_checkm = "/pasteur/scratch/amine/tmp_chkm"
params.min_bin_size = 1000000
params.conta_threshold = 0.9 // at max 10% contamination
params.qos_assembly = "normal" //should be either fast, normal (default) or long
params.filt_readsDir = ""
params.plot_width = "750"
params.plot_heigth = "1000"
params.metagen_src = ""
params.min_contig_length = 1000
params.sim_data = "T"
params.blast_ncbi = "F"
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
    // Map reads on the genome(s) of organisms that could have been sequenced along the wanted genome.
    // Only reads that don't map on it/them will be kept 
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
    
    // reads cleaning
    process trimming {
        //publishDir "$cleanDir", mode: 'copy'
        
        input:
        set pair_id, file(forward), file(reverse) from unmappedChannel
        
        output:
        set pair_id, file("*_1.fastq"), file("*_2.fastq") into trimChannel
        file("*_1.fastq") into R1Channel
        file("*_2.fastq") into R2Channel
        file("Ech*.fastq") into mappingChannel mode flatten
        
        script:
        """
        AlienTrimmer -if ${forward} -ir ${reverse} -of ${pair_id}_alien_1.fastq \
        -or ${pair_id}_alien_2.fastq -os un-conc-mate_sgl.fastq -c ${params.alienseq} \
        -l ${params.minlength}
        """
    }
    mappingChannel.subscribe { it.copyTo(cleanDir) }

}
else if( params.skiptrim == "F" ) {
    // reads cleaning
    process trimming {
        //publishDir "$myDir", mode: 'copy'
        
        input:
        set pair_id, file(reads) from readChannel
        
        output:
        set pair_id, file("*_alien_1.fastq"), file("*_alien_2.fastq") into trimChannel
        file("*_1.fastq") into R1Channel
        file("*_2.fastq") into R2Channel
        file("Ech*.fastq") into mappingChannel mode flatten
        //~ file("Ech*.fastq") into clean_readsChannel
        
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
    R1Channel = Channel.fromPath("${params.cleaned_reads}/*_1.fastq").ifEmpty { exit 1, "No forward clean reads were found"}
    R2Channel = Channel.fromPath("${params.cleaned_reads}/*_2.fastq").ifEmpty { exit 1, "No reverse clean reads were found"}
    //~ clean_readsChannel = Channel.fromPath("${params.cleaned_reads}/Ech*.fastq").ifEmpty { exit 1, "No clean reads file were found"}
}

// Perform a co-assembly
if( params.contigs == "" && params.multi_assembly == "F" ) {

    if( params.skiptrim == "F" && params.filt_readsDir == "" ) {
        // filter reads redoudancy and even their aboudance
        process khmer {
            cpus params.cpus
            memory "100000"
            clusterOptions='--qos=normal -p common'
            
            input:
            set pair_id, file(fw), file(rv) from trimChannel
            
            output:
            //~ set pair_id, file("*t_1.fastq"), file("*t_2.fastq") into khmerChannel
            file("*_filt_1.fastq") into R1_Channel
            file("*_filt_2.fastq") into R2_Channel
            //~ file("interleaved.pe") into conv_fastqChannel
            
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
    else if( params.filt_readsDir == "" ) {
        // filter reads redoudancy and even their aboudance
        process khmer {
            cpus params.cpus
            memory "100000"
            clusterOptions='--qos=normal -p common'
            
            input:
            set pair_id, file(clean_reads) from skiptrimChannel
            
            output:
            //~ set pair_id, file("*t_1.fastq"), file("*t_2.fastq") into khmerChannel
            file("*_filt_1.fastq") into R1_Channel
            file("*_filt_2.fastq") into R2_Channel
            //~ file("interleaved.pe") into conv_fastqChannel
            
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
    else {
        R1_Channel = Channel.fromPath("${params.filt_readsDir}/*filt_1.fastq").ifEmpty { exit 1, "No forward filt reads were found"}
        R2_Channel = Channel.fromPath("${params.filt_readsDir}/*filt_2.fastq").ifEmpty { exit 1, "No reverse filt reads were found"}
    }
    
    // Pull all the forward/reverse reads together to make a co-assembly
    process merge {
        
        input:
        file R1 from R1_Channel.toList()
        file R2 from R2_Channel.toList()
        
        output:
        set val("mergefstq"), file("mergefstq_1.fastq"), file("mergefstq_2.fastq") into mergeChannel
        
        shell:
        """
        cat !{R1} > mergefstq_1.fastq
        cat !{R2} > mergefstq_2.fastq
        """
    }
    
    process assembly {
        // permet de copier dans myDir les fichiers de l'output
        publishDir "$myDir", mode: 'copy'
        
        if(params.mode != "ray") {
            cpus params.cpusassembly
            memory "150G"
        }
        else {
            cpus 1
            clusterOptions="-p dedicated --qos=hubbioit"
            exitReadTimeout = '604800sec'
        }
        
        if(params.mode == "clc") {
            clusterOptions="--qos=clcbio"
        }
        else if(params.mode != "ray") {
            clusterOptions="--qos=${params.qos_assembly}"
        }
        //~ else if(params.qos_assembly == "l") {
            //~ clusterOptions="--qos=long"
        //~ }
        //~ else if(params.mode == "minia" || params.mode == "megahit" || params.mode == "ray"){
            //~ clusterOptions='--qos=long'
        //~ }
        
        input:
        set pair_id, file(forward), file(reverse) from mergeChannel
        
        output:
        file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel
        file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_2
        file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_3
        //~ file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_4
        //~ file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_5
        //~ file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_6
        //~ file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_7
        file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_8
        file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_9
        //~ file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_10
        file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_11
        file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_12
        //~ file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_13
        //~ file("assembly/*_{spades,clc,minia,megahit,ray}.fasta") into assemblyChannel_14
        
        shell:
        """
        #!/bin/bash
        if [ !{params.mode} == "spades" ];then
            mkdir assembly
            spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpusassembly} -o assembly/
            #mv assembly/scaffolds.fasta assembly/!{pair_id}_spades.fasta
            mv assembly/contigs.fasta assembly/contigs_spades.fasta
        elif [ !{params.mode} ==  "clc" ];then
            mkdir assembly
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpusassembly}
            #clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpusassembly} -n
            mv assembly/contigs.fasta assembly/!{pair_id}_clc.fasta
        elif [ !{params.mode} == "minia" ];then
            mkdir assembly
            /pasteur/projets/policy01/Matrix/metagenomics/entvirus/bin/gatb-minia-pipeline/gatb -1 !{forward} -2 !{reverse} -o assembly/!{pair_id}_minia.fasta
        elif [ !{params.mode} == "megahit" ];then
            megahit -m 0.5 --mem-flag 0 -t !{params.cpusassembly} -1 !{forward} -2 !{reverse} -o assembly --out-prefix !{pair_id}_megahit
            mv assembly/!{pair_id}_megahit.contigs.fa assembly/!{pair_id}_megahit.fasta
            sed -i "s/ /_/g" !{pair_id}_megahit.fasta
        elif [ !{params.mode} == "ray" ];then
            sbatch -n !{params.cpusassembly} --mem-per-cpu 30000 --qos=!{params.qos_assembly} -J ray --wait -o ray.log --wrap="mpirun -n !{params.cpusassembly} Ray -k 31 -p !{forward} !{reverse} -o assembly"
            #mv assembly/Scaffolds.fasta assembly/!{pair_id}_ray.fasta
            mv assembly/Contigs.fasta assembly/!{pair_id}_ray.fasta
            sed -i "s/ /_/g" !{pair_id}_ray.fasta
        else
            echo "ERROR : the mode should be either spades, clc, minia, megahit or ray"
            exit 1
        fi
        """
    }
    
    //~ if( params.sim_data == "T" && params.blast_db != "" && params.lcs != "") {
    
        //~ process contigs_annotation {
            //~ cpus params.cpus
            
            //~ input:
            //~ file assembly from assemblyChannel_14
            
            //~ output:
            //~ file("contigs_annotation.txt") into annotationChannel
            
            //~ shell:
            //~ """
            //~ blastn -query !{assembly} -out blast_contigs_on_refs.txt -outfmt
                   //~ "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore"
                   //~ -db !{params.blast_db}
                   //~ -evalue !{params.evalue} -num_threads !{params.cpus}
                   //~ -max_target_seqs !{params.hit}
            
            //~ !{params.scripts}/ExtractAnnotation.py -f blast_contigs_on_refs.txt -a !{params.lcs} -o contigs_annotation.txt -nb 1 -r . -fc !{params.coverage} -fi !{params.identity}
            
            //~ Rscript !{params.scripts}/sum_contigs_length_per_annotation.R contigs_annotation.txt sum_contigs_length_per_annotation.tsv
            
            //~ cp contigs_annotation sum_contigs_length_per_annotation.tsv !{params.out}/assembly
            //~ """
        //~ }
    //~ }
    //~ else {
        //~ exit 1, "If you are running the pipeline on simulated data, make sure you have specified --blast_db path to the BLAST database files and --lcs a file containing the association between sequence ids and the genome they come from and a reference number"
    //~ }
}

// perform a multi-assembly
else if( params.contigs == "" && params.multi_assembly == "T" ) {

    if( params.skiptrim == "F" && params.filt_readsDir == "" ) {
        // filter reads redoudancy and even their aboudance
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
    else if( params.filt_readsDir == "" ) {
        // filter reads redoudancy and even their aboudance
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
    else {
        R1_Channel = Channel.fromPath("${params.filt_readsDir}/*filt_1.fastq").ifEmpty { exit 1, "No forward filt reads were found"}
        R2_Channel = Channel.fromPath("${params.filt_readsDir}/*filt_2.fastq").ifEmpty { exit 1, "No reverse filt reads were found"}
    }
    
    // make an assembly per sample
    process assembly {
        // permet de copier dans myDir les fichiers de l'output
        publishDir "$myDir", mode: 'copy'
        
        if(params.mode != "ray") {
            cpus params.cpusassembly
            memory "70G"
        }
        else {
            cpus 1
        }
        
        if(params.mode == "clc") {
            clusterOptions='--qos=${params.qos_assembly} -C clcbio'
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
            spades.py --meta -1 !{reads[0]} -2 !{reads[1]} -t !{params.cpusassembly} -o assembly/
            mv assembly/scaffolds.fasta assembly/!{pair_id}_spades.fasta
        elif [ !{params.mode} ==  "clc" ];then
            clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{reads[0]} !{reads[1]} --cpus !{params.cpusassembly}
            mv assembly/contigs.fasta assembly/!{pair_id}_clc.fasta
        elif [ !{params.mode} == "minia" ];then
            /pasteur/projets/policy01/Matrix/metagenomics/entvirus/bin/gatb-minia-pipeline/gatb -1 !{forward} -2 !{reverse} -o assembly/!{pair_id}_minia.fasta
        elif [ !{params.mode} == "megahit" ];then
            megahit -m 0.5 --mem-flag 0 -t !{params.cpusassembly} -1 !{forward} -2 !{reverse} -o assembly/ --out-prefix !{pair_id}_megahit
            mv assembly/!{pair_id}_megahit.contigs.fa assembly/!{pair_id}_megahit.fasta
            sed -i "s/ /_/g" !{pair_id}_megahit.fasta 
        elif [ !{params.mode} == "ray" ];then
            sbatch -n !{params.cpusassembly} --mem-per-cpu 10000 --qos=!{params.qos_assembly} -J ray --wait -o ray.log --wrap="mpirun -n !{params.cpusassembly} Ray -k 31 -p !{forward} !{reverse} -o assembly"
            mv assembly/Scaffolds.fasta assembly/!{pair_id}_ray.fasta
        else
            echo "ERROR : the mode should be either spades, clc, minia, megahit or ray"
            exit 1
        fi
        """
        //~ #spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpusassembly} -o assembly/
        //~ #clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q -i !{forward} !{reverse} --cpus !{params.cpusassembly}
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
        //~ file("cata_contig_nr.fasta") into assemblyChannel_4
        //~ file("cata_contig_nr.fasta") into assemblyChannel_5
        //~ file("cata_contig_nr.fasta") into assemblyChannel_6
        //~ file("cata_contig_nr.fasta") into assemblyChannel_7
        file("cata_contig_nr.fasta") into assemblyChannel_8
        file("cata_contig_nr.fasta") into assemblyChannel_9
        //~ file("cata_contig_nr.fasta") into assemblyChannel_10
        file("cata_contig_nr.fasta") into assemblyChannel_11
        file("cata_contig_nr.fasta") into assemblyChannel_12
        //~ file("cata_contig_nr.fasta") into assemblyChannel_13
        
        shell:
        """
        cat !{contigs} > cata_contigs.fasta
        !{params.tools}/cdhit-master/cd-hit-est -i cata_contigs.fasta -o cata_contig_nr.fasta -c 0.95 -T !{params.cpus} -aS 0.9 -d 0 -g 1 -M 0
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
    //~ assemblyChannel_4 = Channel.fromPath("${params.contigs}")
                     //~ .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    //~ assemblyChannel_5 = Channel.fromPath("${params.contigs}")
                     //~ .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    //~ assemblyChannel_6 = Channel.fromPath("${params.contigs}")
                     //~ .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    //~ assemblyChannel_7 = Channel.fromPath("${params.contigs}")
                     //~ .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_8 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_9 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    //~ assemblyChannel_10 = Channel.fromPath("${params.contigs}")
                     //~ .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_11 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    assemblyChannel_12 = Channel.fromPath("${params.contigs}")
                     .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    //~ assemblyChannel_13 = Channel.fromPath("${params.contigs}")
                     //~ .ifEmpty { exit 1, "Cannot find the contigs file: ${params.contigs}" }
    if( params.metacluster5 !=  " " ) {
        conv_fastqChannel = Channel.fromPath("${params.interleaved}")
                         .ifEmpty { exit 1, "You have to give the path to the interleaved.pe file if you want to run metacluster5 without doing the assembly" }
    }
}

bowt_refDir = file(params.bowt_index)
bowt_refDir.mkdirs()

//file(params.mappingDir).mkdirs()

if( params.bamDir == "" && params.index_prefix != "" && ! file("${params.bowt_index}/${params.index_prefix}.1.bt2").exists() ) {
    
    // generate bowtie2 index of the contigs
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
    
    // map cleaned reads on the contigs
    process mapping_count {
        
        queue = 'dedicated'
        clusterOptions = '--qos=hubbioit'
        module 'Python/2.7.8:samtools/1.3:mbma/nextflow'
        
        input:
        file R1 from R1Channel.toList()
        file R2 from R2Channel.toList()
        file idx from indexChannel.first()
        
        output:
        //~ file("mapping/bam/Ech*.bam") into bamChannel mode flatten
        file("mapping/bam/sorted*.bam") into sortedChannel mode flatten
        file("mapping/bam/sorted*.bam") into sortedChannel_2 mode flatten
        file("mapping/bam/sorted*.bam") into sortedChannel_3
        file("mapping/bam/sorted*.bam") into sortedChannel_4
        file("mapping/bam/sorted*.bam") into sortedChannel_5
        file("mapping/bam/sorted*.bam") into sortedChannel_6
        file("mapping/comptage/count_matrix.txt") into countChannel
        
        shell:
        """
        #!/bin/bash
        mbma.py mapping -o mapping -db !{params.bowt_index}/!{params.index_prefix} -t !{params.cpusmapping} -q normal --bowtie2 --best -e !{params.email} --r1 !{R1} --r2 !{R2} -p common
        bash !{params.scripts}/summarise_mapping_PE.sh mapping mapping/stats_mapping.tsv
        rm -r mapping/sam mapping/bam/filtered*
        cp -r mapping/ !{params.out}/
        """
    }
    
    //~ process sort_bam {
        
        //~ cpus params.cpus

        //~ input:
        //~ file bam from bamChannel

        //~ output:
        //~ file("*.sorted") into sortedChannel
        //~ file("*.sorted") into sortedChannel_2
        //~ file("*.sorted") into sortedChannel_3
        //~ file("*.sorted") into sortedChannel_4

        //~ shell:
        //~ """
        //~ samtools sort -o !{bam}.sorted -@ !{params.cpus} -O bam !{bam}
        //~ """
    //~ }
}
else if (params.bamDir == "" && params.index_prefix != "") {
    
    // map cleaned reads on the contigs
    process mapping_count {
    
        queue = 'dedicated'
        clusterOptions = '--qos=hubbioit'
        module 'Python/2.7.8:samtools/1.3:mbma/nextflow'
        
        input:
        file R1 from R1Channel.toList()
        file R2 from R2Channel.toList()
        
        output:
        //~ file("mapping/bam/Ech*.bam") into bamChannel mode flatten
        file("mapping/bam/sorted*.bam") into sortedChannel mode flatten
        file("mapping/bam/sorted*.bam") into sortedChannel_2 mode flatten
        file("mapping/bam/sorted*.bam") into sortedChannel_3
        file("mapping/bam/sorted*.bam") into sortedChannel_4
        file("mapping/bam/sorted*.bam") into sortedChannel_5
        file("mapping/bam/sorted*.bam") into sortedChannel_6
        file("mapping/comptage/count_matrix.txt") into countChannel
        
        shell:
        """
        #!/bin/bash
        mbma.py mapping -o mapping -db !{params.bowt_index}/!{params.index_prefix} -t !{params.cpusmapping} -q normal --bowtie2 --best -e !{params.email} --r1 !{R1} --r2 !{R2} -p common
        bash !{params.scripts}/summarise_mapping_PE.sh mapping mapping/stats_mapping.tsv
        rm -r mapping/sam mapping/bam/filtered*
        cp -r mapping/ !{params.out}/
        """
    }

    //~ if( params.concoct != " " || params.cocacola != " " || params.groopm != " " || params.maxbin!= " " || params.metabat != " " || params.mycc != " " || params.all == "T" ) {
        
        //~ process sort_bam {

            //~ cpus params.cpus

            //~ input:
            //~ file bam from bamChannel

            //~ output:
            //~ file("*.sorted") into sortedChannel
            //~ file("*.sorted") into sortedChannel_2
            //~ file("*.sorted") into sortedChannel_3
            //~ file("*.sorted") into sortedChannel_4

            //~ shell:
            //~ """
            //~ samtools sort -o !{bam}.sorted -@ !{params.cpus} -O bam !{bam}
            //~ """
        //~ }
    //~ }
}
else if ( params.bamDir != "" && params.count_matrix != "" && (params.concoct != " " || params.cocacola != " " || params.maxbin!= " " || params.metabat != " " || params.mycc != " " || params.metagen || params.binsanity || params.all == "T") ) {

    //~ bamChannel = Channel.fromPath("${params.bamDir}/Ech*.bam").ifEmpty { exit 1, "Cannot find any BAM files in the directory : ${params.bamDir}" }
    sortedChannel = Channel.fromPath("${params.bamDir}/sorted*.bam").ifEmpty { exit 1, "Cannot find any sorted BAM files in the directory : ${params.bamDir}" }
    sortedChannel_2 = Channel.fromPath("${params.bamDir}/sorted*.bam").ifEmpty { exit 1, "Cannot find any sorted BAM files in the directory : ${params.bamDir}" }
    sortedChannel_3 = Channel.fromPath("${params.bamDir}/sorted*.bam").ifEmpty { exit 1, "Cannot find any sorted BAM files in the directory : ${params.bamDir}" }
    sortedChannel_4 = Channel.fromPath("${params.bamDir}/sorted*.bam").ifEmpty { exit 1, "Cannot find any sorted BAM files in the directory : ${params.bamDir}" }
    sortedChannel_5 = Channel.fromPath("${params.bamDir}/sorted*.bam").ifEmpty { exit 1, "Cannot find any sorted BAM files in the directory : ${params.bamDir}" }
    sortedChannel_6 = Channel.fromPath("${params.bamDir}/sorted*.bam").ifEmpty { exit 1, "Cannot find any sorted BAM files in the directory : ${params.bamDir}" }
    countChannel = Channel.fromPath("${params.count_matrix}").ifEmpty { exit 1, "Cannot find the count_matrix in the directory : ${params.count_matrix}" }
    
    //~ process sort_bam {

        //~ cpus params.cpus

        //~ input:
        //~ file bam from bamChannel

        //~ output:
        //~ file("*.sorted") into sortedChannel
        //~ file("*.sorted") into sortedChannel_2
        //~ file("*.sorted") into sortedChannel_3
        //~ file("*.sorted") into sortedChannel_4

        //~ shell:
        //~ """
        //~ samtools sort -o !{bam}.sorted -@ !{params.cpus} -O bam !{bam}
        //~ """
    //~ }
}
else {
    exit 1, "If the mapping of the reads on the contigs as already been done. You have to give the path to the BAM files directory and sorted BAM files and the path to the count_matrix.txt file and be sure to have specified the wanted binning softwares to run"
}

if( params.groopm != " " || params.metagen != " " || params.binsanity != " " || params.all == "T" ) {
    
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
        """
    }
}

if( params.metagen != " " ) {
    
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
    
    // make abundance tables based on the coverage files and link tables from alignment files
    process abun_and_link_profile {
        
        cpus params.cpus
        
        input:
        file bed_file from gen_cov_bedChannel.toList()
        file assembly from assemblyChannel_2
        
        output:
        //~ file("covTable.tsv") into abundanceProfileChannel
        //~ file("covTable.tsv") into abundanceProfileChannel_2
        file("covTable.tsv") into abundanceProfileChannel_3
        //file("linkTable.tsv") into linkageChannel
        
        shell:
        """
        python !{params.tools}/CONCOCT-0.4.0/scripts/gen_input_table.py !{assembly} !{bed_file} --isbedfiles > covTable.tsv
        """
    }
}

//~ if( params.metacluster5 != " " ) {
    //~ process fastq_to_fasta {

        //~ input:
        //~ file fastq from conv_fastqChannel.toList()

        //~ output:
        //~ file("merged.fa") into metaclustinputChannel

        //~ shell:
        //~ """
        //~ if [ `awk '{print NF; exit}' \$args` -gt 1 ];then
            //~ cat !{fastq} > merged.fastq
            //~ egrep "^@[a-z]+" -A 1 merged.fastq | sed -e 's/^@/>/' -e '/^--/d' > merged.fa
        //~ else
            //~ egrep "^@[a-z]+" -A 1 !{fastq} | sed -e 's/^@/>/' -e '/^--/d' > merged.fa
        //~ fi
        //~ """
    //~ }
//~ }

if( params.metabat != " " || params.mycc != " " || params.all == "T" ) {
    
    // make abundance file for metabat
    process jgi_summa_depths {
        
        input:
        file bams from sortedChannel_3.toList()
        
        output:
        file("depth.txt") into abunChannel
        //~ file("depth.txt") into abunChannel_2
        
        shell:
        """
        !{params.tools}/metabat/jgi_summarize_bam_contig_depths --outputDepth depth.txt !{bams}
        """
    }
}

//~ if( params.likelybin != " " ) {

    //~ process LikelyBin {

        //~ input:
        //~ file assembly from assemblyChannel_10

        //~ output:
        //~ file("*") into likelybinChannel

        //~ shell:
        //~ """
        //~ perl !{params.tools}/likelybin-master/mcmc.pl !{assembly} -num_source=!{params.nb_cluster}
        //~ """
    //~ }
//~ }
//~ else { likelybinChannel = Channel.from("a") }

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
        val("SCIMM") into scimmChannel_4
        
        shell:
        """
        samtools faidx !{assembly}
        cut -f1-2 !{assembly}.fai > contigs_length.tsv
        awk -v min_size=!{params.min_contig_length} '\$2 >= min_size {print \$1}' contigs_length.tsv > headers_contigs_gt1000.txt
        python !{params.scripts}/extract_fasta_from_list.py !{assembly} headers_contigs_gt1000.txt assembly_gt1000.fasta
        
        python !{params.tools}/scimm/bin/scimm.py -s assembly_gt1000.fasta -k !{params.nb_cluster} -p !{params.cpusscimm} --cs 0 --lt !{params.cpusscimm}
        
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
    scimmChannel_4 = Channel.empty()
}

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
        val("Canopy") into canopyChannel_4
        
        shell:
        '''
        sed -r s/'(\t[0-9]+)\\.[0-9]+'/'\\1'/g !{count_mat} > tmp.txt
        awk -v min_size=!{params.min_contig_length} '\$2 >= min_size {print \$0}' tmp.txt > tmp2.txt
        tail -n +2 tmp2.txt | cut -f 1,3- -d $'\t' > count_matrix_gt1000.tsv
        #tail -n +2 tmp.txt | cut -f 1,3- -d $'\t' > count_matrix.tsv
        
        cc.bin -n !{params.cpus} -i count_matrix_gt1000.tsv -o cluster.out
        #cc.bin -n !{params.cpus} -i count_matrix.tsv -o cluster.out
        
        sed -i "s/>//" cluster.out
        nb_clust=`cut -f 1 cluster.out | uniq | wc -l`
        for clust in `seq $nb_clust`;do
            grep -P "CAG0*$clust\t" cluster.out | cut -f 2 > list.txt
            python !{params.scripts}/extract_fasta_from_list.py !{assembly} list.txt bin.$clust.fa
        done
        cp bin* !{params.canopyDir}
        '''
    }
}
else {
    canopyChannel = Channel.from("a")
    canopyChannel_2 = Channel.empty()
    canopyChannel_3 = Channel.from("a")
    canopyChannel_4 = Channel.empty()
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
        val("MaxBin") into maxbinChannel_4
        
        shell:
        """
        nb_col=`awk '{print NF; exit}' !{abun}`
        bash !{params.scripts}/split_abun_file.sh !{abun} \$nb_col
        ls cov_mean* > list_abun_files.txt
        
        run_MaxBin.pl -contig !{assembly} -out out -abund_list list_abun_files.txt -thread !{params.cpus} -min_contig_length !{params.min_contig_length}
        
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
    maxbinChannel_4 = Channel.empty()
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
        val("Metabat") into metabatChannel_4
        
        shell:
        """
        if [ !{params.min_contig_length} -le 1500 ];then
            !{params.tools}/metabat/metabat -i !{assembly} -a depth.txt -o bin -t !{params.cpus} --minSamples 5
        else
            !{params.tools}/metabat/metabat -i !{assembly} -a depth.txt -o bin -t !{params.cpus} --minSamples 5 -m !{params.min_contig_length}
        fi
        """
    }
}
else {
    metabatChannel = Channel.from("a")
    metabatChannel_2 = Channel.empty()
    metabatChannel_3 = Channel.from("a")
    metabatChannel_4 = Channel.empty()
}

//~ if( params.concoct != " " || params.all == "T" ) {
    //~ file(params.concoctDir).mkdirs()
    
    //~ process CONCOCT {
        //~ cpus params.cpus
        
        //~ input:
        //~ file assembly from assemblyChannel_4
        //~ file abun_prof from abundanceProfileChannel
        
        //~ output:
        //~ file("bin*") into concoctChannel
        //~ set val("Concoct"), file("bin*") into concoctChannel_2
        //~ val("a") into concoctChannel_3
        
        //~ shell:
        //~ """
        //~ cut -f1,3- !{abun_prof} > covTableR.tsv
        //~ concoct -c !{params.nb_cluster} --coverage_file covTableR.tsv --composition_file !{assembly} -b Concoct/ -l !{params.min_contig_length}
        
        //~ nb_clust=`sort -t "," -rnk2,2 Concoct/clustering_gt1000.csv | head -n 1 | cut -f 2 -d ","`
        //~ for clust in `seq 0 \$nb_clust`;do
            //~ grep ",\$clust\$" Concoct/clustering_gt1000.csv | cut -f 1 -d "," > list.txt
            //~ python !{params.scripts}/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        //~ done
        //~ cp bin* !{params.concoctDir}
        //~ """
    //~ }
//~ }
//~ else {
    //~ concoctChannel = Channel.from("a")
    //~ concoctChannel_2 = Channel.empty()
    //~ concoctChannel_3 = Channel.from("a")
//~ }

//~ if( params.cocacola != " " || params.all == "T" ) {
    //~ file(params.cocacolaDir).mkdirs()

    //~ process COCACOLA {
        //~ cpus params.cpus
        
        //~ input:
        //~ file abun_table from abundanceProfileChannel_2
        //~ file assembly from assemblyChannel_5
        
        //~ output:
        //~ file("bin*") into cocacolaChannel
        //~ set val("Cocacola"), file("bin*") into cocacolaChannel_2
        //~ val("a") into cocacolaChannel_3
        
        //~ shell:
        //~ """
        //~ awk -v min_size=!{params.min_contig_length} '\$2 >= min_size {print \$0}' !{abun_table} > covtable_gt1000.tsv
        //~ tail -n +2 covtable_gt1000.tsv | awk '{print \$1}' > headers.txt
        
        //~ nb_contigs=`grep "^>" !{assembly} | wc -l`
        //~ python !{params.tools}/CONCOCT-0.4.0/scripts/fasta_to_features.py !{assembly} \$nb_contigs 4 kmer_4.csv
        //~ head -n1 kmer_4.csv > kmer_4_gt1000.csv
        //~ grep -Fwf headers.txt kmer_4.csv >> kmer_4_gt1000.csv
        
        //~ python !{params.tools}/COCACOLA-python/cocacola.py --contig_file !{assembly} --abundance_profiles covtable_gt1000.tsv --composition_profiles kmer_4_gt1000.csv --out result_gt1000.csv --aux_dir !{params.tools}/COCACOLA-python --threads !{params.cpus}
        
        //~ num_clust=`sort -unk2,2 -t "," result_gt1000.csv | cut -f 2 -d ","`
        //~ for clust in \$num_clust;do
            //~ grep ",\$clust\$" result_gt1000.csv | cut -f 1 -d "," > list.txt
            //~ python !{params.scripts}/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        //~ done
        //~ #cp bin* !{params.cocacolaDir}
        //~ """
    //~ }
//~ }
//~ else {
    //~ cocacolaChannel = Channel.from("a")
    //~ cocacolaChannel_2 = Channel.empty()
    //~ cocacolaChannel_3 = Channel.from("a")
//~ }

//~ if( params.mycc != " " || params.all == "T" ) {
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
//~ }
//~ else {
    //~ myccChannel = Channel.from("a")
    //~ myccChannel_2 = Channel.empty()
    //~ myccChannel_3 = Channel.from("a")
//~ }

if( params.metagen != " " || params.all == "T" ) {
    
    file(params.metagenDir).mkdirs()
    
    process MetaGen {
        
        cpus params.cpus
        module = 'R/3.3.2:bowtie2/2.2.6:samtools/1.3:MetaGen/1.1.0'
        
        input:
        file assembly from assemblyChannel_12
        file count from statChannel.toList()
        //~ file reads from clean_readsChannel.toList()
        
        output:
        file("bin*") into metagenChannel
        set val("MetaGen"), file("bin*") into metagenChannel_2
        val("a") into metagenChannel_3
        val("MetaGen") into metagenChannel_4
        
        shell:
        """
        mkdir map output contigs reads
        ln !{assembly} contigs/Contigs.fasta
        
        for f in !{count};do
            bn=`basename \$f .bam.stat`
            mkdir map/\$bn
            ln \$f map/\$bn/count.dat
        done
        
        !{params.tools}/MetaGen-master/script/combine-counts.sh -p .
        !{params.tools}/MetaGen-master/script/sum-reads.sh -p !{params.cleaned_reads} .
        Rscript !{params.tools}/MetaGen-master/R/metagen.R -w . -n !{params.cpus} -l !{params.min_contig_length} -m !{params.tools}/MetaGen-master
        
        sed -i -e 's/"//g' -e 's/>//' -e 's/ /,/' output/segs.txt
        num_clust=`sort -unk2,2 -t "," output/segs.txt | cut -f 2 -d ","`
        for clust in \$num_clust;do
            grep ",\$clust\$" output/segs.txt | cut -f 1 -d "," > list.txt
            python !{params.scripts}/extract_fasta_from_list.py !{assembly} list.txt bin.\$clust.fa
        done
        cp bin* !{params.metagenDir}
        """
    }
}
else {
    metagenChannel = Channel.from("a")
    metagenChannel_2 = Channel.empty()
    metagenChannel_3 = Channel.from("a")
    metagenChannel_4 = Channel.empty()
}

//~ if( params.binsanity != " " || params.all == "T" ) {
    
    //~ file(params.binsanityDir).mkdirs()
    
    //~ process Binsanity {
        //~ cpus params.cpus
        //~ memory "100G"
        
        //~ input:
        //~ file assembly from assemblyChannel_13
        //~ file sorted_bam from sortedChannel_5.toList()
        //~ file bam_index from indexedChannel_3.toList()
        
        //~ output:
        //~ file("out.txt") into binsanityChannel
        //~ set val("BinSanity"), file("out.txt") into binsanityChannel_2
        //~ val("a") into binsanityChannel_3
        
        //~ shell:
        //~ """
        //~ python !{params.tools}/BinSanity-master/utils/get-ids -f . -l !{assembly} -o ids.txt -x !{params.min_contig_length}
        //~ python !{params.tools}/BinSanity-master/bin/Binsanity-profile -i !{assembly} -s . --ids ids.txt -c matrix -T !{params.cpus}
        //~ python !{params.tools}/BinSanity-master/bin/Binsanity-refine -f . -c matrix.cov -l !{assembly} -o . -x !{params.min_contig_length}
        //~ """
    //~ }
//~ }
//~ else {
    //~ binsanityChannel = Channel.from("a")
    //~ binsanityChannel_2 = Channel.empty()
    //~ binsanityChannel_3 = Channel.from("a")
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
// voir si ajouter un fichier bidon dans la sortie pour éviter le warning de cardinalité
process binnings_done {

    input:
    val inp1 from scimmChannel_3
    val inp2 from canopyChannel_3
    val inp3 from maxbinChannel_3
    //~ val inp4 from concoctChannel_3
    val inp5 from metabatChannel_3
    //~ val inp6 from cocacolaChannel_3
    //~ val inp7 from myccChannel_3
    val inp8 from metagenChannel_3
    //~ val inp9 from binsanityChannel_3

    output:
    val("") into concatChannel
    val("") into concatChannel_2

    shell:
    """
    echo "plop"
    """
}

checkmInputChannel = concatChannel_2.concat( scimmChannel_4, canopyChannel_4, maxbinChannel_4, metabatChannel_4, metagenChannel_4)

tmp_dir_checkm = file(params.tmp_checkm)
tmp_dir_checkm.mkdirs()

// evaluate the quality of the bins generated by the chosen binning softwares
process checkm {
    cpus params.cpus
    memory "100G"
    clusterOptions='--qos=normal -p common'
    
    input:
    val(soft) from checkmInputChannel
    
    output:
    file("chkm_res/tree_qa+qa.tsv") into checkmChannel
    val("plup") into resumeChannel_2
    
    when:
    soft != ""
    
    shell:
    '''
    if ! mkdir chkm_res 2>/dev/null ; then
        rm -r chkm_res
        mkdir chkm_res
    fi
    
    checkm tree -t !{params.cpus} -x fa --tmpdir !{params.tmp_checkm} !{params.binDir}/!{soft} chkm_res
    checkm tree_qa -f chkm_res/tree_qa.tsv --tab_table --tmpdir !{params.tmp_checkm} chkm_res
    checkm lineage_set --tmpdir !{params.tmp_checkm} chkm_res chkm_res/lineage.ms
    checkm analyze -t !{params.cpus} --tmpdir !{params.tmp_checkm} -x fa chkm_res/lineage.ms !{params.binDir}/!{soft} chkm_res
    checkm qa -t !{params.cpus} -f chkm_res/qa_res.tsv --tab_table --tmpdir !{params.tmp_checkm} chkm_res/lineage.ms chkm_res

    echo -e 'Bin id\tTaxonomy\tMarker lineage\t# genomes\t# markers\tmarker sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain heterogeneity' > chkm_res/tree_qa+qa.tsv
    join -t $'\t' -1 1 -2 1 -o 1.1,1.4,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 <(tail -n +2 chkm_res/tree_qa.tsv | sort -k1,1) <(tail -n +2 chkm_res/qa_res.tsv | sort -k1,1) >> chkm_res/tree_qa+qa.tsv

    cp -r chkm_res !{params.binDir}/!{soft}
    '''
}

blastInputChannel = concatChannel.concat( scimmChannel_2, canopyChannel_2, maxbinChannel_2, metabatChannel_2, metagenChannel_2) //, binsanityChannel_2 , concoctChannel_2, cocacolaChannel_2

process blast {
    //~ publishDir "$params.annotDir", mode: 'copy'
    cpus params.cpus
    memory "20G"
    
    input:
    set val(soft), file(fasta) from blastInputChannel
    
    output:
    set val(soft), file("*_nt.txt") into blastChannel //file(fasta), file("*_catalogue.txt")
    
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
               -db ${params.blast_db} \
               -evalue ${params.evalue} -num_threads ${params.cpus} \
               -max_target_seqs ${params.hit}
        # NCBI
        #blastn -query \$fasta_bn.fa -out \$fasta_bn"_nt.txt" -outfmt \
               "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" \
               -db ${params.blast_db} \
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

if( params.lcs != "" ) {
    process taxonomy {
        //~ publishDir "$params.annotDir", mode: 'copy'
        memory "20G"
        
        //~ clusterOptions='--qos=normal'
        
        input:
        set soft, file(nt) from blastChannel //file(fasta), file(catalogue)
        //~~ file(nt) from blastChannel //, file(catalogue)
        //~~ file(fasta) from binChannel_2
        
        output:
        //file("*_krona.txt") into taxChannel
        //file("*_not_annotated.fasta") into notAnnotatedChannel
        //file("log.txt") into logChannel
        //file("*_{taxonomy,annotation}.txt") into resChannel mode flatten
        set val(soft), file("*_annotation.txt") into taxoChannel
        
        script:
        //fasta_bn = fasta.baseName
        """
        #!/bin/bash
        for infile in ${nt};do
            tax_count=`wc -l < \$infile`
            if [ "\$tax_count" -gt "0" ]; then
                bin_name=`echo \$infile | egrep -o "bin.[0-9]+"`
                ${params.scripts}/ExtractAnnotation.py -f \$infile -a ${params.lcs} -o \$bin_name"_annotation.txt" -nb 1 -r . -fc ${params.coverage} -fi ${params.identity}
            fi
        done
        
        cp *_annotation.txt ${params.binDir}/${soft}/Annotation
        """
        
        //~ # Annot ncbi
        //~ #fasta_bn=`echo \$infile | cut -f 1,2 -d "." | cut -f 1 -d "_"`
        //~ python ${params.scripts}/get_taxonomy.py -i \$infile \
                //~ -o \$fasta_bn"_taxonomy.txt" -t ${params.gitaxidnucl} \
                //~ -n ${params.names} -d ${params.nodes}
        //~ python ${params.scripts}/ExtractNCBIDB.py \
                //~ -f \$infile -g \$fasta_bn"_taxonomy.txt" -fc ${params.coverage} \
                //~ -o \$fasta_bn"_annotation.txt" -nb 1
        //~ # Interest column for krona
        //~ cut -s -f 3-10 \$fasta_bn"_annotation.txt" > \$fasta_bn"_annotation_interest.txt"
        //~ count_reads=\$(grep "^>" -c \$fasta_bn.fa)
    
        //~ # Get sequence not annotated
        //~ if [ -f \$fasta_bn"_catalogue_annotation.txt" ]; then
            //~ cat  \$fasta_bn"_annotation.txt" \$fasta_bn"_catalogue_annotation.txt"\
             //~ > annotated
            //~ python ${params.scripts}/extract_fasta.py -q annotated \
                //~ -t \$fasta_bn.fa -n -o \$fasta_bn"_not_annotated.fasta"
        //~ else
            //~ python ${params.scripts}/extract_fasta.py \
                //~ -q \$fasta_bn"_annotation.txt" -t \$fasta_bn.fa -n \
                //~ -o \$fasta_bn"_not_annotated.fasta"
        //~ fi
        
        
        //~ # Create Krona annotation
        //~ while read line; do
            //~ echo -e "1\t\${line}"
        //~ done < \$fasta_bn"_annotation_interest.txt" > \$fasta_bn"_krona.txt"
        //~ annot=`wc -l \$fasta_bn"_krona.txt" | cut -f 1 -d ' '`
        //~ #echo "\$count_reads" > log.txt
        //~ #echo "\$annot" >>log.txt
        //~ # Count not annoted elements
        //~ if [ "\$count_reads" -gt "\$annot" ]; then
            //~ val=\$(( count_reads - annot))
            
            //~ echo -e "\$val\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> \$fasta_bn"_krona.txt"
        //~ fi
        //~ else
            //~ cat \$fasta_bn.fa > \$fasta_bn"_not_annotated.fasta"
            //~ touch \$fasta_bn"_krona.txt" \$fasta_bn"_annotation.txt"
            
        //~ cp *_krona.txt ${params.binDir}/${soft}/Annotation
        //~ cp *_not_annotated.fasta ${params.binDir}/${soft}/Annotation
        //~ cp *_taxonomy.txt ${params.binDir}/${soft}/Annotation
        //~ cp *_annotation.txt ${params.binDir}/${soft}/Annotation
    }
}
else {
    exit 1, "If you want to extract BLAST best hits and replace contigs name by the name of the species to which they belong specify the --lcs option with the path to the file annotation_num.tsv"
}

process eval_complet_conta {
    
    input:
    set val(soft), file(annotations) from taxoChannel
    
    output:
    file("contamination.png") into evalChannel
    val("plop") into resumeChannel
    
    shell:
    """
    Rscript !{params.scripts}/binning_stats_dominant.R !{params.binDir}/!{soft}/bin !{params.binDir}/!{soft} !{params.binDir}/!{soft}/Annotation !{params.refs_info} !{params.min_bin_size} !{params.plot_width} !{params.plot_heigth}
    cp contamination.png completeness.png !{params.binDir}/!{soft}
    """
}

//~ process resume_res {
    
    //~ input:
    //~ file prec_rec_blast from resumeChannel.toList()
    //~ file prec_rec_checkm from resumeChannel_2.toList()
    
    //~ output:
    //~ file("barplot*.png") into finishChannel
    
    //~ shell:
    //~ """
    //~ Rscript !{params.scripts}/nb_bin_per_threshold.R !{params.binDir} !{params.conta_threshold}
    //~ cp barplot*.png !{params.binDir}
    //~ """
//~ }

