#!/bin/bash
#Check arguments
if [ $# -ne 2  ]
then
    echo "$0 <mapping_dir> <output_file>"
    exit
fi

echo -e "Sample\t# reads\tUnique pair reads\tUnique pair reads (%)\tMultiple pair reads\tMultiple pair reads (%)\tNot mapped pair reads\tNot mapped pair reads (%)\tDiscordant pair read\tDiscordant pair read (%)\t# reads from unmapped pairs\tUnique read\tUnique read (%)\tMultiple read\tMultiple read (%)\tNot mapped reads\tNot mapped reads (%)\tOverall alignment rate (%)" > $2
for logfile in $(ls $1/out/*)
do
    sample=$(basename $logfile)
    sample="${sample%.*}"
    reads=$(grep "reads;" $logfile | cut -f 1 -d " ")
    uniq_pair=$(grep "concordantly exactly" $logfile | awk '{print $1}')
    pct_uniq_pair=$(grep "concordantly exactly" $logfile | grep -o -E "\w+\.\w+")
    multiple_pair=$(grep "concordantly >1" $logfile | awk '{print $1}')
    pct_multiple_pair=$(grep "concordantly >1" $logfile | grep -o -E "\w+\.\w+")
    not_map_pair=$(grep "concordantly 0 times$" $logfile | awk '{print $1}')
    pct_not_map_pair=$(grep "concordantly 0 times$" $logfile | grep -o -E "\w+\.\w+")
    discordant_pair=$(grep "discordantly 1" $logfile | awk '{print $1}')
    pct_discordant_pair=$(grep "discordantly 1" $logfile | grep -o -E "\w+\.\w+")
    reads_unmapped_pairs=$(grep "mates" $logfile | awk '{print $1}')
    unique=$(grep "aligned exactly 1 time" $logfile | awk '{print $1}' )
    pct_unique=$(grep "aligned exactly 1 time" $logfile | grep -o -E "\w+\.\w+")
    multiple=$(grep "aligned >1 times" $logfile | awk '{print $1}')
    pct_multiple=$(grep "aligned >1 times" $logfile | grep -o -E "\w+\.\w+")
    not_mapped=$(grep "aligned 0 times$" $logfile | awk '{print $1}')
    pct_not_mapped=$(grep "aligned 0 times$" $logfile | grep -o -E "\w+\.\w+")
    total=$(grep "overall alignment rate" $logfile | grep -o -E "\w+\.\w+")
    echo -e "$sample\t$reads\t$uniq_pair\t$pct_uniq_pair\t$multiple_pair\t$pct_multiple_pair\t$not_map_pair\t$pct_not_map_pair\t$discordant_pair\t$pct_discordant_pair\t$reads_unmapped_pairs\t$unique\t$pct_unique\t$multiple\t$pct_multiple\t$not_mapped\t$pct_not_mapped\t$total"
done >> $2

