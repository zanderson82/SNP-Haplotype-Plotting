#!/bin/bash
#module load samtools
mkdir -p /n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card
dir="/n/1000g/align-card-2.24-hg38/FIRST_100/"
for bam in $dir/*/*PMDV_FINAL.haplotagged.bam; do
    sample_name=$(basename "$bam" | cut -d '-' -f 1)
    mkdir -p "/n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card/${sample_name}"
    output_bam="/n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card/${sample_name}/${sample_name}_OPSIN_BAM_card.bam"
    echo "Output bam: $output_bam"
    samtools view -b -o "$output_bam" "$bam" chrX:154140724-154196769 
    ls $output_bam
    if [ $? -eq 0 ] ; then
        echo "samtools view completed successfully, proceeding with samtools sort"

#sort and index the extracted bam file
sorted_output_bam="/n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card/${sample_name}/${sample_name}_sorted_OPSIN_BAM_card.bam"
        
        samtools sort -o "$sorted_output_bam" "$output_bam"
        if [ $? -eq 0 ]; then
            echo "samtools sort completed successfully, proceeding with samtools index..."
            samtools index "$sorted_output_bam"
        else
            echo "samtools sort failed."
        fi
    else
        echo "samtools view failed."
    fi
done
#coords for MECP2 to OPN1MW : chrX:154086114-154197630
#coords for OPN1LW to TEX28: chrX:154143264-154292599
#Filter reads that span the coordinates
for sorted_bams in /n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card/*/*_sorted_OPSIN_BAM_card.bam; do
    sample_name=$(basename "$sorted_bams" "_sorted_OPSIN_BAM_card.bam")
read_filtered_bam="/n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card/${sample_name}/${sample_name}_read_filtered_OPSIN_BAM_card.bam"
    samtools view -h "$sorted_bams" | \
    awk '{ 
        if ($1 ~ /^@/) {print $0; next}
        chr="chrX";
        pos1=154140724;
        pos2=154196769;
        split($0, a, "\t");
        start=a[4];
        len=length(a[10]);
        end=start+len;
        if ((start <= pos1 && end >= pos2)) {
            print $0
        }
    }' | \
    samtools view -Sb - > $read_filtered_bam

done
for read_filtered_bams in /n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card/*/*_read_filtered_OPSIN_BAM_card.bam; do
    sample_name=$(basename "$read_filtered_bams" "_read_filtered_OPSIN_BAM_card.bam" )
filtered_sorted_bam="/n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card/${sample_name}/${sample_name}_promoter_to_OPN1MW_filtered_sorted_card.bam"
    samtools sort -o "$filtered_sorted_bam" "$read_filtered_bams"
    samtools index "$filtered_sorted_bam"
done


