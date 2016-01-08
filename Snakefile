# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from collections import defaultdict
file_info = defaultdict(list)

paired_end = False

with open( "meta.csv", "r" ) as meta_fh:
    next(meta_fh)
    for line in meta_fh:
        info = line.strip().split(",")
        right_mate = info[0].replace("_R1_", "_R2_")
        if os.path.isfile("./concat_per_sample_fastq/" + right_mate):
            file_info[info[1]] = [info[0], right_mate]
            paired_end = True
        else:
            file_info[info[1]] = [info[0]]
            

def get_fastq( wildcards ):
    return [ os.path.join( "concat_per_sample_fastq", f ) for f in file_info[wildcards.sample] ]

def trim_output( wildcards ):
    

rule target:
    input:
        expand( "analysis/trimmomatic/{sample}/{sample}.left.paired.trim.fastq.gz", sample=file_info.keys() ),
        expand( "analysis/trimmomatic/{sample}/{sample}.right.paired.trim.fastq.gz", sample=file_info.keys() )

rule run_trim:
    input:
        get_fastq
    output:
        left_paired_trim="analysis/trimmomatic/{sample}/{sample}.left.paired.trim.fastq.gz",
        right_paired_trim="analysis/trimmomatic/{sample}/{sample}.right.paired.trim.fastq.gz",
        left_unpaired_trim="analysis/trimmomatic/{sample}/{sample}.left.unpaired.trim.fastq.gz",
        right_unpaired_trim="analysis/trimmomatic/{sample}/{sample}.right.unpaired.trim.fastq.gz"
    log:
        "analysis/trimmomatic/{sample}/{sample}.trim.log"
    params:
        PE_adapter="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TruSeq3-PE.fa"
    threads: 4
    shell:
        "java -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/trimmomatic-0.32.jar "
        " PE -threads {threads} {input} {output.left_paired_trim} {output.left_unpaired_trim} {output.right_paired_trim} {output.right_unpaired_trim}"
        " ILLUMINACLIP:{params.PE_adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 >&{log}" 
       
 
