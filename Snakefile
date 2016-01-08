# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from collections import defaultdict
file_info = defaultdict(list)

paired_end = False

with open( "metasheet.csv", "r" ) as meta_fh:
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
    trim_out_files = []
    for sample in file_info.keys():
        if( paired_end ):
            trim_out_files.append( "analysis/trimmomatic/" + sample + "/" + sample + ".left.paired.trim.fastq.gz" )
            trim_out_files.append( "analysis/trimmomatic/" + sample + "/" + sample + ".right.paired.trim.fastq.gz" )
        else:
            trim_out_files.append( "analysis/trimmomatic/" + sample + "/" + sample + ".single.trim.fastq.gz" )
    
    return trim_out_files

def trim_report( wildcards ):
    if(paired_end):
        return "analysis/trimmomatic/trim_pe_report.tab"
    else:
        return "analysis/trimmomatic/trim_se_report.tab"

rule target:
    input:
        trim_output,
        expand( "analysis/trimmomatic/{sample}/{sample}.trim.log", sample=file_info.keys() ),
        trim_report

rule run_trim_pe:
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
       
rule run_trim_se:
    input:
        get_fastq
    output:
        "analysis/trimmomatic/{sample}/{sample}.single.trim.fastq.gz"
    log:
        "analysis/trimmomatic/{sample}/{sample}.trim.log"
    params:
        SE_adapter="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TruSeq3-SE.fa"
    threads: 4
    shell:
        "java -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/trimmomatic-0.32.jar "
        " SE -threads {threads} {input} {output}"
        " ILLUMINACLIP:{params.SE_adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 >&{log}"

rule trim_report_pe:
    input:
        trim_log_files = expand( "analysis/trimmomatic/{sample}/{sample}.trim.log", sample=file_info.keys() )
    output:
        trim_report="analysis/trimmomatic/trim_pe_report.tab",
        trim_plot="analysis/trimmomatic/trim_pe_report.png"
    run:
        log_file_list = " -l ".join( input.trim_log_files )
        shell( "perl trim_and_align/scripts/trim_report_pe.pl -f {log_file_list} 1>{output.trim_report}" )
        shell( "Rscript trim_and_align/scripts/trim_plot_pe.R {output.trim_report} {output.trim_plot}" )

rule trim_report_se:
    input:
        trim_log_files = expand( "analysis/trimmomatic/{sample}/{sample}.trim.log", sample=file_info.keys() )
    output:
        trim_report="analysis/trimmomatic/trim_se_report.tab",
        trim_plot="analysis/trimmomatic/trim_se_report.png"
    run:
        log_file_list = " -l ".join( input.trim_log_files )
        shell( "perl trim_and_align/scripts/trim_report_se.pl -f {log_file_list} 1>{output.trim_report}" )
        shell( "Rscript trim_and_align/scripts/trim_plot_se.R {output.trim_report} {output.trim_plot}" )






