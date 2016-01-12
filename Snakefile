# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
#
import os
from collections import defaultdict
file_info = defaultdict(list)

configfile: "config.yaml"
strand_command=""
paired_end = False

if( config["stranded"] ):
    strand_command="--outFilterIntronMotifs RemoveNoncanonical"
else:
    strand_command="--outSAMstrandField intronMotif"

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

def get_trim_fastq( wildcards ):
    if( paired_end ):
        return [ "analysis/trimmomatic/" + wildcards.sample + "/" + wildcards.sample + ".left.paired.trim.fastq.gz",
        "analysis/trimmomatic/" + wildcards.sample + "/" + wildcards.sample + ".right.paired.trim.fastq.gz" ]
    else:
        return [ "analysis/trimmomatic/" + wildcards.sample + "/" + wildcards.sample + ".single.trim.fastq.gz" ]

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
        trim_report,
        "analysis/STAR/STAR_Align_Report.csv",
        "analysis/STAR/STAR_Align_Report.png",
        "analysis/STAR/STAR_Gene_Counts.csv"

rule run_trim_pe:
    input:
        get_fastq
    output:
        left_paired_trim="analysis/trimmomatic/{sample}/{sample}.left.paired.trim.fastq.gz",
        right_paired_trim="analysis/trimmomatic/{sample}/{sample}.right.paired.trim.fastq.gz",
        left_unpaired_trim="analysis/trimmomatic/{sample}/{sample}.left.unpaired.trim.fastq.gz",
        right_unpaired_trim="analysis/trimmomatic/{sample}/{sample}.right.unpaired.trim.fastq.gz",
        trim_log="analysis/trimmomatic/{sample}/{sample}.trim.pe.log"
    params:
        PE_adapter="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TruSeq3-PE.fa"
    threads: 4
    shell:
        "java -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/trimmomatic-0.32.jar "
        " PE -threads {threads} {input} {output.left_paired_trim} {output.left_unpaired_trim} {output.right_paired_trim} {output.right_unpaired_trim}"
        " ILLUMINACLIP:{params.PE_adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 >&{output.trim_log}" 
       
rule run_trim_se:
    input:
        get_fastq
    output:
        "analysis/trimmomatic/{sample}/{sample}.single.trim.fastq.gz",
        "analysis/trimmomatic/{sample}/{sample}.trim.se.log"
    params:
        SE_adapter="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TruSeq3-SE.fa"
    threads: 4
    shell:
        "java -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/trimmomatic-0.32.jar "
        " SE -threads {threads} {input} {output[0]}"
        " ILLUMINACLIP:{params.SE_adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 >&{output[1]}"

rule trim_report_pe:
    input:
        trim_log_files = expand( "analysis/trimmomatic/{sample}/{sample}.trim.pe.log", sample=file_info.keys() )
    output:
        trim_report="analysis/trimmomatic/trim_pe_report.tab",
        trim_plot="analysis/trimmomatic/trim_pe_report.png"
    run:
        log_file_list = " -l ".join( input.trim_log_files )
        shell( "perl trim_and_align/scripts/trim_report_pe.pl -f {log_file_list} 1>{output.trim_report}" )
        shell( "Rscript trim_and_align/scripts/trim_plot_pe.R {output.trim_report} {output.trim_plot}" )

rule trim_report_se:
    input:
        trim_log_files = expand( "analysis/trimmomatic/{sample}/{sample}.trim.se.log", sample=file_info.keys() )
    output:
        trim_report="analysis/trimmomatic/trim_se_report.tab",
        trim_plot="analysis/trimmomatic/trim_se_report.png"
    run:
        log_file_list = " -l ".join( input.trim_log_files )
        shell( "perl trim_and_align/scripts/trim_report_se.pl -f {log_file_list} 1>{output.trim_report}" )
        shell( "Rscript trim_and_align/scripts/trim_plot_se.R {output.trim_report} {output.trim_plot}" )


rule run_STAR:
    input:
        get_trim_fastq
    output:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        counts="analysis/STAR/{sample}/{sample}.counts.tab",
        log_file="analysis/STAR/{sample}/{sample}.Log.final.out"
    params:
        stranded=strand_command,
        prefix=lambda wildcards: "analysis/STAR/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    threads: 8
    shell:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR/STAR-STAR_2.4.2a/source/STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[star_index]}"
        " --readFilesIn {input} --readFilesCommand zcat --outFileNamePrefix {params.prefix}."
        "  --outSAMmode Full --outSAMattributes All {params.stranded} --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000 --quantMode GeneCounts"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"
        " && /usr/bin/samtools index {output.bam}"

rule generate_STAR_report:
    input:
        star_log_files=expand( "analysis/STAR/{sample}/{sample}.Log.final.out", sample=file_info.keys() ),
        star_gene_count_files=expand( "analysis/STAR/{sample}/{sample}.counts.tab", sample=file_info.keys() )
    output:
        csv="analysis/STAR/STAR_Align_Report.csv",
        png="analysis/STAR/STAR_Align_Report.png",
        gene_counts="analysis/STAR/STAR_Gene_Counts.csv"
    run:
        log_files = " -l ".join( input.star_log_files )
        count_files = " -f ".join( input.star_gene_count_files )
        shell( "perl trim_and_align/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript trim_and_align/scripts/map_stats.R {output.csv} {output.png}" )
        shell( "perl trim_and_align/scripts/raw_and_fpkm_count_matrix.pl -f {count_files} 1>{output.gene_counts}" )

rule copy_output:
    input:
        trim_dir="analysis/trimmomatic/",
        star_dir="analysis/STAR/"
    output:
        trim_dir="analysis/final_output/trimmed_fastq/",
        align_dir="analysis/final_output/alignment/",
        sum_dir="analysis/final_output/summary/"
    shell:
        "cp -rf {input.trim_dir}/* {output.trim_dir}/"
        " && find {input.star_dir}/ -type f -name '*sorted.ba*' -exec cp -ft {output.align_dir}/ {{}} \;"
        " && cp {input.trim_dir}/trim_*.png {output.sum_dir}/"
        " && cp {input.trim_dir}/trim_*.tab {output.sum_dir}/"
        " && cp {input.star_dir}/STAR_Align_Report.* {output.sum_dir}/"
