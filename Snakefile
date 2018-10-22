"""
Author: E. Frouin
Affiliation: MIO
Aim: A MetaG pipeline.
Date: Wen Jul 13 14:35 2016
Run: snakemake
"""
from os.path import join
configfile: "config.yaml"

WD= config["wdir"]
SAMPLES, = glob_wildcards(join(WD, '{sample,[^/]+}_1.fastq.gz'))
CPU= config["threads"]

rule all:
    input:
        end1= expand("kegg50_2ndAlign/{sample}_reads_perko.csv",sample=SAMPLES),
        end2="norm_coverage/CovCOG_summary.tsv",
        end3="norm_naive/COG_count.csv",
        end4=expand("prokka/{sample}_ORFs.faa",sample=SAMPLES),
        end5=expand("diamond/{sample}_prot.daa",sample=SAMPLES)


rule pigz :
    input:
        R1= join(WD,'{sample}_1.fastq.gz'),
        R2= join(WD,'{sample}_2.fastq.gz')
    output:
        R1 = temp("deduplicate/{sample}_R1.fastq"),
        R2 = temp("deduplicate/{sample}_R2.fastq")
    resources: ram=4
    shell:
        "pigz -cd {input.R1} > {output.R1}; "
        "pigz -cd {input.R2} > {output.R2} "


rule deduplicate :
    input :
        R1= "deduplicate/{sample}_R1.fastq",
        R2= "deduplicate/{sample}_R2.fastq"
    output:
        R1_dedup="deduplicate/{sample}_R1_dedup.fastq",
        R2_dedup="deduplicate/{sample}_R2_dedup.fastq"
    resources: ram=10
    shell:
        "cd-hit-dup -i {input.R1} -i2 {input.R2} "
        "-o {output.R1_dedup} -o2 {output.R2_dedup} > /dev/null 2>&1"

rule merging:
    input:
        R1= "deduplicate/{sample}_R1_dedup.fastq",
        R2= "deduplicate/{sample}_R2_dedup.fastq"
    output:
        merged=temp("merging/{sample}_merged.fastq"),
        unmerged="merging/{sample}_unaligned.txt"
    threads: CPU
    log:
        "logs/merging/{sample}.log.bz2"
    shell:
        "pandaseq -f  {input.R1}  -r  {input.R2}  -F -w {output.merged} "
        "-U {output.unmerged} -T {threads} -G {log};"



rule trimming :
    input :
        "merging/{sample}_merged.fastq"
    params:
        adapt =config["T_ILLUMINA_ADAPTORS"],
        min_len= config['T_MIN_LEN'],
        score_win4=config['T_QUAL_SCORE_WIN4']
    threads: CPU
    log:
        "logs/trimming/{sample}.log"
    output:
        "trimming/{sample}_trim.fastq.gz"
    shell:
        "trimmomatic SE -threads {threads} -phred33 {input} "
        "{output} ILLUMINACLIP:{params.adapt}:2:30:10 LEADING:20 TRAILING:20 "
        "SLIDINGWINDOW:4:{params.score_win4} MINLEN:{params.min_len}>{log} 2>&1"


rule gunzip :
    input:
        "trimming/{sample}_trim.fastq.gz"
    output:
        temp("trimming/{sample}_trim.fastq")
    resources: ram=4
    shell:
        "gzip -dk {input}"

rule fq2fa:
    input:
        "trimming/{sample}_trim.fastq"
    output:
        fa=temp("trimming/{sample}_trim.fasta")
    shell :
        " awk '{{print \">\" substr($0,2);getline;print;getline;getline}}' "
        "{input} > {output} "


rule assembly :
    input:
        "trimming/{sample}_trim.fasta"
    output:
        dir= "assembly/{sample}_dir",
        conigs = "assembly/{sample}_dir/contig.fa"
    threads: CPU
    log:
        "logs/assembly/{sample}.log"
    shell:
        " idba_ud -r {input}  -o {output.dir} --num_threads {threads} "
        "> {log} 2>&1 ||true"

rule index_contig :
    input:
        contigs="assembly/{sample}_dir/contig.fa"
    output:
        expand("assembly/{{sample}}_dir/contig.fa.{extension}", extension=['amb','ann','bwt','pac','sa'])
    log:
        "logs/index_contigs/{sample}.log"
    shell:
        "bwa index {input.contigs}"
        #>{log} 2>&1"

rule bwa_map :
    input:
        reads="trimming/{sample}_trim.fasta",
        contigs="assembly/{sample}_dir/contig.fa",
        index= expand("assembly/{{sample}}_dir/contig.fa.{extension}", extension=['amb','ann','bwt','pac','sa'])
    output:
        temp("bwa_map/{sample}.bam")
    log:
        "logs/bwa_map/{sample}.log"
    threads:CPU
    shell:
        "(bwa mem -t {threads} {input.contigs}  {input.reads}"
        "| samtools view -Sb -F 2308 - >{output}) 2> {log}"


rule samtools_sort:
    input:
        "bwa_map/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"



rule annotation :
    input :
        contigs="assembly/{sample}_dir/contig.fa"
    output:
        faa="annotation/{sample}_ORFs.faa",
        gff="annotation/{sample}_ORFs.gff"
    log:
        "logs/annotation/{sample}.log"
    shell :
        "prodigal -a {output.faa} -i {input.contigs} -f gff "
        "-p meta > {output.gff} 2> {log}"


rule prokka :
    input :
        contigs="assembly/{sample}_dir/contig.fa"
    output:
        faa="prokka/{sample}_ORFs.faa",
        gff="prokka/{sample}_ORFs.gff"
    log:
        "logs/prokka/{sample}.log"
    threads: CPU
    params:
        dir_inter="prokka/{sample}_dir",
        dir_glob="prokka",
        name="{sample}_ORFs"
    shell :
        "prokka --metagenome --cpus {threads}  --outdir {params.dir_inter} "
        "--centre C --locustag c {input} 2> {log} ; "
        "mv {params.dir_inter}/*.faa {params.dir_glob}/{params.name}.faa; "
        "mv {params.dir_inter}/*.gff {params.dir_glob}/{params.name}.gff; "
        "rm -r {params.dir_inter} "

rule diamond :
    input :
        contigs="annotation/{sample}_ORFs.faa"
    output:
        "diamond/{sample}_prot.daa"
    threads: CPU
    params:
        path_db= config["D_PATH_DB"],
        evalue= config["D_EVALUE"]
    shell:
        "diamond blastp --query {input} --db {params.path_db} --daa {output} "
        " -p {threads} -e {params.evalue}"


rule COGing :
    input:
        "annotation/{sample}_ORFs.faa"
    output:
        "COGing/{sample}_blast.xml"
    threads: CPU
    log:
        "logs/COGing/{sample}.log"
    params:
        path_db= config["C_PATH_DB"],
        evalue= config["C_EVALUE"]
    shell :
        "rpsblast -query {input} -db {params.path_db} -evalue {params.evalue} "
        "-outfmt 5 -out {output} -num_threads {threads} > {log} "

rule KEGGing :
    input:
        "annotation/{sample}_ORFs.faa"
    output:
        "blastp_keeg/{sample}_kegg.txt"
    threads: CPU
    log:
        "logs/KEGG/{sample}.log"
    params:
        path_db= config["K_PATH_DB"],
        evalue= config["K_EVALUE"]
    shell :
        " diamond blastp --threads {threads} --query {input} --out {output} "
        " --outfmt 6 qseqid stitle --db {params.path_db} -k 1 "
        " -e {params.evalue} > {log} "

rule sortKegg :
    input:
        "blastp_keeg/{sample}_kegg.txt"
    output:
        "kegg50_2ndAlign/{sample}_Onlyko.txt"
    shell :
        " grep -P 'K[0-9]{{5}}' {input} | sed  -E "
        " 's/(contig-100_[0-9]*_[0-9]*)\t.*(K[0-9]{{5}}).*/\\1\t\\2/' > {output} "


rule parseXML:
    input :
        "COGing/{sample}_blast.xml"
    output :
        COG="norm_naive/{sample}_COG.csv",
        cat="norm_naive/{sample}_cat.csv"
    params:
            dbCOG=config["S_DBCOG"]
    shell:
        "Parse_RpsBlast.py {input} -COGfunction {output.COG} "
        "-COGcategory {output.cat} -COG_db {params.dbCOG}"


rule rawCOGtable :
    input:
        COG=expand("norm_naive/{sample}_COG.csv",sample=SAMPLES),
        cat=expand("norm_naive/{sample}_cat.csv",sample=SAMPLES)
    output :
        tabCOG="norm_naive/COG_count.csv"
        # tabcat="norm_naive/cat_count.csv"
    params:
        listCOG= config["S_LISTCOG"]
    shell:
        "cp {params.listCOG} {output.tabCOG};"
        "for filename in norm_naive/*_COG.csv; do "
        "join {output.tabCOG} \"${{filename}}\" -t $\'\t\' > norm_naive/temp;"
        "rm -f {output.tabCOG};"
        "mv norm_naive/temp {output.tabCOG};"
        "done"


rule XML2WMGA:
    input :
        "COGing/{sample}_blast.xml"
    output :
        "norm_coverage/{sample}_WMGA.csv"
    shell:
        "rpsBlast_to_WMGA.py -output {output} {input}"

rule compute_cov:
    input:
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai",
        gff="annotation/{sample}_ORFs.gff"
    resources: ram=4
    output:
        "norm_coverage/{sample}_covORFs.hist"
    shell:
        "bedtools coverage -hist -a {input.gff} -b {input.bam} > {output}"



rule coverage_perCOG:
    input:
        faa="annotation/{sample}_ORFs.faa",
        gff="annotation/{sample}_ORFs.gff",
        WMGA="norm_coverage/{sample}_WMGA.csv",
        hist="norm_coverage/{sample}_covORFs.hist"
    output:
        "norm_coverage/{sample}_covCOG.tsv"
    resources: ram=8
    params:
        name="{sample}"
    shell:
        "cov_per_cog.py {input.gff} {input.faa} {input.WMGA} {input.hist} "
        "--samplenames {params.name} > {output} "


rule coverage_perKEGG:
    input:
        ko="kegg50_2ndAlign/{sample}_Onlyko.txt",
        hist="norm_coverage/{sample}_covORFs.hist"
    output:
        "kegg50_2ndAlign/{sample}_reads_perko.csv"
    resources: ram=5
    params:
        name="{sample}"
    shell:
        "coverage_per_Ko.py --sampleK {input.ko} --sampleHist {input.hist} "
	    " --sampleName {params.name} --out {output} "

rule template_COG:
    input :
        expand("norm_coverage/{sample}_covCOG.tsv", sample=SAMPLES)
    output :
        "norm_coverage/CovCOG_summary.tsv"
    params:
        header='cog_hit\tcog_description\tcog_class\tcog_class_description'
    shell:
        "cat {input} |grep COG |sed 's/\t/::/g' | awk '{{$NF=\"\"; print $0}}'"
        " |  sed 's/::/\t/g' | sort |uniq  > {output};"
        "sed -i '1i{params.header}' {output};"
        "for filename in norm_coverage/*_covCOG.tsv; do "
        "sed 's/\t.*\t/\t/' \"${{filename}}\" >\"${{filename}}\"_temp;"
        "join {output} \"${{filename}}\"_temp -t $\'\t\' > norm_coverage/temp;"
        "rm -f {output} \"${{filename}}\"_temp;"
        "mv norm_coverage/temp {output};"
        "done"
