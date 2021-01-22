## Global Variables
vThreads = 40

## Pipeline, outputs
vOutdir              = "/projects/bgmp/shared/groups/2020/algae/personal_jonas/snakemake_mini/output/"
vOutdir_fastp        = vOutdir + "fastp/"
vOutdir_samtools_1   = vOutdir + "samtools/1_oceana/"
vOutdir_samtools_2   = vOutdir + "samtools/2_salina/"
vOutdir_metaspades   = vOutdir + "metaspades/"
vOutdir_quast        = vOutdir + "quast/"
vOutdir_maxbin2      = vOutdir + "maxbin2/"
vOutdir_mapping      = vOutdir + "mapping/"
vOutdir_metabat2     = vOutdir + "metabat2/"
vOutdir_concoct      = vOutdir + "concoct/"
vOutdir_bin3c        = vOutdir + "bin3c/"
vOutdir_dastool      = vOutdir + "dastool/"

## Pipeline, inputs
vInput_fastp      = "/projects/bgmp/shared/groups/2020/algae/personal_jonas/data_raw_mini/"

## Pipeline, Resources
vRsrc_ref_oceana = "/projects/bgmp/shared/groups/2020/algae/personal_jonas/nano_chlor_ref/ncbi_dataset_oceana/data/GCA_004519485.1"
vRsrc_ref_salina = "/projects/bgmp/shared/groups/2020/algae/personal_jonas/nano_chlor_ref/ncbi_dataset_salina/data/GCA_004565275.1"



## Input file organization, fastp
vPond = ["algae52", "algae53", "algae114", "algae115"]
vType = ["hic", "sho"]
vRead = ["R1.fastq.gz", "R2.fastq.gz"]

## Input file organization, bwa
vIndex = [".amb", ".ann", ".bwt", ".pac", ".sa"]

## Initialize outputs
rule all:
    """
    > Function: Define all target output files
    """
    input:
        ## fastp
        expand(vOutdir_fastp + "{pond}/{type}/{pond}_{type}_{read}", pond=vPond, type=vType, read=vRead),
        ## bwa (reference index)
        expand(vRsrc_ref_oceana + "ref_seqs.fa{index}", index=vIndex),   ## index, oceana
        expand(vRsrc_ref_salina + "ref_seqs.fa{index}", index=vIndex),   ## index, salina
        ## samtools
        expand(vOutdir_samtools_1 + "{pond}/{type}/unmapped_{pond}_{type}_{read}", pond=vPond, type=vType, read=vRead), ## N. oceana
        expand(vOutdir_samtools_2 + "{pond}/{type}/unmapped_{pond}_{type}_{read}", pond=vPond, type=vType, read=vRead), ## N. salina
        ## spades
        # expand(vOutdir_metaspades + "{pond}_sho/", pond=vPond),
        ## QUAST
        expand(vOutdir_quast + "{pond}_sho", pond=vPond),
        ## MaxBin2
        expand(vOutdir_maxbin2 + "{pond}_sho/", pond=vPond),
        ## Mapping (post-assembly mapping files: index, alignment, depth)
        expand(vOutdir_mapping + "index/{pond}_sho/scaffolds.fasta", pond=vPond),
        expand(vOutdir_mapping + "alignment/{pond}_{type}.bam", pond=vPond, type=vType),
        expand(vOutdir_mapping + "depth/{pond}_depth.txt", pond=vPond),
        ## MetaBat2
        expand(vOutdir_metabat2 + "{pond}_sho/", pond=vPond),
        ## CONCOCT
        # expand(vOutdir_concoct + "{pond}_sho/", pond=vPond),
        ## bin3c
        expand(vOutdir_bin3c + "{pond}_sho/", pond=vPond)





## 01 - Read Processing, fastp
rule fastp:
    """
    > Function: Removes adapters and low-quality reads.
    """
    input:
        r1 = vInput_fastp + "{pond}/{type}/{pond}_{type}_R1.fastq.gz",      ## R1 Input
        r2 = vInput_fastp + "{pond}/{type}/{pond}_{type}_R2.fastq.gz"       ## R2 Input
    output:
        r1 = vOutdir_fastp + "{pond}/{type}/{pond}_{type}_R1.fastq.gz",      ## R1 Output
        r2 = vOutdir_fastp + "{pond}/{type}/{pond}_{type}_R2.fastq.gz",      ## R2 Output
        vH = vOutdir_fastp + "summary/{pond}_{type}_summary.html",           ## HTML Summary
        vJ = vOutdir_fastp + "summary/{pond}_{type}_summary.json"            ## JSON Summary
    conda:
        "envs/env_1.yaml"
    shell:
        "/usr/bin/time -v fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.vH} -j {output.vJ}"



## 01 - Read Processing, bwa Index
rule bwa_index:
    """
    > Before BWA can use the reference sequence to align reads, it must generate index files.
    > Index files are generated based on the reference sequence file.
    """
    input:
        vRsrc_ref_salina = vRsrc_ref_oceana + "ref_seqs.fa",
        vRsrc_ref_oceana = vRsrc_ref_salina + "ref_seqs.fa"
    output:
        expand(vRsrc_ref_oceana + "ref_seqs.fa{index}", index=vIndex),   ## index, oceana
        expand(vRsrc_ref_salina + "ref_seqs.fa{index}", index=vIndex),   ## index, salina
    conda:
        "envs/env_1.yaml"
    threads:
        vThreads
    shell:
        """
        bwa index {input.vRsrc_ref_salina}
        bwa index {input.vRsrc_ref_oceana}
        """



## 01 - Read Processing, bwa mem & samtools (2X filtering)
rule bwa_samtools_1:
    """
    > Function: Align reads to reference genomes N. oceana and N. salina then extract unmapped reads with samtools
    """
    input:
        ## Reference Genomes
        ref_oceana = vRsrc_ref_oceana + "ref_seqs.fa",
        ## Shotgun
        sho_r1 = vOutdir_fastp + "{pond}/sho/{pond}_sho_R1.fastq.gz",
        sho_r2 = vOutdir_fastp + "{pond}/sho/{pond}_sho_R2.fastq.gz",
        ## Hi-C
        hic_r1 = vOutdir_fastp + "{pond}/hic/{pond}_hic_R1.fastq.gz",
        hic_r2 = vOutdir_fastp + "{pond}/hic/{pond}_hic_R2.fastq.gz"
    output:
        ## Shotgun
        sho_r1 = vOutdir_samtools_1 + "{pond}/sho/unmapped_{pond}_sho_R1.fastq.gz",
        sho_r2 = vOutdir_samtools_1 + "{pond}/sho/unmapped_{pond}_sho_R2.fastq.gz",
        ## Hi-C
        hic_r1 = vOutdir_samtools_1 + "{pond}/hic/unmapped_{pond}_hic_R1.fastq.gz",
        hic_r2 = vOutdir_samtools_1 + "{pond}/hic/unmapped_{pond}_hic_R2.fastq.gz"
    conda:
        "envs/env_1.yaml"
    threads:
        vThreads
    shell:
        """
        /usr/bin/time -v bwa mem -t {threads} {input.ref_oceana} {input.sho_r1} {input.sho_r2} | \
        /usr/bin/time -v samtools view -@ {threads} -S -b - | \
        /usr/bin/time -v samtools view -@ {threads} -b -f 4 - | \
        /usr/bin/time -v samtools sort -@ {threads} -n - | \
        /usr/bin/time -v samtools fastq -@ {threads} -1 {output.sho_r1} -2 {output.sho_r2} -0 /dev/null -s /dev/null -n -
        
        /usr/bin/time -v bwa mem -5SP -t {threads} {input.ref_oceana} {input.hic_r1} {input.hic_r2} | \
        /usr/bin/time -v samtools view -@ {threads} -S -b - | \
        /usr/bin/time -v samtools view -@ {threads} -b -f 4 - | \
        /usr/bin/time -v samtools sort -@ {threads} -n - | \
        /usr/bin/time -v samtools fastq -@ {threads} -1 {output.hic_r1} -2 {output.hic_r2} -0 /dev/null -s /dev/null -n -
        """



## 01 - Read Processing, bwa mem & samtools (2X filtering)
rule bwa_samtools_2:
    """
    > Function: Align reads to reference genomes N. oceana and N. salina then extract unmapped reads with samtools
    """
    input:
        ## Reference Genomes
        ref_salina = vRsrc_ref_salina + "ref_seqs.fa",
        ## Shotgun
        sho_r1 = vOutdir_samtools_1 + "{pond}/sho/unmapped_{pond}_sho_R1.fastq.gz",
        sho_r2 = vOutdir_samtools_1 + "{pond}/sho/unmapped_{pond}_sho_R2.fastq.gz",
        ## Hi-C
        hic_r1 = vOutdir_samtools_1 + "{pond}/hic/unmapped_{pond}_hic_R1.fastq.gz",
        hic_r2 = vOutdir_samtools_1 + "{pond}/hic/unmapped_{pond}_hic_R2.fastq.gz"
    output:
        ## Shotgun
        sho_r1 = vOutdir_samtools_2 + "{pond}/sho/unmapped_{pond}_sho_R1.fastq.gz",
        sho_r2 = vOutdir_samtools_2 + "{pond}/sho/unmapped_{pond}_sho_R2.fastq.gz",
        ## Hi-C
        hic_r1 = vOutdir_samtools_2 + "{pond}/hic/unmapped_{pond}_hic_R1.fastq.gz",
        hic_r2 = vOutdir_samtools_2 + "{pond}/hic/unmapped_{pond}_hic_R2.fastq.gz"
    conda:
        "envs/env_1.yaml"
    threads:
        vThreads
    shell:
        """
        /usr/bin/time -v bwa mem -t {threads} {input.ref_salina} {input.sho_r1} {input.sho_r2} | \
        /usr/bin/time -v samtools view -@ {threads} -S -b - | \
        /usr/bin/time -v samtools view -@ {threads} -b -f 4 - | \
        /usr/bin/time -v samtools sort -@ {threads} -n - | \
        /usr/bin/time -v samtools fastq -@ {threads} -1 {output.sho_r1} -2 {output.sho_r2} -0 /dev/null -s /dev/null -n -

        /usr/bin/time -v bwa mem -5SP -t {threads} {input.ref_salina} {input.hic_r1} {input.hic_r2} | \
        /usr/bin/time -v samtools view -@ {threads} -S -b - | \
        /usr/bin/time -v samtools view -@ {threads} -b -f 4 - | \
        /usr/bin/time -v samtools sort -@ {threads} -n - | \
        /usr/bin/time -v samtools fastq -@ {threads} -1 {output.hic_r1} -2 {output.hic_r2} -0 /dev/null -s /dev/null -n -
        """



# 02 - Assembly, metaspades
rule metaspades:
    """
    Function:   Produce metagenome assemblies of all shotgun populations
    Parameters:
        > --meta        Indicates that spades should be ran as metaspades
        > -k            Specify kmer lengths to perform assembly with
        > -1            Input Read 1
        > -2            Input Read 2
        > --threads     Number of threads to use
        > -o            Output Directory
    Input:      Nannochloropsis-free reads (.fastq.gz)
    Output:     Shotgun assembly (format = scaffolds.fasta)
    """
    input:
        r1 = vOutdir_samtools_2 + "{pond}/sho/unmapped_{pond}_sho_R1.fastq.gz",
        r2 = vOutdir_samtools_2 + "{pond}/sho/unmapped_{pond}_sho_R2.fastq.gz"
    output:
        vOutdir_metaspades + "{pond}_sho/"
    conda:
        "envs/env_1.yaml"
    threads:
        vThreads
    shell:
        """
        /usr/bin/time -v spades.py --meta \
        -k 21,33,55,77,99,127 \
        --threads {threads} \
        -1 {input.r1} \
        -2 {input.r2} \
        -o {output}
        """



## 02 - Assembly, QUAST
rule quast:
    input:
        vOutdir_metaspades + "{pond}_sho/scaffolds.fasta"
    output:
        vOutdir_quast + "{pond}_sho"
    conda:
        "envs/env_2.yaml"
    threads:
        vThreads
    shell:
        """
        /usr/bin/time -v metaquast \
        --threads {threads} \
        -o {output} \
        {input}
        """



## 03 - Binning, MaxBin2
rule maxbin2:
    """
    Function:   Take MAG and group genomes using _.
    Resources:
        > https://sourceforge.net/projects/maxbin2/files/
    """
    input:
        contigs = vOutdir_metaspades + "{pond}_sho/scaffolds.fasta",
        r1      = vOutdir_metaspades + "{pond}_sho/corrected/unmapped_s_{pond}_sho_R1.fastq.00.0_0.cor.fastq.gz",
        r2      = vOutdir_metaspades + "{pond}_sho/corrected/unmapped_s_{pond}_sho_R2.fastq.00.0_0.cor.fastq.gz"
    output:
        directory(vOutdir_maxbin2 + "{pond}_sho/")
    conda:
        "envs/env_1.yaml"
    threads:
        vThreads
    shell:
        """
        /usr/bin/time -v run_MaxBin.pl \
        -contig {input.contigs} \
        -out {output} \
        -reads {input.r1} \
        -reads2 {input.r2} \
        -thread {threads}
        """



# 03 - Binning, Generate alignment (used by MetaBat2 and CONCOCT)
rule align_assembly_sho:
    """
    Function:
        1. Copy metagenome assemblies into a new directory, where they are indexed.
        2. Generate sho alignment BAM files for use with MetaBat2 and CONCOCT
        3. The MetaBat2 "depth" file is created
    Resources:
        > https://onestopdataanalysis.com/binning/
    """
    input:
        vRef     = vOutdir_metaspades + "{pond}_sho/scaffolds.fasta",
        vR1_sho  = vOutdir_metaspades + "{pond}_sho/corrected/unmapped_s_{pond}_sho_R1.fastq.00.0_0.cor.fastq.gz",
        vR2_sho  = vOutdir_metaspades + "{pond}_sho/corrected/unmapped_s_{pond}_sho_R2.fastq.00.0_0.cor.fastq.gz",
    output:
        vDir_index     = directory(vOutdir_mapping + "index/{pond}_sho/"),
        vIndex         =           vOutdir_mapping + "index/{pond}_sho/scaffolds.fasta",
        vAlignment_sho =           vOutdir_mapping + "alignment/{pond}_sho.bam",
        vDepth         =           vOutdir_mapping + "depth/{pond}_depth.txt"
    conda:
        "envs/env_1.yaml"
    threads:
        vThreads
    shell:
        """
        ## 1.
        mkdir -p {output.vDir_index}
        cp {input.vRef} {output.vIndex}
        /usr/bin/time -v bwa index {output.vIndex}

        ## 2.
        /usr/bin/time -v bwa mem -t {threads} {output.vIndex} {input.vR1_sho} {input.vR2_sho} \
            | /usr/bin/time -v samtools sort -@ {threads} -o {output.vAlignment_sho}
        /usr/bin/time -v samtools index -@ {threads} {output.vAlignment_sho}

        ## 3.
        /usr/bin/time -v jgi_summarize_bam_contig_depths --outputDepth {output.vDepth} {output.vAlignment_sho}
        """



# 03 - Binning, Generate alignment (used by Bin3c)
rule align_assembly_hic:
    """
    Function: Generate hic alignment BAM files for use with Bin3c
    """
    input:
        vRef     = vOutdir_metaspades + "{pond}_sho/scaffolds.fasta",
        vR1_hic  = vOutdir_fastp      + "{pond}/hic/{pond}_hic_R1.fastq.gz",
        vR2_hic  = vOutdir_fastp      + "{pond}/hic/{pond}_hic_R2.fastq.gz",
        vIndex   = vOutdir_mapping    + "index/{pond}_sho/scaffolds.fasta",
    output:
        vAlignment_hic = vOutdir_mapping + "alignment/{pond}_hic.bam",
    conda:
        "envs/env_1.yaml"
    threads:
        vThreads
    shell:
        """
        /usr/bin/time -v bwa mem -t {threads} -5SP {input.vIndex} {input.vR1_hic} {input.vR2_hic} \
            | /usr/bin/time -v samtools view -@ {threads} -Sb - \
            | /usr/bin/time -v samtools sort -@ {threads} -o {output.vAlignment_hic} -
        """
            


## 03 - Binning, MetaBat2
rule metabat2:
    """
    Function: 
    """
    input:
        vAssembly  = vOutdir_metaspades + "{pond}_sho/scaffolds.fasta",
        vDepth     = vOutdir_mapping + "depth/{pond}_depth.txt",
    output:
        directory(vOutdir_metabat2 + "{pond}_sho/"),
    conda:
        "envs/env_1.yaml"
    threads:
        vThreads
    shell:
        """
        /usr/bin/time -v metabat2 -t {threads} -i {input.vAssembly} -a {input.vDepth} -o {output}
        """



## 03 - Binning, CONCOCT
rule concoct:
    """
    Resources:
        > https://onestopdataanalysis.com/binning/
    """
    input:
        vScaffold  = vOutdir_metaspades + "{pond}_sho/scaffolds.fasta",
        vAlignment = vOutdir_mapping + "alignment/{pond}_sho.bam",
    output:
        directory(vOutdir_concoct + "{pond}_sho/")
    conda:
        "envs/env_3.yaml"
    threads:
        vThreads
    shell:
        """
        mkdir {output}/fasta_bins

        /usr/bin/time -v cut_up_fasta.py {input.vScaffold} -c 10000 -o 0 --merge_last -b {output}/contigs_10K.bed \
            > {output}/contigs_10K.fa

        /usr/bin/time -v concoct_coverage_table.py {output}/contigs_10K.bed {input.vAlignment} \
            > {output}/coverage_table.tsv

        /usr/bin/time -v concoct --threads {threads} --composition_file {output}/contigs_10K.fa --coverage_file {output}/coverage_table.tsv -b {output}

        /usr/bin/time -v merge_cutup_clustering.py {output}/clustering_gt1000.csv \
            > {output}/clustering_merged.csv

        /usr/bin/time -v extract_fasta_bins.py {input.vScaffold} {output}/clustering_merged.csv --output_path {output}/fasta_bins
        """



## 03 - Binning, Bin3c
rule bin3c:
    input:
        vScaffold      = vOutdir_metaspades + "{pond}_sho/scaffolds.fasta",
        vAlignment_hic = vOutdir_mapping + "alignment/{pond}_hic.bam",
    output:
        vDir_mkmap = directory(vOutdir_bin3c + "{pond}/mkmap/")
        vDir_clust = directory(vOutdir_bin3c + "{pond}/cluster/")
    conda:
        "envs/env_4.yaml"
    threads:
        vThreads
    shell:
        """
        ## Install bin3c
        git clone --recursive https://github.com/cerebis/bin3C

        ## Run Bin3c, 
        /usr/bin/time -v python ./bin3C mkmap -v {input.vScaffold} {input.vAlignment_hic} {output.vDir_mkmap}
        
        /usr/bin/time -v python2 ./bin3C cluster -v {output.vDir_mkmap}/contact_map.p.gz {output.vDir_clust}
        """