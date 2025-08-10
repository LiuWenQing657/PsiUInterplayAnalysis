# PsiUInterplayAnalysis

This is the bioinformatics guide and scripts for PRAISE performed on small RNAs.   
For all the codes shown in this document, the parameters are the same as our data GSE299274. You can change suitable parameters for yourself.  
The original bioinformatics pipeline of PRAISE is in https://github.com/Zhe-jiang/PRAISE.  

## 1. Requirements
### Softwares:
cutadapt  
umi_tools  
seqkit  
 hisat2  
samtools  
python3  
packages of python3  
pysam  

## 2. Bioinformatics pipeline of PRAISE performed on small RNAs.
The eCLIP Library is suitable for quantitative purposes of RNA with short length (e.g. tRNA sequencing).
### get cleaned reads
Cut adapters
>cutadapt -j {CORES} --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 30 \
>-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
>-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
>-o {output.fix_R1} \
>-p {output.fix_R2} \
> {input.raw_R1} \
> {input.raw_R2}

>seqkit rmdup {input.fix_R2} -s -j {CORES} -o {output.dedup_R2}

>umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNN \
>-I {input} -S {output}

### Mapping cleaned reads.
>hisat2 -p {CORES} \
>-x {INDEX_HISAT2} \
>-q --repeat --no-spliced-alignment --very-sensitive \
>-U {input.processed_R2} \
>-S {output.SAM}

>samtools sort -O BAM -o {output.bam} -@ {CORES} -m {MEM} -T {output.bam_temp} {input.sam}

>samtools index {input.bam} {output.bai}

>samtools view -h -@ {CORES} -bF 4 {input.bam_sorted} -o {output.bam_mapped}

>samtools sort -@ {CORES} -m {} -O BAM -n -o {output.bam_name_sorted} -T {output.bam_name_sorted.temp} {input.bam_mapped}

>python {Realign_script} --fast -t {CORES} -ms 4.8 \
>-x {Reference} -i {input.bam_name_sorted} -o {output.bam_realigned} -f {output.bam_filtered}

## 4. Count the BAM file and call psiU signals.

>samtools mpileup -d {depth} -BQ0 -f {REF} {input.bam} -o {output.mpileup} --ff UNMAP,QCFAIL -aa

>python {Count_script} -p {CORES} -i {input.mpileup} -o {output.bmat}

## 5. R functions for analysis on the interplay of PUS enzymes.







