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
cutadapt -j {CORES} --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 30 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o {output.fix_R1} \
-p {output.fix_R2} \
 {input.raw_R1} \
 {input.raw_R2}
 # {CORES}: core number used
 # {output.fix_R1}: Output read 1 file (.fq.gz)
 # {output.fix_R2}: Output read 2 file (.fq.gz)
 # {input.raw_R1}: Input read 1 file (.fq.gz)
 # {input.raw_R2}: Input read 2 file (.fq.gz)

seqkit rmdup {input.fix_R2} -s -j {CORES} -o {output.dedup_R2}
# {input.fix_2}: read 2 file after cutadapt (.fq.gz)
# {CORES}: core number used
# {output.dedup_R2}: read 2 file after deduplication (.fq.gz)

umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNN \
-I {input} -S {output}

### Mapping cleaned reads.
hisat2 -p {CORES} \
-x {INDEX_HISAT2} \
-q --repeat --no-spliced-alignment --very-sensitive \
-U {input.processed_R2} \
-S {output.SAM}
# {CORES}: core number used
# {INDEX_HISAT2}: index needed for hisat2
# {input.processed_R2}: pre-processed read 2 file (.fq.gz)
# {output.SAM}: output SAM file (.SAM)

samtools sort -O BAM -o {output.bam} -@ {CORES} -m {MEM} -T {output.bam_temp} {input.sam}
# {CORES}: core number used
# {MEM}: Memory used per core (e.g. 2G)
# {output.bam}: output sorted BAM file (.BAM)
# {output.bam.temp}: temperory file
# {input.sam}: input SAM file (.SAM)

samtools index {input.bam} {output.bai}
# {input.bam}: Input sorted BAM file (.BAM)
# {output.bai}: Output index file (.BAM.BAI)

samtools view -h -@ {CORES} -bF 4 {input.bam_sorted} -o {output.bam_mapped}
# {CORES}: core number used
# {input.bam_sorted}: Input sorted BAM file (.BAM)
# {output.bam_mapped}: Output BAM file without unmapped reads (.BAM)

samtools sort -@ {CORES} -m {} -O BAM -n -o {output.bam_name_sorted} -T {output.bam_name_sorted.temp} {input.bam_mapped}
# {CORES}: core number used
# {MEM}: Memory used per core (e.g. 2G)
# {output.bam_name_sorted}: output name sorted BAM file (.BAM)
# {output.bam_name_sorted.temp}: temporary file
# {input.bam_mapped}: input BAM file without unmapped reads (.BAM)

python {Realign_script} --fast -t {CORES} -ms 4.8 \
-x {Reference} -i {input.bam_name_sorted} -o {output.bam_realigned} -f {output.bam_filtered}
# {Realign_script}: use realignment_forward.py for reads aligned to the forward sequences, use realignment_reverse.py for reads aligned to the reverse sequences
# {CORES}: core number used
# {Reference}: Reference used (.fa, same reference as the hisat2 index created from)
# {input.bam_name_sorted}: Output name sorted BAM file (.BAM)
# {output.bam_realigned}: Result BAM file of the realignment (.BAM)
# {output.bam_filtered}: BAM file contains all filtered reads (.BAM)

## 4. Count the BAM file and call psiU signals.

samtools mpileup -d {depth} -BQ0 -f {REF} {input.bam} -o {output.mpileup} --ff UNMAP,QCFAIL -aa
# {depth}: max depth, 15000000 was used
# {REF}: reference, use the reference you used in mapping (.fasta)
# {input.bam}: input BAM file after mapping (.bam, need to be sorted)
# {output.mpileup}: output MPILEUP file (.mpilup)

python {Count_script} -p {CORES} -i {input.mpileup} -o {output.bmat}
# {Count_script}: use parse-mpileup.py
# {CORES}: core number used
# {input.mpileup}: input MPILEUP file (.mpileup)
# {output.bmat}: output BMAT file (.bmat)

## 5. R functions for analysis on the interplay of PUS enzymes.







