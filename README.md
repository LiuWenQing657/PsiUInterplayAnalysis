# PsiUInterplayAnalysis

This is the bioinformatics guide and scripts for PRAISE performed on small RNAs.   
For all the codes shown in this document, the parameters are the same as our data GSE299274. You can change suitable parameters for yourself.  
The original bioinformatics pipeline of PRAISE is in https://github.com/Zhe-jiang/PRAISE.  
For detailed information on the principles and experiments of PRAISE， please refer to the paper: https://www.nature.com/articles/s41589-023-01304-7.

## 1. Requirements
### Softwares:
cutadapt
umi_tools  
seqkit  
hisat2  
samtools  
Python3  
## The following Python packages are needed:  
pysam, sys, os, multiprocessing, shutil, random, string, time.

## 2. Bioinformatics pipeline of PRAISE performed on small RNAs.
Our analysis pipeline is customized for the eCLIP library construction, which is suitable for quantitative purposes of RNA with short length (e.g. tRNA sequencing).

### 1) getting cleaned reads
- Cut the adapters of raw sequencing reads
> cutadapt -j {cores} --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 30 \
> -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
> -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
> -o {output.fix_R1} \
> -p {output.fix_R2} \
> {input.raw_R1} \  
> {input.raw_R2}
  
The meanings of the parameters above are as follows:  
{cores}: core number used  
{output.fix_R1}: Output read 1 file (.fq.gz)  
{output.fix_R2}: Output read 2 file (.fq.gz)  
{input.raw_R1}: Input read 1 file (.fq.gz)  
{input.raw_R2}: Input read 2 file (.fq.gz)  

- Remove PCR duplication
> seqkit rmdup {input.fix_R2} -s -j {cores} -o {output.dedup_R2}  
  
The meanings of the parameters above are as follows:  
{cores}: core number used  
{input.fix_R2}: read 2 file after cutadapt (.fq.gz)  
{output.dedup_R2}:  read 2 file after deduplication (.fq.gz)  

- Cut 10 mer of the 5' end of the read
  + the length of UMI is 8 mer.
  + cut additional 2 mer of possible template switch to increase the confidence of the mapping result.
> umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNN \  
> -I {input} -S {output}
  
The meanings of the parameters above are as follows:  
{input.dedup_R2}: read 2 file after deduplication (.fq.gz)  
{output.cut5prime_R2}: read file after cutting 14 mer of 5' end of read (.fq.gz)   

This is the final step of reads pre-processing of eCLIP library.

### 2) Mapping cleaned reads.
- Map cleaned reads to the reference of small RNAs using hisat2 or bwa.
> hisat2 -p {cores} \
> -x {index_hisat2} \
> -q --repeat --no-spliced-alignment --very-sensitive \
> -U {input.processed_R2} \
> -S {output.sam}  

The meanings of the parameters above are as follows:  
{cores}: core number used  
{index_hisat2}: index needed for hisat2  
{input.processed_R2}: pre-processed read 2 file (.fq.gz)  
{output.sam}: output sam file (.sam)  

- Use samtools to process the alignment file to make it meet the input format for the next step.
> samtools sort -O BAM -o {output.bam} -@ {cores} -m {mem} -T {output.bam_temp} {input.sam}
> samtools index {input.bam} {output.bai}
> samtools view -h -@ {cores} -bF 4 {input.bam_sorted} -o {output.bam_mapped}
> samtools sort -@ {cores} -m {} -O BAM -n -o {output.bam_name_sorted} -T {output.bam_name_sorted.temp} {input.bam_mapped}

The meanings of the parameters above are as follows:  
{cores}: core number used  
{mem}: Memory used per core (e.g. 2G)  
{output.bam_name_sorted}: output name sorted BAM file (.bam)  
{output.bam_name_sorted.temp}: temporary file  
{input.bam_mapped}: input BAM file without unmapped reads (.bam)  

## 3. realignment of bam files.
- realign the bam files using the script realignment.py.
> python {realign_script} --fast -t {cores} -ms 4.8 \
> -x {reference} -i {input.bam_name_sorted} -o {output.bam_realigned} -f {output.bam_filtered}

The meanings of the parameters above are as follows:   
{realign_script}: use realignment_forward.py for reads aligned to the forward sequences, use realignment_reverse.py for reads aligned to the reverse sequences  
{cores}: core number used  
{reference}: Reference used (.fa, same reference as the hisat2 index created from)  
{input.bam_name_sorted}: Output name sorted bam file (.bam)  
{output.bam_realigned}: Result bam file of the realignment (.bam)  
{output.bam_filtered}: Bam file contains all filtered reads (.bam)  

- Filter multiple alignments to retain the highest confidence realignment results.
> samtools sort -O BAM -n -o {output.bam} -@ {cores} -m 2G -T {output.temp} {input.bam}
> python remove_multi_mapping.py -rm keep -i {input.bam} -o {output.bam}

The meanings of the parameters above are as follows:   
{cores}: core number used  
{input.bam}: Bam file contains multiple mapped reads
{output.bam}: Output BAM files (.bam)

- Remove the deletion signals at the end of reads, which may cause false positive signals.
> python remove_end_signal.py -t {cores} -i {input.bam} -o {output.bam}
> samtools sort -O BAM -o {output.bam} -@ {cores} -m 2G -T {output.temp} {input.bam}

The meanings of the parameters above are as follows:  
{cores}: core number used  
{input.bam}: Bam file contains the deletion signals at the end of reads
{output.bam}: Output BAM files (.bam)

## 4. Count the BAM file and call psiU signals.

- Count MPILEUP file, and make it into a BMAT file using the script parse-mpileup_2.py. The script here originally from [Howard MENG](https://github.com/MengHoward). To ensure code clarity, it is provided directly in the repository.  
> samtools mpileup -d {depth} -BQ0 -f {ref} {input.bam} -o {output.mpileup} --ff UNMAP,QCFAIL -aa  
> python parse-mpileup_2.py -p {cores} -i {input.mpileup} -o {output.bmat}  

The meanings of the parameters above are as follows:   
{count_script}: use parse-mpileup.py  
{cores}: core number used  
{input.mpileup}: input MPILEUP file (.mpileup)  
{output.bmat}: output BMAT file (.bmat)  

- Call pseudouridine signals.
> python call_signal.py -pc {p_value} -u {input.untreated_bmat} -i {input.treated_bmat} -o {output.signals_csv} -o {output.signals_bed}

The meanings of the parameters above are as follows:  
{p_value}: The p-value, which is used for pseudouridine signal calls. Here we use 0.0001.
{input.untreated_bmat}: input BMAT file of the untreated group (.bmat)  
{input.treated_bmat}: input BMAT file of the treated group (.bmat)  
{output.signals_csv}: output csv file (.csv)  
{output.signals_bed}: output bed file (.bed)  

## 5. R functions for analysis on the interplay of PUS enzymes.
After the pseudouridine signals are detected, we use R scripts for downstream data analysis and graphic processing.






