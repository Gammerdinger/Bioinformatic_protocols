# Alignment

## Short-read Alignment

### bwa

For several downstream software packages, read groups are sometimes required. As such, I generally add them as I do the alignment in bwa-mem however you can also add them in Picard. That method can be found at:

Before carrying out an alignment with `bwa` you are going to need to index your reference genome using the following command:

`bwa index <reference.fasta>`

The input I provide is:

><reference.fasta>: This is the path to the reference FASTA file.

After the indexing has run, you will find several files with your FASTA file name along different sufficies

#### mem

For most paired-end short-read alignment, I use `bwa mem`. The syntax that I generally use is:

`bwa mem -t <Threads> -R '@RG\tID:<RGID>\tLB:<RGLB>\tPL:<RGPL>\tPU:<RGPU>\tSM:<RGSM>' -M <reference.fasta> <Left_reads.fastq> <Right_reads.fastq> > <Alignment.sam>`

For single-end short read-alignment, I generally use:

`bwa mem -t <Threads> -R '@RG\tID:<RGID>\tLB:<RGLB>\tPL:<RGPL>\tPU:<RGPU>\tSM:<RGSM>' -M <reference.fasta> <Reads.fastq> > <Alignment.sam>`

The options and input that I provide are:

>-t \<Threads> Number of threads you'd like to use to process the data.
>  
>-R 'Read group information': There are several parts to the read group information field. Different software packages using require different ones, but ID and SM are the most common.
>>
>>ID: This **NEEDS** to be unique. This is the ID for this batch of reads.
>>
>>LB: This is not used much, but the idea is if you ran, MarkDuplicates in Picard and you had run the same DNA library on multiple lanes. I usually just use the same tag as I use for the SM tag.
>>
>>PL: This is the platform that the sequencing was run on. For aligning Illumina reads, you should use *ILLUMINA* here.
>>
>>PU: This is the platform unit and it is ideally supposed to hold <FLOWCELL_BARCODE>.<LANE>.<SAMPLE_BARCODE>, where <FLOWCELL_BARCODE> is the barcode of the flowcell, <LANE> is the lane the data was run on and <SAMPLE_BARCODE> is supposed to be a library/sample specific identifer. That being said, you may not have this information and I have not found it to ever matter. For most practices, anything can go in this field.
>>
>>SM: This is to mark which sample your reads are coming from. Note, this **does not** need to be unique like the ID field since you may have multiple read groups coming from a single sample.
>
>-M <reference.fasta> This is the path to the reference FASTA file. Your index file should also be located here.
>
><Left_reads.fastq> <Right_reads.fastq>/<Reads.fastq> This is the path the the FASTQ reads you would like to align to the reference genome.
>
><Alignment.sam> This is the path and file that you would like to write the alignments to. Note that this is a SAM file which are suaully quite large.


## bowtie2

An alternative alignment software package you may be interested in running for short-read data is `bowtie2`. Similarly to bwa, before we run `bowtie2`, we will need to index our reference genome. In order to index our reference genome with `bowtie2` we will need to run this command:

`bowtie2-build <reference.fasta> <reference_base>`

These inputs are: 
><reference.fasta> This is the path to the reference FASTA file.
>
><reference_base> This is the full path to where you want the index (Ideally in the same directory as your reference FASTA file) and what you want the index to be called.

The output from this should be several files with your reference_base followed by .bt2.

Now to run the alignment, you will need to use the following command:

`bowtie2 -x <reference_base> -1 <Left_reads.fastq> -2 <Right_reads.fastq> --very-sensitive -p <Threads> --rg-id <RGID> --rg LB:<RGLB> --rg PL:<RGPL> --rg PU:<RGPU> --rg SM:<RGSM> -S <Alignment.sam>`

>-x This is the full path to the reference base that you made during the indexing process
>
>-1 <Left_reads.fastq> The full path to the left FASTQ reads.
>
>-2 <Right_reads.fastq> The full path to the right FASTQ reads.
>
>--very-sensitive I generally align with this setting. It translates into using the following options -D 20 -R 3 -N 0 -L 20 -i S,1,0.50. 
>
>-p \<Threads> Number of threads you'd like to use to process the data.
>
>--rg-id <RGID> This is your read group ID and it **NEEDS** to be unique. This is the ID for this batch of reads. Note: I have not tested this option.
>
>--rg LB:<RGLB> This is not used much, but the idea is if you ran, MarkDuplicates in Picard and you had run the same DNA library on multiple lanes. I usually just use the same tag as I use for the SM tag. Note: I have not tested this option.
>
>--rg PL:<RGPL> This is the platform that the sequencing was run on. For aligning Illumina reads, you should use *ILLUMINA* here. Note: I have not tested this option.
>
>--rg PU:<RGPU> This is the platform unit and it is ideally supposed to hold <FLOWCELL_BARCODE>.<LANE>.<SAMPLE_BARCODE>, where <FLOWCELL_BARCODE> is the barcode of the flowcell, <LANE> is the lane the data was run on and <SAMPLE_BARCODE> is supposed to be a library/sample specific identifer. That being said, you may not have this information and I have not found it to ever matter. For most practices, anything can go in this field. Note: I have not tested this option.
>
>--rg SM:<RGSM> This is to mark which sample your reads are coming from. Note, this **does not** need to be unique like the ID field since you may have multiple read groups coming from a single sample. Note: I have not tested this option.
>
>-S <Alignment.sam> This is the path and file that you would like to write the alignments to. Note that this is a SAM file which are suaully quite large.

## HISAT2

From my understanding, the developers of `bowtie2` for the most part went on to develop `HISAT` and `HISAT2` and encourage its use. Most of its options are thus similar to `bowtie2`. Before running HISAT2, you will need to index your reference genome using the following command:

`hisat2-build <Reference.fasta> <reference_base>`

These inputs are: 
><Reference.fasta> This is the path to the reference FASTA file.
>
><reference_base> This is the full path to where you want the index (Ideally in the same directory as your reference FASTA file) and what you want the index to be called.

The output from this should be several files with your reference_base followed by .ht2.

Then you can carry the alignment using the same options as `bowtie2`:

`hisat2 -x <reference_base> -1 <Left_reads.fastq> -2 <Right_reads.fastq>-p <Threads> --rg-id <RGID> --rg LB:<RGLB> --rg PL:<RGPL> --rg PU:<RGPU> --rg SM:<RGSM> -S <Alignment.sam>`

>-x This is the full path to the reference base that you made during the indexing process
>
>-1 <Left_reads.fastq> The full path to the left FASTQ reads.
>
>-2 <Right_reads.fastq> The full path to the right FASTQ reads.
>
>-p \<Threads> Number of threads you'd like to use to process the data.
>
>--rg-id <RGID> This is your read group ID and it **NEEDS** to be unique. This is the ID for this batch of reads. Note: I have not tested this option.
>
>--rg LB:<RGLB> This is not used much, but the idea is if you ran, MarkDuplicates in Picard and you had run the same DNA library on multiple lanes. I usually just use the same tag as I use for the SM tag. Note: I have not tested this option.
>
>--rg PL:<RGPL> This is the platform that the sequencing was run on. For aligning Illumina reads, you should use *ILLUMINA* here. Note: I have not tested this option.
>
>--rg PU:<RGPU> This is the platform unit and it is ideally supposed to hold <FLOWCELL_BARCODE>.<LANE>.<SAMPLE_BARCODE>, where <FLOWCELL_BARCODE> is the barcode of the flowcell, <LANE> is the lane the data was run on and <SAMPLE_BARCODE> is supposed to be a library/sample specific identifer. That being said, you may not have this information and I have not found it to ever matter. For most practices, anything can go in this field. Note: I have not tested this option.
>
>--rg SM:<RGSM> This is to mark which sample your reads are coming from. Note, this **does not** need to be unique like the ID field since you may have multiple read groups coming from a single sample. Note: I have not tested this option.
>
>-S <Alignment.sam> This is the path and file that you would like to write the alignments to. Note that this is a SAM file which are suaully quite large.

### NextGenMap

Another alignment option is `NextGenMap`. This software package indexes the genome when you carry out the alignment. The syntax for alignment is: 

`ngm -r <Reference.fasta> -b -1 <Left_reads.fastq> -2 <Right_reads.fastq> -o <Alignment.bam> -t <Threads>`

>-r <Reference.fasta> This is the path to the reference FASTA file.
>
>-b Output in BAM format.
>
>-1 <Left_reads.fastq> Full path to the left read FASTQ files.
>
>-2 <Right_reads.fastq> Full path to the right read FASTQ files.
> 
>-o <Alignment.bam> Full path to the output BAM file.
>
>-t \<Threads> Number of threads you'd like to use.

## RNA-Seq Specific Aligners

### Kallisto

### Sailfish

## STAR

Another alignment package is `STAR`. Like the other alignment packages you first need to index the genome. You can index your genome using the following command:

`STAR --runThreadN <Threads> --runMode genomeGenerate --genomeFastaFiles <Reference.fasta> --genomeDir <Directory_for_index> --limitGenomeGenerateRAM <RAM_allocation>`

>--runThreadN \<Threads> Number of threads you'd like to use.
>
>--runMode genomeGenerate This is the option to make the index
>
>--genomeFastaFiles <Reference.fasta> This is the path to the reference FASTA file.
>
>--genomeDir <Directory_for_index> A directory to hold the index in.
>
>--limitGenomeGenerateRAM <RAM_allocation> This is in bytes. I think the only time I needed to run this I needed to increase this from the default of 31000000000 to 600000000000. This may not be needed for your genome. 

Now you can carry out the alignment with:

`STAR --runThreadN <Threads> --genomeDir <Directory_for_index> --readFilesIn <Left_reads.fastq> <Right_reads.fastq> --outSAMtype BAM SortedByCoordinate --outFileNamePrefix <Output_prefix>`

>--runThreadN \<Threads> Number of threads you'd like to use.
>
>--genomeDir <Directory_for_index> A directory with the index.
>
>--readFilesIn <Left_reads.fastq> <Right_reads.fastq> This is the path the the FASTQ reads you would like to align to the reference genome.
>
>--outSAMtype BAM This software package allows you to output directly to BAM, which is a compressed version of your SAM files. However, `Samtools` allows for SAM<=>BAM conversion very easily.
>
>SortedByCoordinate This is a nice option so that it pre-sorts your alignment output. This ia a set that is usually done after alignment in Samtools.
>
>--outFileNamePrefix <Output_prefix> The full path to a the prefix you'd like to see the output files come out with.

# Pac-Bio Long-read Alignment

## bwa

### mem

`bwa mem -t <Threads> -R '@RG\tID:<RGID>\tLB:<RGLB>\tPL:<RGPL>\tPU:<RGPU>\tSM:<RGSM>' -M <Reference.fasta> -x pacbio <Reads.fastq>  > <Alignment.sam>`

The options and input that I provide are:

>-t \<Threads> Number of threads you'd like to use to process the data.
>  
>-R 'Read group information': There are several parts to the read group information field. Different software packages using require different ones, but ID and SM are the most common.
>>
>>ID: This **NEEDS** to be unique. This is the ID for this batch of reads.
>>
>>LB: This is not used much, but the idea is if you ran, MarkDuplicates in Picard and you had run the same DNA library on multiple lanes. I usually just use the same tag as I use for the SM tag.
>>
>>PL: This is the platform that the sequencing was run on. For aligning Illumina reads, you should use *ILLUMINA* here.
>>
>>PU: This is the platform unit and it is ideally supposed to hold <FLOWCELL_BARCODE>.<LANE>.<SAMPLE_BARCODE>, where <FLOWCELL_BARCODE> is the barcode of the flowcell, <LANE> is the lane the data was run on and <SAMPLE_BARCODE> is supposed to be a library/sample specific identifer. That being said, you may not have this information and I have not found it to ever matter. For most practices, anything can go in this field.
>>
>>SM: This is to mark which sample your reads are coming from. Note, this **does not** need to be unique like the ID field since you may have multiple read groups coming from a single sample.
>
>-M <reference.fasta> This is the path to the reference FASTA file. Your index file should also be located here.
>
> -x pacbio This lets bwa mem know that you are aligning Pac-Bio reads.
><Reads.fastq> This is the path the the FASTQ reads you would like to align to the reference genome.
>
><Alignment.sam> This is the path and file that you would like to write the alignments to. Note that this is a SAM file which are suaully quite large.

