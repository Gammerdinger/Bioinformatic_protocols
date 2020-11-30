# Alignment

## Short-read Alignment

### bwa

For several downstream software packages, read groups are sometimes required. As such, I generally add them as I do the alignment in bwa-mem however you can also add them in Picard. That method can be found at:

Before carrying out an alignment with `bwa` you are going to need to index your reference genome using the following command:

`bwa index <reference.fasta>`

><reference.fasta>: This is the path to the reference FASTA file.

After the indexing has run, you will find several files with your FASTA file name along different sufficies

#### mem

For most paired-end short-read alignment, I use `bwa mem`. The syntax that I generally use is:

`bwa mem -t <threads> -R '@RG\tID:<RGID>\tLB:<RGLB>\tPL:<RGPL>\tPU:<RGPU>\tSM:<RGSM>' -M <reference.fasta> <Left_reads.fastq> <Right_reads.fastq> > <Alignment.sam>`

For single-end short read-alignment, I generally use:

`bwa mem -t <threads> -R '@RG\tID:<RGID>\tLB:<RGLB>\tPL:<RGPL>\tPU:<RGPU>\tSM:<RGSM>' -M <reference.fasta> <Reads.fastq> > <Alignment.sam>`

Options that I employ:

>-t <threads>: How many threads you'd like to use to process the data.
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
>>SM: This is to mark which sample your reads are coming from. Note, this **does not** need to be unique like the ID field since you may have multiple read groups coming from a single sample
>
>-M <reference.fasta>: This is the path to the reference FASTA file. Your index file should also be located here.
>
><Left_reads.fastq> <Right_reads.fastq>/<Reads.fastq>: This is the path the the reads you would like to align to the reference genome.
>
><Alignment.sam>: This is the path and file that you would like to write the alignments to. Note that this is a SAM file which are suaully quite large.
