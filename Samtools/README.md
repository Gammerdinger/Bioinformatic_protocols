# Samtools

## view

This is the command used to convert between your SAM and BAM files. The syntax for using it is:

`samtools view -bS <Input.sam> > <Output.bam>`

>-bS Input format is SAM and the output format will be BAM
>
><Input.sam> Full path to the input SAM file.
>
><Output.bam> Full path to the output BAM file.

## sort

You will almost always need to sort your BAM file by coordiante. To do this, you need to use:

`samtools sort <input.bam> -o <sorted.bam>`

><input.bam> Full path to the input BAM file
>
><sorted.bam> Full path to the sorted BAM file

## index

You will oftentimes need to also index your sorted BAM files. IN order to accomplish this, you will use this command:

`samtools index <sorted_input.bam> <sorted_input.bam.bai>`

><sorted_input.bam> This is the full path to your sorted input BAM file.
>
><sorted_output.bam.bai> This is the full path to your index of that BAM file. Note, the only difference here is that it has a `.bai` extension. It isn't really an input as I have labelled it, but at the same time, the name need to be the same, so I wasn't sure how to best show that. I hope this is clear.

## idxstats

You can get information on your sorted BAM file once you have indexed it using the `idxstats` command. The syntax is:

`samtools idxstats <sorted_input.bam> > <output.idxstats>`

><sorted_input.bam> This is the full path to your sorted input BAM file.
>
><output.idxstats> This is the full path to your output file with the stats concerning your index.

## faidx

Sometimes software packages ask that you index your FASTA reference file. To index a FASTA file, you will use this syntax:

`samtools faidx <Reference.fasta>`

><Reference.fasta> Full path to the reference FASTA file

## merge

If you need to merge multiple BAM files together (like you have multiple read groups from the same sample), then you would use the `merge` function. The syntax is as follows:

`samtools merge <Merged.bam> <Read_group_1.bam> <Read_group_2.bam>`

## mpileup

In order to take your BAM files and convert them to mpileup format you can use `samtools mpileup`. You can include as many bam files in your mpileup file. The syntax is as follows:

`samtools mpileup -f <Reference.fasta> <Sample_1.bam> ... <Sample_N.bam>  > <Output.mpileup>`

>-f <Reference.fasta> Full path to the reference FASTA.
>
><Sample_1.bam> ... <Sample_N.bam> Full paths to each of your BAM files.
>
><Output.mpileup> Output mpileup file.
