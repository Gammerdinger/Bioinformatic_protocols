# Mpileup

In order to take your BAM files and convert them to mpileup format you can use `samtools mpileup`. You can include as many bam files in your mpileup file. The syntax is as follows:

`samtools mpileup -f <Reference.fasta> <Sample_1.bam> ... <Sample_N.bam>  > <Output.mpileup>`

>-f <Reference.fasta> Full path to the reference FASTA.
>
><Sample_1.bam> ... <Sample_N.bam> Full paths to each of your BAM files.
>
><Output.mpileup> Output mpileup file.
