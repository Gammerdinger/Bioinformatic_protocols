# Popoolation2

## mpileup2sync

The first step to using `Popoolation2` is to convert your .mpileup file into a .sync file. To do this, you will need to:

`java -ea -Xmx110g -jar popoolation2_1201/mpileup2sync.jar --input <input.mpileup>  --output <output.sync>  --fastq-type sanger --min-qual <Minimum_quality> --threads <Threads>`

><input.mpileup> Full path to your .mpileup input file from Samtools
>
><output.sync> Full path to your .sync output file
>
>--fastq-type sanger This is the phred encoding to use. Most the the current Illumina phred scores are equivalent to sanger
>
>--min-qual <Minimum_quality> This is minimum quality to use. I usually use 20. 
>
>--threads <Threads> Number of threads to use.
