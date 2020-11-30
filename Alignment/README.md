# Alignment

## Short-read Alignment

### bwa

For several downstream software packages, read groups are sometimes required. As such, I generally add them as I do the alignment in bwa-mem however you can also add them in Picard. That method can be found at:

#### bwa mem

For most short-read alignment, I use `bwa mem`. The syntax that I generally use is:

`bwa mem -t 20 -R '@RG\tID:S_melanotheron_females\tLB:S_melanotheron_females\tPL:illumina\tPU:WJG\tSM:S_melanotheron_females' -M /lustre/cichlid-labs/reference_assemblies/O_niloticus_UMD1/O_niloticus_UMD1.fasta /lustre/wgammerd/fastq_files/S_melanotheron/S_melanotheronFemales_R1_combined.fastq /lustre/wgammerd/fastq_files/S_melanotheron/S_melanotheronFemales_R2_combined.fastq > /lustre/wgammerd/alignments/bwa/S_melanotheron_PBA/S_melanotheron_female_bwa_PBA.sam`
