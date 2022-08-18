# Hapbam-G:

### Input:
The input for Hapbam-G should be 4 .bam files. Using <a href="https://github.com/fenderglass/hapdup">hapdup</a> two phased .fasta files should be mapped to a diploid reference sequence using <a href="https://github.com/lh3/minimap2">minimap2</a>. The output of the sequence alignment map format files should then be converted to .bam using <a href="https://github.com/samtools/samtools">samtools</a>. When generating .sam files in <a href="https://github.com/lh3/minimap2">minimap2</a>, use the --eqx flag to distinguish between mismatches and matches in the alignment file. 

<img src="https://i.postimg.cc/kXkYFvTq/Hapdup-Phased-to-Ref-1.jpg" alt="mapping" style="width:500px;height:auto;">


### General Usage:

Input Format: 
```
Maternal Phased 1.bam Maternal Phased 2.bam Paternal Phased 1.bam Paternal Phased 2.bam -g -flag
```

Flags: 

#### -g: generates a .csv file that you can use for a graphing comparison of each phaseblock.
#### -mm: calculates the mismatch frequency of each phaseblock
#### -ma: calculates the match frequency of each phaseblock 

Frequency Calculation:

```python
(sum_of_mismatches)/(query_sequence_length - soft_clip)
```

### Overview:
Hapbam-G parses through the CIGAR string in .bam files generated from <a href="https://github.com/fenderglass/hapdup">hapdup</a> assemblies. As <a href="https://github.com/fenderglass/hapdup">hapdup</a> does not achieve chromosomal-level phasing, and the contigs in one file might not be in sync, Hapbam-G is a tool that determine the precision of the hapdup assemblies. Aligning each one of the phased reads and mapping the reads to each diploid reference we can determine whether the reads belong to either maternal or paternal chromosome. 
