# Hapbam-G: Getting Started:

### User Guide:

Input Format: 
```
Maternal Phased 1.bam Maternal Phased 2.bam Paternal Phased 1.bam Paternal Phased 2.bam -g -flag
```

Flags: 
-g: generates a .csv file that you can use for a graphing comparison of each phaseblock.
-mm: calculates the mismatch frequency of each phaseblock
-ma: calculates the match frequency of each phaseblock 

Frequency Calculation:

```python
sum of mismatches/(query sequence length - soft clip)
```
