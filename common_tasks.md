## `fastq` to `fasta` conversion
```bash
#!/bin/bash
ls all_trimmed/trimmed_*.fastq | parallel -j32 seqtk seq -A {} '>' {.}.fasta
```
