## `fastq` to `fasta` conversion
```bash
#!/bin/bash
ls all_trimmed/trimmed_*.fastq | parallel -j32 seqtk seq -A {} '>' {.}.fasta
```
## Killing a process with pid (job id)
```bash
kill <pid>
```
Then if we run `jobs`, it will show that the specific job is terminated. `jobs` can be used to see all running processes as well and we can select which job to stop.
