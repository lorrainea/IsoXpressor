IsoXpressor: A tool to assess transcriptional activity within isochores
===

<b>Description</b>: Given a genome, read directory and conditions file, IsoXpressor statistically analyses the TPM and RPKM scores of the number of reads identified in each isochore for each condition.

<b>Installation</b>: To install IsoXpressor, please follow the instructions given in file INSTALL.
```
usage: IsoXpressor.py [-h] -g GENOME -r READS_DIR -o OUTPUT_DIR -c CONDITIONS
                      [-w WINDOW_SIZE] [-s SEED_ERRORS] [-e TOTAL_ERRORS]
                      [-l SEED_LENGTH] [-a STATISTICAL_ANALYSIS] [-t THREADS]

Assessing transcriptional activity within isochores

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Genome file name
  -r READS_DIR, --reads_dir READS_DIR
                        The directory name where the read files exist.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory.
  -c CONDITIONS, --conditions CONDITIONS
                        The conditions file identifying which read belongs to
                        which condition.
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        Set window size of isochores in bp for IsoSegmenter
                        (default: 100,000).
  -s SEED_ERRORS, --seed_errors SEED_ERRORS
                        Maximum number of seed errors for REAL (default: 2).
  -e TOTAL_ERRORS, --total_errors TOTAL_ERRORS
                        Total number of errors for REAL (default: 5).
  -l SEED_LENGTH, --seed_length SEED_LENGTH
                        Length of seed in bp for REAL (default: 32).
  -a STATISTICAL_ANALYSIS, --statistical_analysis STATISTICAL_ANALYSIS
                        'TPM' (Transcripts Per Kilobase Million) or 'RPKM'
                        (Reads Per Kilobase Million) (default: 'TPM').
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 1).

```

<b>See https://github.com/lorrainea/IsoXpressor/wiki for more help.</b>

<b>License</b>: GNU GPLv3 License; Copyright (C) 2020 Stilianos Arhondakis, Lorraine A. K. Ayad, Athanasia-Maria Dourou, Solon P. Pissis


