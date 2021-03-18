# Teleman Lab software used for analysis of Ribosome  Footprinting Data

This sofware was used in publications from the Teleman Lab involving 40S and 80S Ribosomes Footprinting. It is written in C and can be compiled and run on any platform using a C compiler such as GCC.

This 2021 update includes all the software used for publications from the Teleman lab in 2020, plus the ones from 2021.



## Installing & Running

Download the C source files from the subdirectory /Source-Code-in-C .



Compile them with a C compiler. For instance:

```
gcc count_features_nt_range_v2.c -o count_features_nt_range_v2
```



Execute the program. Details about command line arguments and input files are below in the section "Description of the software".

```
count_features_nt_range_v2 SAM.input
```





## Description of the software



### count_features_nt_range_v2.c ###

*Description*: Takes as input a SAM file with mapped reads, and a file containing a list of transcripts including a nucleotide range for each transcript, and the software counts the number of reads that fall within the nucleotide range of each transcript. Secondary mappings are excluded.

 

*Inputs*:

•the name of the SAM file to be processed should be indicated as a command line argument

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: nucleotide range start. Third column: nucleotide range stop.

 

*Output*: a tab-delimited file with read counts per transcript. Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

 

### metagene_plot_5UTR_CDS_3UTR_v2.c ###

*Description:* Performs the counting to generate a combined metagene plot of reads in 5’UTR, ORF and 3’UTR of transcripts. Read counts were normalized for the length of each of these features to make them comparable (yielding a measurement of ‘read density’).

 

*Inputs*:

•the name of the SAM file to be processed should be indicated as a command line argument

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: total length of the transcript in nucleotides. Third column: position of the first nucleotide of the start codon of the main Open Reading Frame on the transcript. Fourth column: position of the last nucleotide of the stop codon of the main Open Reading Frame on the transcript.

 

*Output*: a tab-delimited file with read counts per position, first for 5’UTRs, then ORFs, then 3’UTRs. Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

 

### metagene_plot_interval_v6.c ###

*Description:* Performs the counting to generate a 2-dimensional metagene plot of reads within an interval, with position in that interval normalized to 0-100% (e.g. all transcript 5’UTRs), resolved for footprint length (1nt-100nt). Secondary mappings are not counted to avoid biasing genes with many transcript isoforms or reads with low sequence complexity. Only one interval per transcript is allowed.

 

*Inputs*:

•SAM file entitled “input.sam” to be processed

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: interval start position in nucleotides. Third column: interval end position in nucleotides.

 

*Output*: a tab-delimited file ‘output.txt’ with read counts per position percentile (in rows), resolved for footprint length (in columns). Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

 

### metagene_plot_onePoint_v7.c ###

*Description:* Performs the counting to generate a 2-dimensional metagene plot of reads relative to fixed point (e.g. start codon), from 100nt upstream of the fixed point to 100nt downstream, resolved for footprint length (1nt-100nt). Secondary mappings are not counted to avoid biasing genes with many transcript isoforms or reads with low sequence complexity. Only one interval per transcript is allowed.

 

*Inputs*:

•SAM file entitled “input.sam” to be processed

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: position of the fixed point, in nucleotides, for that transcript.

 

*Output*: a tab-delimited file ‘output.txt’ with read counts per position (in rows), resolved for footprint length (in columns). Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

 

### single_gene_plot_v2.c ###

*Description:* Counts the number of reads as a function of position on a single transcript, resolved for footprint length (1nt-100nt). Secondary mappings are counted.

 

*Inputs*:

•the name of the SAM file to be processed should be indicated as the first command line argument

•the name of the transcript to be analysed should be indicated as the second command line argument

•a tab-delimited file called “in_features.txt” containing information for all transcripts, as for the program “metagene_plot_5UTR_CDS_3UTR_v2.c”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: total length of the transcript in nucleotides. Third column: position of the first nucleotide of the start codon of the main Open Reading Frame on the transcript. Fourth column: position of the last nucleotide of the stop codon of the main Open Reading Frame on the transcript.

 

*Output*: a tab-delimited file ‘output.txt’ with read counts per position (in rows), resolved for footprint length (in columns). Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.



### transcripts_with_artefacts_v4_20_buckets.c ###

*Description:* Takes as input two files - a SAM alignment file with footprint reads mapped to the transcriptome, and a list of 'features' with nucleotide coordinates for each transcript (for instance the location of the start and stop codons on the transcript). This program with divide each feature into 20 consecutive pieces ("buckets") and will count the number of reads in each bucket for each feature. If more than 80% of the reads fall into 2 buckets, it will flag the transcript as having a problem (in the last column of the output file). 

This program can be used for different purposes: 1) to generate a low-resolution map of reads per 'bucket', or 2) to flag transcripts that might have PCR amplification artefacts because the reads are not distributed uniformly throughout the entire Open Reading Frame, but rather cluster entirely into one or two regions. 



*Inputs*:

•the name of the SAM file to be processed should be indicated as the first command line argument

•a tab-delimited file called “in_features.txt” containing information for all transcripts, as for the program “metagene_plot_5UTR_CDS_3UTR_v2.c”. It should have no header row. First column: transcript name (to match the transcript name in the SAM file). Second column: position where the feature starts on the transcript (e.g. for an ORF, the first nucleotide of the start codon). Third column: position where the feature ends on the transcipt (e.g. for an ORF, the 3rd nucleotide of the stop codon).

 

*Output*: a tab-delimited file called ‘output_[input file name].txt’ with read counts. Each row is a different transcript as listed in the in_features.txt file. Each column is the number of reads per position bucket. In  the last column, if >80% of all reads in that transcript map to 2 of the 20 buckets, the transcript is flagged as potentially containing a PCR artefact.  Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed in the SAM file.

 

## Authors

* **Aurelio Teleman** - a.teleman@dkfz.de



## License

This project is licensed under the GNU General Public License, version 3 (GPLv3) - see the [LICENSE](LICENSE) file for details

