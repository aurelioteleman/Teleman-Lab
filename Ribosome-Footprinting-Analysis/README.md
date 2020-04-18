# Teleman Lab software used for analysis of Ribosome  Footprinting Data

This sofware was used in publications from the Teleman Lab involving 40S and 80S Ribosomes Footprinting. It is written in C and can be compiled and run on any platform using a C compiler such as GCC.



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



**count_features_nt_range_v2.c**

*Description*: Takes as input a SAM file with mapped reads, and a file containing a list of transcripts including a nucleotide range for each transcript, and the software counts the number of reads that fall within the nucleotide range of each transcript. Secondary mappings are excluded.

 

*Inputs*:

•the name of the SAM file to be processed should be indicated as a command line argument

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: nucleotide range start. Third column: nucleotide range stop.

 

*Output*: a tab-delimited file with read counts per transcript. Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

 

**metagene_plot_5UTR_CDS_3UTR_v2.c**

*Description:* Performs the counting to generate a combined metagene plot of reads in 5’UTR, ORF and 3’UTR of transcripts. Read counts were normalized for the length of each of these features to make them comparable (yielding a measurement of ‘read density’).

 

*Inputs*:

•the name of the SAM file to be processed should be indicated as a command line argument

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: total length of the transcript in nucleotides. Third column: position of the first nucleotide of the start codon of the main Open Reading Frame on the transcript. Fourth column: position of the last nucleotide of the stop codon of the main Open Reading Frame on the transcript.

 

*Output*: a tab-delimited file with read counts per position, first for 5’UTRs, then ORFs, then 3’UTRs. Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

 

**metagene_plot_interval_v6.c**

*Description:* Performs the counting to generate a 2-dimensional metagene plot of reads within an interval, with position in that interval normalized to 0-100% (e.g. all transcript 5’UTRs), resolved for footprint length (1nt-100nt). Secondary mappings are not counted to avoid biasing genes with many transcript isoforms or reads with low sequence complexity. Only one interval per transcript is allowed.

 

*Inputs*:

•SAM file entitled “input.sam” to be processed

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: interval start position in nucleotides. Third column: interval end position in nucleotides.

 

*Output*: a tab-delimited file ‘output.txt’ with read counts per position percentile (in rows), resolved for footprint length (in columns). Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

 

**metagene_plot_onePoint_v7.c**

*Description:* Performs the counting to generate a 2-dimensional metagene plot of reads relative to fixed point (e.g. start codon), from 100nt upstream of the fixed point to 100nt downstream, resolved for footprint length (1nt-100nt). Secondary mappings are not counted to avoid biasing genes with many transcript isoforms or reads with low sequence complexity. Only one interval per transcript is allowed.

 

*Inputs*:

•SAM file entitled “input.sam” to be processed

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: position of the fixed point, in nucleotides, for that transcript.

 

*Output*: a tab-delimited file ‘output.txt’ with read counts per position (in rows), resolved for footprint length (in columns). Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

 

**single_gene_plot_v2.c**

*Description:* Counts the number of reads as a function of position on a single transcript, resolved for footprint length (1nt-100nt). Secondary mappings are counted.

 

*Inputs*:

•the name of the SAM file to be processed should be indicated as the first command line argument

•the name of the transcript to be analysed should be indicated as the second command line argument

•a tab-delimited file called “in_features.txt” containing information for all transcripts, as for the program “metagene_plot_5UTR_CDS_3UTR_v2.c”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: total length of the transcript in nucleotides. Third column: position of the first nucleotide of the start codon of the main Open Reading Frame on the transcript. Fourth column: position of the last nucleotide of the stop codon of the main Open Reading Frame on the transcript.

 

*Output*: a tab-delimited file ‘output.txt’ with read counts per position (in rows), resolved for footprint length (in columns). Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10000 lines processed.

 

## Authors

* **Aurelio Teleman** - a.teleman@dkfz.de



## License

This project is licensed under the GNU General Public License, version 3 (GPLv3) - see the [LICENSE](LICENSE) file for details

