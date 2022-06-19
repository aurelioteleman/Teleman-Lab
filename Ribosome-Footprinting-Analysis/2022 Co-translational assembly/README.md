# Teleman Lab software used for detecting co-translational assembly events in ribosome footprinting data

This sofware is written in C and can be compiled and run on any platform using a C compiler such as GCC.



## Installing & Running

Download the C source files from the subdirectory /Source_code_in_C .



Compile them with a C compiler. For instance:

```unix
gcc extract_density_curves_all_transcripts_from_sam_v1.c -o extract_density_curves_all_transcripts_from_sam_v1
```



Execute the program. Details about command line arguments and input files are below in the section "Description of the software" below.

```unix
extract_density_curves_all_transcripts_from_sam_v1 input.SAM
```





## General outline of analysis for translational co-assembly events

1. Perform total 80S footprinting and selective 80S footprinting for your factor of interest, and sequence the reads.
2. Trim the reads and align them to the transcriptome. This software is not included here. The transcriptome should consist of a FASTA file with all transcripts in the genome.
3. Run the software "extract_density_curves_all_transcripts_from_sam_v1" (included in this package) separately on both the total-80S and selective-80S datasets. This software will divide up the CDS for each transcript into 100 buckets and count the number of reads in each bucket. You will get two output files: one with counts for the selective-80S dataset, and one for the total-80S dataset.
4. Use the excel sheet provided in the "Sample_input_files" folder called "excel_template.xlsx" to calculate the ratio of selective-80S / total-80S for each bucket for each transcript. Paste the total-80S data into the "total" sheet and the selective-80S data into the "selective" sheet. The "ratio" sheet contains the formula to calculate the ratio. Copy and paste it down to cover your entire dataset. Note that if the total footprints are 0, the ratio is set to 1. Then save the "ratio" sheet as tab-delimited text file called "input.txt" (sample file provided in "Sample_input_files" folder). Alternatively, use any other method to calculate this ratio and save it to a tab-delimited"input.txt" file where the first column is the transcript name, and subsequent columns are the ratio of selective-80S / total-80S for each position bucket.
5. Run the software "fit_2step_v3_remove_two_outliers" on the input file. The output will contain information about the best-fitting step function, including the size of the step (the larger it is, the more clear the translational co-assembly) and how well the data fit to the step function (sum of the squares of the residuals). 



## Description of the software



### extract_density_curves_all_transcripts_from_sam_v1.c ###

*Description*: Takes as input a SAM file containing reads that have been mapped to a transcriptome (not genome), and a file called "in_features.txt" containing a list of transcripts and a nucleotide range for each transcript. The software counts the number of reads that fall within the nucleotide range of each transcript.  The nucleotide range is divided up into a fixed number of buckets for all transcripts. This can be set on line 6 of the program, but by default is 100 (ie percentiles). For analyzing translational co-assembly events, this should be done once for a total-80S library and once for a selective-80S library.

 

*Inputs*:

•the name of the SAM file to be processed should be indicated as a command line argument

•a tab-delimited file called “in_features.txt”. No header row. First column: transcript name (to match the transcript name in the SAM file). Second column: nucleotide range start. Third column: nucleotide range stop. For analyzing translational co-assembly events, the nucleotide range given here should be the position of the main Coding Sequence / ORF (ie the 2nd column should be the position of the first nucleotide of the start codon; the 3rd column should be the position of the last nucleotide of the stop codon).

 

*Output*: a tab-delimited file with read counts per transcript. First column: transcript name. Next 100 columns: number of reads for that transcript and percentile bucket. Note, if you change the number of buckets that each transcript is divided into, the number of columns will change. Progress of the program is shown to the user by displaying a dot (‘.’) on the screen every 10.000 lines processed.

 

 

### fit_2step_v3_remove_two_outliers.c ###

*Description:* Returns a co-translational assembly score for each transcript based on the output of the previous program "extract_density_curves_all_transcripts_from_sam_v1"

 

*Input:* a file called "input.txt". This file should be tab-delimited. No header row. First column: transcript name. Next 100 columns: ratio of selective-80S / total-80S reads for each position bucket along the transcript CDS. A sample input file is provided in the folder "Sample_input_files". See step 4 in the section above "General outline of analysis for translational co-assembly events".

 

*Output*: a tab-delimited file called "output.txt" with the following columns:

- column 1: transcript name (from input.txt)
- columns 2-101: the data from the input file input.txt containing the ratio of selective-80S / total-80S for the transcript
- column 102: the position where the step-function breaks
- column 103: the average value of the step-function to the left of the break-point
- column 104: the average value of the step-function to the right of the break-point
- column 105: the difference between columns 103 and 104 (ie how large the 'step' is)
- column 106: the sum of the squares of the residuals between the datapoints and the step function (ie how good the fit is)

 

 

 

## Authors

* **Aurelio Teleman** - a.teleman@dkfz.de



## License

This project is licensed under the GNU General Public License, version 3 (GPLv3) - see the [LICENSE](LICENSE) file for details

